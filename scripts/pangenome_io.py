#!/usr/bin/env python3
"""Unified input handling for PIRATE and Panaroo pangenome outputs.

The analytical parts of PANOPTICON should consume the same objects regardless
of which pangenome program produced the source files.

Public API
----------
load_pangenome(tool, outdir, ...) -> PangenomeData
load_presence_absence(tool, outdir, ...) -> pandas.DataFrame
load_family_metadata(tool, outdir, ...) -> pandas.DataFrame
load_pangenome_summary(tool, outdir, ...) -> PangenomeSummary

Canonical orientation
---------------------
presence_absence:
    rows    = sample/genome IDs
    columns = gene-family IDs
    values  = uint8 0/1

family_metadata:
    index = gene-family ID
    standard columns include, where available:
        gene_family, gene_name, annotation, n_genomes, prevalence,
        source_tool

The loaders retain source-specific columns in family_metadata where practical.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Literal, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ToolName = Literal["pirate", "panaroo"]


@dataclass
class PangenomeSummary:
    """Tool-neutral pangenome summary."""

    tool: str
    n_genomes: Optional[int] = None
    total_families: Optional[int] = None
    core: Optional[int] = None
    soft_core: Optional[int] = None
    shell: Optional[int] = None
    cloud: Optional[int] = None
    raw: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PangenomeData:
    """Normalised pangenome inputs used by downstream PANOPTICON scripts."""

    tool: str
    outdir: Path
    presence_absence: pd.DataFrame
    family_metadata: pd.DataFrame
    summary: PangenomeSummary

    def validate(self) -> None:
        """Raise a clear error if the normalised objects are inconsistent."""
        if self.presence_absence.index.has_duplicates:
            dup = self.presence_absence.index[self.presence_absence.index.duplicated()].unique()
            raise ValueError(f"Duplicate sample IDs: {list(dup[:10])}")
        if self.presence_absence.columns.has_duplicates:
            dup = self.presence_absence.columns[
                self.presence_absence.columns.duplicated()
            ].unique()
            raise ValueError(f"Duplicate gene-family IDs: {list(dup[:10])}")
        if not set(self.presence_absence.columns).issubset(self.family_metadata.index):
            missing = [
                x for x in self.presence_absence.columns
                if x not in self.family_metadata.index
            ]
            raise ValueError(
                "Family metadata is missing presence/absence columns: "
                f"{missing[:10]}"
            )
        vals = self.presence_absence.to_numpy(copy=False)
        if vals.size and not np.isin(vals, [0, 1]).all():
            raise ValueError("Presence/absence matrix contains values other than 0 and 1")


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def _normalise_tool(tool: str) -> ToolName:
    value = str(tool).strip().lower()
    aliases = {"pirate": "pirate", "panaroo": "panaroo"}
    if value not in aliases:
        raise ValueError("--pangenome-tool must be either 'pirate' or 'panaroo'")
    return aliases[value]  # type: ignore[return-value]


def _first_existing(outdir: Path, names: Sequence[str], required: bool = True) -> Optional[Path]:
    for name in names:
        path = outdir / name
        if path.is_file() and path.stat().st_size > 0:
            return path
    if required:
        raise FileNotFoundError(
            f"None of the expected files were found in {outdir}: {', '.join(names)}"
        )
    return None


def _safe_fraction(counts: pd.Series, denominator: int) -> pd.Series:
    if denominator <= 0:
        return pd.Series(np.nan, index=counts.index, dtype=float)
    return counts.astype(float) / float(denominator)


def _frequency_bins(counts: pd.Series, n_genomes: int) -> Dict[str, int]:
    """Use standard 99/95/15% core, soft-core, shell and cloud bins."""
    if n_genomes <= 0:
        return {}
    freq = counts.astype(float) / float(n_genomes)
    return {
        "core": int((freq >= 0.99).sum()),
        "soft_core": int(((freq >= 0.95) & (freq < 0.99)).sum()),
        "shell": int(((freq >= 0.15) & (freq < 0.95)).sum()),
        "cloud": int((freq < 0.15).sum()),
    }


def _make_unique(values: Iterable[str], prefix: str = "family") -> List[str]:
    """Make IDs unique without silently discarding rows."""
    seen: Dict[str, int] = {}
    result: List[str] = []
    for i, raw in enumerate(values, start=1):
        base = str(raw).strip() or f"{prefix}_{i}"
        n = seen.get(base, 0)
        seen[base] = n + 1
        result.append(base if n == 0 else f"{base}__dup{n + 1}")
    return result


# ---------------------------------------------------------------------------
# Panaroo
# ---------------------------------------------------------------------------

PANAROO_ANNOTATION_COLUMNS = ("Gene", "Non-unique Gene name", "Annotation")


def read_panaroo_rtab(path: Path) -> pd.DataFrame:
    """Read Panaroo gene_presence_absence.Rtab.

    Panaroo writes rows as gene clusters and columns as genomes. This function
    returns the canonical transpose: genomes x gene clusters.
    """
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if raw.empty or raw.shape[1] < 2:
        raise ValueError(f"Panaroo Rtab is empty or malformed: {path}")

    family_col = raw.columns[0]
    families = _make_unique(raw[family_col].astype(str), prefix="panaroo_family")
    matrix = raw.drop(columns=[family_col]).apply(pd.to_numeric, errors="coerce")
    matrix = matrix.fillna(0).ne(0).astype(np.uint8)
    matrix.index = families
    matrix.index.name = "gene_family"
    out = matrix.T
    out.index = out.index.astype(str)
    out.index.name = "sample"
    out.columns.name = "gene_family"
    return out


def _panaroo_csv_header(path: Path) -> Tuple[List[str], List[str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        header = next(csv.reader(fh))
    annotation_cols = [c for c in PANAROO_ANNOTATION_COLUMNS if c in header]
    sample_cols = [c for c in header if c not in annotation_cols]
    return annotation_cols, sample_cols


def read_panaroo_csv_presence_absence(
    path: Path,
    chunksize: int = 500,
) -> pd.DataFrame:
    """Build a binary matrix directly from gene_presence_absence.csv.

    This is the fallback when gene_presence_absence.Rtab is unavailable. The
    file can be hundreds of MB, so it is processed in row chunks.
    """
    annotation_cols, sample_cols = _panaroo_csv_header(path)
    if "Gene" not in annotation_cols:
        raise ValueError(f"Panaroo CSV has no 'Gene' column: {path}")
    if not sample_cols:
        raise ValueError(f"Panaroo CSV has no genome columns: {path}")

    blocks: List[pd.DataFrame] = []
    family_ids: List[str] = []
    usecols = ["Gene", *sample_cols]
    for chunk in pd.read_csv(
        path,
        usecols=usecols,
        dtype=str,
        keep_default_na=False,
        chunksize=chunksize,
        low_memory=False,
    ):
        family_ids.extend(chunk["Gene"].astype(str).tolist())
        present = chunk[sample_cols].ne("").astype(np.uint8)
        blocks.append(present)

    if not blocks:
        raise ValueError(f"No gene-family rows found in {path}")
    matrix = pd.concat(blocks, axis=0, ignore_index=True)
    matrix.index = _make_unique(family_ids, prefix="panaroo_family")
    matrix.index.name = "gene_family"
    out = matrix.T
    out.index = out.index.astype(str)
    out.index.name = "sample"
    out.columns.name = "gene_family"
    return out


def read_panaroo_family_metadata(
    csv_path: Path,
    presence_absence: pd.DataFrame,
) -> pd.DataFrame:
    """Read Panaroo family names/annotations without loading genome cells."""
    meta = pd.read_csv(
        csv_path,
        usecols=lambda c: c in PANAROO_ANNOTATION_COLUMNS,
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )
    if "Gene" not in meta.columns:
        raise ValueError(f"Panaroo CSV has no 'Gene' column: {csv_path}")

    ids = _make_unique(meta["Gene"].astype(str), prefix="panaroo_family")
    meta = meta.copy()
    meta["gene_family"] = ids
    meta = meta.rename(
        columns={
            "Non-unique Gene name": "gene_name",
            "Annotation": "annotation",
        }
    )
    meta = meta.set_index("gene_family", drop=False)

    counts = presence_absence.sum(axis=0).astype(int)
    meta = meta.reindex(presence_absence.columns)
    meta["gene_family"] = meta.index
    meta["n_genomes"] = counts.reindex(meta.index).fillna(0).astype(int)
    meta["prevalence"] = _safe_fraction(meta["n_genomes"], presence_absence.shape[0])
    meta["source_tool"] = "panaroo"
    return meta



def make_minimal_panaroo_family_metadata(
    presence_absence: pd.DataFrame,
) -> pd.DataFrame:
    """Create family metadata when Panaroo CSV annotations are unavailable."""
    counts = presence_absence.sum(axis=0).astype(int)
    meta = pd.DataFrame(index=presence_absence.columns.copy())
    meta.index.name = "gene_family"
    meta["gene_family"] = meta.index.astype(str)
    meta["gene_name"] = pd.NA
    meta["annotation"] = pd.NA
    meta["n_genomes"] = counts.reindex(meta.index).fillna(0).astype(int)
    meta["prevalence"] = _safe_fraction(
        meta["n_genomes"],
        presence_absence.shape[0],
    )
    meta["source_tool"] = "panaroo"
    return meta


def parse_panaroo_summary(
    path: Optional[Path],
    presence_absence: pd.DataFrame,
) -> PangenomeSummary:
    counts = presence_absence.sum(axis=0).astype(int)
    bins = _frequency_bins(counts, presence_absence.shape[0])
    raw: Dict[str, Any] = {}

    if path is not None:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if not line:
                continue
            if "\t" in line:
                key, value = line.split("\t", 1)
            elif ":" in line:
                key, value = line.split(":", 1)
            else:
                continue
            raw[key.strip()] = value.strip()

    return PangenomeSummary(
        tool="panaroo",
        n_genomes=presence_absence.shape[0],
        total_families=presence_absence.shape[1],
        core=bins.get("core"),
        soft_core=bins.get("soft_core"),
        shell=bins.get("shell"),
        cloud=bins.get("cloud"),
        raw=raw,
    )


def load_panaroo(
    outdir: Path,
    prefer_rtab: bool = True,
    csv_chunksize: int = 500,
) -> PangenomeData:
    """Load Panaroo output, allowing Rtab-only partial output directories."""
    rtab_path = _first_existing(
        outdir,
        ["gene_presence_absence.Rtab", "gene_presence_absence.rtab"],
        required=False,
    )
    csv_path = _first_existing(
        outdir,
        ["gene_presence_absence.csv"],
        required=False,
    )

    if prefer_rtab and rtab_path is not None:
        pa = read_panaroo_rtab(rtab_path)
        matrix_source = rtab_path.name
    elif csv_path is not None:
        pa = read_panaroo_csv_presence_absence(
            csv_path,
            chunksize=csv_chunksize,
        )
        matrix_source = csv_path.name
    elif rtab_path is not None:
        pa = read_panaroo_rtab(rtab_path)
        matrix_source = rtab_path.name
        print(
            "WARNING: gene_presence_absence.csv is unavailable or empty; "
            "falling back to gene_presence_absence.Rtab.",
            file=sys.stderr,
        )
    else:
        raise FileNotFoundError(
            "No usable Panaroo matrix found. Expected a non-empty "
            "gene_presence_absence.Rtab or gene_presence_absence.csv "
            f"in {outdir}."
        )

    if csv_path is not None:
        try:
            metadata = read_panaroo_family_metadata(csv_path, pa)
            annotation_source = csv_path.name
        except (ValueError, pd.errors.EmptyDataError, OSError) as exc:
            print(
                "WARNING: Panaroo CSV annotations could not be read "
                f"({exc}). Continuing with minimal family metadata.",
                file=sys.stderr,
            )
            metadata = make_minimal_panaroo_family_metadata(pa)
            annotation_source = "minimal_from_presence_absence"
    else:
        metadata = make_minimal_panaroo_family_metadata(pa)
        annotation_source = "minimal_from_presence_absence"
        print(
            "WARNING: gene_presence_absence.csv is unavailable or empty. "
            "Gene names and annotations will be omitted.",
            file=sys.stderr,
        )

    summary_path = _first_existing(
        outdir,
        ["summary_statistics.txt"],
        required=False,
    )
    summary = parse_panaroo_summary(summary_path, pa)
    summary.raw["matrix_source"] = matrix_source
    summary.raw["annotation_source"] = annotation_source

    data = PangenomeData("panaroo", outdir, pa, metadata, summary)
    data.validate()
    return data


# ---------------------------------------------------------------------------
# PIRATE
# ---------------------------------------------------------------------------


def _read_fasta_records(path: Path) -> Tuple[List[str], List[str]]:
    ids: List[str] = []
    seqs: List[str] = []
    current_id: Optional[str] = None
    current: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    ids.append(current_id)
                    seqs.append("".join(current))
                current_id = line[1:].split()[0]
                current = []
            else:
                current.append(line)
    if current_id is not None:
        ids.append(current_id)
        seqs.append("".join(current))
    if not ids:
        raise ValueError(f"No FASTA records found: {path}")
    return ids, seqs


def read_pirate_binary_fasta(path: Path) -> pd.DataFrame:
    """Read PIRATE binary_presence_absence.fasta.

    Supports A/C (A=absent, C=present) and 0/1 encodings. PIRATE normally
    writes one record per genome.
    """
    ids, seqs = _read_fasta_records(path)
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"PIRATE binary FASTA is not rectangular: {path}")

    def decode(seq: str) -> np.ndarray:
        return np.fromiter(
            (1 if char.upper() in {"C", "1"} else 0 for char in seq),
            dtype=np.uint8,
            count=len(seq),
        )

    matrix = np.vstack([decode(s) for s in seqs])
    columns = [f"gf_{i + 1}" for i in range(matrix.shape[1])]
    out = pd.DataFrame(matrix, index=ids, columns=columns, dtype=np.uint8)
    out.index.name = "sample"
    out.columns.name = "gene_family"
    return out


def _pirate_family_id_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "gene_family",
        "gene_name",
        "family",
        "cluster",
        "allele_name",
        "g_name",
    ]
    return next((c for c in candidates if c in df.columns), None)


def _pirate_count_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "number_genomes",
        "number_genome",
        "No. isolates",
        "No. isolates ",
        "No_isolates",
        "No. genomes",
        "No. Genomes",
    ]
    return next((c for c in candidates if c in df.columns), None)


def read_pirate_family_metadata(
    path: Path,
    presence_absence: pd.DataFrame,
    threshold: Optional[float] = 95,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read PIRATE family table and attach IDs to binary FASTA columns.

    PIRATE versions and settings produce slightly different tables. If a
    multi-threshold table is supplied, the requested threshold is selected
    when possible. If the selected row count matches the binary matrix width,
    its family IDs replace generic gf_1...gf_n labels.
    """
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    selected = table.copy()
    if threshold is not None and "threshold" in selected.columns:
        numeric = pd.to_numeric(selected["threshold"], errors="coerce")
        mask = numeric.eq(float(threshold))
        if mask.any():
            selected = selected.loc[mask].copy()

    id_col = _pirate_family_id_column(selected)
    if id_col is not None:
        ids = _make_unique(selected[id_col].astype(str), prefix="pirate_family")
    else:
        ids = [f"gf_{i + 1}" for i in range(len(selected))]

    pa = presence_absence.copy()
    if len(ids) == pa.shape[1]:
        pa.columns = ids
        pa.columns.name = "gene_family"
    else:
        # Keep safe positional IDs rather than risking a wrong mapping.
        ids = list(pa.columns)
        selected = selected.iloc[: len(ids)].copy()
        if len(selected) < len(ids):
            selected = selected.reindex(range(len(ids)))

    selected = selected.copy()
    selected["gene_family"] = ids
    selected = selected.set_index("gene_family", drop=False)
    selected = selected.reindex(pa.columns)

    counts = pa.sum(axis=0).astype(int)
    count_col = _pirate_count_column(selected)
    if count_col is not None:
        parsed_counts = pd.to_numeric(selected[count_col], errors="coerce")
        selected["n_genomes_source"] = parsed_counts
    selected["n_genomes"] = counts.reindex(selected.index).fillna(0).astype(int)
    selected["prevalence"] = _safe_fraction(selected["n_genomes"], pa.shape[0])
    selected["source_tool"] = "pirate"

    # Tool-neutral aliases where the source table provides likely fields.
    for candidate in ["gene_product", "product", "annotation"]:
        if candidate in selected.columns:
            selected["annotation"] = selected[candidate]
            break
    if "gene_name" in selected.columns:
        selected["gene_name"] = selected["gene_name"]

    return pa, selected


def parse_pirate_summary(
    path: Optional[Path],
    presence_absence: pd.DataFrame,
) -> PangenomeSummary:
    counts = presence_absence.sum(axis=0).astype(int)
    bins = _frequency_bins(counts, presence_absence.shape[0])
    raw: Dict[str, Any] = {}

    if path is not None:
        text = path.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if "\t" in line:
                key, value = line.split("\t", 1)
                raw[key.strip()] = value.strip()
            elif ":" in line:
                key, value = line.split(":", 1)
                raw[key.strip()] = value.strip()

    return PangenomeSummary(
        tool="pirate",
        n_genomes=presence_absence.shape[0],
        total_families=presence_absence.shape[1],
        core=bins.get("core"),
        soft_core=bins.get("soft_core"),
        shell=bins.get("shell"),
        cloud=bins.get("cloud"),
        raw=raw,
    )


def load_pirate(
    outdir: Path,
    threshold: Optional[float] = 95,
) -> PangenomeData:
    binary_path = _first_existing(
        outdir,
        ["binary_presence_absence.fasta", "binary_presence_absence.fa"],
    )
    family_path = _first_existing(
        outdir,
        [
            "PIRATE.gene_families.ordered.tsv",
            "PIRATE.gene_families.tsv",
            "PIRATE.gene_families.tab",
        ],
    )
    pa = read_pirate_binary_fasta(binary_path)
    pa, metadata = read_pirate_family_metadata(family_path, pa, threshold=threshold)
    summary_path = _first_existing(
        outdir,
        ["PIRATE.pangenome_summary.txt", "pangenome_summary.txt"],
        required=False,
    )
    summary = parse_pirate_summary(summary_path, pa)
    data = PangenomeData("pirate", outdir, pa, metadata, summary)
    data.validate()
    return data


# ---------------------------------------------------------------------------
# Public dispatch functions
# ---------------------------------------------------------------------------


def load_pangenome(
    tool: str,
    outdir: str | os.PathLike[str],
    *,
    pirate_threshold: Optional[float] = 95,
    prefer_panaroo_rtab: bool = True,
    panaroo_from_csv: Optional[bool] = None,
    panaroo_csv_chunksize: int = 500,
) -> PangenomeData:
    """Load and normalise a PIRATE or Panaroo output directory."""
    normalised = _normalise_tool(tool)

    if panaroo_from_csv is not None:
        prefer_panaroo_rtab = not bool(panaroo_from_csv)

    root = Path(outdir).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Pangenome output directory does not exist: {root}")

    if normalised == "panaroo":
        return load_panaroo(
            root,
            prefer_rtab=prefer_panaroo_rtab,
            csv_chunksize=panaroo_csv_chunksize,
        )
    return load_pirate(root, threshold=pirate_threshold)


def load_presence_absence(
    tool: str,
    outdir: str | os.PathLike[str],
    **kwargs: Any,
) -> pd.DataFrame:
    return load_pangenome(tool, outdir, **kwargs).presence_absence


def load_family_metadata(
    tool: str,
    outdir: str | os.PathLike[str],
    **kwargs: Any,
) -> pd.DataFrame:
    return load_pangenome(tool, outdir, **kwargs).family_metadata


def load_pangenome_summary(
    tool: str,
    outdir: str | os.PathLike[str],
    **kwargs: Any,
) -> PangenomeSummary:
    return load_pangenome(tool, outdir, **kwargs).summary


# ---------------------------------------------------------------------------
# Small command-line inspector
# ---------------------------------------------------------------------------


def _main() -> None:
    parser = argparse.ArgumentParser(
        description="Load PIRATE/Panaroo outputs and report their normalised dimensions."
    )
    parser.add_argument("--pangenome-tool", required=True, choices=["pirate", "panaroo"])
    parser.add_argument("--pangenome-out", required=True)
    parser.add_argument("--pirate-threshold", type=float, default=95)
    parser.add_argument(
        "--panaroo-from-csv",
        action="store_true",
        help="Ignore gene_presence_absence.Rtab and derive the matrix from the CSV.",
    )
    parser.add_argument(
        "--write-prefix",
        help="Optionally write normalised .presence_absence.tsv and .family_metadata.tsv files.",
    )
    args = parser.parse_args()

    data = load_pangenome(
        args.pangenome_tool,
        args.pangenome_out,
        pirate_threshold=args.pirate_threshold,
        prefer_panaroo_rtab=not args.panaroo_from_csv,
    )
    print(f"tool: {data.tool}")
    print(f"samples: {data.presence_absence.shape[0]}")
    print(f"gene families: {data.presence_absence.shape[1]}")
    print(f"core (>=99%): {data.summary.core}")
    print(f"soft-core (95-99%): {data.summary.soft_core}")
    print(f"shell (15-95%): {data.summary.shell}")
    print(f"cloud (<15%): {data.summary.cloud}")

    if args.write_prefix:
        prefix = Path(args.write_prefix)
        data.presence_absence.to_csv(
            str(prefix) + ".presence_absence.tsv", sep="\t", index=True
        )
        data.family_metadata.to_csv(
            str(prefix) + ".family_metadata.tsv", sep="\t", index=False
        )


if __name__ == "__main__":
    _main()
