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
    """Read PIRATE ``binary_presence_absence.fasta``.

    PIRATE A/C encoding uses ``A = present`` and ``C = absent``. 0/1 input is
    also supported with ``1 = present``. Records are genomes and sequence
    positions are anonymous family positions.

    This function is retained as a backwards-compatible fallback. The primary
    PIRATE loader uses the per-genome columns in ``PIRATE.gene_families*.tsv``
    because those columns carry explicit family IDs and can be validated
    against ``number_genomes``.
    """
    ids, seqs = _read_fasta_records(path)
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"PIRATE binary FASTA is not rectangular: {path}")

    alphabet = set("".join(seqs).upper())
    if alphabet.issubset({"A", "C"}):
        presence = {"A"}
    elif alphabet.issubset({"0", "1"}):
        presence = {"1"}
    else:
        raise ValueError(
            "Unsupported PIRATE binary FASTA alphabet "
            f"{sorted(alphabet)} in {path}; expected A/C or 0/1."
        )

    matrix = np.vstack([
        np.fromiter(
            (1 if char.upper() in presence else 0 for char in seq),
            dtype=np.uint8,
            count=len(seq),
        )
        for seq in seqs
    ])
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



def read_pirate_family_table_matrix(
    path: Path,
    sample_ids: Optional[Sequence[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Build the canonical PIRATE matrix directly from the family table.

    PIRATE ``PIRATE.gene_families.ordered.tsv`` contains one row per gene
    family and one column per genome. A non-empty genome cell denotes presence.

    Using this table avoids two unsafe assumptions about
    ``binary_presence_absence.fasta``:
      1. A/C polarity (PIRATE uses A=present, C=absent).
      2. Positional family order, which need not match the ordered family TSV.

    When ``sample_ids`` are provided (normally from the binary FASTA headers),
    those exact columns are used and their order is preserved. Otherwise,
    columns following ``cluster_order`` are treated as genome columns.

    The resulting row sums are validated against PIRATE ``number_genomes`` when
    that field is available. Any discrepancy is treated as a hard error.
    """
    table = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )
    if table.empty:
        raise ValueError(f"PIRATE family table is empty: {path}")

    id_col = _pirate_family_id_column(table)
    if id_col is None:
        raise ValueError(
            f"Could not identify a PIRATE family-ID column in {path}. "
            f"Found: {list(table.columns[:25])}"
        )

    requested = [str(x) for x in sample_ids] if sample_ids is not None else []
    if requested:
        missing = [x for x in requested if x not in table.columns]
        if missing:
            raise ValueError(
                "PIRATE family table is missing genome columns found in the "
                f"binary FASTA, e.g. {missing[:10]}"
            )
        sample_cols = requested
    elif "cluster_order" in table.columns:
        start = table.columns.get_loc("cluster_order") + 1
        sample_cols = [str(c) for c in table.columns[start:]]
    else:
        raise ValueError(
            "Could not infer PIRATE genome columns. Supply a binary FASTA with "
            "matching genome IDs or use a family table containing cluster_order."
        )

    if not sample_cols:
        raise ValueError(f"No PIRATE genome columns detected in {path}")

    family_ids = _make_unique(table[id_col].astype(str), prefix="pirate_family")

    # PIRATE genome cells contain locus/allele identifiers when present and are
    # blank when absent.
    family_by_sample = table[sample_cols].apply(
        lambda col: col.astype(str).str.strip().ne("")
    ).astype(np.uint8)
    family_by_sample.index = family_ids
    family_by_sample.index.name = "gene_family"

    observed = family_by_sample.sum(axis=1).astype(int)

    count_col = _pirate_count_column(table)
    source_counts: Optional[pd.Series] = None
    if count_col is not None:
        source_counts = pd.to_numeric(table[count_col], errors="coerce")
        comparable = source_counts.notna()
        mismatch = comparable & (
            observed.to_numpy() != source_counts.fillna(-1).astype(int).to_numpy()
        )
        if mismatch.any():
            bad_idx = np.flatnonzero(mismatch)[:10]
            examples = [
                (
                    family_ids[i],
                    int(observed.iloc[i]),
                    int(source_counts.iloc[i]),
                )
                for i in bad_idx
            ]
            raise ValueError(
                "PIRATE family-table presence/absence does not agree with "
                f"'{count_col}'. Example (family, observed, expected): {examples}"
            )

    pa = family_by_sample.T.copy()
    pa.index = pd.Index(sample_cols, name="sample")
    pa.columns.name = "gene_family"

    metadata = table.copy()
    metadata["gene_family"] = family_ids
    metadata = metadata.set_index("gene_family", drop=False)
    metadata["n_genomes"] = observed.to_numpy()
    if source_counts is not None:
        metadata["n_genomes_source"] = source_counts.to_numpy()
    metadata["prevalence"] = _safe_fraction(metadata["n_genomes"], len(sample_cols))
    metadata["source_tool"] = "pirate"

    # Tool-neutral aliases while retaining all original PIRATE fields.
    if "consensus_gene_name" in metadata.columns:
        metadata["gene_name"] = metadata["consensus_gene_name"]
    elif "gene_name" not in metadata.columns:
        metadata["gene_name"] = pd.NA

    if "consensus_product" in metadata.columns:
        metadata["annotation"] = metadata["consensus_product"]
    else:
        for candidate in ["gene_product", "product", "annotation"]:
            if candidate in metadata.columns:
                metadata["annotation"] = metadata[candidate]
                break
        else:
            metadata["annotation"] = pd.NA

    return pa, metadata

def read_pirate_family_metadata(
    path: Path,
    presence_absence: pd.DataFrame,
    threshold: Optional[float] = 95,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Backwards-compatible wrapper using explicit PIRATE genome columns.

    ``threshold`` is retained for API compatibility but is not used to filter
    ``PIRATE.gene_families.ordered.tsv``: the ordered table already represents
    the final family set, and its ``threshold`` field is family metadata rather
    than a request to select only rows equal to 95.
    """
    del threshold
    return read_pirate_family_table_matrix(
        path,
        sample_ids=[str(x) for x in presence_absence.index],
    )

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
    """Load PIRATE using the family table as the authoritative matrix source."""
    del threshold  # retained in public API for backwards compatibility

    family_path = _first_existing(
        outdir,
        [
            "PIRATE.gene_families.ordered.tsv",
            "PIRATE.gene_families.tsv",
            "PIRATE.gene_families.tab",
        ],
    )
    binary_path = _first_existing(
        outdir,
        ["binary_presence_absence.fasta", "binary_presence_absence.fa"],
        required=False,
    )

    sample_ids: Optional[List[str]] = None
    binary_encoding: Optional[str] = None
    if binary_path is not None:
        ids, seqs = _read_fasta_records(binary_path)
        sample_ids = [str(x) for x in ids]
        alphabet = set("".join(seqs).upper())
        if alphabet.issubset({"A", "C"}):
            binary_encoding = "A=present,C=absent"
        elif alphabet.issubset({"0", "1"}):
            binary_encoding = "1=present,0=absent"
        else:
            binary_encoding = f"unrecognised:{''.join(sorted(alphabet))}"

    pa, metadata = read_pirate_family_table_matrix(
        family_path,
        sample_ids=sample_ids,
    )

    summary_path = _first_existing(
        outdir,
        ["PIRATE.pangenome_summary.txt", "pangenome_summary.txt"],
        required=False,
    )
    summary = parse_pirate_summary(summary_path, pa)
    summary.raw["matrix_source"] = family_path.name
    summary.raw["family_count_validation"] = "passed"
    if binary_path is not None:
        summary.raw["binary_fasta"] = binary_path.name
        summary.raw["binary_encoding"] = binary_encoding
        summary.raw["binary_fasta_used_for_matrix"] = False

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
