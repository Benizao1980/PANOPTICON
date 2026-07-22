#!/usr/bin/env python3
"""Quantify gene-family and exact sequence-allele sharing between host groups.

PANOPTICON module
-----------------
This module uses the tool-neutral ``pangenome_io.load_pangenome`` interface and
therefore supports both Panaroo and PIRATE gene presence/absence outputs.

Main analyses
-------------
1. Observed gene-family sharing between groups.
2. Sample-size-standardised (rarefied) gene-family sharing.
3. Family-by-group prevalence and sharing classifications.
4. Pairwise heatmaps, an UpSet plot, and a host-sharing network.
5. Optional exact sequence-allele sharing from family-specific FASTA files.

Important interpretation
------------------------
Raw overlap is strongly affected by the number of genomes sampled from each
group. Rarefied metrics should normally be preferred for biological comparison.
High sharing can also reflect shared clonal lineages rather than recent HGT.
Use ``--structure-col`` to generate lineage-collapsed sensitivity analyses.

Expected canonical matrix
-------------------------
Rows = isolates; columns = gene families; values = 0/1.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

try:
    from pangenome_io import PangenomeData, load_pangenome
except ImportError as exc:
    raise SystemExit(
        "Could not import pangenome_io.py. Place it beside this script or add "
        "its directory to PYTHONPATH."
    ) from exc


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------


def ensure_outdir(path: str | os.PathLike[str]) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def set_plot_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "font.size": 11,
        "axes.titlesize": 15,
        "axes.labelsize": 12,
        "legend.fontsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def load_metadata(path: str, id_col: str) -> pd.DataFrame:
    kwargs = dict(dtype=str, encoding="utf-8", encoding_errors="replace")
    if path.lower().endswith(".csv"):
        meta = pd.read_csv(path, **kwargs)
    else:
        meta = pd.read_csv(path, sep="\t", **kwargs)

    if id_col not in meta.columns:
        raise ValueError(
            f"Metadata must contain ID column '{id_col}'. Found: {list(meta.columns)}"
        )
    if meta[id_col].duplicated().any():
        examples = meta.loc[meta[id_col].duplicated(), id_col].head(10).tolist()
        raise ValueError(f"Duplicate metadata sample IDs, e.g. {examples}")

    meta[id_col] = meta[id_col].astype(str)
    return meta.set_index(id_col)


def align_matrix_metadata(
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    minimum_samples: int = 2,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    matrix = matrix.copy()
    matrix.index = matrix.index.astype(str)
    meta = meta.copy()
    meta.index = meta.index.astype(str)

    if group_col not in meta.columns:
        raise ValueError(f"Metadata lacks group column '{group_col}'.")

    shared = matrix.index[matrix.index.isin(meta.index)]
    if len(shared) < minimum_samples:
        raise ValueError(
            f"Only {len(shared)} pangenome samples match metadata. Check IDs."
        )

    matrix = matrix.loc[shared]
    meta = meta.reindex(shared)
    valid = meta[group_col].notna() & meta[group_col].astype(str).str.strip().ne("")
    matrix = matrix.loc[valid]
    meta = meta.loc[valid]

    if matrix.empty:
        raise ValueError(f"No samples have a non-empty '{group_col}' value.")

    print(
        f"Using {matrix.shape[0]} matched samples, "
        f"{matrix.shape[1]} gene families, "
        f"{meta[group_col].nunique()} groups.",
        file=sys.stderr,
    )
    return matrix, meta


def safe_name(value: str) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return out.strip("_") or "group"


# ---------------------------------------------------------------------------
# Gene-family sharing
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GroupSets:
    sets: Dict[str, Set[str]]
    sizes: Dict[str, int]
    counts: pd.DataFrame
    prevalence: pd.DataFrame


def build_group_family_sets(
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    min_group_n: int,
    min_host_count: int,
    min_host_prevalence: float,
) -> GroupSets:
    if min_group_n < 1:
        raise ValueError("--min-group-n must be >= 1.")
    if min_host_count < 1:
        raise ValueError("--min-host-count must be >= 1.")
    if not 0 <= min_host_prevalence <= 1:
        raise ValueError("--min-host-prevalence must be between 0 and 1.")

    labels = meta[group_col].astype(str)
    group_sizes = labels.value_counts()
    groups = sorted(group_sizes[group_sizes >= min_group_n].index.astype(str))
    if len(groups) < 2:
        raise ValueError(
            f"Need at least two groups with n >= {min_group_n}; found {groups}."
        )

    count_cols: Dict[str, pd.Series] = {}
    prev_cols: Dict[str, pd.Series] = {}
    sets: Dict[str, Set[str]] = {}

    for group in groups:
        ids = labels.index[labels == group]
        counts = matrix.loc[ids].sum(axis=0).astype(int)
        prevalence = counts / float(len(ids))
        keep = (counts >= min_host_count) & (prevalence >= min_host_prevalence)

        count_cols[group] = counts
        prev_cols[group] = prevalence
        sets[group] = set(matrix.columns[keep])

    count_df = pd.DataFrame(count_cols).fillna(0).astype(int)
    prev_df = pd.DataFrame(prev_cols).fillna(0.0).astype(float)
    return GroupSets(
        sets=sets,
        sizes={g: int(group_sizes[g]) for g in groups},
        counts=count_df,
        prevalence=prev_df,
    )


def pairwise_set_metrics(
    sets: Mapping[str, Set[str]],
    sizes: Optional[Mapping[str, int]] = None,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for a, b in combinations(sorted(sets), 2):
        A, B = sets[a], sets[b]
        inter = A & B
        union = A | B
        minimum = min(len(A), len(B))

        jaccard = len(inter) / len(union) if union else np.nan
        dice = 2 * len(inter) / (len(A) + len(B)) if (len(A) + len(B)) else np.nan
        overlap = len(inter) / minimum if minimum else np.nan
        contain_a_b = len(inter) / len(A) if A else np.nan
        contain_b_a = len(inter) / len(B) if B else np.nan

        rows.append({
            "group_a": a,
            "group_b": b,
            "n_genomes_a": int(sizes[a]) if sizes else np.nan,
            "n_genomes_b": int(sizes[b]) if sizes else np.nan,
            "families_a": len(A),
            "families_b": len(B),
            "shared_families": len(inter),
            "union_families": len(union),
            "unique_a": len(A - B),
            "unique_b": len(B - A),
            "jaccard": jaccard,
            "dice": dice,
            "overlap_coefficient": overlap,
            "containment_a_in_b": contain_a_b,
            "containment_b_in_a": contain_b_a,
        })
    return pd.DataFrame(rows)


def classify_family_sharing(
    group_sets: GroupSets,
    family_metadata: pd.DataFrame,
) -> pd.DataFrame:
    groups = sorted(group_sets.sets)
    all_families = group_sets.counts.index.astype(str)
    out = pd.DataFrame(index=all_families)
    out.index.name = "gene_family"
    out["gene_family"] = out.index

    for group in groups:
        key = safe_name(group)
        out[f"{key}_count"] = group_sets.counts[group].reindex(out.index).fillna(0).astype(int)
        out[f"{key}_prevalence"] = (
            group_sets.prevalence[group].reindex(out.index).fillna(0.0).astype(float)
        )
        out[f"{key}_included"] = out.index.isin(group_sets.sets[group])

    included_cols = [f"{safe_name(g)}_included" for g in groups]
    included = out[included_cols].astype(bool)
    out["n_groups"] = included.sum(axis=1).astype(int)
    out["groups_present"] = included.apply(
        lambda row: ";".join(
            group for group, col in zip(groups, included_cols) if bool(row[col])
        ),
        axis=1,
    )

    n_groups_total = len(groups)
    out["sharing_class"] = np.select(
        [
            out["n_groups"].eq(0),
            out["n_groups"].eq(1),
            out["n_groups"].eq(2),
            out["n_groups"].eq(n_groups_total),
        ],
        [
            "below_threshold_all_groups",
            "group_specific",
            "shared_between_two_groups",
            "ubiquitous_all_groups",
        ],
        default="shared_between_multiple_groups",
    )

    if not family_metadata.empty:
        cols = [
            c for c in [
                "gene_name", "annotation", "n_genomes", "prevalence", "source_tool"
            ]
            if c in family_metadata.columns
        ]
        if cols:
            out = out.join(family_metadata[cols].reindex(out.index), how="left")

    front = [
        "gene_family", "gene_name", "annotation", "n_genomes", "prevalence",
        "source_tool", "n_groups", "groups_present", "sharing_class",
    ]
    ordered = [c for c in front if c in out.columns]
    ordered += [c for c in out.columns if c not in ordered]
    return out[ordered].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Rarefaction
# ---------------------------------------------------------------------------


def rarefied_pairwise_metrics(
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    groups: Sequence[str],
    min_host_count: int,
    min_host_prevalence: float,
    rarefy_n: Optional[int],
    reps: int,
    seed: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if reps < 1:
        raise ValueError("--rarefaction-reps must be >= 1.")

    labels = meta[group_col].astype(str)
    rng = np.random.default_rng(seed)
    replicate_rows: List[Dict[str, object]] = []

    for pair_index, (a, b) in enumerate(combinations(sorted(groups), 2)):
        ids_a = labels.index[labels == a].to_numpy()
        ids_b = labels.index[labels == b].to_numpy()
        n = min(len(ids_a), len(ids_b)) if rarefy_n is None else int(rarefy_n)
        if n < 1:
            continue
        if n > len(ids_a) or n > len(ids_b):
            raise ValueError(
                f"Requested rarefaction n={n} exceeds sample size for {a} "
                f"(n={len(ids_a)}) or {b} (n={len(ids_b)})."
            )

        for rep in range(reps):
            sub_a = rng.choice(ids_a, size=n, replace=False)
            sub_b = rng.choice(ids_b, size=n, replace=False)

            counts_a = matrix.loc[sub_a].sum(axis=0)
            counts_b = matrix.loc[sub_b].sum(axis=0)
            prev_a = counts_a / float(n)
            prev_b = counts_b / float(n)

            A = set(matrix.columns[
                (counts_a >= min_host_count) & (prev_a >= min_host_prevalence)
            ])
            B = set(matrix.columns[
                (counts_b >= min_host_count) & (prev_b >= min_host_prevalence)
            ])

            metrics = pairwise_set_metrics({a: A, b: B}, {a: n, b: n}).iloc[0].to_dict()
            metrics.update({
                "replicate": rep + 1,
                "rarefied_n_per_group": n,
                "seed": seed + pair_index,
            })
            replicate_rows.append(metrics)

    reps_df = pd.DataFrame(replicate_rows)
    if reps_df.empty:
        return reps_df, reps_df

    metric_cols = [
        "families_a", "families_b", "shared_families", "union_families",
        "unique_a", "unique_b", "jaccard", "dice", "overlap_coefficient",
        "containment_a_in_b", "containment_b_in_a",
    ]
    summary_rows: List[Dict[str, object]] = []
    for (a, b), block in reps_df.groupby(["group_a", "group_b"], sort=True):
        row: Dict[str, object] = {
            "group_a": a,
            "group_b": b,
            "rarefied_n_per_group": int(block["rarefied_n_per_group"].iloc[0]),
            "replicates": len(block),
        }
        for metric in metric_cols:
            vals = pd.to_numeric(block[metric], errors="coerce")
            row[f"{metric}_mean"] = vals.mean()
            row[f"{metric}_lo"] = vals.quantile(0.025)
            row[f"{metric}_hi"] = vals.quantile(0.975)
        summary_rows.append(row)

    return reps_df, pd.DataFrame(summary_rows)


# ---------------------------------------------------------------------------
# Lineage-collapsed sensitivity analysis
# ---------------------------------------------------------------------------


def collapse_by_structure(
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    structure_col: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if structure_col not in meta.columns:
        raise ValueError(f"Metadata lacks structure column '{structure_col}'.")

    valid = (
        meta[group_col].notna()
        & meta[group_col].astype(str).str.strip().ne("")
        & meta[structure_col].notna()
        & meta[structure_col].astype(str).str.strip().ne("")
    )
    m = matrix.loc[valid]
    md = meta.loc[valid, [group_col, structure_col]].copy()
    md["_unit"] = (
        md[group_col].astype(str)
        + "|||"
        + md[structure_col].astype(str)
    )

    collapsed = m.groupby(md["_unit"], sort=False).max().astype(np.uint8)
    unit_meta = md.drop_duplicates("_unit").set_index("_unit")
    unit_meta.index = unit_meta.index.astype(str)
    collapsed.index = collapsed.index.astype(str)
    return collapsed, unit_meta


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def metric_matrix(
    pairwise: pd.DataFrame,
    groups: Sequence[str],
    metric: str,
    diagonal: float = 1.0,
) -> pd.DataFrame:
    mat = pd.DataFrame(np.nan, index=groups, columns=groups, dtype=float)
    np.fill_diagonal(mat.values, diagonal)
    for row in pairwise.itertuples(index=False):
        value = getattr(row, metric)
        mat.loc[row.group_a, row.group_b] = value
        mat.loc[row.group_b, row.group_a] = value
    return mat


def plot_heatmap(
    matrix: pd.DataFrame,
    title: str,
    label: str,
    outbase: Path,
    fmt: str,
) -> None:
    fig, ax = plt.subplots(figsize=(max(6, 0.75 * len(matrix)), max(5, 0.65 * len(matrix))))
    image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto")
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index)
    ax.set_title(title)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iat[i, j]
            if np.isfinite(value):
                text = format(value, fmt)
                ax.text(j, i, text, ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(label)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_pairwise_bars(
    pairwise: pd.DataFrame,
    metric: str,
    title: str,
    ylabel: str,
    outbase: Path,
) -> None:
    if pairwise.empty:
        return
    d = pairwise.copy()
    d["pair"] = d["group_a"].astype(str) + " vs " + d["group_b"].astype(str)
    d = d.sort_values(metric, ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.38 * len(d))))
    ax.barh(d["pair"], d[metric])
    ax.set_xlabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_upset(
    family_table: pd.DataFrame,
    groups: Sequence[str],
    outbase: Path,
    top_n: int,
) -> bool:
    try:
        from upsetplot import UpSet, from_memberships
    except ImportError:
        print(
            "Skipping UpSet plot: optional dependency 'upsetplot' is not installed. "
            "Install with: conda install -c conda-forge upsetplot",
            file=sys.stderr,
        )
        return False

    memberships = [
        tuple(x.split(";"))
        for x in family_table.loc[
            family_table["groups_present"].astype(str).ne(""), "groups_present"
        ]
    ]
    if not memberships:
        return False

    data = from_memberships(memberships)
    fig = plt.figure(figsize=(11, 7))
    upset = UpSet(
        data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        max_subset_rank=top_n,
    )
    upset.plot(fig=fig)
    fig.suptitle("Gene-family sharing among groups", y=1.02)
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return True


def plot_network(
    pairwise: pd.DataFrame,
    group_sizes: Mapping[str, int],
    metric: str,
    minimum_edge: float,
    outbase: Path,
    seed: int,
) -> bool:
    try:
        import networkx as nx
    except ImportError:
        print(
            "Skipping network plot: optional dependency 'networkx' is not installed. "
            "Install with: conda install -c conda-forge networkx",
            file=sys.stderr,
        )
        return False

    graph = nx.Graph()
    for group, n in group_sizes.items():
        graph.add_node(group, size=n)

    for row in pairwise.itertuples(index=False):
        value = float(getattr(row, metric))
        if np.isfinite(value) and value >= minimum_edge:
            graph.add_edge(row.group_a, row.group_b, weight=value)

    if graph.number_of_nodes() == 0:
        return False

    pos = nx.spring_layout(graph, seed=seed, weight="weight")
    node_sizes = [
        350 + 1000 * math.sqrt(graph.nodes[node]["size"] / max(group_sizes.values()))
        for node in graph.nodes
    ]
    edge_widths = [
        0.5 + 7.0 * graph.edges[edge]["weight"] for edge in graph.edges
    ]

    fig, ax = plt.subplots(figsize=(8, 7))
    nx.draw_networkx_nodes(graph, pos, node_size=node_sizes, ax=ax)
    nx.draw_networkx_labels(graph, pos, font_size=10, ax=ax)
    if graph.number_of_edges():
        nx.draw_networkx_edges(graph, pos, width=edge_widths, alpha=0.65, ax=ax)
        edge_labels = {
            edge: f"{graph.edges[edge]['weight']:.2f}" for edge in graph.edges
        }
        nx.draw_networkx_edge_labels(
            graph, pos, edge_labels=edge_labels, font_size=8, ax=ax
        )
    ax.set_title(f"Group sharing network ({metric.replace('_', ' ')})")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# Exact sequence-allele sharing
# ---------------------------------------------------------------------------


FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna", ".ffn", ".aln", ".faa"}


def normalise_sequence(sequence: str, protein: bool) -> str:
    seq = str(sequence).upper().replace("-", "").replace(".", "").replace(" ", "")
    if protein:
        return re.sub(r"[^A-Z*]", "", seq)
    return re.sub(r"[^ACGTN]", "", seq)


def hash_allele(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:20]


def infer_sample_id(record_id: str, metadata_ids: Set[str]) -> Optional[str]:
    rid = str(record_id).split()[0]
    if rid in metadata_ids:
        return rid

    candidates = [
        rid.split("|")[0],
        rid.split(";")[0],
        rid.rsplit("_", 1)[0],
    ]
    for candidate in candidates:
        if candidate in metadata_ids:
            return candidate

    # Conservative prefix matching: only accept an unambiguous metadata prefix.
    matches = [sample for sample in metadata_ids if rid.startswith(sample)]
    if len(matches) == 1:
        return matches[0]
    return None


def read_exact_alleles_from_fasta_dir(
    seq_dir: Path,
    metadata_ids: Set[str],
    protein: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    files = sorted(
        p for p in seq_dir.rglob("*")
        if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS
    )
    if not files:
        raise FileNotFoundError(f"No FASTA-like files found under {seq_dir}")

    rows: List[Dict[str, object]] = []
    unmatched: List[Dict[str, str]] = []

    for path in files:
        family = path.stem
        for record in SeqIO.parse(str(path), "fasta"):
            sample = infer_sample_id(record.id, metadata_ids)
            if sample is None:
                unmatched.append({
                    "file": str(path),
                    "record_id": str(record.id),
                })
                continue
            sequence = normalise_sequence(str(record.seq), protein=protein)
            if not sequence:
                continue
            rows.append({
                "sample": sample,
                "gene_family": family,
                "allele": hash_allele(sequence),
                "sequence_length": len(sequence),
                "source_file": str(path),
            })

    return pd.DataFrame(rows), pd.DataFrame(unmatched)


def load_allele_table(
    path: Path,
    sample_col: str,
    family_col: str,
    allele_col: str,
) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        table = pd.read_csv(path, dtype=str)
    else:
        table = pd.read_csv(path, sep="\t", dtype=str)

    missing = [
        c for c in [sample_col, family_col, allele_col] if c not in table.columns
    ]
    if missing:
        raise ValueError(
            f"Allele table lacks columns {missing}. Found: {list(table.columns)}"
        )
    out = table[[sample_col, family_col, allele_col]].rename(columns={
        sample_col: "sample",
        family_col: "gene_family",
        allele_col: "allele",
    })
    return out.dropna().drop_duplicates()


def allele_sharing_tables(
    alleles: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    min_group_n: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    data = alleles.copy()
    data["sample"] = data["sample"].astype(str)
    joined = data.merge(
        meta[[group_col]],
        left_on="sample",
        right_index=True,
        how="inner",
    )
    joined = joined.dropna(subset=[group_col])
    joined[group_col] = joined[group_col].astype(str)

    group_sizes = meta[group_col].astype(str).value_counts()
    valid_groups = set(group_sizes[group_sizes >= min_group_n].index.astype(str))
    joined = joined[joined[group_col].isin(valid_groups)]

    family_group = (
        joined.groupby(["gene_family", group_col])["allele"]
        .agg(lambda x: set(map(str, x)))
    )

    rows: List[Dict[str, object]] = []
    for family in sorted(joined["gene_family"].astype(str).unique()):
        block = family_group.loc[family] if family in family_group.index.levels[0] else None
        if block is None:
            continue
        groups = sorted(block.index.astype(str))
        for a, b in combinations(groups, 2):
            A = set(block.loc[a])
            B = set(block.loc[b])
            inter = A & B
            union = A | B
            rows.append({
                "gene_family": family,
                "group_a": a,
                "group_b": b,
                "alleles_a": len(A),
                "alleles_b": len(B),
                "shared_alleles": len(inter),
                "allele_union": len(union),
                "allele_jaccard": len(inter) / len(union) if union else np.nan,
                "same_family_no_shared_allele": bool(A and B and not inter),
            })

    family_pair = pd.DataFrame(rows)
    if family_pair.empty:
        return family_pair, family_pair

    pair_summary = (
        family_pair.groupby(["group_a", "group_b"], as_index=False)
        .agg(
            shared_gene_families=("gene_family", "nunique"),
            families_with_shared_allele=("shared_alleles", lambda x: int((x > 0).sum())),
            families_without_shared_allele=(
                "same_family_no_shared_allele", "sum"
            ),
            total_shared_alleles=("shared_alleles", "sum"),
            mean_allele_jaccard=("allele_jaccard", "mean"),
            median_allele_jaccard=("allele_jaccard", "median"),
        )
    )
    pair_summary["proportion_families_sharing_allele"] = (
        pair_summary["families_with_shared_allele"]
        / pair_summary["shared_gene_families"].replace(0, np.nan)
    )
    return family_pair, pair_summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Quantify observed and rarefied gene-family sharing, and optional "
            "exact sequence-allele sharing, among metadata groups."
        )
    )
    parser.add_argument(
        "--pangenome-tool", required=True, choices=["panaroo", "pirate"]
    )
    parser.add_argument("--pangenome-out", required=True)
    parser.add_argument("--panaroo-from-csv", action="store_true")
    parser.add_argument("--meta", required=True)
    parser.add_argument("--meta-id-col", default="sample")
    parser.add_argument("--group-col", default="host")
    parser.add_argument("--structure-col")
    parser.add_argument("--outdir", required=True)

    parser.add_argument("--min-group-n", type=int, default=10)
    parser.add_argument("--min-host-count", type=int, default=2)
    parser.add_argument("--min-host-prevalence", type=float, default=0.01)

    parser.add_argument("--rarefy", action="store_true")
    parser.add_argument(
        "--rarefy-n",
        type=int,
        help="Genomes sampled per group; default is the smaller group in each pair.",
    )
    parser.add_argument("--rarefaction-reps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=1)

    parser.add_argument("--heatmaps", action="store_true")
    parser.add_argument("--pairwise-bars", action="store_true")
    parser.add_argument("--upset", action="store_true")
    parser.add_argument("--upset-top-n", type=int, default=30)
    parser.add_argument("--network", action="store_true")
    parser.add_argument(
        "--network-metric",
        default="jaccard",
        choices=["jaccard", "dice", "overlap_coefficient"],
    )
    parser.add_argument("--network-min-edge", type=float, default=0.0)

    allele = parser.add_argument_group("optional exact allele sharing")
    allele.add_argument("--allele-table")
    allele.add_argument("--allele-sample-col", default="sample")
    allele.add_argument("--allele-family-col", default="gene_family")
    allele.add_argument("--allele-col", default="allele")
    allele.add_argument(
        "--allele-seq-dir",
        help=(
            "Directory containing family-specific FASTA/alignment files. "
            "Alleles are exact hashes of ungapped sequences."
        ),
    )
    allele.add_argument(
        "--allele-sequence-type",
        choices=["nucleotide", "protein"],
        default="nucleotide",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    set_plot_style()
    outdir = ensure_outdir(args.outdir)
    tables_dir = ensure_outdir(outdir / "tables")
    plots_dir = ensure_outdir(outdir / "plots")

    pg: PangenomeData = load_pangenome(
        tool=args.pangenome_tool,
        outdir=args.pangenome_out,
        prefer_panaroo_rtab=not args.panaroo_from_csv,
    )
    meta = load_metadata(args.meta, args.meta_id_col)
    matrix, meta = align_matrix_metadata(
        pg.presence_absence,
        meta,
        args.group_col,
    )

    group_sets = build_group_family_sets(
        matrix=matrix,
        meta=meta,
        group_col=args.group_col,
        min_group_n=args.min_group_n,
        min_host_count=args.min_host_count,
        min_host_prevalence=args.min_host_prevalence,
    )
    groups = sorted(group_sets.sets)

    observed = pairwise_set_metrics(group_sets.sets, group_sets.sizes)
    observed.to_csv(
        tables_dir / "pairwise_gene_family_sharing_observed.tsv",
        sep="\t",
        index=False,
    )

    family_table = classify_family_sharing(group_sets, pg.family_metadata)
    family_table.to_csv(
        tables_dir / "gene_family_sharing_by_group.tsv",
        sep="\t",
        index=False,
    )

    # Matrices are useful even when plots are not requested.
    for metric in ["shared_families", "jaccard", "dice", "overlap_coefficient"]:
        diagonal = (
            float(group_sets.sizes[g]) if metric == "shared_families" else 1.0
            for g in groups
        )
        if metric == "shared_families":
            mat = pd.DataFrame(
                0, index=groups, columns=groups, dtype=float
            )
            for g in groups:
                mat.loc[g, g] = len(group_sets.sets[g])
            for row in observed.itertuples(index=False):
                mat.loc[row.group_a, row.group_b] = row.shared_families
                mat.loc[row.group_b, row.group_a] = row.shared_families
        else:
            mat = metric_matrix(observed, groups, metric)
        mat.to_csv(tables_dir / f"matrix_{metric}_observed.tsv", sep="\t")

    rarefied_reps = pd.DataFrame()
    rarefied_summary = pd.DataFrame()
    if args.rarefy:
        rarefied_reps, rarefied_summary = rarefied_pairwise_metrics(
            matrix=matrix,
            meta=meta,
            group_col=args.group_col,
            groups=groups,
            min_host_count=args.min_host_count,
            min_host_prevalence=args.min_host_prevalence,
            rarefy_n=args.rarefy_n,
            reps=args.rarefaction_reps,
            seed=args.seed,
        )
        rarefied_reps.to_csv(
            tables_dir / "pairwise_gene_family_sharing_rarefied_replicates.tsv",
            sep="\t",
            index=False,
        )
        rarefied_summary.to_csv(
            tables_dir / "pairwise_gene_family_sharing_rarefied_summary.tsv",
            sep="\t",
            index=False,
        )

        for metric in ["shared_families", "jaccard", "dice", "overlap_coefficient"]:
            col = f"{metric}_mean"
            mat = metric_matrix(
                rarefied_summary.rename(columns={col: metric}),
                groups,
                metric,
            )
            mat.to_csv(
                tables_dir / f"matrix_{metric}_rarefied.tsv",
                sep="\t",
            )

    # Optional lineage-collapsed analysis.
    if args.structure_col:
        collapsed, unit_meta = collapse_by_structure(
            matrix, meta, args.group_col, args.structure_col
        )
        lineage_sets = build_group_family_sets(
            matrix=collapsed,
            meta=unit_meta,
            group_col=args.group_col,
            min_group_n=1,
            min_host_count=1,
            min_host_prevalence=args.min_host_prevalence,
        )
        lineage_pairwise = pairwise_set_metrics(
            lineage_sets.sets, lineage_sets.sizes
        )
        lineage_pairwise.to_csv(
            tables_dir
            / f"pairwise_gene_family_sharing_collapsed_by_{safe_name(args.structure_col)}.tsv",
            sep="\t",
            index=False,
        )

    # Plots.
    if args.heatmaps:
        shared_mat = pd.read_csv(
            tables_dir / "matrix_shared_families_observed.tsv",
            sep="\t",
            index_col=0,
        )
        jac_mat = pd.read_csv(
            tables_dir / "matrix_jaccard_observed.tsv",
            sep="\t",
            index_col=0,
        )
        plot_heatmap(
            shared_mat,
            "Observed shared gene families",
            "Shared gene families",
            plots_dir / "heatmap_shared_families_observed",
            ".0f",
        )
        plot_heatmap(
            jac_mat,
            "Observed gene-family Jaccard similarity",
            "Jaccard similarity",
            plots_dir / "heatmap_jaccard_observed",
            ".2f",
        )

        if args.rarefy and not rarefied_summary.empty:
            rare_jac = pd.read_csv(
                tables_dir / "matrix_jaccard_rarefied.tsv",
                sep="\t",
                index_col=0,
            )
            rare_shared = pd.read_csv(
                tables_dir / "matrix_shared_families_rarefied.tsv",
                sep="\t",
                index_col=0,
            )
            plot_heatmap(
                rare_jac,
                "Rarefied gene-family Jaccard similarity",
                "Mean Jaccard similarity",
                plots_dir / "heatmap_jaccard_rarefied",
                ".2f",
            )
            plot_heatmap(
                rare_shared,
                "Rarefied shared gene families",
                "Mean shared families",
                plots_dir / "heatmap_shared_families_rarefied",
                ".0f",
            )

    if args.pairwise_bars:
        plot_pairwise_bars(
            observed,
            metric="jaccard",
            title="Observed pairwise gene-family similarity",
            ylabel="Jaccard similarity",
            outbase=plots_dir / "pairwise_jaccard_observed",
        )
        if args.rarefy and not rarefied_summary.empty:
            plot_pairwise_bars(
                rarefied_summary,
                metric="jaccard_mean",
                title="Rarefied pairwise gene-family similarity",
                ylabel="Mean Jaccard similarity",
                outbase=plots_dir / "pairwise_jaccard_rarefied",
            )

    if args.upset:
        plot_upset(
            family_table=family_table,
            groups=groups,
            outbase=plots_dir / "upset_gene_family_sharing",
            top_n=args.upset_top_n,
        )

    network_source = observed
    network_metric = args.network_metric
    if args.rarefy and not rarefied_summary.empty:
        network_source = rarefied_summary.rename(
            columns={f"{network_metric}_mean": network_metric}
        )
    if args.network:
        plot_network(
            pairwise=network_source,
            group_sizes=group_sets.sizes,
            metric=network_metric,
            minimum_edge=args.network_min_edge,
            outbase=plots_dir / f"network_{network_metric}",
            seed=args.seed,
        )

    # Optional allele sharing.
    allele_table: Optional[pd.DataFrame] = None
    unmatched = pd.DataFrame()
    if args.allele_table:
        allele_table = load_allele_table(
            Path(args.allele_table),
            sample_col=args.allele_sample_col,
            family_col=args.allele_family_col,
            allele_col=args.allele_col,
        )
    elif args.allele_seq_dir:
        allele_table, unmatched = read_exact_alleles_from_fasta_dir(
            Path(args.allele_seq_dir),
            metadata_ids=set(meta.index.astype(str)),
            protein=args.allele_sequence_type == "protein",
        )
        allele_table.to_csv(
            tables_dir / "exact_sequence_alleles.tsv",
            sep="\t",
            index=False,
        )
        if not unmatched.empty:
            unmatched.to_csv(
                tables_dir / "unmatched_allele_sequence_headers.tsv",
                sep="\t",
                index=False,
            )

    if allele_table is not None:
        family_pair, allele_pair = allele_sharing_tables(
            allele_table,
            meta,
            args.group_col,
            args.min_group_n,
        )
        family_pair.to_csv(
            tables_dir / "allele_sharing_by_family_and_group_pair.tsv",
            sep="\t",
            index=False,
        )
        allele_pair.to_csv(
            tables_dir / "pairwise_allele_sharing_summary.tsv",
            sep="\t",
            index=False,
        )
        if not allele_pair.empty and args.heatmaps:
            allele_mat = metric_matrix(
                allele_pair,
                groups,
                "proportion_families_sharing_allele",
            )
            allele_mat.to_csv(
                tables_dir / "matrix_proportion_families_sharing_allele.tsv",
                sep="\t",
            )
            plot_heatmap(
                allele_mat,
                "Proportion of shared families containing a shared exact allele",
                "Proportion",
                plots_dir / "heatmap_shared_exact_alleles",
                ".2f",
            )

    manifest = {
        "pangenome_tool": args.pangenome_tool,
        "pangenome_out": str(Path(args.pangenome_out).resolve()),
        "metadata": str(Path(args.meta).resolve()),
        "group_col": args.group_col,
        "structure_col": args.structure_col,
        "n_samples": int(matrix.shape[0]),
        "n_gene_families": int(matrix.shape[1]),
        "groups": group_sets.sizes,
        "min_group_n": args.min_group_n,
        "min_host_count": args.min_host_count,
        "min_host_prevalence": args.min_host_prevalence,
        "rarefied": bool(args.rarefy),
        "rarefaction_reps": args.rarefaction_reps if args.rarefy else 0,
        "rarefy_n": args.rarefy_n,
        "allele_analysis": bool(allele_table is not None),
        "allele_definition": (
            f"exact_{args.allele_sequence_type}_sequence"
            if args.allele_seq_dir
            else "provided_allele_table"
            if args.allele_table
            else None
        ),
    }
    with (outdir / "analysis_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2)

    print(f"Done. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
