#!/usr/bin/env python3
"""
pirate_plotting.py (version 2.5)

Plotting helpers for PIRATE outputs (Campylobacter-scale datasets).
Generates publication-ready summaries without requiring core alignments.

Expected in --pirate-out:
  - binary_presence_absence.fasta
Optional (if present / requested):
  - PIRATE.gene_families.tsv
  - PIRATE.pangenome_summary.txt

Dependencies:
  - pandas, numpy, matplotlib
  - scikit-learn
  - biopython
Optional:
  - umap-learn

Examples:
  python pirate_plotting.py --pirate-out PIRATE_out --outdir plots --all
  python pirate_plotting.py --pirate-out PIRATE_out --outdir plots --pca --meta meta.txt --meta-col host
  python pirate_plotting.py --pirate-out PIRATE_out --outdir plots --pca --meta meta.txt --meta-col host \
    --group-colors "human_gastroenteritis=#9C2C2C,human_asymptomatic=#003B5C"
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors


# -----------------------------
# Palettes
# -----------------------------
WES_PALETTES: Dict[str, List[str]] = {
    "life_aquatic_blues": [
        "#003B5C",  # deep navy
        "#005F73",  # teal
        "#0A9396",  # sea green
        "#94D2BD",  # pale aqua
        "#E9D8A6",  # sand
        "#1D3557",  # dark blue
        "#457B9D",  # steel blue
        "#A8DADC",  # light blue
    ],
    "grand_budapest_soft": [
        "#5B1A1A", "#9C2C2C", "#D04E4E", "#F2B5B5", "#F7EDE2", "#3E4E50", "#BFD7EA"
    ],
}

DEFAULT_UNKNOWN_COLOR = "#BDBDBD"  # also used for REFERENCE unless overridden


def get_palette(name: str) -> List[str]:
    if name not in WES_PALETTES:
        raise ValueError(f"Unknown palette '{name}'. Available: {', '.join(WES_PALETTES.keys())}")
    return WES_PALETTES[name]


# -----------------------------
# IO helpers
# -----------------------------
def ensure_outdir(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)


def warn(msg: str) -> None:
    print(f"[WARN] {msg}")


def read_binary_presence_absence_fasta(path: str) -> pd.DataFrame:
    """
    Read PIRATE binary presence/absence FASTA into a dataframe:
      rows = genomes, cols = gene-family positions, values = 0/1 integers

    Supports PIRATE encodings:
      - A/C (A=0 absent, C=1 present)
      - 0/1
      - '-' treated as 0
    """
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"No records found in {path}")

    ids: List[str] = []
    arrs: List[np.ndarray] = []
    for r in records:
        seq = str(r.seq).strip().upper()
        # Map: C/1 -> 1, everything else -> 0 (A/0/-/N)
        arr = np.fromiter((1 if c in ("C", "1") else 0 for c in seq), dtype=np.int8)
        ids.append(r.id)
        arrs.append(arr)

    mat = np.vstack(arrs)
    df = pd.DataFrame(mat, index=ids)
    df.index.name = "sample"
    return df


def read_gene_families_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def filter_variable_columns(df01: pd.DataFrame, min_var: float = 0.0) -> pd.DataFrame:
    """Drop gene families with zero variance (all 0 or all 1)."""
    v = df01.var(axis=0)
    keep = v > min_var
    return df01.loc[:, keep]


def load_metadata(meta_path: str, id_col: str = "sample") -> pd.DataFrame:
    """
    Loads a metadata TSV/CSV.
    - .csv -> comma-separated
    - otherwise -> tab-separated (meta.txt typical)
    Values are loaded as strings to avoid mixed-type warnings.
    """
    if not os.path.exists(meta_path):
        raise FileNotFoundError(meta_path)

    if meta_path.endswith(".csv"):
        df = pd.read_csv(meta_path, dtype=str)
    else:
        df = pd.read_csv(meta_path, sep="\t", dtype=str)

    if id_col not in df.columns:
        raise ValueError(f"Metadata file must contain column '{id_col}'. Found: {list(df.columns)}")

    df = df.set_index(id_col)
    return df


# -----------------------------
# Styling + legend utilities
# -----------------------------
def apply_matplotlib_style() -> None:
    plt.rcParams.update({
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "font.size": 12,
        "axes.titlesize": 18,
        "axes.labelsize": 14,
        "legend.fontsize": 12,
        "legend.framealpha": 0.95,
        "savefig.bbox": "tight",
    })


def parse_group_colors(s: Optional[str]) -> Dict[str, str]:
    """
    Parse "a=#RRGGBB,b=#RRGGBB" into dict.
    Accepts commas and optional whitespace.
    """
    if not s:
        return {}
    out: Dict[str, str] = {}
    parts = [p.strip() for p in s.split(",") if p.strip()]
    for p in parts:
        if "=" not in p:
            continue
        k, v = p.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def build_color_map(
    labels: pd.Series,
    palette: List[str],
    group_colors: Dict[str, str],
    unknown_label: str = "Unknown",
) -> Tuple[pd.Series, Dict[str, str], List[str]]:
    lab = labels.copy()
    lab = lab.astype(str)
    lab = lab.replace({"nan": np.nan, "None": np.nan})
    lab = lab.where(lab.notna(), other=unknown_label)

    levels = sorted(pd.unique(lab))
    def _key(x: str) -> Tuple[int, str]:
        if x == "REFERENCE":
            return (0, x)
        if x == unknown_label:
            return (1, x)
        return (2, x)

    levels = sorted(levels, key=_key)

    color_map: Dict[str, str] = {}
    pal_i = 0
    for lvl in levels:
        if lvl in group_colors:
            color_map[lvl] = group_colors[lvl]
        elif lvl == "REFERENCE":
            color_map[lvl] = group_colors.get("REFERENCE", DEFAULT_UNKNOWN_COLOR)
        elif lvl == unknown_label:
            color_map[lvl] = group_colors.get(unknown_label, DEFAULT_UNKNOWN_COLOR)
        else:
            color_map[lvl] = palette[pal_i % len(palette)]
            pal_i += 1

    return lab, color_map, levels


def export_legend_only(
    levels: List[str],
    color_map: Dict[str, str],
    title: str,
    out_png: str,
    marker: str = "o",
    markersize: float = 10.0,
) -> None:
    fig = plt.figure(figsize=(6, max(2.5, 0.5 * len(levels))))
    ax = fig.add_subplot(111)
    ax.axis("off")

    handles = []
    for lvl in levels:
        handles.append(
            plt.Line2D(
                [0], [0],
                marker=marker, linestyle="",
                label=str(lvl),
                markerfacecolor=color_map[lvl],
                markeredgecolor="white",
                markeredgewidth=0.9,
                markersize=markersize,
            )
        )
    ax.legend(handles=handles, title=title, loc="center left", frameon=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# -----------------------------
# Core plots
# -----------------------------
def save_coords_csv(coords: np.ndarray, samples: List[str], out_csv: str, extra: Optional[pd.DataFrame] = None) -> None:
    df = pd.DataFrame(coords, index=samples, columns=["dim1", "dim2"])
    if extra is not None:
        df = pd.concat([df, extra.reindex(df.index)], axis=1)
    df.to_csv(out_csv)


def add_knn_edges(ax: plt.Axes, coords: np.ndarray, k: int = 8, alpha: float = 0.18, lw: float = 1.0) -> None:
    if k is None or k <= 0:
        return
    if coords.shape[0] < 3:
        return
    nn = NearestNeighbors(n_neighbors=min(k + 1, coords.shape[0]), metric="euclidean")
    nn.fit(coords)
    inds = nn.kneighbors(coords, return_distance=False)
    for i in range(coords.shape[0]):
        for j in inds[i, 1:]:
            ax.plot(
                [coords[i, 0], coords[j, 0]],
                [coords[i, 1], coords[j, 1]],
                color="#CFCFCF",
                alpha=alpha,
                linewidth=lw,
                zorder=1,
            )


def plot_scatter_pretty(
    coords: np.ndarray,
    labels: Optional[pd.Series],
    palette: List[str],
    title: str,
    out_png: str,
    xlabel: str = "Dim 1",
    ylabel: str = "Dim 2",
    group_colors: Optional[Dict[str, str]] = None,
    legend_outside: bool = False,
    knn_edges: int = 0,
    point_size: float = 38.0,
    point_alpha: float = 0.95,
    legend_title: Optional[str] = None,
    legend_out_png: Optional[str] = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 6.8))

    if knn_edges and knn_edges > 0:
        add_knn_edges(ax, coords, k=knn_edges)

    if labels is None:
        ax.scatter(coords[:, 0], coords[:, 1], s=point_size, c=palette[0],
                   edgecolors="white", linewidths=0.6, alpha=point_alpha, zorder=2)
    else:
        group_colors = group_colors or {}
        lab, cmap, levels = build_color_map(labels, palette, group_colors)

        for lvl in levels:
            idx = (lab == lvl).values
            ax.scatter(
                coords[idx, 0], coords[idx, 1],
                s=point_size,
                c=cmap[lvl],
                edgecolors="white",
                linewidths=0.6,
                alpha=point_alpha,
                label=str(lvl),
                zorder=3,
            )

        if legend_outside:
            ax.legend(
                title=legend_title or (labels.name if labels is not None else "group"),
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=True,
            )
        else:
            ax.legend(
                title=legend_title or (labels.name if labels is not None else "group"),
                loc="best",
                frameon=True,
            )

        if legend_out_png:
            export_legend_only(levels, cmap, title=legend_title or (labels.name or "group"), out_png=legend_out_png)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_gene_frequency_hist(gene_families: pd.DataFrame, out_png: str) -> None:
    col_candidates = [
        "number_genomes",
        "number_genome",
        "No. isolates", "No. isolates ", "No_isolates",
        "No. genomes", "No. Genomes",
    ]
    col = None
    for c in col_candidates:
        if c in gene_families.columns:
            col = c
            break
    if col is None:
        raise ValueError(f"Could not find isolates/genomes count column. Columns: {list(gene_families.columns)}")

    if "threshold" in gene_families.columns:
        gf95 = gene_families[gene_families["threshold"].astype(str) == "95"]
        if len(gf95) > 0:
            gene_families = gf95

    freq = pd.to_numeric(gene_families[col], errors="coerce").dropna().astype(int)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.hist(freq, bins=60)
    ax.set_xlabel("Number of genomes containing gene family")
    ax.set_ylabel("Number of gene families")
    ax.set_title("Gene family frequency distribution")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def compute_pangenome_bins_from_gene_families(
    gene_families: pd.DataFrame,
    n_genomes: int,
    core: float = 0.99,
    soft_core: float = 0.95,
    shell: float = 0.15,
) -> Dict[str, int]:
    col_candidates = [
        "number_genomes",
        "number_genome",
        "No. isolates", "No. isolates ", "No_isolates",
        "No. genomes", "No. Genomes",
    ]
    col = None
    for c in col_candidates:
        if c in gene_families.columns:
            col = c
            break
    if col is None:
        raise ValueError(f"Could not find isolates/genomes count column. Columns: {list(gene_families.columns)}")

    if "threshold" in gene_families.columns:
        gf95 = gene_families[gene_families["threshold"].astype(str) == "95"]
        if len(gf95) > 0:
            gene_families = gf95

    freq = pd.to_numeric(gene_families[col], errors="coerce").dropna().astype(int)

    core_thr = int(np.ceil(core * n_genomes))
    soft_thr = int(np.ceil(soft_core * n_genomes))
    shell_thr = int(np.ceil(shell * n_genomes))

    return {
        "Core (≥99%)": int((freq >= core_thr).sum()),
        "Soft-core (95–99%)": int(((freq >= soft_thr) & (freq < core_thr)).sum()),
        "Shell (15–95%)": int(((freq >= shell_thr) & (freq < soft_thr)).sum()),
        "Cloud (<15%)": int((freq < shell_thr).sum()),
    }


def plot_pangenome_bars(bins: Dict[str, int], n_genomes: int, out_png: str, palette: List[str]) -> None:
    labels_bar = list(bins.keys())
    values_bar = list(bins.values())
    colors = [palette[i % len(palette)] for i in range(len(values_bar))]

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.bar(labels_bar, values_bar, color=colors)
    ax.set_ylabel("Number of gene families")
    ax.set_title(f"Pangenome composition (n={n_genomes} genomes)")
    ax.set_xticklabels(labels_bar, rotation=20, ha="right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# -----------------------------
# Rarefaction
# -----------------------------
def compute_rarefaction_curves(
    df01: pd.DataFrame,
    steps: int = 60,
    reps: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n = df01.shape[0]

    if steps >= n:
        ks = np.arange(1, n + 1)
    else:
        ks = np.unique(np.round(np.linspace(1, n, steps)).astype(int))

    X = df01.values.astype(np.uint8)

    pan_mat = np.zeros((reps, len(ks)), dtype=np.int32)
    core_mat = np.zeros((reps, len(ks)), dtype=np.int32)

    for r in range(reps):
        order = rng.permutation(n)
        Xp = X[order, :]
        csum = np.cumsum(Xp, axis=0)
        for j, k in enumerate(ks):
            pan_mat[r, j] = int((csum[k - 1, :] > 0).sum())
            core_mat[r, j] = int((csum[k - 1, :] == k).sum())

    def summarise(mat: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        mean = mat.mean(axis=0)
        lo = np.percentile(mat, 2.5, axis=0)
        hi = np.percentile(mat, 97.5, axis=0)
        return mean, lo, hi

    pan_mean, pan_lo, pan_hi = summarise(pan_mat)
    core_mean, core_lo, core_hi = summarise(core_mat)

    return pd.DataFrame({
        "n_genomes": ks,
        "pan_mean": pan_mean, "pan_lo": pan_lo, "pan_hi": pan_hi,
        "core_mean": core_mean, "core_lo": core_lo, "core_hi": core_hi,
    })


def compute_rarefaction_curves_by_group(
    df01: pd.DataFrame,
    groups: pd.Series,
    steps: int = 60,
    reps: int = 20,
    seed: int = 0,
    min_group_size: int = 30,
    max_groups: Optional[int] = None,
) -> Dict[str, pd.DataFrame]:
    g = groups.reindex(df01.index)
    keep = g.notna()
    df01 = df01.loc[keep]
    g = g.loc[keep].astype(str)

    sizes = g.value_counts()
    sizes = sizes[sizes >= min_group_size]
    if max_groups is not None:
        sizes = sizes.iloc[:max_groups]

    out: Dict[str, pd.DataFrame] = {}
    for grp in sizes.index:
        sub = df01.loc[g == grp]
        out[grp] = compute_rarefaction_curves(sub, steps=steps, reps=reps, seed=seed)
    return out


def plot_rarefaction(
    curves: pd.DataFrame,
    out_png: str,
    palette: List[str],
    title: str = "Gene family accumulation (PIRATE)",
    legend_outside: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 6.2))
    x = curves["n_genomes"].values
    pan_c = palette[1]
    core_c = palette[0]

    ax.plot(x, curves["pan_mean"].values, linewidth=3, color=pan_c, label="Pangenome (total)")
    ax.fill_between(x, curves["pan_lo"].values, curves["pan_hi"].values, alpha=0.20, color=pan_c)

    ax.plot(x, curves["core_mean"].values, linewidth=3, color=core_c, label="Core (present in all)")
    ax.fill_between(x, curves["core_lo"].values, curves["core_hi"].values, alpha=0.20, color=core_c)

    ax.set_xlabel("Number of genomes")
    ax.set_ylabel("Number of gene families")
    ax.set_title(title)

    if legend_outside:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
    else:
        ax.legend(frameon=True, loc="best")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_rarefaction_by_group(
    curves_by_group: Dict[str, pd.DataFrame],
    out_png: str,
    palette: List[str],
    group_colors: Optional[Dict[str, str]] = None,
    title: str = "Gene family accumulation by group",
    include_core: bool = False,
    legend_outside: bool = False,
) -> None:
    group_colors = group_colors or {}
    fig, ax = plt.subplots(figsize=(8.8, 6.4))

    ordered = list(curves_by_group.keys())

    for i, grp in enumerate(ordered):
        dfc = curves_by_group[grp]
        c = group_colors.get(grp, palette[i % len(palette)])
        x = dfc["n_genomes"].values

        ax.plot(x, dfc["pan_mean"].values, linewidth=3.5, color=c, label=grp)
        ax.fill_between(x, dfc["pan_lo"].values, dfc["pan_hi"].values, alpha=0.18, color=c)

        if include_core:
            ax.plot(x, dfc["core_mean"].values, linewidth=2.5, color=c, linestyle="--", alpha=0.9)

    ax.set_xlabel("Number of isolates")
    ax.set_ylabel("Number of gene families")
    ax.set_title(title)

    if legend_outside:
        ax.legend(title="Group", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
    else:
        ax.legend(title="Group", frameon=True, loc="best")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# -----------------------------
# Embeddings + diagnostics
# -----------------------------
def run_pca_accessory(df01: pd.DataFrame, n_components: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    pca = PCA(n_components=n_components, random_state=0)
    coords = pca.fit_transform(df01.values)
    return coords, pca.explained_variance_ratio_


def run_mds_jaccard(df01: pd.DataFrame) -> np.ndarray:
    X = df01.values.astype(bool)
    dist = pairwise_distances(X, metric="jaccard")
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0, n_init=4, max_iter=300)
    return mds.fit_transform(dist)


def run_umap(df01: pd.DataFrame, n_neighbors: int = 15, min_dist: float = 0.1) -> np.ndarray:
    try:
        import umap  # type: ignore
    except ImportError as e:
        raise ImportError("UMAP requested but umap-learn is not installed. Install with: conda install -c conda-forge umap-learn") from e
    reducer = umap.UMAP(n_components=2, n_neighbors=n_neighbors, min_dist=min_dist, metric="jaccard", random_state=0)
    return reducer.fit_transform(df01.values.astype(bool))


def gene_counts_per_genome(df01: pd.DataFrame) -> pd.Series:
    counts = df01.sum(axis=1)
    counts.name = "gene_count"
    return counts


def plot_gene_count_distribution(gene_counts: pd.Series, out_png: str) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.hist(gene_counts, bins=50)
    ax.set_xlabel("Number of gene families per genome")
    ax.set_ylabel("Number of genomes")
    ax.set_title("Gene content per genome")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_pca_vs_gene_count(coords: np.ndarray, gene_counts: pd.Series, out_png: str, var: Tuple[float, float]) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=gene_counts.values, cmap="viridis", s=14)
    ax.set_xlabel(f"PC1 ({var[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]*100:.1f}%)")
    ax.set_title("Accessory PCA coloured by gene count")
    fig.colorbar(sc, ax=ax, label="Genes per genome")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_top_variance_heatmap(df01: pd.DataFrame, top_n: int, out_png: str) -> None:
    variances = df01.var(axis=0)
    top_cols = variances.sort_values(ascending=False).head(top_n).index
    mat = df01[top_cols]

    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(mat.T, aspect="auto", interpolation="nearest", cmap="viridis")
    fig.colorbar(im, ax=ax, label="Presence (0/1)")
    ax.set_xlabel("Genomes")
    ax.set_ylabel(f"Top {top_n} variable gene families")
    ax.set_title(f"Top {top_n} variable gene families (presence/absence)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# -----------------------------
# CLI runner with skip/summary
# -----------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="Plot PIRATE summaries + accessory structure (Campy-friendly).")
    ap.add_argument("--pirate-out", required=True, help="Path to PIRATE_out directory")
    ap.add_argument("--outdir", default="pirate_plots", help="Output directory for plots")
    ap.add_argument("--palette", default="life_aquatic_blues", help=f"Palette: {', '.join(WES_PALETTES.keys())}")

    ap.add_argument("--all", action="store_true",
                    help="Run the standard suite of outputs (matches historical default outputs).")

    ap.add_argument("--pca", action="store_true")
    ap.add_argument("--mds", action="store_true")
    ap.add_argument("--umap", action="store_true")

    ap.add_argument("--gene-freq", action="store_true")
    ap.add_argument("--pangenome-bars", action="store_true")

    ap.add_argument("--gene-counts", action="store_true")
    ap.add_argument("--pca-gene-count", action="store_true")
    ap.add_argument("--heatmap-top", type=int, default=None)

    ap.add_argument("--rarefaction", action="store_true")
    ap.add_argument("--rare-steps", type=int, default=60)
    ap.add_argument("--rare-reps", type=int, default=20)
    ap.add_argument("--rare-seed", type=int, default=0)

    ap.add_argument("--rarefaction-by", default=None,
                    help="Metadata column to stratify rarefaction curves (e.g. host). Requires --meta.")
    ap.add_argument("--min-group-size", type=int, default=30)
    ap.add_argument("--max-groups", type=int, default=None)
    ap.add_argument("--rarefaction-core", action="store_true")

    ap.add_argument("--meta", default=None)
    ap.add_argument("--meta-id-col", default="sample")
    ap.add_argument("--meta-col", default=None)

    ap.add_argument("--umap-n-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.1)

    ap.add_argument("--group-colors", default=None)
    ap.add_argument("--legend-outside", action="store_true")
    ap.add_argument("--legend-out", default=None)
    ap.add_argument("--knn-edges", type=int, default=0)

    args = ap.parse_args()

    apply_matplotlib_style()

    pirate_out = args.pirate_out.rstrip("/")
    outdir = args.outdir
    ensure_outdir(outdir)

    palette = get_palette(args.palette)
    group_colors = parse_group_colors(args.group_colors)

    if args.all:
        args.pangenome_bars = True
        args.gene_freq = True
        args.pca = True
        args.mds = True
        args.gene_counts = True
        args.pca_gene_count = True
        args.heatmap_top = 200
        args.rarefaction = True

    ran: List[str] = []
    skipped: List[str] = []

    meta_df: Optional[pd.DataFrame] = None
    labels: Optional[pd.Series] = None
    extra: Optional[pd.DataFrame] = None

    if args.meta:
        try:
            meta_df = load_metadata(args.meta, id_col=args.meta_id_col)
            extra = meta_df
            if args.meta_col:
                if args.meta_col not in meta_df.columns:
                    warn(f"Metadata column '{args.meta_col}' not found; available: {list(meta_df.columns)}")
                else:
                    labels = meta_df[args.meta_col]
                    labels.name = args.meta_col
        except Exception as e:
            warn(f"Failed to load metadata '{args.meta}': {e}")
            meta_df = None
            labels = None
            extra = None

    gf: Optional[pd.DataFrame] = None
    gf_path = os.path.join(pirate_out, "PIRATE.gene_families.tsv")
    if args.pangenome_bars or args.gene_freq:
        if os.path.exists(gf_path):
            try:
                gf = read_gene_families_tsv(gf_path)
            except Exception as e:
                warn(f"Could not read {gf_path}: {e}")
                gf = None
        else:
            warn(f"Missing {gf_path}; pangenome bars / gene frequency will be skipped.")
            gf = None

    need_df01 = (
        args.pca or args.mds or args.umap or args.gene_counts or args.pca_gene_count
        or (args.heatmap_top is not None) or args.rarefaction
    )

    df01: Optional[pd.DataFrame] = None
    bin_path = os.path.join(pirate_out, "binary_presence_absence.fasta")
    if need_df01:
        if not os.path.exists(bin_path):
            warn(f"Missing {bin_path}; PCA/MDS/UMAP/rarefaction/diagnostics will be skipped.")
            df01 = None
        else:
            try:
                df01 = read_binary_presence_absence_fasta(bin_path)
                df01 = filter_variable_columns(df01)
                if df01.shape[1] < 2:
                    warn("After filtering, <2 variable gene families remain; skipping embeddings.")
                    df01 = None
            except Exception as e:
                warn(f"Failed to load binary presence/absence FASTA: {e}")
                df01 = None

    if df01 is not None and meta_df is not None:
        extra = meta_df.reindex(df01.index)
        if labels is not None:
            labels = labels.reindex(df01.index)

    if args.rarefaction:
        if df01 is None:
            skipped.append("rarefaction (overall): missing binary_presence_absence.fasta")
        else:
            curves = compute_rarefaction_curves(df01, steps=args.rare_steps, reps=args.rare_reps, seed=args.rare_seed)
            out_csv = os.path.join(outdir, "pangenome_rarefaction.csv")
            out_png = os.path.join(outdir, "pangenome_rarefaction.png")
            curves.to_csv(out_csv, index=False)
            plot_rarefaction(curves, out_png, palette,
                             title="Gene family accumulation (Campylobacter PIRATE)",
                             legend_outside=args.legend_outside)
            ran.append(f"rarefaction (overall): {os.path.basename(out_csv)}, {os.path.basename(out_png)}")

            if args.rarefaction_by:
                if meta_df is None:
                    skipped.append(f"rarefaction by '{args.rarefaction_by}': missing/invalid --meta")
                elif args.rarefaction_by not in meta_df.columns:
                    skipped.append(f"rarefaction by '{args.rarefaction_by}': column not in metadata")
                else:
                    grp_series = meta_df[args.rarefaction_by].reindex(df01.index)
                    curves_by = compute_rarefaction_curves_by_group(
                        df01,
                        grp_series,
                        steps=args.rare_steps,
                        reps=args.rare_reps,
                        seed=args.rare_seed,
                        min_group_size=args.min_group_size,
                        max_groups=args.max_groups,
                    )
                    out_png2 = os.path.join(outdir, f"pangenome_rarefaction_by_{args.rarefaction_by}.png")
                    plot_rarefaction_by_group(
                        curves_by,
                        out_png2,
                        palette,
                        group_colors=group_colors,
                        title=f"Gene family accumulation by {args.rarefaction_by}",
                        include_core=args.rarefaction_core,
                        legend_outside=args.legend_outside,
                    )
                    ran.append(f"rarefaction by '{args.rarefaction_by}': {os.path.basename(out_png2)}")
    else:
        skipped.append("rarefaction (overall): flag not set")

    if args.pangenome_bars:
        if gf is None:
            skipped.append("pangenome composition bars: missing PIRATE.gene_families.tsv")
        else:
            n_genomes = df01.shape[0] if df01 is not None else 0
            try:
                bins = compute_pangenome_bins_from_gene_families(gf, n_genomes)
                out_png = os.path.join(outdir, "pangenome_composition.png")
                out_counts = os.path.join(outdir, "pangenome_composition_counts.csv")
                plot_pangenome_bars(bins, n_genomes, out_png, palette)
                pd.Series(bins).to_csv(out_counts)
                ran.append(f"pangenome composition bars: {os.path.basename(out_png)}, {os.path.basename(out_counts)}")
            except Exception as e:
                warn(f"Failed to compute pangenome bins: {e}")
                skipped.append("pangenome composition bars: parse/compute error")
    else:
        skipped.append("pangenome composition bars: flag not set")

    if args.gene_freq:
        if gf is None:
            skipped.append("gene frequency histogram: missing PIRATE.gene_families.tsv")
        else:
            try:
                out_png = os.path.join(outdir, "gene_frequency_hist.png")
                plot_gene_frequency_hist(gf, out_png)
                ran.append(f"gene frequency histogram: {os.path.basename(out_png)}")
            except Exception as e:
                warn(f"Failed to plot gene frequency histogram: {e}")
                skipped.append("gene frequency histogram: plot error")
    else:
        skipped.append("gene frequency histogram: flag not set")

    if df01 is None:
        if args.gene_counts:
            skipped.append("gene count per genome: missing binary_presence_absence.fasta")
        if args.heatmap_top is not None:
            skipped.append("top variable genes heatmap: missing binary_presence_absence.fasta")
        if args.pca:
            skipped.append("accessory PCA: missing binary_presence_absence.fasta")
        if args.mds:
            skipped.append("accessory MDS: missing binary_presence_absence.fasta")
        if args.umap:
            skipped.append("accessory UMAP: missing binary_presence_absence.fasta")
        if args.pca_gene_count:
            skipped.append("accessory PCA coloured by gene count: missing binary_presence_absence.fasta")
    else:
        gene_counts = gene_counts_per_genome(df01)

        if args.gene_counts:
            out_png = os.path.join(outdir, "gene_count_per_genome.png")
            plot_gene_count_distribution(gene_counts, out_png)
            ran.append(f"gene count per genome: {os.path.basename(out_png)}")
        else:
            skipped.append("gene count per genome: flag not set")

        if args.heatmap_top is not None:
            out_png = os.path.join(outdir, f"top_{args.heatmap_top}_variable_genes_heatmap.png")
            plot_top_variance_heatmap(df01, top_n=args.heatmap_top, out_png=out_png)
            ran.append(f"top variable genes heatmap: {os.path.basename(out_png)}")
        else:
            skipped.append("top variable genes heatmap: flag not set")

        if args.pca:
            coords, var = run_pca_accessory(df01)
            out_csv = os.path.join(outdir, "accessory_pca_coords.csv")
            save_coords_csv(coords, df01.index.tolist(), out_csv, extra=extra)
            out_png = os.path.join(outdir, "accessory_pca.png")
            plot_scatter_pretty(
                coords,
                labels=labels,
                palette=palette,
                title=f"Accessory genome PCA (PC1 {var[0]*100:.1f}%, PC2 {var[1]*100:.1f}%)",
                out_png=out_png,
                xlabel="PC1",
                ylabel="PC2",
                group_colors=group_colors,
                legend_outside=args.legend_outside,
                knn_edges=args.knn_edges,
                legend_out_png=args.legend_out,
            )
            ran.append(f"accessory PCA: {os.path.basename(out_csv)}, {os.path.basename(out_png)}")

            if args.pca_gene_count:
                out_png2 = os.path.join(outdir, "accessory_pca_gene_count.png")
                plot_pca_vs_gene_count(coords, gene_counts, out_png2, (float(var[0]), float(var[1])))
                ran.append(f"accessory PCA coloured by gene count: {os.path.basename(out_png2)}")
            else:
                skipped.append("accessory PCA coloured by gene count: flag not set")
        else:
            skipped.append("accessory PCA: flag not set")
            if args.pca_gene_count:
                skipped.append("accessory PCA coloured by gene count: requires --pca")

        if args.mds:
            coords = run_mds_jaccard(df01)
            out_csv = os.path.join(outdir, "accessory_mds_coords.csv")
            save_coords_csv(coords, df01.index.tolist(), out_csv, extra=extra)
            out_png = os.path.join(outdir, "accessory_mds.png")
            plot_scatter_pretty(
                coords,
                labels=labels,
                palette=palette,
                title="Accessory genome MDS (Jaccard distance)",
                out_png=out_png,
                xlabel="MDS1",
                ylabel="MDS2",
                group_colors=group_colors,
                legend_outside=args.legend_outside,
                knn_edges=args.knn_edges,
                legend_out_png=args.legend_out,
            )
            ran.append(f"accessory MDS: {os.path.basename(out_csv)}, {os.path.basename(out_png)}")
        else:
            skipped.append("accessory MDS: flag not set")

        if args.umap:
            try:
                coords = run_umap(df01, n_neighbors=args.umap_n_neighbors, min_dist=args.umap_min_dist)
                out_csv = os.path.join(outdir, "accessory_umap_coords.csv")
                save_coords_csv(coords, df01.index.tolist(), out_csv, extra=extra)
                out_png = os.path.join(outdir, "accessory_umap.png")
                plot_scatter_pretty(
                    coords,
                    labels=labels,
                    palette=palette,
                    title=f"Accessory genome UMAP (Jaccard; n={args.umap_n_neighbors}, min_dist={args.umap_min_dist})",
                    out_png=out_png,
                    xlabel="UMAP1",
                    ylabel="UMAP2",
                    group_colors=group_colors,
                    legend_outside=args.legend_outside,
                    knn_edges=args.knn_edges,
                    legend_out_png=args.legend_out,
                )
                ran.append(f"accessory UMAP: {os.path.basename(out_csv)}, {os.path.basename(out_png)}")
            except Exception as e:
                warn(str(e))
                skipped.append("accessory UMAP: missing dependency (umap-learn) or runtime error")
        else:
            skipped.append("accessory UMAP: flag not set")

    print("\n" + "=" * 72)
    print("PIRATE plotting summary")
    print("=" * 72)

    print(f"Ran: {len(ran)}")
    for x in ran:
        print(f"  - {x}")

    seen = set()
    skipped_u = []
    for s in skipped:
        if s not in seen:
            skipped_u.append(s)
            seen.add(s)

    print(f"\nSkipped: {len(skipped_u)}")
    for x in skipped_u:
        print(f"  - {x}")

    print("=" * 72)
    print(f"Done. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
