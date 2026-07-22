#!/usr/bin/env python3
"""Tool-neutral pangenome plotting for PIRATE and Panaroo outputs.

All source-format parsing is delegated to :mod:`pangenome_io`. Downstream
plotting functions consume a canonical samples-by-gene-family binary matrix,
standardised family metadata, and a tool-neutral pangenome summary.

Examples
--------
Panaroo::

    python pangenome_plotting.py \
      --pangenome-tool panaroo \
      --pangenome-out panaroo_output \
      --outdir plots --pca --gene-freq --pangenome-bars

PIRATE::

    python pangenome_plotting.py \
      --pangenome-tool pirate \
      --pangenome-out PIRATE_out \
      --outdir plots --all

The legacy ``--pirate-out`` option remains available for backwards
compatibility and implies ``--pangenome-tool pirate``.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances

# Allow execution from outside the source directory while keeping a simple
# single-directory installation model for PANOPTICON scripts.
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from pangenome_io import PangenomeData, load_pangenome

# -----------------------------
# Wes Anderson-inspired palettes
# -----------------------------
# Note: These are "Life Aquatic"-inspired blues/teals (not a pixel-perfect extraction).
# If you want exact hexes from a specific palette source, tell me which package/source you prefer.
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


def get_palette(name: str) -> List[str]:
    if name not in WES_PALETTES:
        raise ValueError(f"Unknown palette '{name}'. Available: {', '.join(WES_PALETTES.keys())}")
    return WES_PALETTES[name]


# -----------------------------
# IO helpers
# -----------------------------
def ensure_outdir(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)


def filter_by_counts(df01: pd.DataFrame, min_count: int = 2, max_count: Optional[int] = None) -> pd.DataFrame:
    n = df01.shape[0]
    if max_count is None:
        max_count = n - 2
    counts = df01.sum(axis=0)
    keep = (counts >= min_count) & (counts <= max_count)
    return df01.loc[:, keep]

def filter_variable_columns(df01: pd.DataFrame, min_var: float = 0.0) -> pd.DataFrame:
    """Drop gene families with zero variance (all 0 or all 1)."""
    v = df01.var(axis=0)
    keep = v > min_var
    return df01.loc[:, keep]

def load_metadata(meta_path: str, id_col: str = "sample") -> pd.DataFrame:
    """
    Loads a metadata TSV/CSV (auto-detected by extension).
    Must contain a column matching id_col that matches sample IDs in PIRATE outputs.
    """
    if meta_path.endswith(".csv"):
        df = pd.read_csv(meta_path)
    else:
        df = pd.read_csv(meta_path, sep="\t")
    if id_col not in df.columns:
        raise ValueError(f"Metadata file must contain column '{id_col}'. Found: {list(df.columns)}")
    df = df.set_index(id_col)
    return df


# -----------------------------
# Plot helpers
# -----------------------------
def save_coords_csv(coords: np.ndarray, samples: List[str], out_csv: str, extra: Optional[pd.DataFrame] = None) -> None:
    df = pd.DataFrame(coords, index=samples, columns=["dim1", "dim2"])
    if extra is not None:
        df = pd.concat([df, extra.reindex(df.index)], axis=1)
    df.to_csv(out_csv)


def plot_scatter(
    x: np.ndarray,
    y: np.ndarray,
    labels: Optional[pd.Series],
    palette: List[str],
    title: str,
    out_png: str,
    xlabel: str = "Dim 1",
    ylabel: str = "Dim 2",
    point_size: float = 8.0,
) -> None:
    plt.figure(figsize=(7, 6))

    if labels is None:
        plt.scatter(x, y, s=point_size, c=palette[0])
    else:
        cats = labels.astype("category")
        levels = list(cats.cat.categories)
        color_map = {lvl: palette[i % len(palette)] for i, lvl in enumerate(levels)}
        colors = cats.map(color_map)
        plt.scatter(x, y, s=point_size, c=colors)

        # Legend
        handles = []
        for lvl in levels:
            handles.append(plt.Line2D([0], [0], marker="o", linestyle="", label=str(lvl),
                                      markerfacecolor=color_map[lvl], markeredgecolor=color_map[lvl], markersize=6))
        plt.legend(handles=handles, title=labels.name or "group", loc="best", frameon=True)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_gene_frequency_hist(family_metadata: pd.DataFrame, out_png: str) -> None:
    """Plot the number of genomes containing each normalised gene family."""
    if "n_genomes" not in family_metadata.columns:
        raise ValueError(
            "Normalised family metadata has no 'n_genomes' column. "
            "Check pangenome_io.py parsing for this input."
        )
    freq = pd.to_numeric(family_metadata["n_genomes"], errors="coerce").dropna()
    if freq.empty:
        raise ValueError("No usable gene-family frequency values were found.")

    plt.figure(figsize=(7, 4.5))
    plt.hist(freq.astype(int), bins=60)
    plt.xlabel("Number of genomes containing gene family")
    plt.ylabel("Number of gene families")
    plt.title("Gene family frequency distribution")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_pangenome_bars(summary: PangenomeSummary, out_png: str, palette: List[str]) -> None:
    labels = []
    values = []

    for name, val in [("Core", summary.core), ("Soft-core", summary.soft_core), ("Shell", summary.shell), ("Cloud", summary.cloud)]:
        if val is not None:
            labels.append(name)
            values.append(val)

    if not values:
        raise ValueError("Could not parse core/soft-core/shell/cloud from PIRATE.pangenome_summary.txt. "
                         "If you paste the top ~60 lines, I can add a parser for your exact format.")

    plt.figure(figsize=(7, 4.5))
    colors = [palette[i % len(palette)] for i in range(len(values))]
    plt.bar(labels, values, color=colors)
    plt.ylabel("Number of gene families")
    plt.title("Pangenome composition")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


from typing import Dict, Optional

def compute_rarefaction_curves_by_group(
    df01: pd.DataFrame,
    groups: pd.Series,
    steps: int = 60,
    reps: int = 20,
    seed: int = 0,
    min_group_size: int = 30,
    max_groups: Optional[int] = None,
) -> Dict[str, pd.DataFrame]:
    
    """
    Compute rarefaction curves separately for each group in `groups`.
    Returns dict: {group_name: curves_df}
    """
    
    g = groups.reindex(df01.index)

    # drop missing group labels
    keep = g.notna()
    df01 = df01.loc[keep]
    g = g.loc[keep].astype(str)

    # order groups by size (desc), filter by min size
    sizes = g.value_counts()
    sizes = sizes[sizes >= min_group_size]
    if max_groups is not None:
        sizes = sizes.iloc[:max_groups]

    out: Dict[str, pd.DataFrame] = {}
    for grp in sizes.index:
        sub = df01.loc[g == grp]
        out[grp] = compute_rarefaction_curves(
            sub, steps=steps, reps=reps, seed=seed
        )
    return out

def compute_rarefaction_curves(
    df01: pd.DataFrame,
    steps: int = 60,
    reps: int = 20,
    seed: int = 0,
) -> pd.DataFrame:

    rng = np.random.default_rng(seed)
    n = df01.shape[0]

    # choose genome counts to evaluate (dense early, coarser later)
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

        # cumulative sums across genomes
        csum = np.cumsum(Xp, axis=0)  # shape: (n, genes)

        for j, k in enumerate(ks):
            # pan: genes seen at least once
            pan_mat[r, j] = int((csum[k - 1, :] > 0).sum())
            # core: genes present in all k genomes so far
            core_mat[r, j] = int((csum[k - 1, :] == k).sum())

    def summarise(mat: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        mean = mat.mean(axis=0)
        lo = np.percentile(mat, 2.5, axis=0)
        hi = np.percentile(mat, 97.5, axis=0)
        return mean, lo, hi

    pan_mean, pan_lo, pan_hi = summarise(pan_mat)
    core_mean, core_lo, core_hi = summarise(core_mat)

    out = pd.DataFrame({
        "n_genomes": ks,
        "pan_mean": pan_mean, "pan_lo": pan_lo, "pan_hi": pan_hi,
        "core_mean": core_mean, "core_lo": core_lo, "core_hi": core_hi,
    })
    return out


def plot_rarefaction(
    curves: pd.DataFrame,
    out_png: str,
    palette: List[str],
    title: str = "Gene family accumulation (PIRATE)",
) -> None:
    """
    Plot pan (increasing) and core (decreasing) rarefaction curves with 95% CI.
    """
    plt.figure(figsize=(7.5, 5.5))

    x = curves["n_genomes"].values

    # pick two blues from the palette
    pan_c = palette[1]   # teal-ish
    core_c = palette[0]  # deep navy

    # pangenome
    plt.plot(x, curves["pan_mean"].values, linewidth=2, color=pan_c, label="Pangenome (total)")
    plt.fill_between(x, curves["pan_lo"].values, curves["pan_hi"].values, alpha=0.20, color=pan_c)

    # core genome
    plt.plot(x, curves["core_mean"].values, linewidth=2, color=core_c, label="Core (present in all)")
    plt.fill_between(x, curves["core_lo"].values, curves["core_hi"].values, alpha=0.20, color=core_c)

    plt.xlabel("Number of genomes")
    plt.ylabel("Number of gene families")
    plt.title(title)
    plt.legend(frameon=True, loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_rarefaction_by_group(
    curves_by_group: Dict[str, pd.DataFrame],
    out_png: str,
    palette: List[str],
    title: str = "Gene family accumulation by group",
    include_core: bool = False,
) -> None:
    plt.figure(figsize=(8.5, 6.0))

    for i, (grp, dfc) in enumerate(curves_by_group.items()):
        c = palette[i % len(palette)]
        x = dfc["n_genomes"].values

        # pangenome
        plt.plot(x, dfc["pan_mean"].values, linewidth=3.0, color=c, label=grp)
        plt.fill_between(x, dfc["pan_lo"].values, dfc["pan_hi"].values, alpha=0.18, color=c)

        # optional core (dashed)
        if include_core:
            plt.plot(x, dfc["core_mean"].values, linewidth=2.2, color=c, linestyle="--", alpha=0.9)

    plt.xlabel("Number of isolates")
    plt.ylabel("Number of gene families")
    plt.title(title)
    plt.legend(title="Group", frameon=True, loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def filter_accessory(df01: pd.DataFrame, min_freq: float = 0.01, max_freq: float = 0.99) -> pd.DataFrame:
    # freq of presence per gene family
    freq = df01.mean(axis=0)
    keep = (freq >= min_freq) & (freq <= max_freq)
    out = df01.loc[:, keep]
    # also drop any remaining zero-variance cols (belt + braces)
    out = out.loc[:, out.var(axis=0) > 0]
    return out


# -----------------------------
# Analyses
# -----------------------------
def run_pca_accessory(df01: pd.DataFrame, n_components: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """
    PCA on presence/absence matrix.
    Returns coords (n x 2) and explained variance ratios.
    """
    pca = PCA(n_components=n_components, random_state=0)
    coords = pca.fit_transform(df01.values)
    return coords, pca.explained_variance_ratio_


def run_mds_jaccard(df01: pd.DataFrame) -> np.ndarray:
    """
    Metric MDS on Jaccard distances of presence/absence (good for accessory structure).
    """
    # Jaccard distance expects boolean
    X = df01.values.astype(bool)
    dist = pairwise_distances(X, metric="jaccard")
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0, n_init=4, max_iter=300)
    coords = mds.fit_transform(dist)
    return coords


def run_umap(df01: pd.DataFrame, n_neighbors: int = 15, min_dist: float = 0.1) -> np.ndarray:
    try:
        import umap  # type: ignore
    except ImportError as e:
        raise ImportError("UMAP requested but umap-learn is not installed. Install with: conda install -c conda-forge umap-learn") from e

    reducer = umap.UMAP(n_components=2, n_neighbors=n_neighbors, min_dist=min_dist, metric="jaccard", random_state=0)
    coords = reducer.fit_transform(df01.values.astype(bool))
    return coords


# -----------------------------
# Extra diagnostics
# -----------------------------

def gene_counts_per_genome(df01: pd.DataFrame) -> pd.Series:
    """
    Number of gene families present per genome.
    """
    counts = df01.sum(axis=1)
    counts.name = "gene_count"
    return counts


def plot_gene_count_distribution(gene_counts: pd.Series, out_png: str) -> None:
    plt.figure(figsize=(7, 4.5))
    plt.hist(gene_counts, bins=50)
    plt.xlabel("Number of gene families per genome")
    plt.ylabel("Number of genomes")
    plt.title("Gene content per genome")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_pca_vs_gene_count(
    coords: np.ndarray,
    gene_counts: pd.Series,
    out_png: str,
    var: Tuple[float, float],
) -> None:
    plt.figure(figsize=(7, 6))
    sc = plt.scatter(
        coords[:, 0],
        coords[:, 1],
        c=gene_counts.values,
        cmap="viridis",
        s=10
    )
    plt.xlabel(f"PC1 ({var[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({var[1]*100:.1f}%)")
    plt.title("Accessory PCA coloured by gene count")
    plt.colorbar(sc, label="Genes per genome")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_top_variance_heatmap(
    df01: pd.DataFrame,
    top_n: int,
    out_png: str,
) -> None:
    """
    Heatmap of top N most variable gene families.
    """
    variances = df01.var(axis=0)
    top_cols = variances.sort_values(ascending=False).head(top_n).index
    mat = df01[top_cols]

    plt.figure(figsize=(12, 6))
    plt.imshow(mat.T, aspect="auto", interpolation="nearest", cmap="viridis")
    plt.colorbar(label="Presence (0/1)")
    plt.xlabel("Genomes")
    plt.ylabel(f"Top {top_n} variable gene families")
    plt.title(f"Top {top_n} variable gene families (presence/absence)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def fst_per_gene_nei(df01: pd.DataFrame, groups: pd.Series, min_group_size: int = 20) -> pd.DataFrame:
    """
    df01: rows=samples, cols=genes, values 0/1
    groups: index aligned to df01.index, values=group labels
    Returns: DataFrame with columns: gene, fst, ht, hs, n_groups, n_total
    """
    g = groups.reindex(df01.index)
    keep = g.notna()
    X = df01.loc[keep]
    g = g.loc[keep].astype(str)

    # keep only groups with enough samples
    sizes = g.value_counts()
    good = sizes[sizes >= min_group_size].index
    keep2 = g.isin(good)
    X = X.loc[keep2]
    g = g.loc[keep2]

    if len(good) < 2:
        raise ValueError("Need >=2 groups passing --fst-min-group-size to compute FST.")

    # weights
    sizes = g.value_counts().sort_index()
    w = (sizes / sizes.sum()).to_dict()

    # per-group frequencies p_i (vector per gene)
    # build p_i matrix: rows=groups, cols=genes
    p = []
    for grp in sizes.index:
        p.append(X.loc[g == grp].mean(axis=0).values)
    P = np.vstack(p)  # shape (G, M)
    W = np.array([w[grp] for grp in sizes.index]).reshape(-1, 1)

    # Hs: weighted mean within-group heterozygosity
    H_i = 2 * P * (1 - P)
    Hs = (W * H_i).sum(axis=0)

    # Ht: heterozygosity of pooled frequency
    pbar = (W * P).sum(axis=0)
    Ht = 2 * pbar * (1 - pbar)

    # fst
    fst = np.zeros_like(Ht)
    mask = Ht > 0
    fst[mask] = (Ht[mask] - Hs[mask]) / Ht[mask]
    fst = np.clip(fst, 0, 1)

    out = pd.DataFrame({
        "gene": X.columns,
        "fst": fst,
        "ht": Ht,
        "hs": Hs,
        "n_groups": len(sizes),
        "n_total": int(sizes.sum()),
    })
    return out.sort_values("fst", ascending=False).reset_index(drop=True)


def fst_pairwise_per_gene_nei(df01: pd.DataFrame, groups: pd.Series, min_group_size: int = 20) -> pd.DataFrame:
    """
    Returns a long table: gene, group_a, group_b, fst
    """
    g = groups.reindex(df01.index)
    keep = g.notna()
    X = df01.loc[keep]
    g = g.loc[keep].astype(str)

    sizes = g.value_counts()
    good = sizes[sizes >= min_group_size].index
    X = X.loc[g.isin(good)]
    g = g.loc[g.isin(good)]

    labs = sorted(g.unique())
    rows = []
    for i in range(len(labs)):
        for j in range(i + 1, len(labs)):
            a, b = labs[i], labs[j]
            Xa = X.loc[g == a]
            Xb = X.loc[g == b]
            pa = Xa.mean(axis=0).values
            pb = Xb.mean(axis=0).values

            wa = Xa.shape[0] / (Xa.shape[0] + Xb.shape[0])
            wb = 1 - wa

            Hs = wa * (2 * pa * (1 - pa)) + wb * (2 * pb * (1 - pb))
            pbar = wa * pa + wb * pb
            Ht = 2 * pbar * (1 - pbar)

            fst = np.zeros_like(Ht)
            mask = Ht > 0
            fst[mask] = (Ht[mask] - Hs[mask]) / Ht[mask]
            fst = np.clip(fst, 0, 1)

            rows.append(pd.DataFrame({"gene": X.columns, "group_a": a, "group_b": b, "fst": fst}))
    return pd.concat(rows, ignore_index=True)



# -----------------------------
# CLI
# -----------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Plot pangenome summaries and accessory-genome structure from PIRATE or Panaroo."
    )
    ap.add_argument("--pangenome-tool", choices=["pirate", "panaroo"],
                    help="Pangenome source format.")
    ap.add_argument("--pangenome-out",
                    help="Path to the PIRATE or Panaroo output directory.")
    ap.add_argument("--pirate-out",
                    help="Deprecated alias for --pangenome-tool pirate --pangenome-out PATH.")
    ap.add_argument("--pirate-threshold", type=float, default=95,
                    help="PIRATE clustering threshold to select when outputs contain multiple thresholds (default 95).")
    ap.add_argument("--panaroo-from-csv", action="store_true",
                    help="Derive Panaroo presence/absence from gene_presence_absence.csv even when Rtab exists.")

    ap.add_argument("--outdir", default="pangenome_plots", help="Output directory for plots")
    ap.add_argument("--palette", default="life_aquatic_blues",
                    help=f"Palette name: {', '.join(WES_PALETTES.keys())}")

    ap.add_argument("--all", action="store_true")
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
                    help="Metadata column used to stratify rarefaction curves. Requires --meta.")
    ap.add_argument("--min-group-size", type=int, default=30)
    ap.add_argument("--max-groups", type=int, default=None)
    ap.add_argument("--rarefaction-core", action="store_true")

    ap.add_argument("--meta", default=None)
    ap.add_argument("--meta-id-col", default="sample")
    ap.add_argument("--meta-col", default=None)

    ap.add_argument("--umap-n-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.1)

    ap.add_argument("--fst", action="store_true")
    ap.add_argument("--fst-by", default=None)
    ap.add_argument("--fst-min-group-size", type=int, default=20)
    ap.add_argument("--fst-min-freq", type=float, default=0.01)
    ap.add_argument("--fst-max-freq", type=float, default=0.99)
    ap.add_argument("--fst-pairwise", action="store_true")
    args = ap.parse_args()

    # Backwards-compatible input resolution.
    if args.pirate_out:
        if args.pangenome_out and Path(args.pangenome_out).resolve() != Path(args.pirate_out).resolve():
            ap.error("Use either --pirate-out or --pangenome-out, not two different directories.")
        if args.pangenome_tool and args.pangenome_tool != "pirate":
            ap.error("--pirate-out cannot be combined with --pangenome-tool panaroo.")
        args.pangenome_tool = "pirate"
        args.pangenome_out = args.pirate_out
        print("WARNING: --pirate-out is deprecated; use --pangenome-tool pirate --pangenome-out PATH.",
              file=sys.stderr)

    if not args.pangenome_tool or not args.pangenome_out:
        ap.error("--pangenome-tool and --pangenome-out are required (or use legacy --pirate-out).")

    outdir = args.outdir
    ensure_outdir(outdir)
    palette = get_palette(args.palette)

    if args.all:
        args.pangenome_bars = True
        args.gene_freq = True
        args.pca = True
        args.mds = True
        args.gene_counts = True
        args.pca_gene_count = True
        args.heatmap_top = 200
        args.rarefaction = True

    pg: PangenomeData = load_pangenome(
        args.pangenome_tool,
        args.pangenome_out,
        pirate_threshold=args.pirate_threshold,
        prefer_panaroo_rtab=not args.panaroo_from_csv,
    )
    family_metadata = pg.family_metadata
    summary = pg.summary

    print(
        f"Loaded {pg.tool}: {pg.presence_absence.shape[0]} genomes x "
        f"{pg.presence_absence.shape[1]} gene families",
        file=sys.stderr,
    )

    meta_df = None
    labels = None
    extra = None
    if args.meta:
        meta_df = load_metadata(args.meta, id_col=args.meta_id_col)
    if args.meta_col:
        if meta_df is None:
            raise ValueError("--meta-col requires --meta")
        if args.meta_col not in meta_df.columns:
            raise ValueError(f"Metadata column '{args.meta_col}' not found. Available: {list(meta_df.columns)}")
        labels = meta_df[args.meta_col]
        labels.name = args.meta_col

    need_df01 = any([
        args.pca, args.mds, args.umap, args.gene_counts, args.pca_gene_count,
        args.heatmap_top is not None, args.rarefaction, args.rarefaction_by is not None,
        args.fst, args.fst_pairwise,
    ])
    df01 = pg.presence_absence if need_df01 else None
    if df01 is not None:
        if meta_df is not None:
            extra = meta_df.reindex(df01.index)
            if labels is not None:
                labels = labels.reindex(df01.index)
        df01 = filter_variable_columns(df01)
        if df01.shape[1] < 2:
            raise ValueError("After filtering, fewer than two variable gene families remain.")

    plt.rcParams.update({
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "font.size": 12,
        "axes.titlesize": 16,
        "axes.labelsize": 13,
        "legend.fontsize": 11,
        "legend.framealpha": 0.95,
        "savefig.bbox": "tight",
    })

    if args.rarefaction:
        assert df01 is not None
        curves = compute_rarefaction_curves(
            df01, steps=args.rare_steps, reps=args.rare_reps, seed=args.rare_seed
        )
        curves.to_csv(os.path.join(outdir, "pangenome_rarefaction.csv"), index=False)
        plot_rarefaction(
            curves,
            os.path.join(outdir, "pangenome_rarefaction.png"),
            palette,
            title=f"Gene family accumulation ({pg.tool})",
        )

    if args.rarefaction_by:
        assert df01 is not None
        if meta_df is None:
            raise ValueError("--rarefaction-by requires --meta")
        if args.rarefaction_by not in meta_df.columns:
            raise ValueError(f"--rarefaction-by '{args.rarefaction_by}' not found in metadata columns")
        grp_series = meta_df[args.rarefaction_by].reindex(df01.index)
        curves_by = compute_rarefaction_curves_by_group(
            df01, grp_series, steps=args.rare_steps, reps=args.rare_reps,
            seed=args.rare_seed, min_group_size=args.min_group_size,
            max_groups=args.max_groups,
        )
        curves_by.to_csv(
            os.path.join(outdir, f"pangenome_rarefaction_by_{args.rarefaction_by}.csv"),
            index=False,
        )
        plot_rarefaction_by_group(
            curves_by,
            os.path.join(outdir, f"pangenome_rarefaction_by_{args.rarefaction_by}.png"),
            palette,
            title=f"Gene family accumulation by {args.rarefaction_by}",
            include_core=args.rarefaction_core,
        )

    if args.fst or args.fst_pairwise:
        assert df01 is not None
        if meta_df is None or args.fst_by is None:
            raise ValueError("--fst/--fst-pairwise require --meta and --fst-by")
        if args.fst_by not in meta_df.columns:
            raise ValueError(f"--fst-by '{args.fst_by}' not found in metadata columns")
        grp = meta_df[args.fst_by].reindex(df01.index)
        freq = df01.mean(axis=0)
        df_fst = df01.loc[:, (freq >= args.fst_min_freq) & (freq <= args.fst_max_freq)]
        if args.fst:
            fst_table = fst_per_gene_nei(df_fst, grp, min_group_size=args.fst_min_group_size)
            fst_table.to_csv(
                os.path.join(outdir, f"fst_per_gene_by_{args.fst_by}.tsv"), sep="\t", index=False
            )
        if args.fst_pairwise:
            fst_pw = fst_pairwise_per_gene_nei(df_fst, grp, min_group_size=args.fst_min_group_size)
            fst_pw.to_csv(
                os.path.join(outdir, f"fst_pairwise_per_gene_by_{args.fst_by}.tsv"), sep="\t", index=False
            )

    if args.pangenome_bars:
        bins = {
            "Core (≥99%)": summary.core,
            "Soft-core (95–99%)": summary.soft_core,
            "Shell (15–95%)": summary.shell,
            "Cloud (<15%)": summary.cloud,
        }
        bins = {k: int(v) for k, v in bins.items() if v is not None}
        if not bins:
            raise ValueError("No pangenome frequency-bin counts were available.")
        plt.figure(figsize=(7, 4.5))
        colors = [palette[i % len(palette)] for i in range(len(bins))]
        plt.bar(list(bins), list(bins.values()), color=colors)
        plt.ylabel("Number of gene families")
        plt.title(f"Pangenome composition (n={summary.n_genomes} genomes)")
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "pangenome_composition.png"), dpi=300)
        plt.close()
        pd.Series(bins, name="n_gene_families").to_csv(
            os.path.join(outdir, "pangenome_composition_counts.csv"), header=True
        )

    if args.gene_freq:
        plot_gene_frequency_hist(family_metadata, os.path.join(outdir, "gene_frequency_hist.png"))

    if df01 is not None:
        gene_counts = gene_counts_per_genome(df01)
        if args.gene_counts:
            plot_gene_count_distribution(gene_counts, os.path.join(outdir, "gene_count_per_genome.png"))
        if args.heatmap_top is not None:
            plot_top_variance_heatmap(
                df01, top_n=args.heatmap_top,
                out_png=os.path.join(outdir, f"top_{args.heatmap_top}_variable_genes_heatmap.png"),
            )
        if args.pca:
            coords, var = run_pca_accessory(df01)
            save_coords_csv(coords, df01.index.tolist(), os.path.join(outdir, "accessory_pca_coords.csv"), extra=extra)
            plot_scatter(
                coords[:, 0], coords[:, 1], labels, palette,
                title=f"Accessory genome PCA (PC1 {var[0]*100:.1f}%, PC2 {var[1]*100:.1f}%)",
                out_png=os.path.join(outdir, "accessory_pca.png"), xlabel="PC1", ylabel="PC2",
            )
            if args.pca_gene_count:
                plot_pca_vs_gene_count(
                    coords, gene_counts, os.path.join(outdir, "accessory_pca_gene_count.png"),
                    (float(var[0]), float(var[1])),
                )
        if args.mds:
            coords = run_mds_jaccard(df01)
            save_coords_csv(coords, df01.index.tolist(), os.path.join(outdir, "accessory_mds_coords.csv"), extra=extra)
            plot_scatter(
                coords[:, 0], coords[:, 1], labels, palette,
                title="Accessory genome MDS (Jaccard distance)",
                out_png=os.path.join(outdir, "accessory_mds.png"), xlabel="MDS1", ylabel="MDS2",
            )
        if args.umap:
            coords = run_umap(df01, n_neighbors=args.umap_n_neighbors, min_dist=args.umap_min_dist)
            save_coords_csv(coords, df01.index.tolist(), os.path.join(outdir, "accessory_umap_coords.csv"), extra=extra)
            plot_scatter(
                coords[:, 0], coords[:, 1], labels, palette,
                title=f"Accessory genome UMAP (Jaccard; n={args.umap_n_neighbors}, min_dist={args.umap_min_dist})",
                out_png=os.path.join(outdir, "accessory_umap.png"), xlabel="UMAP1", ylabel="UMAP2",
            )

    print(f"Done. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
