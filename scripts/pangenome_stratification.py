#!/usr/bin/env python3
"""Group-aware pangenome analyses for PIRATE and Panaroo outputs.

This is the tool-neutral successor to ``pirate_host_association.py``.
All pangenome parsing is delegated to ``pangenome_io.load_pangenome``.

Features
--------
1. Per-group pangenome/core-genome rarefaction curves with 95% intervals.
2. Gene-family association testing using a Cochran-Mantel-Haenszel (CMH)
   comparison of each target group against all other groups, stratified by a
   population-structure variable such as clonal complex.

Canonical input after loading
-----------------------------
rows = samples, columns = gene families, values = 0/1.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from math import erf, sqrt
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from pangenome_io import PangenomeData, load_pangenome
except ImportError as exc:
    raise SystemExit(
        "Could not import pangenome_io.py. Place it beside this script or add "
        "its directory to PYTHONPATH."
    ) from exc


WES_SUPERSET = [
    "#EE9B00", "#9C2C2C", "#003B5C", "#005F73",
    "#0A9396", "#94D2BD", "#A8DADC", "#457B9D",
    "#E9D8A6", "#CA6702", "#BB3E03", "#AE2012",
    "#5B1A1A", "#D04E4E", "#F2B5B5", "#F7EDE2",
    "#3E4E50", "#BFD7EA", "#2A9D8F", "#264653",
]


def set_plot_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "font.size": 12,
        "axes.titlesize": 16,
        "axes.labelsize": 13,
        "legend.fontsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "lines.linewidth": 2.8,
    })


def ensure_outdir(outdir: str | os.PathLike[str]) -> Path:
    path = Path(outdir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_metadata(meta_path: str, id_col: str = "sample") -> pd.DataFrame:
    kwargs = dict(dtype=str, encoding="utf-8", encoding_errors="replace")
    if meta_path.lower().endswith(".csv"):
        meta = pd.read_csv(meta_path, **kwargs)
    else:
        meta = pd.read_csv(meta_path, sep="\t", **kwargs)

    if id_col not in meta.columns:
        raise ValueError(
            f"Metadata must contain ID column '{id_col}'. Found: {list(meta.columns)}"
        )
    if meta[id_col].duplicated().any():
        dup = meta.loc[meta[id_col].duplicated(), id_col].head(10).tolist()
        raise ValueError(f"Metadata contains duplicate sample IDs, e.g. {dup}")

    meta[id_col] = meta[id_col].astype(str)
    return meta.set_index(id_col)


def align_pangenome_and_metadata(
    df01: pd.DataFrame,
    meta: pd.DataFrame,
    minimum_samples: int = 10,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Intersect sample IDs while preserving pangenome matrix order."""
    df01 = df01.copy()
    df01.index = df01.index.astype(str)
    meta = meta.copy()
    meta.index = meta.index.astype(str)

    shared = df01.index[df01.index.isin(meta.index)]
    if len(shared) < minimum_samples:
        raise ValueError(
            f"Too few shared samples between pangenome and metadata: n={len(shared)}. "
            "Check ID formatting and --meta-id-col."
        )

    missing_meta = len(df01.index) - len(shared)
    missing_pg = len(meta.index) - len(shared)
    if missing_meta or missing_pg:
        print(
            f"Sample join: {len(shared)} shared; "
            f"{missing_meta} pangenome samples lack metadata; "
            f"{missing_pg} metadata samples lack pangenome data.",
            file=sys.stderr,
        )

    return df01.loc[shared], meta.reindex(shared)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR, retaining NaN positions."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out

    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    restored = np.empty_like(q)
    restored[order] = q
    out[valid] = restored
    return out


def filter_by_prevalence(
    df01: pd.DataFrame,
    min_prev: float,
    max_prev: float,
) -> pd.DataFrame:
    if not 0 <= min_prev <= max_prev <= 1:
        raise ValueError("Prevalence thresholds must satisfy 0 <= min <= max <= 1.")
    prevalence = df01.mean(axis=0)
    keep = prevalence.between(min_prev, max_prev, inclusive="both")
    out = df01.loc[:, keep]
    if out.shape[1] < 1:
        raise ValueError(
            "No gene families remain after prevalence filtering. "
            "Relax --min-prevalence/--max-prevalence."
        )
    return out


def rarefaction_curve(
    X: np.ndarray,
    steps: int = 60,
    reps: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """Compute pan/core accumulation without storing a full cumulative matrix.

    This implementation keeps only one feature-count vector in memory, which is
    important for pangenomes containing thousands of genomes and tens of
    thousands of families.
    """
    Xb = np.asarray(X, dtype=np.uint8)
    if Xb.ndim != 2 or Xb.shape[0] == 0:
        raise ValueError("Rarefaction requires a non-empty 2D matrix.")
    if reps < 1 or steps < 1:
        raise ValueError("Rarefaction reps and steps must be >= 1.")

    rng = np.random.default_rng(seed)
    n, n_features = Xb.shape
    ks = np.unique(np.round(np.linspace(1, n, min(steps, n))).astype(int))
    k_to_col = {int(k): i for i, k in enumerate(ks)}

    pan = np.zeros((reps, len(ks)), dtype=np.int32)
    core = np.zeros((reps, len(ks)), dtype=np.int32)
    count_dtype = np.uint16 if n <= np.iinfo(np.uint16).max else np.uint32

    for rep in range(reps):
        order = rng.permutation(n)
        counts = np.zeros(n_features, dtype=count_dtype)
        for position, row_index in enumerate(order, start=1):
            counts += Xb[row_index]
            col = k_to_col.get(position)
            if col is not None:
                pan[rep, col] = np.count_nonzero(counts)
                core[rep, col] = np.count_nonzero(counts == position)

    def summarise(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (
            matrix.mean(axis=0),
            np.percentile(matrix, 2.5, axis=0),
            np.percentile(matrix, 97.5, axis=0),
        )

    pan_mean, pan_lo, pan_hi = summarise(pan)
    core_mean, core_lo, core_hi = summarise(core)
    return pd.DataFrame({
        "n_genomes": ks,
        "pan_mean": pan_mean,
        "pan_lo": pan_lo,
        "pan_hi": pan_hi,
        "core_mean": core_mean,
        "core_lo": core_lo,
        "core_hi": core_hi,
    })


def rarefaction_per_group(
    df01: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    reps: int,
    steps: int,
    seed: int,
    min_group_n: int,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, int]]:
    if group_col not in meta.columns:
        raise ValueError(f"Metadata lacks group column '{group_col}'.")

    valid = meta[group_col].notna() & (meta[group_col].astype(str).str.strip() != "")
    matrix = df01.loc[valid]
    metadata = meta.loc[valid]

    counts = metadata[group_col].astype(str).value_counts()
    groups = counts[counts >= min_group_n].index.tolist()
    if not groups:
        raise ValueError(f"No groups in '{group_col}' have n >= {min_group_n}.")

    curves: Dict[str, pd.DataFrame] = {}
    sizes: Dict[str, int] = {}
    for offset, group in enumerate(groups):
        ids = metadata.index[metadata[group_col].astype(str) == str(group)]
        sizes[str(group)] = len(ids)
        curves[str(group)] = rarefaction_curve(
            matrix.loc[ids].to_numpy(copy=False),
            steps=steps,
            reps=reps,
            seed=seed + offset,
        )
    return curves, sizes


def write_and_plot_rarefaction(
    curves: Dict[str, pd.DataFrame],
    sizes: Dict[str, int],
    group_col: str,
    out_png: Path,
    out_tsv: Path,
    truncate_to_min: bool = False,
) -> None:
    set_plot_style()
    x_max: Optional[int] = min(sizes.values()) if truncate_to_min else None

    tables = []
    plt.figure(figsize=(8.5, 6.2))
    for i, (group, curve) in enumerate(curves.items()):
        c = WES_SUPERSET[i % len(WES_SUPERSET)]
        current = curve.copy()
        current.insert(0, "group_level", group)
        current.insert(0, "group_col", group_col)
        current["group_n"] = sizes[group]
        tables.append(current)

        plot_data = curve if x_max is None else curve[curve["n_genomes"] <= x_max]
        x = plot_data["n_genomes"].to_numpy()
        plt.plot(x, plot_data["pan_mean"], color=c, label=f"{group} (pan)", alpha=0.95)
        plt.fill_between(x, plot_data["pan_lo"], plot_data["pan_hi"], color=c, alpha=0.15)
        plt.plot(
            x,
            plot_data["core_mean"],
            color=c,
            linestyle="--",
            label=f"{group} (core)",
            alpha=0.95,
        )
        plt.fill_between(x, plot_data["core_lo"], plot_data["core_hi"], color=c, alpha=0.10)

    plt.xlabel(f"Number of genomes (within {group_col})")
    plt.ylabel("Number of gene families")
    title = f"Gene family accumulation by {group_col}"
    if x_max is not None:
        title += f" (truncated to n={x_max})"
    plt.title(title)
    plt.legend(ncol=2, frameon=True, loc="best")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

    pd.concat(tables, ignore_index=True).to_csv(out_tsv, sep="\t", index=False)


def cmh_for_feature(
    x: np.ndarray,
    y: np.ndarray,
    strata: np.ndarray,
) -> Tuple[float, float, int]:
    """Return common odds ratio, two-sided CMH p-value and informative strata."""
    or_num = 0.0
    or_den = 0.0
    sum_a_minus_e = 0.0
    var_sum = 0.0
    informative = 0

    for stratum in pd.unique(strata):
        mask = strata == stratum
        if mask.sum() < 2:
            continue
        xs = x[mask]
        ys = y[mask]

        a = float(np.count_nonzero((xs == 1) & (ys == 1)))
        b = float(np.count_nonzero((xs == 1) & (ys == 0)))
        c = float(np.count_nonzero((xs == 0) & (ys == 1)))
        d = float(np.count_nonzero((xs == 0) & (ys == 0)))
        n = a + b + c + d

        if n <= 1:
            continue
        if (a + c) == 0 or (b + d) == 0 or (a + b) == 0 or (c + d) == 0:
            continue

        informative += 1
        or_num += (a * d) / n
        or_den += (b * c) / n

        row1, row0 = a + b, c + d
        col1, col0 = a + c, b + d
        expected_a = (row1 * col1) / n
        variance_a = (row1 * row0 * col1 * col0) / (n**2 * (n - 1.0))
        sum_a_minus_e += a - expected_a
        var_sum += variance_a

    if or_den == 0.0:
        common_or = np.inf if or_num > 0 else np.nan
    else:
        common_or = or_num / or_den

    if var_sum <= 0 or informative == 0:
        return common_or, 1.0, informative

    corrected = max(abs(sum_a_minus_e) - 0.5, 0.0)
    chi2 = (corrected**2) / var_sum
    z = sqrt(max(chi2, 0.0))
    normal_cdf = 0.5 * (1.0 + erf(z / sqrt(2.0)))
    p_value = float(np.clip(2.0 * (1.0 - normal_cdf), 0.0, 1.0))
    return common_or, p_value, informative


def run_group_association(
    df01: pd.DataFrame,
    meta: pd.DataFrame,
    family_metadata: pd.DataFrame,
    group_col: str,
    structure_col: str,
    min_group_n: int,
    min_prev: float,
    max_prev: float,
) -> pd.DataFrame:
    for column in (group_col, structure_col):
        if column not in meta.columns:
            raise ValueError(f"Metadata lacks required column '{column}'.")

    valid = (
        meta[group_col].notna()
        & meta[structure_col].notna()
        & (meta[group_col].astype(str).str.strip() != "")
        & (meta[structure_col].astype(str).str.strip() != "")
    )
    matrix = filter_by_prevalence(df01.loc[valid], min_prev, max_prev)
    metadata = meta.loc[valid]

    group_counts = metadata[group_col].astype(str).value_counts()
    groups = group_counts[group_counts >= min_group_n].index.tolist()
    if not groups:
        raise ValueError(f"No groups in '{group_col}' have n >= {min_group_n}.")

    X = matrix.to_numpy(dtype=np.uint8, copy=False)
    features = matrix.columns.astype(str).tolist()
    strata = metadata[structure_col].astype(str).to_numpy()
    group_values = metadata[group_col].astype(str).to_numpy()
    prevalence = matrix.mean(axis=0).to_numpy()

    results = []
    for group in groups:
        y = (group_values == str(group)).astype(np.uint8)
        pvals = np.ones(len(features), dtype=float)
        odds_ratios = np.full(len(features), np.nan, dtype=float)
        n_strata = np.zeros(len(features), dtype=int)

        for j in range(len(features)):
            odds_ratios[j], pvals[j], n_strata[j] = cmh_for_feature(X[:, j], y, strata)

        result = pd.DataFrame({
            "group_col": group_col,
            "group_level": group,
            "structure_col": structure_col,
            "gene_family": features,
            "or_cmh": odds_ratios,
            "p_value": pvals,
            "q_value_bh": bh_fdr(pvals),
            "prevalence": prevalence,
            "n_informative_strata": n_strata,
            "n_total": len(y),
            "n_in_group": int(y.sum()),
        })
        results.append(result)

    output = pd.concat(results, ignore_index=True)

    annotation_columns = [
        c for c in ["gene_name", "annotation", "n_genomes", "source_tool"]
        if c in family_metadata.columns
    ]
    if annotation_columns:
        annotations = family_metadata.copy()
        if "gene_family" in annotations.columns:
            annotations = annotations.drop(columns=["gene_family"])
        annotations.index = annotations.index.astype(str)
        annotations.index.name = "gene_family"
        output = output.merge(
            annotations[annotation_columns].reset_index(),
            on="gene_family",
            how="left",
        )

    front = [
        "group_col", "group_level", "structure_col", "gene_family",
        "gene_name", "annotation", "or_cmh", "p_value", "q_value_bh",
        "prevalence", "n_genomes", "n_informative_strata", "n_total",
        "n_in_group", "source_tool",
    ]
    return output[[c for c in front if c in output.columns] + [c for c in output.columns if c not in front]]


def resolve_cli_input(args: argparse.Namespace, parser: argparse.ArgumentParser) -> Tuple[str, str]:
    """Resolve modern arguments and the legacy --pirate-out option."""
    if args.pirate_out:
        if args.pangenome_out or args.pangenome_tool:
            parser.error("Use either --pirate-out or --pangenome-tool/--pangenome-out, not both.")
        warnings.warn(
            "--pirate-out is deprecated; use --pangenome-tool pirate --pangenome-out PATH.",
            DeprecationWarning,
            stacklevel=2,
        )
        return "pirate", args.pirate_out

    if not args.pangenome_tool or not args.pangenome_out:
        parser.error("--pangenome-tool and --pangenome-out are required (or use legacy --pirate-out).")
    return args.pangenome_tool, args.pangenome_out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Group-aware rarefaction and CMH association for PIRATE/Panaroo pangenomes."
    )
    parser.add_argument("--pangenome-tool", choices=["pirate", "panaroo"])
    parser.add_argument("--pangenome-out")
    parser.add_argument("--pirate-out", help=argparse.SUPPRESS)
    parser.add_argument("--pirate-threshold", type=float, default=95.0)
    parser.add_argument(
        "--panaroo-from-csv",
        action="store_true",
        help="Ignore gene_presence_absence.Rtab and derive Panaroo matrix from CSV.",
    )

    parser.add_argument("--meta", required=True)
    parser.add_argument("--meta-id-col", default="sample")
    parser.add_argument("--outdir", default="pangenome_host_assoc")
    parser.add_argument("--group-col", default="host", help="Grouping variable, e.g. host or country.")
    parser.add_argument(
        "--structure-col",
        default="CC",
        help="Population-structure stratum, e.g. CC, ST, cgMLST or PopPUNK cluster.",
    )
    parser.add_argument("--min-group-n", type=int, default=30)
    parser.add_argument("--min-shared-samples", type=int, default=10,
                        help="Minimum pangenome/metadata ID overlap required (default: 10).")
    parser.add_argument("--min-prevalence", type=float, default=0.01)
    parser.add_argument("--max-prevalence", type=float, default=0.99)

    parser.add_argument("--rarefaction-per-group", action="store_true")
    parser.add_argument("--truncate-to-min", action="store_true")
    parser.add_argument("--rare-steps", type=int, default=60)
    parser.add_argument("--rare-reps", type=int, default=20)
    parser.add_argument("--rare-seed", type=int, default=0)
    parser.add_argument(
        "--assoc",
        action="store_true",
        help="Run CMH group-vs-rest association stratified by --structure-col.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    tool, pangenome_out = resolve_cli_input(args, parser)

    if not args.rarefaction_per_group and not args.assoc:
        parser.error("Select at least one analysis: --rarefaction-per-group and/or --assoc.")

    outdir = ensure_outdir(args.outdir)
    pangenome: PangenomeData = load_pangenome(
        tool,
        pangenome_out,
        pirate_threshold=args.pirate_threshold,
        prefer_panaroo_rtab=not args.panaroo_from_csv,
    )
    metadata = load_metadata(args.meta, id_col=args.meta_id_col)
    matrix, metadata = align_pangenome_and_metadata(
        pangenome.presence_absence,
        metadata,
        minimum_samples=args.min_shared_samples,
    )

    print(
        f"Loaded {matrix.shape[0]} samples and {matrix.shape[1]} gene families "
        f"from {pangenome.tool}.",
        file=sys.stderr,
    )

    if args.rarefaction_per_group:
        curves, sizes = rarefaction_per_group(
            matrix,
            metadata,
            group_col=args.group_col,
            reps=args.rare_reps,
            steps=args.rare_steps,
            seed=args.rare_seed,
            min_group_n=args.min_group_n,
        )
        write_and_plot_rarefaction(
            curves,
            sizes,
            group_col=args.group_col,
            out_png=outdir / f"rarefaction_by_{args.group_col}.png",
            out_tsv=outdir / f"rarefaction_by_{args.group_col}.tsv",
            truncate_to_min=args.truncate_to_min,
        )

    if args.assoc:
        association = run_group_association(
            matrix,
            metadata,
            pangenome.family_metadata,
            group_col=args.group_col,
            structure_col=args.structure_col,
            min_group_n=args.min_group_n,
            min_prev=args.min_prevalence,
            max_prev=args.max_prevalence,
        )
        association.to_csv(
            outdir / f"cmh_assoc_{args.group_col}_strat_{args.structure_col}.tsv",
            sep="\t",
            index=False,
        )

    print(f"Done. Outputs in: {outdir}")


if __name__ == "__main__":
    main()
