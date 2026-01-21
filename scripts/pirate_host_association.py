#!/usr/bin/env python3
"""
pirate_host_association.py

Host-/country-aware analyses for PIRATE presence/absence matrices.

Features (v1):
  1) Rarefaction / accumulation curves per group (e.g. Host) with 95% CI
  2) Host association testing for gene families using CMH stratified by structure (default: CC)
     - outputs odds ratio + p + FDR per gene-family
     - controls for population structure by stratifying within CC (or ST/cgMLST cluster)

Inputs:
  - PIRATE_out/binary_presence_absence.fasta
  - meta.tsv (must include: sample, host, country, CC at minimum)

Example:
  python pirate_host_association.py \
    --pirate-out PIRATE_out \
    --meta meta.tsv \
    --outdir pirate_host_assoc \
    --group-col host \
    --structure-col CC \
    --rarefaction-per-group \
    --assoc \
    --min-prevalence 0.01 --max-prevalence 0.99 \
    --min-group-n 30

Notes:
  - binary_presence_absence.fasta contains A/C only (A=0, C=1). This script decodes that.
  - CMH is run as: "group == target_level" vs "all other groups"
"""

from __future__ import annotations

import argparse
import os
import re
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO


# -----------------------------
# Styling + palette
# -----------------------------
WES_SUPERSET = [
    # A broad “Wes-ish” mix (blues/teals + warm accents + muted neutrals)
    "#EE9B00", "#9C2C2C", "#003B5C", "#005F73",
    "#0A9396", "#94D2BD", "#A8DADC", "#457B9D",
    "#E9D8A6", "#CA6702", "#BB3E03", "#AE2012",
    "#5B1A1A", "#D04E4E", "#F2B5B5", "#F7EDE2",
    "#3E4E50", "#BFD7EA", "#2A9D8F", "#264653",
]

HEX_RE = re.compile(r"^#(?:[0-9a-fA-F]{6})$")


def set_plot_style():
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
        "lines.linewidth": 2.8,  # thicker lines globally
    })


def parse_group_colors(spec: Optional[str]) -> Dict[str, str]:
    """
    Parse 'LEVEL=#RRGGBB,LEVEL2=#RRGGBB' into a dict.
    Whitespace is tolerated around tokens.

    Example:
      --group-colors "Human (gastroenteritis)=#E66101,Chicken (commercial)=#E6AB02"
    """
    if not spec:
        return {}
    out: Dict[str, str] = {}
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(
                f"Invalid --group-colors entry '{chunk}'. Expected format LEVEL=#RRGGBB"
            )
        level, hexcol = chunk.split("=", 1)
        level = level.strip()
        hexcol = hexcol.strip()
        if not level:
            raise ValueError(f"Invalid --group-colors entry '{chunk}': empty level name")
        if not HEX_RE.match(hexcol):
            raise ValueError(
                f"Invalid HEX colour '{hexcol}' for level '{level}'. Expected #RRGGBB"
            )
        out[level] = hexcol
    return out


def pick_group_color(group_level: str, i: int, group_color_map: Dict[str, str]) -> str:
    """
    If an explicit mapping exists for this group level, use it; otherwise use WES_SUPERSET.
    """
    if group_level in group_color_map:
        return group_color_map[group_level]
    return WES_SUPERSET[i % len(WES_SUPERSET)]


# -----------------------------
# IO
# -----------------------------
def ensure_outdir(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)


def load_metadata(meta_path: str, id_col: str = "sample") -> pd.DataFrame:
    read_kwargs = dict(dtype=str)  # keep IDs stable
    if meta_path.endswith(".csv"):
        meta = pd.read_csv(meta_path, encoding="utf-8", encoding_errors="replace", **read_kwargs)
    else:
        # .tsv, .txt, etc assumed tab-delimited unless user provides csv
        meta = pd.read_csv(meta_path, sep="\t", encoding="utf-8", encoding_errors="replace", **read_kwargs)

    if id_col not in meta.columns:
        raise ValueError(f"Metadata must contain '{id_col}'. Found: {list(meta.columns)}")
    return meta.set_index(id_col)


def read_binary_presence_absence_fasta(path: str) -> pd.DataFrame:
    """
    PIRATE binary_presence_absence.fasta:
      - records are genomes
      - sequence alphabet: A/C (A=0, C=1)

    Returns:
      df01 rows=samples, cols=positions (gene families), values 0/1
      (Note: without PIRATE mapping file, columns are just positions; still fine for stats/plots.)
    """
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found: {path}")

    ids = []
    mat = []
    seqlen = None
    for r in records:
        s = str(r.seq).upper().strip()
        if seqlen is None:
            seqlen = len(s)
        elif len(s) != seqlen:
            raise ValueError("FASTA sequences have inconsistent lengths; cannot build matrix.")
        # PIRATE: A=0, C=1
        arr = np.fromiter((1 if c == "C" else 0 for c in s), dtype=np.uint8)
        ids.append(r.id)
        mat.append(arr)

    X = np.vstack(mat)  # (n_samples, n_features)
    cols = [f"gf_{i+1}" for i in range(X.shape[1])]
    df01 = pd.DataFrame(X, index=ids, columns=cols)
    df01.index.name = "sample"
    return df01


# -----------------------------
# Helpers
# -----------------------------
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR (q-values)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def filter_by_prevalence(df01: pd.DataFrame, min_prev: float, max_prev: float) -> pd.DataFrame:
    prev = df01.mean(axis=0).values
    keep = (prev >= min_prev) & (prev <= max_prev)
    out = df01.loc[:, keep]
    if out.shape[1] < 2:
        raise ValueError(
            "After prevalence filtering, <2 features remain. "
            "Try relaxing --min-prevalence/--max-prevalence."
        )
    return out


# -----------------------------
# Rarefaction per group
# -----------------------------
def rarefaction_curve(
    X: np.ndarray,
    steps: int = 60,
    reps: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """
    X: boolean/0-1 matrix (n_genomes, n_features)

    Returns: n_genomes, pan_mean/lo/hi, core_mean/lo/hi
    """
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    ks = np.unique(np.round(np.linspace(1, n, min(steps, n))).astype(int))

    pan = np.zeros((reps, len(ks)), dtype=np.int32)
    core = np.zeros((reps, len(ks)), dtype=np.int32)

    Xb = X.astype(np.uint8)

    for r in range(reps):
        order = rng.permutation(n)
        Xp = Xb[order, :]
        csum = np.cumsum(Xp, axis=0)
        for j, k in enumerate(ks):
            pan[r, j] = int((csum[k - 1, :] > 0).sum())
            core[r, j] = int((csum[k - 1, :] == k).sum())

    def summarise(M):
        return M.mean(axis=0), np.percentile(M, 2.5, axis=0), np.percentile(M, 97.5, axis=0)

    pan_mean, pan_lo, pan_hi = summarise(pan)
    core_mean, core_lo, core_hi = summarise(core)

    return pd.DataFrame({
        "n_genomes": ks,
        "pan_mean": pan_mean, "pan_lo": pan_lo, "pan_hi": pan_hi,
        "core_mean": core_mean, "core_lo": core_lo, "core_hi": core_hi,
    })


def plot_rarefaction_per_group(
    df01: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    out_png: str,
    reps: int,
    steps: int,
    seed: int,
    min_group_n: int,
    truncate_to_min: bool = False,
    group_color_map: Optional[Dict[str, str]] = None,
) -> None:
    set_plot_style()
    group_color_map = group_color_map or {}

    # align
    meta2 = meta.reindex(df01.index)
    ok = meta2[group_col].notna()
    df01 = df01.loc[ok]
    meta2 = meta2.loc[ok]

    groups = meta2[group_col].astype("category").cat.categories.tolist()
    groups = [g for g in groups if (meta2[group_col] == g).sum() >= min_group_n]
    if not groups:
        raise ValueError(f"No groups in '{group_col}' with n >= {min_group_n}")

    # compute curves
    curves: Dict[str, pd.DataFrame] = {}
    ns = []
    for g in groups:
        idx = meta2.index[meta2[group_col] == g]
        X = df01.loc[idx].values
        curves[str(g)] = rarefaction_curve(X, steps=steps, reps=reps, seed=seed)
        ns.append(len(idx))

    x_max = min(ns) if truncate_to_min else None

    plt.figure(figsize=(8.5, 6.2))
    for i, g in enumerate(groups):
        g_str = str(g)
        c = pick_group_color(g_str, i, group_color_map)
        cur = curves[g_str].copy()

        if x_max is not None:
            cur = cur[cur["n_genomes"] <= x_max]

        x = cur["n_genomes"].values
        # Pangenome line
        plt.plot(x, cur["pan_mean"].values, color=c, label=f"{g_str} (pan)", alpha=0.95)
        plt.fill_between(x, cur["pan_lo"].values, cur["pan_hi"].values, color=c, alpha=0.15)

        # Core-within-host line (same colour, dashed)
        plt.plot(x, cur["core_mean"].values, color=c, linestyle="--", label=f"{g_str} (core)", alpha=0.95)
        plt.fill_between(x, cur["core_lo"].values, cur["core_hi"].values, color=c, alpha=0.10)

    plt.xlabel(f"Number of genomes (within {group_col})")
    plt.ylabel("Number of gene families")
    title = f"Gene family accumulation by {group_col}"
    if truncate_to_min:
        title += f" (truncated to n={min(ns)})"
    plt.title(title)
    plt.legend(ncol=2, frameon=True, loc="best")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


# -----------------------------
# CMH host association (stratified by structure)
# -----------------------------
def cmh_for_feature(
    x: np.ndarray,               # 0/1 feature vector (n,)
    y: np.ndarray,               # 0/1 target host vs rest (n,)
    strata: np.ndarray,          # stratum labels (n,)
) -> Tuple[float, float]:
    """
    Returns: (OR_cmh, p_value) for a single feature.

    Implements CMH chi-square test with continuity correction.
    Assumes 2x2 within each stratum.
    """
    uniq = pd.unique(strata)
    or_num = 0.0
    or_den = 0.0
    sum_a_minus_e = 0.0
    var_sum = 0.0

    for s in uniq:
        m = (strata == s)
        if m.sum() < 2:
            continue
        xs = x[m]
        ys = y[m]

        # 2x2 table:
        #           y=1   y=0
        # x=1        a     b
        # x=0        c     d
        a = float(((xs == 1) & (ys == 1)).sum())
        b = float(((xs == 1) & (ys == 0)).sum())
        c = float(((xs == 0) & (ys == 1)).sum())
        d = float(((xs == 0) & (ys == 0)).sum())
        n = a + b + c + d
        if n <= 1:
            continue

        # Skip strata with no variation in y or x
        if (a + c) == 0 or (b + d) == 0:
            continue
        if (a + b) == 0 or (c + d) == 0:
            continue

        # CMH OR components
        or_num += (a * d) / n
        or_den += (b * c) / n

        # CMH chi-square components
        row1 = a + b
        col1 = a + c
        e_a = (row1 * col1) / n

        row0 = c + d
        col0 = b + d
        var_a = (row1 * row0 * col1 * col0) / (n**2 * (n - 1.0)) if n > 1 else 0.0

        sum_a_minus_e += (a - e_a)
        var_sum += var_a

    # OR
    if or_den == 0.0:
        or_cmh = np.inf if or_num > 0 else np.nan
    else:
        or_cmh = or_num / or_den

    # p-value from chi-square(1), via normal approximation
    if var_sum <= 0:
        return or_cmh, 1.0

    num = (abs(sum_a_minus_e) - 0.5)**2  # continuity correction
    chi2 = num / var_sum

    from math import erf, sqrt
    z = np.sqrt(max(chi2, 0.0))
    phi = 0.5 * (1.0 + erf(z / sqrt(2.0)))
    p = 2.0 * (1.0 - phi)
    p = float(np.clip(p, 0.0, 1.0))
    return or_cmh, p


def run_host_association(
    df01: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    structure_col: str,
    min_group_n: int,
    out_csv: str,
    min_prev: float,
    max_prev: float,
) -> None:
    # align + filter
    meta2 = meta.reindex(df01.index)
    keep = meta2[group_col].notna() & meta2[structure_col].notna()
    df01 = df01.loc[keep]
    meta2 = meta2.loc[keep]

    # prevalence filter
    df01 = filter_by_prevalence(df01, min_prev=min_prev, max_prev=max_prev)

    # eligible target groups
    groups = meta2[group_col].astype("category").cat.categories.tolist()
    groups = [g for g in groups if (meta2[group_col] == g).sum() >= min_group_n]
    if not groups:
        raise ValueError(f"No groups in '{group_col}' with n >= {min_group_n}")

    strata = meta2[structure_col].astype(str).values

    results = []
    X = df01.values.astype(np.uint8)
    features = df01.columns.tolist()

    for g in groups:
        y = (meta2[group_col].astype(str).values == str(g)).astype(np.uint8)

        if y.sum() < min_group_n:
            continue

        pvals = np.ones(len(features), dtype=float)
        ors = np.full(len(features), np.nan, dtype=float)

        for j in range(len(features)):
            or_cmh, p = cmh_for_feature(X[:, j], y, strata)
            ors[j] = or_cmh
            pvals[j] = p

        qvals = bh_fdr(pvals)

        res = pd.DataFrame({
            "group_col": group_col,
            "group_level": str(g),
            "structure_col": structure_col,
            "feature": features,
            "or_cmh": ors,
            "p_value": pvals,
            "q_value_bh": qvals,
            "n_total": int(len(y)),
            "n_in_group": int(y.sum()),
        })
        results.append(res)

    out = pd.concat(results, axis=0, ignore_index=True)
    out.to_csv(out_csv, index=False)


# -----------------------------
# CLI
# -----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Host association + per-group accumulation curves for PIRATE presence/absence."
    )
    ap.add_argument("--pirate-out", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--meta-id-col", default="sample")
    ap.add_argument("--outdir", default="pirate_host_assoc")

    ap.add_argument("--group-col", default="host", help="e.g. host or country")
    ap.add_argument("--structure-col", default="CC", help="e.g. CC (default), ST, cgMLST cluster, PopPUNK cluster")

    ap.add_argument("--min-group-n", type=int, default=30)

    ap.add_argument("--min-prevalence", type=float, default=0.01)
    ap.add_argument("--max-prevalence", type=float, default=0.99)

    ap.add_argument("--rarefaction-per-group", action="store_true")
    ap.add_argument("--truncate-to-min", action="store_true")
    ap.add_argument("--rare-steps", type=int, default=60)
    ap.add_argument("--rare-reps", type=int, default=20)
    ap.add_argument("--rare-seed", type=int, default=0)

    ap.add_argument(
        "--assoc",
        action="store_true",
        help="Run CMH host association (group vs rest), stratified by structure",
    )

    ap.add_argument(
        "--group-colors",
        default=None,
        help='Optional: pin colours for specific group levels. Format: "LEVEL=#RRGGBB,LEVEL2=#RRGGBB"',
    )

    args = ap.parse_args()
    ensure_outdir(args.outdir)

    bin_fa = os.path.join(args.pirate_out, "binary_presence_absence.fasta")
    if not os.path.exists(bin_fa):
        raise FileNotFoundError(f"Missing: {bin_fa}")

    meta = load_metadata(args.meta, id_col=args.meta_id_col)

    # Optional colour map for group levels
    group_color_map = parse_group_colors(args.group_colors)

    # Load PIRATE matrix and align to metadata
    df01 = read_binary_presence_absence_fasta(bin_fa)
    df01 = df01.loc[df01.index.intersection(meta.index)]

    if df01.shape[0] < 10:
        raise ValueError(f"Too few samples after joining meta to PIRATE: n={df01.shape[0]}")

    if args.rarefaction_per_group:
        out_png = os.path.join(args.outdir, f"rarefaction_by_{args.group_col}.png")
        plot_rarefaction_per_group(
            df01, meta,
            group_col=args.group_col,
            out_png=out_png,
            reps=args.rare_reps,
            steps=args.rare_steps,
            seed=args.rare_seed,
            min_group_n=args.min_group_n,
            truncate_to_min=args.truncate_to_min,
            group_color_map=group_color_map,
        )

    if args.assoc:
        out_csv = os.path.join(args.outdir, f"cmh_assoc_{args.group_col}_strat_{args.structure_col}.csv")
        run_host_association(
            df01, meta,
            group_col=args.group_col,
            structure_col=args.structure_col,
            min_group_n=args.min_group_n,
            out_csv=out_csv,
            min_prev=args.min_prevalence,
            max_prev=args.max_prevalence,
        )

    print(f"Done. Outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
    
