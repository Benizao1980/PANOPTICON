#!/usr/bin/env python3
"""
Mask recombination from ClonalFrameML output and optionally draw an SVG.

Upgrades vs original:
- Robust ID normalization (strip/first token)
- Node lookup dict + graceful missing-node handling
- Optional ancestral masking
- Merge overlapping intervals per isolate before masking
- Stream FASTA output (memory efficient)
- Optional SVG without hard dependency unless requested
- TSV metrics export per isolate
- Alignment rectangularity check
"""

import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import csv
from collections import defaultdict, Counter

# -------------------------
# CLI
# -------------------------
parser = argparse.ArgumentParser(
    formatter_class=RawTextHelpFormatter,
    description="Mask recombination from ClonalFrameML output and (optionally) draw SVG of recombinant regions",
    usage="\n  %(prog)s --aln FASTA --out OUTPREFIX <CFML_PREFIX>",
)

parser.add_argument("prefix", metavar="CFML_PREFIX",
                    help="prefix used for CFML output files (required)")
parser.add_argument("--aln", metavar="FASTA", required=True,
                    help="multiFASTA alignment to mask (required; must be rectangular)")
parser.add_argument("--out", metavar="OUTPREFIX", default="maskrc",
                    help='output prefix for masked alignment (default="maskrc")')
parser.add_argument("--symbol", metavar="SYMBOL", default="N",
                    help='symbol to use for masking (default="N")')

parser.add_argument("--regions", metavar="FILE",
                    help="write recombinant regions table (TSV) to file")
parser.add_argument("--metrics", metavar="FILE",
                    help="write per-isolate masking metrics (TSV) to file")

parser.add_argument("--svg", metavar="FILE",
                    help="draw SVG output of recombinant regions and save as specified file")
parser.add_argument("--svgsize", metavar="WIDExHIGH", default="800x600",
                    help='specify width and height of SVG in pixels (default="800x600")')
parser.add_argument("--svgorder", metavar="FILE",
                    help="file containing list of taxa (1 per line) in desired order")
parser.add_argument("--svgcolour", metavar="COLOUR", default="black",
                    help='colour of extant recombination regions (default="black"; can be HEX)')
parser.add_argument("--consensus", action="store_true",
                    help="add consensus row of recombination hotspots to SVG")

parser.add_argument("--mask-ancestral", action="store_true",
                    help="also mask ancestral segments (default: extant only)")

parser.add_argument("--version", action="version", version="v0.4")
args = parser.parse_args()

# -------------------------
# Helpers
# -------------------------
def msg(*a, **k):
    print(*a, file=sys.stderr, **k)

def err(*a, **k):
    msg(*a, **k)
    sys.exit(1)

def check_file(path: str, label: str = "file"):
    if not os.path.isfile(path):
        err(f'ERROR: Cannot find {label} "{path}". Check CFML output files exist here or prefix is correct.')

def norm_id(x: str) -> str:
    """Normalize IDs so tree/FASTA/CFML match: strip + first whitespace-delimited token."""
    return x.strip().split()[0]

def merge_intervals(intervals):
    """Merge overlapping/adjacent intervals. Input: list of (start,end,type). Output merged per type."""
    by_type = defaultdict(list)
    for s, e, z in intervals:
        by_type[z].append((s, e))
    out = []
    for z, ints in by_type.items():
        ints = sorted(ints)
        merged = []
        for s, e in ints:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        out.extend([(s, e, z) for s, e in merged])
    return sorted(out, key=lambda t: (t[2], t[0], t[1]))

def wrap_write_fasta(out_handle, header: str, seq_bytes: bytes, width: int = 80):
    out_handle.write(f">{header}\n")
    for i in range(0, len(seq_bytes), width):
        out_handle.write(seq_bytes[i:i+width].decode() + "\n")

def parse_cfml_em(em_path: str):
    """
    Best-effort parse of global CFML parameters from .em.txt.
    CFML output formats vary; we extract float-like tokens keyed by common names.
    """
    params = {}
    if not os.path.isfile(em_path):
        return params
    txt = open(em_path, "r", encoding="utf-8", errors="replace").read().splitlines()
    # Heuristic: look for lines like "rho/theta = ..." or "r/m = ..."
    keys = [
        ("rho/theta", ["rho/theta", "rho/theta="]),
        ("r/m", ["r/m", "r/m="]),
        ("delta", ["delta", "delta="]),
        ("nu", ["nu", "nu="]),
        ("theta", ["theta", "theta="]),
        ("rho", ["rho", "rho="]),
    ]
    for line in txt:
        low = line.lower().replace(" ", "")
        for k, patterns in keys:
            if any(p in low for p in patterns):
                # grab first float-ish token in the original line
                tokens = line.replace("=", " ").replace("\t", " ").split()
                for tok in tokens[::-1]:
                    try:
                        val = float(tok)
                        params[k] = val
                        break
                    except ValueError:
                        continue
    return params

# -------------------------
# CFML input files
# -------------------------
cfmlPREFIX = str(args.prefix)
cfmlTREE = cfmlPREFIX + ".labelled_tree.newick"
cfmlRECOMB = cfmlPREFIX + ".importation_status.txt"
cfmlEM = cfmlPREFIX + ".em.txt"

check_file(cfmlTREE, "tree")
check_file(cfmlRECOMB, "importation status")

# -------------------------
# Imports (delay optional)
# -------------------------
try:
    from ete3 import Tree
except Exception as e:
    err(f"ERROR: Failed to import ete3. Install into this environment (conda install -c conda-forge ete3).\n{e}")

# Only import svgwrite if we need SVG
if args.svg:
    try:
        import svgwrite
    except Exception as e:
        err(f"ERROR: Failed to import svgwrite. Install into this environment (conda install -c conda-forge svgwrite).\n{e}")

# -------------------------
# Load CFML tree + build node lookup
# -------------------------
t = Tree(cfmlTREE, format=1)
node_by_name = {norm_id(n.name): n for n in t.traverse() if getattr(n, "name", None)}

# Determine leaf order for SVG and for metrics (optional)
leafLIST = []
if args.svgorder:
    with open(args.svgorder, "r", encoding="utf-8", errors="replace") as fh:
        leafLIST = [norm_id(x) for x in fh.read().splitlines() if x.strip()]
else:
    leafLIST = [norm_id(leaf.name) for leaf in t.iter_leaves()]

leafSET = set(leafLIST)

# -------------------------
# Read recombinant regions and assign to leaves
# d[leaf_id] = list of (start, stop, type)
# type in {"extant","ancestral"}
# -------------------------
d = defaultdict(list)

# The importation_status format is expected to have at least:
# node_id \t start \t stop ...
with open(cfmlRECOMB, "r", encoding="utf-8", errors="replace") as fh:
    rdr = csv.reader(fh, delimiter="\t")
    # Skip header line robustly
    header = next(rdr, None)

    missing_nodes = 0
    for row in rdr:
        if not row or len(row) < 3:
            continue
        seq = norm_id(row[0])
        try:
            RCstart = int(row[1])
            RCstop = int(row[2])
        except ValueError:
            continue

        if RCstop < RCstart:
            RCstart, RCstop = RCstop, RCstart
        if RCstop <= 0:
            continue

        node = node_by_name.get(seq)
        if node is None:
            missing_nodes += 1
            if missing_nodes <= 20:
                msg(f'WARNING: node "{seq}" not found in tree labels; skipping this record')
            continue

        if node.is_leaf():
            # extant segment
            d[norm_id(node.name)].append((RCstart, RCstop, "extant"))
        else:
            # ancestral segment applies to all descendant leaves
            if not args.mask_ancestral and not args.svg:
                # if we aren't masking ancestral and no SVG requested, we can skip distributing ancestral blocks
                # (but we still might want metrics; keep extant only)
                continue

            for leaf in node.get_leaves():
                lid = norm_id(leaf.name)
                d[lid].append((RCstart, RCstop, "ancestral"))

if missing_nodes > 20:
    msg(f"WARNING: {missing_nodes} CFML records had node IDs not found in tree labels (only first 20 shown).")

# Merge intervals per isolate to reduce masking operations
for k in list(d.keys()):
    d[k] = merge_intervals(d[k])

# -------------------------
# Mask alignment (streaming output)
# -------------------------
masked_fasta = args.out if args.out.endswith((".fa", ".fasta", ".fas", ".aln")) else (args.out + ".fasta")
symbol = args.symbol
if len(symbol) != 1:
    err("ERROR: --symbol must be a single character (e.g., N, ?, -)")

# We also collect metrics while streaming
metrics = {}
# We'll need alignment length; ensure rectangular FASTA
expected_len = None
n_records = 0

try:
    from Bio import SeqIO
except Exception as e:
    err(f"ERROR: Biopython not available in this environment. Install it (conda install -c conda-forge biopython).\n{e}")

msg(f'Writing masked alignment to "{masked_fasta}" ...')

with open(masked_fasta, "w", encoding="utf-8") as out:
    for record in SeqIO.parse(args.aln, "fasta"):
        n_records += 1
        rid = norm_id(record.id)
        seq_str = str(record.seq)
        seqlen = len(seq_str)

        if expected_len is None:
            expected_len = seqlen
        elif seqlen != expected_len:
            err(f'ERROR: Alignment is not rectangular: "{record.id}" has length {seqlen}, expected {expected_len}')

        # mask
        seq = bytearray(seq_str.encode())
        reg = d.get(rid, [])

        # metrics
        ext_bp = anc_bp = 0
        ext_n = anc_n = 0

        if reg:
            for (start, stop, z) in reg:
                # clamp to alignment bounds (1-based inclusive in CFML)
                s = max(1, start)
                e = min(expected_len, stop)
                if e <= s:
                    continue
                # decide whether to mask ancestral
                if z == "ancestral" and not args.mask_ancestral:
                    continue

                # convert to 0-based half-open
                s0 = s - 1
                e0 = e
                L = e0 - s0
                seq[s0:e0] = (symbol.encode() * L)

                if z == "extant":
                    ext_bp += L
                    ext_n += 1
                else:
                    anc_bp += L
                    anc_n += 1

        wrap_write_fasta(out, record.id, seq, width=80)

        metrics[rid] = {
            "seqlen": expected_len,
            "n_segments_extant": ext_n,
            "bp_recomb_extant": ext_bp,
            "n_segments_ancestral": anc_n,
            "bp_recomb_ancestral": anc_bp,
            "bp_recomb_total": ext_bp + anc_bp,
        }

msg(f"Done masking {n_records} sequences. Alignment length={expected_len} bp.")

# -------------------------
# Optional: write recombinant regions table (TSV)
# -------------------------
if args.regions:
    msg(f'Writing recombinant regions table to "{args.regions}" ...')
    with open(args.regions, "w", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["taxon", "start", "stop", "type"])
        for taxon, intervals in d.items():
            for s, e, z in intervals:
                w.writerow([taxon, s, e, z])

# -------------------------
# Optional: write per-isolate metrics (TSV)
# -------------------------
if args.metrics:
    msg(f'Writing per-isolate metrics to "{args.metrics}" ...')
    global_params = parse_cfml_em(cfmlEM)

    # Ensure we output all taxa (even those not in alignment) if they exist in leafLIST
    all_taxa = set(metrics.keys()) | leafSET

    with open(args.metrics, "w", encoding="utf-8") as fh:
        cols = [
            "taxon", "seqlen",
            "n_segments_extant", "bp_recomb_extant",
            "n_segments_ancestral", "bp_recomb_ancestral",
            "bp_recomb_total", "frac_recomb_total",
        ]
        # include global params as repeated columns (handy for joins)
        for k in ["rho/theta", "r/m", "delta", "nu", "theta", "rho"]:
            if k in global_params:
                cols.append(f"cfml_{k.replace('/','_')}")
        fh.write("\t".join(cols) + "\n")

        for taxon in sorted(all_taxa):
            m = metrics.get(taxon, {
                "seqlen": expected_len if expected_len is not None else "",
                "n_segments_extant": 0,
                "bp_recomb_extant": 0,
                "n_segments_ancestral": 0,
                "bp_recomb_ancestral": 0,
                "bp_recomb_total": 0,
            })
            seql = m["seqlen"]
            frac = (m["bp_recomb_total"] / float(seql)) if seql else 0.0
            row = [
                taxon, str(seql),
                str(m["n_segments_extant"]), str(m["bp_recomb_extant"]),
                str(m["n_segments_ancestral"]), str(m["bp_recomb_ancestral"]),
                str(m["bp_recomb_total"]), f"{frac:.6f}",
            ]
            for k in ["rho/theta", "r/m", "delta", "nu", "theta", "rho"]:
                if k in global_params:
                    row.append(str(global_params[k]))
            fh.write("\t".join(row) + "\n")

# -------------------------
# Optional: SVG drawing
# -------------------------
if args.svg:
    msg(f'Drawing SVG to "{args.svg}" ...')

    svgsize = args.svgsize.split("x", 1)
    width = int(svgsize[0])
    height = int(svgsize[1])

    numseqs = max(1, len(leafLIST))
    h = height / float(numseqs + (1 if args.consensus else 0))
    s = width / float(expected_len if expected_len else 1)
    fsize = max(6, int(h * 0.9))
    font = f"font-size:{fsize}px; font-family:Arial"
    colour = args.svgcolour
    main_colour = "black"
    ancestral_colour = "#A9A9A9"
    interval = 500000  # tick interval in bp

    dwg = svgwrite.Drawing(args.svg)

    def rect(x_bp, row_idx, w_px, fill):
        dwg.add(dwg.rect(insert=((x_bp * s) + 5, (row_idx * h) + 5),
                         size=(w_px, h * 0.8), fill=fill))
        if args.consensus:
            dwg.add(dwg.rect(insert=((x_bp * s) + 5, 5),
                             size=(w_px, h * 0.8), fill=main_colour))

    def label(n, row_idx):
        dwg.add(dwg.text(n, insert=((expected_len * s) + 15, (((row_idx + 1) * h) - (0.2 * h)) + 5),
                         fill=main_colour, style=font))

    def ticks(q):
        count = interval
        while expected_len and count < expected_len:
            dwg.add(dwg.line((count * s, (q * h)), (count * s, (q * h) - 3), stroke=main_colour))
            dwg.add(dwg.text((count / 1_000_000.0), insert=((count * s) - 6, (q * h) + 13),
                             fill=main_colour, style=font))
            count += interval

    def box(w, hh, p_last):
        dwg.add(dwg.line((0, 0), ((w + 10), 0), stroke=main_colour))
        dwg.add(dwg.line(((w + 10), 0), ((w + 10), ((p_last + 1) * hh)), stroke=main_colour))
        dwg.add(dwg.line(((w + 10), ((p_last + 1) * hh)), (0, ((p_last + 1) * hh)), stroke=main_colour))
        dwg.add(dwg.line((0, ((p_last + 1) * hh)), (0, 0), stroke=main_colour))

    # draw
    p = 1 if args.consensus else 0
    if args.consensus:
        dwg.add(dwg.text("Consensus", insert=((expected_len * s) + 15, (0.8 * h) + 5),
                         fill=main_colour, style=font))

    for seq in leafLIST:
        v = d.get(seq, [])
        for (x, y, z) in v:
            if z == "ancestral" and not args.mask_ancestral:
                continue
            # x,y are bp coords 1-based inclusive
            # convert to bp for rendering (just use x)
            w_bp = max(1, (y - x))
            w_px = w_bp * s
            if w_px < 1:
                w_px = 1
            if z == "extant":
                rect(x, p, w_px, colour)
            else:
                rect(x, p, w_px, ancestral_colour)

        label(seq, p)
        p += 1

    box(width, h, p)
    ticks(p + 1)
    dwg.save()

msg("Done.")
sys.exit(0)
