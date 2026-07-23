#!/usr/bin/env python3
"""Standardise host metadata and add broader biological groupings.

This script is intended for PANOPTICON metadata tables. It:

1. preserves the original host label in ``host_original``;
2. normalises case, whitespace, and punctuation in ``host``;
3. adds ``host_group`` for biologically similar host taxa;
4. adds ``ecology_group`` for broad ecological categories;
5. writes audit tables showing label counts and any unmapped values.

The mappings below are explicit and deliberately conservative. In particular,
pig and wild boar remain separate because that comparison is central to the
current S. aureus project.

Example
-------
python prepare_panopticon_metadata.py \
  --input metadata_final.txt \
  --output metadata_grouped.tsv \
  --host-col host
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Explicit host mappings
# ---------------------------------------------------------------------------

# Normalised host label -> (host_group, ecology_group)
HOST_MAP: Dict[str, Tuple[str, str]] = {
    # Humans
    "human": ("human", "human"),

    # Pigs: keep wild boar separate
    "pig": ("pig", "livestock"),
    "pig carcass": ("pig", "food"),
    "wild boar": ("wild_boar", "wildlife"),

    # Cattle and bovids
    "cow": ("cattle", "livestock"),
    "cow milk": ("cattle", "food"),
    "buffalo milk": ("other_bovid", "food"),
    "yak": ("other_bovid", "livestock"),

    # Small ruminants
    "sheep": ("small_ruminant", "livestock"),
    "sheep milk": ("small_ruminant", "food"),
    "goat": ("small_ruminant", "livestock"),
    "goat milk": ("small_ruminant", "food"),
    "wild goat milk": ("small_ruminant", "food"),
    "mountain goat": ("small_ruminant", "wildlife"),
    "european muflon": ("small_ruminant", "wildlife"),
    "iberian ibex": ("small_ruminant", "wildlife"),

    # Camelids
    "camel": ("camelid", "livestock"),
    "camel milk": ("camelid", "food"),

    # Equids
    "horse": ("equine", "livestock"),
    "donkey": ("equine", "livestock"),
    "mule": ("equine", "livestock"),

    # Poultry and domestic birds
    "chicken": ("poultry", "livestock"),
    "chicken carcass": ("poultry", "food"),
    "turkey": ("poultry", "livestock"),
    "duck": ("poultry", "livestock"),
    "goose": ("poultry", "livestock"),

    # Game and wild birds
    "wild turkey": ("game_bird", "wildlife"),
    "pheasant": ("game_bird", "wildlife"),
    "partridge": ("game_bird", "wildlife"),
    "grouse": ("game_bird", "wildlife"),
    "magpie": ("wild_bird", "wildlife"),
    "great tit": ("wild_bird", "wildlife"),
    "turkey vulture": ("wild_bird", "wildlife"),
    "parrot": ("wild_bird", "wildlife"),
    "lesser yellowlegs": ("wild_bird", "wildlife"),
    "screech owl": ("wild_bird", "wildlife"),
    "glaucous blue grosbeak": ("wild_bird", "wildlife"),
    "yellow eyed penguin": ("wild_bird", "wildlife"),
    "white bearded manakin": ("wild_bird", "wildlife"),

    # Companion animals
    "dog": ("companion_animal", "companion"),
    "cat": ("companion_animal", "companion"),
    "ferret": ("companion_animal", "companion"),

    # Rodents and related small mammals
    "mouse": ("rodent", "wildlife"),
    "yellow necked mouse": ("rodent", "wildlife"),
    "rat": ("rodent", "wildlife"),
    "common vole": ("rodent", "wildlife"),
    "bank vole": ("rodent", "wildlife"),
    "field vole": ("rodent", "wildlife"),
    "hamster": ("rodent", "companion"),
    "cloudrat": ("rodent", "wildlife"),
    "blank vole": ("rodent", "wildlife"),  # likely source typo for bank vole
    "laboratory mouse": ("rodent", "laboratory"),
    "cavy": ("rodent", "other"),

    # Squirrels
    "red squirrel": ("squirrel", "wildlife"),

    # Hedgehogs and shrews
    "hedgehog": ("hedgehog", "wildlife"),
    "european hedgehog": ("hedgehog", "wildlife"),
    "common shrew": ("shrew", "wildlife"),

    # Rabbits
    "rabbit": ("rabbit", "wildlife"),
    "european rabbit": ("rabbit", "wildlife"),

    # Cervids
    "red deer": ("cervid", "wildlife"),
    "roe deer": ("cervid", "wildlife"),
    "fallow deer": ("cervid", "wildlife"),

    # Bats
    "bat": ("bat", "wildlife"),
    "livingstone's fruit bat": ("bat", "wildlife"),

    # Primates
    "monkey": ("primate", "wildlife"),
    "macaque": ("primate", "wildlife"),
    "chimpanzee": ("primate", "wildlife"),
    "geoffroy's spider monkey": ("primate", "wildlife"),

    # Carnivores and mustelids
    "mink": ("mustelid", "wildlife"),
    "beech marten": ("mustelid", "wildlife"),
    "fox": ("wild_carnivore", "wildlife"),
    "mongoose": ("wild_carnivore", "wildlife"),
    "meerkat": ("wild_carnivore", "wildlife"),
    "common genet": ("wild_carnivore", "wildlife"),

    # Marine mammals
    "dolphin": ("marine_mammal", "wildlife"),
    "porpoise": ("marine_mammal", "wildlife"),
    "seal": ("marine_mammal", "wildlife"),
    "kangaroo": ("marsupial", "wildlife"),

    # Insects and invertebrates
    "housefly": ("insect", "environment"),
    "fly": ("insect", "environment"),
    "bee": ("insect", "environment"),
    "shrimp": ("aquatic_invertebrate", "food"),

    # Fish
    "fish": ("fish", "food"),
    "largehead hairtail": ("fish", "food"),
}


def normalise_label(value: object) -> str:
    """Normalise host labels while retaining biologically meaningful words."""
    if pd.isna(value):
        return "unknown"

    text = str(value).strip().lower()
    if not text:
        return "unknown"

    # Standardise separators and punctuation.
    text = text.replace("_", " ")
    text = re.sub(r"\s*,\s*", " ", text)
    text = re.sub(r"\s+", " ", text)

    # Normalise apostrophes and hyphens.
    text = text.replace("’", "'")
    text = text.replace("-", " ")

    # Remove punctuation except apostrophes.
    text = re.sub(r"[^\w\s']", "", text)
    text = re.sub(r"\s+", " ", text).strip()

    return text or "unknown"


def load_table(path: Path) -> pd.DataFrame:
    """Load CSV or tab-delimited metadata as strings."""
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, dtype=str)
    return pd.read_csv(path, sep="\t", dtype=str)


def prepare_metadata(
    data: pd.DataFrame,
    host_col: str,
    unknown_group: str,
) -> pd.DataFrame:
    """Add standardised host and ecological grouping columns."""
    if host_col not in data.columns:
        raise ValueError(
            f"Host column '{host_col}' was not found. "
            f"Available columns: {list(data.columns)}"
        )

    out = data.copy()
    out["host_original"] = out[host_col]
    out["host"] = out[host_col].map(normalise_label)

    mapped = out["host"].map(HOST_MAP)
    out["host_group"] = mapped.map(
        lambda value: value[0] if isinstance(value, tuple) else unknown_group
    )
    out["ecology_group"] = mapped.map(
        lambda value: value[1] if isinstance(value, tuple) else "unknown"
    )

    # Keep a clear flag for manual review.
    out["host_mapping_status"] = out["host"].map(
        lambda value: "mapped" if value in HOST_MAP else "unmapped"
    )

    return out


def write_audit_tables(data: pd.DataFrame, output: Path) -> None:
    """Write label and grouping summaries beside the prepared metadata."""
    prefix = output.with_suffix("")

    original_summary = (
        data.groupby(
            ["host_original", "host", "host_group", "ecology_group",
             "host_mapping_status"],
            dropna=False,
        )
        .size()
        .reset_index(name="n")
        .sort_values(["host_mapping_status", "n"], ascending=[True, False])
    )
    original_summary.to_csv(
        Path(str(prefix) + ".host_mapping_summary.tsv"),
        sep="\t",
        index=False,
    )

    host_group_summary = (
        data["host_group"]
        .value_counts(dropna=False)
        .rename_axis("host_group")
        .reset_index(name="n")
    )
    host_group_summary.to_csv(
        Path(str(prefix) + ".host_group_counts.tsv"),
        sep="\t",
        index=False,
    )

    ecology_summary = (
        data["ecology_group"]
        .value_counts(dropna=False)
        .rename_axis("ecology_group")
        .reset_index(name="n")
    )
    ecology_summary.to_csv(
        Path(str(prefix) + ".ecology_group_counts.tsv"),
        sep="\t",
        index=False,
    )

    unmapped = original_summary[
        original_summary["host_mapping_status"] == "unmapped"
    ].copy()
    unmapped.to_csv(
        Path(str(prefix) + ".unmapped_hosts.tsv"),
        sep="\t",
        index=False,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Normalise host labels and add host_group and ecology_group "
            "columns for PANOPTICON analyses."
        )
    )
    parser.add_argument("--input", required=True, help="Input metadata TSV/CSV.")
    parser.add_argument("--output", required=True, help="Output metadata TSV.")
    parser.add_argument(
        "--host-col",
        default="host",
        help="Input host column name [default: host].",
    )
    parser.add_argument(
        "--unknown-group",
        default="other",
        help="host_group assigned to unmapped labels [default: other].",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = load_table(input_path)
    prepared = prepare_metadata(
        data=data,
        host_col=args.host_col,
        unknown_group=args.unknown_group,
    )
    prepared.to_csv(output_path, sep="\t", index=False)
    write_audit_tables(prepared, output_path)

    print(f"Rows: {len(prepared)}")
    print(f"Unique normalised hosts: {prepared['host'].nunique()}")
    print("\nHost-group counts:")
    print(prepared["host_group"].value_counts().to_string())
    print("\nEcology-group counts:")
    print(prepared["ecology_group"].value_counts().to_string())

    unmapped = prepared.loc[
        prepared["host_mapping_status"] == "unmapped", "host"
    ].value_counts()

    if len(unmapped):
        print("\nWARNING: Unmapped host labels remain:", file=sys.stderr)
        print(unmapped.to_string(), file=sys.stderr)
        print(
            "\nReview the .unmapped_hosts.tsv audit file before analysis.",
            file=sys.stderr,
        )
    else:
        print("\nAll host labels were mapped.")

    print(f"\nPrepared metadata written to: {output_path}")


if __name__ == "__main__":
    main()
