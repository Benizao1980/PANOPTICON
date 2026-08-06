#!/usr/bin/env python3
"""Regression test for PANOPTICON PIRATE input parsing."""

from pathlib import Path
import tempfile
import pandas as pd

from pangenome_io import load_pangenome, read_pirate_binary_fasta


def main():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)

        # Family table order: famA, famB, famC.
        table = pd.DataFrame({
            "allele_name": ["a1", "b1", "c1"],
            "gene_family": ["famA", "famB", "famC"],
            "consensus_gene_name": ["A", "B", "C"],
            "consensus_product": ["prodA", "prodB", "prodC"],
            "threshold": ["95", "95", "95"],
            "number_genomes": ["2", "1", "2"],
            "cluster_order": ["1", "2", "3"],
            "s1": ["l1", "", "l3"],
            "s2": ["l1b", "l2", ""],
            "s3": ["", "", "l3b"],
        })
        table.to_csv(root / "PIRATE.gene_families.ordered.tsv", sep="\t", index=False)

        # Binary FASTA positions deliberately use a DIFFERENT family order:
        # [famC, famA, famB]. PIRATE A=present, C=absent.
        # s1: C present, A present, B absent -> A A C
        # s2: C absent,  A present, B present -> C A A
        # s3: C present, A absent,  B absent -> A C C
        (root / "binary_presence_absence.fasta").write_text(
            ">s1\nAAC\n>s2\nCAA\n>s3\nACC\n"
        )

        pg = load_pangenome("pirate", root)

        assert pg.presence_absence.shape == (3, 3)
        assert list(pg.presence_absence.columns) == ["famA", "famB", "famC"]
        assert list(pg.presence_absence.index) == ["s1", "s2", "s3"]

        expected = pd.DataFrame(
            [[1, 0, 1], [1, 1, 0], [0, 0, 1]],
            index=["s1", "s2", "s3"],
            columns=["famA", "famB", "famC"],
            dtype="uint8",
        )
        expected.index.name = "sample"
        expected.columns.name = "gene_family"
        pd.testing.assert_frame_equal(pg.presence_absence, expected)

        # Standalone binary fallback must also decode A as presence.
        binary = read_pirate_binary_fasta(root / "binary_presence_absence.fasta")
        assert binary.iloc[0].tolist() == [1, 1, 0]

        print("PASS: PIRATE ordered-table matrix and A/C decoding are correct.")


if __name__ == "__main__":
    main()
