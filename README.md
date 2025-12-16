# Pangenome analysis of the Zang et al. dataset

## Overview

To characterise gene content diversity, structural variation, and host-associated genomic features across the *Zang et al.* genome collection, we performed a comprehensive pangenome analysis using **PIRATE** (Pangenome Iterative Refinement and Threshold Evaluation). PIRATE is specifically designed for bacterial pangenome analysis across multiple sequence identity thresholds and explicitly accounts for allelic diversity, gene fragmentation, fission/fusion events, and gene duplication.

All analyses described below were performed on the full dataset and form the basis for downstream comparative and host-association analyses presented in this study.

---

## Genome dataset

The pangenome analysis comprised **2,327 whole-genome assemblies** derived from the *Zang et al.* dataset. Assemblies were processed in a consistent manner prior to pangenome reconstruction to minimise technical artefacts in gene clustering and annotation.

---

## Pangenome reconstruction

Pangenome reconstruction was performed using **PIRATE**, which clusters homologous genes iteratively across a descending series of amino-acid identity thresholds. This approach enables robust delineation of orthologous gene families while capturing biologically meaningful variation arising from sequence divergence, gene fragmentation, and paralogy.

PIRATE was executed on the complete genome set using default parameters unless otherwise specified:

```slurm
#!/bin/bash
#SBATCH --job-name=pirate_zang
#SBATCH --output=slurm_logs/pirate_%A.out
#SBATCH --error=slurm_logs/pirate_%A.err
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=72:00:00
#SBATCH --mem=250G

# --------------------
# Environment
# --------------------
source ~/.bashrc
conda activate pirate

# --------------------
# Directories
# --------------------
GFF_DIR="/groups/cooperma/bpascoe/Zang_cdtB/pirate_gff"
OUTDIR="/groups/cooperma/bpascoe/Zang_cdtB/PIRATE_out"

mkdir -p "$OUTDIR"

# --------------------
# PIRATE Command (Campylobacter thresholds)
# --------------------
PIRATE \
    -i "$GFF_DIR" \
    -o "$OUTDIR" \
    -t 16 \
    -a \
    -s "80,85,90,92,94,95,96,97,98"
```

During execution, PIRATE:

* Identifies homologous gene clusters across genomes
* Resolves allelic variants within clusters
* Detects gene fission and fusion events arising from fragmentation or recombination
* Identifies gene duplication and copy-number variation

---

## Global pangenome properties

The resulting pangenome comprised **5,387 gene families** across the 2,327 genomes analysed. Of these:

* **682 gene families** contained more than one allele at the identity thresholds analysed, indicating substantial allelic diversity.
* **1,701 gene families** showed evidence of gene fission or fusion events.
* **998 gene families** exhibited gene duplication or loss, reflecting copy-number variation across isolates.

These features indicate a highly dynamic genome architecture with extensive structural and allelic variability beyond simple presence–absence differences.

---

## Gene frequency and structural variation analysis

Gene families were stratified according to their frequency across isolates to assess how allelic diversity and structural variation scale with prevalence in the population.

| % isolates | # gene clusters | >1 allele | fission/fusion | multicopy |
| ---------: | --------------: | --------: | -------------: | --------: |
|      0–10% |           3,518 |       249 |            325 |        83 |
|     10–25% |             228 |       118 |            143 |       105 |
|     25–50% |             128 |        60 |             92 |        68 |
|     50–75% |              68 |        29 |             55 |        41 |
|     75–90% |              50 |        21 |             41 |        35 |
|     90–95% |              29 |         6 |             25 |        17 |
|    95–100% |           1,366 |       199 |          1,020 |       649 |

Rare gene families (≤10% prevalence) accounted for the majority of the pangenome, consistent with a large accessory gene pool. Notably, even highly prevalent genes frequently exhibited fission/fusion and multicopy signals, indicating ongoing structural plasticity within the core genome.

---

## Visualisation of pangenome structure

To summarise pangenome composition and structural variation, we generated a series of visualisations from the PIRATE outputs using a custom plotting workflow:

```bash
python3 pirate_plotting.py \
  --pirate-out PIRATE_out \
  --outdir pirate_plots \
  --all \
  --palette life_aquatic_blues
```

These plots provide an overview of core and accessory genome structure, allelic diversity, and the distribution of structural variation across gene families.

---

## Host-associated gene analysis

To investigate host-specific patterns of gene presence and absence, we conducted host-association analyses using PIRATE’s association testing framework. Genome metadata describing host origin were supplied, and both rarefaction and association testing were enabled to account for unequal sampling across host groups.

```bash
python3 pirate_host_association.py \
  --pirate-out PIRATE_out \
  --meta meta.tsv \
  --outdir pirate_host_assoc \
  --group-col host \
  --rarefaction-per-group \
  --assoc \
  --min-group-n 1 \
  --min-prevalence 0.01 \
  --max-prevalence 0.99
```

This analysis framework enables identification of accessory genes whose distribution is significantly associated with host category while controlling for gene prevalence and sampling depth.

---

## Summary

Together, these analyses demonstrate that the *Zang et al.* dataset harbours a large and highly variable pangenome characterised by extensive accessory gene diversity, widespread structural variation, and detectable host-associated gene content patterns. This pangenome framework provides the foundation for subsequent functional and evolutionary analyses presented in this study.
