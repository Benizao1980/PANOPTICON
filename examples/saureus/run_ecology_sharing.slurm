#!/bin/bash
#SBATCH --job-name=saureus_ecology
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=saureus_ecology_%j.out
#SBATCH --error=saureus_ecology_%j.err

source ~/.bashrc
conda activate panopticon

ROOT=/groups/cooperma/bpascoe

python "$ROOT/PANOPTICON/scripts/pangenome_ecology.py" \
  --pangenome-tool panaroo \
  --pangenome-out "$ROOT/Saureus" \
  --meta "$ROOT/Saureus/metadata_grouped.tsv" \
  --meta-id-col id \
  --group-col ecology_group \
  --outdir "$ROOT/Saureus/ecology_gene_sharing" \
  --min-group-n 20 \
  --min-group-family-count 2 \
  --min-group-family-prevalence 0.01 \
  --rarefy \
  --rarefy-n 40 \
  --rarefaction-reps 100 \
  --group-core-threshold 0.99 \
  --all-plots \
  --summary-figure \
  --stratification-plot "$ROOT/Saureus/ecology_stratified/rarefaction_by_ecology_group.png"
