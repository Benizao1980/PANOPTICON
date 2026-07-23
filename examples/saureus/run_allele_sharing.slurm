#!/bin/bash
#SBATCH --job-name=saureus_alleles
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=saureus_alleles_%j.out
#SBATCH --error=saureus_alleles_%j.err

source ~/.bashrc
conda activate panopticon

ROOT=/groups/cooperma/bpascoe

python "$ROOT/PANOPTICON/scripts/pangenome_ecology.py" \
  --pangenome-tool panaroo \
  --pangenome-out "$ROOT/Saureus" \
  --meta "$ROOT/Saureus/metadata_grouped.tsv" \
  --meta-id-col id \
  --group-col ecology_group \
  --outdir "$ROOT/Saureus/ecology_allele_sharing" \
  --min-group-n 20 \
  --allele-seq-dir "$ROOT/Saureus/aligned_gene_sequences" \
  --allele-sequence-type nucleotide \
  --heatmaps
