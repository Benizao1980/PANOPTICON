#!/bin/bash
#SBATCH --job-name=saureus_strat
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=saureus_strat_%j.out
#SBATCH --error=saureus_strat_%j.err

source ~/.bashrc
conda activate panopticon

ROOT=/groups/cooperma/bpascoe

python "$ROOT/PANOPTICON/scripts/pangenome_stratification.py" \
  --pangenome-tool panaroo \
  --pangenome-out "$ROOT/Saureus" \
  --meta "$ROOT/Saureus/metadata_grouped.tsv" \
  --meta-id-col id \
  --group-col ecology_group \
  --outdir "$ROOT/Saureus/ecology_stratified" \
  --rarefaction-per-group
