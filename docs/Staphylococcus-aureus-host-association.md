# Animal-associated *Staphylococcus aureus*: reproducible PANOPTICON workflow

This tutorial documents the workflow used for the animal-associated *S. aureus* collection. It is written so that the analysis can be repeated from annotated assemblies through pangenome exploration, ecology sharing, allele sharing, phylogenetics, and preparation for GWAS.

## 1. Questions

The main biological questions are:

1. Is genomic sharing greater between physiologically similar hosts, such as pigs and wild boar, or among hosts linked by agricultural ecology, such as pigs, cattle, and poultry?
2. Does observed gene-family sharing remain after equalising sample size?
3. Do ecologies carrying the same gene families also carry identical alleles?
4. Which associations persist after controlling for clonal population structure?

Do not infer recent HGT from gene-family overlap alone. Population structure, shared ancestry, geography, and study composition must be evaluated.

## 2. Project layout

```text
Saureus/
├── assemblies/
├── gff/
├── panaroo_output/
├── metadata/
├── phylogeny/
├── panopticon_qc/
├── ecology_stratified/
├── ecology_gene_sharing/
├── ecology_allele_sharing/
└── slurm_jobs/
```

## 3. Annotation

Use one annotation program and version for all assemblies. Bakta is preferred; Prokka is acceptable.

```bash
bakta \
  --db /path/to/bakta_db \
  --output annotations/sample_001 \
  --prefix sample_001 \
  --threads 8 \
  assemblies/sample_001.fasta
```

Collect the GFF3 files and verify that every expected assembly has one annotation.

```bash
find annotations -name '*.gff3' | wc -l
```

Record the Bakta/Prokka version and database version.

## 4. Construct the Panaroo pangenome

```bash
panaroo \
  -i gff/*.gff \
  -o panaroo_output \
  --clean-mode strict \
  --remove-invalid-genes \
  -t 16
```

For several thousand genomes, run Panaroo as a scheduled HPC job with generous memory and walltime. Do not delete partial output until the run is confirmed complete.

### Required files for PANOPTICON

```text
gene_presence_absence.Rtab
gene_presence_absence.csv
```

`gene_presence_absence.Rtab` is the preferred binary matrix. `gene_presence_absence.csv` supplies gene names, descriptions, and locus tags. `summary_statistics.txt` is useful but optional.

Check that files are non-empty:

```bash
ls -lh panaroo_output/gene_presence_absence.Rtab \
       panaroo_output/gene_presence_absence.csv
```

PANOPTICON can continue from an Rtab alone, but annotation columns will be blank.

## 5. Core-genome alignment

Generate the alignment separately if the original Panaroo run omitted it:

```bash
panaroo-msa \
  -a core \
  --core_threshold 0.98 \
  -o panaroo_output \
  --aligner mafft \
  -t 16
```

Retain:

```text
core_gene_alignment.aln
core_gene_alignment_filtered.aln
```

Use the filtered alignment for the primary maximum-likelihood tree where appropriate. Keep the unfiltered alignment for sensitivity analyses and troubleshooting.

## 6. Install PANOPTICON

```bash
git clone https://github.com/Benizao1980/PANOPTICON.git
cd PANOPTICON
mamba env create -f environment.yml
conda activate panopticon
```

Check the scripts:

```bash
python -m py_compile scripts/*.py
```

## 7. Prepare metadata

The minimum metadata table contains one row per genome and an ID matching the Panaroo column name exactly.

```text
id  host  host_group  ecology_group  country  collection_date  ST  CC  study
```

Retain the original host label. Add grouped variables rather than overwriting it.

Recommended distinctions:

- pig and wild boar remain separate
- cattle-derived milk/meat can be retained as `food` rather than silently merged with living cattle
- poultry and wild/game birds remain separate where possible
- broad ecology categories may include `livestock`, `wildlife`, `companion`, `food`, and `environment`

Prepare and audit labels:

```bash
python scripts/prepare_panopticon_metadata.py \
  --input metadata_final.txt \
  --output metadata_grouped.tsv \
  --host-col host
```

Review:

```text
metadata_grouped.host_mapping_summary.tsv
metadata_grouped.host_group_counts.tsv
metadata_grouped.ecology_group_counts.tsv
metadata_grouped.unmapped_hosts.tsv
```

Do not proceed until unmapped labels and unexpected tiny categories have been reviewed.

## 8. Validate the pangenome

```bash
python scripts/pangenome_io.py \
  --pangenome-tool panaroo \
  --pangenome-out panaroo_output
```

For the current collection, the expected values are approximately:

```text
samples: 8574
gene families: 15655
core (>=99%): 1250
soft-core (95-99%): 596
shell (15-95%): 1252
cloud (<15%): 12557
```

A discrepancy usually indicates a different dataset, incomplete copy, or parsing problem.

## 9. Initial pangenome exploration

```bash
python scripts/pangenome_plotting.py \
  --pangenome-tool panaroo \
  --pangenome-out panaroo_output \
  --outdir panopticon_qc \
  --pangenome-bars \
  --gene-freq \
  --gene-counts \
  --pca
```

Inspect:

- pangenome category totals sum to all families
- per-genome family counts are biologically plausible
- PCA clusters are not driven by obvious low-quality assemblies
- colour PCA by CC, ST, host, ecology, country, study, and gene-family count

## 10. Ecology-stratified accumulation

```bash
python scripts/pangenome_stratification.py \
  --pangenome-tool panaroo \
  --pangenome-out panaroo_output \
  --meta metadata_grouped.tsv \
  --meta-id-col id \
  --group-col ecology_group \
  --outdir ecology_stratified \
  --rarefaction-per-group
```

The raw endpoints reflect different sample sizes. Use the curves to inspect saturation and sampling depth; do not interpret the largest observed pangenome as intrinsically most diverse without rarefaction.

## 11. Gene-family sharing

Use a fixed rarefaction size that is below the smallest retained ecology group. In the current analysis, `n=40` allows the environment group to be included.

```bash
python scripts/pangenome_ecology.py \
  --pangenome-tool panaroo \
  --pangenome-out panaroo_output \
  --meta metadata_grouped.tsv \
  --meta-id-col id \
  --group-col ecology_group \
  --outdir ecology_gene_sharing \
  --min-group-n 20 \
  --min-group-family-count 2 \
  --min-group-family-prevalence 0.01 \
  --rarefy \
  --rarefy-n 40 \
  --rarefaction-reps 100 \
  --group-core-threshold 0.99 \
  --all-plots \
  --summary-figure \
  --stratification-plot ecology_stratified/rarefaction_by_ecology_group.png
```

### Main outputs

```text
tables/pairwise_gene_family_sharing_observed.tsv
tables/pairwise_gene_family_sharing_rarefied_summary.tsv
tables/gene_family_sharing_by_group.tsv
tables/private_gene_family_counts.tsv
tables/pairwise_group_core_sharing.tsv
plots/heatmap_jaccard_rarefied.png
plots/heatmap_group_core_jaccard.png
plots/private_gene_family_counts.png
plots/upset_gene_family_sharing.png
plots/network_jaccard.png
ecology_summary_figure.png
```

Prefer the rarefied Jaccard/Dice results for biological comparison. Raw shared counts are descriptive QC.

## 12. Allele sharing

The concatenated core alignment is not sufficient for assigning alleles to individual families. Supply either:

1. family-specific FASTA/alignment files; or
2. a table containing `sample`, `gene_family`, and `allele`.

With family alignments:

```bash
python scripts/pangenome_ecology.py \
  --pangenome-tool panaroo \
  --pangenome-out panaroo_output \
  --meta metadata_grouped.tsv \
  --meta-id-col id \
  --group-col ecology_group \
  --outdir ecology_allele_sharing \
  --min-group-n 20 \
  --allele-seq-dir aligned_gene_sequences \
  --allele-sequence-type nucleotide \
  --heatmaps
```

The current implementation defines an allele as an exact ungapped sequence. This is reproducible but strict. Report this definition explicitly. A later sensitivity analysis can cluster alleles by SNP distance or sequence identity.

Interpret:

- shared family + shared allele: compatible with recent sharing, common ancestry, or conserved sequence
- shared family + no shared exact allele: common function but structured sequence variation
- ecology-specific allele in a widespread family: potential host/ecology-associated diversification

## 13. Phylogeny

```bash
iqtree2 \
  -s panaroo_output/core_gene_alignment_filtered.aln \
  --seqtype DNA \
  -m MFP \
  -B 1000 \
  -T AUTO \
  --prefix phylogeny/saureus_core
```

Overlay host, ecology, country, study, ST, and CC. Quantify whether ecology groups are concentrated in particular lineages.

Gubbins or ClonalFrameML can be considered, but validate runtime and biological assumptions on a collection this large. A representative or lineage-specific approach may be more tractable than one global recombination analysis.

## 14. Is GWAS next?

Yes, after the structure audit and allele-sharing analysis.

Start with hypothesis-led binary comparisons rather than a single five-class test:

- pig versus wild boar
- pig versus cattle
- pig versus poultry
- livestock versus wildlife
- livestock versus companion

Control for:

- core-genome relatedness/kinship
- major CC or lineage
- country
- study/source dataset
- sampling imbalance

Use discovery/validation where possible: public genomes for discovery and the newly sequenced Polish isolates for validation, or vice versa depending on group sizes.

Top hits should be evaluated for:

- effect size
- lineage breadth
- within-CC replication
- annotation confidence
- mobile-element context
- consistency between gene-family and unitig analyses

## 15. Resource requirements

A default interactive allocation of 4 GB was insufficient. For 8,574 genomes × 15,655 families, use scheduled jobs.

Recommended starting allocations:

- validation and simple plots: 8–16 GB
- ecology sharing with 100 rarefaction replicates: 32 GB
- allele processing: 32 GB or more depending on alignment count and size

After a run:

```bash
sacct -j JOBID --format=JobID,State,Elapsed,MaxRSS,ReqMem,AllocCPUS
```

## 16. Troubleshooting

### Process prints only `Killed`

Check:

```bash
dmesg | tail -20
```

A cgroup memory message means Slurm killed the process for exceeding requested RAM.

### CSV exists but cannot be read

Check file size:

```bash
ls -lh gene_presence_absence.csv
```

A zero-byte CSV is treated as unavailable by the robust loader. Matrix analyses continue from the Rtab, but annotations are absent.

### Very few metadata matches

Compare metadata `id` values with Panaroo sample columns. Do not use display names or accession aliases unless they match exactly.

### UpSet or network plot missing

Install optional dependencies:

```bash
conda install -c conda-forge upsetplot networkx
```

### Summary figure missing

Install Pillow and provide the stratification PNG:

```bash
conda install -c conda-forge pillow
```

## 17. Paper-level interpretation

The present outputs suggest extensive sharing among livestock, wildlife, companion, and food-associated populations, with environment-associated isolates more distinct. These patterns should be described as repertoire similarity, not direct evidence of HGT. The central test is whether gene and allele associations remain after accounting for shared lineages and sampling provenance.
