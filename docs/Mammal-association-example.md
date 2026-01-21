## Overview

To characterise gene content diversity, structural variation, and host-associated genomic features across the *Zang et al.* genome collection, we performed a comprehensive pangenome analysis using **PIRATE** (Pangenome Iterative Refinement and Threshold Evaluation). PIRATE is designed for large bacterial genome collections and explicitly models allelic diversity, gene fragmentation, fission/fusion events, and gene duplication by clustering genes across multiple amino-acid identity thresholds.

This document combines:
- a **step-by-step computational walkthrough**,
- **inline commentary on observed results**, and
- **reproducibility metadata**, including software versions and references.

All analyses were performed on the full dataset and form the basis for downstream comparative, evolutionary, and host-association analyses, including recombination-aware phylogenetics and genome-wide association studies.

---

## Genome dataset

The analysis comprised **2,327 whole-genome assemblies** derived from the *Zang et al.* dataset. These genomes represent isolates sampled across multiple host species and epidemiological contexts. All assemblies were processed using a uniform annotation and pangenome workflow to minimise technical artefacts.

contigs are shared online at FigShare (doi: add) and can be downloaded: 

```bash
mkdir -p contigs
cd contigs
wget -O zang_contigs.tar.gz "TODO:FIGSHARE_DIRECT_DOWNLOAD_URL"
tar -xzf zang_contigs.tar.gz
```

Genomes and isolate data are also grouped as a collection on [pubMLST](add link)

*Commentary:*  
The large sample size and host diversity provide sufficient power to characterise both rare accessory genes and structural variation within highly prevalent (core) gene families.

---

## Genome annotation with Prokka

Consistent gene annotation is critical for pangenome analysis. Differences in gene calling or annotation standards can artificially inflate accessory gene counts or introduce spurious gene fragmentation. All genomes were therefore annotated de novo using **Prokka**, ensuring consistent gene prediction and functional annotation across the dataset.

### Input preparation

```bash
mkdir -p input
```

```bash
find contigs -maxdepth 1 -type f \
  \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fas" \) \
  -printf "%f\n" | sort > input/InputFiles
```

```bash
wc -l input/InputFiles
head input/InputFiles
```

### High-throughput annotation

Prokka was executed as a SLURM array job, processing one genome per task (up to 500 concurrent jobs):

```slurm
#!/bin/bash
#SBATCH --job-name=prokka_batch
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --array=1-500

source ~/.bashrc
conda activate prokka

WORKDIR="/home/u12/bpascoe/Zang_cdtB"
cd "$WORKDIR"

INPUT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input/InputFiles)
[ -z "$INPUT_FILE" ] && exit 0

PREFIX="${INPUT_FILE%.*}"
OUTDIR="output/${PREFIX}"
mkdir -p "$OUTDIR"

prokka --outdir "$OUTDIR" --prefix "$PREFIX" --cpus 4 --force "contigs/${INPUT_FILE}"
```

*Commentary:*  
This produced one annotated GFF per genome. These GFFs form the direct input for PIRATE and preserve consistent gene boundaries across the dataset.

---

## Pangenome reconstruction with PIRATE

Unlike single-threshold clustering methods, PIRATE iteratively clusters genes across descending amino-acid identity thresholds. This allows separation of true orthologues from divergent alleles and paralogues while retaining information on structural variation.

### Run PIRATE

First, collect GFFs into a single input directory:

```
mkdir -p pirate_gff
find output -name "*.gff" -type f -print -exec ln -sf {} pirate_gff/ \;
```

Then run PIRATE:

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

source ~/.bashrc
conda activate pirate

PIRATE -i pirate_gff -o PIRATE_out -t 16 -a -s "80,85,90,92,94,95,96,97,98"
```

*Commentary:*  
Identity thresholds spanning 80–98% are appropriate for *Campylobacter*, capturing both deep divergence and fine-scale allelic variation.

---

## Visualisation

Additional plots summarising core/accessory structure, allelic diversity, and structural variation:

```bash
python3 scripts/pirate_plotting.py \
  --pirate-out PIRATE_out \
  --outdir pirate_plots \
  --all \
  --palette life_aquatic_blues
```

---

## Global pangenome structure

The overall pangenome structure was summarised using rarefaction and gene-frequency analyses derived from PIRATE outputs.

![Pangenome rarefaction](figures/pangenome_rarefaction.png)

**Figure 1. Pangenome rarefaction curve.**
*Total pangenome size (upper curve) and core genome size (lower curve) as a function of the number of genomes sampled. Shaded regions indicate variability across random genome orderings. The continued increase in total gene families indicates an open pangenome, while the gradual decline in core genes reflects increasing genetic diversity.*

![Pangenome composition](figures/pangenome_composition.png)

**Figure 2. Pangenome composition across prevalence categories.**
*Gene families were classified into core (≥99%), soft-core (95–99%), shell (15–95%), and cloud (<15%) categories. The pangenome is dominated by low-frequency accessory genes, consistent with extensive horizontal gene transfer and lineage-specific gene content.*

The inferred pangenome comprised **5,387 gene families** across 2,327 genomes.
- **682 gene families** exhibited allelic diversity (>1 allele)
- **1,701 gene families** showed evidence of gene fission or fusion
- **998 gene families** exhibited gene duplication or loss

*Interpretation:*  
These results indicate a highly dynamic genome architecture. Structural variation is not restricted to rare accessory genes but is also prevalent among widely distributed gene families.

---

## Gene frequency and structural variation

Gene families were stratified by prevalence to assess how structural variation scales with frequency.

| % isolates | # gene clusters | >1 allele | fission/fusion | multicopy |
|-----------:|----------------:|----------:|---------------:|----------:|
| 0–10%      | 3,518           | 249       | 325            | 83        |
| 10–25%     | 228             | 118       | 143            | 105       |
| 25–50%     | 128             | 60        | 92             | 68        |
| 50–75%     | 68              | 29        | 55             | 41        |
| 75–90%     | 50              | 21        | 41             | 35        |
| 90–95%     | 29              | 6         | 25             | 17        |
| 95–100%    | 1,366           | 199       | 1,020          | 649       |

*Interpretation:*  
Rare gene families dominate the accessory genome, consistent with an open pangenome. However, even near-core genes frequently show fission/fusion and copy-number variation, indicating ongoing structural evolution within the core genome.

---

### Accessory genome structure

The structure of the accessory genome was explored using principal component analysis (PCA) based on presence–absence of accessory gene families.

![Accessory genome PCA](figures/accessory_pca.png)

**Figure 3. Accessory genome PCA.**
*Principal component analysis of accessory gene presence–absence reveals clear clustering, indicating non-random structure in accessory gene content consistent with lineage- and host-associated gene pools.*

---

## Host-stratified pangenome accumulation analysis

To investigate whether pangenome structure differs systematically by host, we extended the PIRATE-based analysis to incorporate host metadata and performed **host-stratified gene family accumulation (rarefaction) analyses**. This allows direct comparison of pangenome openness and core genome stability across host-associated populations while controlling for unequal sampling depth.

Genomes were grouped by host category (human, mammal, bird) using curated metadata, and gene family accumulation curves were estimated independently for each group.
For each host category, genomes were randomly subsampled without replacement across increasing sample sizes. At each subsampling depth, we calculated:

- the total number of gene families observed (**pangenome size**)
- the number of gene families present in all sampled genomes (**core genome size**)

This procedure was repeated across multiple random permutations to estimate mean trends and variability.
The analysis was implemented using the PIRATE-derived gene presence–absence matrix and host metadata, and integrated into the existing host-association workflow.

```bash
python3 pirate_host_association.py \
  --pirate-out PIRATE_out \
  --meta meta.tsv \
  --outdir pirate_host_assoc \
  --group-col host \
  --rarefaction-per-group \
  --assoc \
  --min-prevalence 0.01 \
  --max-prevalence 0.99
```

*Commentary:*
Performing rarefaction independently within each host category controls for differences in sample size and avoids conflating biological signal with uneven sampling.

Host-stratified accumulation curves reveal clear differences in pangenome structure between host-associated populations.

![Host-stratified gene family accumulation](figures/rarefaction_by_host.png)

**Figure 4. Host-stratified gene family accumulation curves.**
*Pangenome size (solid lines) and core genome size (dashed lines) plotted as a function of the number of genomes sampled within each host category. Shaded regions indicate variability across permutations.*

### Key observations:
- All host groups exhibit open pangenomes, with the number of gene families continuing to increase as additional genomes are sampled.
- Human-associated isolates show the largest and most rapidly expanding pangenome, even at large sample sizes (>1,000 genomes), with no evidence of saturation.
- Bird-associated isolates show the most rapid core genome erosion, with the number of core gene families declining steeply as additional genomes are included.
- Mammal-associated isolates display intermediate behaviour, with both pangenome growth and core genome contraction falling between human- and bird-associated populations.

*Interpretation:*
These host-specific accumulation patterns suggest distinct evolutionary regimes. The expansive pangenome in human-associated isolates likely reflects ecological heterogeneity and repeated host switching, while the rapid core erosion in birds is consistent with elevated recombination and population structure in avian hosts.

---

## Core-genome phylogeny with IQ-TREE

To infer a high-quality maximum-likelihood (ML) phylogeny from the PIRATE core alignment, we used IQ-TREE2 on the nucleotide core alignment (PIRATE_out/core_alignment.fasta) under a GTR+F+I+G4 model, with branch support from both ultrafast bootstrap and SH-aLRT. The resulting tree provides the fixed topology required by ClonalFrameML.

```bash
iqtree2 \
  -s PIRATE_out/core_alignment.fasta \
  -m GTR+F+I+G4 \
  -T 40 \
  -B 1000 \
  --alrt 1000 \
  --prefix zang_core \
  --seed 12345 \
  --safe
```

*Commentary:*
This run used IQ-TREE2 v2.3.6, inferred a tree from 2,327 sequences with an alignment length of 1,371,024 bp, and used a fixed seed for reproducibility. 

### Key outputs

- zang_core.treefile — **ML tree** (Newick; used for ClonalFrameML)
- zang_core.iqtree / zang_core.log — **run record and model details**
- zang_core.ckp.gz — **checkpoint file for safe/restart mode**
- zang_core.mldist — **pairwise ML distances** (optional downstream use)

## Recombination-aware inference with ClonalFrameML (CF-ML)

We next used ClonalFrameML to infer recombinant tracts on the fixed IQ-TREE topology, enabling downstream analyses on a recombination-masked alignment (or recombination-aware branch lengths, depending on the application).

### Inputs

- Tree: zang_core.treefile (IQ-TREE output)
- Alignment: PIRATE_out/core_alignment.fasta (core alignment)

```bash
ClonalFrameML zang_core.treefile PIRATE_out/core_alignment.fasta cfml_core
```

Example SLURM job
```slurm
#!/bin/bash
#SBATCH --job-name=cfml_core
#SBATCH --output=slurm_logs/cfml_%A.out
#SBATCH --error=slurm_logs/cfml_%A.err
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
#SBATCH --mem=160G

set -euo pipefail
source ~/.bashrc
conda activate clonalframeml

cd /groups/cooperma/bpascoe/Zang_cdtB
mkdir -p slurm_logs cfml_snps_out

TREE="iqtree/zang_core.treefile"
ALN="PIRATE_out/core_alignment.fasta"
OUT="cfml_core_out/cfml_core"

test -s "$TREE"
test -s "$ALN"

ClonalFrameML "$TREE" "$ALN" "$OUT"
```

Given OUT="cfml_core_out/cfml_core", ClonalFrameML will produce (among others):
- cfml_core_out/cfml_core.labelled_tree.newick
- cfml_core_out/cfml_core.importation_status.txt
- cfml_core_out/cfml_core.em.txt
- 
*ClonalFrameML inferred:*

- Genome-wide recombination rate relative to mutation (r/m ≈ 0.62)
- Low recombination initiation rate (ρ/θ ≈ 0.003)
- Recombination tract lengths on the order of a few hundred base pairs

Together, these indicate frequent but relatively short homologous recombination events, consistent with prior observations in Campylobacter populations.

## Mask recombination tracts (CF-ML → masked alignment)

To generate a recombination-masked core alignment for downstream phylogenetic inference and association testing, we applied a masking script that converts CF-ML inferred recombinant segments into per-isolate masked sites (default mask symbol: N).

We used the updated script cfml-maskrc_updated.py, which adds robustness and useful outputs (ID normalization, interval merging, optional ancestral masking, per-isolate masking metrics, optional SVG plotting, and rectangular alignment checks). 

### Masking approach

- Recombinant tracts were extracted from cfml_core.importation_status.txt
- Masking was applied to the original core genome alignment
- By default, only extant recombination events were masked (ancestral recombination retained)
- Masked sites were replaced with N, preserving alignment structure

This approach removes recent horizontal signal while retaining deep phylogenetic structure.

Minimal masking run (extant recombination only)

```bash
python3 cfml-maskrc_updated.py \
  cfml_core_out/cfml_core \
  --aln PIRATE_out/core_alignment.fasta \
  --out cfml_core_out/core_alignment.masked \
  --metrics cfml_core_out/recomb_masking_metrics.tsv
```

Outputs included:
- core_alignment.masked.fasta – recombination-masked alignment
- recomb_masking_metrics.tsv – per-isolate recombination burden
- recomb_regions.tsv – genomic coordinates of inferred recombination

*Optional SVG visualisations were generated to illustrate recombination density along the genome.*
*Optional:* also mask ancestral segments inferred on internal branches (not just leaf-specific “extant” imports), add:

```
  --mask-ancestral
```

*Commentary:*
Masking ancestral segments is more aggressive and may be appropriate if your downstream method is sensitive to deep recombination signal; for many applications, extant-only masking is a good default.

- Alignment length preserved across all isolates
- Recombination tracts replaced by ambiguity characters (N)
- Substantial heterogeneity in recombination burden across isolates
- Masking removed large recombinant blocks while preserving phylogenetic signal

This masked alignment represents a recombination-aware substrate suitable for downstream evolutionary inference.

## Downstream: phylogeny on recombination-masked alignment

Once masking is complete, re-infer an ML tree on the masked alignment to obtain branch lengths/topology less driven by homologous recombination:

```bash
iqtree2 \
  -s cfml_core_out/core_alignment.masked.fasta \
  -m GTR+F+I+G4 \
  -T 40 \
  -B 1000 \
  --alrt 1000 \
  --prefix zang_core.masked \
  --seed 12345 \
  --safe
```

*Interpretation:*
This “masked-core” tree is typically the best default for host-association, trait mapping, and selection/structure analyses where recombination can inflate apparent homoplasy or distort branch lengths.

---

## Genome-wide asociation with PYSEER (mammal isoaltes vs bird isolates)

This section describes genome-wide association analyses performed using pyseer, correcting for population structure using a recombination-masked core-genome phylogeny.

Unless otherwise stated, GWAS were performed using fixed effects (MDS covariates derived from patristic distances). This approach is widely used in bacterial GWAS and avoids instability observed with LMMs on large, strongly structured datasets.

### Prepare PYSEER input files

#### Patristic distance matrix from the masked tree
Build patristic distance matrix from the masked tree:

```bash
python - <<'EOF' > zang_core_masked.dist.tsv
from ete3 import Tree
import sys

t = Tree("zang_core_masked.treefile", format=1)

leaves = list(t.iter_leaves())
names  = [leaf.name for leaf in leaves]
n = len(leaves)

print("\t" + "\t".join(names))

for i, a in enumerate(leaves, 1):
    if i % 50 == 0 or i in (1, n):
        print(f"[progress] {i}/{n} rows", file=sys.stderr, flush=True)
    row = [a.name]
    for b in leaves:
        row.append(str(t.get_distance(a, b)))
    print("\t".join(row))
EOF
```

OR convert the IQ-TREE ML distance matrix:

```bash
python - <<'EOF'
import sys

infile  = "zang_core_masked.mldist"
outfile = "zang_core_masked.dist.tsv"

# Read PHYLIP distance matrix (first line n, then n rows: name + n floats)
with open(infile) as f:
    first = f.readline().strip()
    n = int(first)
    rows = []
    names = []

    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        name = parts[0]
        vals = parts[1:]
        # Some PHYLIP writers wrap rows across lines; IQ-TREE usually doesn't, but be safe:
        while len(vals) < n:
            nxt = f.readline()
            if not nxt:
                break
            vals += nxt.strip().split()
        if len(vals) != n:
            raise ValueError(f"{name}: expected {n} distances, got {len(vals)}")
        names.append(name)
        rows.append(vals)

if len(rows) != n:
    raise ValueError(f"Expected {n} rows, got {len(rows)}")

with open(outfile, "w") as out:
    out.write("\t" + "\t".join(names) + "\n")
    for name, vals in zip(names, rows):
        out.write(name + "\t" + "\t".join(vals) + "\n")

print(f"Wrote {outfile} ({n}x{n})", file=sys.stderr)
EOF
```

This file is used to derive population structure covariates.

#### Phenotype file (bird vs mammal)

```bash
python - <<'EOF'
import pandas as pd

df = pd.read_csv("phenotypes.csv")

df["Host"] = df["Host"].astype(str).str.strip().str.lower()
df = df[df["Host"].isin(["birds", "mammal"])].copy()

df["phenotype"] = (df["Host"] == "birds").astype(int)

df[["id", "phenotype"]].to_csv(
    "pyseer_birds_vs_mammal.pheno.tsv",
    sep="\t",
    index=False,
    header=False
)

print(df["phenotype"].value_counts())
EOF
```

Phenotype encoding:
- 1 = bird
- 0 = mammal

### Run GWAS with genes, SNPs and unitigs

#### 1. Genes (presence/absence)

Gene presence/absence GWAS (PIRATE)
Generate a binary Rtab file (required by pyseer)

We convert PIRATE allele calls to a binary presence/absence matrix using the [PIRATE script](https://github.com/SionBayliss/PIRATE/blob/master/tools/convert_format/PIRATE_to_Rtab.pl).

```
PIRATE_to_Rtab.pl \
  --input  PIRATE_out/PIRATE.unique_alleles.tsv \
  --output pyseer/genes.Rtab \
  --samples pyseer/gwas.samples \
  --low  0.01 \
  --high 0.99 \
  --all
```

As not all isolate genomes from the dataset used to construct the pangenome will be used in the GWAS, we can subset the .Rtab to just include those that will be used in the GWAS (xxx mammal isolates vs xxx bird isoaltes): 

Create the sample list from the phenotype file
```bash
cut -f1 pyseer/pyseer_birds_vs_mammal.pheno.in_masked.tsv > pyseer/gwas.samples
wc -l pyseer/gwas.samples
head pyseer/gwas.samples
```

(That wc -l should match the “Read 962 phenotypes / Analysing 962 samples…”)

Subset genes.Rtab down to those samples

```bash
python - <<'EOF'
keep = set(x.strip() for x in open("pyseer/gwas.samples") if x.strip())

inp  = "pyseer/genes.Rtab"
outp = "pyseer/genes.subset.Rtab"

with open(inp) as fin:
    header = fin.readline().rstrip("\n").split("\t")
    idx = [0] + [i for i, s in enumerate(header[1:], start=1) if s in keep]
    new_header = [header[i] for i in idx]

    with open(outp, "w") as fout:
        fout.write("\t".join(new_header) + "\n")
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            fout.write("\t".join(parts[i] for i in idx) + "\n")

print("Wrote", outp, "with", len(new_header)-1, "samples")
EOF
```

#### GWAS using fixed effects

Gene presence/absence association was tested using pyseer under a fixed-effects logistic regression model, with population structure controlled using the first 10 multidimensional scaling components derived from patristic distances on a recombination-masked core genome phylogeny.

**Inputs required**
- pyseer/pyseer_birds_vs_mammal.pheno.in_masked.tsv (2-column TSV, no header: sample_id phenotype)
- pyseer/genes.Rtab (binary Rtab; rows=variants, columns=samples)
- pyseer/zang_core_masked.dist.tsv (square, tab-delimited distance matrix with header)

*Important:* ``--pres`` must be a binary Rtab for fixed-effects GWAS. The two-line “variant–sample” list format will trigger ``ValueError: Rtab file not binary.``

#### Step 1 — Build MDS cache

```slurm
#!/bin/bash
#SBATCH --job-name=pyseer_mds
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=60G

set -euo pipefail
export BASHRCSOURCED=1
source ~/.bashrc
conda activate pyseer

cd /groups/cooperma/bpascoe/Zang_cdtB
mkdir -p slurm_logs pyseer

pyseer \
  --phenotypes pyseer/pyseer_birds_vs_mammal.pheno.in_masked.tsv \
  --pres pyseer/dummy.pres \
  --distances pyseer/zang_core_masked.dist.tsv \
  --mds classic --max-dimensions 10 \
  --save-m pyseer/zang_core_masked.mds.pkl \
  --cpu 1 \
  --output pyseer/mds_dummy.tsv
```

#### Step 2 — Run gene GWAS using the saved MDS (lower memory)

```bash
pyseer \
  --phenotypes pyseer/pyseer_birds_vs_mammal.pheno.in_masked.tsv \
  --pres pyseer/genes.subset.Rtab \
  --load-m pyseer/mds_components.tsv.pkl \
  --min-af 0.01 --max-af 0.99 \
  --cpu 1 \
  --uncompressed \
  > pyseer/gwas_genes_birds_vs_mammal.mds.tsv
```

#### Step 3 - Multiple-testing correction (Bonferroni)

```bash
python - <<'EOF'
import pandas as pd

df = pd.read_csv("pyseer/gwas_genes_birds_vs_mammal.mds.tsv", sep="\t")

pcol = "lrt-pvalue" if "lrt-pvalue" in df.columns else "pvalue"
bonf = 0.05 / len(df)

print("rows:", len(df))
print("Bonferroni:", bonf)
print("Significant hits:", (df[pcol] < bonf).sum())

print("\nTop hits:")
print(df.sort_values(pcol).head(10)[["variant","af","beta",pcol]].to_string(index=False))
EOF

```

**Expected outputs**

Cached MDS components derived from the distance matrix (reuse this for SNP and unitig runs on the same sample set).
``pyseer/zang_core_masked.mds.pkl``

Main GWAS results table (one row per tested variant/allele), including effect size, SE, and p-values.
``pyseer/gwas_genes_birds_vs_mammal.mds.tsv``

Bonferroni-significant hits only (sorted by p-value).
``pyseer/gwas_genes_birds_vs_mammal.mds.bonf_sig.tsv``

- 961 samples analysed (birds vs mammals, masked, harmonised)
- MDS cache loaded: Loaded projection with dimension (961, 479) → confirms you’re using the same population structure across all variants
- 8,309 gene alleles tested
- No silent filtering
- Fixed-effects model with 10 MDS covariates

**Results summary**
- Multiple testing
- Tests: 8,309
- Bonferroni threshold: 6.02 × 10⁻⁶
- 211 significant gene alleles

Gene-level GWAS identified 211 PIRATE alleles significantly associated with host (bird vs mammal) after Bonferroni correction (α = 0.05), including multiple alleles from recurrent gene families.
- Strong effects (|β| ≈ 2–13)
- A mix of high- and low-frequency alleles
- Both bird- and mammal-associated directions

*Interpretation of columns*

variant | af | filter-pvalue | lrt-pvalue | beta | beta-std-err | intercept | PC1..PC10 | notes
- variant → PIRATE allele ID
- af → allele frequency in analysed samples
- beta → log-odds effect (positive = bird-associated; negative = mammal-associated)
- lrt-pvalue → the one you report
- filter-pvalue → pre-filtering check (ignore for reporting)
- PC1..PC10 → population structure covariates
- notes → usually empty or flags
- **Use lrt-pvalue for significance.**

#### 2. SNP-based (VCF) GWAS

**Required inputs**

- You need a VCF aligned to the same sample IDs as: pyseer/gwas.samples (962 IDs)
- pyseer/mds_components.tsv.pkl (projection for 961 samples)

Files we already have: 
- pyseer/gwas.samples
- pyseer/mds_components.tsv.pkl
- (optional) pyseer/zang_core_masked.dist.tsv if you ever rebuild MDS

File you need to locate / create:
- snps.vcf.gz (and its .tbi)

*Where are SNPs likely to come from?*
- PIRATE_out/core_alignment.snps.fasta exists, but that’s FASTA, not VCF

We can generate one from the core alignment we have (PIRATE_out/core_alignment.fasta; or the masked version) using snp-sites:

```bash
snp-sites -v -o pyseer/core_snps.vcf PIRATE_out/core_alignment.fasta
bgzip -c pyseer/core_snps.vcf > pyseer/core_snps.vcf.gz
tabix -p vcf pyseer/core_snps.vcf.gz
```

*Important:*
If you want SNPs from the masked alignment instead, use cfml_core_out/core_alignment.masked.fasta — just be aware masking inserts Ns, which is fine but changes missingness patterns.

### Step 1 - Subset VCF to your samples

Using the VCF we just created:
```
bcftools view -S pyseer/gwas.samples -Oz -o pyseer/core_snps.birds_vs_mammal.vcf.gz pyseer/core_snps.vcf.gz
tabix -f -p vcf pyseer/core_snps.birds_vs_mammal.vcf.gz
```

### Step 2 -  Run pyseer SNP GWAS with the same MDS cache
Use the same model as genes: fixed effects + cached MDS.

```slurm
#!/bin/bash
#SBATCH --job-name=pyseer_snps_BvM
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=120G

set -euo pipefail
export BASHRCSOURCED=1
source ~/.bashrc
conda activate pyseer

cd /groups/cooperma/bpascoe/Zang_cdtB
mkdir -p slurm_logs pyseer

pyseer \
  --phenotypes pyseer/pyseer_birds_vs_mammal.pheno.in_masked.tsv \
  --vcf pyseer/snps.birds_vs_mammal.vcf.gz \
  --load-m pyseer/mds_components.tsv.pkl \
  --min-af 0.01 --max-af 0.99 \
  --cpu ${SLURM_CPUS_PER_TASK} \
  --uncompressed \
  > pyseer/gwas_snps_birds_vs_mammal.mds.tsv
```

### Step 3 - Bonferroni + top hits

```bash
python - <<'EOF'
import pandas as pd
fn="pyseer/gwas_snps_birds_vs_mammal.mds.tsv"
df=pd.read_csv(fn, sep="\t")

pcol = "lrt-pvalue" if "lrt-pvalue" in df.columns else ("pvalue" if "pvalue" in df.columns else None)
print("rows:", len(df))
print("cols:", df.columns.tolist()[:20])
if pcol is None:
    raise SystemExit("No p-value column found")

bonf = 0.05/len(df)
print("Bonferroni:", bonf)
print("n_sig:", (df[pcol] < bonf).sum())

print("\nTop hits:")
print(df.sort_values(pcol).head(10)[["variant","af","beta",pcol]].to_string(index=False))
EOF
```

#### 3. Unitig-based GWAS

unitig presence/absence GWAS using the same fixed-effects structure correction as gene + SNP GWAS.

#### Required inputs
- `pyseer/pyseer_birds_vs_mammal.pheno.in_masked.tsv` (2-col, no header)
- `pyseer/mds_components.tsv.pkl` (cached MDS covariates)
- `pyseer/gwas.samples` (sample IDs used in GWAS)
- assemblies: `contigs/<ID>.fas`

#### Build paths list for assemblies
- ``pyseer/gwas.samples``
- ``contigs in contigs/<ID>.fas``

```bash
cd /groups/cooperma/bpascoe/Zang_cdtB

python - <<'EOF'
import os
samples=[x.strip() for x in open("pyseer/gwas.samples") if x.strip()]
paths=[]
missing=[]
for s in samples:
    p=f"contigs/{s}.fas"
    (paths if os.path.exists(p) else missing).append(p)
open("pyseer/unitigs/paths.txt","w").write("\n".join(paths)+"\n")
print("paths:",len(paths))
print("missing:",len(missing))
if missing: print("example missing:", missing[:5])
EOF
```

#### Build unitigs + presence/absence (unitig-caller)
*Note:* this step is RAM-intensive because it builds a coloured compacted de Bruijn graph (Bifrost). Expect large memory requirements.

```slurm
#!/bin/bash
#SBATCH --job-name=unitigs_call_k31
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=450G

set -euo pipefail
source ~/.bashrc
conda activate unitig-counter

cd /groups/cooperma/bpascoe/Zang_cdtB
mkdir -p slurm_logs pyseer/unitigs_call

unitig-caller \
  --call \
  --refs pyseer/unitigs/refs_962.txt \
  --kmer 31 \
  --out pyseer/unitigs_call/k31 \
  --rtab
```

*Output:* This should produce an Rtab (presence/absence) suitable for pyseer --pres / --kmers input (depending on chosen output).

#### Run unitig GWAS

Run using the same cached MDS as genes/SNPs:

```bash
pyseer \
  --phenotypes pyseer/pyseer_birds_vs_mammal.pheno.in_masked.tsv \
  --kmers pyseer/unitigs_call/k31.unitigs.Rtab \
  --load-m pyseer/mds_components.tsv.pkl \
  --min-af 0.01 --max-af 0.99 \
  --cpu 8 \
  --uncompressed \
  > pyseer/gwas_unitigs_birds_vs_mammal.mds.tsv
```

(Adjust --kmers/--pres path to the actual unitig-caller output filename.)

#### 4. Combining GWAS results across variant types

Gene-, SNP-, and (optionally) unitig-based GWAS results were combined into a single harmonised table using a custom Python script. This ensures consistent multiple-testing correction, effect-size reporting, and downstream visualisation across variant classes.

```
python scripts/combine_pyseer_results.py \
  --genes pyseer/gwas_genes_birds_vs_mammal.mds.tsv \
  --snps  pyseer/gwas_snps_birds_vs_mammal.mds.tsv \
  --out-prefix pyseer/birds_vs_mammal.combined
```

Outputs:
- birds_vs_mammal.combined.all.tsv
- birds_vs_mammal.combined.sig_q0.05.tsv
- birds_vs_mammal.combined.top50.tsv

#### 5. Mapping GWAS hits to the PIRATE pangenome

Significant gene-level GWAS hits were mapped back to PIRATE allele and gene-family definitions using PIRATE.unique_alleles.tsv.

```
python scripts/annotate_gene_hits_pirate.py \
  --gwas pyseer/gwas_genes_birds_vs_mammal.mds.tsv \
  --pirate-alleles PIRATE_out/PIRATE.unique_alleles.tsv \
  --out-prefix pyseer/birds_vs_mammal.genes \
  --sig-q 0.05
```

#### 6. Collapsing GWAS results to PIRATE gene families

Allele-level associations were collapsed to PIRATE gene families by retaining the most significant allele per family and recording the number of alleles tested.

```
python scripts/collapse_gwas_to_families.py \
  PIRATE_out/PIRATE.unique_alleles.tsv \
  pyseer/birds_vs_mammal.genes.annotated_hits.tsv \
  pyseer/birds_vs_mammal.family_level.tsv
```

#### 7. Rationalising gene annotations using reference genomes

To provide biologically interpretable gene names, PIRATE representative protein sequences were aligned against published reference proteomes (NCTC11168, 81-176, 81116). The best hit per gene family was retained.

```
makeblastdb -in pyseer/refs/campy_refs.faa -dbtype prot

blastp \
  -query PIRATE_out/representative_sequences.faa \
  -db pyseer/refs/campy_refs \
  -evalue 1e-20 \
  -max_target_seqs 5 \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  > pyseer/pirate_reps.vs_refs.blast.tsv
```

```
python scripts/make_pirate_family_annotations.py \
  --blast pyseer/pirate_reps.vs_refs.blast.tsv \
  --ref-faa reference_genomes/*.protein.faa \
  --ref-label NCTC11168 81_176 81116 \
  --out pyseer/pirate_family_annotations.tsv
```

These annotations were merged into the family-level GWAS results, with PIRATE representative gene names used as fallback when no reference match was available.

```
python scripts/add_pirate_rep_annotations.py \
  --rep-faa PIRATE_out/representative_sequences.faa \
  --families pyseer/birds_vs_mammal.family_level.annotated.clean.tsv \
  --out pyseer/birds_vs_mammal.family_level.annotated.clean.plus_pirate.tsv \
  --make-bestnames
```

**Key fields:**
- best_gene
- best_product
- pirate_family
- n_alleles
- p_lrt

#### 8. Final GWAS reporting and visualisation

All plots and summary tables were regenerated using the fully annotated, family-level GWAS results.

```
python scripts/gwas_report.py \
  --combined pyseer/birds_vs_mammal.combined.all.tsv \
  --families pyseer/birds_vs_mammal.family_level.annotated.clean.plus_pirate.tsv \
  --outdir pyseer/report_pack_birds_vs_mammal \
  --sig-q 0.05 \
  --top-n 50
```

#### Results

Genome-wide association analysis identified a limited number of robust host-associated loci after correction for population structure. The strongest gene-family signal corresponded to the cytolethal distending toxin locus (cdtA/B), with multiple low-frequency alleles showing large effect sizes, consistent with host-restricted functional adaptation. In addition, significant associations were concentrated in surface glycan biosynthesis pathways, including LOS and CPS-associated genes, often occurring as multi-gene clusters. In contrast, few housekeeping genes remained significant after correction, indicating that residual lineage structure did not dominate the analysis. Together, these results suggest that mammal–bird adaptation in C. jejuni is driven by selective tuning of host-interaction loci rather than broad metabolic divergence.

**cdtB**

- very strong effect size (β ≈ −6.5)
- low allele frequency (~4–5%)
- multiple alleles (n≈5) within the family

This is exactly the pattern you expect if:
- cdtB alleles are not lineage markers
- but instead rare, host-restricted functional variants
- maintained under negative frequency–dependent selection

**Surface glycan biology dominates the rest of the signal**

Outside cdtB, the GWAS signal concentrates in:
- LOS / CPS / glycosylation pathways
- genes involved in: sugar transfer, outer-membrane structure, host–immune interface

These hits often appear as **clusters of adjacent gene families**
- rather than isolated single genes
- consistent with locus-level selection

*Interpretation:*

This strongly supports:
- host-specific selection on cell surface architecture
- rather than generic metabolic adaptation
- and aligns with host immune recognition, gut niche differences & known Campylobacter host tropism biology

**Few housekeeping genes survive correction**

Genes like:
- purB
- dxr / dxs
- gmhA

either:
- drop out after correction
- or appear weak and isolated

👉 This argues against residual lineage confounding.

**No explosion of high-frequency “easy” hits**

We do not see:
- hundreds of ubiquitous genes
- or clear CC markers dominating the results

👉 This shows the MDS-based structure correction worked.

---
## Reproducibility, software versions, and resources

All analyses were performed on an HPC cluster using SLURM.

### Software
- Prokka v1.14.x — [https://github.com/tseemann/prokka](https://github.com/tseemann/prokka)
- PIRATE v1.0.x — [https://github.com/SionBayliss/PIRATE](https://github.com/SionBayliss/PIRATE)
- PIRATE additional analyses scripts — [pirate_host_association.py](https://github.com/Benizao1980/Mammal-adaptation-in-Campylobacter-jejuni/blob/main/scripts/pirate_host_association.py); [pirate_plotting.py](https://github.com/Benizao1980/Mammal-adaptation-in-Campylobacter-jejuni/blob/main/scripts/pirate_plotting.py)
- IQ-TREE2 v2.3.6 — [https://github.com/iqtree/iqtree2](https://github.com/iqtree/iqtree2)
- ClonalFrameML v1.0.x — [https://github.com/xavierdidelot/ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML)
- ClonalFrameML masking script — [cfml-maskrc_updated.py](https://github.com/Benizao1980/Mammal-adaptation-in-Campylobacter-jejuni/blob/main/scripts/cfml-maskrc_updated.py)
- PYSEER post-processing scripts: 
- Python v3.9+ — [https://www.python.org](https://www.python.org)
- Conda — [https://docs.conda.io](https://docs.conda.io)

### Reproducibility notes
- All genomes processed with identical parameters
- FigShare link to genomes -
- PubMLST shared proejct - 
- microreact [tree](https://microreact.org/project/gy6keE6LWQ4As8neqw9Zwz-zang-and-pascoe-et-almammal-adaptation-in-c-jejuni)
- Modular Conda environments used throughout, e.g.
```
conda create -n prokka -c conda-forge -c bioconda prokka=1.14
conda create -n pirate -c conda-forge -c bioconda pirate
conda create -n iqtree -c conda-forge -c bioconda iqtree=2.3.6
conda create -n clonalframeml -c conda-forge -c bioconda clonalframeml
conda create -n pyseer -c conda-forge -c bioconda \
    pyseer ete3 pandas numpy scipy statsmodels scikit-learn
```

### Please cite: 
TBC <link to manuscript>
