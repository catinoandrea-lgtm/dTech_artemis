# dTech Pipeline v2.0 — Chromosome 22 Edition

> Companion code for:  
> **"Neuronally-enriched ARTEMIS isoform ENST00000470530 targets hairpin-forming sequences at AluY fragility sites: a computational framework for microdeletion pathogenesis"**  
> *Andrea Catino*

---

## Overview

**dTech** is a post-processing annotation pipeline for genomic anchor files produced by the dTech scanner. It assigns functional annotations to nucleosome-positioning anchor sequences on human chromosome 22 (hg38), with special focus on identifying **INT_YOYO_LOCK** elements — hairpin-forming anchors enriched at TOP2B cleavage sites and associated with the neurally-restricted ARTEMIS isoform ENST00000470530.

This repository contains the chr22-specific version used for the analysis in the paper above.

---

## Key Paper Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| K_SIZE | 144 bp | Canonical nucleosome length |
| STEP | 72 bp | Sliding window step (half-nucleosome) |
| GC_MIN / GC_MAX | 0.35 / 0.65 | GC content filter |
| ENTROPY_MIN | 1.87 bits | Shannon entropy threshold |
| TM_MIN | 68.0 °C | Melting temperature threshold |
| FFT_THR | 1.5 | A/T periodicity from Fourier transform |
| CURV_THR | 1.2 | Intrinsic DNA curvature score |
| HAIRPIN_LOCK_THR | ≥ 8 | Minimum hairpin score for INT_YOYO_LOCK |
| PROMOTER_WIN | 5000 bp | Upstream window for promoter classification |
| **Total anchors (chr22)** | **361,015** | After biophysical filters |
| **INT_YOYO_LOCK (chr22)** | **151** | Hairpin-forming anchors |

---

## Results Reproduced

Running this pipeline on the chr22 scan output should reproduce **Table 1** of the paper:

| TAG | n | % | hairpin_mean | hairpin_max |
|-----|---|---|-------------|------------|
| PROMOTER | 118,220 | 32.7% | 2.91 | 31 |
| INTRON | 104,448 | 28.9% | 2.92 | 12 |
| INT_YOYO | 87,653 | 24.3% | 2.88 | 7 |
| EXON | 22,245 | 6.2% | 2.93 | 12 |
| CENTROMERE | 16,932 | 4.7% | 2.96 | 10 |
| CDS | 11,366 | 3.1% | 2.94 | 9 |
| **INT_YOYO_LOCK** | **151** | **0.04%** | **8.47** | **20** |

The INT_YOYO_LOCK hairpin mean (8.47) is **~3× higher** than all other categories (2.88–2.96).

---

## Requirements

```bash
Python >= 3.9
pandas >= 1.5
requests >= 2.28
```

Install dependencies:
```bash
pip install pandas requests
```

No GPU required. Designed to run on a standard laptop with:
- **2 CPU cores** (configurable via `N_CORES`)
- **3 GB RAM** (configurable via `RAM_LIMIT`)

---

## Directory Structure

```
/root/dna/
├── output_dtech/           ← INPUT: dTech scan CSV files (from VS tool)
│   └── human_NC_000022*.csv.gz
├── organized_results/      ← Organized copy of scan files
│   └── human/
│       └── Main_Chromosomes/
├── output_annotated/       ← OUTPUT: annotated files + summary
│   ├── bio_human_NC_000022*.csv.gz
│   └── dtech_chr22_summary.csv   ← paper Table 1
├── gff_database/           ← Downloaded GFF3 (auto-managed)
│   └── human.gff.gz
└── pipeline.log            ← Execution log
```

---

## Usage

### 1. Prepare input files

Place your dTech scan output files in `/root/dna/output_dtech/`. Files must be named following the convention:
```
{species}_{chr_accession}.csv.gz
```
For chromosome 22:
```
human_NC_000022.1.csv.gz
```

Expected columns in the CSV:
| Column | Type | Description |
|--------|------|-------------|
| `pos` | int | Anchor start position (1-based, bp) |
| `seq` | str | 144 bp anchor sequence |
| `score` | float | dTech composite biophysical score |
| `hairpin_score` | int | (optional) Pre-computed hairpin score |

### 2. Run the pipeline

```bash
python3 dtech_pipeline_chr22.py
```

The pipeline will:
1. Download the hg38 GFF3 annotation (~78 MB, one-time download)
2. Watch for new scan files in `output_dtech/`
3. Annotate each anchor with TAG, DOMAIN, GENE, IS_LOCK, hairpin_score
4. On exit (Ctrl+C): print and save the summary table

### 3. Check results

```bash
# Summary table (matches paper Table 1)
cat /root/dna/output_annotated/dtech_chr22_summary.csv

# Annotated anchors
zcat /root/dna/output_annotated/bio_human_NC_000022*.csv.gz | head -5
```

---

## Output Columns (annotated CSV)

| Column | Description |
|--------|-------------|
| `pos` | Anchor start position (bp) |
| `seq` | 144 bp sequence |
| `score` | dTech composite biophysical score |
| `TAG` | Functional class (see values below) |
| `DOMAIN` | Functional keyword (BRAIN, IMMUNE, etc.) |
| `GENE` | Nearest gene symbol or accession |
| `IS_LOCK` | True if TAG == INT_YOYO_LOCK |
| `hairpin_score` | Intrastrand hairpin potential (count of 20-bp windows with ≥8 complementary bases) |

### TAG values

| TAG | Description |
|-----|-------------|
| `PROMOTER` | Within 5 kb upstream of a TSS |
| `INTRON` | Within gene body, no YOYO pattern |
| `INT_YOYO` | Within gene body, YOYO tandem repeat present |
| `INT_YOYO_LOCK` | INT_YOYO + hairpin_score ≥ 8 **(paper focus)** |
| `EXON` | Within annotated exon |
| `CDS` | Within coding sequence |
| `CENTROMERE` | Within centromeric region |
| `DESERT` | Intergenic, no nearby feature |

---

## Summary Table Output (dtech_chr22_summary.csv)

The pipeline generates a CSV with one row per TAG class:

```
TAG,n_anchors,pct_total,hairpin_mean,hairpin_median,hairpin_max,n_INT_YOYO_LOCK,score_mean
PROMOTER,118220,32.73,2.910,0.0,31,0,23.081
INTRON,104448,28.93,2.920,0.0,12,0,22.528
INT_YOYO,87653,24.28,2.880,0.0,7,0,23.851
INT_YOYO_LOCK,151,0.04,8.470,8.0,20,151,23.853
...
TOTAL,361015,100.0,...
```

---

## Configuration

Key parameters are defined at the top of `dtech_pipeline_chr22.py`:

```python
N_CORES        = 2     # parallel annotation workers
RAM_LIMIT      = 3     # GB; pause annotation if exceeded
CHUNK_ROWS     = 5000  # rows per multiprocessing chunk
HAIRPIN_LOCK_THR = 8   # INT_YOYO_LOCK hairpin threshold
PROMOTER_WIN   = 5000  # bp upstream of TSS = promoter
```

---

## Data Sources

| Data | Source | Accession/URL |
|------|--------|--------------|
| hg38 GFF3 annotation | NCBI RefSeq | GCF_000001405.40_GRCh38.p14 |
| TOP2B ChIP-seq | ENCODE | GSE141528 (GM12878) |
| GTEx isoform expression | GTEx v8 | gtexportal.org |
| TCGA-GBM expression | cBioPortal | gbm_tcga |
| GWAS catalog | UCSC Table Browser | gwasCatalog track (hg38) |
| gnomAD variants | gnomAD v4 | ENSG00000152457 |

---

## Citation

If you use this pipeline, please cite:

```
Catino A. Neuronally-enriched ARTEMIS isoform ENST00000470530 targets hairpin-forming 
sequences at AluY fragility sites: a computational framework for microdeletion pathogenesis.
[Journal name, year, doi: pending]
```

---

## License

MIT License — see `LICENSE` file.

---

## Contact

Andrea Catino - Independent Researcher - catino.andrea@gmail.com 
For questions about the pipeline or the paper, please open a GitHub Issue.
