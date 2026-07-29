# Data Preprocessing

This page covers how to produce the `titan_<study_id>.rds` file that TITAN loads for a given study. The entry point is `scripts/prepare_titan_inputs.R`.

---

## Overview

`prepare_titan_inputs.R` takes a study config YAML and assembles six data sections into a single RDS object:

| Section | Script step | Source |
|---|---|---|
| ncORF backbone | 1/6 | Pipeline ncORF table (CSV) |
| Target translation | 2/6 | Study ribo-seq PPM + P-site matrices |
| RiboCrypt external | 3/6 | External RiboCrypt PPM matrix |
| Target expression + GTEx | 4/6 | RNA-seq TPM + GTEx quant + DESeq2 DE results |
| TCGA | 5/6 | TCGA Salmon TPM |
| Output | 6/6 | Joined table + matrices → RDS + CSV |

---

## Usage

```bash
# On an HPC login node — submit via SLURM so it doesn't run on the submit node
sbatch --mem=64G --cpus-per-task=4 --wrap \
  "Rscript scripts/prepare_titan_inputs.R scripts/configs/my_study.yaml"
```

Or in an interactive session on a compute node:
```r
Rscript scripts/prepare_titan_inputs.R scripts/configs/my_study.yaml
```

The script outputs:
- `<output.dir>/titan_<study_id>.rds` — app input (large; typically 500 MB–2 GB)
- `<output.dir>/titan_<study_id>.csv` — flat CSV for inspection

---

## Input Files

### 1. ncORF Candidates (`paths.ncorfs`)

A CSV/TSV from the upstream ribo-seq / ncORF-discovery pipeline. Required columns:

| Column | Type | Description |
|---|---|---|
| `orf_id` | character | Unique ORF identifier (hash or stable ID) |
| `gene_id` | character | Ensembl gene ID (versioned OK; version stripped internally) |
| `gene_name` | character | Gene symbol |
| `gene_biotype` | character | Ensembl gene biotype |
| `orf_biotype_single` | character | TITAN biotype label (see legend in App Guide) |
| `protein_seq` | character | Amino acid sequence |
| `protein_length` | integer | Length in aa |
| `start_codon` | character | ATG / CTG / etc. |
| `chr` | character | Chromosome (without "chr" prefix) |
| `orf_start` | integer | Genomic start (0-based) |
| `orf_end` | integer | Genomic end (0-based) |
| `strand` | character | `+` or `-` |

### 2. Ribo-seq PPM matrix (`paths.ribo_ppm`)

CSV/TSV with:
- First column: `orf_id` (matching the ncORF table)
- Remaining columns: one per sample (column names = sample IDs)
- Values: P-sites per million (PPM) — non-negative

### 3. Ribo-seq P-site matrix (`paths.ribo_psites`)

Same format as the PPM matrix but with raw P-site counts. Used to compute `target_translation_median_psites`.

### 4. RiboCrypt external matrix (`paths.ribocrypt_ext`)

CSV/TSV of ORF × cross-study sample PPM values from the RiboCrypt database.

- First column: `orf_id`
- Column names encode sample identity: prefix-based convention determines primary tissue vs cell-line classification (see [Catalog Setup → Sample Classification](Catalog-Setup)).

### 5. GTEx quantification (`paths.gtex_quant`)

Gene × sample TPM matrix. Accepts:
- `.rds` — a tximport-style list with `$abundance` (gene × sample matrix)
- `.csv` / `.tsv` — first column = gene IDs (rownames), remaining = samples

Two usage modes:

**Mode A — combined object** (default): one file contains both tumor samples and GTEx normal samples. `paths.gtex_coldata` (see below) specifies which is which via `tissue_type`.

**Mode B — separate files**: set `paths.tumor_quant` to a separate tumor-only quantification file. `paths.gtex_quant` then contains only GTEx normals.

### 6. GTEx column data (`paths.gtex_coldata`)

RDS file: a data.frame with at minimum:

| Column | Description |
|---|---|
| `sample_id` | Must match column names of the expression matrix |
| `tissue_type` | GTEx tissue label (e.g., `Brain`, `Lung`) or cancer label (e.g., `RMS`) |

Tumor samples are identified by `tissue_type == cfg$tumor_type`. GTEx normals are all other rows with non-NA `tissue_type`.

### 7. TCGA quantification (`paths.tcga_quant`)

Same format as GTEx quant (`.rds`, `.csv`, or `.tsv`). Should include:
- Primary tumor samples: TCGA barcodes ending in `-0[1-9][A-Z]`
- Normal/peritumoral samples: TCGA barcodes ending in `-1[0-1][A-Z]`

The script auto-classifies tumor vs normal by barcode suffix — no coldata needed.

### 8. TCGA column data (`paths.tcga_coldata`)

RDS file: a data.frame with at minimum:

| Column | Description |
|---|---|
| `sample_id` | TCGA barcode |
| `tissue_type` | Cancer type code (e.g., `RMS`, `EWS`) |
| `sample_type` | `Tumor` or `Normal` |
| `group` | Display group string: e.g., `RMS Tumor`, `RMS Normal` |

### 9. GTEx DE table (`paths.de_sig_all`) — optional but strongly recommended

A tab-delimited file with pre-computed DESeq2 results from a GTEx differential expression analysis (tumor vs all 28 GTEx normal tissues). Required columns:

| Column | Description |
|---|---|
| `gene_id` | Ensembl gene ID (versioned OK) |
| `sig_in_all` | `TRUE` if DE in all 28 tissues |
| `low_all_tissues` | `TRUE` if the gene is tumor-only (low in all GTEx tissues) |
| `Q3_GTEx` | Q3 GTEx TPM across tissues (used to classify tumor-enriched: `Q3_GTEx < 1`) |

Without this file, `GTEX_tumor_only`, `GTEX_tumor_enriched`, and `GTEX_DE_sig_in_all` are all `NA` for every candidate in this study, and the Specificity column in TITAN will show "Unavailable".

### 10. Tumor quantification (`paths.tumor_quant`) — optional

Only needed in Mode B. A gene × sample TPM matrix (`.rds`, `.csv`, `.tsv`) containing only tumor samples from the study. All columns are used as tumor samples; no coldata is required.

---

## What `prepare_titan_inputs.R` Computes

### Per-ORF columns added to the backbone table

**Target translation (ribo-seq):**
- `target_translation_num_samples` — number of samples with PPM ≥ threshold
- `target_translation_pct_samples` — percentage
- `target_translation_median_PPM` / `target_translation_max_PPM`
- `target_translation_median_psites`

**RiboCrypt external:**
- `ribocrypt_primary_{num,pct,median,max}_*` — primary tissue metrics
- `ribocrypt_cell-line_{num,pct,median,max}_*` — cell-line metrics

**Target expression (RNA-seq, tumor samples):**
- `target_expression_num_samples`, `target_expression_pct_samples`
- `target_expression_median_TPM`, `target_expression_max_TPM`

**GTEx metrics (gene-level):**
- `GTEX_max_median_TPM` — maximum per-tissue median TPM across all 28 tissues
- `GTEX_median_TPM` — overall median TPM
- `GTEX_tissues_q3_gt1` — `|`-separated string of `tissue=Q3TPM` pairs where Q3 > threshold
- `GTEX_DE_sig_in_all`, `GTEX_tumor_only`, `GTEX_tumor_enriched` (from DE table)

**TCGA metrics (gene-level):**
- `TCGA_tumor_{num,pct,median,max}_samples/TPM`
- `TCGA_normal_{num,pct,median,max}_samples/TPM`

### Per-sample matrices stored in the RDS

| Element | Rows | Columns | Values |
|---|---|---|---|
| `$ribo_ppm_samples` | orf_id | ribo-seq samples | PPM |
| `$rna_tpm_mat` | gene_id (no version) | tumor RNA-seq samples | TPM |
| `$gtex_tpm_mat` | gene_id | GTEx normal samples | TPM |
| `$tcga_tpm_mat` | gene_id | TCGA tumor + normal samples | TPM |
| `$ribocrypt_mat` | orf_id | all RiboCrypt samples | PPM |

### Metadata data.frames

| Element | Columns |
|---|---|
| `$ribo_sample_meta` | `sample_id`, `condition` |
| `$rna_sample_meta` | `sample_id`, `tissue_type`, `condition` |
| `$gtex_sample_meta` | `sample_id`, `tissue_type` |
| `$tcga_sample_meta` | `sample_id`, `tissue_type`, `sample_type`, `group` |
| `$ribocrypt_sample_meta` | `sample_id`, `group` (`"Primary"` or `"Cell-line"`) |
| `$ribocrypt_meta` | list with `$primary_samples` and `$cell_line_samples` character vectors |

---

## Pipelines Used for Each Data Type

> **Note:** The specific Nextflow pipelines, module versions, and workflow configurations used to generate the input files vary by study. Fill in the gaps below with your study's pipeline details.

### ncORF Candidates
- **Pipeline:** [VanHeeschTools Ribo-seq / ncORF discovery pipeline] _(link TBD)_
- **Output format:** CSV with columns as listed above
- **Key steps:** Ribo-seq alignment → ORF calling → filtering for immunopeptidome candidates (≥ 8 aa, no canonical overlap, etc.)

### Ribo-seq PPM / P-site matrices
- **Pipeline:** Same as above; PPM normalisation applied per sample
- **Output format:** ORF × sample CSV (first column = `orf_id`)

### RiboCrypt external
- **Source:** [RiboCrypt database](https://ribocrypt.org/)
- **Version used:** _(record version/date of download)_
- **Output format:** ORF × sample CSV

### GTEx quantification
- **Source:** GTEx portal v8 or v9 _(confirm version)_
- **Processing:** Salmon TPM aligned to GRCh38, tximport aggregation to gene level
- **Output format:** tximport RDS (`$abundance` matrix) or gene × sample CSV

### GTEx DE classification
- **Tool:** DESeq2 _(version TBD)_
- **Design:** tumor vs each of 28 GTEx tissue groups; summary table of sig_in_all / low_all / Q3
- **Output format:** TSV with columns as listed above

### TCGA quantification
- **Source:** GDC data portal _(project IDs TBD)_
- **Processing:** Salmon TPM, tximport gene-level aggregation
- **Samples:** 257 tumor + 257 peritumoral/normal _(study-specific)_
- **Output format:** gene × sample CSV or tximport RDS
