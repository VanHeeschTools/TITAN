<img src="app/www/titan_logo_blue.svg" alt="TITAN" height="72">

**Tumor Immunopeptidomics Target Atlas of Non‑canonical ORFs**

A Shiny app for prioritising non-canonical ORF (ncORF) peptide candidates identified by immunopeptidomics mass spectrometry. TITAN integrates ribo-seq translation evidence, RNA-seq tumour expression, GTEx normal-tissue specificity, TCGA pan-cancer coverage and Ribocrypt external translation data into a single scored, interactive candidate table. The ORF Detail tab provides per-candidate safety checks: canonical cross-reactivity and allergen cross-reactivity (Biostrings exact/1-mismatch matching) and BLAST homology (blastp), all run offline against pre-built Ensembl 114 reference databases.

---

## Overview

```
scripts/
  prepare_titan_inputs.R     # Builds the app RDS from upstream analyses
  prepare_ms_fragpipe.R      # Converts a FragPipe PSM file to TITAN-ready TSV
app/
  app.R                      # Shiny app (UI + server)
  global.R                   # Package loading, reference data, constants
  www/                       # CSS, JS, SVG assets
  data/                      # Runtime data (gitignored; generated locally)
  ref/                       # Reference databases (gitignored; built by prep scripts)
  bin/                       # Bundled blastp binary
  scripts/                   # One-time reference preparation scripts (sbatch)
documentation/
  input_fields               # Column dictionary for the orf_table
```

---

## Requirements

R packages (install once, preferably in a Singularity image):

```r
install.packages(c("shiny", "bslib", "DT", "plotly", "dplyr", "tidyr",
                   "stringr", "shinyWidgets", "data.table", "digest", "yaml"))
# Bioconductor
BiocManager::install(c("tximport", "matrixStats", "Biostrings", "rBLAST"))
```

> **Note**: `Biostrings` and `rBLAST` are required for the cross-reactivity and BLAST homology checks in the ORF Detail tab. `rBLAST` calls the bundled `blastp` binary in `app/bin/`; no system BLAST installation or Singularity image is needed at runtime.

---

## Step 0: Build reference databases (one-time, HPC)

The ORF Detail safety checks require three reference files built from public databases. Run these once before launching the app:

```bash
cd /path/to/titan
mkdir -p logs

# ~1 h, 16 GB — downloads Ensembl 114 pep, deduplicates, builds BLAST db + index
sbatch app/scripts/01_prep_ensembl_pep.sbatch

# ~30 min, 4 GB — downloads UniProt allergen set (KW-0020, reviewed)
sbatch app/scripts/02_prep_allergen.sbatch

# ~30 min, 8 GB — offline biomaRt query for gene annotations
sbatch app/scripts/03_prep_annotation.sbatch
```

Outputs written to `app/ref/` (gitignored):

| File | Used for |
|---|---|
| `ref/ensembl114_pep/ensembl114_pep_index.rds` | Canonical cross-reactivity (Biostrings) |
| `ref/ensembl114_pep/ensembl114_pep.*` | BLAST homology database |
| `ref/ensembl114_pep/ensembl_gene_annotation.rds` | BLAST hit annotation |
| `ref/allergen_uniprot/allergen_index.rds` | Allergen cross-reactivity (Biostrings) |

The app starts without these files (safety check cards show an informational prompt); build them to activate the full ORF Detail functionality.

---

## Step 1: Prepare the app data (HPC)

Edit the `CONFIGURATION` block at the top of `scripts/prepare_titan_inputs.R` to point to your project paths:

| Variable | Description |
|---|---|
| `BASE_TARGET` | Root of your tumour project directory |
| `BASE_TITAN` | Output directory (usually `app/data/`) |
| `TARGET_TUMOR_TYPE` | The `tissue_type` label used for tumour samples in the GTEx coldata |
| `PATHS` | Individual file paths for ncORF table, ribo-seq, RNA-seq, GTEx and TCGA inputs |

Then run on an interactive RStudio session or submit to the HPC:

```bash
#!/bin/bash
#SBATCH --job-name=titan_prep
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00

singularity exec /hpc/local/Rocky8/pmc_vanheesch/singularity_images/<image>.sif \
  Rscript /hpc/pmc_oatv/projects/tools_dev/titan/scripts/prepare_titan_inputs.R
```

This writes `app/data/titan_orf_table.rds` (the app's primary input).

### Required upstream files

| Path key | Content |
|---|---|
| `ncorfs` | Harmonised ncORF candidate table (CSV, one row per ORF) from nf_riboseq_pipeline |
| `ribo_ppm` | ORF × sample ribo-seq PPM matrix (CSV) from nf_riboseq_pipeline |
| `ribo_psites` | ORF × sample p-site count matrix (CSV) from nf_riboseq_pipeline |
| `ribocrypt_ext` | ORF × sample Ribocrypt PPM matrix (CSV) |
| `de_sig_all` | GTEx DE result table with `sig_in_all`, `low_all_tissues`, `Q3_GTEx` columns |
| `gtex_txi` | tximport object (RDS) — GTEx + tumour samples combined |
| `gtex_coldata` | Sample metadata (RDS) with `sample_id` and `tissue_type` columns |
| `tcga_txi` | tximport object (RDS) — TCGA matched tumour/normal |

---

## Step 2: Prepare MS peptides (optional, FragPipe)
MS data should have at least one column with peptide identifications (sequences). Any additional columns will just be displayed as complementary information. 

If your immunopeptidomics data comes from FragPipe, convert the PSM output before uploading to TITAN:

```bash
singularity exec <image>.sif \
  Rscript scripts/prepare_ms_fragpipe.R \
    --input /path/to/psm.tsv \
    --output /path/to/titan_ms.tsv \
    --qval 0.01
```

The script filters decoys, contaminants and poor-quality PSMs, deduplicates to one PSM per peptide sequence, and drops redundant columns. Upload the resulting TSV directly in the TITAN app.

---

## Step 3: Run the app (HPC)

Launch from a compute node:

```bash
singularity exec <image>.sif \
  Rscript -e "shiny::runApp('app/', port=3838, host='0.0.0.0')"
```

Then open a port-forwarded connection from your local machine:

```bash
ssh -L 3838:localhost:3838 <hpc-login-node>
```

Navigate to `http://localhost:3838` in your browser.

---

## App workflow

1. **Data** tab — load the prepared RDS (`titan_orf_table.rds`) or use the built-in demo dataset
2. **Overview** tab — filter by biotype, ribo-seq and RNA-seq thresholds; inspect translation/expression distributions
3. **MS** tab — upload a TITAN-ready peptide file; the app matches peptides to ORF protein sequences. Peptides that match an `ORF-annotated` or `NC-variant` ORF are restricted to those canonical biotypes and excluded from ncORF evidence.
4. **Prioritisation** tab — candidates are displayed one row per gene (grouped by gene × biotype × peptide set), with the highest-scoring ORF as the representative; click **+** to expand and inspect co-identified ORFs within the same group. Rows are ranked by a weighted priority score; filter by tumour specificity category or GTEx expression; click a row for the candidate detail panel (tissue Q3 flags, score profile radar, dimension breakdown).
5. **ORF Detail** tab — per-ORF safety assessment and expression context. The protein sequence card shows the full predicted protein with matched MS peptides highlighted. The **Cross-reactivity card** runs two checks on selection: (a) canonical self-protein matching against the Ensembl 114 human proteome (exact and 1-mismatch, gene-level collapse, isoform count) and (b) allergen matching against UniProt reviewed allergen proteins (KW-0020). The **BLAST homology card** runs `blastp` against the same Ensembl 114 proteome (debounced 750 ms), reporting hits at ≥ 50 % identity and ≥ 30 % alignment coverage annotated with gene description; click a hit row to view its expression profile. All checks are session-cached per ORF. Clicking any gene row opens a three-panel expression modal (target tumour / GTEx / TCGA).

---

## Scoring

The priority score is a weighted sum of eight dimensions (0–100 scale):

| Dimension | Signal | Default direction |
|---|---|---|
| % expressed (RNA-seq) | Tumour sample coverage | + |
| % translated (Ribo-seq) | Tumour sample coverage | + |
| Tumour specificity (GTEx) | Tumor-only=1, enriched=0.5, other=0 | + |
| GTEx expression level | Max tissue median TPM | − |
| TCGA tumour coverage | % TCGA tumour samples ≥ 1 TPM | + |
| TCGA normal expression | % TCGA peritumoral samples ≥ 1 TPM | − |
| RC primary tissue | Ribocrypt normal primary PPM | − |
| RC cell-line | Ribocrypt cell-line PPM | + |

Two presets are provided (**Cancer-specific**, **Pan-cancer**); weights can be freely adjusted via sliders.

---

## Input field dictionary

See `documentation/input_fields` for a full description of every column expected in the `orf_table`.
