<img src="app/www/titan_logo_blue.svg" alt="TITAN" height="72">

**Tumor Immunopeptidomics Target Atlas of Non‑canonical ORFs**

A Shiny app for prioritizing non-canonical ORF (ncORF) peptide candidates identified by immunopeptidomics mass spectrometry. TITAN integrates ribo-seq translation evidence, RNA-seq tumour expression, GTEx normal-tissue specificity, TCGA pan-cancer coverage and Ribocrypt external translation data into a single scored, interactive candidate table. The ORF Detail tab provides per-candidate safety checks: canonical cross-reactivity (Biostrings exact/1-mismatch/2-mismatch matching against the Ensembl 114 proteome) and BLAST homology (blastp), all run offline against pre-built Ensembl 114 reference databases.

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
                   "stringr", "shinyWidgets", "data.table", "digest", "yaml",
                   "polished", "blastula"))
# Bioconductor
BiocManager::install(c("tximport", "matrixStats", "Biostrings", "rBLAST"))
```

> **Note**: `Biostrings` and `rBLAST` are required for the cross-reactivity and BLAST homology checks in the ORF Detail tab. `rBLAST` calls the bundled `blastp` binary in `app/bin/`; no system BLAST installation or Singularity image is needed at runtime. `polished` provides account login/sign-up; `blastula` will send catalog-access-request notification emails once that feature is wired up (see [Authentication (Polished)](#authentication-polished)).

---

## Authentication (Polished)

TITAN uses [Polished](https://polished.tech) for account creation and login (`app/global.R`, `app/app.R`). Configuration:

- **Sign-up**: open self-service registration (email + password + confirm password), no invite or domain restriction — the app is already IP-restricted at the network level. Set via `is_invite_required = FALSE` in `polished::global_sessions_config()`.
- **Default role**: new sign-ups are granted the `"general"` role automatically. This is a **Polished dashboard setting** for this app (Roles → default role on registration), not something set in app code.
- **Password policy**: enforced centrally by Polished, not by TITAN's code. Current dashboard setting for this app: **minimum 10 characters, must include both letters and numbers**. If this is changed in the dashboard, update this line rather than adding validation in the app.
- **Forgot / reset password**: built into Polished's default sign-in page (`sign_in_ui_default()`) — no extra wiring needed in-app. The reset email is sent by Polished itself, so its delivery depends on the **Email/SMTP settings configured in the Polished dashboard** for this app. Point those settings at the same SMTP server used by `blastula` for catalog-access notifications so both flows send from a consistent, deliverable address.
- **Required environment variables**: `POLISHED_API_KEY` (from the Polished dashboard for this app; not yet provisioned as of this branch) and optionally `POLISHED_APP_NAME` (defaults to `"titan"`).

Role-request / catalog-access-specific behavior (`app/R/catalog_access.R`) is scaffolded separately and not yet wired into sign-up.

---

## Step 0: Build reference databases (one-time, HPC)

The ORF Detail safety checks require three reference files built from public databases. Run these once before launching the app:

```bash
cd /path/to/titan
mkdir -p logs

# ~1 h, 16 GB — downloads Ensembl 114 pep, deduplicates, builds BLAST db + index
sbatch app/scripts/01_prep_ensembl_pep.sbatch

# ~30 min, 8 GB — offline biomaRt query for gene annotations
sbatch app/scripts/03_prep_annotation.sbatch
```

Outputs written to `app/ref/` (gitignored):

| File | Used for |
|---|---|
| `ref/ensembl114_pep/ensembl114_pep_index.rds` | Canonical cross-reactivity (Biostrings) |
| `ref/ensembl114_pep/ensembl114_pep.*` | BLAST homology database |
| `ref/ensembl114_pep/ensembl_gene_annotation.rds` | BLAST hit annotation |

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

1. **Data** tab — load the prepared RDS (`titan_orf_table.rds`) or use the built-in catalog data.
2. **Overview** tab — filter by biotype, ribo-seq and RNA-seq thresholds; inspect translation/expression distributions
3. **MS** tab — upload a TITAN-ready peptide file; the app matches peptides to ORF protein sequences. Peptides that match an `ORF-annotated` or `NC-variant` ORF are restricted to those canonical biotypes and excluded from ncORF evidence.
4. **Prioritization** tab — candidates are displayed one row per gene (grouped by gene × biotype × peptide set), with the highest-scoring ORF as the representative; click **+** to expand and inspect co-identified ORFs within the same group. Rows are ranked by a weighted priority score; filter by tumour specificity category or GTEx expression; click a row for the candidate detail panel (tissue Q3 flags, score profile radar, dimension breakdown).
5. **ORF Detail** tab — per-ORF safety assessment and expression context. The left column shows the protein sequence card (full predicted protein with matched MS peptides highlighted) and the cross-reactivity card (canonical self-protein matching against the Ensembl 114 human proteome — exact and 1-mismatch, gene-level collapse, isoform count). The right column shows the **BLAST homology card**, which runs `blastp` against the same Ensembl 114 proteome (debounced 750 ms), reporting hits at ≥ 50 % identity and ≥ 30 % alignment coverage annotated with gene description; click a hit row to view its expression profile. All checks are session-cached per ORF. Clicking any gene row opens a three-panel expression modal (target tumour / GTEx / TCGA).

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
