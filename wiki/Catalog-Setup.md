# Catalog Setup

This page covers everything needed to make a study available in TITAN's Study Library: the per-study config YAML, the RDS preparation run, catalog registration, and optional MS peptide bundling.

---

## Step 1: Write a Study Config YAML

Create a file under `scripts/configs/<study_id>.yaml`. All path values must be absolute or resolvable from the working directory where you run the script.

### Minimal config

```yaml
# ── Identity ──────────────────────────────────────────────────────────────────
study_id:     rms_organoids         # snake_case, unique across catalog
display_name: "RMS organoids"       # shown in the Study Library
cancer_type:  RMS                   # used for filtering
cohort:       Organoids             # used for filtering
tumor_type:   RMS                   # must match tissue_type values in gtex_coldata

# ── Input paths ───────────────────────────────────────────────────────────────
paths:
  ncorfs:         /path/to/ncorfs.csv
  ribo_ppm:       /path/to/ribo_ppm.csv
  ribo_psites:    /path/to/ribo_psites.csv
  ribocrypt_ext:  /path/to/ribocrypt_ppm.csv
  gtex_quant:     /path/to/combined_quant.rds      # Mode A: tumor + GTEx in one file
  gtex_coldata:   /path/to/gtex_coldata.rds

  # Required in Mode B only (separate tumor quantification):
  # tumor_quant:  /path/to/tumor_quant.rds

  tcga_quant:     /path/to/tcga_quant.rds
  tcga_coldata:   /path/to/tcga_coldata.rds

  # Optional but strongly recommended:
  de_sig_all:     /path/to/gtex_de_summary.tsv

# ── Output ────────────────────────────────────────────────────────────────────
output:
  dir: /path/to/output/             # directory must exist or be creatable
  # rds: titan_rms_organoids.rds   # optional override; default = titan_<study_id>.rds
  # csv: titan_rms_organoids.csv   # optional override
```

### Full config reference

```yaml
study_id:     rms_organoids
display_name: "RMS organoids"
cancer_type:  RMS
cohort:       Organoids
tumor_type:   RMS
target_label: "RMS"                 # Optional: display label for the tumor type in plots

paths:
  ncorfs:         /path/to/ncorfs.csv
  ribo_ppm:       /path/to/ribo_ppm.csv
  ribo_psites:    /path/to/ribo_psites.csv
  ribocrypt_ext:  /path/to/ribocrypt_ppm.csv
  gtex_quant:     /path/to/combined_quant.rds
  gtex_coldata:   /path/to/gtex_coldata.rds
  tcga_quant:     /path/to/tcga_quant.rds
  tcga_coldata:   /path/to/tcga_coldata.rds
  de_sig_all:     /path/to/gtex_de_summary.tsv
  # tumor_quant:  /path/to/separate_tumor_quant.rds  # Mode B only

thresholds:
  expression: 1     # TPM/PPM threshold for "expressed" (default: 1)
  gtex_q3:    1     # GTEx Q3 threshold for off-tissue risk (default: 1)

# Regex applied to ribo-seq sample column names to shorten display labels.
# E.g., remove a common prefix like "sample_RMS_" to leave just sample IDs.
sample_id_regex: "^sample_RMS_"

# Optional: classify ribo-seq samples into conditions by name pattern.
# Without this, all samples get the tumor_type label.
riboseq_condition:
  pattern:      "_T_"
  match_label:   "Tumor"
  nomatch_label: "Organoid"

output:
  dir: /hpc/pmc_oatv/projects/tools_dev/titan/app/data/rms_organoids/
```

---

## Step 2: Run `prepare_titan_inputs.R`

```bash
# On the HPC — submit via SLURM from the project root
sbatch --mem=64G --cpus-per-task=4 \
  --output=logs/prepare_rms_organoids_%j.out \
  --wrap="Rscript scripts/prepare_titan_inputs.R scripts/configs/rms_organoids.yaml"
```

The script prints progress in 6 numbered steps:
```
=== TITAN input preparation ===
    Study : RMS organoids (rms_organoids)
    Config: scripts/configs/rms_organoids.yaml

[1/6] Loading ncORF candidate table...
[2/6] Computing target translation metrics (RMS ribo-seq)...
[3/6] Computing RiboCrypt external metrics...
[4/6] Loading GTEx and tumour quantification data...
[5/6] Computing TCGA expression metrics...
[6/6] Assembling final object and writing outputs...
```

On completion you will have `titan_rms_organoids.rds` (and `.csv`) in `output.dir`.

---

## Step 3: Register the Study

```bash
Rscript scripts/catalog/register_study.R \
  scripts/configs/rms_organoids.yaml \
  /path/to/titan_rms_organoids.rds
```

This script:
1. Validates the RDS using `validate_titan_rds()`.
2. Reads `study_id`, `display_name`, `cancer_type`, `cohort` from the config.
3. Extracts `n_orfs`, `n_ribo_samples`, `n_rna_samples`, `prepared_on` from the RDS.
4. Computes the RDS path **relative to the `app/` directory** (what gets stored in the catalog).
5. Appends or updates the entry in `app/data/catalog.yaml`.

Prints a diff when updating an existing entry so you can review the change before it is written.

Custom catalog path (e.g., for staging):
```bash
Rscript scripts/catalog/register_study.R \
  scripts/configs/rms_organoids.yaml \
  /path/to/titan_rms_organoids.rds \
  --catalog /path/to/custom/catalog.yaml
```

---

## The `catalog.yaml` Format

`app/data/catalog.yaml` is the app's study library. It is read once at app startup. Do not edit it by hand — use `register_study.R`.

```yaml
studies:
  - study_id:       rms_organoids
    display_name:   "RMS organoids"
    cancer_type:    RMS
    cohort:         Organoids
    rds_path:       data/rms_organoids/titan_rms_organoids.rds  # relative to app/
    n_orfs:         12847
    n_ribo_samples: 6
    n_rna_samples:  6
    prepared_on:    "2025-11-14"

  - study_id:       ews_primary
    display_name:   "EWS primary tumours"
    cancer_type:    EWS
    cohort:         Primary
    rds_path:       data/ews_primary/titan_ews_primary.rds
    n_orfs:         9234
    n_ribo_samples: 8
    n_rna_samples:  8
    prepared_on:    "2025-10-01"
```

> **Important:** `rds_path` is resolved relative to the `app/` directory at runtime. The app calls `readRDS(entry$rds_path)` from its working directory (`/srv/titan/` in Docker).

---

## Step 4: Bundle MS Peptides (Optional)

If a peptide file is placed at:
```
app/data/<study_id>/peptides_<study_id>.csv   (or .tsv or .txt)
```

TITAN will auto-detect and auto-load it when the study is loaded from the Study Library, without requiring the user to upload the file manually. A notification confirms the auto-load.

The file must contain at least one column of peptide sequences (≥ 8 aa). Standard column names auto-detected: `Peptide`, `Sequence`, `Annotated Sequence`, `Modified Sequence`. Additional columns (intensity, probability scores, etc.) are carried through to the MS metadata panel.

---

## RiboCrypt Sample Classification

Sample names in `paths.ribocrypt_ext` are classified as **primary tissue** or **cell-line** by prefix:

| Prefix | Class |
|---|---|
| `primary_` | Primary tissue |
| `hepatocyte_` | Primary tissue |
| `Myoblast_` | Primary tissue |
| `HSPC_` | Primary tissue |
| `huvec_` | Primary tissue |
| Any other | Cell-line |

If the external RiboCrypt matrix uses a different naming convention, the patterns must be updated in `scripts/prepare_titan_inputs.R` (`classify_ribocrypt_samples` function).

---

## Directory Layout After Setup

```
app/
├── data/
│   ├── catalog.yaml
│   ├── rms_organoids/
│   │   ├── titan_rms_organoids.rds      # produced by prepare_titan_inputs.R
│   │   └── peptides_rms_organoids.csv   # optional MS bundling
│   └── ews_primary/
│       ├── titan_ews_primary.rds
│       └── peptides_ews_primary.csv
└── ref/
    └── ensembl114_pep/
        ├── ensembl114_pep_index.rds
        ├── ensembl_gene_annotation.rds
        ├── ensembl114_pep.pdb           # BLAST db files
        ├── ensembl114_pep.phr
        ├── ensembl114_pep.pin
        ├── ensembl114_pep.pjs
        ├── ensembl114_pep.pot
        └── ensembl114_pep.psq
```

None of these files should be committed to git (all covered by `.gitignore`). In Docker, they are provided via volume mounts.
