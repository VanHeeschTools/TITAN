# TITAN Wiki

**Tumor Immunopeptidomics Target Atlas of Non-canonical ORFs**

TITAN is a Shiny application for prioritising non-canonical ORF (ncORF)-derived peptide candidates for cancer immunotherapy. It integrates study-specific ribo-seq and RNA-seq with publicly available databases (GTEx, TCGA, RiboCrypt) and runs offline canonical cross-reactivity and BLAST homology checks.

---

## Pages

| Page | Audience | Description |
|---|---|---|
| [App Guide](App-Guide) | End users | How to use each tab of the deployed app |
| [Data Preprocessing](Data-Preprocessing) | Bioinformaticians | How to prepare study data with `prepare_titan_inputs.R` |
| [Catalog Setup](Catalog-Setup) | Admins | Config YAML, catalog.yaml, and registering studies |
| [Reference Data](Reference-Data) | Admins | Building the Ensembl pep index and BLAST database |
| [Docker Deployment](Docker-Deployment) | DevOps / Admins | Building and running the TITAN Docker container |

---

## Quick Start

### For app users
Load a study on the **Data** tab → upload or auto-load MS peptides → click **EXPLORE TARGETS**.

### For admins deploying a new study
1. Run `scripts/prepare_titan_inputs.R <config.yaml>` to build the RDS.
2. Run `scripts/catalog/register_study.R <config.yaml> <rds_path>` to add it to the catalog.
3. Place `peptides_<study_id>.csv` in `app/data/<study_id>/` for auto-loading.
4. Restart the app (or rebuild the Docker image).

### For reference data setup (one-time)
```bash
# Run these once, from the app/ directory
Rscript scripts/database_prep/01_prep_ensembl_pep.R    # Ensembl pep index + BLAST db
Rscript scripts/database_prep/02_prep_allergen.R        # Allergen database
Rscript scripts/database_prep/03_prep_annotation.R      # Gene annotation table
```

---

## Repository Structure

```
titan/
├── app/
│   ├── app.R             # Main Shiny app
│   ├── global.R          # Libraries, reference data, constants
│   ├── R/                # Modular functions
│   │   ├── fct_data_loading.R
│   │   ├── fct_peptides.R
│   │   ├── fct_scoring.R
│   │   ├── fct_plots.R
│   │   ├── fct_table.R
│   │   ├── fct_utils.R
│   │   ├── report_html.R
│   │   └── report_plots.R
│   ├── www/              # CSS, JS, SVGs
│   ├── data/             # catalog.yaml + per-study data (NOT in git)
│   └── ref/              # Reference databases (NOT in git)
│       └── ensembl114_pep/
├── scripts/
│   ├── prepare_titan_inputs.R   # Study data preparation
│   ├── catalog/
│   │   └── register_study.R     # Catalog registration
│   ├── database_prep/
│   │   ├── 01_prep_ensembl_pep.R
│   │   ├── 02_prep_allergen.R
│   │   └── 03_prep_annotation.R
│   └── configs/                 # Per-study config YAMLs
├── Dockerfile
├── .dockerignore
└── .gitignore
```

---

## Data Security

Study data (patient-derived RDS objects, MS peptide files, catalog) **must never be committed to git** or included in the Docker image. All paths are provided to the container via volume mounts at runtime. See [Docker Deployment](Docker-Deployment) for details.
