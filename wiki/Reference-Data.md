# Reference Data

TITAN requires three reference datasets to run the ORF Detail safety checks. These are built once and stored in `app/ref/`. They are **not** included in the git repository or Docker image — they are provided via volume mounts at runtime.

---

## Overview

| File | Script | Size (approx.) | Purpose |
|---|---|---|---|
| `ref/ensembl114_pep/ensembl114_pep_index.rds` | `01_prep_ensembl_pep.R` | ~300 MB | Exact/near-exact cross-reactivity matching |
| `ref/ensembl114_pep/ensembl114_pep.*` (BLAST db) | `01_prep_ensembl_pep.R` | ~150 MB | BLAST homology search |
| `ref/ensembl114_pep/ensembl_gene_annotation.rds` | `03_prep_annotation.R` | ~5 MB | Gene symbol + description for BLAST hits |
| `ref/gencode_orfs_phase2.csv` | External download | ~13 MB | GENCODE TransCode Phase 2 ncORF cross-matching |

All scripts are run from the `app/` directory.

---

## Script 01 — Ensembl 114 Pep Index + BLAST Database

**File:** `scripts/database_prep/01_prep_ensembl_pep.R`

### What it does

1. Downloads the Ensembl release 114 human proteome FASTA from the Ensembl FTP.
2. Parses all ~118,000 protein sequences (ENSP → ENSG → gene symbol mapping).
3. Deduplicates by MD5 of amino acid sequence: multiple ENSP isoforms with identical sequences are collapsed into one representative entry (`ensembl114_pep_dedup.fa`).
4. Builds a lookup index RDS (`ensembl114_pep_index.rds`) for fast cross-reactivity queries.
5. Runs `makeblastdb` to build the blastp database from the deduplicated FASTA.

### Outputs

```
app/ref/ensembl114_pep/
├── Homo_sapiens.GRCh38.pep.all.fa         # full FASTA (~60 MB, gitignored)
├── Homo_sapiens.GRCh38.pep.all.fa.gz      # compressed original (gitignored)
├── ensembl114_pep_dedup.fa                 # deduplicated sequences (~30 MB, gitignored)
├── ensembl114_pep_index.rds                # in-app index (gitignored)
└── ensembl114_pep.*                        # BLAST database files (gitignored)
```

### Index structure (`ensembl114_pep_index.rds`)

An R list with:

| Element | Type | Description |
|---|---|---|
| `$seqs` | named character | MD5 → AA sequence (60–80 K unique sequences) |
| `$md5_to_ensg` | named list | MD5 → character vector of ENSG IDs |
| `$md5_to_sym` | named list | MD5 → character vector of gene symbols |
| `$md5_to_ensp` | named list | MD5 → character vector of ENSP IDs |
| `$ensp_to_md5` | named character | ENSP → MD5 (for BLAST sseqid back-lookup) |
| `$rep_table` | data.frame | One row per unique sequence (md5, representative ENSP/ENSG, etc.) |

### How cross-reactivity uses this index

For each matched peptide against the selected ORF, TITAN calls:
```r
Biostrings::vmatchPattern(pep_aa, AAStringSet(index$seqs), max.mismatch = 0|1|2)
```
Results at each mismatch level are back-mapped via `md5_to_ensg` and `md5_to_sym`, then collapsed to one row per gene (keeping the lowest mismatch level). The ORF's own gene is excluded from results.

### BLAST database

Built from `ensembl114_pep_dedup.fa` using:
```bash
makeblastdb -in ensembl114_pep_dedup.fa -dbtype prot \
            -out ensembl114_pep \
            -title "Ensembl114_pep_dedup"
```

The app uses `rBLAST::blast(db = "ref/ensembl114_pep/ensembl114_pep", type = "blastp")` and then filters hits at ≥50% identity AND ≥30% alignment coverage.

### Running the script

Requires internet access and `makeblastdb` in PATH (provided by the bundled NCBI BLAST+ binary at `app/bin/`).

```bash
# From the app/ directory
Rscript ../scripts/database_prep/01_prep_ensembl_pep.R
```

Or via SLURM from the project root:
```bash
sbatch --mem=16G \
  --output=logs/prep_ensembl_%j.out \
  --wrap="cd app && Rscript ../scripts/database_prep/01_prep_ensembl_pep.R"
```

**Runtime:** ~10–15 min (download + parse + dedup + makeblastdb).

---

## Script 02 — Allergen Database

**File:** `scripts/database_prep/02_prep_allergen.R`

> Details for this script are TBD. It prepares an allergenicity reference database used for an optional safety check. Run it from the `app/` directory in the same way as script 01.

---

## Script 03 — Gene Annotation Table

**File:** `scripts/database_prep/03_prep_annotation.R`

### What it does

Fetches offline gene-level annotation from Ensembl 114 biomaRt (or from a cached local biomaRt dump) and saves a data.frame to `ref/ensembl114_pep/ensembl_gene_annotation.rds`.

### Output columns

| Column | Description |
|---|---|
| `ensembl_gene_id` | ENSG ID (no version suffix) |
| `external_gene_name` | Gene symbol |
| `uniprotswissprot` | UniProt/SwissProt accession |
| `description` | Gene description string |

### Used by

BLAST hit annotation: when a hit ENSP is back-mapped to an ENSG, TITAN looks up the gene symbol and description in this table for display in the BLAST hits table and expression modal.

### Running the script

Requires biomaRt access (internet or local Ensembl mirror). Run from `app/`:
```bash
Rscript ../scripts/database_prep/03_prep_annotation.R
```

Or via SLURM:
```bash
sbatch --mem=8G \
  --output=logs/prep_annotation_%j.out \
  --wrap="cd app && Rscript ../scripts/database_prep/03_prep_annotation.R"
```

---

## GENCODE ORFs — TransCode Phase 2

**File:** `app/ref/gencode_orfs_phase2.csv`

A reference table of GENCODE-annotated ncORFs from the TransCode Phase 2 study. Used for cross-matching: if an MS peptide also matches a TransCode ORF (but not the study's own in-house ncORF table), TITAN flags it as a "Gencode cross-match" in the ORF metadata sidebar.

### Required columns

| Column in CSV | Mapped to app column |
|---|---|
| `releasev45_id` | `orf_id` |
| `gene_id` | `gene_id` + `gene_id_clean` (version stripped) |
| `gene_name` | `gene_name` |
| `gene_biotype` | `gene_biotype` |
| `orf_type` | `orf_biotype_single` (PT → Processed_transcript_ORF) |
| `sequence_aa` | `protein_seq` + `protein_length` |
| `chrm` | `chr` |
| `starts (0-based)` | `orf_start` |
| `ends (0-based)` | `orf_end` |
| `strand` | `strand` |
| `initiation_codon` | `start_codon` |

### Source

> _(Download URL, DOI, or internal path to be filled in by the team.)_

### Behaviour if missing

If `app/ref/gencode_orfs_phase2.csv` is absent, TITAN loads normally but silently disables Gencode cross-matching. No error is shown.

---

## BLAST Binary (`app/bin/`)

The app ships `blastp.REAL` — the `blastp` executable from NCBI BLAST+ 2.14.0 compiled for Linux x86-64. At Docker build time, a symlink `app/bin/blastp → app/bin/blastp.REAL` is created so rBLAST can find it.

The tarball `app/bin/ncbi-blast-2.14.0+-x64-linux.tar.gz` and expanded directory `app/bin/ncbi-blast-2.14.0+/` are gitignored and excluded from the Docker build context.

On the HPC, the NCBI BLAST+ binaries can also be accessed via the Singularity image:
```
/hpc/local/Rocky8/pmc_vanheesch/singularity_images/ncbi_blast-2.14.0.sif
```

---

## Complete Reference Prep — Run Order

```bash
# From the project root; all scripts expect CWD = app/

# 1. Ensembl pep index + BLAST database (needs internet + ~15 min)
cd app && Rscript ../scripts/database_prep/01_prep_ensembl_pep.R

# 2. Allergen database (details TBD)
Rscript ../scripts/database_prep/02_prep_allergen.R

# 3. Gene annotation (needs biomaRt internet access)
Rscript ../scripts/database_prep/03_prep_annotation.R

# 4. Copy GENCODE Phase 2 CSV to ref/
cp /path/to/gencode_orfs_phase2.csv app/ref/gencode_orfs_phase2.csv
```

After these steps, `app/ref/` should contain:
```
ref/
├── ensembl114_pep/
│   ├── ensembl114_pep_index.rds        ✓
│   ├── ensembl_gene_annotation.rds     ✓
│   ├── ensembl114_pep.pdb              ✓
│   ├── ensembl114_pep.phr              ✓
│   ├── ensembl114_pep.pin              ✓
│   ├── ensembl114_pep.pjs              ✓
│   ├── ensembl114_pep.pot              ✓
│   └── ensembl114_pep.psq             ✓
└── gencode_orfs_phase2.csv             ✓
```
