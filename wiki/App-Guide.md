# App Guide

This page describes the TITAN Shiny interface for end users — researchers using an already-deployed instance of the app.

---

## Data Tab

The entry point. You must load both an ORF candidate dataset and an MS peptides file before analysis begins.

### Study Library

A searchable list of pre-configured studies. Click **Load** to read the corresponding RDS from disk. The currently loaded study shows an expanded card with an accent border and a trash icon to clear it.

You can filter the library by cancer type or cohort using the dropdowns, or type in the search box.

### Upload Data

Click **Upload Data** (toggle button) to load an arbitrary RDS that was prepared with `prepare_titan_inputs.R` but not registered in the catalog.

### MS Peptides

After loading a study, if a bundled `peptides_<study_id>.csv` exists in `app/data/<study_id>/`, it is loaded automatically and shown as a green "Ready" badge. Otherwise, upload a CSV or TSV with at least one column of peptide sequences (≥8 aa). TITAN auto-detects the peptide column from common names (`Peptide`, `Sequence`, `Annotated Sequence`, `Modified Sequence`); you can override this with the dropdown.

### Starting Analysis

When both sources are ready, a readiness summary appears with a count of ORFs × peptides. Click **EXPLORE TARGETS** to run the peptide-to-ORF matching and navigate to the Overview tab.

---

## Overview Tab

Summary statistics and four overview plots for the filtered candidate set.

### Stats row
| Stat | Meaning |
|---|---|
| Total candidates | ORFs passing the current filter (sidebar) |
| Unique genes | Distinct genes encoding those ORFs |
| Matched ORFs | ORFs with at least one MS peptide hit |
| MS peptide hits | Total peptide–ORF match rows |

### Sidebar — Filtering

| Control | Effect |
|---|---|
| PPM threshold | Minimum ribo-seq PPM per sample |
| Min. samples ≥ threshold | How many samples must meet the PPM threshold |
| TPM threshold | Minimum RNA-seq TPM per sample |
| Min. samples ≥ threshold | How many samples must meet the TPM threshold |
| ORF biotypes | Restrict to selected biotype classes |
| Reset filters | Restore all sliders/pickers to defaults |

Filters apply to all tabs simultaneously.

### Plots

**ORF biotype distribution** — stacked bar for all candidates and for peptide-matched ORFs only.

**Translation vs expression** — scatter of log₁₀(Median PPM+0.1) vs log₁₀(Median TPM+0.1), coloured by biotype. Peptide-matched ORFs appear as larger, fully-opaque points. Toggle "Matched only" to hide the background cloud.

**PPM distribution by biotype** — violin/box of log₁₀(Median PPM+0.1) split by biotype. When peptides are loaded, each biotype is shown as a split violin: grey left half = unmatched, coloured right half = matched.

---

## Prioritisation Tab

A ranked, filterable table of peptide-matched ORF candidates, grouped by gene.

### Scoring

Candidates are scored across 9 weighted dimensions (see the Scoring section below). Use the **Sidebar → Scoring** panel to adjust weights or select a preset.

**Presets:**

| Preset | Strategy |
|---|---|
| Cancer-specific | Maximises tumor-only candidates; heavily penalises any normal tissue signal |
| Pan-cancer | Broad coverage; tolerates tumor-enriched but not tumor-only candidates |

**Custom weights:** Each dimension slider accepts values from −1 (penalise) to +1 (reward), with 0 meaning "ignore". A bipolar slider means negative values are meaningful. The total weight display shows the sum of absolute weights; scores are normalised to 0–100.

#### Scoring Dimensions

| Dimension | Signal |
|---|---|
| % expressed (RNA-seq) | Fraction of tumor samples with TPM ≥ threshold |
| % translated (Ribo-seq) | Fraction of tumor samples with PPM ≥ threshold |
| Tumor specificity (GTEx) | tumor-only=1, enriched=0.5, non-specific=0 |
| Off-tissue risk (GTEx) | Safe=1, Acceptable=0.75, Borderline=0.5, Critical=0.25 |
| GTEx expression level | Normal tissue signal (use negative weight to penalise) |
| TCGA tumor coverage | % TCGA tumor samples expressing the gene |
| TCGA normal expression | % TCGA peritumoral samples (use negative weight to penalise) |
| RC primary tissue | Ribocrypt primary tissue translation (penalise if unwanted) |
| RC cell-line | Ribocrypt cell-line translation |

### Table Columns

| Column | Description |
|---|---|
| Gene | Gene symbol (click to open detail panel); expands to child ORF rows if multiple ORFs share the same peptide evidence |
| ORF-biotype | Coloured biotype badge |
| Peptides | Matched peptide sequences (expand with "…N more") |
| ORF-id | ORF hash identifier |
| Location | Genomic coordinates + strand + start codon |
| Specificity | GTEx tumor-only / enriched / non-specific badge |
| Off-tissue risk | Safety tier (Safe / Acceptable / Borderline / Critical) |
| Score | Visual score bar (0–100) |
| Transl. % / PPM | % tumor samples translated / median PPM |
| Expr. % / TPM | % tumor samples expressed / median TPM |
| TCGA columns | Tumor and normal coverage and expression |
| RC prim / CL | Ribocrypt primary tissue and cell-line metrics |

### Candidate Detail Panel

Click a gene name to open a slide-in panel below the table with:

- Per-sample expression boxplots: Target tumor · GTEx tissues · TCGA tumour vs normal
- Per-sample translation lollipop plots: Target tumor · Ribocrypt primary · Ribocrypt cell-line
- Score radar chart
- Dimension contribution bars

### Selecting for Report

Check the checkbox in the leftmost column (or use the header checkbox for the current page) to mark candidates. Selected candidates are exported with the **Export selection** button or included in the report.

### Filters

**Tumor specificity filter** and **Off-tissue risk filter** in the scoring sidebar further restrict which genes appear in the table.

### Downloads

| Button | Contents |
|---|---|
| Export parameters | CSV of current preset + weight values |
| Export ranked | CSV of all ranked candidates with all numeric columns |
| Export selection | CSV of checked candidates only |

---

## ORF Detail Tab

Per-ORF view with protein sequence and two safety checks. Navigate via the left sidebar.

### Gene / ORF Selector

1. Type a gene name in the **Gene** selectize dropdown.
2. If the gene has more than one matched ORF, a second **ORF** dropdown appears.
3. The counter "ORF N of M for GENE" shows position in the ranked list.
4. The metadata card below shows gene type, ORF biotype, coordinates (UCSC link), protein length, start codon, and ORF ID — all with copy buttons.
5. Gencode cross-matches and co-identified sibling ORFs (same peptide evidence) are listed below.

### Protein Sequence & Peptide Matches

The full protein sequence with matched peptides highlighted in the text. Hover a highlighted region to see MS metadata from the uploaded file.

### Cross-Reactivity (Canonical Proteome)

Automatically checked against the Ensembl 114 pep index whenever you navigate to an ORF.

- Uses `Biostrings::vmatchPattern` at 0, 1, and 2 mismatches.
- Results are collapsed to gene level (ENSG); isoform-level duplicates are merged.
- Self (the ORF's own gene) is excluded.
- If the peptide matches a canonical protein exactly → flagged as a warning.
- Click a row to open an expression modal showing Target / GTEx / TCGA expression for that canonical gene.

**Status messages:**

| Message | Meaning |
|---|---|
| ✓ No canonical matches (up to 2 mismatches) | Peptide(s) unique to the non-canonical proteome |
| ⚠ N exact gene(s), M 1-mismatch gene(s)… | Cross-reactive hits found — review before advancing the candidate |
| ✗ Ensembl 114 pep index not loaded | Reference data missing; run `scripts/database_prep/01_prep_ensembl_pep.R` |

### BLAST Homology

Run automatically 750 ms after you stop changing ORFs, using the full-protein sequence vs the Ensembl 114 pep BLAST database.

- Filter: ≥50% identity AND ≥30% alignment coverage.
- Results deduplicated to one row per gene (highest-identity hit kept).
- Click a row to open the expression modal with the pairwise alignment text.

**Status messages follow the same green/amber/red pattern as cross-reactivity.**

---

## Report Tab

Generates a self-contained HTML report for all checked candidates.

1. Check candidates in the Prioritisation table.
2. Navigate to Report; the green alert confirms how many are selected.
3. Choose Y-axis scale (log or raw).
4. Click **Generate report** to download.

Open the downloaded `.html` in any browser. Use **File → Print → Save as PDF** (A3 landscape, no margins, fit to page) for a clean printable version.

Each candidate page includes:
- Gene metadata, score, specificity and risk badges
- Per-sample expression plots (Target / GTEx / TCGA)
- Per-sample translation plots (Target / Ribocrypt)
- Score radar + dimension bars
- Protein sequence with highlighted peptides

---

## About Tab

Lists all data sources, scoring dimensions, and the offline reference pipeline scripts. The sidebar shows installed R package versions.

---

## ORF Biotype Legend

| Biotype | Description |
|---|---|
| ORF-annotated | Falls within an annotated CDS |
| Processed_transcript_ORF | ORF in a processed transcript |
| NC-variant | Non-coding variant ORF |
| uORF | Upstream ORF |
| uoORF | Upstream overlapping ORF |
| dORF | Downstream ORF |
| doORF | Downstream overlapping ORF |
| intORF | Intronic ORF |
| lncRNA-ORF | ORF in a lncRNA |
| pseudogene-ORF | ORF in a pseudogene |

---

## Off-Tissue Risk Tiers

Risk tiers are derived from GTEx expression in 28 normal adult tissues. The risk level of a gene reflects the **most critical** tissue where it is expressed above the GTEx Q3 threshold.

| Tier | Adult tissues in this category |
|---|---|
| Critical | Brain, Heart, Lung, Liver, Kidney, Colon, Esophagus, Stomach, Small Intestine, Adrenal Gland, Artery, Pancreas, Pituitary, Nerve, Whole Blood |
| Borderline | Bladder, Breast, Muscle, Ovary, Thyroid |
| Acceptable | Adipose, Cervix, Fallopian Tube, Prostate, Salivary, Skin, Spleen, Uterus, Vagina |
| Safe | No GTEx tissue above threshold |
| Unavailable | GTEx data not available for this study |

> Pediatric patients have a stricter risk table — several tissues classified as "Acceptable" in adults (e.g., Breast, Ovary, Cervix) are reclassified as "Critical" or "Borderline". This is handled internally by the scoring module.
