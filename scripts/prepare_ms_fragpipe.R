#!/usr/bin/env Rscript
# prepare_ms_fragpipe.R
#
# Converts a FragPipe PSM CSV to a TITAN-ready TSV.
#
# Steps:
#   1. Remove decoys / contaminants (if columns are present)
#   2. Apply Q-value cutoff (default 0.01, if column is present)
#   3. Deduplicate to one PSM per peptide (best score kept)
#   4. Write a flat TSV — upload directly to TITAN
#
# Usage:
#   Rscript prepare_ms_fragpipe.R --input psm.csv [--output titan_ms.tsv] [--qval 0.01]
#
# SBATCH example:
#   #!/bin/bash
#   #SBATCH --job-name=titan_prep
#   #SBATCH --cpus-per-task=4
#   #SBATCH --mem=16G
#   #SBATCH --time=00:30:00
#   singularity exec /hpc/local/Rocky8/pmc_vanheesch/singularity_images/<image>.sif \
#     Rscript /hpc/pmc_oatv/projects/tools_dev/titan/prepare_ms_fragpipe.R \
#     --input psm.csv --output titan_ms.tsv

suppressPackageStartupMessages({
  library(data.table)
})

# ── Parse arguments ────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[idx + 1L]
}

input_file  <- get_arg("--input")
output_file <- get_arg("--output")
qval_cutoff <- as.numeric(get_arg("--qval", "0.01"))

if (is.null(input_file)) {
  cat("Usage: Rscript prepare_ms_fragpipe.R --input psm.csv [--output titan_ms.tsv] [--qval 0.01]\n")
  quit(status = 1)
}

if (is.null(output_file)) {
  output_file <- sub("\\.[^.]+$", "_titan.tsv", input_file)
  if (output_file == input_file) output_file <- paste0(input_file, "_titan.tsv")
}

# ── Read file ──────────────────────────────────────────────────────────────────
cat(sprintf("Reading: %s\n", input_file))
ext <- tolower(tools::file_ext(input_file))
sep <- if (ext == "csv") "," else "\t"
psm <- fread(input_file, sep = sep, header = TRUE, data.table = FALSE,
             stringsAsFactors = FALSE, check.names = FALSE)
# Strip UTF-8 BOM that FragPipe sometimes writes into the first column name
colnames(psm)[1] <- gsub("^\xef\xbb\xbf", "", colnames(psm)[1], useBytes = TRUE)
cat(sprintf("  %s PSMs, %s columns\n", formatC(nrow(psm), big.mark = ","), ncols <- ncol(psm)))
cat(sprintf("  Columns: %s\n", paste(head(colnames(psm), 20), collapse = ", "),
            if (ncols > 20) "…" else ""))

# ── Identify key columns (case-insensitive) ────────────────────────────────────
find_col <- function(candidates, df) {
  candidates[candidates %in% colnames(df)][1]
}

pep_col       <- find_col(c("Peptide", "peptide", "PEPTIDE", "Sequence", "Modified Sequence"), psm)
prob_col      <- find_col(c("Probability", "PEP.score", "Score"), psm)
qval_col      <- find_col(c("Qvalue", "Q-value", "q-value", "FDR"), psm)
decoy_col     <- find_col(c("Is Decoy", "IsDecoy", "Decoy"), psm)
contam_col    <- find_col(c("Is Contaminant", "IsContaminant", "Contaminant"), psm)

if (is.na(pep_col)) stop("Cannot find a peptide sequence column. Check column names.")
cat(sprintf("  Peptide column   : %s\n", pep_col))
cat(sprintf("  Quality column   : %s\n", if (is.na(prob_col))  "(none found)" else prob_col))
cat(sprintf("  Q-value column   : %s\n", if (is.na(qval_col))  "(none found)" else qval_col))
cat(sprintf("  Decoy column     : %s\n", if (is.na(decoy_col)) "(none found)" else decoy_col))
cat(sprintf("  Contaminant col  : %s\n", if (is.na(contam_col))"(none found)" else contam_col))

n_start <- nrow(psm)

# ── Filter decoys ──────────────────────────────────────────────────────────────
if (!is.na(decoy_col)) {
  before <- nrow(psm)
  val    <- psm[[decoy_col]]
  keep   <- if (is.logical(val)) !val else !(tolower(as.character(val)) %in% c("true", "1", "yes"))
  psm    <- psm[keep, , drop = FALSE]
  cat(sprintf("  Removed %s decoys → %s remain\n",
              formatC(before - nrow(psm), big.mark = ","),
              formatC(nrow(psm), big.mark = ",")))
} else {
  cat("  [skip] No decoy column found — not filtering\n")
}

# ── Filter contaminants ────────────────────────────────────────────────────────
if (!is.na(contam_col)) {
  before <- nrow(psm)
  val    <- psm[[contam_col]]
  keep   <- if (is.logical(val)) !val else !(tolower(as.character(val)) %in% c("true", "1", "yes"))
  psm    <- psm[keep, , drop = FALSE]
  cat(sprintf("  Removed %s contaminants → %s remain\n",
              formatC(before - nrow(psm), big.mark = ","),
              formatC(nrow(psm), big.mark = ",")))
} else {
  cat("  [skip] No contaminant column found — not filtering\n")
}

# ── Filter by Q-value ──────────────────────────────────────────────────────────
if (!is.na(qval_col)) {
  before <- nrow(psm)
  psm    <- psm[as.numeric(psm[[qval_col]]) <= qval_cutoff, , drop = FALSE]
  cat(sprintf("  Q-value ≤ %s: removed %s → %s remain\n",
              qval_cutoff,
              formatC(before - nrow(psm), big.mark = ","),
              formatC(nrow(psm), big.mark = ",")))
} else {
  cat(sprintf("  [skip] No Q-value column found — not applying cutoff of %s\n", qval_cutoff))
}

# ── Deduplicate: best PSM per peptide ─────────────────────────────────────────
n_before_dedup <- nrow(psm)
n_unique_peps  <- length(unique(psm[[pep_col]]))

if (n_unique_peps < n_before_dedup) {
  # Pick ordering: Probability (high good), else Qvalue (low good), else first row
  if (!is.na(prob_col)) {
    psm <- psm[order(psm[[pep_col]], -as.numeric(psm[[prob_col]])), , drop = FALSE]
    cat(sprintf("  Deduplicating by highest %s\n", prob_col))
  } else if (!is.na(qval_col)) {
    psm <- psm[order(psm[[pep_col]], as.numeric(psm[[qval_col]])), , drop = FALSE]
    cat(sprintf("  Deduplicating by lowest %s\n", qval_col))
  } else {
    psm <- psm[order(psm[[pep_col]]), , drop = FALSE]
    cat("  Deduplicating by first occurrence (no quality column found)\n")
  }
  psm <- psm[!duplicated(psm[[pep_col]]), , drop = FALSE]
  cat(sprintf("  Deduplicated: %s → %s unique peptides\n",
              formatC(n_before_dedup, big.mark = ","),
              formatC(nrow(psm), big.mark = ",")))
} else {
  cat(sprintf("  Already unique: %s peptides\n", formatC(nrow(psm), big.mark = ",")))
}

# ── Drop columns not used in TITAN ────────────────────────────────────────────
# Technical / redundant columns that add no value in the app.
# Any column NOT in this list that is present in the file is kept by default,
# so novel columns from other FragPipe versions are preserved automatically.
DROP_COLS <- c(
  # Spectrum / file identifiers
  "Spectrum", "Spectrum File",
  # Redundant sequence context (Prev AA / Next AA already capture this)
  "Extended Peptide",
  # Detailed mass accuracy (not displayed in app)
  "Observed Mass", "Calibrated Observed Mass",
  "Observed M/Z", "Calibrated Observed M/Z",
  "Calculated Peptide Mass", "Calculated M/Z", "Delta Mass",
  # Secondary / internal scores
  "SpectralSim", "RTScore", "IMScore", "Nextscore",
  # FragPipe protein coordinates (TITAN re-maps peptides from scratch)
  "Protein Start", "Protein End",
  # Boolean filter flags (always FALSE after filtering)
  "Is Decoy", "Is Contaminant", "Is Unique",
  # Redundant protein identifiers (Protein column already has sp|ACC|NAME)
  "Protein ID", "Entry Name",
  # Ion mobility detail columns
  "Ion Mobility", "Purity",
  "Apex Ion Mobility", "Ion Mobility Start", "Ion Mobility End", "Ion Mobility FWHM",
  # Scan-level metadata
  "Parent Scan Number", "Apex Scan Number", "Apex Retention Time",
  "Retention Time Start", "Retention Time End", "Retention Time FWHM", "Traced Scans",
  # Enzymatic detail (not shown in app)
  "Number of Enzymatic Termini"
)

to_drop  <- intersect(DROP_COLS, colnames(psm))
to_keep  <- setdiff(colnames(psm), to_drop)
psm      <- psm[, to_keep, drop = FALSE]
cat(sprintf("  Dropped %d technical columns → %d kept\n", length(to_drop), ncol(psm)))

# ── Write output ───────────────────────────────────────────────────────────────
cat(sprintf("Writing: %s\n", output_file))
write.table(psm, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Done. %s peptides × %s columns\n",
            formatC(nrow(psm), big.mark = ","), ncol(psm)))
cat(sprintf("\nNext step: upload '%s' to TITAN and select '%s' as the peptide column.\n",
            basename(output_file), pep_col))

