#!/usr/bin/env Rscript
#
# TITAN — Tumour Immunopeptidomics Target ANnotation
# Data preparation script: computes all per-ORF metrics and saves the integrated
# table for the Shiny app.
#
# Usage:
#   Rscript prepare_titan_inputs.R <config.yaml>
#   Rscript prepare_titan_inputs.R configs/tnbc_navarra.yaml
#
# Outputs (paths resolved from config):
#   titan_<study_id>.rds  — full table + per-sample matrices (app input)
#   titan_<study_id>.csv  — flat CSV for inspection

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(matrixStats)
  library(yaml)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L)
  stop("Usage: Rscript prepare_titan_inputs.R <config.yaml>", call. = FALSE)

config_path <- args[1L]
if (!file.exists(config_path))
  stop(sprintf("Config file not found: %s", config_path), call. = FALSE)

cfg <- yaml::read_yaml(config_path)

# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION  (fail fast before any heavy data loading)
# ─────────────────────────────────────────────────────────────────────────────

required_top <- c("study_id", "display_name", "cancer_type", "cohort",
                  "tumor_type", "paths", "output")
missing_top <- setdiff(required_top, names(cfg))
if (length(missing_top))
  stop(sprintf("Config missing required top-level fields: %s",
               paste(missing_top, collapse = ", ")), call. = FALSE)

if (is.null(cfg$output$dir))
  stop("Config missing: output.dir", call. = FALSE)

required_path_keys <- c("ncorfs", "ribocrypt_ext",
                        "gtex_quant", "gtex_coldata", "tcga_quant", "tcga_coldata")
missing_keys <- setdiff(required_path_keys, names(cfg$paths))
if (length(missing_keys))
  stop(sprintf("Config missing required paths: %s",
               paste(missing_keys, collapse = ", ")), call. = FALSE)

path_errors <- Filter(nchar, sapply(required_path_keys, function(p) {
  f <- cfg$paths[[p]]
  if (!file.exists(f)) sprintf("  paths.%s: %s", p, f) else ""
}))

# de_sig_all is optional; if provided, the path must exist
has_de_sig <- !is.null(cfg$paths$de_sig_all)
if (has_de_sig && !file.exists(cfg$paths$de_sig_all)) {
  path_errors <- c(path_errors,
    sprintf("  paths.de_sig_all (specified but not found): %s", cfg$paths$de_sig_all))
  has_de_sig <- FALSE
}

# ribo_ppm/ribo_psites (internal ribo-seq) are optional together; if either is
# provided, both must be provided and both paths must exist. Studies with no
# internal ribo-seq (translation evidence from RiboCrypt instead) omit both.
has_riboseq <- !is.null(cfg$paths$ribo_ppm) || !is.null(cfg$paths$ribo_psites)
if (has_riboseq) {
  missing_ribo <- setdiff(c("ribo_ppm", "ribo_psites"), names(cfg$paths))
  if (length(missing_ribo))
    stop(sprintf(
      "paths.ribo_ppm/ribo_psites: provide both or neither (missing: %s)",
      paste(missing_ribo, collapse = ", ")), call. = FALSE)
  for (p in c("ribo_ppm", "ribo_psites")) {
    if (!file.exists(cfg$paths[[p]]))
      path_errors <- c(path_errors,
        sprintf("  paths.%s (specified but not found): %s", p, cfg$paths[[p]]))
  }
}

# tumor_quant is optional; if provided, the path must exist
if (!is.null(cfg$paths$tumor_quant) && !file.exists(cfg$paths$tumor_quant))
  path_errors <- c(path_errors,
    sprintf("  paths.tumor_quant (specified but not found): %s", cfg$paths$tumor_quant))

if (length(path_errors))
  stop(sprintf("The following configured paths do not exist on disk:\n%s",
               paste(path_errors, collapse = "\n")), call. = FALSE)

if (!is.null(cfg$condition) && !is.null(cfg$condition$pattern)) {
  missing_rc <- setdiff(c("match_label", "nomatch_label"), names(cfg$condition))
  if (length(missing_rc))
    stop(sprintf("condition.pattern is set but missing: %s",
                 paste(missing_rc, collapse = ", ")), call. = FALSE)
}

for (th_name in c("expression", "gtex_q3")) {
  v <- cfg$thresholds[[th_name]]
  if (!is.null(v) && (!is.numeric(v) || v <= 0))
    stop(sprintf("thresholds.%s must be a positive number (got: %s)", th_name, v), call. = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# RESOLVE CONFIG VALUES
# ─────────────────────────────────────────────────────────────────────────────

TARGET_TUMOR_TYPE <- cfg$tumor_type
target_label      <- cfg$target_label %||% cfg$tumor_type
expr_threshold    <- cfg$thresholds$expression %||% 1
gtex_q3_threshold <- cfg$thresholds$gtex_q3   %||% 1

output_rds <- file.path(cfg$output$dir,
                        cfg$output$rds %||% paste0("titan_", cfg$study_id, ".rds"))
output_csv <- file.path(cfg$output$dir,
                        cfg$output$csv %||% paste0("titan_", cfg$study_id, ".csv"))

dir.create(cfg$output$dir, recursive = TRUE, showWarnings = FALSE)

if (!has_de_sig)
  warning(sprintf(
    paste0("de_sig_all not configured for study '%s'.\n",
           "  GTEX_DE_sig_in_all, GTEX_tumor_only, and GTEX_tumor_enriched\n",
           "  will be NA for ALL candidates in this study's output."),
    cfg$study_id), call. = FALSE)

if (!has_riboseq)
  warning(sprintf(
    paste0("ribo_ppm/ribo_psites not configured for study '%s'.\n",
           "  target_translation_* columns will be NA for ALL candidates in this\n",
           "  study's output (RiboCrypt data, if configured, is unaffected)."),
    cfg$study_id), call. = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

strip_ensg_version <- function(x) sub("\\..*", "", x)

# Computes num_samples, pct_samples, median, max from a gene/ORF × sample matrix.
compute_expression_metrics <- function(mat, threshold = 1) {
  mat <- as.matrix(mat)
  data.frame(
    num_samples  = rowSums(mat >= threshold, na.rm = TRUE),
    pct_samples  = 100 * rowMeans(mat >= threshold, na.rm = TRUE),
    median_value = apply(mat, 1, median, na.rm = TRUE),
    max_value    = apply(mat, 1, max, na.rm = TRUE),
    row.names    = rownames(mat)
  )
}

# Classifies external RiboCrypt sample names as primary tissue or cell line.
# These prefixes reflect a fixed convention of the external RiboCrypt database.
classify_ribocrypt_samples <- function(sample_names) {
  primary_patterns <- c("^primary_", "^hepatocyte_", "^Myoblast_", "^HSPC_", "^huvec_")
  is_primary <- sapply(sample_names, function(s) {
    any(sapply(primary_patterns, function(p) grepl(p, s)))
  })
  list(primary = sample_names[is_primary], cell_line = sample_names[!is_primary])
}

# Infers sample condition from sample name using the study config pattern.
# Applied to both RNA-seq and ribo-seq target tumor samples.
get_sample_condition <- function(sample_names) {
  rc <- cfg$condition
  if (is.null(rc) || is.null(rc$pattern))
    return(rep(TARGET_TUMOR_TYPE, length(sample_names)))
  ifelse(grepl(rc$pattern, sample_names, ignore.case = TRUE),
         rc$match_label, rc$nomatch_label)
}

# Loads an expression matrix from an RDS txi object or a flat CSV/TSV file.
# Returns a numeric matrix with genes as rows and samples as columns.
# For .rds: expects a txi-style list with an $abundance element (matrix).
# For .csv/.tsv: expects gene IDs in the first column (used as rownames).
load_expression_data <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!is.list(obj) || is.null(obj$abundance) || !is.matrix(obj$abundance))
      stop(sprintf(
        "RDS file does not contain a txi-style list with a matrix $abundance element:\n  %s",
        path), call. = FALSE)
    obj$abundance
  } else if (ext %in% c("csv", "tsv")) {
    sep <- if (ext == "csv") "," else "\t"
    df  <- read.delim(path, sep = sep, check.names = FALSE, row.names = 1)
    as.matrix(df)
  } else {
    stop(sprintf(
      "Unrecognized expression file extension '.%s' — expected .rds, .csv, or .tsv:\n  %s",
      ext, path), call. = FALSE)
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 — ORF BACKBONE
# ─────────────────────────────────────────────────────────────────────────────

cat("=== TITAN input preparation ===\n")
cat(sprintf("    Study : %s (%s)\n", cfg$display_name, cfg$study_id))
cat(sprintf("    Config: %s\n\n", config_path))
cat("[1/6] Loading ncORF candidate table...\n")

ncorfs <- fread(cfg$paths$ncorfs, data.table = FALSE)
ncorfs$gene_id_clean <- strip_ensg_version(ncorfs$gene_id)

cat(sprintf("      %d ncORF candidates, %d unique genes\n",
            nrow(ncorfs), n_distinct(ncorfs$gene_id_clean)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2 — TARGET TRANSLATION (ribo-seq)
# ─────────────────────────────────────────────────────────────────────────────

candidate_ids <- ncorfs$orf_id

if (has_riboseq) {
  cat(sprintf("[2/6] Computing target translation metrics (%s ribo-seq)...\n", target_label))

  ribo_ppm    <- fread(cfg$paths$ribo_ppm,    data.table = FALSE)
  ribo_psites <- fread(cfg$paths$ribo_psites, data.table = FALSE)

  rownames(ribo_ppm)    <- ribo_ppm$orf_id;    ribo_ppm$orf_id    <- NULL
  rownames(ribo_psites) <- ribo_psites$orf_id; ribo_psites$orf_id <- NULL

  # Shorten sample IDs for display using study-specific regex (null = skip)
  if (!is.null(cfg$sample_id_regex)) {
    ribo_ppm    <- rename_with(ribo_ppm,    ~ sub(cfg$sample_id_regex, "", .x))
    ribo_psites <- rename_with(ribo_psites, ~ sub(cfg$sample_id_regex, "", .x))
  }

  common_ribo <- intersect(candidate_ids, rownames(ribo_ppm))

  ribo_ppm_mat    <- as.matrix(ribo_ppm[common_ribo, ])
  ribo_psites_mat <- as.matrix(ribo_psites[common_ribo, ])

  transl_metrics <- compute_expression_metrics(ribo_ppm_mat, threshold = expr_threshold) %>%
    rename(
      target_translation_num_samples = num_samples,
      target_translation_pct_samples = pct_samples,
      target_translation_median_PPM  = median_value,
      target_translation_max_PPM     = max_value
    ) %>%
    mutate(
      target_translation_median_psites = apply(ribo_psites_mat, 1, median, na.rm = TRUE),
      orf_id = rownames(.)
    )

  ribo_sample_meta <- data.frame(
    sample_id = colnames(ribo_ppm_mat),
    condition = get_sample_condition(colnames(ribo_ppm_mat))
  )

  cat(sprintf("      %d ORFs with ribo-seq data, %d samples\n",
              nrow(transl_metrics), nrow(ribo_sample_meta)))
} else {
  cat("[2/6] No internal ribo-seq configured for this study — target_translation_* set to NA...\n")

  ribo_ppm_mat <- matrix(nrow = length(candidate_ids), ncol = 0,
                         dimnames = list(candidate_ids, NULL))

  transl_metrics <- data.frame(
    orf_id                            = candidate_ids,
    target_translation_num_samples   = NA_integer_,
    target_translation_pct_samples   = NA_real_,
    target_translation_median_PPM    = NA_real_,
    target_translation_max_PPM       = NA_real_,
    target_translation_median_psites = NA_real_,
    stringsAsFactors = FALSE
  )

  ribo_sample_meta <- data.frame(
    sample_id = character(0),
    condition = character(0),
    stringsAsFactors = FALSE
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 — RIBOCRYPT EXTERNAL (primary tissue / cell line)
# ─────────────────────────────────────────────────────────────────────────────

cat("[3/6] Computing RiboCrypt external metrics...\n")

ribocrypt_ext <- fread(cfg$paths$ribocrypt_ext, data.table = FALSE)
rownames(ribocrypt_ext) <- ribocrypt_ext$orf_id; ribocrypt_ext$orf_id <- NULL

sample_classes <- classify_ribocrypt_samples(colnames(ribocrypt_ext))
cat(sprintf("      %d primary tissue samples: %s...\n",
            length(sample_classes$primary),
            paste(head(sample_classes$primary, 3), collapse = ", ")))
cat(sprintf("      %d cell-line samples: %s...\n",
            length(sample_classes$cell_line),
            paste(head(sample_classes$cell_line, 3), collapse = ", ")))

common_rc <- intersect(candidate_ids, rownames(ribocrypt_ext))
rc_mat    <- as.matrix(ribocrypt_ext[common_rc, ])

primary_mat   <- rc_mat[, sample_classes$primary,   drop = FALSE]
cell_line_mat <- rc_mat[, sample_classes$cell_line, drop = FALSE]

primary_metrics <- compute_expression_metrics(primary_mat, threshold = expr_threshold) %>%
  rename(
    ribocrypt_primary_num_samples = num_samples,
    ribocrypt_primary_pct_samples = pct_samples,
    ribocrypt_primary_median_PPM  = median_value,
    ribocrypt_primary_max_PPM     = max_value
  ) %>%
  mutate(orf_id = rownames(.))

cell_line_metrics <- compute_expression_metrics(cell_line_mat, threshold = expr_threshold) %>%
  rename(
    `ribocrypt_cell-line_num_samples` = num_samples,
    `ribocrypt_cell-line_pct_samples` = pct_samples,
    `ribocrypt_cell-line_median_PPM`  = median_value,
    `ribocrypt_cell-line_max_PPM`     = max_value
  ) %>%
  mutate(orf_id = rownames(.))

ribocrypt_mat         <- rc_mat
ribocrypt_sample_meta <- data.frame(
  sample_id = colnames(ribocrypt_mat),
  group     = ifelse(colnames(ribocrypt_mat) %in% sample_classes$primary,
                     "Primary", "Cell-line"),
  stringsAsFactors = FALSE
)
rm(ribocrypt_ext, primary_mat, cell_line_mat); gc()
cat("      Done.\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4 — TARGET EXPRESSION + GTEX METRICS
#
# paths.gtex_quant may be either:
#   (a) A combined object containing both tumour (tissue_type == TARGET_TUMOR_TYPE)
#       and GTEx normal samples — the default when paths.tumor_quant is absent.
#   (b) A GTEx-only object, paired with a separate paths.tumor_quant file that
#       contains only tumour samples (all columns used, no coldata needed).
# ─────────────────────────────────────────────────────────────────────────────

has_separate_tumor <- !is.null(cfg$paths$tumor_quant)

cat("[4/6] Loading GTEx and tumour quantification data...\n")

coldata_gtex <- readRDS(cfg$paths$gtex_coldata)
stopifnot("tissue_type" %in% colnames(coldata_gtex))

if (has_separate_tumor) {
  # (b) Standalone GTEx — all coldata rows with a non-NA tissue_type are GTEx normals.
  # Intersect with matrix columns: the coldata may be a filtered subset of a larger
  # quantification object, so not all coldata sample IDs are guaranteed to be present.
  gtex_full_m           <- load_expression_data(cfg$paths$gtex_quant)
  rownames(gtex_full_m) <- strip_ensg_version(rownames(gtex_full_m))
  gtex_mask             <- !is.na(coldata_gtex$tissue_type)
  gtex_ids              <- intersect(coldata_gtex$sample_id[gtex_mask], colnames(gtex_full_m))
  gtex_tpm_m            <- gtex_full_m[, gtex_ids, drop = FALSE]
  rm(gtex_full_m)

  tumor_tpm_m           <- load_expression_data(cfg$paths$tumor_quant)
  rownames(tumor_tpm_m) <- strip_ensg_version(rownames(tumor_tpm_m))
  tumor_ids             <- colnames(tumor_tpm_m)

  cat(sprintf("      Tumour samples (%s): %d (separate quantification)\n",
              TARGET_TUMOR_TYPE, length(tumor_ids)))
  cat(sprintf("      GTEx normal samples: %d across %d tissues\n",
              length(gtex_ids),
              n_distinct(coldata_gtex$tissue_type[gtex_mask])))
} else {
  # (a) Combined — split by tissue_type in coldata
  combined_m           <- load_expression_data(cfg$paths$gtex_quant)
  rownames(combined_m) <- strip_ensg_version(rownames(combined_m))

  tumor_mask <- coldata_gtex$tissue_type %in% TARGET_TUMOR_TYPE
  gtex_mask  <- coldata_gtex$tissue_type != TARGET_TUMOR_TYPE & !is.na(coldata_gtex$tissue_type)

  cat(sprintf("      Tumour samples (%s): %d\n", TARGET_TUMOR_TYPE, sum(tumor_mask)))
  cat(sprintf("      GTEx normal samples: %d across %d tissues\n",
              sum(gtex_mask),
              n_distinct(coldata_gtex$tissue_type[gtex_mask])))

  tumor_ids   <- coldata_gtex$sample_id[tumor_mask]
  tumor_tpm_m <- combined_m[, tumor_ids, drop = FALSE]
  gtex_ids    <- coldata_gtex$sample_id[gtex_mask]
  gtex_tpm_m  <- combined_m[, gtex_ids, drop = FALSE]
  rm(combined_m)
}

# ── 4a. Target expression (tumour samples) ───────────────────────────────────

rna_tpm_sub <- tumor_tpm_m

expr_metrics <- compute_expression_metrics(rna_tpm_sub, threshold = expr_threshold) %>%
  dplyr::rename(
    target_expression_num_samples = num_samples,
    target_expression_pct_samples = pct_samples,
    target_expression_median_TPM  = median_value,
    target_expression_max_TPM     = max_value
  ) %>%
  mutate(gene_id_clean = rownames(.))

rna_sample_meta <- data.frame(
  sample_id   = tumor_ids,
  tissue_type = TARGET_TUMOR_TYPE,
  condition   = get_sample_condition(tumor_ids)
)

cat(sprintf("      %d genes with target expression data (%d samples)\n",
            nrow(expr_metrics), length(tumor_ids)))

# ── 4b. GTEx DE classification (from pre-computed DE table) ──────────────────

if (has_de_sig) {
  de_sig <- read.delim(cfg$paths$de_sig_all, check.names = FALSE)
  de_sig$gene_id_clean <- strip_ensg_version(de_sig$gene_id)
  gtex_de_metrics <- de_sig %>%
    transmute(
      gene_id_clean,
      GTEX_DE_sig_in_all  = sig_in_all,
      GTEX_tumor_only     = low_all_tissues,
      GTEX_tumor_enriched = Q3_GTEx < 1
    )
  cat(sprintf("      %d genes with GTEx DE classification\n", nrow(gtex_de_metrics)))
} else {
  gtex_de_metrics <- NULL
  cat("      WARNING: de_sig_all not provided — GTEX_DE_sig_in_all, GTEX_tumor_only,\n")
  cat("               GTEX_tumor_enriched will be NA for all candidates.\n")
}

# ── 4c. GTEx per-tissue median TPM (normal samples only) ─────────────────────

gtex_tissues <- unique(coldata_gtex$tissue_type[gtex_mask])

tissue_med_mat <- sapply(gtex_tissues, function(tt) {
  samp <- intersect(coldata_gtex$sample_id[coldata_gtex$tissue_type %in% tt], colnames(gtex_tpm_m))
  rowMedians(as.matrix(gtex_tpm_m[, samp, drop = FALSE]), na.rm = TRUE)
})
rownames(tissue_med_mat) <- rownames(gtex_tpm_m)

q3_mat <- sapply(gtex_tissues, function(tt) {
  samp <- intersect(coldata_gtex$sample_id[coldata_gtex$tissue_type %in% tt], colnames(gtex_tpm_m))
  rowQuantiles(as.matrix(gtex_tpm_m[, samp, drop = FALSE]), probs = 0.75, na.rm = TRUE)
})
rownames(q3_mat) <- rownames(gtex_tpm_m)

# Per-gene: tissue=Q3TPM pairs ("|"-separated) where Q3 > gtex_q3_threshold
gtex_tissue_q3_gt1 <- data.frame(
  gene_id_clean = rownames(q3_mat),
  GTEX_tissues_q3_gt1 = apply(q3_mat, 1, function(x) {
    hits <- x[!is.na(x) & x > gtex_q3_threshold]
    if (length(hits) == 0L) NA_character_
    else paste(paste0(names(hits), "=", round(hits, 1)), collapse = "|")
  }),
  stringsAsFactors = FALSE
)

gtex_tissue_stats <- data.frame(
  gene_id_clean       = rownames(gtex_tpm_m),
  GTEX_max_median_TPM = apply(tissue_med_mat, 1, max, na.rm = TRUE),
  GTEX_median_TPM     = rowMedians(gtex_tpm_m, na.rm = TRUE)
)

# Assemble gtex_data; the DE columns are handled differently depending on whether
# de_sig_all was available for this study:
#   - Present: join, then apply replace_na(FALSE) for per-gene join misses only
#   - Absent:  set all three DE columns to NA for every row (not replace_na'd to FALSE)
gtex_data <- gtex_tissue_stats %>%
  left_join(gtex_tissue_q3_gt1, by = "gene_id_clean") %>%
  mutate(across(c(GTEX_max_median_TPM, GTEX_median_TPM), ~ replace_na(.x, 0)))

if (!is.null(gtex_de_metrics)) {
  gtex_data <- gtex_data %>%
    left_join(gtex_de_metrics, by = "gene_id_clean") %>%
    mutate(across(c(GTEX_DE_sig_in_all, GTEX_tumor_only, GTEX_tumor_enriched),
                  ~ replace_na(.x, FALSE)))
} else {
  gtex_data <- gtex_data %>%
    mutate(GTEX_DE_sig_in_all  = NA,
           GTEX_tumor_only     = NA,
           GTEX_tumor_enriched = NA)
}

# Subset GTEx matrix to candidate genes for per-sample plots
gtex_tpm_sub <- gtex_tpm_m[, gtex_ids, drop = FALSE]
gtex_sample_meta    <- data.frame(
  sample_id   = gtex_ids,
  tissue_type = coldata_gtex$tissue_type[match(gtex_ids, coldata_gtex$sample_id)],
  stringsAsFactors = FALSE
)
cat(sprintf("      GTEx sub-matrix: %d genes × %d samples\n",
            nrow(gtex_tpm_sub), ncol(gtex_tpm_sub)))

rm(tumor_tpm_m, gtex_tpm_m, tissue_med_mat, gtex_de_metrics, gtex_tissue_q3_gt1,
   gtex_tissue_stats); gc()
cat("      GTEx computation complete.\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5 — TCGA QUANTIFICATION
# ─────────────────────────────────────────────────────────────────────────────

cat("[5/6] Computing TCGA expression metrics...\n")

tcga_tpm <- load_expression_data(cfg$paths$tcga_quant)
rownames(tcga_tpm) <- strip_ensg_version(rownames(tcga_tpm))
ids <- colnames(tcga_tpm)

# TCGA barcode: -01x = primary tumour, -11x/-10x = normal/peritumoral
# These are fixed TCGA standard barcode conventions, not study-specific.
is_tumor  <- grepl("-0[1-9][A-Z]$", ids)
is_normal <- grepl("-1[0-1][A-Z]$", ids)
cat(sprintf("      Tumour: %d  Normal: %d\n", sum(is_tumor), sum(is_normal)))

tcga_tumor_metrics <- compute_expression_metrics(tcga_tpm[, is_tumor,  drop = FALSE],
                                                  threshold = expr_threshold) %>%
  dplyr::rename(TCGA_tumor_num_samples = num_samples, TCGA_tumor_pct_samples = pct_samples,
                TCGA_tumor_median_TPM  = median_value, TCGA_tumor_max_TPM    = max_value) %>%
  mutate(gene_id_clean = rownames(.))

tcga_normal_metrics <- compute_expression_metrics(tcga_tpm[, is_normal, drop = FALSE],
                                                   threshold = expr_threshold) %>%
  dplyr::rename(TCGA_normal_num_samples = num_samples, TCGA_normal_pct_samples = pct_samples,
                TCGA_normal_median_TPM  = median_value, TCGA_normal_max_TPM    = max_value) %>%
  mutate(gene_id_clean = rownames(.))

tcga_coldata        <- readRDS(cfg$paths$tcga_coldata)
keep_tcga    <- is_tumor | is_normal
tcga_tpm_sub <- tcga_tpm[, keep_tcga, drop = FALSE]

ids_kept    <- colnames(tcga_tpm)[keep_tcga]
cd_idx      <- match(ids_kept, tcga_coldata$sample_id)
cancer_type <- ifelse(!is.na(cd_idx), tcga_coldata$tissue_type[cd_idx], "Unknown")
sample_type <- ifelse(is_tumor[keep_tcga], "Tumor", "Normal")
tcga_sample_meta <- data.frame(
  sample_id   = ids_kept,
  tissue_type = cancer_type,
  sample_type = sample_type,
  group       = paste(cancer_type, sample_type),
  stringsAsFactors = FALSE
)
cat(sprintf("      TCGA sub-matrix: %d genes × %d samples (%d tumour, %d normal, %d cancer types)\n",
            nrow(tcga_tpm_sub), ncol(tcga_tpm_sub), sum(is_tumor), sum(is_normal),
            n_distinct(cancer_type)))

rm(tcga_tpm); gc()
cat("      TCGA computation complete.\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6 — ASSEMBLE & SAVE
# ─────────────────────────────────────────────────────────────────────────────

cat("[6/6] Assembling final table and saving...\n")

titan_table <- ncorfs %>%
  select(
    orf_id, summary_id, gene_id, gene_name, gene_biotype,
    orf_biotype_single, protein_seq, protein_length,
    start_codon, stop_codon, chr, orf_start, orf_end, strand,
    tx_id, orf_biotypes_all, caller_count, gene_id_clean
  ) %>%
  left_join(transl_metrics,       by = "orf_id") %>%
  left_join(primary_metrics,      by = "orf_id") %>%
  left_join(cell_line_metrics,    by = "orf_id") %>%
  left_join(expr_metrics,         by = "gene_id_clean") %>%
  left_join(gtex_data,            by = "gene_id_clean") %>%
  left_join(tcga_tumor_metrics,   by = "gene_id_clean") %>%
  left_join(tcga_normal_metrics,  by = "gene_id_clean")

# ─── ORF ID normalisation ────────────────────────────────────────────────────
# Target format:
#   orf_id     = {versioned_gene_id}_{8-char md5 hex}
#   summary_id = {orf_id}_{chr}_{orf_start}_{orf_end}_{strand}_{orf_biotype_single}
# Hash inputs: protein_seq | starts | ends | chr | strand
#   where starts/ends are pipeline multi-exon coordinate columns in ncorfs.
# Rows already matching the target regex are left unchanged (source_orf_id = NA).
# Collision resolution: all rows sharing the same candidate id get _1, _2, ...
cat("\nNormalising orf_id and summary_id...\n")

CONFORMING_RE <- "^ENSG[0-9]+\\.[0-9]+_[0-9a-f]{8}$"
conforming    <- grepl(CONFORMING_RE, ncorfs$orf_id)

hash_na <- is.na(ncorfs$protein_seq) | is.na(ncorfs$starts) |
           is.na(ncorfs$ends)        | is.na(ncorfs$chr)     |
           is.na(ncorfs$strand)
to_hash <- !conforming & !hash_na

if (any(hash_na & !conforming)) {
  warning(sprintf(
    "%d non-conforming orf(s) have NA in hash input column(s) — orf_id unchanged:\n  %s",
    sum(hash_na & !conforming),
    paste(head(ncorfs$orf_id[hash_na & !conforming], 5), collapse = ", ")
  ), call. = FALSE)
}

new_orf_id <- ncorfs$orf_id
if (any(to_hash)) {
  hash_str <- paste0(
    ncorfs$protein_seq[to_hash], "|",
    ncorfs$starts[to_hash],      "|",
    ncorfs$ends[to_hash],        "|",
    ncorfs$chr[to_hash],         "|",
    ncorfs$strand[to_hash]
  )
  h8 <- substr(
    vapply(hash_str, digest::digest, character(1L), algo = "md5", USE.NAMES = FALSE),
    1L, 8L
  )
  new_orf_id[to_hash] <- paste0(ncorfs$gene_id[to_hash], "_", h8)
}

# Collision resolution — every row in a collision group gets _1, _2, ...
dupes <- unique(new_orf_id[duplicated(new_orf_id)])
if (length(dupes) > 0L) {
  cat(sprintf("  Resolving %d collision group(s) with _N suffix.\n", length(dupes)))
  oid_count <- list()
  for (i in seq_along(new_orf_id)) {
    oid <- new_orf_id[i]
    if (oid %in% dupes) {
      oid_count[[oid]] <- (oid_count[[oid]] %||% 0L) + 1L
      new_orf_id[i]    <- paste0(oid, "_", oid_count[[oid]])
    }
  }
}

source_orf_id <- ifelse(new_orf_id != ncorfs$orf_id, ncorfs$orf_id, NA_character_)
id_map        <- setNames(new_orf_id, ncorfs$orf_id)

# summary_id recomputed unconditionally (after collision resolution)
new_summary_id <- paste0(
  new_orf_id,               "_",
  ncorfs$chr,               "_",
  ncorfs$orf_start,         "_",
  ncorfs$orf_end,           "_",
  ncorfs$strand,            "_",
  ncorfs$orf_biotype_single
)
sm_map <- setNames(new_summary_id, ncorfs$orf_id)

# Apply to titan_table using the original orf_id as the lookup key
orig_ids                  <- titan_table$orf_id
titan_table$source_orf_id <- source_orf_id[match(orig_ids, ncorfs$orf_id)]
titan_table$orf_id        <- id_map[orig_ids]
titan_table$summary_id    <- sm_map[orig_ids]
titan_table               <- dplyr::relocate(titan_table, source_orf_id, .after = summary_id)

# Keep matrix rownames in sync so the app's per-ORF lookups stay valid
rownames(ribo_ppm_mat)  <- id_map[rownames(ribo_ppm_mat)]
rownames(ribocrypt_mat) <- id_map[rownames(ribocrypt_mat)]

cat(sprintf("  %d orf_ids normalised | %d already conforming | %d skipped (NA inputs)\n",
            sum(!is.na(titan_table$source_orf_id)), sum(conforming), sum(hash_na & !conforming)))

cat(sprintf("      Final table: %d ORFs × %d columns\n",
            nrow(titan_table), ncol(titan_table)))

# ─── Coverage summary ────────────────────────────────────────────────────────
cat("\nCoverage summary:\n")
pct <- function(x) sprintf("%.1f%%", 100 * mean(!is.na(x)))
cat(sprintf("  target_expression      : %s have data\n", pct(titan_table$target_expression_median_TPM)))
if (is.null(cfg$paths$de_sig_all)) {
  cat("  GTEX_DE_sig_in_all     : de_sig_all not provided for this study — all NA\n")
  cat("  GTEX_tumor_only        : de_sig_all not provided for this study — all NA\n")
  cat("  GTEX_tumor_enriched    : de_sig_all not provided for this study — all NA\n")
} else {
  cat(sprintf("  GTEX_DE_sig_in_all     : %s have data\n", pct(titan_table$GTEX_DE_sig_in_all)))
  cat(sprintf("  GTEX_tumor_only        : %s have data\n", pct(titan_table$GTEX_tumor_only)))
}
cat(sprintf("  GTEX_max_median_TPM    : %s have data\n", pct(titan_table$GTEX_max_median_TPM)))
cat(sprintf("  TCGA_tumor_median_TPM  : %s have data\n", pct(titan_table$TCGA_tumor_median_TPM)))
if (!has_riboseq) {
  cat("  target_translation_PPM : ribo_ppm/ribo_psites not provided for this study — all NA\n")
} else {
  cat(sprintf("  target_translation_PPM : %s have data\n", pct(titan_table$target_translation_median_PPM)))
}
cat(sprintf("  ribocrypt_primary_PPM  : %s have data\n", pct(titan_table$ribocrypt_primary_median_PPM)))

# ─── Save ────────────────────────────────────────────────────────────────────
# Convert expression matrices to 32-bit float before saving.
# Halves their in-memory footprint; the app coerces back to double at the point
# of use (as.numeric / arithmetic), so no downstream code changes are needed.
if (!requireNamespace("float", quietly = TRUE))
  stop("Package 'float' is required. Install with: install.packages('float')")
to_fl <- function(m) if (!is.null(m) && is.matrix(m)) float::fl(m) else m

app_data <- list(
  orf_table             = titan_table,
  ribo_ppm_samples      = to_fl(ribo_ppm_mat),
  rna_tpm_mat           = to_fl(rna_tpm_sub),
  ribo_sample_meta      = ribo_sample_meta,
  rna_sample_meta       = rna_sample_meta,
  gtex_tpm_mat          = to_fl(gtex_tpm_sub),
  gtex_sample_meta      = gtex_sample_meta,
  tcga_tpm_mat          = to_fl(tcga_tpm_sub),
  tcga_sample_meta      = tcga_sample_meta,
  ribocrypt_mat         = to_fl(ribocrypt_mat),
  ribocrypt_sample_meta = ribocrypt_sample_meta,
  ribocrypt_meta        = list(
    primary_samples   = sample_classes$primary,
    cell_line_samples = sample_classes$cell_line
  ),
  study_id    = cfg$study_id,
  prepared_on = Sys.time()
)

saveRDS(app_data, output_rds)
write.csv(titan_table, output_csv, row.names = FALSE)

cat(sprintf("\nDone!\n  %s\n  %s\n", output_rds, output_csv))
