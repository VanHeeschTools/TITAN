#!/usr/bin/env Rscript
#
# TITAN — Tumour Immunopeptidomics Target ANnotation
# Data preparation script: computes all 31 input_fields per ORF candidate
# and saves the integrated table for the Shiny app.
#
# Run once on HPC before launching the app:
#   Rscript prepare_titan_inputs.R
#
# Outputs:
#   app/data/titan_orf_table.rds   — full table + per-sample matrices (app input)
#   app/data/titan_orf_table.csv   — flat CSV version for inspection

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(tximport)
  library(matrixStats)
})

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

BASE_TARGET   <- "/hpc/pmc_oatv/projects/tnbc_navarra/"
BASE_TITAN    <- "/hpc/pmc_oatv/projects/tools_dev/titan/app/data/"

# Value of coldata$tissue_type that identifies the target tumour samples.
# All other tissue_type values in the coldata are treated as GTEx normal tissues.
TARGET_TUMOR_TYPE <- "TNBC"

PATHS <- list(
  # Ribo-seq inputs
  ncorfs        = file.path(BASE_TARGET, "analysis/navarra/riboseq_pipeline/harmonise_orfs/harmonised_orf_table.csv"),
  ribo_ppm      = file.path(BASE_TARGET, "analysis/navarra/riboseq_pipeline/orf_expression/orf_table_psites_permillion.csv"),
  ribo_psites   = file.path(BASE_TARGET, "analysis/navarra/riboseq_pipeline/orf_expression/orf_table_psites.csv"),
  ribocrypt_ext = file.path(BASE_TARGET, "analysis/navarra/riboseq_database_quantification/ribocrypt/ribocrypt_tnbc_quantification_psites_permillion.csv"),

  # Combined GTEx + tumour txi (coldata$tissue_type distinguishes tumour from GTEx tissues)
  de_sig_all    = file.path(BASE_TARGET, "results/DE/TNBC_polyA_DE_sig_all.tsv"),
  gtex_txi      = file.path(BASE_TARGET, "analysis/GTEx_DE/rds/TNBC_polyA_GTEx_txi_filtered_raw.RDS"),
  gtex_coldata  = file.path(BASE_TARGET, "analysis/GTEx_DE/rds/TNBC_polyA_GTEx_coldata_filtered.RDS"),

  # TCGA inputs
  tcga_txi      = file.path("/hpc/pmc_vanheesch/shared_resources/quantification/TCGA_matched_TN_quantification/gencode_48_quantification/rds_objects/TCGA_matched_TN_gencode48_genes.RDS"),

  # Outputs
  output_dir    = file.path(BASE_TITAN, ""),
  output_rds    = file.path(BASE_TITAN, "titan_orf_table.rds"),
  output_csv    = file.path(BASE_TITAN, "titan_orf_table.csv")
)

dir.create(PATHS$output_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

strip_ensg_version <- function(x) sub("\\..*", "", x)

# Computes num_samples, pct_samples, median, max from a gene/ORF × sample matrix
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

# Classifies external ribocrypt sample names as primary tissue or cell line.
# Primary: samples prefixed with "primary_", "hepatocyte_", "Myoblast_",
#          "HSPC_", or "huvec_" — all others are treated as cell lines.
classify_ribocrypt_samples <- function(sample_names) {
  primary_patterns <- c("^primary_", "^hepatocyte_", "^Myoblast_", "^HSPC_", "^huvec_")
  is_primary <- sapply(sample_names, function(s) {
    any(sapply(primary_patterns, function(p) grepl(p, s)))
  })
  list(primary = sample_names[is_primary], cell_line = sample_names[!is_primary])
}

# Infer RMS ribo-seq sample condition from sample name
get_riboseq_condition <- function(sample_names) {
  ifelse(grepl("IFNg|IFN", sample_names, ignore.case = TRUE), "IFNg", "Unstimulated")
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 — ORF BACKBONE
# ─────────────────────────────────────────────────────────────────────────────

cat("=== TITAN input preparation ===\n\n")
cat("[1/6] Loading ncORF candidate table...\n")

ncorfs <- fread(PATHS$ncorfs, data.table = FALSE)
ncorfs$gene_id_clean <- strip_ensg_version(ncorfs$gene_id)

cat(sprintf("      %d ncORF candidates, %d unique genes\n",
            nrow(ncorfs), n_distinct(ncorfs$gene_id_clean)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2 — TARGET TRANSLATION (RMS ribo-seq)
# ─────────────────────────────────────────────────────────────────────────────

cat("[2/6] Computing target translation metrics (RMS ribo-seq)...\n")

ribo_ppm    <- fread(PATHS$ribo_ppm,    data.table = FALSE)
ribo_psites <- fread(PATHS$ribo_psites, data.table = FALSE)

rownames(ribo_ppm)    <- ribo_ppm$orf_id;    ribo_ppm$orf_id    <- NULL
rownames(ribo_psites) <- ribo_psites$orf_id; ribo_psites$orf_id <- NULL

# Shorten sample IDs for display (keep leading numeric prefix only)
ribo_ppm    <- rename_with(ribo_ppm,    ~ sub("-SL_.*|-EW_.*|.*TIS-|.*ORG-", "", .x))
ribo_psites <- rename_with(ribo_psites, ~ sub("-SL_.*|-EW_.*|.*TIS-|.*ORG-", "", .x))

candidate_ids <- ncorfs$orf_id
common_ribo   <- intersect(candidate_ids, rownames(ribo_ppm))

ribo_ppm_mat    <- as.matrix(ribo_ppm[common_ribo, ])
ribo_psites_mat <- as.matrix(ribo_psites[common_ribo, ])

transl_metrics <- compute_expression_metrics(ribo_ppm_mat, threshold = 1) %>%
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
  condition = get_riboseq_condition(colnames(ribo_ppm_mat))
)

cat(sprintf("      %d ORFs with ribo-seq data, %d samples\n",
            nrow(transl_metrics),
            nrow(ribo_sample_meta)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 — RIBOCRYPT EXTERNAL (primary tissue / cell line)
# ─────────────────────────────────────────────────────────────────────────────

cat("[3/6] Computing ribocrypt external metrics...\n")

ribocrypt_ext <- fread(PATHS$ribocrypt_ext, data.table = FALSE)
rownames(ribocrypt_ext) <- ribocrypt_ext$orf_id; ribocrypt_ext$orf_id <- NULL

sample_classes <- classify_ribocrypt_samples(colnames(ribocrypt_ext))
cat(sprintf("      %d primary tissue samples: %s...\n",
            length(sample_classes$primary),
            paste(head(sample_classes$primary, 3), collapse = ", ")))
cat(sprintf("      %d cell-line samples: %s...\n",
            length(sample_classes$cell_line),
            paste(head(sample_classes$cell_line, 3), collapse = ", ")))

common_rc  <- intersect(candidate_ids, rownames(ribocrypt_ext))
rc_mat     <- as.matrix(ribocrypt_ext[common_rc, ])

primary_mat   <- rc_mat[, sample_classes$primary,   drop = FALSE]
cell_line_mat <- rc_mat[, sample_classes$cell_line, drop = FALSE]

primary_metrics <- compute_expression_metrics(primary_mat, threshold = 1) %>%
  rename(
    ribocrypt_primary_num_samples = num_samples,
    ribocrypt_primary_pct_samples = pct_samples,
    ribocrypt_primary_median_PPM  = median_value,
    ribocrypt_primary_max_PPM     = max_value
  ) %>%
  mutate(orf_id = rownames(.))

cell_line_metrics <- compute_expression_metrics(cell_line_mat, threshold = 1) %>%
  rename(
    `ribocrypt_cell-line_num_samples` = num_samples,
    `ribocrypt_cell-line_pct_samples` = pct_samples,
    `ribocrypt_cell-line_median_PPM`  = median_value,
    `ribocrypt_cell-line_max_PPM`     = max_value
  ) %>%
  mutate(orf_id = rownames(.))

rm(ribocrypt_ext, rc_mat, primary_mat, cell_line_mat); gc()
cat("      Done.\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4 — TARGET EXPRESSION + GTEX METRICS (single combined txi)
#
# The txi object contains both tumour samples (tissue_type == TARGET_TUMOR_TYPE)
# and GTEx normal tissue samples (all other tissue_type values).
# ─────────────────────────────────────────────────────────────────────────────

cat("[4/6] Loading combined txi (~1.2 GB, may take a few minutes)...\n")

txi_gtex     <- readRDS(PATHS$gtex_txi)
coldata_gtex <- readRDS(PATHS$gtex_coldata)

# Validate that the required column is present
stopifnot("tissue_type" %in% colnames(coldata_gtex))

tumor_mask <- coldata_gtex$tissue_type %in% TARGET_TUMOR_TYPE
gtex_mask  <- coldata_gtex$tissue_type != TARGET_TUMOR_TYPE & !is.na(coldata_gtex$tissue_type)

cat(sprintf("      Tumour samples (%s): %d\n", TARGET_TUMOR_TYPE, sum(tumor_mask)))
cat(sprintf("      GTEx normal samples: %d across %d tissues\n",
            sum(gtex_mask),
            n_distinct(coldata_gtex$tissue_type[gtex_mask])))

# Strip ENSG version from txi row names once (applies to both subsets)
txi_gene_ids <- strip_ensg_version(rownames(txi_gtex$abundance))

# ── 4a. Target expression (tumour samples) ───────────────────────────────────

tumor_ids   <- coldata_gtex$sample_id[tumor_mask]
tumor_tpm_m <- txi_gtex$abundance[, tumor_ids, drop = FALSE]
rownames(tumor_tpm_m) <- txi_gene_ids

gene_ids_available <- intersect(ncorfs$gene_id_clean, rownames(tumor_tpm_m))
rna_tpm_sub <- tumor_tpm_m[gene_ids_available, , drop = FALSE]

expr_metrics <- compute_expression_metrics(rna_tpm_sub, threshold = 1) %>%
  dplyr::rename(
    target_expression_num_samples = num_samples,
    target_expression_pct_samples = pct_samples,
    target_expression_median_TPM  = median_value,
    target_expression_max_TPM     = max_value
  ) %>%
  mutate(gene_id_clean = rownames(.))

rna_sample_meta <- data.frame(
  sample_id   = tumor_ids,
  tissue_type = TARGET_TUMOR_TYPE
)

cat(sprintf("      %d genes with target expression data (%d samples)\n",
            nrow(expr_metrics), length(tumor_ids)))

# ── 4b. GTEx DE classification (from pre-computed DE table) ──────────────────

de_sig <- read.delim(PATHS$de_sig_all, check.names = FALSE)
de_sig$gene_id_clean <- strip_ensg_version(de_sig$gene_id)

gtex_de_metrics <- de_sig %>%
  transmute(
    gene_id_clean,
    GTEX_DE_sig_in_all  = sig_in_all,
    GTEX_tumor_only     = low_all_tissues,
    GTEX_tumor_enriched = Q3_GTEx < 1
  )

cat(sprintf("      %d genes with GTEx DE classification\n", nrow(gtex_de_metrics)))

# ── 4c. GTEx per-tissue median TPM (normal samples only) ─────────────────────

gtex_ids   <- coldata_gtex$sample_id[gtex_mask]
gtex_tpm_m <- txi_gtex$abundance[, gtex_ids, drop = FALSE]
rownames(gtex_tpm_m) <- txi_gene_ids

gtex_tissues <- unique(coldata_gtex$tissue_type[gtex_mask])

tissue_med_mat <- sapply(gtex_tissues, function(tt) {
  samp <- coldata_gtex$sample_id[coldata_gtex$tissue_type %in% tt]
  rowMedians(as.matrix(gtex_tpm_m[, samp, drop = FALSE]), na.rm = TRUE)
})
rownames(tissue_med_mat) <- rownames(gtex_tpm_m)


q3_mat <- sapply(gtex_tissues, function(tt) {
  samp <- coldata_gtex$sample_id[coldata_gtex$tissue_type %in% tt]
  rowQuantiles(as.matrix(gtex_tpm_m[, samp, drop = FALSE]), probs = 0.75, na.rm = TRUE)
})
rownames(q3_mat) <- rownames(gtex_tpm_m)

# Per-gene: tissue=Q3TPM pairs ("|"-separated) where Q3 > 1 TPM
gtex_tissue_q3_gt1 <- data.frame(
  gene_id_clean = rownames(q3_mat),
  GTEX_tissues_q3_gt1 = apply(q3_mat, 1, function(x) {
    hits <- x[!is.na(x) & x > 1]
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

gtex_data <- gtex_tissue_stats %>%
  left_join(gtex_de_metrics,     by = "gene_id_clean") %>%
  left_join(gtex_tissue_q3_gt1,  by = "gene_id_clean") %>%
  # Only replace NA for logical / numeric columns — leave GTEX_tissues_q3_gt1 as NA
  mutate(across(c(GTEX_max_median_TPM, GTEX_median_TPM), ~ replace_na(.x, 0)),
         across(c(GTEX_DE_sig_in_all, GTEX_tumor_only, GTEX_tumor_enriched),
                ~ replace_na(.x, FALSE)))

rm(txi_gtex, tumor_tpm_m, gtex_tpm_m, tissue_med_mat, gtex_de_metrics, gtex_tissue_q3_gt1,
   gtex_tissue_stats); gc()
cat("      GTEx computation complete.\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6 — TCGA QUANTIFICATION
# ─────────────────────────────────────────────────────────────────────────────

cat("[5/6] Computing TCGA expression metrics (514 samples via tximport)...\n")


txi_tcga <- readRDS(PATHS$tcga_txi)
tcga_tpm <- txi_tcga$abundance   # gene × sample (no-version ENSG IDs)
rownames(tcga_tpm) <- strip_ensg_version(rownames(tcga_tpm))
ids       <- colnames(tcga_tpm)

# TCGA barcode: -01x = primary tumour, -11x/-10x = normal/peritumoral
is_tumor  <- grepl("-0[1-9][A-Z]$", ids)
is_normal <- grepl("-1[0-1][A-Z]$", ids)
cat(sprintf("      Tumour: %d  Normal: %d\n", sum(is_tumor), sum(is_normal)))

tcga_tumor_metrics <- compute_expression_metrics(tcga_tpm[, is_tumor,  drop = FALSE]) %>%
  dplyr::rename(TCGA_tumor_num_samples = num_samples, TCGA_tumor_pct_samples = pct_samples,
         TCGA_tumor_median_TPM  = median_value, TCGA_tumor_max_TPM    = max_value) %>%
  mutate(gene_id_clean = rownames(.))

tcga_normal_metrics <- compute_expression_metrics(tcga_tpm[, is_normal, drop = FALSE]) %>%
  dplyr::rename(TCGA_normal_num_samples = num_samples, TCGA_normal_pct_samples = pct_samples,
         TCGA_normal_median_TPM  = median_value, TCGA_normal_max_TPM    = max_value) %>%
  mutate(gene_id_clean = rownames(.))

rm(txi_tcga, tcga_tpm); gc()
cat("      TCGA computation complete.\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7 — ASSEMBLE & SAVE
# ─────────────────────────────────────────────────────────────────────────────

cat("[6/6] Assembling final table and saving...\n")

titan_table <- ncorfs %>%
  select(
    orf_id, summary_id, gene_id, gene_name, gene_biotype,
    orf_biotype_single, protein_seq, protein_length,
    start_codon, stop_codon, chr, orf_start, orf_end, strand,
    tx_id, orf_biotypes_all, caller_count, gene_id_clean
  ) %>%
  # ORF-level joins (ribo-seq)
  left_join(transl_metrics,    by = "orf_id") %>%
  left_join(primary_metrics,   by = "orf_id") %>%
  left_join(cell_line_metrics, by = "orf_id") %>%
  # Gene-level joins (RNA-seq)
  left_join(expr_metrics,          by = "gene_id_clean") %>%
  left_join(gtex_data,     by = "gene_id_clean") %>%
  left_join(tcga_tumor_metrics,    by = "gene_id_clean") %>%
  left_join(tcga_normal_metrics,   by = "gene_id_clean")

cat(sprintf("      Final table: %d ORFs × %d columns\n",
            nrow(titan_table), ncol(titan_table)))

# ─── Coverage summary ────────────────────────────────────────────────────────
cat("\nCoverage summary:\n")
pct <- function(x) sprintf("%.1f%%", 100 * mean(!is.na(x)))
cat(sprintf("  target_expression      : %s have data\n", pct(titan_table$target_expression_median_TPM)))
cat(sprintf("  GTEX_DE_sig_in_all     : %s have data\n", pct(titan_table$GTEX_DE_sig_in_all)))
cat(sprintf("  GTEX_max_median_TPM    : %s have data\n", pct(titan_table$GTEX_max_median_TPM)))
cat(sprintf("  TCGA_tumor_median_TPM  : %s have data\n", pct(titan_table$TCGA_tumor_median_TPM)))
cat(sprintf("  target_translation_PPM : %s have data\n", pct(titan_table$target_translation_median_PPM)))
cat(sprintf("  ribocrypt_primary_PPM  : %s have data\n", pct(titan_table$ribocrypt_primary_median_PPM)))

# ─── Save ────────────────────────────────────────────────────────────────────
app_data <- list(
  orf_table             = titan_table,
  ribo_ppm_samples      = ribo_ppm_mat,      # ORF × RMS-sample matrix (for per-sample plots)
  rna_tpm_mat           = rna_tpm_sub,       # gene × RMS-sample matrix (for dynamic TPM threshold)
  ribo_sample_meta      = ribo_sample_meta,
  rna_sample_meta       = rna_sample_meta,
  ribocrypt_meta        = list(
    primary_samples   = sample_classes$primary,
    cell_line_samples = sample_classes$cell_line
  ),
  prepared_on = Sys.time()
)

saveRDS(app_data, PATHS$output_rds)
write.csv(titan_table, PATHS$output_csv, row.names = FALSE)

cat(sprintf("\nDone!\n  %s\n  %s\n", PATHS$output_rds, PATHS$output_csv))

