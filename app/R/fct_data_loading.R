## Shared validation for TITAN RDS objects.
## Called by register_study.R (catalog gating) and app.R upload handler.

validate_titan_rds <- function(obj) {
  required_elements <- c(
    "orf_table", "ribo_ppm_samples", "rna_tpm_mat",
    "ribo_sample_meta", "rna_sample_meta",
    "gtex_tpm_mat", "gtex_sample_meta",
    "tcga_tpm_mat", "tcga_sample_meta",
    "ribocrypt_mat", "ribocrypt_sample_meta", "ribocrypt_meta"
  )
  required_orf_cols <- c(
    "orf_id", "gene_id", "gene_name", "orf_biotype_single",
    "protein_seq", "protein_length", "chr", "orf_start", "orf_end", "strand"
  )

  if (!is.list(obj) || is.data.frame(obj))
    stop("Expected a named list; got: ", paste(class(obj), collapse = ", "), call. = FALSE)

  missing_els <- setdiff(required_elements, names(obj))
  if (length(missing_els) > 0)
    stop("RDS is missing required elements: ", paste(missing_els, collapse = ", "), call. = FALSE)

  if (!is.data.frame(obj$orf_table))
    stop("$orf_table must be a data.frame", call. = FALSE)

  missing_cols <- setdiff(required_orf_cols, names(obj$orf_table))
  if (length(missing_cols) > 0)
    stop("$orf_table is missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)

  is_mat <- function(x) is.matrix(x) || inherits(x, "float32")

  n_orfs <- nrow(obj$orf_table)
  if (!is_mat(obj$ribo_ppm_samples) || nrow(obj$ribo_ppm_samples) != n_orfs)
    stop(sprintf(
      "$ribo_ppm_samples must be a matrix with %d rows (orf_table has %d rows); got %s with %s rows",
      n_orfs, n_orfs,
      paste(class(obj$ribo_ppm_samples), collapse = "/"),
      if (is_mat(obj$ribo_ppm_samples)) nrow(obj$ribo_ppm_samples) else "N/A"),
      call. = FALSE)

  if (!is_mat(obj$ribocrypt_mat) || nrow(obj$ribocrypt_mat) != n_orfs)
    stop(sprintf(
      "$ribocrypt_mat must be a matrix with %d rows; got %s rows",
      n_orfs, if (is_mat(obj$ribocrypt_mat)) nrow(obj$ribocrypt_mat) else "N/A"),
      call. = FALSE)

  invisible(TRUE)
}
