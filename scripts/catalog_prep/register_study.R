#!/usr/bin/env Rscript
#
# TITAN — register a prepared study into the app study catalog.
#
# Usage:
#   Rscript register_study.R <config.yaml> <rds_path> [--catalog <catalog.yaml>]
#
# Arguments:
#   config.yaml   The prepare_titan_inputs.R config used for this study.
#                 Supplies: study_id, display_name, cancer_type, cohort.
#   rds_path      Path to the titan_<study_id>.rds output file.
#   --catalog     Path to catalog.yaml. Defaults to ../app/data/catalog.yaml
#                 relative to this script.
#
# If study_id already exists in the catalog, its entry is updated in-place.
# Prints a diff before writing so you can review the change.

suppressPackageStartupMessages({
  library(yaml)
})

# sys.frame(1)$ofile only resolves when this script is source()'d — it errors
# ("not that many frames on the stack") when run directly via `Rscript
# register_study.R ...` (this script's own documented usage). Fall back to
# parsing --file= from commandArgs() for the direct-invocation case.
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0)
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  frame_ofile <- sys.frame(1)$ofile
  if (!is.null(frame_ofile)) return(dirname(frame_ofile))
  getwd()
}

script_dir <- get_script_dir()

source(file.path(script_dir, "../../app/R/fct_data_loading.R"))

# ─── CLI ─────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

catalog_flag <- which(args == "--catalog")
if (length(catalog_flag) > 0) {
  catalog_path <- args[catalog_flag + 1L]
  args <- args[-c(catalog_flag, catalog_flag + 1L)]
} else {
  catalog_path <- normalizePath(
    file.path(script_dir, "../../app/data/catalog.yaml"),
    mustWork = FALSE
  )
}

if (length(args) != 2L)
  stop("Usage: Rscript register_study.R <config.yaml> <rds_path> [--catalog <path>]",
       call. = FALSE)

config_path <- args[1L]
rds_path    <- args[2L]

if (!file.exists(config_path))
  stop("Config not found: ", config_path, call. = FALSE)
if (!file.exists(rds_path))
  stop("RDS not found: ", rds_path, call. = FALSE)

# ─── Read metadata from config ───────────────────────────────────────────────

cfg <- yaml::read_yaml(config_path)

required_meta <- c("study_id", "display_name", "cancer_type", "cohort")
missing_meta  <- setdiff(required_meta, names(cfg))
if (length(missing_meta) > 0)
  stop("Config is missing required fields: ", paste(missing_meta, collapse = ", "), call. = FALSE)

# ─── Validate and inspect RDS ─────────────────────────────────────────────────

cat(sprintf("Loading %s …\n", rds_path))
obj <- tryCatch(readRDS(rds_path), error = function(e)
  stop("Failed to read RDS: ", e$message, call. = FALSE))

cat("Validating structure …\n")
tryCatch(
  validate_titan_rds(obj),
  error = function(e) stop("Validation failed: ", e$message, call. = FALSE)
)

n_orfs         <- nrow(obj$orf_table)
n_ribo_samples <- ncol(obj$ribo_ppm_samples)
n_rna_samples  <- if (!is.null(obj$rna_tpm_mat)) ncol(obj$rna_tpm_mat) else
                    nrow(obj$rna_sample_meta)
prepared_on    <- if (!is.null(obj$prepared_on))
                    format(as.Date(obj$prepared_on), "%Y-%m-%d") else
                    format(Sys.Date(), "%Y-%m-%d")

# ─── Compute rds_path relative to app/ ───────────────────────────────────────

app_dir    <- normalizePath(file.path(dirname(catalog_path), ".."), mustWork = FALSE)
# The catalog lives in app/data/ so its parent is app/
app_dir    <- normalizePath(dirname(catalog_path), mustWork = FALSE)
# Use path relative to app_dir (the directory containing catalog.yaml's parent data/)
app_dir    <- normalizePath(file.path(dirname(catalog_path), ".."), mustWork = FALSE)
rds_abs    <- normalizePath(rds_path, mustWork = TRUE)
rds_rel    <- sub(paste0("^", app_dir, "/"), "", rds_abs)

# Sanity check: relative path must not be the same as absolute
if (rds_rel == rds_abs)
  warning("Could not make rds_path relative to app dir. Storing absolute path.")

# ─── Build new entry ─────────────────────────────────────────────────────────

new_entry <- list(
  study_id       = cfg$study_id,
  display_name   = cfg$display_name,
  cancer_type    = cfg$cancer_type,
  cohort         = cfg$cohort,
  rds_path       = rds_rel,
  n_orfs         = n_orfs,
  n_ribo_samples = n_ribo_samples,
  n_rna_samples  = n_rna_samples,
  prepared_on    = prepared_on
)

# ─── Read existing catalog ───────────────────────────────────────────────────

if (file.exists(catalog_path)) {
  cat_yaml <- yaml::read_yaml(catalog_path)
  studies  <- cat_yaml$studies %||% list()
} else {
  cat_yaml <- list()
  studies  <- list()
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Check for existing entry
existing_idx <- which(sapply(studies, `[[`, "study_id") == cfg$study_id)

if (length(existing_idx) > 0) {
  old <- studies[[existing_idx[1L]]]
  cat(sprintf("\nUpdating existing entry for '%s':\n", cfg$study_id))
  for (k in names(new_entry)) {
    ov <- old[[k]]; nv <- new_entry[[k]]
    if (!identical(ov, nv))
      cat(sprintf("  %-18s  %s  →  %s\n", k,
                  if (is.null(ov)) "(null)" else as.character(ov),
                  as.character(nv)))
  }
  studies[[existing_idx[1L]]] <- new_entry
} else {
  cat(sprintf("\nAdding new entry for '%s':\n", cfg$study_id))
  for (k in names(new_entry))
    cat(sprintf("  %-18s  %s\n", k, as.character(new_entry[[k]])))
  studies <- c(studies, list(new_entry))
}

# ─── Write catalog ───────────────────────────────────────────────────────────

cat_yaml$studies <- studies
cat(sprintf("\nWriting catalog: %s\n", catalog_path))
yaml::write_yaml(cat_yaml, catalog_path)
cat("Done.\n")
