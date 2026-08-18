#!/usr/bin/env Rscript
# 03_prep_gencode_kmer_index.R
# Precompute the k-mer seed index that match_peptides() (app/R/fct_peptides.R)
# builds internally from an ORF table's protein_seq column.
#
# gencode_orf_tbl (app/ref/gencode_orfs_phase2.csv, ~28,300 ORFs) is static -
# it never changes between MS-peptide uploads - but match_peptides() rebuilds
# its k-mer index (hashing + radix-sorting every protein_seq) from scratch on
# every call. That rebuild is the "Running Gencode cross-match…" step users see
# on each upload/study switch. Building the index once here and loading it as
# gencode_kmer_index in global.R lets match_peptides() skip straight to the
# peptide-side lookup instead.
#
# Run via sbatch (containerised R):
#   singularity exec \
#     --env "R_LIBS_USER=/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/Rstudio_4.4.0_3.19_libs" \
#     /hpc/local/Rocky8/pmc_vanheesch/singularity_images/bioconductor_docker_RELEASE_3_19.sif \
#     Rscript scripts/reference_prep/03_prep_gencode_kmer_index.R
#
# Output:
#   app/ref/gencode_kmer_index.rds
#
# Rerun whenever app/ref/gencode_orfs_phase2.csv changes — match_peptides()
# detects a stale index (orf_id mismatch) and falls back to a live rebuild,
# but that silently loses the speedup, so keep this in sync.

.args     <- commandArgs(trailingOnly = FALSE)
.this_file <- sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])
APP_DIR   <- normalizePath(file.path(dirname(.this_file), "..", "..", "app"), mustWork = TRUE)
GENCODE_CSV <- file.path(APP_DIR, "ref", "gencode_orfs_phase2.csv")
OUT_FILE  <- file.path(APP_DIR, "ref", "gencode_kmer_index.rds")

source(file.path(APP_DIR, "R", "fct_gencode_orf.R"))
source(file.path(APP_DIR, "R", "fct_peptides.R"))

message("Loading Gencode ORF table: ", GENCODE_CSV)
gencode_orf_tbl <- load_gencode_orf_table(GENCODE_CSV)
if (is.null(gencode_orf_tbl))
  stop("Failed to load Gencode ORF table from ", GENCODE_CSV, call. = FALSE)
message(sprintf("  %s ORFs loaded", formatC(nrow(gencode_orf_tbl), big.mark = ",")))

message("Building k-mer seed index (k = 6)…")
t0 <- proc.time()
index <- build_pep_kmer_index(gencode_orf_tbl, k = 6L)
elapsed <- (proc.time() - t0)["elapsed"]
message(sprintf("  built in %.1f s — %s unique k-mers", elapsed, formatC(length(index$uniq_hash), big.mark = ",")))

message("Saving index RDS…")
saveRDS(index, OUT_FILE, compress = "xz")
message("  → ", OUT_FILE)

message("\n--- DONE ---")
message("Restart the Shiny app to load the new gencode_kmer_index.")
