#!/usr/bin/env Rscript
# 02_prep_annotation.R
# Build offline gene-annotation lookup: Ensembl 114 gene ID → gene name,
# UniProt/SwissProt accession, and gene description.
# Used by the app to annotate BLAST homology hits without live biomaRt calls.
#
# Run via sbatch (needs internet access for one-time biomaRt query):
#   sbatch app/scripts/02_prep_annotation.sbatch
#
# Output:
#   ref/ensembl114_pep/ensembl_gene_annotation.rds

suppressPackageStartupMessages(library(biomaRt))

APP_DIR  <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = FALSE)
OUT_FILE <- file.path(APP_DIR, "ref", "ensembl114_pep", "ensembl_gene_annotation.rds")

# ── Connect to Ensembl 114 archive ────────────────────────────────────────────
message("Connecting to Ensembl 114 BioMart archive…")
mart <- tryCatch(
  useEnsembl("genes", dataset = "hsapiens_gene_ensembl", version = 114),
  error = function(e) {
    message("useEnsembl() failed, trying archive mirror…")
    useMart("ensembl",
            dataset = "hsapiens_gene_ensembl",
            host    = "https://oct2024.archive.ensembl.org")
  }
)
message("  Connected: ", mart@host)

# ── Fetch attributes ──────────────────────────────────────────────────────────
message("Fetching gene annotations (this takes a few minutes)…")
attrs <- c(
  "ensembl_gene_id",      # ENSG ID
  "external_gene_name",   # HGNC gene symbol
  "uniprotswissprot",     # UniProt/SwissProt accession (may be empty)
  "description"           # Gene description from Ensembl
)

annot <- getBM(attributes = attrs, mart = mart)

# Clean: description sometimes has " [Source:…]" suffix — strip it
annot$description <- sub("\\s*\\[Source:.*\\]$", "", annot$description)

# Remove empty UniProt rows (keep one row per ENSG even without UniProt)
annot <- annot[!duplicated(paste0(annot$ensembl_gene_id, "|", annot$uniprotswissprot)), ]

message(sprintf("  %d gene–UniProt pairs fetched (%d unique ENSG IDs)",
                nrow(annot), length(unique(annot$ensembl_gene_id))))

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(annot, OUT_FILE, compress = "xz")
message("  → ", OUT_FILE)
message("\n--- DONE ---")
message("All three prep steps complete. Restart the Shiny app to load new references.")
