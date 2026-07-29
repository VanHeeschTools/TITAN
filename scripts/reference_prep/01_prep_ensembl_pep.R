#!/usr/bin/env Rscript
# 01_prep_ensembl_pep.R
# Download Ensembl 114 pep FASTA, deduplicate by protein sequence,
# build in-app index (RDS) and blastp database.
#
# Run via sbatch (needs internet access + makeblastdb in PATH):
#   sbatch scripts/01_prep_ensembl_pep.sbatch
#
# Outputs (relative to app/):
#   ref/ensembl114_pep/Homo_sapiens.GRCh38.pep.all.fa   (full FASTA)
#   ref/ensembl114_pep/ensembl114_pep_dedup.fa           (deduplicated FASTA, for blastp db)
#   ref/ensembl114_pep/ensembl114_pep_index.rds          (in-app lookup index)
#   ref/ensembl114_pep/ensembl114_pep.*                  (blastp database files)

suppressPackageStartupMessages(library(digest))

APP_DIR  <- getwd()
OUT_DIR  <- file.path(APP_DIR, "ref", "ensembl114_pep")
MAKEBLASTDB <- file.path(APP_DIR, "bin", "ncbi-blast-2.14.0+", "bin", "makeblastdb")
FASTA_GZ <- file.path(OUT_DIR, "Homo_sapiens.GRCh38.pep.all.fa.gz")
FASTA    <- file.path(OUT_DIR, "Homo_sapiens.GRCh38.pep.all.fa")
DEDUP_FA <- file.path(OUT_DIR, "ensembl114_pep_dedup.fa")
INDEX_RDS <- file.path(OUT_DIR, "ensembl114_pep_index.rds")
BLAST_OUT <- file.path(OUT_DIR, "ensembl114_pep")

# ── 1. Download ───────────────────────────────────────────────────────────────
URL <- "https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

if (!file.exists(FASTA)) {
  if (!file.exists(FASTA_GZ)) {
    message("Downloading Ensembl 114 pep FASTA…")
    download.file(URL, FASTA_GZ, method = "wget", quiet = FALSE)
  }
  message("Decompressing…")
  system2("gunzip", args = c("-k", FASTA_GZ))   # -k keeps the .gz
}
message("FASTA at: ", FASTA)

# ── 2. Parse headers + sequences ─────────────────────────────────────────────
# Ensembl 114 header format:
# >ENSP00000493376.2 pep chromosome:GRCh38:... gene:ENSG00000277248.1
#   transcript:ENST00000631435.1 gene_biotype:TR_V_gene
#   transcript_biotype:TR_V_gene gene_symbol:TRAV6 description:T cell ...
message("Parsing FASTA (this takes ~1 min for 118K sequences)…")
lines   <- readLines(FASTA, warn = FALSE)
hdr_idx <- which(startsWith(lines, ">"))
n       <- length(hdr_idx)
message("  Total sequences: ", n)

ensp_ids  <- character(n)
ensg_ids  <- character(n)
enst_ids  <- character(n)
gene_syms <- character(n)
seqs      <- character(n)

for (i in seq_len(n)) {
  end        <- if (i < n) hdr_idx[i + 1L] - 1L else length(lines)
  seqs[i]    <- paste(lines[(hdr_idx[i] + 1L):end], collapse = "")
  hdr        <- lines[hdr_idx[i]]

  # Protein ID (first token after >)
  ensp_ids[i] <- sub("^>([^ .]+).*", "\\1", hdr)   # strip .version suffix

  # gene:ENSG…
  m <- regmatches(hdr, regexpr("gene:ENSG[0-9]+", hdr))
  ensg_ids[i] <- if (length(m)) sub("gene:", "", m) else NA_character_

  # transcript:ENST…
  m <- regmatches(hdr, regexpr("transcript:ENST[0-9]+", hdr))
  enst_ids[i] <- if (length(m)) sub("transcript:", "", m) else NA_character_

  # gene_symbol:XXXX
  m <- regmatches(hdr, regexpr("gene_symbol:[^ ]+", hdr))
  gene_syms[i] <- if (length(m)) sub("gene_symbol:", "", m) else NA_character_
}

# ── 3. Deduplicate by exact protein sequence ──────────────────────────────────
message("Computing sequence hashes…")
seq_md5   <- vapply(seqs, digest, character(1), algo = "md5", USE.NAMES = FALSE)
dedup_idx <- !duplicated(seq_md5)

n_total <- n
n_dedup <- sum(dedup_idx)
message(sprintf("  %d total → %d unique sequences after dedup (%.0f%% reduction)",
                n_total, n_dedup, 100 * (1 - n_dedup / n_total)))

# ── 4. Build lookup tables ────────────────────────────────────────────────────
message("Building lookup tables…")

# Per-sequence metadata (all isoforms)
meta_df <- data.frame(
  ensp     = ensp_ids,
  ensg     = ensg_ids,
  enst     = enst_ids,
  gene_sym = gene_syms,
  seq_md5  = seq_md5,
  stringsAsFactors = FALSE
)

# md5 → all unique ENSG IDs (a single sequence may be shared by multiple genes)
md5_to_ensg <- tapply(
  meta_df$ensg, meta_df$seq_md5,
  function(x) unique(x[!is.na(x)]),
  simplify = FALSE
)
# md5 → all unique gene symbols
md5_to_sym <- tapply(
  meta_df$gene_sym, meta_df$seq_md5,
  function(x) unique(x[!is.na(x)]),
  simplify = FALSE
)
# md5 → all ENSP (for BLAST sseqid back-mapping)
md5_to_ensp <- tapply(
  meta_df$ensp, meta_df$seq_md5,
  function(x) unique(x[!is.na(x)]),
  simplify = FALSE
)

# Representative ENSP for each unique sequence (first encountered)
rep_ensp <- ensp_ids[dedup_idx]
rep_ensg <- ensg_ids[dedup_idx]
rep_sym  <- gene_syms[dedup_idx]
rep_md5  <- seq_md5[dedup_idx]
rep_seq  <- seqs[dedup_idx]

# Also build ENSP → md5 reverse map (for BLAST hit back-lookup)
ensp_to_md5 <- setNames(meta_df$seq_md5, meta_df$ensp)

index <- list(
  # Named character vector: md5 → sequence (deduplicated)
  seqs         = setNames(rep_seq, rep_md5),
  # Lookup tables keyed by md5
  md5_to_ensg  = md5_to_ensg,
  md5_to_sym   = md5_to_sym,
  md5_to_ensp  = md5_to_ensp,
  # Reverse: ENSP → md5 (for BLAST results where sseqid = ENSP)
  ensp_to_md5  = ensp_to_md5,
  # Per-ENSP representative table for BLAST database header construction
  rep_table    = data.frame(
    md5      = rep_md5,
    ensp     = rep_ensp,
    ensg     = rep_ensg,
    gene_sym = rep_sym,
    stringsAsFactors = FALSE
  )
)

message("Saving index RDS…")
saveRDS(index, INDEX_RDS, compress = "xz")
message("  → ", INDEX_RDS)

# ── 5. Write deduplicated FASTA for blastp database ───────────────────────────
message("Writing deduplicated FASTA (", n_dedup, " sequences)…")
con <- file(DEDUP_FA, "w")
for (i in seq_len(n_dedup)) {
  # sseqid = ENSP (no version), so BLAST hits carry ENSP → app looks up via ensp_to_md5
  header <- paste0(">", rep_ensp[i],
                   " gene:", rep_ensg[i] %||% "unknown",
                   " sym:", rep_sym[i]  %||% "unknown")
  writeLines(header, con)
  s      <- rep_seq[i]
  starts <- seq(1L, nchar(s), 60L)
  for (st in starts) writeLines(substr(s, st, min(st + 59L, nchar(s))), con)
}
close(con)
message("  → ", DEDUP_FA)

# ── 6. Build blastp database ──────────────────────────────────────────────────
message("Building blastp database…")
ret <- system2(MAKEBLASTDB, args = c(
  "-in",    DEDUP_FA,
  "-dbtype", "prot",
  "-out",   BLAST_OUT,
  "-title", "Ensembl_114_pep_dedup"
))
if (ret != 0L) stop("makeblastdb failed with exit code ", ret)
message("  → ", BLAST_OUT, ".*")

# ── 7. Quick runtime benchmark ────────────────────────────────────────────────
message("\nBenchmarking Biostrings vmatchPattern on dedup seqs…")
requireNamespace("Biostrings", quietly = TRUE)
ref_set <- Biostrings::AAStringSet(setNames(rep_seq, rep_md5))
test_pep <- Biostrings::AAString("SIINFEKL")   # canonical 9-mer test
t0 <- proc.time()
m0 <- Biostrings::vmatchPattern(test_pep, ref_set, max.mismatch = 0L, fixed = TRUE)
m1 <- Biostrings::vmatchPattern(test_pep, ref_set, max.mismatch = 1L, fixed = TRUE)
elapsed <- (proc.time() - t0)["elapsed"]
n_exact <- sum(lengths(m0) > 0L)
n_near  <- sum(lengths(m1) > 0L)
message(sprintf(
  "  1 peptide × %d sequences: %.2f s  (exact hits: %d, 1-mm hits: %d)",
  n_dedup, elapsed, n_exact, n_near
))
message(sprintf(
  "  Estimated time for 3 peptides: %.2f s", elapsed * 3))

message("\n--- DONE ---")
message("Run 02_prep_allergen.R next, then 03_prep_annotation.R")
