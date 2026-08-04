## Shared parser for the TransCode Phase 2 Gencode ORF reference table.
## Used by: app/global.R (MS-peptide cross-matching against gencode_orf_tbl) and
##          scripts/catalog_prep/build_candidate_backbone.R (TITAN-FORGE ncorfs backbone
##          for studies with no internal ribo-seq ORF calling, e.g. TCGA PAAD).

load_gencode_orf_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(
    data.table::fread(path, data.table = FALSE, showProgress = FALSE),
    error = function(e) { message("Gencode ORF table not loaded: ", e$message); NULL }
  )
  if (is.null(df)) return(NULL)

  starts_raw <- df[["starts (0-based)"]]
  ends_raw   <- df[["ends (0-based)"]]
  # Multi-exon rows are semicolon-joined; overall genomic span is min(starts)/max(ends).
  span <- function(x, f) vapply(strsplit(x, ";", fixed = TRUE),
                                 function(v) f(as.integer(v)), numeric(1))

  aa <- sub("\\*$", "", df$sequence_aa_MS_ms)
  nt <- df$sequence_nt

  df$orf_id             <- df$releasev45_id
  df$orf_biotype_single <- ifelse(df$orf_type == "PT", "Processed_transcript_ORF", df$orf_type)
  df$orf_biotypes_all   <- df$orf_biotype_single
  df$protein_seq        <- aa
  df$protein_length     <- nchar(aa)
  df$chr                <- df$chrm
  df$starts              <- starts_raw
  df$ends                <- ends_raw
  df$orf_start           <- span(starts_raw, min)
  df$orf_end             <- span(ends_raw,   max)
  df$start_codon         <- df$initiation_codon
  df$stop_codon          <- substr(nt, nchar(nt) - 2L, nchar(nt))
  df$tx_id               <- df$transcript
  df$caller_count        <- NA_integer_
  df$gene_id_clean       <- sub("\\..*", "", df$gene_id)
  df$summary_id          <- paste(df$orf_id, df$chr, df$orf_start, df$orf_end,
                                   df$strand, df$orf_biotype_single, sep = "_")

  df
}
