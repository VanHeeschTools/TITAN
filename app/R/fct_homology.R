## Shared peptide-vs-Ensembl-114-proteome homology scan.
## Used by both the ORF Detail cross-reactivity card (app.R) and the Peptide
## tab's per-row homology check (app.R). Depends on the static reference
## objects built once in global.R: ensembl_pep_index, ensembl_pep_seqs_lazy(),
## ensembl_gene_annot — safe to reference directly here since none of them
## are reactive.
##
## Restricted to 0- and 1-mismatch hits only (a 2-mismatch hit was judged too
## weak a signal to act on for cross-reactivity purposes) — this is the one
## piece of logic both call sites share, so changing it here changes it
## everywhere at once.

# Builds the Biostrings reference (AAStringSet + md5 names) once, so a batch
# of independent single-peptide scan_peptide_homology() calls (e.g. the
# Peptide tab precomputing homology for every row on tab entry) can reuse it
# instead of rebuilding it from ~60-80K sequences on every call. Returns
# list(ref_set=, ref_md5=) or data.frame(Error=...) if unavailable - same
# "Error" convention as scan_peptide_homology() so callers can check either
# result the same way before looping.
build_pep_reference <- function() {
  if (is.null(ensembl_pep_index)) {
    return(data.frame(Error = "Ensembl 114 pep index not loaded — run scripts/reference_prep/01_prep_ensembl_pep.R first."))
  }
  pep_seqs <- ensembl_pep_seqs_lazy()
  if (is.null(pep_seqs) || length(pep_seqs) == 0L) {
    return(data.frame(Error = "Ensembl 114 pep sequences could not be loaded."))
  }
  list(ref_set = Biostrings::AAStringSet(pep_seqs), ref_md5 = names(pep_seqs))
}

# Scans `peptides` (character vector) against the Ensembl 114 proteome.
# Gene-level dedup: keeps the best (lowest) Mismatches per ENSG, across all
# input peptides combined. Does NOT self-exclude the caller's own gene(s) -
# callers filter that out afterwards, since "self" differs per call site
# (the ORF being viewed vs. the gene a Peptide-tab row is attributed to).
# `rna_mat` (optional, e.g. rna_tpm_rv()) disambiguates gene symbols that
# resolve to multiple ENSGs, keeping only ENSGs present in that matrix.
# `ref` (optional) is a pre-built build_pep_reference() result, to skip
# rebuilding the AAStringSet on every call in a batch loop - if omitted, it's
# built fresh internally (unchanged behaviour for existing single-shot callers).
# Returns a data.frame(Peptide, Query_html, Target_html, Gene_sym, ENSG,
# Mismatches) — possibly zero rows — or data.frame(Error = <message>) if the
# reference index/sequences aren't available.
scan_peptide_homology <- function(peptides, rna_mat = NULL, ref = NULL) {
  peptides <- unique(peptides)
  if (length(peptides) == 0L) return(data.frame())
  if (is.null(ref)) ref <- build_pep_reference()
  if ("Error" %in% names(ref)) return(ref)
  ref_set <- ref$ref_set
  ref_md5 <- ref$ref_md5

  # Render a peptide as HTML, bolding the given (1-based) mismatch positions.
  pep_html <- function(chars, mm_pos) {
    paste(ifelse(seq_along(chars) %in% mm_pos,
                 paste0("<strong>", chars, "</strong>"), chars),
          collapse = "")
  }

  hits <- do.call(rbind, Filter(Negate(is.null), lapply(peptides, function(pep) {
    pep_aa    <- tryCatch(Biostrings::AAString(pep), error = function(e) NULL)
    if (is.null(pep_aa)) return(NULL)
    pep_len   <- nchar(pep)
    pep_chars <- strsplit(pep, "")[[1]]

    m0 <- tryCatch(Biostrings::vmatchPattern(pep_aa, ref_set, max.mismatch = 0L, fixed = TRUE),
                   error = function(e) NULL)
    m1 <- tryCatch(Biostrings::vmatchPattern(pep_aa, ref_set, max.mismatch = 1L, fixed = TRUE),
                   error = function(e) NULL)

    # Strict per-level indices (parallel to ref_set rows)
    idx0 <- if (!is.null(m0)) which(lengths(m0) > 0L) else integer(0)
    idx1 <- if (!is.null(m1)) setdiff(which(lengths(m1) > 0L), idx0) else integer(0)
    if (!length(idx0) && !length(idx1)) return(NULL)

    # For each matched sequence index, extract the target subsequence and
    # build HTML-rendered query/target strings with mismatches bolded.
    build_rows <- function(seq_idx, views_list, mm) {
      if (!length(seq_idx)) return(NULL)
      do.call(rbind, Filter(Negate(is.null), lapply(seq_idx, function(i) {
        views <- views_list[[i]]
        if (!length(views)) return(NULL)
        s <- IRanges::start(views)[1L]
        tgt <- tryCatch(
          as.character(Biostrings::subseq(ref_set[[i]], start = s, width = pep_len)),
          error = function(e) NULL
        )
        if (is.null(tgt) || nchar(tgt) != pep_len) return(NULL)
        tgt_chars   <- strsplit(tgt, "")[[1]]
        mm_pos      <- which(pep_chars != tgt_chars)
        query_html  <- pep_html(pep_chars, mm_pos)
        target_html <- pep_html(tgt_chars, mm_pos)
        md5         <- ref_md5[i]
        ensg_vec    <- ensembl_pep_index$md5_to_ensg[[md5]] %||% NA_character_
        sym_vec     <- ensembl_pep_index$md5_to_sym[[md5]]  %||% "unknown"
        # Iterate over ENSGs explicitly to avoid length-mismatch in data.frame()
        # when md5_to_ensg and md5_to_sym have different lengths.
        do.call(rbind, lapply(seq_along(ensg_vec), function(j) {
          ensg_j <- as.character(ensg_vec[[j]])[1L]
          gene_sym_j <- if (!is.null(ensembl_gene_annot) && !is.na(ensg_j)) {
            idx_a <- match(ensg_j, ensembl_gene_annot$ensembl_gene_id)
            if (!is.na(idx_a)) ensembl_gene_annot$external_gene_name[idx_a]
            else if (length(sym_vec) > 0L) as.character(sym_vec[1L]) else NA_character_
          } else {
            if (length(sym_vec) > 0L) as.character(sym_vec[1L]) else NA_character_
          }
          data.frame(
            Peptide     = pep,
            Query_html  = query_html,
            Target_html = target_html,
            Gene_sym    = gene_sym_j,
            ENSG        = ensg_j,
            Mismatches  = mm,
            stringsAsFactors = FALSE
          )
        }))
      })))
    }

    rbind(build_rows(idx0, m0, 0L),
          build_rows(idx1, m1, 1L))
  })))

  if (is.null(hits) || nrow(hits) == 0L) return(data.frame())

  # Gene-level dedup: keep best (lowest) Mismatches per gene, across all
  # input peptides combined.
  hits <- hits[order(hits$ENSG, hits$Mismatches), ]
  hits <- hits[!duplicated(hits$ENSG), ]

  # For gene names that resolve to multiple ENSGs, keep only those present
  # in the tumor quantification matrix.
  if (nrow(hits) > 0L && !is.null(rna_mat)) {
    gene_n     <- table(hits$Gene_sym)
    multi_syms <- names(gene_n)[gene_n > 1L]
    if (length(multi_syms) > 0L) {
      in_mat   <- hits$ENSG %in% rownames(rna_mat)
      is_multi <- hits$Gene_sym %in% multi_syms
      hits     <- hits[!is_multi | in_mat, , drop = FALSE]
    }
  }

  hits
}

# Drops rows in `hits` (as returned by scan_peptide_homology()) whose ENSG or
# Gene_sym matches the given self gene - shared by both call sites so the
# exact self-exclusion rule (match by either ENSG or gene name, since some
# hits carry NA ENSG) can't drift between them.
exclude_self_gene <- function(hits, self_ensg, self_gene_name) {
  if (is.null(hits) || nrow(hits) == 0L || "Error" %in% names(hits)) return(hits)
  if (!is.na(self_ensg) && nzchar(self_ensg %||% ""))
    hits <- hits[is.na(hits$ENSG) | hits$ENSG != self_ensg, , drop = FALSE]
  if (!is.na(self_gene_name) && nzchar(self_gene_name %||% ""))
    hits <- hits[is.na(hits$Gene_sym) | hits$Gene_sym != self_gene_name, , drop = FALSE]
  hits
}
