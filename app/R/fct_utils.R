## Formatting utility functions — no Shiny or reactive dependencies.

# Single-row extraction from a matrix as a plain numeric vector.
mat_row <- function(mat, i) as.numeric(mat[i, , drop = FALSE])

# Compute off-tissue risk for a canonical gene from the loaded GTEx TPM matrix.
# Returns "Safe", "Acceptable", "Borderline", "Critical", or "Unavailable".
gtex_ensg_risk <- function(ensg, gtex_mat, gtex_meta, tpm_threshold = 1) {
  if (is.null(gtex_mat) || is.null(gtex_meta) || is.na(ensg) || !nzchar(ensg %||% ""))
    return("Unavailable")
  if (!ensg %in% rownames(gtex_mat)) return("Unavailable")
  expr    <- mat_row(gtex_mat, ensg)
  samples <- colnames(gtex_mat)
  tissue  <- gtex_meta$tissue_type[match(samples, gtex_meta$sample_id)]
  q3_tis  <- Filter(function(tis) {
    vals <- expr[!is.na(tissue) & tissue == tis]
    q3   <- quantile(vals, 0.75, na.rm = TRUE)
    !is.na(q3) && q3 > tpm_threshold
  }, unique(tissue[!is.na(tissue)]))
  # No tissues above threshold → Safe for cross-reactivity purposes.
  # (off_tissue_risk() with tumor_only=FALSE and empty tissues returns "Acceptable"
  #  by design for ncORF candidates, but canonical genes with no Q3>1 tissue are Safe.)
  if (length(q3_tis) == 0L) return("Safe")
  off_tissue_risk(tumor_only = FALSE, tissues_q3_gt1 = paste(q3_tis, collapse = "|"))
}

# Colored circle-exclamation icon for off-tissue risk; "" for Safe / Unavailable.
risk_icon_html <- function(risk) {
  col <- switch(risk,
    "Critical"   = "#b33e3e",
    "Borderline" = "#f3c677",
    "Acceptable" = "#37a4a2",
    NULL
  )
  if (is.null(col)) return("")
  sprintf(
    '<i class="fa-solid fa-circle-exclamation ms-1" style="color:%s" title="%s off-tissue risk (GTEx Q3 &gt; 1 TPM)"></i>',
    col, risk
  )
}

fmt1     <- function(x) if (is.na(x) || length(x) == 0) "—" else sprintf("%.1f%%", x)
fmt2     <- function(x) if (is.na(x) || length(x) == 0) "—" else sprintf("%.3f", x)
bool_fmt <- function(x) if (is.na(x) || length(x) == 0) "—" else if (isTRUE(x)) "Yes ✓" else "No"

# Annotate the Query lines in a blastp -outfmt 0 alignment text so that amino
# acid positions covered by MS peptide matches are rendered bold.
# Returns an HTML string safe for use with htmltools::HTML() inside <pre>.
highlight_blast_query <- function(aln_text, peptides, protein_seq) {
  if (is.null(aln_text) || !nzchar(aln_text %||% "")) return("")
  esc <- function(x) htmltools::htmlEscape(x)

  # Bail out early if no peptide/sequence info — return plain-escaped text
  if (!length(peptides) || !nzchar(protein_seq %||% ""))
    return(esc(aln_text))

  # Build a logical vector: covered[i] == TRUE if position i (1-based) in
  # protein_seq is spanned by at least one matched peptide.
  n_aa    <- nchar(protein_seq)
  covered <- logical(n_aa)
  for (pep in peptides) {
    pep <- trimws(pep)
    if (!nzchar(pep)) next
    m <- gregexpr(pep, protein_seq, fixed = TRUE)[[1]]
    if (m[1L] != -1L) {
      plen <- nchar(pep)
      for (s in as.integer(m)) covered[s:(s + plen - 1L)] <- TRUE
    }
  }

  if (!any(covered)) return(esc(aln_text))

  # Process each line; only Query lines get special treatment.
  lines <- strsplit(aln_text, "\n", fixed = TRUE)[[1]]
  processed <- vapply(lines, function(line) {
    parts <- regmatches(line,
      regexec("^(Query\\s+)(\\d+)(\\s+)([A-Z*-]+)(\\s*)([0-9]*)\\s*$", line))[[1]]
    if (!length(parts)) return(esc(line))

    start <- as.integer(parts[3L])
    chars <- strsplit(parts[5L], "")[[1]]
    pos   <- start
    seq_html <- paste(vapply(chars, function(ch) {
      if (ch == "-") return(esc(ch))
      bold <- if (pos <= n_aa) covered[pos] else FALSE
      pos <<- pos + 1L
      if (bold) paste0("<strong>", esc(ch), "</strong>") else esc(ch)
    }, character(1L)), collapse = "")

    paste0(esc(parts[2L]), esc(parts[3L]), esc(parts[4L]),
           seq_html, esc(parts[6L]), esc(parts[7L]))
  }, character(1L))

  paste(processed, collapse = "\n")
}
