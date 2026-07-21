## Priority table cell HTML builders.
## Depends on: biotype_badge_html, spec_badge_html, score_bar_html, pct_bar_html (fct_scoring.R).

make_expand_cell <- function(count, items_str) {
  items <- trimws(strsplit(as.character(items_str), ",\\s*")[[1]])
  if (length(items) <= 1L) return(as.character(count))
  items_html <- paste(items, collapse = "<br>")
  paste0(
    count,
    ' <span class="titan-expand-btn">+</span>',
    '<div class="titan-xcontent" style="display:none">',
    items_html,
    '</div>'
  )
}

make_peptide_cell <- function(items_str) {
  items <- trimws(strsplit(as.character(items_str), ",\\s*")[[1]])
  mono  <- function(s) sprintf('<span class="font-monospace" style="font-size:10px">%s</span>', s)
  if (length(items) <= 1L) return(mono(items[1]))
  n_more   <- length(items) - 1L
  all_html <- paste(vapply(items[-1], mono, character(1)), collapse = "<br>")
  paste0(
    mono(items[1]),
    ' <span class="titan-pep-more">and ', n_more, ' more...</span>',
    '<span class="titan-pep-less" style="display:none">less</span>',
    '<div class="titan-pep-extra" style="display:none; margin-top:3px; line-height:1.7">',
    all_html,
    '</div>'
  )
}

make_child_html <- function(orfs_df) {
  # Returns concatenated <tr class="titan-child-row"> strings, one per ORF.
  # orfs_df: non-best ORF rows from the group (already sorted desc by score),
  # with orf_biotype_single and matched_peptides restored from the grouping key.
  # Cell layout must match prio_table_df() transmute column order (27 cols total):
  # Sel(0) Gene(1) ORF-biotype(2) Peptides(3) ORF-id(4) Location(5) Spec(6) Score(7)
  # Transl%(8) TranslPPM(9) Expr%(10) ExprTPM(11) GTEx(12) TCGAT%(13) TCGATPM(14)
  # TCGAN%(15) TCANPM(16) RCprim%(17) RCprimPPM(18) RCCL%(19) RCCLPPM(20)
  # .biotype_sort(21) .spec_sort(22) .score_sort(23) .transl_sort(24) .expr_sort(25) .child_rows(26)
  r2 <- function(x) if (is.na(x) || !is.finite(x)) "&mdash;" else sprintf("%.2f", x)
  r1 <- function(x) if (is.na(x) || !is.finite(x)) "&mdash;" else sprintf("%.1f", x)
  r3 <- function(x) if (is.na(x) || !is.finite(x)) "&mdash;" else sprintf("%.3f", x)
  rows <- vapply(seq_len(nrow(orfs_df)), function(i) {
    r <- orfs_df[i, ]
    paste0(
      '<tr class="titan-child-row">',
      '<td class="dt-center titan-sel-col"></td>',
      '<td></td>',
      '<td>', biotype_badge_html(r$orf_biotype_single), '</td>',
      '<td class="titan-pep-cell">', make_peptide_cell(r$matched_peptides), '</td>',
      '<td><span class="font-monospace" style="font-size:10px;word-break:break-all">',
        r$orf_id, '</span></td>',
      '<td style="font-size:11px;white-space:nowrap">',
        r$chr, ':', formatC(r$orf_start, format = "d", big.mark = ","),
        '&ndash;', formatC(r$orf_end, format = "d", big.mark = ","),
        ' ', r$strand, ' ', r$start_codon, '</td>',
      '<td>', spec_badge_html(r$GTEX_tumor_only, r$GTEX_tumor_enriched), '</td>',
      '<td>', score_bar_html(r$priority_score), '</td>',
      '<td>', pct_bar_html(r$target_translation_pct_samples, "#28646E"), '</td>',
      '<td style="font-size:12px">', r2(r$target_translation_median_PPM), '</td>',
      '<td>', pct_bar_html(r$target_expression_pct_samples, "#7EB8BF"), '</td>',
      '<td style="font-size:12px">', r2(r$target_expression_median_TPM), '</td>',
      '<td style="font-size:12px">', r3(r$GTEX_max_median_TPM), '</td>',
      '<td style="font-size:12px">', r1(r$TCGA_tumor_pct_samples), '</td>',
      '<td style="font-size:12px">', r2(r$TCGA_tumor_median_TPM), '</td>',
      '<td style="font-size:12px">', r1(r$TCGA_normal_pct_samples), '</td>',
      '<td style="font-size:12px">', r2(r$TCGA_normal_median_TPM), '</td>',
      '<td style="font-size:12px">', r1(r$ribocrypt_primary_pct_samples), '</td>',
      '<td style="font-size:12px">', r2(r$ribocrypt_primary_median_PPM), '</td>',
      '<td style="font-size:12px">', r1(r$`ribocrypt_cell-line_pct_samples`), '</td>',
      '<td style="font-size:12px">', r2(r$`ribocrypt_cell-line_median_PPM`), '</td>',
      '<td></td><td></td><td></td><td></td><td></td><td></td>',
      '</tr>'
    )
  }, character(1))
  paste(rows, collapse = "")
}
