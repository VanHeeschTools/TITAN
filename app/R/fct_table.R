## Priority table cell HTML builders.
## Depends on: biotype_badge_html, spec_badge_html, score_bar_html, pct_bar_html (fct_scoring.R).
## Depends on: risk_icon_html, gtex_ensg_risk (fct_utils.R); exclude_self_gene (fct_homology.R).

# KPI stat box (Overview / Prioritization stat rows). Built as plain HTML
# rather than bslib::value_box()'s showcase/showcase_layout mechanism -
# that mechanism's default layout appears to vary across bslib versions and
# wasn't reliably producing the icon-left horizontal layout on this app's
# installed version, so this sidesteps it entirely with markup we fully
# control (see titan.css ".titan-kpi-box" for styling).
titan_kpi_box <- function(icon_name, label, value_output_id, unit = NULL) {
  tags$div(class = "titan-kpi-box",
    tags$div(class = "titan-kpi-icon", icon(icon_name)),
    tags$div(class = "titan-kpi-text",
      tags$div(class = "titan-kpi-label", label),
      tags$div(class = "titan-kpi-value",
        textOutput(value_output_id, inline = TRUE),
        if (!is.null(unit)) tags$span(class = "titan-kpi-unit", unit))
    )
  )
}

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

# Vectorised over the WHOLE items_str column in one pass instead of one R call
# (strsplit + nested vapply) per row - at thousands of matched-peptide-group
# rows the old per-row make_peptide_cell() was the one column in this file
# that never got the same treatment as score_bar_html/spec_badge_html/etc.
# (see prio_table_df() comment in app.R), and dominated table build time.
make_peptide_cell <- function(items_str) {
  items_str <- as.character(items_str)
  # "" / NA marks an unmatched ORF (no peptide evidence at all - see the
  # prio_include_unmatched path in prioritised_data(), app.R). strsplit("")
  # behavior isn't relied on here; handled explicitly below instead.
  is_empty <- is.na(items_str) | !nzchar(trimws(items_str))
  mono <- function(s) sprintf('<span class="font-monospace" style="font-size:10px">%s</span>', s)

  split_list <- strsplit(items_str, ",\\s*")
  n_items    <- lengths(split_list)
  first_item <- vapply(split_list, function(x) if (length(x)) x[1] else NA_character_, character(1))
  out        <- mono(trimws(first_item))

  multi <- which(n_items > 1L)
  if (length(multi) > 0L) {
    extra_items <- lapply(split_list[multi], function(x) trimws(x[-1]))
    extra_rep   <- rep(seq_along(multi), lengths(extra_items))
    extra_flat  <- unlist(extra_items, use.names = FALSE)
    extra_html  <- mono(extra_flat)
    # factor(levels=) keeps rows in `multi` order - split()'s default
    # character-sort on unfactored group keys would otherwise scramble
    # order past 9 groups (e.g. "10" sorting before "2").
    by_row  <- split(extra_html, factor(extra_rep, levels = seq_along(multi)))
    all_html <- vapply(by_row, paste, character(1), collapse = "<br>")
    n_more   <- n_items[multi] - 1L
    out[multi] <- paste0(
      out[multi],
      ' <span class="titan-pep-more">and ', n_more, ' more...</span>',
      '<span class="titan-pep-less" style="display:none">less</span>',
      '<div class="titan-pep-extra" style="display:none; margin-top:3px; line-height:1.7">',
      all_html,
      '</div>'
    )
  }
  out[is_empty] <- '<span class="text-muted fst-italic" style="font-size:11px">No peptide evidence</span>'
  out
}

make_child_rows_html <- function(orfs_df) {
  # Returns a character VECTOR, one <tr class="titan-child-row"> string per
  # row of orfs_df (NOT concatenated/grouped - caller splits+collapses per
  # gene/biotype/peptide group). orfs_df: non-best ORF rows across ALL
  # multi-ORF groups at once (already sorted desc by score within each
  # group), with orf_biotype_single and matched_peptides restored from the
  # grouping key.
  #
  # Vectorised across the WHOLE input in one pass instead of looping per row
  # via orfs_df[i, ] + scalar badge/bar calls (the original approach): at
  # ~10k child rows (rms_organoids, real data) the old per-row version - each
  # iteration paying a data.frame row-extraction plus a full sprintf/ifelse
  # call for what are otherwise vectorised helpers - measured at 23s
  # standalone (profiling/child_html_isolate.R). This version does the exact
  # same badge/bar/number-formatting helpers but ONCE, over the full column,
  # matching how prio_table_df()'s top-level columns already work.
  #
  # Cell layout must match prio_table_df() transmute column order (28 cols total):
  # Sel(0) Gene(1) ORF-biotype(2) Peptides(3) ORF-id(4) Location(5) Spec(6) Off-tissue(7) Score(8)
  # Transl%(9) TranslPPM(10) Expr%(11) ExprTPM(12) GTEx(13) TCGAT%(14) TCGATPM(15)
  # TCGAN%(16) TCANPM(17) RCprim%(18) RCprimPPM(19) RCCL%(20) RCCLPPM(21)
  # .biotype_sort(22) .spec_sort(23) .off_tissue_sort(24) .score_sort(25) .transl_sort(26) .expr_sort(27) .child_rows(28)
  if (nrow(orfs_df) == 0L) return(character(0))
  r2 <- function(x) ifelse(is.na(x) | !is.finite(x), "&mdash;", sprintf("%.2f", x))
  r1 <- function(x) ifelse(is.na(x) | !is.finite(x), "&mdash;", sprintf("%.1f", x))
  r3 <- function(x) ifelse(is.na(x) | !is.finite(x), "&mdash;", sprintf("%.3f", x))
  pep_html <- make_peptide_cell(orfs_df$matched_peptides)
  loc_html <- paste0(
    orfs_df$chr, ':', formatC(orfs_df$orf_start, format = "d", big.mark = ","),
    '&ndash;', formatC(orfs_df$orf_end, format = "d", big.mark = ","),
    ' ', orfs_df$strand, ' ', orfs_df$start_codon
  )
  paste0(
    '<tr class="titan-child-row">',
    '<td class="dt-center titan-sel-col"></td>',
    '<td></td>',
    '<td>', biotype_badge_html(orfs_df$orf_biotype_single), '</td>',
    '<td class="titan-pep-cell">', pep_html, '</td>',
    '<td><span class="font-monospace" style="font-size:10px;word-break:break-all">',
      orfs_df$orf_id, '</span></td>',
    '<td style="font-size:11px;white-space:nowrap">', loc_html, '</td>',
    '<td>', spec_badge_html(orfs_df$GTEX_tumor_only, orfs_df$GTEX_tumor_enriched), '</td>',
    '<td>', off_tissue_risk_html(orfs_df$off_tissue_label), '</td>',
    '<td>', score_bar_html(orfs_df$priority_score), '</td>',
    '<td>', pct_bar_html(orfs_df$target_translation_pct_samples, "#2F3D46"), '</td>',
    '<td style="font-size:12px">', r2(orfs_df$target_translation_median_PPM), '</td>',
    '<td>', pct_bar_html(orfs_df$target_expression_pct_samples, "#7EB8BF"), '</td>',
    '<td style="font-size:12px">', r2(orfs_df$target_expression_median_TPM), '</td>',
    '<td style="font-size:12px">', r3(orfs_df$GTEX_max_median_TPM), '</td>',
    '<td style="font-size:12px">', r1(orfs_df$TCGA_tumor_pct_samples), '</td>',
    '<td style="font-size:12px">', r2(orfs_df$TCGA_tumor_median_TPM), '</td>',
    '<td style="font-size:12px">', r1(orfs_df$TCGA_normal_pct_samples), '</td>',
    '<td style="font-size:12px">', r2(orfs_df$TCGA_normal_median_TPM), '</td>',
    '<td style="font-size:12px">', r1(orfs_df$ribocrypt_primary_pct_samples), '</td>',
    '<td style="font-size:12px">', r2(orfs_df$ribocrypt_primary_median_PPM), '</td>',
    '<td style="font-size:12px">', r1(orfs_df$`ribocrypt_cell-line_pct_samples`), '</td>',
    '<td style="font-size:12px">', r2(orfs_df$`ribocrypt_cell-line_median_PPM`), '</td>',
    '<td></td><td></td><td></td><td></td><td></td><td></td>',
    '</tr>'
  )
}

# Peptide tab's "Homology (0/1mm)" cell — one row at a time (this table is
# scoped to just the currently-selected candidates, so it's small; not worth
# the vectorisation effort the priority-table cells need at thousands of
# rows). `cache` is pep_xreact_cache_rv()'s list (keyed by peptide sequence).
# `rowid` is the row's .pep_row_id, used to route the "Check" click back to
# the right peptide via data attributes.
make_pep_homology_cell <- function(peptide, gene_id_clean, gene_name, rowid, cache,
                                    gtex_mat, gtex_meta) {
  hits <- cache[[peptide]]
  if (is.null(hits)) {
    return(sprintf(
      '<a href="#" class="titan-pep-homology-check" data-peptide="%s" data-rowid="%s">Check &#8981;</a>',
      htmltools::htmlEscape(peptide), htmltools::htmlEscape(rowid)
    ))
  }
  if ("Error" %in% names(hits))
    return(sprintf('<span class="text-danger small">%s</span>', htmltools::htmlEscape(hits$Error[1L])))

  hits <- exclude_self_gene(hits, gene_id_clean, gene_name)
  if (is.null(hits) || nrow(hits) == 0L)
    return('<span class="text-success small"><i class="fa-solid fa-circle-check"></i> None</span>')

  hits  <- hits[order(hits$Mismatches, hits$Gene_sym), ]
  icons <- vapply(hits$ENSG, function(e) risk_icon_html(gtex_ensg_risk(e, gtex_mat, gtex_meta)), character(1))
  items <- sprintf('%s (%dmm)%s', hits$Gene_sym, hits$Mismatches, icons)

  # "first item + N more..." toggle, same convention as make_peptide_cell()'s
  # .titan-pep-more/.titan-pep-less/.titan-pep-extra (tbl_peptide's own JS
  # callback binds the same handler, scoped to #tbl_peptide).
  out <- items[1L]
  if (length(items) > 1L) {
    out <- paste0(
      out,
      ' <span class="titan-pep-more">and ', length(items) - 1L, ' more...</span>',
      '<span class="titan-pep-less" style="display:none">less</span>',
      '<div class="titan-pep-extra" style="display:none; margin-top:3px; line-height:1.7">',
      paste(items[-1L], collapse = "<br>"),
      '</div>'
    )
  }
  out
}
