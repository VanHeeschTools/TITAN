## Report HTML page builders — pure functions, no Shiny reactive context required.
## Depends on: report_plots.R, fct_utils.R, BIOTYPE_COLORS (global.R).

.rpt_seq_html <- function(seq, peps) {
  # Simplified version of render_protein_seq_html with no Bootstrap/popover
  if (is.na(seq) || !nzchar(seq)) return("<em style='color:black'>No sequence available.</em>")
  PEP_COLS <- c("#28646E", "#D4850A", "#8E44AD", "#C0392B", "#0097A7")
  peps     <- unique(peps[!is.na(peps) & nzchar(peps)])
  n        <- nchar(seq)
  chars    <- strsplit(seq, "")[[1]]
  coverage <- integer(n)
  for (pi in seq_along(peps)) {
    m <- gregexpr(peps[[pi]], seq, fixed = TRUE)[[1]]
    if (m[1L] > 0L) {
      plen <- nchar(peps[[pi]])
      for (s in m) for (j in s:min(s + plen - 1L, n)) if (!coverage[j]) coverage[j] <- pi
    }
  }
  BLOCK <- 59L
  blocks <- vapply(seq_len(ceiling(n / BLOCK)), function(b) {
    i0 <- (b - 1L) * BLOCK + 1L; i1 <- min(b * BLOCK, n)
    cv <- coverage[i0:i1]; sv <- chars[i0:i1]
    rl  <- rle(cv); r_len <- rl$lengths; r_val <- rl$values
    r_end <- cumsum(r_len); r_st <- c(1L, r_end[-length(r_end)] + 1L)
    parts <- mapply(function(rs, re, cvi) {
      txt <- paste(sv[rs:re], collapse = "")
      if (cvi == 0L) return(sprintf('<span style="color:#444">%s</span>', txt))
      col <- PEP_COLS[(cvi - 1L) %% 5L + 1L]
      sprintf('<span style="background:%s33;color:%s;font-weight:bold" title="%s">%s</span>',
              col, col, peps[[cvi]], txt)
    }, r_st, r_end, r_val, SIMPLIFY = TRUE)
    sprintf('<div><span style="color:#aaa;margin-right:4px">%4d</span>%s</div>',
            i0, paste(parts, collapse = ""))
  }, character(1))
  legend <- if (length(peps))
    paste(vapply(seq_along(peps), function(pi) {
      col <- PEP_COLS[(pi - 1L) %% 5L + 1L]
      sprintf('<code style="background:%s22;color:%s;border:1px solid %s55;border-radius:3px;padding:0 3px;font-size:14px">%s</code>',
              col, col, col, peps[[pi]])
    }, character(1)), collapse = " ")
  else '<span style="color:black;font-size:14px">No MS peptides matched.</span>'
  paste0('<div style="font-size:14px;color:black;margin-bottom:4px">', legend, '</div>',
         '<div style="font-family:monospace;font-size:13px;line-height:1.55">',
         paste(blocks, collapse = ""), '</div>')
}

.rpt_badge <- function(tumor_only, tumor_enriched) {
  if (isTRUE(as.logical(tumor_only)))
    '<span style="background:#d4edda;color:#155724;border-radius:4px;padding:2px 6px;font-size:11px;font-weight:600">Tumor-only</span>'
  else if (isTRUE(as.logical(tumor_enriched)))
    '<span style="background:#cce5ff;color:#004085;border-radius:4px;padding:2px 6px;font-size:11px;font-weight:600">Tumor-enriched</span>'
  else
    '<span style="background:#e2e3e5;color:#383d41;border-radius:4px;padding:2px 6px;font-size:11px;font-weight:600">Non-specific</span>'
}

.rpt_risk_badge <- function(tumor_only, tissues_q3_gt1) {
  label <- off_tissue_risk(tumor_only, tissues_q3_gt1, off_tissue_risk_adult)
  col <- switch(label,
    "Safe"        = "#21ae7f",
    "Acceptable"  = "#37a4a2",
    "Borderline"  = "#f3c677",
    "Critical"    = "#b33e3e",
    "Unavailable" = "#adb5bd",
    "#adb5bd"
  )
  txt <- if (identical(label, "Borderline")) "#333333" else "#ffffff"
  sprintf('<span style="background:%s;color:%s;border-radius:4px;padding:2px 6px;font-size:11px;font-weight:600">%s</span>',
          col, txt, label)
}

.rpt_build_page <- function(row, pep_list,
                             rna_mat, rna_meta, gtex_mat, gtex_meta,
                             tcga_mat, tcga_meta, ribo_m, ribo_sm,
                             rc_mat,  rc_meta,
                             logo_uri, vh_logo_uri, gen_date, log_scale) {
  LW <- 10.5; PH_E <- 3.72; PH_T <- 3.36
  apsc    <- if (log_scale) function(x) log(pmax(as.numeric(x), 0) + 1) else function(x) pmax(as.numeric(x), 0)
  ylabel_e <- if (log_scale) "log(TPM+1)" else "TPM"
  ylabel_t <- if (log_scale) "log(PPM+1)" else "PPM"
  y_ref   <- if (log_scale) log(2) else 1

  gid <- row$gene_id_clean
  oid <- row$orf_id

  # Shared y-range for expression group (target TPM + GTEx TPM + TCGA TPM)
  expr_ylim <- .rpt_ylim(c(
    if (!is.null(rna_mat)  && isTRUE(gid %in% rownames(rna_mat)))  list(apsc(as.numeric(rna_mat[gid, ])))  else list(),
    if (!is.null(gtex_mat) && isTRUE(gid %in% rownames(gtex_mat))) list(apsc(as.numeric(gtex_mat[gid, ]))) else list(),
    if (!is.null(tcga_mat) && isTRUE(gid %in% rownames(tcga_mat))) list(apsc(as.numeric(tcga_mat[gid, ]))) else list()
  ))

  # Shared y-range for translation group (ribo PPM + RC PPM)
  transl_ylim <- .rpt_ylim(c(
    if (!is.null(ribo_m) && isTRUE(oid %in% rownames(ribo_m))) list(apsc(as.numeric(ribo_m[oid, ]))) else list(),
    if (!is.null(rc_mat) && isTRUE(oid %in% rownames(rc_mat))) list(apsc(as.numeric(rc_mat[oid, ]))) else list()
  ))

  p_et  <- .rpt_expr_target(row, rna_mat, rna_meta, apsc, ylabel_e, y_ref = y_ref, ylim = expr_ylim)
  p_eg  <- .rpt_expr_gtex(row, gtex_mat, gtex_meta, apsc,           y_ref = y_ref, ylim = expr_ylim)
  p_etc <- .rpt_expr_tcga(row, tcga_mat, tcga_meta, apsc,           y_ref = y_ref, ylim = expr_ylim)
  p_tt  <- .rpt_transl_target(row, ribo_m, ribo_sm, apsc, ylabel_t, y_ref = y_ref, ylim = transl_ylim)
  p_tp  <- .rpt_transl_rc(row, rc_mat, rc_meta, "RC Primary",
                           apsc, label_fn = function(x) sub("^primary_", "", x),
                           y_ref = y_ref, ylim = transl_ylim)
  p_tcl <- .rpt_transl_rc(row, rc_mat, rc_meta, "RC Cell-line", apsc,
                           y_ref = y_ref, ylim = transl_ylim)

  # Align panel boundaries within each row so varying label lengths don't shift panels.
  # Skip alignment if any plot in the row has no data (void placeholder).
  .is_void <- function(p) isTRUE(attr(p, "is_void"))
  .save <- function(p, w, h) .rpt_plot_uri(cowplot::ggdraw(p), w, h)

  al_e <- if (!.is_void(p_et) && !.is_void(p_eg) && !.is_void(p_etc))
    cowplot::align_plots(p_et, p_eg, p_etc, align = "h", axis = "tb")
  else list(p_et, p_eg, p_etc)

  al_t <- if (!.is_void(p_tt) && !.is_void(p_tp) && !.is_void(p_tcl))
    cowplot::align_plots(p_tt, p_tp, p_tcl, align = "h", axis = "tb")
  else list(p_tt, p_tp, p_tcl)

  uri_et  <- .save(al_e[[1]], LW * 0.12, PH_E)
  uri_eg  <- .save(al_e[[2]], LW * 0.53, PH_E)
  uri_etc <- .save(al_e[[3]], LW * 0.35, PH_E)
  uri_tt  <- .save(al_t[[1]], LW * 0.10, PH_T)
  uri_tp  <- .save(al_t[[2]], LW * 0.20, PH_T)
  uri_tcl <- .save(al_t[[3]], LW * 0.70, PH_T)

  seq_html <- .rpt_seq_html(row$protein_seq %||% NA_character_, pep_list)

  bio_col <- unname(BIOTYPE_COLORS[row$orf_biotype_single])
  if (is.na(bio_col)) bio_col <- "#95A5A6"

  tiles <- list(
    list("Expr. %",      fmt1(row$target_expression_pct_samples)),
    list("Median TPM",   fmt2(row$target_expression_median_TPM)),
    list("Transl. %",    fmt1(row$target_translation_pct_samples)),
    list("Median PPM",   fmt2(row$target_translation_median_PPM)),
    list("GTEx max TPM", fmt2(row$GTEX_max_median_TPM)),
    list("TCGA tumor %", fmt1(row$TCGA_tumor_pct_samples)),
    list("TCGA normal %",fmt1(row$TCGA_normal_pct_samples)),
    list("RC primary %", fmt1(row$ribocrypt_primary_pct_samples)),
    list("RC CL %",      fmt1(row$`ribocrypt_cell-line_pct_samples`))
  )
  tile_html <- paste(vapply(tiles, function(t)
    sprintf('<div style="flex:1;border:1px solid #ddd;border-radius:4px;padding:4px 5px;text-align:center;min-width:0">
               <div style="font-size:10px;color:black;white-space:nowrap;overflow:hidden;text-overflow:ellipsis">%s</div>
               <div style="font-size:13px;font-weight:600;color:#1a1a1a">%s</div></div>', t[[1]], t[[2]]),
    character(1)), collapse = "")

  n_orfs_str <- sprintf('%s ORF%s: %s', row$n_orfs, if (isTRUE(row$n_orfs == 1L)) "" else "s",
                        row$orf_ids)
  n_pep_str  <- sprintf('%s peptide%s: %s', row$n_peptides, if (isTRUE(row$n_peptides == 1L)) "" else "s",
                        row$matched_peptides)

  sprintf('
<div class="rpt-slide">
  <div class="rpt-inner">
    <div style="display:flex;justify-content:space-between;align-items:flex-start;margin-bottom:5pt">
      <div>
        <div style="display:flex;align-items:baseline;gap:6pt;flex-wrap:wrap">
          <span style="font-size:17pt;font-weight:700;color:#1a1a1a;line-height:1.1">%s</span>
          <span style="background:%s22;color:#111;border:1px solid %s55;border-radius:4px;padding:1px 6px;font-size:10px;font-weight:600">%s</span>
          %s
        </div>
        <div style="font-size:9px;color:#555;margin-top:3px;word-break:break-all">%s</div>
        <div style="font-size:9px;color:#555;word-break:break-all">%s</div>
      </div>
      <div style="text-align:right;flex-shrink:0;margin-left:12pt">
        <div style="font-size:20pt;font-weight:700;color:#28646E;line-height:1">%.1f</div>
        <div style="font-size:9px;color:black">/ 100</div>
      </div>
    </div>
    <div style="border-top:1px solid #e0e0e0;margin-bottom:5pt"></div>
    <div style="display:flex;gap:4pt;margin-bottom:6pt">%s</div>
    <div style="display:flex;gap:0.18in">
      <div style="flex:0 0 62%%;min-width:0">
        <div style="font-size:10.5pt;font-weight:700;color:#444;margin-bottom:2pt">
          Expression (%s)
        </div>
        <div style="display:flex;width:100%%;height:%.2fin">
          <img src="%s" style="width:12%%;height:100%%;object-fit:fill">
          <img src="%s" style="width:53%%;height:100%%;object-fit:fill">
          <img src="%s" style="width:35%%;height:100%%;object-fit:fill">
        </div>
        <div style="border-top:1px solid #eee;margin:4pt 0 3pt"></div>
        <div style="font-size:10.5pt;font-weight:700;color:#444;margin-bottom:2pt">
          Translation (%s)
        </div>
        <div style="display:flex;width:100%%;height:%.2fin">
          <img src="%s" style="width:10%%;height:100%%;object-fit:fill">
          <img src="%s" style="width:20%%;height:100%%;object-fit:fill">
          <img src="%s" style="width:70%%;height:100%%;object-fit:fill">
        </div>
      </div>
      <div style="flex:1;min-width:0;padding-left:8pt;border-left:1px solid #eee">
        <div style="font-size:10.5pt;font-weight:700;color:#444;margin-bottom:4pt">
          Protein sequence &amp; MS peptide matches
        </div>
        <div style="overflow:hidden">%s</div>
      </div>
    </div>
  </div>
  <div class="rpt-footer">
    <img src="%s" style="height:24px">
    <span>Generated: %s</span>
    <img src="%s" style="height:24px">
  </div>
</div>',
    row$gene_name, bio_col, bio_col, row$orf_biotype_single,
    paste0(.rpt_badge(row$GTEX_tumor_only, row$GTEX_tumor_enriched), " ",
           .rpt_risk_badge(row$GTEX_tumor_only, row$GTEX_tissues_q3_gt1)),
    n_orfs_str, n_pep_str, row$priority_score,
    tile_html,
    ylabel_e, PH_E, uri_et, uri_eg, uri_etc,
    ylabel_t, PH_T, uri_tt, uri_tp, uri_tcl,
    seq_html, logo_uri, gen_date, vh_logo_uri
  )
}

.rpt_build_summary_page <- function(params_list, data_info, logo_uri, vh_logo_uri, gen_date) {
  weight_rows <- paste(vapply(params_list$weights, function(r)
    sprintf('<tr><td style="padding:2px 8px">%s</td><td style="padding:2px 8px;text-align:right;color:%s">%+.2f</td></tr>',
            r$label, if (r$weight > 0) "#00A555" else if (r$weight < 0) "#C0392B" else "black", r$weight),
    character(1)), collapse = "")

  filter_rows <- paste(vapply(seq_along(params_list$filters), function(i) {
    f <- params_list$filters[[i]]
    sprintf('<tr><td style="padding:2px 8px">%s</td><td style="padding:2px 8px">%s</td></tr>', f[[1]], f[[2]])
  }, character(1)), collapse = "")

  data_rows <- paste(vapply(seq_along(data_info), function(i) {
    r <- data_info[[i]]
    sprintf('<tr><td style="padding:2px 8px">%s</td><td style="padding:2px 8px">%s</td></tr>', r[[1]], r[[2]])
  }, character(1)), collapse = "")

  tbl_style <- 'style="border-collapse:collapse;font-size:8.5px;width:100%"'
  th_style  <- 'style="padding:3px 8px;background:#f5f5f5;text-align:left;font-size:9px;font-weight:700;border-bottom:2px solid #ddd"'
  td_style  <- 'style="padding:2px 8px;vertical-align:top"'

  biotype_rows <- paste(vapply(list(
    c("ORF-annotated",           "Matches the structure of a canonical protein"),
    c("NC-variant",              "Overlaps a canonical protein with N- or C-terminal truncation/extension, or both"),
    c("ovCDS_lncRNA-ORF",        "On a lncRNA transcript, but with in-frame overlap with a CDS"),
    c("lncRNA-ORF",              "On a transcript with biotype lncRNA"),
    c("pseudogene-ORF",          "On a transcript with biotype pseudogene"),
    c("intORF",                  "Within exon boundaries of a transcript but not in the same frame"),
    c("uoORF",                   "Upstream overlapping ORF, not in-frame with reference CDS"),
    c("doORF",                   "Downstream overlapping ORF, not in-frame with reference CDS"),
    c("uORF",                    "Upstream ORF, not in-frame with reference CDS"),
    c("dORF",                    "Downstream ORF, not in-frame with reference CDS"),
    c("Processed_transcript_ORF","Within exon boundaries of a transcript with no defined CDS"),
    c("novel-ORF",               "Transcript created by StringTie"),
    c("undefined",               "No transcript assignment")
  ), function(r)
    sprintf('<tr><td %s><code style="font-size:8px">%s</code></td><td %s>%s</td></tr>',
            td_style, r[1], td_style, r[2]),
  character(1)), collapse = "")

  spec_rows <- paste(vapply(list(
    c("Safe",        "Tumor-only (logFC &gt; 2 vs every GTEx tissue; Q3 &lt; 1 TPM per tissue). No appreciable expression in adult healthy tissues."),
    c("Acceptable",  "Q3 &gt; 1 TPM only in tissues classified as low-concern by the project safety framework (e.g. Adipose, Skin, Spleen)."),
    c("Borderline",  "Q3 &gt; 1 TPM in Borderline-classified tissues (e.g. Bladder, Muscle, Thyroid, Ovary, Breast)."),
    c("Critical",    "Q3 &gt; 1 TPM in at least one Critical-classified tissue (e.g. Brain, Heart, Liver, Kidney, Lung, Nerve)."),
    c("Unavailable", "GTEx tumor-specificity data not available for this candidate. Off-tissue risk could not be assessed.")
  ), function(r)
    sprintf('<tr><td %s><b>%s</b></td><td %s>%s</td></tr>', td_style, r[1], td_style, r[2]),
  character(1)), collapse = "")

  sprintf('
<div class="rpt-slide">
  <div class="rpt-inner">
    <div style="font-size:14pt;font-weight:700;color:#1a1a1a;margin-bottom:4pt">
      Report settings &amp; inputs
    </div>
    <div style="border-top:1px solid #e0e0e0;margin-bottom:10pt"></div>
    <div style="display:flex;gap:0.4in;margin-bottom:10pt">
      <div style="flex:1">
        <div style="font-size:9px;font-weight:700;color:#444;margin-bottom:6pt">
          Scoring — Preset: <span style="color:#28646E">%s</span>
        </div>
        <table %s>
          <thead><tr><th %s>Dimension</th><th %s>Weight</th></tr></thead>
          <tbody>%s</tbody>
        </table>
      </div>
      <div style="flex:1">
        <div style="font-size:9px;font-weight:700;color:#444;margin-bottom:6pt">
          Active filters
        </div>
        <table %s>
          <thead><tr><th %s>Parameter</th><th %s>Value</th></tr></thead>
          <tbody>%s</tbody>
        </table>
        <div style="font-size:9px;font-weight:700;color:#444;margin-top:14pt;margin-bottom:6pt">
          Data inputs
        </div>
        <table %s>
          <thead><tr><th %s>Item</th><th %s>Value</th></tr></thead>
          <tbody>%s</tbody>
        </table>
      </div>
    </div>
    <div style="border-top:1px solid #e0e0e0;margin-bottom:8pt"></div>
    <div style="display:flex;gap:0.4in">
      <div style="flex:3">
        <div style="font-size:9px;font-weight:700;color:#444;margin-bottom:6pt">ORF biotype glossary</div>
        <table %s>
          <thead><tr><th %s>Biotype</th><th %s>Definition</th></tr></thead>
          <tbody>%s</tbody>
        </table>
      </div>
      <div style="flex:2">
        <div style="font-size:9px;font-weight:700;color:#444;margin-bottom:6pt">Off-tissue risk tiers</div>
        <table %s>
          <thead><tr><th %s>Tier</th><th %s>Criterion</th></tr></thead>
          <tbody>%s</tbody>
        </table>
        <div style="font-size:7.5px;color:black;margin-top:6pt">
          Low expression threshold: Q3 &lt; 1 TPM. Tissue criticality classifications reflect project-specific safety assumptions, not universal clinical standards.
        </div>
      </div>
    </div>
  </div>
  <div class="rpt-footer">
    <img src="%s" style="height:24px">
    <span>Generated: %s</span>
    <img src="%s" style="height:24px">
  </div>
</div>',
    params_list$preset %||% "Custom",
    tbl_style, th_style, th_style, weight_rows,
    tbl_style, th_style, th_style, filter_rows,
    tbl_style, th_style, th_style, data_rows,
    tbl_style, th_style, th_style, biotype_rows,
    tbl_style, th_style, th_style, spec_rows,
    logo_uri, gen_date, vh_logo_uri
  )
}

.rpt_wrap_html <- function(pages) {
  css <- '
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: Arial, Helvetica, sans-serif; background: #ccc; }
@page { size: A3 landscape; margin: 0; }
@media print {
  body { background: white; }
  .no-print { display: none !important; }
  .rpt-slide { page-break-after: always; break-after: page; }
}
.rpt-slide {
  width: 408mm; height: 284mm;
  position: relative; overflow: hidden;
  background: white;
  margin: 0 auto 12px;
  box-shadow: 0 2px 8px rgba(0,0,0,0.18);
}
.rpt-inner {
  position: absolute;
  left: 0.39in; right: 0.39in; top: 0.39in; bottom: 0.52in;
}
.rpt-footer {
  position: absolute;
  left: 0.39in; right: 0.39in; bottom: 0.07in; height: 0.38in;
  display: flex; align-items: center; justify-content: space-between;
  border-top: 1px solid #e0e0e0; padding-top: 4px;
  font-size: 12px; color: black;
}
'
  print_btn <- '
<div class="no-print" style="text-align:center;padding:16px;background:#f0f0f0;margin-bottom:0">
  <button onclick="window.print()"
          style="background:#28646E;color:white;border:none;border-radius:6px;padding:8px 20px;font-size:14px;cursor:pointer">
    &#128438; Print / Save as PDF
  </button>
  <span style="margin-left:16px;font-size:12px;color:#555">
    In the print dialog set: <b>A3 landscape</b> &middot; No margins &middot; Fit to page
  </span>
</div>'

  paste0('<!DOCTYPE html>\n<html>\n<head>\n<meta charset="UTF-8">',
         '\n<title>TITAN Report</title>\n<style>', css, '</style>\n</head>\n<body>\n',
         print_btn, '\n',
         paste(pages, collapse = "\n"),
         '\n</body>\n</html>')
}
