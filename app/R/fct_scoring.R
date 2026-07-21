## Scoring functions and HTML badge renderers.
## Depends on: BIOTYPE_COLORS, WEIGHT_META (global.R).

preset_weights <- function(preset_name) {
  field <- switch(preset_name,
    "Pan-cancer"      = "pancancer",
    "Cancer-specific" = "specific"
  )
  setNames(sapply(WEIGHT_META, `[[`, field), sapply(WEIGHT_META, `[[`, "id"))
}

dim_col <- function(wid) sub("^w_", "dim_", wid)

score_candidates <- function(df, w) {
  if (nrow(df) == 0) return(df)

  safe_max <- function(x) { v <- max(x, na.rm = TRUE); if (is.finite(v) && v > 0) v else 1 }
  max_gtex      <- safe_max(df$GTEX_max_median_TPM)
  max_ribo_prim <- safe_max(df$ribocrypt_primary_max_PPM)
  max_ribo_cell <- safe_max(df$`ribocrypt_cell-line_max_PPM`)

  wv <- setNames(as.numeric(unlist(w)), names(w))

  # Signals are always in [0, 1]; contribution = signal × weight.
  # max_possible = sum of positive weights (all signals at 1 for those dims)
  # min_possible = sum of negative weights (all signals at 1 for those dims)
  max_possible <- sum(pmax(0, wv))
  min_possible <- sum(pmin(0, wv))
  score_range  <- max_possible - min_possible

  df %>%
    mutate(
      dim_pct_samples  = pmin(replace_na(target_expression_pct_samples, 0) / 100, 1) * wv["w_pct_samples"],
      dim_pct_transl   = pmin(replace_na(target_translation_pct_samples, 0) / 100, 1) * wv["w_pct_transl"],
      dim_tumor_spec   = case_when(
                           GTEX_tumor_only %in% TRUE     ~ 1,
                           GTEX_tumor_enriched %in% TRUE ~ 0.5,
                           TRUE                          ~ 0
                         ) * wv["w_tumor_spec"],
      dim_gtex_penalty = pmin(replace_na(GTEX_max_median_TPM, 0) / max_gtex, 1) * wv["w_gtex_penalty"],
      dim_tcga_cov     = pmin(replace_na(TCGA_tumor_pct_samples, 0) / 100, 1) * wv["w_tcga_cov"],
      dim_peri_penalty = pmin(replace_na(TCGA_normal_pct_samples, 0) / 100, 1) * wv["w_peri_penalty"],
      dim_ribo_primary = pmin(replace_na(ribocrypt_primary_median_PPM, 0) / (max_ribo_prim * 0.5), 1) * wv["w_ribo_primary"],
      dim_ribo_cell    = pmin(replace_na(`ribocrypt_cell-line_median_PPM`, 0) / (max_ribo_cell * 0.5), 1) * wv["w_ribo_cell"],
      raw_score        = dim_pct_samples + dim_pct_transl + dim_tumor_spec + dim_gtex_penalty +
                         dim_tcga_cov + dim_peri_penalty + dim_ribo_primary + dim_ribo_cell,
      priority_score   = if (score_range == 0) 50 else
                           pmax(0, pmin(100, (raw_score - min_possible) / score_range * 100))
    )
}

score_bar_html <- function(s) {
  col <- if (s > 60) "#00A555" else if (s > 35) "#D4850A" else "#C0392B"
  sprintf(
    '<div class="titan-score-wrap"><div class="titan-score-track"><div class="titan-score-fill" style="width:%.0f%%;background:%s"></div></div><span class="titan-score-val">%.1f</span></div>',
    s, col, s
  )
}

pct_bar_html <- function(pct, color) {
  v <- if (is.na(pct) || !is.finite(pct)) 0 else max(0, min(100, pct))
  sprintf(
    '<div class="titan-score-wrap"><div class="titan-score-track"><div class="titan-score-fill" style="width:%.0f%%;background:%s"></div></div><span class="titan-score-val">%.0f%%</span></div>',
    v, color, v
  )
}

spec_badge_html <- function(tumor_only, tumor_enriched) {
  if (isTRUE(as.logical(tumor_only)))          '<span class="titan-badge titan-badge-specific">Tumor-only</span>'
  else if (isTRUE(as.logical(tumor_enriched))) '<span class="titan-badge titan-badge-enriched">Tumor-enriched</span>'
  else                                         '<span class="titan-badge titan-badge-other">Non-specific</span>'
}

biotype_badge_html <- function(biotype) {
  col <- unname(BIOTYPE_COLORS[biotype])
  if (is.na(col)) col <- "#95A5A6"
  sprintf('<span class="titan-badge titan-badge-biotype" style="background:%s22;color:#111111;border:1px solid %s55;">%s</span>',
          col, col, biotype)
}
