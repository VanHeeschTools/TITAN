## Report ggplot2 helpers — pure functions, no Shiny reactive context required.
## Depends on: ggplot2, gtex_colors, gtex_colors_subtissue, tcga_colors_normal,
##             tcga_colors_tumor, rc_color_map (all global.R).

.rpt_theme <- function(bs = 9) {
  theme_bw(base_size = bs, base_family = "sans") + theme(
    panel.grid.major.x = element_blank(),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 8,  family = "sans"),
    axis.text.y   = element_text(size = bs - 1,                     family = "sans"),
    axis.title    = element_text(size = bs, face = "bold",           family = "sans"),
    plot.margin   = margin(3, 4, 2, 3, "pt"),
    legend.position = "none",
    text          = element_text(family = "sans")
  )
}

.rpt_void <- function(xtitle, msg) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             color = "grey65", size = 3.0, hjust = 0.5, vjust = 0.5) +
    labs(x = xtitle, y = "") + theme_void(base_size = 9) +
    theme(axis.title.x = element_text(size = 12, face = "bold", color = "#555"),
          plot.margin   = margin(3, 4, 2, 3, "pt"))
  attr(p, "is_void") <- TRUE
  p
}

.rpt_plot_uri <- function(p, w, h) {
  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp))
  ggplot2::ggsave(tmp, plot = p, width = w, height = h, device = "svg",
                  bg = "white")
  # Cairo substitutes "sans" with the container's local font (e.g. "TeX Gyre Heros"),
  # which is not available in browsers and causes fallback to serif (Times).
  # Replace any quoted font-family in the SVG with a web-safe stack before encoding.
  svg_txt <- rawToChar(readBin(tmp, "raw", n = file.info(tmp)$size))
  svg_txt <- gsub('font-family: "[^"]*"',
                  'font-family: Arial, Helvetica, sans-serif',
                  svg_txt, perl = TRUE)
  paste0("data:image/svg+xml;base64,", jsonlite::base64_enc(charToRaw(svg_txt)))
}

.rpt_expr_target <- function(row, rna_mat, rna_meta, apsc, ylabel, y_ref = NULL, ylim = NULL) {
  gid <- row$gene_id_clean
  if (is.null(rna_mat) || !isTRUE(gid %in% rownames(rna_mat)))
    return(.rpt_void("Target tumor", "No RNA-seq data"))
  tpm  <- apsc(mat_row(rna_mat, gid))
  cond <- if (!is.null(rna_meta)) {
    label_col <- if ("condition"   %in% colnames(rna_meta)) rna_meta$condition
                 else if ("tissue_type" %in% colnames(rna_meta)) rna_meta$tissue_type
                 else rep("Tumor", nrow(rna_meta))
    label_col[match(colnames(rna_mat), rna_meta$sample_id)]
  } else rep("Tumor", length(tpm))
  cond[is.na(cond)] <- "Tumor"
  df <- data.frame(g = factor(cond), y = tpm)
  p <- ggplot(df, aes(x = g, y = y)) +
    geom_boxplot(fill = "#2F3D4622", color = "#2F3D46", width = 0.5,
                 outlier.shape = 1, outlier.size = 0.7) +
    geom_jitter(width = 0.14, size = 0.6, color = "#2F3D46", alpha = 0.5) +
    labs(x = NULL, y = ylabel) + .rpt_theme()
  if (!is.null(y_ref)) p <- p + geom_hline(yintercept = y_ref, linetype = "dashed", color = "#BBBBBB", linewidth = 0.35)
  if (!is.null(ylim))  p <- p + coord_cartesian(ylim = ylim)
  p
}

.rpt_expr_gtex <- function(row, gtex_mat, gtex_meta, apsc, y_ref = NULL, ylim = NULL) {
  gid <- row$gene_id_clean
  if (is.null(gtex_mat) || !isTRUE(gid %in% rownames(gtex_mat)))
    return(.rpt_void("GTEx", "GTEx matrix not loaded"))
  tpm     <- apsc(mat_row(gtex_mat, gid))
  tis_raw <- gtex_meta$tissue_type[match(colnames(gtex_mat), gtex_meta$sample_id)]
  tis_raw[is.na(tis_raw)] <- "Unknown"
  tis_lab <- gsub("_", " ", tis_raw)
  df  <- data.frame(tis_raw = tis_raw, tissue = tis_lab, y = tpm)
  med <- tapply(df$y, df$tissue, median, na.rm = TRUE)
  df$tissue <- factor(df$tissue, levels = names(sort(med, decreasing = TRUE)))
  col_map <- setNames(
    sapply(unique(tis_raw), function(t) {
      idx  <- match(t, gtex_colors_subtissue$Tissue)
      if (!is.na(idx)) return(gtex_colors_subtissue$ColorHex[idx])
      hits <- names(gtex_colors)[startsWith(t, names(gtex_colors))]
      if (length(hits)) gtex_colors[[hits[which.max(nchar(hits))]]] else "#BBBBBB"
    }),
    gsub("_", " ", unique(tis_raw))
  )
  p <- ggplot(df, aes(x = tissue, y = y, fill = tissue, color = tissue)) +
    geom_boxplot(width = 0.7, alpha = 0.18, outlier.shape = 1, outlier.size = 0.3) +
    scale_fill_manual(values = col_map)  +
    scale_color_manual(values = col_map) +
    labs(x = "GTEx", y = "") + .rpt_theme()
  if (!is.null(y_ref)) p <- p + geom_hline(yintercept = y_ref, linetype = "dashed", color = "#BBBBBB", linewidth = 0.35)
  if (!is.null(ylim))  p <- p + coord_cartesian(ylim = ylim)
  p
}

.rpt_expr_tcga <- function(row, tcga_mat, tcga_meta, apsc, y_ref = NULL, ylim = NULL) {
  gid <- row$gene_id_clean
  if (is.null(tcga_mat) || !isTRUE(gid %in% rownames(tcga_mat)))
    return(.rpt_void("TCGA", "TCGA matrix not loaded"))
  tpm     <- apsc(mat_row(tcga_mat, gid))
  grp_raw <- tcga_meta$group[match(colnames(tcga_mat), tcga_meta$sample_id)]
  grp_raw[is.na(grp_raw)] <- "Unknown"
  df       <- data.frame(grp = grp_raw, y = tpm)
  gmeds    <- tapply(df$y, df$grp, median, na.rm = TRUE)
  ct_codes <- sub(" .*", "", unique(grp_raw))
  pk       <- vapply(unique(ct_codes), function(ct) {
    m <- gmeds[sub(" .*", "", names(gmeds)) == ct]
    if (!length(m)) -Inf else max(m, na.rm = TRUE)
  }, numeric(1))
  ct_sorted <- unique(ct_codes)[order(pk, decreasing = TRUE)]
  go <- unlist(lapply(ct_sorted, function(ct) {
    gs <- unique(grp_raw)[sub(" .*", "", unique(grp_raw)) == ct]
    gs[order(match(sub(".* ", "", gs), c("Tumor", "Normal")))]
  }), use.names = FALSE)
  df$grp <- factor(df$grp, levels = go)
  col_map <- sapply(go, function(lvl) {
    ct <- sub(" .*", "", lvl)
    if (grepl("Normal", lvl, ignore.case = TRUE)) tcga_colors_normal[[ct]] %||% "#87C8D4"
    else tcga_colors_tumor[[ct]] %||% "#3B95A5"
  })
  label_map <- setNames(gsub(" Tumor$", " T", gsub(" Normal$", " N", go)), go)
  p <- ggplot(df, aes(x = grp, y = y, fill = grp, color = grp)) +
    geom_boxplot(width = 0.7, alpha = 0.35, outlier.shape = 1, outlier.size = 0.3) +
    scale_fill_manual(values = col_map)  +
    scale_color_manual(values = col_map) +
    scale_x_discrete(labels = label_map) +
    labs(x = "TCGA", y = "") + .rpt_theme()
  if (!is.null(y_ref)) p <- p + geom_hline(yintercept = y_ref, linetype = "dashed", color = "#BBBBBB", linewidth = 0.35)
  if (!is.null(ylim))  p <- p + coord_cartesian(ylim = ylim)
  p
}

.rpt_transl_target <- function(row, ribo_m, ribo_sm, apsc, ylabel, y_ref = NULL, ylim = NULL) {
  oid <- row$orf_id
  if (is.null(ribo_m) || !isTRUE(oid %in% rownames(ribo_m)))
    return(.rpt_void("Target tumor", "No ribo-seq data\nfor this ORF"))
  ppm   <- apsc(mat_row(ribo_m, oid))
  sids  <- colnames(ribo_m)
  cond  <- if (!is.null(ribo_sm) && "condition" %in% colnames(ribo_sm))
    ribo_sm$condition[match(sids, ribo_sm$sample_id)]
  else rep("Tumor", length(ppm))
  cond[is.na(cond)] <- "Tumor"
  df <- data.frame(g = factor(cond), y = ppm)
  p <- ggplot(df, aes(x = g, y = y)) +
    geom_boxplot(fill = "#2F3D4622", color = "#2F3D46", width = 0.5,
                 outlier.shape = 1, outlier.size = 0.7) +
    geom_jitter(width = 0.14, size = 0.6, color = "#2F3D46", alpha = 0.5) +
    labs(x = NULL, y = ylabel) + .rpt_theme()
  if (!is.null(y_ref)) p <- p + geom_hline(yintercept = y_ref, linetype = "dashed", color = "#BBBBBB", linewidth = 0.35)
  if (!is.null(ylim))  p <- p + coord_cartesian(ylim = ylim)
  p
}

.rpt_transl_rc <- function(row, rc_mat, rc_meta, grp_label, apsc, label_fn = identity, y_ref = NULL, ylim = NULL) {
  oid        <- row$orf_id
  tgt_grp    <- if (grepl("Primary", grp_label, ignore.case = TRUE)) "Primary" else "Cell-line"
  if (is.null(rc_mat) || !isTRUE(oid %in% rownames(rc_mat)))
    return(.rpt_void(grp_label, "Ribocrypt matrix not loaded"))
  rc_raw  <- mat_row(rc_mat, oid)
  sids    <- colnames(rc_mat)
  grp     <- rc_meta$group[match(sids, rc_meta$sample_id)]
  grp[is.na(grp)] <- "Unknown"
  mask <- grp == tgt_grp
  if (!any(mask)) return(.rpt_void(grp_label, paste("No", grp_label, "data")))
  s_ids <- sids[mask]; s_y <- apsc(rc_raw[mask])
  ord   <- order(s_y, decreasing = TRUE)
  s_ids <- s_ids[ord]; s_y <- s_y[ord]
  s_lab <- factor(label_fn(s_ids), levels = label_fn(s_ids))
  s_col <- ifelse(s_ids %in% names(rc_color_map), rc_color_map[s_ids], "#AAAAAA")
  n     <- length(s_lab)
  df_st <- data.frame(x = rep(s_lab, each = 3),
                      y = as.vector(rbind(rep(0, n), s_y, NA_real_)),
                      id = rep(seq_len(n), each = 3))
  df_pt <- data.frame(x = s_lab, y = s_y, fill = s_col)
  p <- ggplot() +
    geom_line(data = df_st, aes(x = x, y = y, group = id),
              color = "#CCCCCC", linewidth = 0.45, na.rm = TRUE) +
    geom_point(data = df_pt, aes(x = x, y = y, fill = I(fill), color = I(fill)),
               shape = 21, size = 2.2, stroke = 0) +
    scale_x_discrete(limits = levels(s_lab)) +
    labs(x = grp_label, y = "") + .rpt_theme() +
    theme(axis.text.x = element_text(size = 8))
  if (!is.null(y_ref)) p <- p + geom_hline(yintercept = y_ref, linetype = "dashed", color = "#BBBBBB", linewidth = 0.35)
  if (!is.null(ylim))  p <- p + coord_cartesian(ylim = ylim)
  p
}

.rpt_ylim <- function(vals_list) {
  v <- unlist(vals_list); v <- v[is.finite(v)]
  if (!length(v)) return(NULL)
  rng  <- range(v, na.rm = TRUE)
  span <- max(diff(rng), 1e-6)
  c(rng[1] - span * 0.04, rng[2] + span * 0.04)
}
