## TITAN - Tumor Immunopeptidomics Target ANnotation
## Shiny app for prioritising ncORF-derived peptide candidates
##
## Usage:
##   shiny::runApp("/hpc/pmc_oatv/projects/tools_dev/titan/app")

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(plotly)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(shinyWidgets)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

options(shiny.maxRequestSize = 1000 * 1024^2)  # 1GB upload limit

# ─────────────────────────────────────────────────────────────────────────────
# DATA  (demo loaded at startup for initial UI state; reactive in server)
# ─────────────────────────────────────────────────────────────────────────────

app_data         <- readRDS("data/titan_orf_table.rds")
orf_table        <- app_data$orf_table
ribo_ppm_samples <- app_data$ribo_ppm_samples
ribo_sample_meta <- app_data$ribo_sample_meta
rna_tpm_mat      <- app_data$rna_tpm_mat

if (!"gene_id_clean" %in% colnames(orf_table))
  orf_table$gene_id_clean <- sub("\\..*", "", orf_table$gene_id)

n_ribo_samples <- ncol(ribo_ppm_samples)
n_rna_samples  <- if (!is.null(rna_tpm_mat)) ncol(rna_tpm_mat) else nrow(app_data$rna_sample_meta)

biotypes <- sort(unique(orf_table$orf_biotype_single))

# ─────────────────────────────────────────────────────────────────────────────
# GENCODE ORF TABLE  (optional; enables cross-matching against TransCode Phase 2)
# ─────────────────────────────────────────────────────────────────────────────
gencode_orf_tbl <- local({
  f <- "data/gencode_orfs_phase2.csv"
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(
    data.table::fread(f, data.table = FALSE, showProgress = FALSE),
    error = function(e) { message("gencode_orfs_phase2.csv not loaded: ", e$message); NULL }
  )
  if (is.null(df)) return(NULL)
  # Map TransCode orf_type to app biotype labels (PT is the only rename)
  df$orf_biotype_single <- ifelse(df$orf_type == "PT", "Processed_transcript_ORF", df$orf_type)
  df$protein_seq    <- df$sequence_aa
  df$protein_length <- nchar(df$sequence_aa)
  df$chr            <- df$chrm
  df$orf_start      <- as.integer(df[["starts (0-based)"]])
  df$orf_end        <- as.integer(df[["ends (0-based)"]])
  df$start_codon    <- df$initiation_codon
  df$orf_id         <- df$releasev45_id
  df$gene_id_clean  <- sub("\\..*", "", df$gene_id)
  df[, c("orf_id", "gene_id", "gene_name", "gene_biotype", "orf_biotype_single",
         "protein_seq", "protein_length", "chr", "orf_start", "orf_end",
         "strand", "start_codon", "gene_id_clean")]
})

# ─────────────────────────────────────────────────────────────────────────────
# COLOUR PALETTES  (GTEx tissue groups / TCGA studies)
# ─────────────────────────────────────────────────────────────────────────────

.darken_hex <- function(hex, amount = 0.35) {
  # Darken hex color(s) by reducing HLS lightness (no colorspace dependency).
  sapply(hex, function(h) {
    v  <- col2rgb(h) / 255
    r  <- v[1]; g <- v[2]; b <- v[3]
    mx <- max(r, g, b); mn <- min(r, g, b); d <- mx - mn
    l  <- (mx + mn) / 2
    s  <- if (d == 0) 0 else d / (1 - abs(2 * l - 1))
    hh <- if (d == 0)    0
          else if (mx == r) 60 * (((g - b) / d) %% 6)
          else if (mx == g) 60 * ((b - r) / d + 2)
          else              60 * ((r - g) / d + 4)
    l2 <- max(0, l * (1 - amount))
    cv <- (1 - abs(2 * l2 - 1)) * s
    xv <- cv * (1 - abs((hh / 60) %% 2 - 1))
    m  <- l2 - cv / 2
    rn <- if      (hh < 60)  c(cv, xv,  0)
          else if (hh < 120) c(xv, cv,  0)
          else if (hh < 180) c(0,  cv,  xv)
          else if (hh < 240) c(0,  xv,  cv)
          else if (hh < 300) c(xv, 0,   cv)
          else               c(cv, 0,   xv)
    out <- pmin(1, pmax(0, rn + m))
    rgb(out[1], out[2], out[3])
  }, USE.NAMES = FALSE)
}

gtex_colors <- c(
  Adipose         = "#D1B9A5",
  Adrenal_Gland   = "#BB80B1",
  Artery          = "#FFB08E",
  Bladder         = "#B6A1D9",
  Brain           = "#90AFD4",
  Breast          = "#F0C2B5",
  Cervix          = "#E5B4AE",
  Colon           = "#FACC78",
  Esophagus       = "#FE9C62",
  Fallopian_Tube  = "#DBA6A8",
  Heart           = "#EC8675",
  Kidney          = "#BCB9EB",
  Liver           = "#ECF3A8",
  Lung            = "#DCD4AD",
  Muscle          = "#C4D4AC",
  Nerve           = "#BEBBB8",
  Ovary           = "#D099A1",
  Pancreas        = "#CFEBB7",
  Pituitary       = "#CD9EC7",
  Prostate        = "#B390CD",
  Salivary        = "#C794C0",
  Skin            = "#E1CEA9",
  Small_Intestine = "#FBBA70",
  Spleen          = "#F1E496",
  Stomach         = "#FDAC6A",
  Thyroid         = "#DFBCDD",
  Uterus          = "#C68B9B",
  Vagina          = "#BB7D94",
  Whole_Blood     = "#DD6661"
)

gtex_colors_subtissue <- data.frame(
  Tissue = c(
    "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c_1", "Brain_Substantia_nigra",
    "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis",
    "Stomach", "Small_Intestine_Terminal_Ileum",
    "Colon_Sigmoid", "Colon_Transverse",
    "Spleen", "Liver", "Pancreas", "Lung",
    "Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
    "Heart_Atrial_Appendage", "Heart_Left_Ventricle",
    "Whole_Blood",
    "Kidney_Cortex", "Kidney_Medulla",
    "Bladder", "Prostate", "Nerve_Tibial", "Muscle_Skeletal", "Breast_Mammary_Tissue",
    "Cervix_Ectocervix", "Cervix_Endocervix",
    "Fallopian_Tube", "Ovary", "Uterus", "Vagina", "Adrenal_Gland",
    "Minor_Salivary_Gland", "Pituitary", "Thyroid",
    "Adipose_Subcutaneous", "Adipose_Visceral_Omentum",
    "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg"
  ),
  Group = c(
    rep("Brain", 13),
    rep("Esophagus", 3),
    "Stomach", "Small_Intestine",
    "Colon", "Colon",
    "Spleen", "Liver", "Pancreas", "Lung",
    rep("Artery", 3),
    "Heart", "Heart",
    "Whole_Blood",
    "Kidney", "Kidney",
    "Bladder", "Prostate", "Nerve", "Muscle", "Breast",
    "Cervix", "Cervix",
    "Fallopian_Tube", "Ovary", "Uterus", "Vagina", "Adrenal_Gland",
    "Salivary", "Pituitary", "Thyroid",
    "Adipose", "Adipose",
    "Skin", "Skin"
  ),
  ColorHex = c(
    "#90AFD4", "#88A6CC", "#819DC5", "#7994BD", "#7994BD", "#718BB5", "#6A82AE",
    "#627AA6", "#5B719F", "#536897", "#4B5F8F", "#445688", "#3C4D80",
    "#FE9C62", "#FE945F", "#FF8C5B",
    "#FDAC6A", "#FBBA70",
    "#FACC78", "#FBC374",
    "#F1E496", "#ECF3A8", "#CFEBB7", "#DCD4AD",
    "#FFB08E", "#FAA588", "#F59B81",
    "#EC8675", "#E77B6E",
    "#DD6661",
    "#BCB9EB", "#BAB1E5",
    "#B6A1D9", "#B390CD", "#BEBBB8", "#C4D4AC", "#F0C2B5",
    "#EBBBB2", "#E5B4AE",
    "#DBA6A8", "#D099A1", "#C68B9B", "#BB7D94", "#BB80B1",
    "#C794C0", "#CD9EC7", "#DFBCDD",
    "#CDB599", "#D1B9A5",
    "#E1CEA9", "#ECCBA8"
  ),
  stringsAsFactors = FALSE
)

# TCGA study base ("normal") colors: GTEx group rep or sub-tissue alternate for collisions.
# Collisions with only 1 sub-tissue shade share a color: ACC=PCPG, CHOL=LIHC, LUAD=LUSC, UCEC=UCS.
# KIRC=KICH (Kidney_Cortex); KIRP uses Kidney_Medulla. Brain resolved via gradient extremes.
tcga_colors_normal <- c(
  ACC  = "#BB80B1",   # Adrenal_Gland group rep  (= PCPG, 1 shade)
  BLCA = "#B6A1D9",
  BRCA = "#F0C2B5",
  CESC = "#E5B4AE",
  CHOL = "#ECF3A8",   # Liver group rep  (= LIHC, 1 shade)
  COAD = "#FACC78",   # Colon_Sigmoid
  DLBC = "#B8AADC",   # custom -- lymphoid, no GTEx match
  ESCA = "#FE9C62",
  GBM  = "#90AFD4",   # Brain_Amygdala (lightest end of gradient)
  HNSC = "#E8C0A0",   # custom -- head/neck, no GTEx match
  KICH = "#BCB9EB",   # Kidney_Cortex  (= KIRC, only 2 shades for 3 studies)
  KIRC = "#BCB9EB",   # Kidney_Cortex
  KIRP = "#BAB1E5",   # Kidney_Medulla
  LAML = "#DD6661",
  LGG  = "#3C4D80",   # Brain_Substantia_nigra (darkest end of gradient)
  LIHC = "#ECF3A8",   # Liver group rep  (= CHOL, 1 shade)
  LUAD = "#DCD4AD",   # Lung group rep  (= LUSC, 1 shade)
  LUSC = "#DCD4AD",
  MESO = "#A8C4B0",   # custom -- pleura, no GTEx match
  OV   = "#D099A1",
  PAAD = "#CFEBB7",
  PCPG = "#BB80B1",   # Adrenal_Gland group rep  (= ACC, 1 shade)
  PRAD = "#B390CD",
  READ = "#FBC374",   # Colon_Transverse
  SARC = "#C4D4AC",
  SKCM = "#E1CEA9",
  STAD = "#FDAC6A",
  TGCT = "#C4D4E4",   # custom -- testicular, no GTEx match
  THCA = "#DFBCDD",
  THYM = "#D4C0D4",   # custom -- thymic, no GTEx match
  UCEC = "#C68B9B",   # Uterus group rep  (= UCS, 1 shade)
  UCS  = "#C68B9B"
)

tcga_colors_tumor <- setNames(
  .darken_hex(tcga_colors_normal, amount = 0.35),
  names(tcga_colors_normal)
)

# Ribocrypt sample → GTEx-palette color (GroupColorHex or custom where noted)
rc_color_map <- c(
  # Primary tissues
  hepatocyte_liver              = "#ECF3A8",  # Liver
  HSPC_blood                    = "#DD6661",  # Whole_Blood
  huvec_Umbilical               = "#FFB08E",  # Artery (vascular endothelium)
  Myoblast_muscle               = "#C4D4AC",  # Muscle
  primary_brain                 = "#90AFD4",  # Brain
  primary_corneal_eye           = "#AAAAAA",  # no GTEx match
  primary_fibroblast_connective = "#E1CEA9",  # Skin (dermal fibroblasts)
  primary_heart                 = "#EC8675",  # Heart
  primary_kidney                = "#BCB9EB",  # Kidney
  primary_liver                 = "#ECF3A8",  # Liver
  primary_muscle                = "#C4D4AC",  # Muscle
  primary_skin                  = "#E1CEA9",  # Skin
  # Cell lines
  A549_lung                     = "#DCD4AD",  # Lung
  BJ_skin                       = "#E1CEA9",  # Skin
  Calu3_lung                    = "#DCD4AD",  # Lung
  CTC_blood                     = "#DD6661",  # Whole_Blood
  embryonic_brain               = "#90AFD4",  # Brain
  ENDOC_pancreas                = "#CFEBB7",  # Pancreas
  fibroblast_heart              = "#EC8675",  # Heart
  fibroblast_skin               = "#E1CEA9",  # Skin
  glioblastoma_brain            = "#90AFD4",  # Brain
  H1_embryonic                  = "#AAAAAA",  # no tissue equivalent (hESC)
  H1299_lung                    = "#DCD4AD",  # Lung
  H9_embryonic                  = "#AAAAAA",  # no tissue equivalent (hESC)
  HAP1_bone                     = "#DD6661",  # Whole_Blood (bone-marrow CML origin)
  HBEC_lung                     = "#DCD4AD",  # Lung
  HCT116_colon                  = "#FACC78",  # Colon
  HDF_skin                      = "#E1CEA9",  # Skin
  HEK293_kidney                 = "#BCB9EB",  # Kidney
  HeLa_ovary                    = "#D099A1",  # Ovary
  HepG2_liver                   = "#ECF3A8",  # Liver
  hesc_embryonic                = "#AAAAAA",  # no tissue equivalent (hESC)
  HFF_skin                      = "#E1CEA9",  # Skin
  hTERT_RPE1_eye                = "#AAAAAA",  # no GTEx match (RPE)
  Huh7_liver                    = "#ECF3A8",  # Liver
  IMR90_lung                    = "#DCD4AD",  # Lung
  iPSC_brain                    = "#90AFD4",  # Brain
  iPSC_none                     = "#AAAAAA",  # no tissue specified
  K562_bone                     = "#DD6661",  # Whole_Blood (bone-marrow CML origin)
  LCL_blood                     = "#DD6661",  # Whole_Blood
  LN299_brain                   = "#90AFD4",  # Brain
  LN308_brain                   = "#90AFD4",  # Brain
  MCF10A_breast                 = "#F0C2B5",  # Breast
  MCF7_breast                   = "#F0C2B5",  # Breast
  MDA_breast                    = "#F0C2B5",  # Breast
  MM1_blood                     = "#DD6661",  # Whole_Blood
  MOLM13_blood                  = "#DD6661",  # Whole_Blood
  MRC5_lung                     = "#DCD4AD",  # Lung
  Neuroblast_brain              = "#90AFD4",  # Brain
  NPC_throat                    = "#FE9C62",  # Esophagus (upper aerodigestive)
  PANC1_pancreas                = "#CFEBB7",  # Pancreas
  PC3_prostate                  = "#B390CD",  # Prostate
  RD_muscle                     = "#E4E6AB",  # custom (sarcoma)
  SHSY5Y_brain                  = "#90AFD4",  # Brain
  SUM159_breast                 = "#F0C2B5",  # Breast
  THP1_blood                    = "#DD6661",  # Whole_Blood
  U2392_blood                   = "#DD6661",  # Whole_Blood
  U2OS_bone                     = "#C4D4AC",  # Muscle (mesenchymal/bone)
  UOK262_kidney                 = "#BCB9EB",  # Kidney
  VSMC_muscle                   = "#F59B81",  # custom (Artery_Tibial, vascular SM)
  Wilmstumor_kidney             = "#BCB9EB"   # Kidney
)

# ─────────────────────────────────────────────────────────────────────────────
# SCORING SYSTEM
# ─────────────────────────────────────────────────────────────────────────────

# Previous scoring system (v1) - kept for reference:
# WEIGHT_META_V1 <- list(
#   list(id="w_pct_samples",  label="% samples (expr.)",       dir="+", balanced=18, pancancer=20, specific=15),
#   list(id="w_expr_level",   label="Median expression level",  dir="+", balanced=12, pancancer=14, specific=10),
#   list(id="w_tumor_spec",   label="Tumor specificity (GTEx)",dir="+", balanced=18, pancancer=10, specific=25),
#   list(id="w_gtex_penalty", label="GTEx expression penalty",  dir="−", balanced=14, pancancer=8,  specific=20),
#   list(id="w_tcga_cov",     label="TCGA tumor coverage",     dir="+", balanced=14, pancancer=18, specific=12),
#   list(id="w_peri_penalty", label="Peritumor penalty",       dir="−", balanced=12, pancancer=10, specific=15),
#   list(id="w_ribo_primary", label="Ribocrypt primary tissue", dir="−", balanced=14, pancancer=12, specific=18),
#   list(id="w_ribo_cell",    label="Ribocrypt cell-line",      dir="+", balanced=7,  pancancer=8,  specific=5)
# )
# dim_pct_samples:  pmin(expr_pct/100,1)*w          [0,w]
# dim_expr_level:   pmin(TPM/(max*0.5),1)*w          [0,w]
# dim_tumor_spec:   case(only→1, enriched→0.5, else→0)*w [0,w]
# dim_gtex_penalty: -pmin(GTEx/max,1)*w              [-w,0]
# priority_score:   raw/total_w * 100                [0,100]

# Each dimension exposes a signal in [0, 1].
# Weight = +1 → high signal adds to score; −1 → high signal penalises; 0 → ignored.
# hint describes what "high signal" means for each dimension.
WEIGHT_META <- list(
  list(id="w_pct_samples",  label="% expressed (RNA-seq)",    hint="many tumor samples express the ORF",             radar="Expr. %",    group="Tumor coverage",           balanced=0.6,  pancancer=0.9, specific=0.4),
  list(id="w_pct_transl",   label="% translated (Ribo-seq)",  hint="many tumor samples translate the ORF",           radar="Transl. %",  group="Tumor coverage",           balanced=0.6,  pancancer=0.9, specific=0.4),
  list(id="w_tumor_spec",   label="Tumor specificity (GTEx)", hint="GTEx: tumor-only=1, enriched=0.5, non-specific=0",radar="Specificity",group="Specificity",              balanced=0.5,  pancancer=0.2, specific=0.9),
  list(id="w_gtex_penalty", label="GTEx expression level",    hint="expressed in normal tissues (GTEx)",             radar="GTEx",       group="Specificity",              balanced=-0.4, pancancer=-0.1,specific=-0.8),
  list(id="w_tcga_cov",     label="TCGA tumor coverage",      hint="expressed in many TCGA tumor samples",           radar="TCGA T",     group="TCGA validation",          balanced=0.4,  pancancer=0.7, specific=0.3),
  list(id="w_peri_penalty", label="TCGA normal expression",   hint="expressed in peritumoral / normal TCGA samples", radar="TCGA N",     group="TCGA validation",          balanced=-0.3, pancancer=-0.1,specific=-0.6),
  list(id="w_ribo_primary", label="RC primary tissue",        hint="translated in normal primary tissues (Ribocrypt)",radar="RC primary", group="Normal tissue (Ribocrypt)",balanced=-0.4, pancancer=-0.1,specific=-0.7),
  list(id="w_ribo_cell",    label="RC cell-line",             hint="translated in normal cell lines (Ribocrypt)",    radar="RC CL",      group="Normal tissue (Ribocrypt)",balanced=0.0,  pancancer=0.1, specific=0.0)
)

PRESETS <- list(
  "Cancer-specific" = list(label = "Strict tumor specificity, penalises normal tissue", color = "#FFBEFF"),
  "Pan-cancer"      = list(label = "Broad coverage, tolerates enriched targets",        color = "#28646E")
)

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

# ─────────────────────────────────────────────────────────────────────────────
# FORMATTING HELPERS
# ─────────────────────────────────────────────────────────────────────────────

fmt1     <- function(x) if (is.na(x) || length(x) == 0) "—" else sprintf("%.1f%%", x)
fmt2     <- function(x) if (is.na(x) || length(x) == 0) "—" else sprintf("%.3f", x)
bool_fmt <- function(x) if (is.na(x) || length(x) == 0) "—" else if (isTRUE(x)) "Yes ✓" else "No"

BIOTYPE_COLORS <- c(
  "ORF-annotated"            = "#ED6A5A",
  "Processed_transcript_ORF" = "#F8A598",
  "NC-variant"               = "#ADF7B6",
  "uORF"                     = "#ff7f00",
  "uoORF"                    = "#fdbf6f",
  "dORF"                     = "#F3D264",
  "doORF"                    = "#F8E39C",
  "intORF"                   = "#CDD1E0",
  "lncRNA-ORF"               = "#F6B9FF",
  "pseudogene-ORF"           = "#25CED1"
)

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

match_peptides <- function(peptides, orf_tbl) {
  peptides  <- unique(trimws(peptides))
  peptides  <- peptides[nchar(peptides) >= 8]
  if (length(peptides) == 0) return(NULL)
  canonical <- c("ORF-annotated", "NC-variant")
  results <- lapply(peptides, function(pep) {
    hits <- which(str_detect(orf_tbl$protein_seq, fixed(pep)))
    if (length(hits) == 0) return(NULL)
    matched <- orf_tbl[hits, , drop = FALSE] %>% mutate(matched_peptide = pep)
    # Peptides matching a canonical biotype are not evidence for ncORFs
    if (any(matched$orf_biotype_single %in% canonical))
      matched <- matched[matched$orf_biotype_single %in% canonical, , drop = FALSE]
    matched
  })
  bind_rows(results)
}

pill_badge <- function(text, color = "primary") {
  tags$span(class = paste0("badge rounded-pill bg-", color, " me-1"), text)
}

orf_id_labels <- function(tbl) {
  paste0(tbl$gene_name, "_", tbl$orf_biotype_single, "_",
         tbl$protein_length, "aa_",
         tbl$chr, ":", tbl$orf_start, "-", tbl$orf_end, "_",
         tbl$start_codon)
}

# Escape a string for use inside an HTML attribute value (data-bs-content etc.)
html_attr_escape <- function(s) {
  s <- gsub("&",  "&amp;",  s, fixed = TRUE)
  s <- gsub("<",  "&lt;",   s, fixed = TRUE)
  s <- gsub(">",  "&gt;",   s, fixed = TRUE)
  s <- gsub('"',  "&quot;", s, fixed = TRUE)
  s <- gsub("'",  "&#39;",  s, fixed = TRUE)
  s
}

# Build an HTML table for the Bootstrap popover from a data.frame of MS rows
build_pep_popover <- function(rows_df) {
  if (is.null(rows_df) || nrow(rows_df) == 0L || ncol(rows_df) == 0L) return("")
  make_table <- function(row) {
    cells <- paste(vapply(colnames(row), function(col) {
      sprintf('<tr><td class="pep-tt-key">%s</td><td class="pep-tt-val">%s</td></tr>',
              col, row[[col]])
    }, character(1)), collapse = "")
    sprintf('<table class="pep-tt-table">%s</table>', cells)
  }
  parts <- vapply(seq_len(nrow(rows_df)), function(i) make_table(rows_df[i, , drop = FALSE]),
                  character(1))
  paste(parts, collapse = "<hr class='my-1'>")
}

# Render protein sequence as HTML with per-peptide colour highlights,
# alignment rows, and Bootstrap popover tooltips showing MS data on hover.
# pep_info: named list  peptide → data.frame of MS rows (for popover)
render_protein_seq_html <- function(seq, pep_list, pep_info = list()) {
  PEP_COLS <- c("#28646E", "#D4850A", "#8E44AD", "#C0392B", "#0097A7")
  pep_list <- unique(pep_list[!is.na(pep_list) & nzchar(pep_list)])
  n_chars  <- nchar(seq)
  seq_v    <- strsplit(seq, "")[[1]]

  # Mark which peptide (1-indexed) first covers each position
  coverage   <- integer(n_chars)
  pep_starts <- vector("list", length(pep_list))
  for (pi in seq_along(pep_list)) {
    m <- gregexpr(pep_list[[pi]], seq, fixed = TRUE)[[1]]
    if (m[1L] > 0L) {
      pep_starts[[pi]] <- m
      plen <- nchar(pep_list[[pi]])
      for (s in m) for (j in s:min(s + plen - 1L, n_chars)) if (!coverage[j]) coverage[j] <- pi
    }
  }

  BLOCK  <- 60L
  INDENT <- "     "   # 4-digit line number + 1 space

  blocks <- vapply(seq_len(ceiling(n_chars / BLOCK)), function(b) {
    i0 <- (b - 1L) * BLOCK + 1L
    i1 <- min(b * BLOCK, n_chars)

    # Sequence row - use RLE to group runs of same peptide into one span
    cov_range <- coverage[i0:i1]
    seq_range <- seq_v[i0:i1]
    r_len  <- rle(cov_range)$lengths
    r_val  <- rle(cov_range)$values
    r_end  <- cumsum(r_len)
    r_start <- c(1L, r_end[-length(r_end)] + 1L)

    seq_parts <- mapply(function(rs, re, cv) {
      chars <- paste(seq_range[rs:re], collapse = "")
      if (cv == 0L) return(chars)
      col <- PEP_COLS[(cv - 1L) %% length(PEP_COLS) + 1L]
      pep <- pep_list[[cv]]
      info      <- pep_info[[pep]]
      n_records <- if (!is.null(info) && nrow(info) > 0L) nrow(info) else 0L
      tt_title   <- html_attr_escape(sprintf("MS data (%d record%s)", n_records, if (n_records == 1L) "" else "s"))
      tt_content <- html_attr_escape(build_pep_popover(info))
      sprintf(
        '<span class="pep-hit" style="background:%s33;color:%s;font-weight:bold;" data-bs-toggle="popover" data-bs-html="true" data-bs-placement="top" data-bs-trigger="hover focus" data-bs-title="%s" data-bs-content="%s">%s</span>',
        col, col, tt_title, tt_content, chars
      )
    }, r_start, r_end, r_val, SIMPLIFY = TRUE)
    seq_line <- sprintf('<span class="seq-pos">%4d</span> %s', i0, paste(seq_parts, collapse = ""))

    # One alignment row per peptide overlapping this block
    aln_rows <- vapply(seq_along(pep_list), function(pi) {
      starts <- pep_starts[[pi]]
      if (is.null(starts) || starts[1L] < 0L) return("")
      pv   <- strsplit(pep_list[[pi]], "")[[1]]
      col  <- PEP_COLS[(pi - 1L) %% length(PEP_COLS) + 1L]
      aln  <- rep(" ", i1 - i0 + 1L)
      for (s in starts) {
        if (s + length(pv) - 1L < i0 || s > i1) next
        for (j in seq_along(pv)) {
          ap <- s + j - 1L
          if (ap >= i0 && ap <= i1) aln[ap - i0 + 1L] <- pv[j]
        }
      }
      if (all(aln == " ")) return("")
      spans <- vapply(aln, function(ch) {
        if (ch == " ") return(ch)
        sprintf('<span style="color:%s;font-weight:bold;">%s</span>', col, ch)
      }, character(1))
      paste0(INDENT, paste(spans, collapse = ""))
    }, character(1))
    aln_rows <- aln_rows[nzchar(aln_rows)]

    paste(c(seq_line, aln_rows, ""), collapse = "\n")
  }, character(1))

  legend_html <- if (length(pep_list)) {
    badges <- paste(vapply(seq_along(pep_list), function(pi) {
      col <- PEP_COLS[(pi - 1L) %% length(PEP_COLS) + 1L]
      sprintf('<code class="seq-legend-badge" style="background:%s22;color:%s;border:1px solid %s55;">%s</code>',
              col, col, col, pep_list[[pi]])
    }, character(1)), collapse = " ")
    sprintf('<div class="seq-legend"><span class="text-muted small me-2">%d MS peptide%s - hover to see MS data:</span>%s</div>',
            length(pep_list), if (length(pep_list) > 1L) "s" else "", badges)
  } else {
    '<p class="text-muted small mb-1">No MS peptides identified for this ORF.</p>'
  }

  HTML(paste0(legend_html,
              '<div class="titan-protein-seq">',
              paste(blocks, collapse = "\n"),
              '</div>'))
}

biotype_bar <- function(df) {
  if (nrow(df) == 0) {
    return(
      plot_ly() %>%
        layout(
          annotations = list(list(text = "No data", x = 0.5, y = 0.5,
                                  xref = "paper", yref = "paper",
                                  showarrow = FALSE, font = list(size = 13, color = "#888"))),
          xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
          paper_bgcolor = "white", plot_bgcolor = "white"
        ) %>% config(displayModeBar = FALSE)
    )
  }
  d <- df %>%
    count(orf_biotype_single, name = "n") %>%
    arrange(desc(n)) %>%
    mutate(orf_biotype_single = factor(orf_biotype_single, levels = rev(orf_biotype_single)))
  cols <- coalesce(BIOTYPE_COLORS[as.character(d$orf_biotype_single)], "#95A5A6")
  plot_ly(d, x = ~n, y = ~orf_biotype_single, type = "bar", orientation = "h",
          marker = list(color = cols),
          text  = ~paste0(n, " (", round(100 * n / sum(n), 1), "%)"),
          textposition = "outside",
          cliponaxis = FALSE,
          hovertemplate = "%{y}: %{x} ORFs<extra></extra>") %>%
    layout(
      xaxis  = list(title = "Number of ORFs", showgrid = TRUE, gridcolor = "#EEF2F7",
                    autorange = TRUE),
      yaxis  = list(title = "", automargin = TRUE),
      paper_bgcolor = "white", plot_bgcolor = "white",
      margin = list(l = 5, r = 110, t = 10, b = 30),
      font   = list(family = "Inter", size = 12)
    ) %>% config(displayModeBar = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# THEME
# ─────────────────────────────────────────────────────────────────────────────

titan_theme <- bs_theme(
  version       = 5,
  bg            = "#FFFFFF",
  fg            = "#212529",
  primary       = "#28646E",
  secondary     = "#AADCFF",
  success       = "#00A555",
  warning       = "#D4850A",
  danger        = "#C0392B",
  info          = "#FFBEFF",
  font_scale    = 0.9,
  base_font     = font_google("Inter"),
  heading_font  = font_google("Inter"),
  code_font     = font_google("Fira Code"),
  "navbar-bg"               = "#28646E",
  "navbar-dark-color"       = "rgba(255,255,255,.9)",
  "navbar-dark-hover-color" = "#FFFFFF",
  "sidebar-bg"              = "#EFF7F8",
  "sidebar-border-color"    = "#C8D8DC",
  "card-border-color"       = "#DAE9EC",
  "card-cap-bg"             = "#F0F7F9"
)

# ─────────────────────────────────────────────────────────────────────────────
# UI HELPERS
# ─────────────────────────────────────────────────────────────────────────────

filtering_sidebar_ui <- function() {
  card(
    card_header(tags$span(icon("filter"), " Filtering"), class = "fw-semibold"),
    card_body(
      class = "px-2 pt-2 pb-1",

      tags$p(class = "text-muted small mb-1 fw-semibold", "Translation (Ribo-seq)"),
      numericInput("ppm_threshold", "PPM threshold", value = 1, min = 0, max = 200, step = 0.5),
      sliderInput("ppm_n_samples",
                  label = tags$span("Min. samples ≥ threshold",
                                    tags$small(class = "text-muted fw-normal ms-1",
                                               paste0("(max ", n_ribo_samples, ")"))),
                  min = 1, max = n_ribo_samples, value = floor(n_ribo_samples / 4),
                  step = 1, ticks = FALSE),

      hr(class = "my-2"),

      tags$p(class = "text-muted small mb-1 fw-semibold", "Expression (RNA-seq)"),
      numericInput("tpm_threshold", "TPM threshold", value = 1, min = 0, max = 200, step = 0.5),
      sliderInput("tpm_n_samples",
                  label = tags$span("Min. samples ≥ threshold",
                                    tags$small(class = "text-muted fw-normal ms-1",
                                               paste0("(max ", n_rna_samples, ")"))),
                  min = 1, max = n_rna_samples, value = floor(n_rna_samples / 4),
                  step = 1, ticks = FALSE),

      hr(class = "my-2"),

      tags$p(class = "text-muted small mb-1 fw-semibold", "ORF biotypes"),
      pickerInput(
        "biotype_filter", label = NULL,
        choices = biotypes, selected = biotypes, multiple = TRUE,
        options = pickerOptions(
          actionsBox = TRUE, selectedTextFormat = "count > 2",
          countSelectedText = "{0} / {1} biotypes",
          size = 8, container = "body", dropupAuto = FALSE
        )
      ),

      hr(class = "my-2"),
      actionButton("reset_filters", "Reset filters",
                   class = "btn-sm btn-outline-primary w-100 mb-2",
                   icon  = icon("rotate-left"))
    )
  )
}

scoring_sidebar_ui <- function() {
  groups <- unique(sapply(WEIGHT_META, `[[`, "group"))

  weight_sliders <- lapply(groups, function(grp) {
    metas <- Filter(function(m) m$group == grp, WEIGHT_META)
    tagList(
      tags$span(class = "weight-group-label", grp),
      lapply(metas, function(m) {
        is_pct <- m$id %in% c("w_pct_samples", "w_pct_transl")
        s <- sliderInput(
          m$id,
          label = tags$span(
            m$label,
            tags$br(),
            tags$span(class = "weight-hint", "↑ ", m$hint)
          ),
          min = if (is_pct) 0 else -1,
          max = 1, value = 0, step = 0.1, ticks = FALSE
        )
        if (is_pct) s else tags$div(class = "titan-bipolar-slider", s)
      })
    )
  })

  card(
    card_header(tags$span(icon("sliders"), " Scoring"), class = "fw-semibold"),
    card_body(
      class = "px-2 pt-2 pb-1",

      tags$p(class = "text-muted small mb-2 fw-semibold", "Preset"),

      lapply(names(PRESETS), function(name) {
        p <- PRESETS[[name]]
        btn_id <- paste0("preset_", gsub("[ -]", "_", tolower(name)))
        tags$button(
          id    = btn_id,
          class = "prio-preset-btn",
          style = paste0("border-left: 3px solid ", p$color, " !important;"),
          tags$span(class = "prio-preset-name", name),
          tags$span(class = "prio-preset-desc", p$label),
          onclick = sprintf(
            "Shiny.setInputValue('scoring_preset', '%s', {priority: 'event'});", name
          )
        )
      }),

      hr(class = "my-2"),
      tags$p(class = "text-muted small mb-1 fw-semibold", "Tumor specificity filter"),
      pickerInput(
        "prio_spec_filter", label = NULL,
        choices  = c("Tumor-only", "Tumor-enriched", "Non-specific"),
        selected = c("Tumor-only", "Tumor-enriched", "Non-specific"),
        multiple = TRUE,
        options  = pickerOptions(
          actionsBox = TRUE, container = "body", dropupAuto = FALSE,
          selectedTextFormat = "count > 2", countSelectedText = "{0} / {1} types"
        )
      ),
      checkboxInput("gtex_max_tpm_filter",
        label = tags$span(class = "small", "Exclude GTEx max tissue median TPM > 1"),
        value = FALSE
      ),
      hr(class = "my-2"),

      tags$p(class = "text-muted small mb-1 fw-semibold", "Custom weights",
             tags$span(id = "preset_modified_dot", style = "display:none; margin-left:6px;",
                       tags$span(style = "color:#D4850A; font-size:10px;", "● modified"))),

      weight_sliders,

      div(class = "weight-total-box",
          "Total weight: ", tags$strong(textOutput("total_weight", inline = TRUE)),
          tags$br(),
          tags$span(style = "color:#999; font-size:10px;",
                    "+1 rewards · −1 penalises · 0 neutral. Score normalised 0–100.")
      )
    )
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────

ui <- page_navbar(

  id           = "main_nav",
  title        = tags$span(
    tags$img(src = "titan_logo_white.svg", height = "26px",
             style = "margin-right:8px; vertical-align:middle;"),
    tags$span(" · Tumor Immunopeptidomics Target ANnotation",
              style = "font-weight:400; font-size:.85em; opacity:.85;")
  ),
  theme        = titan_theme,
  window_title = "TITAN",
  lang         = "en",
  bg           = "#28646E",
  header       = tags$head(
    tags$link(rel = "stylesheet", href = "titan.css"),
    tags$script(src = "titan-sliders.js")
  ),

  # ── Global sidebar ──────────────────────────────────────────────────────────
  sidebar = sidebar(
    width = 290,
    style = "overflow-y:auto; height:100%;",

    # Data tab: format descriptions
    conditionalPanel(
      condition = "input.main_nav === 'Data'",
      card(
        card_header(tags$span(icon("circle-info"), " Format guide"), class = "fw-semibold"),
        card_body(
          class = "px-2 pt-2 pb-1",
          tags$p(class = "titan-section-label", "ORF candidate RDS"),
          tags$p(class = "small text-muted mb-1",
                 "A named R list with the following elements:"),
          tags$ul(class = "small text-muted ps-3",
            tags$li(tags$code("$orf_table"),
                    " - data.frame of ncORF candidates with annotation columns"),
            tags$li(tags$code("$ribo_ppm_samples"),
                    " - numeric matrix (ORFs × samples); rownames = orf_id"),
            tags$li(tags$code("$ribo_sample_meta"),
                    " - data.frame with ", tags$code("sample_id"),
                    " and ", tags$code("condition"), " columns"),
            tags$li(tags$code("$rna_tpm_mat"),
                    " - optional matrix (genes × RNA-seq samples); rownames = gene_id ",
                    tags$em("without version suffix")),
            tags$li(tags$code("$prepared_on"),
                    " - optional POSIXct timestamp")
          ),
          hr(class = "my-2"),
          tags$p(class = "titan-section-label", "MS peptides file"),
          tags$p(class = "small text-muted mb-1",
                 "CSV or TSV with one column of peptide sequences (≥ 8 aa)."),
          tags$p(class = "small text-muted mb-1", "Auto-detected column names:"),
          tags$ul(class = "small text-muted ps-3",
            tags$li(tags$code("Peptide")),
            tags$li(tags$code("Sequence")),
            tags$li(tags$code("Annotated Sequence")),
            tags$li(tags$code("Modified Sequence"))
          ),
          tags$p(class = "small text-muted",
                 "Any additional columns (intensity, probability, etc.) are carried through.")
        )
      )
    ),

    # Filtering panel: Overview
    conditionalPanel(
      condition = "input.main_nav === 'Overview'",
      filtering_sidebar_ui()
    ),

    # Scoring panel: Prioritisation
    conditionalPanel(
      condition = "input.main_nav === 'Prioritisation'",
      scoring_sidebar_ui()
    ),

    # Package versions: About
    conditionalPanel(
      condition = "input.main_nav === 'About'",
      card(
        card_header(tags$span(icon("box-open"), " Package versions"), class = "fw-semibold"),
        card_body(
          class = "px-2 pt-1 pb-2",
          tags$ul(
            class = "list-unstyled mb-0 small",
            lapply(
              c("shiny", "bslib", "DT", "plotly",
                "dplyr", "tidyr", "stringr", "ggplot2", "shinyWidgets"),
              function(pkg) {
                tags$li(class = "d-flex justify-content-between py-1 border-bottom",
                  tags$code(pkg),
                  tags$span(class = "text-muted", as.character(packageVersion(pkg)))
                )
              }
            )
          )
        )
      )
    ),

    # ORF selector: ORF Detail
    conditionalPanel(
      condition = "input.main_nav === 'ORF Detail'",
      card(
        card_header(tags$span(icon("search"), " Select ORF"), class = "fw-semibold"),
        card_body(
          class = "px-2 pt-2 pb-1",
          selectizeInput("detail_orf_id", NULL, choices = NULL,
                         options = list(
                           placeholder    = "Type gene, biotype, or position…",
                           maxOptions     = 30,
                           dropdownParent = "body"
                         )),
          uiOutput("detail_orf_meta")
        )
      )
    )
  ),

  # ── Tab 1: Data ─────────────────────────────────────────────────────────────
  nav_panel(
    "Data", icon = icon("database"),

    div(
      class = "titan-data-wrap",

      div(
        class = "titan-hero",
        tags$img(src = "titan_logo_blue.svg", class = "titan-hero-logo", alt = "TITAN"),
        div(
          class = "titan-hero-text",
          tags$h1(class = "titan-hero-title", "Tumor Immunopeptidomics Target ANnotation"),
          tags$p(class = "titan-hero-subtitle", "van Heesch lab")
        )
      ),

      layout_columns(
        col_widths = c(6, 6),
        gap = "1rem",

        card(
          card_header(tags$span(icon("file-arrow-up"), " ORF candidates")),
          card_body(
            fileInput("user_rds_file", "Upload your own RDS",
                      accept = ".rds", buttonLabel = "Browse…",
                      placeholder = "titan_orf_table.rds"),
            div(class = "d-flex align-items-center gap-2 mb-3",
                tags$span(class = "text-muted small", "- or —"),
                actionButton("load_demo_rds", "Load demo dataset",
                             icon = icon("flask"), class = "btn-sm btn-outline-primary")),
            uiOutput("data_load_status")
          )
        ),

        card(
          card_header(tags$span(icon("vials"), " MS peptides")),
          card_body(
            fileInput("ms_file", "Upload MS results file",
                      accept = c(".csv", ".tsv", ".txt"),
                      buttonLabel = "Browse…",
                      placeholder = "peptides.csv / .tsv"),
            div(class = "d-flex align-items-center gap-2 mb-3",
                tags$span(class = "text-muted small", "- or —"),
                actionButton("load_demo_ms", "Load demo peptides",
                             icon = icon("flask"), class = "btn-sm btn-outline-primary")),
            uiOutput("ms_load_status"),
            uiOutput("col_selector"),
            uiOutput("match_summary_badge")
          )
        )
      ),

      uiOutput("start_section_ui")
    )
  ),

  # ── Tab 2: Overview ─────────────────────────────────────────────────────────
  nav_panel(
    "Overview", icon = icon("chart-bar"),
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      value_box("Total candidates",  textOutput("stat_total",        inline = TRUE),
                showcase = icon("microscope"),   theme = "primary"),
      value_box("Unique genes",      textOutput("stat_genes",        inline = TRUE),
                showcase = icon("dna"),          theme = "secondary"),
      value_box("Matched ORFs",      textOutput("stat_matched_orfs", inline = TRUE),
                showcase = icon("vials"),        theme = "success"),
      value_box("MS peptide hits",   textOutput("stat_matches",      inline = TRUE),
                showcase = icon("check-circle"), theme = "info")
    ),

    layout_columns(
      col_widths = c(6, 6),
      card(card_header("ORF biotype distribution (all candidates)"),
           card_body(class = "p-2", plotlyOutput("plot_biotype",         height = "250px"))),
      card(card_header("ORF biotype distribution (peptide-matched ORFs)"),
           card_body(class = "p-2", plotlyOutput("plot_biotype_matched", height = "250px")))
    ),

    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Translation vs expression",
          materialSwitch("transl_expr_matched_only",
                         label  = tags$span(class = "text-muted fw-normal", "Matched only"),
                         value  = FALSE, status = "primary", right = TRUE, inline = TRUE)
        ),
        card_body(class = "p-2", plotlyOutput("plot_transl_expr", height = "350px"))
      ),
      card(card_header("PPM distribution by biotype (matched vs unmatched)"),
           card_body(class = "p-2", plotlyOutput("plot_ppm_dist",        height = "350px")))
    )
  ),

  # ── Tab 3: Prioritisation ───────────────────────────────────────────────────
  nav_panel(
    "Prioritisation", icon = icon("ranking-star"),

    conditionalPanel(
      condition = "output.has_peptides == false",
      div(class = "d-flex flex-column align-items-center justify-content-center mt-5 text-muted",
          icon("upload", class = "fa-3x mb-3"),
          tags$h5("No peptides loaded"),
          tags$p("Go to the Data tab to upload MS results or load the demo peptides."))
    ),

    conditionalPanel(
      condition = "output.has_peptides == true && output.has_started == false",
      div(class = "d-flex flex-column align-items-center justify-content-center mt-5 text-muted",
          icon("play-circle", class = "fa-3x mb-3"),
          tags$h5("Ready to start"),
          tags$p("Press START on the Data tab to begin prioritisation."))
    ),

    conditionalPanel(
      condition = "output.has_started == true",

      layout_columns(
        col_widths = c(3, 3, 3, 3),
        value_box("Candidates ranked", textOutput("stat_prio_total",  inline = TRUE),
                  showcase = icon("list-ol"),  theme = "primary"),
        value_box("MS peptides",       textOutput("stat_prio_pep",    inline = TRUE),
                  showcase = icon("vials"),    theme = "secondary"),
        value_box("Top score",         textOutput("stat_prio_top",    inline = TRUE),
                  showcase = icon("trophy"),   theme = "success"),
        value_box("Active preset",     textOutput("stat_prio_preset", inline = TRUE),
                  showcase = icon("sliders"),  theme = "info")
      ),

      card(
        class = "titan-priority-card",
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Ranked candidates",
          tags$div(
            class = "d-flex gap-1",
            downloadButton("dl_params",    "Export parameters",
                           class = "btn-sm btn-outline-primary"),
            downloadButton("dl_priority",  "Export ranked",
                           class = "btn-sm btn-outline-primary"),
            downloadButton("dl_selected",  "Export selection",
                           class = "btn-sm btn-outline-success")
          )
        ),
        card_body(class = "p-0 titan-priority-body", DTOutput("tbl_priority"))
      ),

      uiOutput("priority_detail_panel")
    )
  ),

  # ── Tab 4: ORF Detail ───────────────────────────────────────────────────────
  nav_panel(
    "ORF Detail", icon = icon("circle-info"),

    conditionalPanel(
      condition = "output.has_started == false",
      div(class = "d-flex flex-column align-items-center justify-content-center mt-5 text-muted",
          icon("play-circle", class = "fa-3x mb-3"),
          tags$h5("Not started"),
          tags$p("Press START on the Data tab to begin."))
    ),

    conditionalPanel(
      condition = "output.has_started == true",

      card(
        card_header("Protein sequence & peptide matches"),
        card_body(class = "p-2", uiOutput("detail_protein_seq_ui"))
      ),

      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Per-sample translation (ribo-seq PPM)"),
             card_body(class = "p-2", plotlyOutput("plot_detail_ribo", height = "300px"))),
        card(
          card_header("Prioritisation metrics"),
          card_body(
            class = "p-2",
            style = "max-height:360px; overflow-y:auto;",
            uiOutput("detail_metric_badges"),
            hr(class = "my-2"),
            tableOutput("detail_metric_table")
          )
        )
      )
    )
  ),

  # ── Tab 5: About ────────────────────────────────────────────────────────────
  nav_panel(
    "About", icon = icon("circle-question"),
    card(
      card_body(
        tags$h4("TITAN - Tumor Immunopeptidomics Target ANnotation"),
        tags$p("TITAN integrates ribo-seq, RNA-seq, and external databases to prioritise",
               " non-canonical ORF-derived peptide candidates for cancer immunotherapy."),
        tags$hr(),
        tags$h6("Data sources"),
        tags$ul(
          tags$li(tags$b("ncORF candidates:"),  " immunopeptidomics-filtered ORFs from the ribo-seq pipeline"),
          tags$li(tags$b("Target translation:"), " ribo-seq (psites per million)"),
          tags$li(tags$b("Ribocrypt external:"), " cross-study ribo-seq (primary tissues and cell lines)"),
          tags$li(tags$b("Target expression:"),  " RNA-seq (Salmon/tximport gene-level TPM)"),
          tags$li(tags$b("GTEx DE:"),            " DESeq2 GTEx (28 tissues); tumor-specificity classification"),
          tags$li(tags$b("TCGA:"),               " 257 tumor + 257 normal/peritumoral samples (Salmon TPM)")
        ),
        tags$hr(),
        tags$h6("Scoring dimensions"),
        tags$ul(lapply(WEIGHT_META, function(m) {
          tags$li(tags$b(m$label), " - ", m$hint, tags$span(class="text-muted small ms-1", paste0("(", m$group, ")")))
        })),
        tags$hr(),
        uiOutput("about_data_info")
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Reactive data (NULL until user loads; replaced on upload) ───────────────
  app_data_rv <- reactiveVal(NULL)

  orf_table_rv <- reactive({
    req(app_data_rv())
    tbl <- app_data_rv()$orf_table
    if (!"gene_id_clean" %in% colnames(tbl))
      tbl$gene_id_clean <- sub("\\..*", "", tbl$gene_id)
    tbl
  })
  ribo_ppm_rv           <- reactive({ req(app_data_rv()); app_data_rv()$ribo_ppm_samples })
  ribo_meta_rv          <- reactive({ req(app_data_rv()); app_data_rv()$ribo_sample_meta })
  rna_tpm_rv            <- reactive({ req(app_data_rv()); app_data_rv()$rna_tpm_mat })
  rna_meta_rv           <- reactive({ req(app_data_rv()); app_data_rv()$rna_sample_meta })
  gtex_tpm_rv           <- reactive({ req(app_data_rv()); app_data_rv()$gtex_tpm_mat })
  gtex_meta_rv          <- reactive({ req(app_data_rv()); app_data_rv()$gtex_sample_meta })
  tcga_tpm_rv           <- reactive({ req(app_data_rv()); app_data_rv()$tcga_tpm_mat })
  tcga_meta_rv          <- reactive({ req(app_data_rv()); app_data_rv()$tcga_sample_meta })
  ribocrypt_mat_rv      <- reactive({ req(app_data_rv()); app_data_rv()$ribocrypt_mat })
  ribocrypt_smeta_rv    <- reactive({ req(app_data_rv()); app_data_rv()$ribocrypt_sample_meta })

  # Update sidebar controls + ORF selector when data changes
  observeEvent(app_data_rv(), {
    bios  <- sort(unique(app_data_rv()$orf_table$orf_biotype_single))
    n_r   <- ncol(app_data_rv()$ribo_ppm_samples)
    n_rna <- if (!is.null(app_data_rv()$rna_tpm_mat))
               ncol(app_data_rv()$rna_tpm_mat)
             else
               nrow(app_data_rv()$rna_sample_meta)
    updateSliderInput(session, "ppm_n_samples", max = n_r,   value = floor(n_r   / 4),
      label = paste0("Min. samples ≥ threshold (max ", n_r,   ")"))
    updateSliderInput(session, "tpm_n_samples", max = n_rna, value = floor(n_rna / 4),
      label = paste0("Min. samples ≥ threshold (max ", n_rna, ")"))
    updatePickerInput(session, "biotype_filter", choices = bios, selected = bios)
  }, ignoreInit = TRUE)

  # ORF Detail selector: only matched ORFs that pass Overview + Prioritisation filters
  observeEvent({
    safe_matched_data()
    input$prio_spec_filter
    input$gtex_max_tpm_filter
  }, {
    md <- safe_matched_data()
    if (is.null(md) || nrow(md) == 0L) {
      # server = FALSE for empty choices - server = TRUE with character(0) sends a
      # malformed AJAX response that triggers "Cannot read properties of undefined"
      updateSelectizeInput(session, "detail_orf_id", choices = character(0), server = FALSE)
      return()
    }
    # Use spec-filtered priority table to restrict choices.
    # orf_ids holds all ORF IDs per group (best + co-identified child ORFs).
    ppd <- tryCatch(prio_table_df(), error = function(e) NULL)
    matched_ids <- if (!is.null(ppd) && nrow(ppd) > 0L) {
      unique(trimws(unlist(strsplit(ppd$orf_ids, ",\\s*", perl = TRUE))))
    } else {
      unique(md$orf_id)
    }
    tbl <- orf_table_rv()[orf_table_rv()$orf_id %in% matched_ids, ]
    if (nrow(tbl) == 0L) {
      updateSelectizeInput(session, "detail_orf_id", choices = character(0), server = FALSE)
      return()
    }
    choices <- setNames(tbl$orf_id, orf_id_labels(tbl))
    top_oid <- if (!is.null(ppd) && nrow(ppd) > 0L) ppd$.orf_id[1L] else choices[[1L]]
    updateSelectizeInput(session, "detail_orf_id",
      choices = choices, selected = top_oid, server = TRUE)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  # ── Data tab: load handlers ─────────────────────────────────────────────────
  observeEvent(input$load_demo_rds, {
    withProgress(message = "Loading demo dataset…", value = 0.2, {
      setProgress(0.6, detail = "Reading RDS…")
      dat <- readRDS("data/titan_orf_table.rds")
      setProgress(0.95)
      app_data_rv(dat)
    })
  })

  observeEvent(input$user_rds_file, {
    req(input$user_rds_file)
    withProgress(message = "Loading dataset…", value = 0.2, {
      setProgress(0.6, detail = "Reading RDS…")
      dat <- tryCatch(readRDS(input$user_rds_file$datapath), error = function(e) NULL)
      setProgress(0.95)
      if (is.null(dat) || is.null(dat$orf_table)) {
        showNotification("Invalid RDS: expected a list with $orf_table.", type = "error")
        return()
      }
      message("[DEBUG] RDS loaded. Top-level names: ", paste(names(dat), collapse = ", "))
      message("[DEBUG] orf_table columns: ", paste(colnames(dat$orf_table), collapse = ", "))
      gtex_col <- dat$orf_table$GTEX_tissues_q3_gt1
      message("[DEBUG] GTEX_tissues_q3_gt1 in orf_table: ", !is.null(gtex_col),
              " non-NA: ", sum(!is.na(gtex_col)),
              " example: ", head(gtex_col[!is.na(gtex_col)], 1))
      app_data_rv(dat)
    })
  })

  output$data_load_status <- renderUI({
    d <- app_data_rv()
    if (is.null(d))
      return(tags$span(class = "text-muted small", "No dataset loaded."))
    n  <- nrow(d$orf_table)
    nr <- ncol(d$ribo_ppm_samples)
    tagList(
      div(class = "d-flex align-items-center gap-2 flex-wrap",
          tags$span(class = "badge bg-secondary",
                    paste0(formatC(n, big.mark = ","), " ORFs")),
          tags$span(class = "badge bg-primary",
                    paste0(nr, " samples")),
          actionButton("clear_rds", "Clear ORF data", icon = icon("trash"),
                       class = "btn-sm btn-outline-danger ms-auto",
                       title = "Clear dataset")),
      if (!is.null(d$prepared_on))
        tags$small(class = "text-muted d-block mt-1",
                   format(d$prepared_on, "%Y-%m-%d"))
    )
  })

  output$ms_load_status <- renderUI({
    ms <- tryCatch(ms_data(), error = function(e) NULL)
    if (is.null(ms)) return(NULL)
    n_pep   <- nrow(ms)
    is_demo <- is.null(user_ms_rv())
    div(class = "d-flex align-items-center gap-2 mb-2",
        tags$span(class = "badge bg-secondary",
                  paste0(formatC(n_pep, big.mark = ","), " peptides")),
        if (is_demo) tags$span(class = "badge bg-info", "demo"),
        actionButton("clear_ms", "Clear MS data", icon = icon("trash"),
                     class = "btn-sm btn-outline-danger ms-auto",
                     title = "Clear MS data"))
  })

  observeEvent(input$clear_rds, {
    app_data_rv(NULL)
    all_matches_rv(NULL)
    started_rv(FALSE)
  })

  observeEvent(input$clear_ms, {
    user_ms_rv(NULL)
    demo_ms_rv(NULL)
    all_matches_rv(NULL)
    started_rv(FALSE)
  })

  started_rv    <- reactiveVal(FALSE)
  all_matches_rv <- reactiveVal(NULL)

  observeEvent(input$start_titan, {
    req(ms_peptides(), orf_table_rv(), ms_meta())
    withProgress(message = "Matching peptides to ORFs…", value = 0.1, {
      hits <- match_peptides(ms_peptides(), orf_table_rv())
      # For peptides that match an ORF-annotated (canonical) entry AND non-canonical
      # entries, discard the non-canonical hits — they are canonical peptide evidence.
      if (!is.null(hits) && nrow(hits) > 0L) {
        hits <- hits %>%
          group_by(matched_peptide) %>%
          filter(
            if (any(orf_biotype_single == "ORF-annotated", na.rm = TRUE))
              orf_biotype_single == "ORF-annotated"
            else
              TRUE
          ) %>%
          ungroup()
      }
      # Gencode cross-match: annotate in-house hits and add Gencode-only rows
      if (!is.null(gencode_orf_tbl)) {
        gc_all <- match_peptides(ms_peptides(), gencode_orf_tbl)
        if (!is.null(gc_all) && nrow(gc_all) > 0L) {
          # Case (a): build per-peptide summary of matching Gencode ORFs
          gc_summary <- gc_all %>%
            group_by(matched_peptide) %>%
            summarise(
              gencode_match_ids = paste(
                sprintf("%s (%s, %s)", orf_id, gene_name, orf_biotype_single),
                collapse = "; "
              ),
              .groups = "drop"
            )
          if (!is.null(hits) && nrow(hits) > 0L) {
            hits <- hits %>%
              left_join(gc_summary, by = "matched_peptide") %>%
              mutate(gencode_match_ids = replace_na(gencode_match_ids, ""),
                     gencode_only       = FALSE)
            # Case (b): peptides with Gencode hits but no in-house hit → new rows
            gc_only_peps <- setdiff(unique(gc_all$matched_peptide), unique(hits$matched_peptide))
          } else {
            gc_only_peps <- unique(gc_all$matched_peptide)
          }
          if (length(gc_only_peps) > 0L) {
            gc_only_rows <- gc_all %>%
              filter(matched_peptide %in% gc_only_peps) %>%
              mutate(gencode_match_ids = "", gencode_only = TRUE)
            # Populate gene-level expression / GTEx / TCGA metrics by borrowing from
            # any in-house ORF of the same gene (these columns are gene-level, not ORF-level)
            expr_cols <- intersect(
              c("target_expression_num_samples", "target_expression_pct_samples",
                "target_expression_median_TPM", "target_expression_max_TPM",
                "GTEX_max_median_TPM", "GTEX_median_TPM", "GTEX_DE_sig_in_all",
                "GTEX_tumor_only", "GTEX_tumor_enriched", "GTEX_tissues_q3_gt1",
                "TCGA_tumor_num_samples", "TCGA_tumor_pct_samples",
                "TCGA_tumor_median_TPM", "TCGA_tumor_max_TPM",
                "TCGA_normal_num_samples", "TCGA_normal_pct_samples",
                "TCGA_normal_median_TPM", "TCGA_normal_max_TPM"),
              colnames(orf_table_rv())
            )
            gene_expr <- orf_table_rv() %>%
              group_by(gene_id_clean) %>%
              summarise(across(all_of(expr_cols), first), .groups = "drop")
            gc_only_rows <- gc_only_rows %>%
              left_join(gene_expr, by = "gene_id_clean")
            hits <- bind_rows(hits, gc_only_rows)
          }
        } else if (!is.null(hits) && nrow(hits) > 0L) {
          hits <- hits %>% mutate(gencode_match_ids = "", gencode_only = FALSE)
        }
      } else if (!is.null(hits) && nrow(hits) > 0L) {
        hits <- hits %>% mutate(gencode_match_ids = "", gencode_only = FALSE)
      }
      setProgress(0.85, detail = "Joining MS metadata…")
      result <- if (!is.null(hits) && nrow(hits) > 0L) {
        left_join(hits, ms_meta(), by = "matched_peptide")
      } else {
        data.frame(orf_id = character(0), matched_peptide = character(0))
      }
      all_matches_rv(result)
    })
    started_rv(TRUE)
    nav_select("main_nav", "Overview")
  })
  output$has_started <- reactive({ isTRUE(started_rv()) })
  outputOptions(output, "has_started", suspendWhenHidden = FALSE)

  output$start_section_ui <- renderUI({
    rds_ok <- !is.null(app_data_rv())
    ms_ok  <- !is.null(user_ms_rv()) || !is.null(demo_ms_rv())
    if (rds_ok && ms_ok && !started_rv()) {
      div(
        class = "mt-4 text-center",
        div(class = "d-flex justify-content-center gap-2 mb-3",
            tags$span(class = "badge bg-primary px-3 py-2",
                      icon("check"), " ORF dataset ready"),
            tags$span(class = "badge bg-primary px-3 py-2",
                      icon("check"), " MS peptides ready")),
        actionButton("start_titan", " START",
                     icon  = icon("play-circle"),
                     class = "btn-lg btn-outline-primary px-5 py-2 fw-bold")
      )
    } else if (started_rv()) {
      div(
        class = "mt-4 text-center text-success",
        icon("circle-check", class = "fa-2x"),
        tags$p(class = "mt-2 mb-0 fw-semibold", "Running - use the tabs above to explore.")
      )
    } else {
      # Show which piece is still missing
      div(
        class = "mt-4 text-center",
        div(class = "d-flex justify-content-center gap-2 mb-2",
            if (rds_ok) tags$span(class = "badge bg-primary px-3 py-2", icon("check"), " ORF dataset")
            else        tags$span(class = "badge bg-secondary px-3 py-2", icon("circle"), " ORF dataset"),
            if (ms_ok)  tags$span(class = "badge bg-primary px-3 py-2", icon("check"), " MS peptides")
            else        tags$span(class = "badge bg-secondary px-3 py-2", icon("circle"), " MS peptides")),
        tags$p(class = "text-muted small mb-0",
               "Load both an ORF dataset and MS peptides to enable prioritisation.")
      )
    }
  })

  # ── MS peptides ─────────────────────────────────────────────────────────────
  demo_ms_rv <- reactiveVal(NULL)
  user_ms_rv <- reactiveVal(NULL)

  observeEvent(input$load_demo_ms, {
    withProgress(message = "Loading demo peptides…", value = 0.5, {
      demo_ms_rv(read.delim("data/demo_ms_peptides.tsv", sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE))
    })
  })

  observeEvent(input$ms_file, {
    req(input$ms_file)
    # Sniff separator from first line - don't trust the file extension
    first_line <- tryCatch(readLines(input$ms_file$datapath, n = 1L, warn = FALSE),
                           error = function(e) "")
    n_tabs   <- nchar(first_line) - nchar(gsub("\t", "", first_line, fixed = TRUE))
    n_commas <- nchar(first_line) - nchar(gsub(",",  "", first_line, fixed = TRUE))
    sep <- if (n_tabs >= n_commas) "\t" else ","
    dat <- tryCatch(
      read.delim(input$ms_file$datapath, sep = sep,
                 stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) {
        showNotification(paste("Could not read MS file:", conditionMessage(e)), type = "error")
        NULL
      }
    )
    if (!is.null(dat)) user_ms_rv(dat)
  })

  ms_data <- reactive({
    if (!is.null(user_ms_rv())) return(user_ms_rv())
    req(demo_ms_rv())
    demo_ms_rv()
  })

  PEPTIDE_COL_CANDIDATES <- c("Peptide", "Sequence", "peptide", "sequence",
                               "Annotated Sequence", "Modified Sequence")
  auto_pep_col <- reactive({
    req(ms_data())
    m <- intersect(PEPTIDE_COL_CANDIDATES, colnames(ms_data()))
    if (length(m)) m[1] else colnames(ms_data())[1]
  })

  output$col_selector <- renderUI({
    req(ms_data())
    selectInput("pep_col", "Peptide sequence column",
                choices = colnames(ms_data()), selected = auto_pep_col())
  })

  ms_peptides <- reactive({
    req(ms_data(), input$pep_col)
    unique(trimws(ms_data()[[input$pep_col]]))
  })

  # Columns tried in order for best-PSM-per-peptide deduplication (higher = better)
  PSM_QUALITY_COLS <- c("Probability", "Hyperscore", "SpectralSim", "Score")

  ms_meta <- reactive({
    req(ms_data(), input$pep_col)
    cols      <- colnames(ms_data())
    score_col <- intersect(PSM_QUALITY_COLS, cols)[1]  # first available, NA if none
    ms_data() %>%
      rename(matched_peptide = !!input$pep_col) %>%
      group_by(matched_peptide) %>%
      slice_max(
        order_by = if (!is.na(score_col)) .data[[score_col]] else row_number(),
        n = 1, with_ties = FALSE
      ) %>%
      ungroup()
  })

  matched_data <- reactive({
    m <- all_matches_rv()
    if (is.null(m) || nrow(m) == 0) return(NULL)
    fd <- filtered_data()
    if (is.null(fd) || nrow(fd) == 0) return(NULL)
    gc_pass <- if ("gencode_only" %in% colnames(m)) m$gencode_only %in% TRUE else FALSE
    m[m$orf_id %in% fd$orf_id | gc_pass, , drop = FALSE]
  })

  safe_matched_data <- reactive({
    tryCatch(matched_data(), error = function(e) NULL)
  })

  output$has_peptides <- reactive({
    !is.null(user_ms_rv()) || !is.null(demo_ms_rv())
  })
  outputOptions(output, "has_peptides", suspendWhenHidden = FALSE)

  output$match_summary_badge <- renderUI({
    req(matched_data())
    n_pep  <- n_distinct(matched_data()$matched_peptide)
    n_orfs <- n_distinct(matched_data()$orf_id)
    div(class = "mt-1",
        pill_badge(paste(n_pep,  "peptides matched"), "success"),
        pill_badge(paste(n_orfs, "ORF hits"),         "primary"))
  })

  # ── Filter reactive (debounced) ─────────────────────────────────────────────
  filtered_data_raw <- reactive({
    df <- orf_table_rv()

    ppm_thr <- max(0, input$ppm_threshold %||% 1)
    ppm_n   <- max(1, input$ppm_n_samples  %||% 1)
    tpm_thr <- max(0, input$tpm_threshold  %||% 1)
    tpm_n   <- max(1, input$tpm_n_samples  %||% 1)

    selected_bios <- input$biotype_filter %||% sort(unique(df$orf_biotype_single))
    if (length(selected_bios) < length(unique(df$orf_biotype_single)))
      df <- filter(df, orf_biotype_single %in% selected_bios)

    ribo_mat <- ribo_ppm_rv()
    in_mat      <- df$orf_id %in% rownames(ribo_mat)
    n_above_ppm <- integer(nrow(df))
    if (any(in_mat))
      n_above_ppm[in_mat] <- as.integer(rowSums(
        ribo_mat[df$orf_id[in_mat], , drop = FALSE] >= ppm_thr, na.rm = TRUE
      ))
    df <- df[n_above_ppm >= ppm_n, ]

    rna_mat <- rna_tpm_rv()
    if (!is.null(rna_mat)) {
      gene_n     <- rowSums(rna_mat >= tpm_thr, na.rm = TRUE)
      pass_genes <- names(gene_n)[gene_n >= tpm_n]
      df <- filter(df, gene_id_clean %in% pass_genes)
    } else {
      df <- filter(df, is.na(target_expression_num_samples) |
                         target_expression_num_samples >= tpm_n)
    }
    df
  })

  filtered_data <- debounce(filtered_data_raw, 500)

  observeEvent(input$reset_filters, {
    d <- app_data_rv()
    if (is.null(d)) return()
    bios  <- sort(unique(d$orf_table$orf_biotype_single))
    n_r   <- ncol(d$ribo_ppm_samples)
    n_rna <- if (!is.null(d$rna_tpm_mat)) ncol(d$rna_tpm_mat) else nrow(d$rna_sample_meta)
    updateNumericInput(session, "ppm_threshold", value = 1)
    updateSliderInput(session,  "ppm_n_samples", value = floor(n_r / 4))
    updateNumericInput(session, "tpm_threshold", value = 1)
    updateSliderInput(session,  "tpm_n_samples", value = floor(n_rna / 4))
    updatePickerInput(session,  "biotype_filter", selected = bios)
  })

  # ── Scoring weights ─────────────────────────────────────────────────────────
  current_weights <- reactive({
    w <- sapply(WEIGHT_META, function(m) {
      v <- input[[m$id]]
      if (is.null(v)) 0 else v
    })
    setNames(as.list(w), sapply(WEIGHT_META, `[[`, "id"))
  })

  active_preset <- reactiveVal(NULL)

  observeEvent(input$scoring_preset, {
    preset <- input$scoring_preset
    if (!preset %in% names(PRESETS)) return()
    pw <- preset_weights(preset)
    for (m in WEIGHT_META) updateSliderInput(session, m$id, value = pw[[m$id]])
    active_preset(preset)
  }, ignoreInit = TRUE)

  lapply(WEIGHT_META, function(m) {
    observeEvent(input[[m$id]], { active_preset(NULL) }, ignoreInit = TRUE)
  })

  output$total_weight <- renderText({ round(sum(abs(unlist(current_weights()))), 1) })
  output$stat_prio_preset <- renderText({
    ap <- active_preset(); if (is.null(ap)) "Custom" else ap
  })

  # ── Scored/ranked data ──────────────────────────────────────────────────────
  prioritised_data <- reactive({
    req(matched_data())
    # Collapse to one row per ORF - all scoring columns are ORF-level (same for
    # all peptides sharing an orf_id), so first() is safe for those columns.
    per_orf <- matched_data() %>%
      group_by(orf_id) %>%
      summarise(
        across(c(gene_name, orf_biotype_single, chr, orf_start, orf_end, strand,
                 protein_length, start_codon, gene_id, gene_biotype, gene_id_clean,
                 starts_with("target_"), starts_with("GTEX_"), starts_with("TCGA_"),
                 starts_with("ribocrypt_")),
               first),
        n_peptides       = n_distinct(matched_peptide),
        matched_peptides = paste(sort(unique(matched_peptide)), collapse = ", "),
        .groups = "drop"
      )
    score_candidates(per_orf, current_weights()) %>%
      arrange(desc(priority_score)) %>%
      mutate(.row_id = row_number())
  })

  gene_prioritised_data <- reactive({
    req(prioritised_data())

    # Group by (gene, biotype, exact peptide set): ORFs that are identified by the
    # same peptides collapse into one row; different peptide evidence = separate rows.
    # prioritised_data() is already score-sorted desc, so group_orfs[1] = best ORF.
    # group_modify strips grouping columns from group_orfs; they come back via key.
    per_group <- prioritised_data() %>%
      group_by(gene_id, orf_biotype_single, matched_peptides) %>%
      group_modify(function(group_orfs, key) {
        n_grp <- nrow(group_orfs)
        child_html <- if (n_grp > 1L) {
          child_orfs <- group_orfs[-1L, , drop = FALSE] %>%
            mutate(orf_biotype_single = key$orf_biotype_single,
                   matched_peptides   = key$matched_peptides)
          make_child_html(child_orfs)
        } else ""
        group_orfs[1L, , drop = FALSE] %>%
          mutate(.child_html = child_html,
                 n_orfs      = n_grp,
                 orf_ids     = paste(group_orfs$orf_id, collapse = ", "))
      }) %>%
      ungroup()

    per_group %>%
      arrange(desc(priority_score)) %>%
      mutate(.row_id = row_number())
  })

  # ── Overview stats ───────────────────────────────────────────────────────────
  output$stat_total <- renderText(formatC(nrow(filtered_data()), big.mark = ","))
  output$stat_genes <- renderText(formatC(n_distinct(filtered_data()$gene_id), big.mark = ","))
  output$stat_matched_orfs <- renderText({
    if (is.null(matched_data())) return("—")
    formatC(n_distinct(matched_data()$orf_id), big.mark = ",")
  })
  output$stat_matches <- renderText({
    if (is.null(matched_data())) return("—")
    formatC(nrow(matched_data()), big.mark = ",")
  })

  # ── Overview plots ───────────────────────────────────────────────────────────
  output$plot_biotype <- renderPlotly({ biotype_bar(filtered_data()) })

  output$plot_biotype_matched <- renderPlotly({
    if (is.null(matched_data()) || nrow(matched_data()) == 0) {
      return(
        plot_ly() %>%
          layout(
            annotations = list(list(
              text = "Load MS peptides to see<br>matched ORF biotype distribution",
              x = 0.5, y = 0.5, xref = "paper", yref = "paper",
              showarrow = FALSE, font = list(size = 13, color = "#888")
            )),
            xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
            paper_bgcolor = "white", plot_bgcolor = "white"
          ) %>% config(displayModeBar = FALSE)
      )
    }
    matched_orfs_df <- filtered_data() %>%
      filter(orf_id %in% unique(matched_data()$orf_id))
    biotype_bar(matched_orfs_df)
  })

  output$plot_transl_expr <- renderPlotly({
    df <- filtered_data() %>%
      filter(!is.na(target_translation_median_PPM), !is.na(target_expression_median_TPM))

    matched_ids <- if (!is.null(safe_matched_data())) unique(safe_matched_data()$orf_id) else character(0)
    if (isTRUE(input$transl_expr_matched_only)) {
      if (length(matched_ids) == 0)
        return(plot_ly() %>%
                 layout(annotations = list(list(
                   text = "No matched peptides loaded",
                   x = 0.5, y = 0.5, xref = "paper", yref = "paper",
                   showarrow = FALSE, font = list(size = 13, color = "#888")
                 )), xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                 paper_bgcolor = "white", plot_bgcolor = "white") %>%
                 config(displayModeBar = FALSE))
      df <- filter(df, orf_id %in% matched_ids)
    }

    df <- df %>%
      mutate(
        is_matched = orf_id %in% matched_ids,
        tip = paste0("<b>", gene_name, "</b><br>", orf_biotype_single, "<br>",
                     "PPM: ", round(target_translation_median_PPM, 2), "<br>",
                     "TPM: ", round(target_expression_median_TPM,  2),
                     if_else(is_matched, "<br><b style='color:#1B4F72'>✓ peptide match</b>", ""))
      )

    has_matches <- length(matched_ids) > 0

    scatter_layout <- function(p)
      p %>%
        layout(
          xaxis  = list(title = "log₁₀(Median TPM + 0.1) - RNA-seq",
                        showgrid = TRUE, gridcolor = "#EEF2F7"),
          yaxis  = list(title = "log₁₀(Median PPM + 0.1) - Ribo-seq",
                        showgrid = TRUE, gridcolor = "#EEF2F7"),
          paper_bgcolor = "white", plot_bgcolor = "white",
          legend = list(title = list(text = "Biotype"), x = 1.01, y = 1),
          margin = list(l = 10, r = 5, t = 10, b = 40),
          font   = list(family = "Inter", size = 12)
        ) %>%
        config(displayModeBar = FALSE)

    if (!has_matches) {
      return(scatter_layout(
        plot_ly(df,
                x = ~log10(target_expression_median_TPM + 0.1),
                y = ~log10(target_translation_median_PPM + 0.1),
                type = "scatter", mode = "markers",
                color = ~orf_biotype_single, colors = BIOTYPE_COLORS,
                marker = list(size = 5, opacity = 0.65, line = list(width = 0)),
                text = ~tip, hovertemplate = "%{text}<extra></extra>")
      ))
    }

    p    <- plot_ly()
    d_no  <- filter(df, !is_matched)
    d_yes <- filter(df, is_matched)

    for (bio in unique(as.character(d_no$orf_biotype_single))) {
      d_bio <- filter(d_no, orf_biotype_single == bio)
      col   <- unname(BIOTYPE_COLORS[bio])
      if (length(col) == 0 || is.na(col)) col <- "#95A5A6"
      p <- p %>% add_trace(
        data = d_bio, type = "scatter", mode = "markers",
        x = ~log10(target_expression_median_TPM + 0.1),
        y = ~log10(target_translation_median_PPM + 0.1),
        name = bio, legendgroup = bio,
        showlegend = !(bio %in% unique(as.character(d_yes$orf_biotype_single))),
        marker = list(color = col, size = 5, opacity = 0.3, line = list(width = 0)),
        text = ~tip, hovertemplate = "%{text}<extra></extra>"
      )
    }

    for (bio in unique(as.character(d_yes$orf_biotype_single))) {
      d_bio <- filter(d_yes, orf_biotype_single == bio)
      col   <- unname(BIOTYPE_COLORS[bio])
      if (length(col) == 0 || is.na(col)) col <- "#95A5A6"
      p <- p %>% add_trace(
        data = d_bio, type = "scatter", mode = "markers",
        x = ~log10(target_expression_median_TPM + 0.1),
        y = ~log10(target_translation_median_PPM + 0.1),
        name = bio, legendgroup = bio, showlegend = TRUE,
        marker = list(color = col, size = 9, opacity = 1,
                      line = list(color = "#ffffff", width = 1.5)),
        text = ~tip, hovertemplate = "%{text}<extra></extra>"
      )
    }

    scatter_layout(p)
  })

  output$plot_ppm_dist <- renderPlotly({
    df_all <- filtered_data() %>%
      filter(!is.na(target_translation_median_PPM)) %>%
      mutate(log_ppm = log10(target_translation_median_PPM + 0.1))

    matched_ids <- if (!is.null(safe_matched_data())) unique(safe_matched_data()$orf_id) else character(0)

    df_all <- df_all %>%
      mutate(match_group = factor(
        if_else(orf_id %in% matched_ids, "Peptide match", "No match"),
        levels = c("No match", "Peptide match")
      ))
    bios <- df_all %>% count(orf_biotype_single) %>% filter(n >= 3) %>% pull(orf_biotype_single)
    df_all <- filter(df_all, orf_biotype_single %in% bios) %>%
      mutate(orf_biotype_single = factor(orf_biotype_single, levels = bios))

    if (nrow(df_all) == 0)
      return(plot_ly() %>% layout(title = "No data", paper_bgcolor = "white") %>%
               config(displayModeBar = FALSE))

    has_matches <- length(matched_ids) > 0
    p <- plot_ly()
    d_no  <- filter(df_all, match_group == "No match")
    d_yes <- filter(df_all, match_group == "Peptide match")

    if (!has_matches) {
      for (bio in bios) {
        d_bio <- filter(d_no, orf_biotype_single == bio)
        if (nrow(d_bio) < 3) next
        col <- coalesce(BIOTYPE_COLORS[bio], "#95A5A6")
        p <- p %>% add_trace(
          type = "violin", x = d_bio$orf_biotype_single, y = d_bio$log_ppm,
          name = bio, legendgroup = "unmatched", showlegend = (bio == bios[1]),
          fillcolor = paste0(col, "99"), line = list(color = col, width = 1.2),
          box = list(visible = TRUE), points = FALSE,
          hovertemplate = "%{x}<br>log₁₀PPM: %{y:.2f}<extra></extra>"
        )
      }
    } else {
      if (nrow(d_no) > 0)
        p <- p %>% add_trace(
          type = "violin", x = d_no$orf_biotype_single, y = d_no$log_ppm,
          name = "No match", side = "negative",
          fillcolor = "rgba(189,195,199,0.55)", line = list(color = "#95A5A6", width = 1),
          box = list(visible = TRUE), points = FALSE,
          hovertemplate = "%{x} - No match<br>log₁₀PPM: %{y:.2f}<extra></extra>"
        )
      matched_bios <- intersect(bios, unique(d_yes$orf_biotype_single))
      for (bio in matched_bios) {
        d_bio <- filter(d_yes, orf_biotype_single == bio)
        if (nrow(d_bio) < 3) next
        col <- coalesce(BIOTYPE_COLORS[bio], "#28646E")
        p <- p %>% add_trace(
          type = "violin", x = d_bio$orf_biotype_single, y = d_bio$log_ppm,
          name = bio, legendgroup = "matched",
          legendgrouptitle = list(text = "Peptide match"),
          showlegend = TRUE, side = "positive",
          fillcolor = paste0(col, "99"), line = list(color = col, width = 1.5),
          box = list(visible = TRUE), points = FALSE,
          hovertemplate = "%{x} - matched<br>log₁₀PPM: %{y:.2f}<extra></extra>"
        )
      }
    }

    subtitle_text <- if (has_matches)
      "Gray (left half): unmatched ORFs of same biotype"
    else ""

    p %>% layout(
      title = list(
        text      = if (nchar(subtitle_text) > 0)
                      paste0("<span style='font-size:10px;color:#8A9CAA;'>", subtitle_text, "</span>")
                    else "",
        x         = 0,
        xanchor   = "left",
        font      = list(size = 10),
        pad       = list(t = 2)
      ),
      xaxis = list(title = "", tickangle = -30, automargin = TRUE),
      yaxis = list(title = "log₁₀(Median PPM + 0.01)", showgrid = TRUE, gridcolor = "#EEF2F7"),
      showlegend = FALSE,
      paper_bgcolor = "white", plot_bgcolor = "white",
      margin = list(l = 10, r = 10, t = 28, b = 40),
      font   = list(family = "Inter", size = 12)
    ) %>% config(displayModeBar = FALSE)
  })

  # ── Prioritisation stats ─────────────────────────────────────────────────────
  output$stat_prio_total <- renderText({
    req(gene_prioritised_data()); formatC(nrow(gene_prioritised_data()), big.mark = ",")
  })
  output$stat_prio_pep <- renderText({
    req(gene_prioritised_data()); formatC(sum(gene_prioritised_data()$n_peptides), big.mark = ",")
  })
  output$stat_prio_top <- renderText({
    req(gene_prioritised_data())
    sprintf("%.1f / 100", max(gene_prioritised_data()$priority_score, na.rm = TRUE))
  })

  # ── Priority table (gene-centric) ────────────────────────────────────────────
  prio_table_df <- reactive({
    req(gene_prioritised_data())
    df <- gene_prioritised_data() %>%
      mutate(
        spec_category = case_when(
          GTEX_tumor_only %in% TRUE     ~ "Tumor-only",
          GTEX_tumor_enriched %in% TRUE ~ "Tumor-enriched",
          TRUE                          ~ "Non-specific"
        ),
        score_html   = vapply(priority_score, score_bar_html, character(1)),
        spec_html    = mapply(spec_badge_html, GTEX_tumor_only, GTEX_tumor_enriched),
        biotype_html = vapply(orf_biotype_single, biotype_badge_html, character(1)),
        transl_html  = vapply(target_translation_pct_samples,
                              function(x) pct_bar_html(x, "#28646E"), character(1)),
        expr_html    = vapply(target_expression_pct_samples,
                              function(x) pct_bar_html(x, "#7EB8BF"), character(1))
      )
    sel_spec    <- input$prio_spec_filter %||% c("Tumor-only", "Tumor-enriched", "Non-specific")
    gtex_cutoff <- if (isTRUE(input$gtex_max_tpm_filter)) 1 else Inf
    df <- filter(df,
                 spec_category %in% sel_spec,
                 replace_na(GTEX_max_median_TPM, 0) <= gtex_cutoff)
    # DT column layout (0-based, after dropping .orf_id + orf_ids before DT):
    # Sel(0) Gene(1) ORF-biotype(2) Peptides(3) ORF-id(4) Location(5)
    # Specificity(6) Score(7) Transl.%(8) Transl.PPM(9) Expr.%(10) Expr.TPM(11)
    # GTEx(12) TCGA T%(13) TCGA T TPM(14) TCGA N%(15) TCGA N TPM(16)
    # RC prim%(17) RC prim PPM(18) RC CL%(19) RC CL PPM(20)
    # .biotype_sort(21) .spec_sort(22) .score_sort(23) .transl_sort(24) .expr_sort(25) .child_rows(26)
    df %>% transmute(
      Sel            = sprintf('<input type="checkbox" class="titan-row-checkbox" data-rowid="%d">', .row_id),
      Gene           = {
        link <- sprintf('<span class="titan-gene-link fw-semibold fst-italic" data-rowid="%d">%s</span>',
                        .row_id, gene_name)
        expand <- ifelse(n_orfs > 1L,
                         ' <span class="titan-expand-btn titan-orf-expand">+</span>', "")
        paste0(link, expand)
      },
      `ORF-biotype`  = biotype_html,
      Peptides       = vapply(matched_peptides, make_peptide_cell, character(1)),
      `ORF-id`       = sprintf('<span class="font-monospace" style="font-size:10px;word-break:break-all">%s</span>',
                               orf_id),
      Location       = sprintf('%s:%s&ndash;%s %s %s',
                               chr,
                               formatC(orf_start, format = "d", big.mark = ","),
                               formatC(orf_end,   format = "d", big.mark = ","),
                               strand, start_codon),
      Specificity    = spec_html,
      Score          = score_html,
      `Transl. %`    = transl_html,
      `Transl. PPM`  = round(target_translation_median_PPM,  2),
      `Expr. %`      = expr_html,
      `Expr. TPM`    = round(target_expression_median_TPM,   2),
      `GTEx max TPM` = round(GTEX_max_median_TPM,            3),
      `TCGA T%`      = round(TCGA_tumor_pct_samples,         1),
      `TCGA T TPM`   = round(TCGA_tumor_median_TPM,          2),
      `TCGA N%`      = round(TCGA_normal_pct_samples,        1),
      `TCGA N TPM`   = round(TCGA_normal_median_TPM,         2),
      `RC prim %`    = round(ribocrypt_primary_pct_samples,  1),
      `RC prim PPM`  = round(ribocrypt_primary_median_PPM,   2),
      `RC CL %`      = round(`ribocrypt_cell-line_pct_samples`, 1),
      `RC CL PPM`    = round(`ribocrypt_cell-line_median_PPM`,  2),
      .biotype_sort  = orf_biotype_single,
      .spec_sort     = spec_category,
      .score_sort    = round(priority_score, 2),
      .transl_sort   = round(target_translation_pct_samples, 1),
      .expr_sort     = round(target_expression_pct_samples,  1),
      .child_rows    = .child_html,
      .orf_id        = orf_id,        # kept for ORF Detail dropdown population
      orf_ids        = orf_ids        # kept for ORF Detail dropdown population
    )
  })

  output$tbl_priority <- renderDT({
    df <- prio_table_df() %>% select(-`.orf_id`, -`orf_ids`)
    # Col layout (0-based): Sel(0) Gene(1) ORF-biotype(2) Peptides(3) ORF-id(4)
    # Location(5) Specificity(6) Score(7) Transl.%(8) PPM(9) Expr.%(10) TPM(11)
    # GTEx(12) TCGAT%(13) TCGATPM(14) TCGAN%(15) TCANPM(16) RCprim%(17) RCprimPPM(18)
    # RCCL%(19) RCCLPPM(20) .biotype_sort(21) .spec_sort(22) .score_sort(23)
    # .transl_sort(24) .expr_sort(25) .child_rows(26)
    n_vis <- ncol(df) - 6L   # 6 hidden cols: 5 sort + .child_rows
    datatable(
      df,
      escape    = FALSE,
      rownames  = FALSE,
      selection = "none",
      class     = "compact hover",
      callback  = JS("
        // Off -> re-add with namespace to avoid stacking on re-render.
        $('#tbl_priority').off('.titanprio');
        $(document).off('.titanpriohdr');

        // Track SELECTED row IDs (starts empty = all unchecked by default).
        window.titanPrioSel = new Set();

        // -- Gene link: flash + open candidate detail panel -----------------
        $('#tbl_priority').on('click.titanprio', '.titan-gene-link', function(e) {
          e.stopPropagation();
          var $el = $(this);
          $el.addClass('titan-gene-active');
          setTimeout(function() { $el.removeClass('titan-gene-active'); }, 450);
          var rid = parseInt($el.attr('data-rowid'));
          if (!isNaN(rid))
            Shiny.setInputValue('prio_gene_click', {rowid: rid, nonce: Math.random()}, {priority: 'event'});
        });

        // -- ORF expand: inject child rows as <tr> siblings -----------------
        $('#tbl_priority').on('click.titanprio', '.titan-orf-expand', function(e) {
          e.stopPropagation();
          var $btn = $(this);
          var $tr  = $btn.closest('tr');
          var ch   = $.data($tr[0], 'child');
          if (!ch) return;
          if ($tr.hasClass('titan-orf-expanded')) {
            $tr.nextUntil(':not(.titan-child-row)').remove();
            $tr.removeClass('titan-orf-expanded');
            $btn.text('+');
          } else {
            $('<tbody>' + ch + '</tbody>').children().insertAfter($tr);
            $tr.addClass('titan-orf-expanded');
            $btn.text('\\u2212');
          }
        });

        // -- Peptide cell: and N more... / less toggle ----------------------
        $('#tbl_priority').on('click.titanprio', '.titan-pep-more, .titan-pep-less', function(e) {
          e.stopPropagation();
          var $cell  = $(this).closest('td');
          var $more  = $cell.find('.titan-pep-more');
          var $less  = $cell.find('.titan-pep-less');
          var $extra = $cell.find('.titan-pep-extra');
          var open   = $extra.css('display') !== 'none';
          $extra.css('display', open ? 'none' : '');
          $more.css('display', open ? '' : 'none');
          $less.css('display', open ? 'none' : '');
        });

        // -- Row checkbox: track selected rows ------------------------------
        $('#tbl_priority').on('change.titanprio', '.titan-row-checkbox', function() {
          var rid = parseInt($(this).data('rowid'));
          this.checked ? window.titanPrioSel.add(rid) : window.titanPrioSel.delete(rid);
          titanSyncHeader();
          Shiny.setInputValue('prio_selected_rowids', Array.from(window.titanPrioSel), {priority: 'event'});
        });

        // -- Header checkbox: select/deselect current page ------------------
        $(document).on('change.titanpriohdr', '#titan-hdr-cb', function() {
          var ok = this.checked;
          $('#tbl_priority tbody .titan-row-checkbox').each(function() {
            var rid = parseInt($(this).data('rowid'));
            $(this).prop('checked', ok);
            ok ? window.titanPrioSel.add(rid) : window.titanPrioSel.delete(rid);
          });
          this.indeterminate = false;
          Shiny.setInputValue('prio_selected_rowids', Array.from(window.titanPrioSel), {priority: 'event'});
        });

        function titanSyncHeader() {
          var cbs = $('#tbl_priority tbody .titan-row-checkbox');
          var n = cbs.length, nc = cbs.filter(':checked').length;
          var h = document.getElementById('titan-hdr-cb');
          if (!h) return;
          h.checked = (nc === n && n > 0); h.indeterminate = (nc > 0 && nc < n);
        }
        window.titanSyncHeader = titanSyncHeader;
      "),
      options   = list(
        pageLength = 20,
        dom        = "Bfrtip",
        scrollX        = TRUE,
        scrollY        = "1px",
        scrollCollapse = TRUE,
        createdRow = JS("function(row, data, index) {
          var ch = data[data.length - 1];
          if (ch) $.data(row, 'child', ch);
        }"),
        headerCallback = JS("function(thead) {
          $(thead).find('th:first').html('<input type=\"checkbox\" id=\"titan-hdr-cb\" style=\"cursor:pointer\" title=\"Select/deselect current page\">');
        }"),
        drawCallback = JS("function() {
          var sel = window.titanPrioSel || new Set();
          $('#tbl_priority tbody .titan-row-checkbox').each(function() {
            $(this).prop('checked', sel.has(parseInt($(this).data('rowid'))));
          });
          if (window.titanSyncHeader) window.titanSyncHeader();
        }"),
        columnDefs = list(
          list(className = "dt-center titan-sel-col", targets = 0L),
          list(className = "titan-pep-cell",          targets = 3L),
          list(visible   = FALSE, targets = seq(n_vis, n_vis + 5L)),
          list(orderData = n_vis,      targets = 2L),   # ORF-biotype
          list(orderData = n_vis + 1L, targets = 6L),   # Specificity
          list(orderData = n_vis + 2L, targets = 7L),   # Score
          list(orderData = n_vis + 3L, targets = 8L),   # Transl. %
          list(orderData = n_vis + 4L, targets = 10L),  # Expr. %
          list(orderable = FALSE, targets = 0L)         # Sel not sortable
        ),
        lengthMenu = list(c(10, 20, 50), c("10", "20", "50"))
      )
    ) %>%
      formatStyle("GTEx max TPM",
                  color      = styleInterval(1, c("inherit", "#C0392B")),
                  fontWeight = styleInterval(1, c("normal", "600")))
  }, server = TRUE)

  prio_row_id          <- reactiveVal(NULL)
  prio_selected_rowids <- reactiveVal(integer(0))

  # Gene name link clicked → open candidate detail panel
  observeEvent(input$prio_gene_click, {
    rid <- as.integer(input$prio_gene_click$rowid)
    prio_row_id(rid)
    # Also pre-populate the ORF Detail tab with the best ORF for this gene
    gdata <- isolate(gene_prioritised_data())
    row   <- filter(gdata, .row_id == rid)
    if (nrow(row) > 0)
      updateSelectizeInput(session, "detail_orf_id", selected = row$orf_id[1])
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(input$prio_selected_rowids, {
    prio_selected_rowids(as.integer(input$prio_selected_rowids %||% integer(0)))
  }, ignoreNULL = FALSE)

  selected_prio_row <- reactive({
    rid <- prio_row_id()
    req(!is.null(rid))
    filter(gene_prioritised_data(), .row_id == rid)
  })

  output$priority_detail_panel <- renderUI({
    req(nrow(selected_prio_row()) > 0)
    row <- selected_prio_row()
    w   <- current_weights()

    metrics <- list(
      list("% samples (expr.)",     fmt1(row$target_expression_pct_samples)),
      list("Median TPM",            fmt2(row$target_expression_median_TPM)),
      list("% samples (transl.)",   fmt1(row$target_translation_pct_samples)),
      list("Median PPM",            fmt2(row$target_translation_median_PPM)),
      list("GTEx max TPM",          fmt2(row$GTEX_max_median_TPM)),
      list("TCGA tumor %",          fmt1(row$TCGA_tumor_pct_samples)),
      list("TCGA normal %",         fmt1(row$TCGA_normal_pct_samples)),
      list("RC primary %",          fmt1(row$ribocrypt_primary_pct_samples)),
      list("RC cell-line %",        fmt1(row$`ribocrypt_cell-line_pct_samples`))
    )
    tile_ui <- lapply(metrics, function(m) {
      tags$div(class = "prio-metric-tile",
        tags$p(class = "prio-metric-label", m[[1]]),
        tags$p(class = "prio-metric-value", m[[2]])
      )
    })

    dim_bars <- lapply(WEIGHT_META, function(m) {
      raw <- as.numeric(row[[dim_col(m$id)]])
      wv  <- as.numeric(w[[m$id]])
      pct     <- if (wv == 0) 0 else min(100, abs(raw) / abs(wv) * 100)
      bar_col <- if (wv == 0) "#CCCCCC" else if (wv > 0) "#00A555" else "#C0392B"
      tags$div(class = "prio-dim-row",
        tags$div(class = "prio-dim-header",
          tags$span(class = "prio-dim-label", m$label),
          tags$span(class = "prio-dim-contribution", sprintf("%.2f", raw))
        ),
        tags$div(class = "prio-dim-track",
          tags$div(class = "prio-dim-fill",
            style = sprintf("width:%.0f%%;background:%s", pct, bar_col)))
      )
    })

    n_orfs_label <- paste0(
      row$n_orfs, " ORF", if (isTRUE(row$n_orfs == 1L)) "" else "s", ": ",
      row$orf_ids
    )
    n_pep_label <- paste0(
      row$n_peptides, " peptide", if (isTRUE(row$n_peptides == 1L)) "" else "s", ": ",
      row$matched_peptides
    )

    div(class = "prio-detail-outer",
      # ── Header ──────────────────────────────────────────────────────────────
      div(class = "d-flex justify-content-between align-items-start mb-2",
        div(class = "flex-grow-1 me-3",
          tags$p(class = "prio-detail-overline", "Candidate detail"),
          tags$p(class = "prio-detail-gene", row$gene_name,
                 HTML(paste0(' ', biotype_badge_html(row$orf_biotype_single)))),
          div(class = "d-flex align-items-center justify-content-between gap-2 mb-0",
            tags$p(class = "prio-detail-orfs text-muted small mb-0",
                   n_orfs_label),
            HTML(spec_badge_html(row$GTEX_tumor_only, row$GTEX_tumor_enriched))
          ),
          tags$p(class = "prio-detail-pep text-muted small mb-0", n_pep_label)
        ),
        div(class = "d-flex flex-column align-items-end gap-1 flex-shrink-0",
          actionButton("close_prio_detail", label = NULL, icon = icon("xmark"),
                       class = "btn-sm btn-outline-primary",
                       title = "Close detail panel"),
          div(class = "text-end mt-1",
            tags$span(class = "prio-detail-score-badge",
                      sprintf("%.1f", row$priority_score)),
            tags$span(class = "prio-detail-score-label", " / 100")
          )
        )
      ),
      # ── Summary metric tiles (9 tiles, one row) ───────────────────────────
      div(class = "prio-metric-grid", tile_ui),
      # ── Expression card ──────────────────────────────────────────────────
      card(
        class = "mt-2",
        card_header(
          div(class = "d-flex justify-content-between align-items-center",
            tags$span("Expression"),
            radioButtons("prio_expr_scale", NULL,
                         choices = c("log(TPM+1)" = "log", "TPM" = "raw"),
                         selected = "log", inline = TRUE)
          )
        ),
        card_body(class = "p-1", plotlyOutput("plot_prio_expr", height = "260px"))
      ),
      # ── Translation card ─────────────────────────────────────────────────
      card(
        class = "mt-2",
        card_header(
          div(class = "d-flex justify-content-between align-items-center",
            tags$span("Translation"),
            radioButtons("prio_transl_scale", NULL,
                         choices = c("log(PPM+1)" = "log", "PPM" = "raw"),
                         selected = "log", inline = TRUE)
          )
        ),
        card_body(class = "p-1", plotlyOutput("plot_prio_transl", height = "200px"))
      ),
      # ── Score profile + dimension contributions ───────────────────────────
      layout_columns(
        col_widths = c(5, 7),
        class = "mt-2",
        card(
          card_header("Score profile"),
          card_body(class = "p-1", plotlyOutput("plot_radar", height = "220px"))
        ),
        card(
          card_header("Dimension contributions"),
          card_body(dim_bars)
        )
      )
    )
  })

  observeEvent(input$close_prio_detail, {
    prio_row_id(NULL)
  })

  output$plot_radar <- renderPlotly({
    req(nrow(selected_prio_row()) > 0)
    row <- selected_prio_row()
    w   <- current_weights()

    radar_labels <- sapply(WEIGHT_META, `[[`, "radar")
    vals <- vapply(WEIGHT_META, function(m) {
      raw <- as.numeric(row[[dim_col(m$id)]])
      wv  <- as.numeric(w[[m$id]])
      if (wv == 0) return(0)
      min(1, max(0, raw / wv))  # recovers signal (∈[0,1]) from contribution = signal × weight
    }, numeric(1))

    r_vals <- c(vals, vals[1])
    t_vals <- c(radar_labels, radar_labels[1])

    plot_ly(
      type = "scatterpolar", r = r_vals, theta = t_vals,
      fill = "toself", fillcolor = "rgba(0,165,85,0.18)",
      line = list(color = "#00A555", width = 1.5),
      mode = "lines+markers",
      marker = list(color = "#00A555", size = 4),
      hovertemplate = "%{theta}: %{r:.2f}<extra></extra>"
    ) %>%
      layout(
        polar = list(
          radialaxis  = list(visible = TRUE, range = c(0, 1),
                             showticklabels = FALSE,
                             gridcolor = "#DAE9EC", linecolor = "#DAE9EC"),
          angularaxis = list(tickfont = list(size = 9, color = "#555"),
                             gridcolor = "#DAE9EC")
        ),
        showlegend = FALSE,
        paper_bgcolor = "white", plot_bgcolor = "white",
        margin = list(l = 45, r = 45, t = 45, b = 45)
      ) %>%
      config(displayModeBar = FALSE)
  })

  # ── Expression card: Target / GTEx / TCGA ─────────────────────────────────
  output$plot_prio_expr <- renderPlotly({
    req(nrow(selected_prio_row()) > 0)
    row       <- selected_prio_row()
    log_scale <- isTRUE((input$prio_expr_scale %||% "log") == "log")
    y_label   <- if (log_scale) "log(TPM+1)" else "TPM"
    y_ref     <- if (log_scale) log(2) else 1   # threshold at TPM = 1
    gid       <- row$gene_id_clean

    apply_scale <- function(x) if (log_scale) log(pmax(as.numeric(x), 0) + 1) else pmax(as.numeric(x), 0)

    make_box_args <- function(fill_hex) list(
      type = "box", boxpoints = "all", jitter = 0.35, pointpos = 0,
      fillcolor    = paste0(fill_hex, "1F"),
      line         = list(color = fill_hex, width = 1.5),
      whiskerwidth = 0.5, showlegend = FALSE, hoveron = "boxes",
      marker       = list(symbol = "circle-open", size = 5, opacity = 0.5, color = fill_hex)
    )
    no_data_plot <- function(x_title, msg) {
      plot_ly(type = "scatter", mode = "markers", x = 0, y = 0,
              marker = list(opacity = 0)) %>%
        layout(xaxis = list(title = x_title, showticklabels = FALSE),
               yaxis = list(title = ""),
               annotations = list(list(
                 text = msg, x = 0.5, y = 0.5,
                 xref = "paper", yref = "paper", showarrow = FALSE,
                 font = list(size = 10, color = "#6C757D"))))
    }

    # --- Plot 1: Target tumor (per-sample) ---
    rna_mat  <- tryCatch(rna_tpm_rv(), error = function(e) NULL)
    rna_meta <- tryCatch(rna_meta_rv(), error = function(e) NULL)

    if (!is.null(rna_mat) && isTRUE(gid %in% rownames(rna_mat))) {
      tpm_raw <- as.numeric(rna_mat[gid, ])
      grp_col <- if (!is.null(rna_meta)) {
        label_col <- if ("condition"   %in% colnames(rna_meta)) rna_meta$condition
                     else if ("tissue_type" %in% colnames(rna_meta)) rna_meta$tissue_type
                     else rep("Tumor", nrow(rna_meta))
        label_col[match(colnames(rna_mat), rna_meta$sample_id)]
      } else rep("Tumor", length(tpm_raw))
      grp_col[is.na(grp_col)] <- "Tumor"
      df_t <- data.frame(g = grp_col, y = apply_scale(tpm_raw))
      p1 <- do.call(plot_ly, c(list(df_t, x = ~g, y = ~y), make_box_args("#28646E"))) %>%
        layout(xaxis = list(title = list(text = "", standoff = 4), automargin = TRUE),
               yaxis = list(title = y_label))
    } else {
      p1 <- no_data_plot("Target tumor", "No RNA-seq matrix")
      p1 <- p1 %>% layout(yaxis = list(title = y_label))
    }

    # --- Plot 2: GTEx – per-sample boxplots by tissue ---
    gtex_mat  <- tryCatch(gtex_tpm_rv(),  error = function(e) NULL)
    gtex_meta <- tryCatch(gtex_meta_rv(), error = function(e) NULL)

    if (!is.null(gtex_mat) && isTRUE(gid %in% rownames(gtex_mat))) {
      gtex_raw   <- as.numeric(gtex_mat[gid, ])
      tissue_raw <- gtex_meta$tissue_type[match(colnames(gtex_mat), gtex_meta$sample_id)]
      tissue_raw[is.na(tissue_raw)] <- "Unknown"
      tissue     <- gsub("_", " ", tissue_raw)
      gtex_y     <- apply_scale(gtex_raw)
      med_by_tis <- tapply(gtex_y, tissue, median, na.rm = TRUE)
      sorted_tis <- names(sort(med_by_tis, decreasing = TRUE))
      # Per-tissue color: exact sub-tissue match first, then group-level prefix fallback
      tis_col_map <- setNames(
        sapply(unique(tissue_raw), function(t) {
          idx <- match(t, gtex_colors_subtissue$Tissue)
          if (!is.na(idx)) return(gtex_colors_subtissue$ColorHex[idx])
          gtex_keys <- names(gtex_colors)
          hits <- gtex_keys[startsWith(t, gtex_keys)]
          if (length(hits)) gtex_colors[[hits[which.max(nchar(hits))]]] else "#BBBBBB"
        }),
        gsub("_", " ", unique(tissue_raw))
      )
      df_g <- data.frame(g = factor(tissue, levels = sorted_tis), y = gtex_y)
      p2 <- plot_ly()
      for (tis in sorted_tis) {
        d <- df_g[as.character(df_g$g) == tis, , drop = FALSE]
        if (nrow(d) == 0L) next
        col <- tis_col_map[[tis]] %||% "#BBBBBB"
        p2 <- do.call(add_trace, c(list(p2, data = d, x = ~g, y = ~y), make_box_args(col)))
      }
      p2 <- p2 %>% layout(
        xaxis = list(title = list(text = "GTEx", standoff = 4),
                     tickangle = -45, automargin = TRUE,
                     categoryorder = "array", categoryarray = sorted_tis),
        yaxis = list(title = "")
      )
    } else {
      p2 <- no_data_plot("GTEx", "GTEx matrix not in app data\n(re-run prepare script)")
    }

    # --- Plot 3: TCGA – tumor (darkened study color) and normal (base study color) traces ---
    tcga_mat  <- tryCatch(tcga_tpm_rv(),  error = function(e) NULL)
    tcga_meta <- tryCatch(tcga_meta_rv(), error = function(e) NULL)

    if (!is.null(tcga_mat) && isTRUE(gid %in% rownames(tcga_mat))) {
      tcga_raw    <- as.numeric(tcga_mat[gid, ])
      grp         <- tcga_meta$group[match(colnames(tcga_mat), tcga_meta$sample_id)]
      grp[is.na(grp)] <- "Unknown"
      tcga_y      <- apply_scale(tcga_raw)
      unique_grps <- unique(grp)

      # Sort studies by max(median_tumor, median_normal), same statistic as GTEx sort.
      # Pairs stay intact (Tumor before Normal within each study).
      grp_medians <- tapply(tcga_y, grp, median, na.rm = TRUE)
      ct_codes    <- unique(sub(" .*", "", unique_grps))
      pair_key    <- vapply(ct_codes, function(ct) {
        meds <- grp_medians[sub(" .*", "", names(grp_medians)) == ct]
        if (length(meds) == 0L) -Inf else max(meds, na.rm = TRUE)
      }, numeric(1))
      ct_sorted   <- ct_codes[order(pair_key, decreasing = TRUE)]
      grp_order   <- unlist(lapply(ct_sorted, function(ct) {
        grps <- unique_grps[sub(" .*", "", unique_grps) == ct]
        grps[order(match(sub(".* ", "", grps), c("Tumor", "Normal")))]
      }), use.names = FALSE)

      grp_factor <- factor(grp, levels = grp_order)
      df_c <- data.frame(g = grp_factor, y = tcga_y)
      p3 <- plot_ly()
      for (lvl in grp_order) {
        d <- df_c[as.character(df_c$g) == lvl, , drop = FALSE]
        if (nrow(d) == 0) next
        ct_code <- sub(" .*", "", lvl)
        col <- if (grepl("Normal", lvl, ignore.case = TRUE))
          tcga_colors_normal[[ct_code]] %||% "#87C8D4"
        else
          tcga_colors_tumor[[ct_code]] %||% "#3B95A5"
        p3 <- do.call(add_trace, c(list(p3, data = d, x = ~g, y = ~y), make_box_args(col)))
      }
      p3 <- p3 %>% layout(
        xaxis = list(title = list(text = "TCGA", standoff = 4),
                     tickangle = -45, automargin = TRUE,
                     categoryorder = "array", categoryarray = grp_order),
        yaxis = list(title = "")
      )
    } else {
      p3 <- no_data_plot("TCGA", "TCGA matrix not in app data\n(re-run prepare script)")
    }

    out <- subplot(p1, p2, p3, nrows = 1, shareY = TRUE, titleX = TRUE,
                   widths = c(0.12, 0.53, 0.35)) %>%
      layout(
        paper_bgcolor = "white", plot_bgcolor = "white",
        margin = list(l = 5, r = 5, t = 10, b = 5),
        font   = list(family = "Inter", size = 11)
      ) %>%
      config(
        toImageButtonOptions = list(format = "svg", filename = "expression"),
        modeBarButtons = list(list("toImage"))
      )
    out$x$layout$shapes <- list(list(
      type = "line",
      xref = "paper", x0 = 0, x1 = 1,
      yref = "y",     y0 = y_ref, y1 = y_ref,
      line = list(color = "#BBBBBB", width = 1.5, dash = "dash")
    ))
    out
  })

  # ── Translation card: Target / Ribocrypt ──────────────────────────────────
  output$plot_prio_transl <- renderPlotly({
    req(nrow(selected_prio_row()) > 0)
    row       <- selected_prio_row()
    log_scale <- isTRUE((input$prio_transl_scale %||% "log") == "log")
    y_label   <- if (log_scale) "log(PPM+1)" else "PPM"
    y_ref     <- if (log_scale) log(2) else 1   # threshold at PPM = 1
    oid       <- row$orf_id

    apply_scale <- function(x) if (log_scale) log(pmax(as.numeric(x), 0) + 1) else pmax(as.numeric(x), 0)

    make_box_args <- function(fill_hex) list(
      type = "box", boxpoints = "all", jitter = 0.35, pointpos = 0,
      fillcolor    = paste0(fill_hex, "1F"),
      line         = list(color = fill_hex, width = 1.5),
      whiskerwidth = 0.5, showlegend = FALSE, hoveron = "boxes",
      marker       = list(symbol = "circle-open", size = 5, opacity = 0.5, color = fill_hex)
    )
    no_data_plot <- function(x_title, msg) {
      plot_ly(type = "scatter", mode = "markers", x = 0, y = 0,
              marker = list(opacity = 0)) %>%
        layout(xaxis = list(title = x_title, showticklabels = FALSE),
               yaxis = list(title = ""),
               annotations = list(list(
                 text = msg, x = 0.5, y = 0.5,
                 xref = "paper", yref = "paper", showarrow = FALSE,
                 font = list(size = 10, color = "#6C757D"))))
    }

    # --- Plot 1: Target tumor (per-sample ribo-seq) ---
    ribo_m  <- tryCatch(ribo_ppm_rv(),  error = function(e) NULL)
    ribo_sm <- tryCatch(ribo_meta_rv(), error = function(e) NULL)

    if (!is.null(ribo_m) && isTRUE(oid %in% rownames(ribo_m))) {
      ppm_raw <- as.numeric(ribo_m[oid, ])
      cond <- if (!is.null(ribo_sm) && "condition" %in% colnames(ribo_sm))
        ribo_sm$condition[match(colnames(ribo_m), ribo_sm$sample_id)]
      else rep("Tumor", length(ppm_raw))
      cond[is.na(cond)] <- "Tumor"
      df_r <- data.frame(g = cond, y = apply_scale(ppm_raw))
      p1 <- do.call(plot_ly, c(list(df_r, x = ~g, y = ~y), make_box_args("#28646E"))) %>%
        layout(xaxis = list(title = list(text = "", standoff = 4), automargin = TRUE),
               yaxis = list(title = y_label))
    } else {
      p1 <- no_data_plot("Target tumor", "No ribo-seq data for this ORF")
      p1 <- p1 %>% layout(yaxis = list(title = y_label))
    }

    # --- Plots 2 & 3: Ribocrypt Primary / Cell-line – per-sample scatter (GTEx palette) ---
    rc_mat  <- tryCatch(ribocrypt_mat_rv(),   error = function(e) NULL)
    rc_meta <- tryCatch(ribocrypt_smeta_rv(), error = function(e) NULL)

    make_rc_panel <- function(s_ids, s_y, x_title, label_fn = identity) {
      if (length(s_ids) == 0L)
        return(no_data_plot(x_title, paste("No", x_title, "data")))
      s_ord    <- order(s_y, decreasing = TRUE)
      s_ids    <- s_ids[s_ord]
      s_y      <- s_y[s_ord]
      s_labels <- label_fn(s_ids)
      s_cols   <- ifelse(s_ids %in% names(rc_color_map), rc_color_map[s_ids], "#AAAAAA")
      n        <- length(s_labels)
      # Lollipop stems: NA-separated segments from y = 0 to each point
      stem_x   <- as.vector(rbind(s_labels, s_labels, NA_character_))
      stem_y   <- as.vector(rbind(rep(0, n), s_y, NA_real_))
      plot_ly() %>%
        add_trace(
          x = stem_x, y = stem_y,
          type = "scatter", mode = "lines",
          line = list(color = "#CCCCCC", width = 1),
          showlegend = FALSE, hoverinfo = "skip"
        ) %>%
        add_trace(
          x = s_labels, y = s_y,
          type = "scatter", mode = "markers",
          marker = list(symbol = "circle", size = 7, color = s_cols,
                        line   = list(color = s_cols, width = 0)),
          showlegend = FALSE,
          hovertext = sprintf("%s<br>%s: %.3f", s_labels, y_label, s_y),
          hoverinfo = "text"
        ) %>%
        layout(
          xaxis = list(title    = list(text = x_title, standoff = 4),
                       tickangle = -45, automargin = TRUE,
                       ticks    = "",
                       categoryorder = "array", categoryarray = s_labels),
          yaxis = list(title = "")
        )
    }

    if (!is.null(rc_mat) && isTRUE(oid %in% rownames(rc_mat))) {
      rc_raw <- as.numeric(rc_mat[oid, ])
      sids   <- colnames(rc_mat)
      grp    <- rc_meta$group[match(sids, rc_meta$sample_id)]
      grp[is.na(grp)] <- "Unknown"
      rc_y   <- apply_scale(rc_raw)

      p_prim <- make_rc_panel(sids[grp == "Primary"],   rc_y[grp == "Primary"],   "RC Primary",
                              label_fn = function(x) sub("^primary_", "", x))
      p_cl   <- make_rc_panel(sids[grp == "Cell-line"], rc_y[grp == "Cell-line"], "RC Cell-line")
    } else {
      p_prim <- no_data_plot("RC Primary",   "Ribocrypt matrix not in app data\n(re-run prepare script)")
      p_cl   <- no_data_plot("RC Cell-line", "Ribocrypt matrix not in app data\n(re-run prepare script)")
    }

    out <- subplot(p1, p_prim, p_cl, nrows = 1, shareY = TRUE, titleX = TRUE,
                   widths = c(0.10, 0.20, 0.70)) %>%
      layout(
        paper_bgcolor = "white", plot_bgcolor = "white",
        margin = list(l = 5, r = 5, t = 10, b = 5),
        font   = list(family = "Inter", size = 11)
      ) %>%
      config(
        toImageButtonOptions = list(format = "svg", filename = "translation"),
        modeBarButtons = list(list("toImage"))
      )
    out$x$layout$shapes <- list(list(
      type = "line",
      xref = "paper", x0 = 0, x1 = 1,
      yref = "y",     y0 = y_ref, y1 = y_ref,
      line = list(color = "#BBBBBB", width = 1.5, dash = "dash")
    ))
    out
  })

  output$dl_priority <- downloadHandler(
    filename = function() paste0("titan_priority_", format(Sys.time(), "%Y-%m-%d_%H%M"), ".csv"),
    content  = function(file) {
      req(prioritised_data())
      df <- prioritised_data() %>%
        select(.row_id, n_peptides, matched_peptides, gene_name, orf_biotype_single,
               chr, orf_start, orf_end, protein_length, start_codon,
               priority_score, starts_with("dim_"),
               target_expression_pct_samples, target_expression_median_TPM,
               target_translation_pct_samples, target_translation_median_PPM,
               GTEX_tumor_only, GTEX_tumor_enriched, GTEX_max_median_TPM,
               TCGA_tumor_pct_samples, TCGA_normal_pct_samples,
               ribocrypt_primary_pct_samples, `ribocrypt_cell-line_pct_samples`,
               orf_id)
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$dl_selected <- downloadHandler(
    filename = function() paste0("titan_selection_", format(Sys.time(), "%Y-%m-%d_%H%M"), ".csv"),
    content  = function(file) {
      req(gene_prioritised_data(), prioritised_data())
      sel       <- prio_selected_rowids()
      sel_genes <- gene_prioritised_data() %>%
        filter(.row_id %in% sel) %>%
        pull(gene_id)
      df <- prioritised_data() %>%
        filter(gene_id %in% sel_genes) %>%
        select(.row_id, n_peptides, matched_peptides, gene_name, orf_biotype_single,
               chr, orf_start, orf_end, protein_length, start_codon,
               priority_score, starts_with("dim_"),
               target_expression_pct_samples, target_expression_median_TPM,
               target_translation_pct_samples, target_translation_median_PPM,
               GTEX_tumor_only, GTEX_tumor_enriched, GTEX_max_median_TPM,
               TCGA_tumor_pct_samples, TCGA_normal_pct_samples,
               ribocrypt_primary_pct_samples, `ribocrypt_cell-line_pct_samples`,
               orf_id)
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$dl_params <- downloadHandler(
    filename = function() paste0("titan_scoring_params_", format(Sys.time(), "%Y-%m-%d_%H%M"), ".csv"),
    content  = function(file) {
      preset <- active_preset()
      w      <- current_weights()
      params <- data.frame(
        preset    = if (is.null(preset)) "Custom" else preset,
        dimension = sapply(WEIGHT_META, `[[`, "label"),
        weight    = sapply(WEIGHT_META, function(m) as.numeric(w[[m$id]])),
        signal    = sapply(WEIGHT_META, `[[`, "hint"),
        stringsAsFactors = FALSE
      )
      write.csv(params, file, row.names = FALSE)
    }
  )

  # ── ORF Detail ───────────────────────────────────────────────────────────────
  detail_orf <- reactive({
    req(input$detail_orf_id)
    oid      <- input$detail_orf_id
    from_tbl <- filter(orf_table_rv(), orf_id == oid)
    if (nrow(from_tbl) > 0L) return(from_tbl)
    # Gencode-only ORF: build one-row summary from matched_data
    md <- tryCatch(matched_data(), error = function(e) NULL)
    if (is.null(md)) return(from_tbl)
    md %>% filter(orf_id == oid) %>% slice(1L)
  })

  # ORFs that share an identical set of matched peptides with the selected ORF.
  detail_orf_siblings <- reactive({
    req(input$detail_orf_id)
    selected <- input$detail_orf_id
    md <- tryCatch(matched_data(), error = function(e) NULL)
    if (is.null(md) || nrow(md) == 0L) return(character(0))

    sel_key <- md %>%
      filter(orf_id == selected) %>%
      pull(matched_peptide) %>%
      { paste(sort(unique(.)), collapse = "|") }

    if (!nzchar(sel_key)) return(character(0))

    md %>%
      group_by(orf_id) %>%
      summarise(key = paste(sort(unique(matched_peptide)), collapse = "|"), .groups = "drop") %>%
      filter(key == sel_key, orf_id != selected) %>%
      pull(orf_id)
  })

  output$detail_orf_meta <- renderUI({
    req(nrow(detail_orf()) > 0)
    o        <- detail_orf()
    siblings <- detail_orf_siblings()
    tbl      <- orf_table_rv()

    main_block <- tags$div(
      tags$p(tags$b("Gene: "), o$gene_name, "  ",
             tags$span(class = "badge bg-secondary", o$gene_biotype)),
      tags$p(tags$b("Biotype: "), HTML(biotype_badge_html(o$orf_biotype_single))),
      tags$p(tags$b("Coordinates: "),
             paste0(o$chr, ":", format(o$orf_start, big.mark = ","), "–",
                    format(o$orf_end, big.mark = ","), " (", o$strand, ")")),
      tags$p(tags$b("Protein length: "), o$protein_length, " aa"),
      tags$p(tags$b("Start codon: "), o$start_codon),
      tags$small(class = "text-muted font-monospace",
                 style = "word-break:break-all;", o$orf_id)
    )

    # Gencode cross-match annotation (case a: in-house ORF also matched a Gencode entry)
    md <- tryCatch(matched_data(), error = function(e) NULL)
    gc_ids <- if (!is.null(md) && "gencode_match_ids" %in% colnames(md)) {
      vals <- unique(md$gencode_match_ids[md$orf_id == o$orf_id])
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals)) paste(vals, collapse = "; ") else ""
    } else ""
    gc_block <- if (nzchar(gc_ids))
      tags$p(tags$b("Gencode cross-match: "),
             tags$small(class = "text-muted font-monospace", gc_ids))
    else NULL

    sib_section <- NULL
    if (length(siblings) > 0L) {
      sib_tbl  <- filter(tbl, orf_id %in% siblings)
      sib_rows <- lapply(seq_len(nrow(sib_tbl)), function(i) {
        s <- sib_tbl[i, ]
        tags$div(class = "d-flex align-items-start gap-2 mb-1",
          HTML(biotype_badge_html(s$orf_biotype_single)),
          tags$small(class = "text-muted font-monospace",
                     paste0(s$protein_length, " aa · ",
                            s$chr, ":", format(s$orf_start, big.mark = ","),
                            "–", format(s$orf_end, big.mark = ","),
                            "  ", s$orf_id))
        )
      })
      sib_section <- tagList(
        tags$hr(class = "my-2"),
        tags$p(class = "text-muted small mb-1",
               paste0(length(siblings), " co-identified ORF",
                      if (length(siblings) > 1L) "s" else "",
                      " (identical peptide evidence):")),
        tagList(sib_rows)
      )
    }

    if (is.null(gc_block) && is.null(sib_section)) return(main_block)

    tagList(
      main_block,
      if (!is.null(gc_block)) tagList(tags$hr(class = "my-2"), gc_block) else NULL,
      sib_section
    )
  })

  output$detail_metric_badges <- renderUI({
    req(nrow(detail_orf()) > 0)
    o <- detail_orf()
    badges <- list()
    if (isTRUE(as.logical(o$GTEX_tumor_only)))
      badges[[length(badges) + 1]] <- pill_badge("GTEx tumor-only",           "success")
    if (isTRUE(as.logical(o$GTEX_tumor_enriched)))
      badges[[length(badges) + 1]] <- pill_badge("GTEx tumor-enriched",       "primary")
    if (isTRUE(as.logical(o$GTEX_DE_sig_in_all)))
      badges[[length(badges) + 1]] <- pill_badge("DE significant (all GTEx)", "info")
    if (length(badges) == 0)
      badges[[1]] <- tags$span(class = "text-muted small", "No tumor-specificity flags")
    div(badges)
  })

  output$detail_metric_table <- renderTable({
    req(nrow(detail_orf()) > 0)
    o <- detail_orf()
    rows <- list(
      c("Translation - % samples (PPM ≥ threshold)",   fmt1(o$target_translation_pct_samples)),
      c("Translation - median PPM",                     fmt2(o$target_translation_median_PPM)),
      c("Translation - max PPM",                        fmt2(o$target_translation_max_PPM)),
      c("Translation - median psites",                  fmt2(o$target_translation_median_psites)),
      c("Expression - % samples (TPM ≥ threshold)",     fmt1(o$target_expression_pct_samples)),
      c("Expression - median TPM",                      fmt2(o$target_expression_median_TPM)),
      c("GTEx - DE sig. in all tissues",                bool_fmt(o$GTEX_DE_sig_in_all)),
      c("GTEx - tumor-enriched",                       bool_fmt(o$GTEX_tumor_enriched)),
      c("GTEx - tumor-only",                           bool_fmt(o$GTEX_tumor_only)),
      c("GTEx - median TPM (all samples)",              fmt2(o$GTEX_median_TPM)),
      c("GTEx - max tissue median TPM",                 fmt2(o$GTEX_max_median_TPM)),
      c("TCGA - % tumor samples (TPM ≥ 1)",           fmt1(o$TCGA_tumor_pct_samples)),
      c("TCGA - median TPM (tumor)",                   fmt2(o$TCGA_tumor_median_TPM)),
      c("TCGA - % normal samples (TPM ≥ 1)",           fmt1(o$TCGA_normal_pct_samples)),
      c("TCGA - median TPM (normal)",                   fmt2(o$TCGA_normal_median_TPM)),
      c("Ribocrypt - % primary samples (PPM ≥ 1)",     fmt1(o$ribocrypt_primary_pct_samples)),
      c("Ribocrypt - median PPM (primary)",             fmt2(o$ribocrypt_primary_median_PPM)),
      c("Ribocrypt - % cell-line samples (PPM ≥ 1)",   fmt1(o$`ribocrypt_cell-line_pct_samples`)),
      c("Ribocrypt - median PPM (cell lines)",          fmt2(o$`ribocrypt_cell-line_median_PPM`))
    )
    df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    colnames(df) <- c("Metric", "Value")
    df
  }, striped = TRUE, hover = TRUE, bordered = FALSE, spacing = "xs",
     width = "100%", colnames = TRUE)

  output$plot_detail_ribo <- renderPlotly({
    req(nrow(detail_orf()) > 0)
    oid     <- input$detail_orf_id
    ribo_m  <- ribo_ppm_rv()
    ribo_sm <- ribo_meta_rv()
    if (!oid %in% rownames(ribo_m))
      return(plot_ly() %>% layout(title = "No ribo-seq data for this ORF",
                                  paper_bgcolor = "white") %>% config(displayModeBar = FALSE))
    ppm_vals <- ribo_m[oid, ]
    df <- data.frame(
      sample    = names(ppm_vals),
      ppm       = log10(as.numeric(ppm_vals)+0.1),
      condition = ribo_sm$condition[match(names(ppm_vals), ribo_sm$sample_id)]
    )
    plot_ly(df, x = ~sample, y = ~ppm, type = "bar",
            marker = list(color = "#317A87",
                          line = list(color =  "#28646E", 
                                      width = 1.5)),
            text = ~round(ppm, 2), textposition = "outside",
            hovertemplate = "%{x}<br>PPM: %{y:.2f}<extra></extra>") %>%
      layout(
        xaxis  = list(title = "", tickangle = -45, automargin = TRUE),
        yaxis  = list(title = "Psites per million (PPM)",
                      showgrid = TRUE, gridcolor = "#EEF2F7"),
        legend = list(title = list(text = "Condition"), orientation = "h", y = -0.3),
        paper_bgcolor = "white", plot_bgcolor = "white",
        margin = list(l = 10, r = 10, t = 10, b = 100),
        font   = list(family = "Inter", size = 11)
      ) %>% config(displayModeBar = FALSE)
  })

  output$detail_protein_seq_ui <- renderUI({
    req(nrow(detail_orf()) > 0)
    seq <- detail_orf()$protein_seq
    if (is.na(seq) || !nzchar(seq))
      return(tags$em(class = "text-muted", "No protein sequence available."))

    md  <- tryCatch(matched_data(), error = function(e) NULL)
    oid <- input$detail_orf_id
    peps <- if (!is.null(md)) unique(md %>% filter(orf_id == oid) %>% pull(matched_peptide)) else character(0)

    # Build per-peptide MS info for hover popovers (all columns from raw upload)
    pep_info <- list()
    if (length(peps) > 0L) {
      raw_ms <- tryCatch(ms_data(), error = function(e) NULL)
      pep_col <- if (!is.null(input$pep_col) && nzchar(input$pep_col)) {
        input$pep_col
      } else if (!is.null(raw_ms)) {
        m <- intersect(PEPTIDE_COL_CANDIDATES, colnames(raw_ms))
        if (length(m)) m[1L] else NULL
      } else NULL
      if (!is.null(raw_ms) && !is.null(pep_col) && pep_col %in% colnames(raw_ms)) {
        pep_info <- setNames(lapply(peps, function(p) {
          rows <- raw_ms[raw_ms[[pep_col]] == p, setdiff(colnames(raw_ms), pep_col), drop = FALSE]
          rows[!duplicated(rows), , drop = FALSE]
        }), peps)
      }
    }

    render_protein_seq_html(seq, peps, pep_info)
  })

  # ── About ─────────────────────────────────────────────────────────────────────
  output$about_data_info <- renderUI({
    lines <- list()
    d <- app_data_rv()
    if (!is.null(d) && !is.null(d$prepared_on))
      lines[[1]] <- tags$p(class = "text-muted small",
                           paste0("Data prepared: ", format(d$prepared_on, "%Y-%m-%d %H:%M")))
    if (!is.null(d) && is.null(tryCatch(rna_tpm_rv(), error = function(e) NULL)))
      lines[[2]] <- tags$p(class = "text-warning small",
                           "Note: RNA-seq TPM matrix not found - re-run prepare_titan_inputs.R to enable dynamic TPM threshold filtering.")
    do.call(tagList, lines)
  })
}

# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
