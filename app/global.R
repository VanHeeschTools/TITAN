## TITAN - global.R
## Sourced by Shiny before app.R, or explicitly from app.R if auto-sourcing is unavailable.
## Contains: libraries, startup data, colour palettes, constants, theme.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(plotly)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(jsonlite)
  library(shinyWidgets)
  library(yaml)
  library(float)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

options(shiny.maxRequestSize = 1000 * 1024^2)  # 1 GB upload limit

# ─────────────────────────────────────────────────────────────────────────────
# STUDY CATALOG  (static at app startup; loaded from data/catalog.yaml)
# ─────────────────────────────────────────────────────────────────────────────

STUDY_CATALOG <- local({
  f <- "data/catalog.yaml"
  empty <- data.frame(
    study_id = character(), display_name = character(),
    cancer_type = character(), cohort = character(),
    rds_path = character(), n_orfs = integer(),
    n_ribo_samples = integer(), n_rna_samples = integer(),
    prepared_on = character(), stringsAsFactors = FALSE
  )
  if (!file.exists(f)) { message("No catalog found at ", f); return(empty) }
  entries <- tryCatch(yaml::read_yaml(f)$studies, error = function(e) {
    message("Failed to read catalog: ", e$message); NULL
  })
  if (is.null(entries) || length(entries) == 0L) return(empty)
  dplyr::bind_rows(lapply(entries, function(e) {
    as.data.frame(lapply(e, function(v) if (is.null(v)) NA else v),
                  stringsAsFactors = FALSE)
  }))
})

# Placeholder values used by UI sliders/pickers before any data is loaded.
# All are updated reactively via observeEvent(app_data_rv(), ...) in app.R.
n_ribo_samples <- 1L
n_rna_samples  <- 1L
biotypes       <- character(0)


# ─────────────────────────────────────────────────────────────────────────────
# CROSS-REACTIVITY REFERENCE DATA
# ─────────────────────────────────────────────────────────────────────────────

# Load Bioconductor packages without attaching to avoid masking dplyr::slice, dplyr::filter etc.
requireNamespace("Biostrings", quietly = TRUE)
requireNamespace("IRanges",   quietly = TRUE)
requireNamespace("rBLAST",    quietly = TRUE)

# Prepend bundled blastp binary so rBLAST finds it without any module or Singularity
local({
  bin <- normalizePath("bin", mustWork = FALSE)
  if (dir.exists(bin))
    Sys.setenv(PATH = paste0(bin, ":", Sys.getenv("PATH")))
})

# BLAST database: Ensembl 114 pep (deduplicated), same release as pipeline
REF_DB_ENSEMBL <- "ref/ensembl114_pep/ensembl114_pep"

# ── Ensembl 114 pep index (self-cross-reactivity + BLAST back-mapping) ────────
# Built by scripts/reference_prep/01_prep_ensembl_pep.R. List with:
#   $seqs         named character: md5 → AA sequence (deduplicated, ~60-80K)
#   $md5_to_ensg  list: md5 → character vector of ENSG IDs
#   $md5_to_sym   list: md5 → character vector of gene symbols
#   $ensp_to_md5  named vector: ENSP → md5 (BLAST sseqid back-lookup)
#   $rep_table    data.frame: one row per unique sequence
ensembl_pep_index <- local({
  f <- "ref/ensembl114_pep/ensembl114_pep_index.rds"
  if (!file.exists(f)) {
    message("Ensembl 114 pep index not found — run scripts/reference_prep/01_prep_ensembl_pep.R first")
    return(NULL)
  }
  idx <- readRDS(f)
  n_seqs <- length(idx$seqs)
  # Drop the sequence strings from global scope — they are only needed for the
  # in-memory exact-match cross-reactivity check. Load lazily on first use via
  # ensembl_pep_seqs_lazy() to avoid holding ~200-400 MB at startup.
  idx$seqs <- NULL
  message(sprintf("Ensembl 114 pep: index loaded (%d unique sequences; AA strings deferred)", n_seqs))
  idx
})

# Convert any float32 matrices in a study data list to plain double.
# Called right after readRDS() so the rest of the app never sees float32 objects.
fl32_to_dbl <- function(dat) {
  mat_names <- c("ribo_ppm_samples", "rna_tpm_mat", "gtex_tpm_mat",
                 "tcga_tpm_mat", "ribocrypt_mat")
  for (nm in mat_names)
    if (!is.null(dat[[nm]]) && inherits(dat[[nm]], "float32"))
      dat[[nm]] <- float::dbl(dat[[nm]])
  dat
}

# Returns the $seqs named character vector, loading from disk on first call.
ensembl_pep_seqs_lazy <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    f <- "ref/ensembl114_pep/ensembl114_pep_index.rds"
    if (!file.exists(f)) return(NULL)
    message("Ensembl 114 pep: loading AA sequences for cross-reactivity check…")
    cache <<- readRDS(f)$seqs
    cache
  }
})

# ── Ensembl gene annotation (offline lookup for BLAST hit annotation) ─────────
# Built by scripts/reference_prep/02_prep_annotation.R. data.frame columns:
#   ensembl_gene_id, external_gene_name, uniprotswissprot, description
ensembl_gene_annot <- local({
  f <- "ref/ensembl114_pep/ensembl_gene_annotation.rds"
  if (!file.exists(f)) {
    message("Ensembl gene annotation not found — run scripts/reference_prep/02_prep_annotation.R first")
    return(NULL)
  }
  df <- readRDS(f)
  message(sprintf("Ensembl gene annotation: %d entries loaded", nrow(df)))
  df
})

# ─────────────────────────────────────────────────────────────────────────────
# GENCODE ORF TABLE  (optional; enables cross-matching against TransCode Phase 2)
# ─────────────────────────────────────────────────────────────────────────────

source("R/fct_gencode_orf.R")

gencode_orf_tbl <- local({
  df <- load_gencode_orf_table("ref/gencode_orfs_phase2.csv")
  if (is.null(df)) return(NULL)
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
  Whole_Blood     = "#DD6661",

  # RETRO-ACTIVE CATALOG STUDIES ONLY: custom aggregate/renamed GTEx groups, not
  # part of the standard xstitchR Group vocabulary. Aliases reuse their
  # constituent group's color; aggregates get a new blended color (see above).
  digestive_system     = "#FDB378",
  female_reproductive  = "#D2A2A3",
  Adipose_Tissue       = "#D1B9A5",  # alias of Adipose
  Blood_Vessel         = "#FFB08E",  # alias of Artery
  Blood                = "#DD6661",  # alias of Whole_Blood
  Salivary_Gland       = "#C794C0"   # alias of Salivary
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
    "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
    # RETRO-ACTIVE CATALOG STUDIES ONLY: self-referential Tissue==Group rows for
    # studies whose DE data was pre-computed at this coarser grouping already
    "digestive_system", "female_reproductive",
    "Adipose_Tissue", "Blood_Vessel", "Blood", "Salivary_Gland"
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
    "Skin", "Skin",
    "digestive_system", "female_reproductive",
    "Adipose_Tissue", "Blood_Vessel", "Blood", "Salivary_Gland"
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
    "#E1CEA9", "#ECCBA8",
    "#FDB378", "#D2A2A3",
    "#D1B9A5", "#FFB08E", "#DD6661", "#C794C0"
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
  primary_corneal_eye           = "#C2B2A1",  # no GTEx match
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
  hTERT_RPE1_eye                = "#C2B2A1",  # no GTEx match (RPE)
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
# BIOTYPE COLOURS
# ─────────────────────────────────────────────────────────────────────────────

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
  "pseudogene-ORF"           = "#25CED1", 
  "novel-ORF"                 = "#2E3D46"
)

# ─────────────────────────────────────────────────────────────────────────────
# RISK CATEGORIES FOR OFF-TISSUE EXPRESSION (GTEx) — constants
# ─────────────────────────────────────────────────────────────────────────────

off_tissue_risk_adult <- c(
  Adipose            = "Acceptable",
  Adrenal_Gland      = "Critical",
  Artery             = "Critical",
  Bladder            = "Borderline",
  Brain              = "Critical",
  Breast             = "Borderline",
  Cervix             = "Acceptable",
  Colon              = "Critical",
  Esophagus          = "Critical",
  Fallopian_Tube     = "Acceptable",
  Heart              = "Critical",
  Kidney             = "Critical",
  Liver              = "Critical",
  Lung               = "Critical",
  Muscle             = "Borderline",
  Nerve              = "Critical",
  Ovary              = "Borderline",
  Pancreas           = "Critical",
  Pituitary          = "Critical",
  Prostate           = "Acceptable",
  Salivary           = "Acceptable",
  Skin               = "Acceptable",
  Small_Intestine    = "Critical",
  Spleen             = "Acceptable",
  Stomach            = "Critical",
  Thyroid            = "Borderline",
  Uterus             = "Acceptable",
  Vagina             = "Acceptable",
  Whole_Blood        = "Critical"
)

off_tissue_risk_pediatric <- c(
  Adipose            = "Acceptable",
  Adrenal_Gland      = "Critical",
  Artery             = "Critical",
  Bladder            = "Borderline",
  Brain              = "Critical",
  Breast             = "Critical",
  Cervix             = "Borderline",
  Colon              = "Critical",
  Esophagus          = "Critical",
  Fallopian_Tube     = "Borderline",
  Heart              = "Critical",
  Kidney             = "Critical",
  Liver              = "Critical",
  Lung               = "Critical",
  Muscle             = "Borderline",
  Nerve              = "Critical",
  Ovary              = "Critical",
  Pancreas           = "Critical",
  Pituitary          = "Critical",
  Prostate           = "Critical",
  Salivary           = "Borderline",
  Skin               = "Borderline",
  Small_Intestine    = "Critical",
  Spleen             = "Borderline",
  Stomach            = "Critical",
  Thyroid            = "Critical",
  Uterus             = "Borderline",
  Vagina             = "Borderline",
  Whole_Blood        = "Critical"
)

# ─────────────────────────────────────────────────────────────────────────────
# SCORING SYSTEM — constants
# ─────────────────────────────────────────────────────────────────────────────

# Each dimension exposes a signal in [0, 1].
# Weight = +1 → high signal adds to score; −1 → high signal penalises; 0 → ignored.
# hint describes what "high signal" means for each dimension.
WEIGHT_META <- list(
  list(id="w_pct_samples",  label="% expressed (RNA-seq)",    hint="many tumor samples express the ORF",             radar="Expr. %",    group="Tumor coverage",  pancancer=0.3, specific=0.3),
  list(id="w_pct_transl",   label="% translated (Ribo-seq)",  hint="many tumor samples translate the ORF",           radar="Transl. %",  group="Tumor coverage",  pancancer=0.3, specific=0.3),
  list(id="w_tumor_spec",   label="Tumor specificity (GTEx)", hint="GTEx: tumor-only=1, enriched=0.5, non-specific=0",                         radar="Specificity",group="Specificity",    pancancer=0.6, specific=1),
  list(id="w_off_tissue",   label="Off-tissue risk (GTEx)",   hint="Safety tier: Safe=1, Acceptable=0.75, Borderline=0.5, Critical=0",         radar="Off-tissue", group="Safety",        pancancer=0.7, specific=0.8),
  list(id="w_gtex_penalty", label="GTEx expression level",    hint="expressed in normal tissues (GTEx)",             radar="GTEx",       group="Specificity",     pancancer=-0.6,specific=-1),
  list(id="w_tcga_cov",     label="TCGA tumor coverage",      hint="expressed in many TCGA tumor samples",           radar="TCGA T",     group="TCGA validation",   pancancer=0.6, specific=-0.3),
  list(id="w_peri_penalty", label="TCGA normal expression",   hint="expressed in peritumoral / normal TCGA samples", radar="TCGA N",     group="TCGA validation",      pancancer=-0.6,specific=-0.6),
  list(id="w_ribo_primary", label="RC primary tissue",        hint="translated in normal primary tissues (Ribocrypt)",radar="RC primary", group="Normal tissue (Ribocrypt)", pancancer=-0.6,specific=-0.7),
  list(id="w_ribo_cell",    label="RC cell-line",             hint="translated in normal cell lines (Ribocrypt)",    radar="RC CL",      group="Normal tissue (Ribocrypt)", pancancer=0.6, specific=-0.6)
)

PRESETS <- list(
  "Cancer-specific" = list(label = "Strict tumor specificity, penalises normal tissue", color = "#FFBEFF"),
  "Pan-cancer"      = list(label = "Broad coverage, tolerates enriched targets",        color = "#2F3D46")
)

# ─────────────────────────────────────────────────────────────────────────────
# SHINY THEME
# ─────────────────────────────────────────────────────────────────────────────

titan_theme <- bs_theme(
  version       = 5,
  bg            = "#f8f9fa",
  fg            = "#212529",
  primary       = "#2F3D46",
  secondary     = "#AADCFF",
  success       = "#41A395",  # verdigris - must match --seed-success in www/titan.css
  warning       = "#F3C677",
  danger        = "#B33E3E",
  info          = "#FFBEFF",
  font_scale    = 0.9,
  base_font     = font_google("Inter"),
  heading_font  = "IBM Plex Sans",
  code_font     = "IBM Plex Mono",
  "card-bg"                 = "#FFFFFF",
  "navbar-light-bg"               = "#2F3D46",
  "navbar-light-color"            = "#F0F7F9",
  "navbar-light-hover-color"      = "#FFFFFF",
  "navbar-light-active-color"     = "#FFFFFF",
  "navbar-light-brand-color"      = "#FFFFFF",
  "navbar-light-brand-hover-color"= "#FFFFFF",
  "input-border-color"            = "#DAE9EC"
)
