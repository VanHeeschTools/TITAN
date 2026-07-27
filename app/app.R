## TITAN - Tumor Immunopeptidomics Target Atlas of Non-canonical ORFs
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
  library(jsonlite)
  library(shinyWidgets)
})

# Explicit source guards — no-op if Shiny already auto-sourced these files.
if (!exists("titan_theme")) source("global.R")
if (!exists("fmt1"))        for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

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
        choices  = c("Tumor-only", "Tumor-enriched", "Non-specific", "Unavailable"),
        selected = c("Tumor-only", "Tumor-enriched", "Non-specific", "Unavailable"),
        multiple = TRUE,
        options  = pickerOptions(
          actionsBox = TRUE, container = "body", dropupAuto = FALSE,
          selectedTextFormat = "count > 3", countSelectedText = "{0} / {1} types"
        )
      ),
      tags$p(class = "text-muted small mb-1 fw-semibold mt-2", "Off-tissue risk filter"),
      pickerInput(
        "prio_risk_filter", label = NULL,
        choices  = c("Safe", "Acceptable", "Borderline", "Critical", "Unavailable"),
        selected = c("Safe", "Acceptable", "Borderline", "Critical", "Unavailable"),
        multiple = TRUE,
        options  = pickerOptions(
          actionsBox = TRUE, container = "body", dropupAuto = FALSE,
          selectedTextFormat = "count > 4", countSelectedText = "{0} / {1} tiers"
        )
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
    style = "cursor:pointer;",
    onclick = "document.querySelector('#main_nav [data-value=\"Data\"]').click();",
    tags$img(src = "titan_logo_white.svg", height = "26px",
             style = "margin-right:8px; vertical-align:middle;")
  ),
  theme        = titan_theme,
  window_title = "TITAN",
  lang         = "en",
  footer       = tags$footer(
    class = "d-flex align-items-center gap-3 px-3 py-2",
    style = "background:#f8f9fa;border-top:1px solid #dee2e6;font-size:12px;color:#6c757d;",
    tags$img(src = "vanHeesch_logo_petrol.svg", height = "25px",
             alt = "Van Heesch Lab", style = "opacity:0.85;"),
    tags$a("vanheeschlab.com", href = "https://vanheeschlab.com", target = "_blank",
           style = "color:#28646E;text-decoration:none;"),
    tags$span("·", style = "color:#ccc;"),
    tags$a("TITAN on GitHub", href = "https://github.com/VanHeeschTools/TITAN",
           target = "_blank", style = "color:#28646E;text-decoration:none;"),
    tags$span("·", style = "color:#ccc;"),
    tags$span("2026")
  ),
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
                 "A named R list produced by ", tags$code("prepare_titan_inputs.R"),
                 ". Required element:"),
          tags$ul(class = "small text-muted ps-3 mb-1",
            tags$li(tags$code("$orf_table"),
                    " — data.frame of ncORF candidates. Backbone columns: ",
                    tags$code("orf_id"), ", ", tags$code("gene_id"), ", ",
                    tags$code("gene_name"), ", ", tags$code("gene_biotype"), ", ",
                    tags$code("orf_biotype_single"), ", ", tags$code("protein_seq"), ", ",
                    tags$code("protein_length"), ", ", tags$code("start_codon"), ", ",
                    tags$code("chr"), ", ", tags$code("orf_start"), ", ",
                    tags$code("orf_end"), ", ", tags$code("strand"), ", ",
                    tags$code("gene_id_clean"), ". Computed columns by source: ",
                    tags$em("target translation"),
                    " (", tags$code("target_translation_*"), "); ",
                    tags$em("GTEx"),
                    " (", tags$code("GTEX_tumor_only"), ", ",
                    tags$code("GTEX_tumor_enriched"), ", ",
                    tags$code("GTEX_max_median_TPM"), ", ",
                    tags$code("GTEX_tissues_q3_gt1"), "); ",
                    tags$em("target expression"),
                    " (", tags$code("target_expression_*"), "); ",
                    tags$em("TCGA"),
                    " (", tags$code("TCGA_tumor_*"), ", ", tags$code("TCGA_normal_*"), "); ",
                    tags$em("Ribocrypt"),
                    " (", tags$code("ribocrypt_primary_*"), ", ",
                    tags$code("ribocrypt_cell-line_*"), ").")
          ),
          tags$p(class = "small text-muted mb-1", "Per-sample matrices (rownames = orf_id or gene_id):"),
          tags$ul(class = "small text-muted ps-3 mb-1",
            tags$li(tags$code("$ribo_ppm_samples"),
                    " — ORF × target ribo-seq sample matrix (PPM)"),
            tags$li(tags$code("$rna_tpm_mat"),
                    " — gene × target RNA-seq sample matrix (TPM); gene_id ",
                    tags$em("without version suffix")),
            tags$li(tags$code("$gtex_tpm_mat"),
                    " — gene × GTEx normal sample matrix (TPM)"),
            tags$li(tags$code("$tcga_tpm_mat"),
                    " — gene × TCGA sample matrix (TPM; tumour + peritumoral)"),
            tags$li(tags$code("$ribocrypt_mat"),
                    " — ORF × Ribocrypt sample matrix (PPM; primary + cell-line)")
          ),
          tags$p(class = "small text-muted mb-1", "Sample metadata data.frames:"),
          tags$ul(class = "small text-muted ps-3 mb-1",
            tags$li(tags$code("$ribo_sample_meta"),
                    " — ", tags$code("sample_id"), ", ", tags$code("condition")),
            tags$li(tags$code("$rna_sample_meta"),
                    " — ", tags$code("sample_id"), ", ", tags$code("tissue_type"),
                    ", ", tags$code("condition")),
            tags$li(tags$code("$gtex_sample_meta"),
                    " — ", tags$code("sample_id"), ", ", tags$code("tissue_type")),
            tags$li(tags$code("$tcga_sample_meta"),
                    " — ", tags$code("sample_id"), ", ", tags$code("tissue_type"),
                    ", ", tags$code("sample_type"), ", ", tags$code("group")),
            tags$li(tags$code("$ribocrypt_sample_meta"),
                    " — ", tags$code("sample_id"), ", ", tags$code("group")),
            tags$li(tags$code("$ribocrypt_meta"),
                    " — list with ", tags$code("$primary_samples"),
                    " and ", tags$code("$cell_line_samples"), " character vectors")
          ),
          tags$ul(class = "small text-muted ps-3",
            tags$li(tags$code("$prepared_on"), " — POSIXct timestamp (optional)")
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
                "dplyr", "tidyr", "stringr", "ggplot2", "shinyWidgets",
                "Biostrings", "rBLAST"),
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
          tags$h1(class = "titan-hero-title", "Tumor Immunopeptidomics Target Atlas of Non‑canonical ORFs")
        )
      ),

      layout_columns(
        col_widths = c(6, 6),
        gap = "1rem",

        card(
          card_header(
            class = "d-flex align-items-center justify-content-between",
            tags$span(icon("file-arrow-up"), " ORF candidates"),
            uiOutput("orf_status_badge", inline = TRUE)
          ),
          card_body(uiOutput("orf_source_ui"))
        ),

        card(
          card_header(
            class = "d-flex align-items-center justify-content-between",
            tags$span(icon("vials"), " MS peptides"),
            uiOutput("ms_status_badge", inline = TRUE)
          ),
          card_body(
            div(class = "d-flex align-items-end gap-2",
              div(class = "flex-grow-1",
                fileInput("ms_file", "Upload MS results file",
                          accept = c(".csv", ".tsv", ".txt"),
                          buttonLabel = "Browse…",
                          placeholder = "peptides.csv / .tsv")
              ),
              div(class = "mb-3", uiOutput("clear_ms_btn"))
            ),
            div(class = "d-flex align-items-center gap-2 mb-3",
                tags$span(class = "text-muted small", "- or —"),
                actionButton("load_demo_ms", "Load demo peptides",
                             icon = icon("flask"), class = "btn-sm btn-outline-primary")),
            uiOutput("ms_load_status"),
            uiOutput("col_selector")
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

      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Protein sequence & peptide matches"),
          card_body(class = "p-2", uiOutput("detail_protein_seq_ui"))
        ),
        card(
          card_header(
            class = "d-flex align-items-center justify-content-between",
            tags$span(icon("shield-halved"), " Cross-reactivity"),
            tags$small(class = "text-muted fst-italic me-1", "Ensembl 114 pep")
          ),
          card_body(
            class = "p-2",
            tags$strong(class = "small text-secondary d-block mb-1", "Canonical self-proteins"),
            uiOutput("xreact_status_ui"),
            DTOutput("xreact_hits_dt"),
            tags$hr(class = "my-2"),
            tags$strong(class = "small text-secondary d-block mb-1", "Allergens (non-self)"),
            uiOutput("allergen_status_ui"),
            DTOutput("allergen_hits_dt")
          )
        )
      ),

      card(
        class = "mt-2",
        card_header(
          class = "d-flex align-items-center justify-content-between",
          tags$span(icon("dna"), " BLAST homology"),
          tags$small(class = "text-muted fst-italic me-1", "Ensembl 114 pep · ≥50% id, ≥30% cov")
        ),
        card_body(
          class = "p-2",
          uiOutput("blast_status_ui"),
          DTOutput("blast_hits_dt")
        )
      )
    )
  ),

  # ── Tab 5: Report ───────────────────────────────────────────────────────────
  nav_panel(
    "Report", icon = icon("file-pdf"),
    div(class = "container-fluid py-3",
      layout_columns(
        col_widths = c(2, 8, 2),
        NULL,
        card(
          card_body(
            uiOutput("report_status_ui"),
            tags$hr(class = "my-3"),
            tags$p(class = "fw-semibold mb-1", "Generate report"),
            tags$p(class = "text-muted small mb-1", "Plot values in:"),
            radioButtons("rpt_scale", NULL,
                         choices  = c("log(TPM+1) / log(PPM+1)" = "log", "Raw values" = "raw"),
                         selected = "log", inline = TRUE),
            div(class = "mb-3",
              downloadButton("dl_report", " Generate report",
                             icon  = icon("download"),
                             class = "btn-primary")
            ),
            tags$p(class = "text-muted small mb-0",
                   icon("circle-info"), " Downloads a self-contained HTML file.",
                   " Open it in any browser and use ", tags$b("File → Print"),
                   " to save as PDF (", tags$b("A3 landscape"), " · no margins · fit to page).")
          )
        ),
        NULL
      )
    )
  ),

  # ── Tab 6: About ────────────────────────────────────────────────────────────
  nav_panel(
    "About", icon = icon("circle-question"),
    card(
      card_body(
        tags$h4("TITAN - Tumor Immunopeptidomics Target Atlas for Non-canonical ORFs"),
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
        tags$h6("ORF Detail — safety checks"),
        tags$ul(
          tags$li(tags$b("Canonical cross-reactivity (Biostrings):"),
                  " exact and 1-mismatch peptide matching against the Ensembl 114 proteome",
                  " (~60–80 K deduplicated sequences, all annotated isoforms); results collapsed",
                  " to gene level (ENSG) with isoform count."),
          tags$li(tags$b("Allergen cross-reactivity (Biostrings):"),
                  " same matching strategy against UniProt KW-0020 reviewed allergen proteins."),
          tags$li(tags$b("BLAST homology (blastp):"),
                  " full-protein search against the Ensembl 114 deduplicated proteome;",
                  " hits filtered to ≥ 50 % identity AND ≥ 30 % alignment coverage;",
                  " annotated with gene description from an offline Ensembl 114 biomaRt lookup.")
        ),
        tags$p(class = "text-muted small mb-0",
               "Reference databases are built once via ",
               tags$code("app/scripts/01_prep_ensembl_pep.sbatch"),
               ", ",
               tags$code("02_prep_allergen.sbatch"),
               ", and ",
               tags$code("03_prep_annotation.sbatch"),
               ". All checks run offline at app runtime — no internet access required."),
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
  app_data_rv    <- reactiveVal(NULL)
  show_upload_rv <- reactiveVal(FALSE)

  # ── Cross-reactivity / BLAST state ───────────────────────────────────────────
  # Per-orf session caches; all keyed by orf_id.
  xreact_cache_rv   <- reactiveVal(list())
  allergen_cache_rv <- reactiveVal(list())
  blast_cache_rv    <- reactiveVal(list())
  modal_gene_rv     <- reactiveVal(NULL)   # list(gid=ENSG, sym=gene_symbol)

  # BLAST fires 750 ms after the user stops changing candidates.
  # xreact + allergen observe input$detail_orf_id directly (fast, no debounce needed).
  detail_orf_id_debounced <- debounce(reactive(input$detail_orf_id), 750)

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
    input$prio_risk_filter
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

  # ── Data tab: study library ──────────────────────────────────────────────────

  filtered_catalog <- reactive({
    df <- STUDY_CATALOG
    if (nrow(df) == 0L) return(df)
    q  <- trimws(input$catalog_search %||% "")
    ct <- input$catalog_ct_filter     %||% "ALL"
    co <- input$catalog_cohort_filter %||% "ALL"
    if (nzchar(q))         df <- df[grepl(q,  df$display_name, ignore.case = TRUE), ]
    if (ct != "ALL")       df <- df[!is.na(df$cancer_type) & df$cancer_type == ct, ]
    if (co != "ALL")       df <- df[!is.na(df$cohort)      & df$cohort      == co, ]
    df
  })

  output$catalog_study_list <- renderUI({
    df <- filtered_catalog()
    if (nrow(df) == 0L)
      return(tags$p(class = "text-muted small mt-2",
                    if (nrow(STUDY_CATALOG) == 0L)
                      "No studies in catalog. Run scripts/register_study.R to add one."
                    else
                      "No studies match the current filter."))
    rows <- lapply(seq_len(nrow(df)), function(i) {
      s   <- df[i, ]
      sid <- s$study_id
      loaded <- !is.null(app_data_rv()) &&
                identical(app_data_rv()$study_id %||% "", sid)

      if (loaded) {
        # ── Active study: expanded card with left accent ──
        div(class = "titan-study-active",
          div(class = "d-flex justify-content-between align-items-start",
            div(
              tags$b(class = "d-block", s$display_name),
              div(class = "d-flex gap-2 mt-1 flex-wrap align-items-center",
                if (!is.na(s$n_orfs))
                  tags$span(class = "badge bg-primary",
                            paste0(formatC(s$n_orfs, big.mark = ",", format = "d"), " ORFs")),
                if (!is.na(s$n_ribo_samples))
                  tags$span(class = "badge bg-secondary",
                            paste0(s$n_ribo_samples, " samples")),
                tags$span(class = "badge rounded-pill text-bg-light border", s$cancer_type),
                if (!is.na(s$prepared_on))
                  tags$small(class = "text-muted", s$prepared_on)
              )
            ),
            actionButton("clear_rds", NULL, icon = icon("trash"),
                         class = "btn-sm btn-outline-danger flex-shrink-0",
                         title = "Clear ORF data")
          )
        )
      } else {
        # ── Available study: compact row ──
        div(class = "titan-study-row d-flex justify-content-between align-items-center",
          div(
            tags$b(class = "d-block", s$display_name),
            tags$small(class = "text-muted",
              paste0(
                if (!is.na(s$n_orfs))
                  paste0(formatC(s$n_orfs, big.mark = ",", format = "d"), " ORFs"),
                if (!is.na(s$n_ribo_samples))
                  paste0(" · ", s$n_ribo_samples, " ribo / ", s$n_rna_samples, " RNA")
              )
            )
          ),
          actionButton(paste0("load_study_", sid), "Load",
                       class = "btn-sm btn-outline-primary flex-shrink-0")
        )
      }
    })
    div(class = "mt-1", tagList(rows))
  })

  output$orf_source_ui <- renderUI({
    show_upload <- show_upload_rv()
    ct_choices     <- c("All cancer types" = "ALL",
                        sort(unique(na.omit(STUDY_CATALOG$cancer_type))))
    cohort_choices <- c("All cohorts" = "ALL",
                        sort(unique(na.omit(STUDY_CATALOG$cohort))))
    div(
      # Toggle buttons: Study Library | Upload Data
      div(class = "d-flex gap-2 mb-3",
        actionButton("show_library_btn", tagList(icon("book-open"), " Study Library"),
                     class = paste("btn-sm",
                                   if (!show_upload) "titan-toggle-active" else "titan-toggle-inactive")),
        actionButton("show_upload_btn", tagList(icon("upload"), " Upload Data"),
                     class = paste("btn-sm",
                                   if (show_upload) "titan-toggle-active" else "titan-toggle-inactive"))
      ),
      if (!show_upload) {
        div(
          textInput("catalog_search", NULL, placeholder = "Search studies…", width = "100%"),
          div(class = "d-flex gap-2 mb-2",
            selectInput("catalog_ct_filter", NULL, width = "150px", choices = ct_choices),
            selectInput("catalog_cohort_filter", NULL, width = "140px", choices = cohort_choices)
          ),
          div(style = "max-height:320px;overflow-y:auto;padding-right:2px;",
            uiOutput("catalog_study_list")
          )
        )
      } else {
        div(class = "mt-1",
          fileInput("user_rds_file", NULL, accept = ".rds",
                    buttonLabel = "Browse…",
                    placeholder = "titan_<study_id>.rds")
        )
      }
    )
  })

  # One observer per catalog entry, registered at session start
  lapply(STUDY_CATALOG$study_id, function(sid) {
    observeEvent(input[[paste0("load_study_", sid)]], {
      entry <- STUDY_CATALOG[STUDY_CATALOG$study_id == sid, ]
      withProgress(message = paste0("Loading ", entry$display_name, "…"), value = 0.2, {
        setProgress(0.6, detail = "Reading RDS (may take ~10 s)…")
        dat <- tryCatch(
          readRDS(entry$rds_path),
          error = function(e) {
            showNotification(paste0("Could not read RDS: ", e$message), type = "error",
                             duration = 10)
            NULL
          }
        )
        if (is.null(dat)) return()
        setProgress(0.85, detail = "Validating…")
        err <- tryCatch({ validate_titan_rds(dat); NULL }, error = function(e) e$message)
        if (!is.null(err)) {
          showNotification(paste0("Invalid study data: ", err), type = "error", duration = 10)
          return()
        }
        if (is.null(dat$study_id)) dat$study_id <- sid
        setProgress(1)
        app_data_rv(dat)
      })
    }, ignoreInit = TRUE)
  })

  # ── Data tab: upload handler ─────────────────────────────────────────────────

  observeEvent(input$user_rds_file, {
    req(input$user_rds_file)
    withProgress(message = "Loading dataset…", value = 0.2, {
      setProgress(0.6, detail = "Reading RDS…")
      dat <- tryCatch(readRDS(input$user_rds_file$datapath), error = function(e) NULL)
      setProgress(0.85, detail = "Validating…")
      if (is.null(dat)) {
        showNotification("Could not read file as RDS.", type = "error")
        return()
      }
      err <- tryCatch({ validate_titan_rds(dat); NULL }, error = function(e) e$message)
      if (!is.null(err)) {
        showNotification(paste0("Invalid RDS: ", err), type = "error", duration = 10)
        return()
      }
      setProgress(1)
      app_data_rv(dat)
    })
  })

  observeEvent(input$show_library_btn, {
    show_upload_rv(FALSE)
  }, ignoreInit = TRUE)

  observeEvent(input$show_upload_btn, {
    show_upload_rv(TRUE)
  }, ignoreInit = TRUE)

  output$ms_load_status <- renderUI({
    ms <- tryCatch(ms_data(), error = function(e) NULL)
    if (is.null(ms)) return(NULL)
    n_pep   <- nrow(ms)
    is_demo <- is.null(user_ms_rv())
    div(class = "d-flex align-items-center gap-2 mb-2",
        tags$span(class = "badge bg-secondary",
                  paste0(formatC(n_pep, big.mark = ","), " peptides")),
        if (is_demo) tags$span(class = "badge bg-info", "demo"))
  })

  output$clear_ms_btn <- renderUI({
    ms <- tryCatch(ms_data(), error = function(e) NULL)
    if (is.null(ms)) return(NULL)
    actionButton("clear_ms", "Clear", icon = icon("trash"),
                 class = "btn-sm btn-outline-danger",
                 title = "Clear MS data")
  })

  output$orf_status_badge <- renderUI({
    if (!is.null(app_data_rv()))
      tags$span(class = "titan-status-badge titan-status-ready", "Ready")
    else
      tags$span(class = "titan-status-badge titan-status-awaiting", "Awaiting data")
  })

  output$ms_status_badge <- renderUI({
    ms_ok <- tryCatch(!is.null(ms_data()), error = function(e) FALSE)
    if (ms_ok)
      tags$span(class = "titan-status-badge titan-status-ready", "Ready")
    else
      tags$span(class = "titan-status-badge titan-status-awaiting", "Awaiting data")
  })

  observeEvent(input$clear_rds, {
    app_data_rv(NULL)
    all_matches_rv(NULL)
    started_rv(FALSE)
    show_upload_rv(FALSE)
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
        actionButton("start_titan", "EXPLORE TARGETS",
                     icon  = icon("play-circle"),
                     class = "btn-lg btn-primary btn-titan-ready px-5 py-2 fw-bold")
      )
    } else if (started_rv()) {
      md <- tryCatch(matched_data(), error = function(e) NULL)
      counts_str <- if (!is.null(md) && nrow(md) > 0L) {
        n_pep  <- n_distinct(md$matched_peptide)
        n_orfs <- n_distinct(md$orf_id)
        paste0(formatC(n_orfs, big.mark = ","), " ORFs and ",
               formatC(n_pep,  big.mark = ","), " peptides matched.")
      } else NULL
      div(
        class = "mt-4 text-center text-success",
        icon("circle-check", class = "fa-2x"),
        tags$p(class = "mt-2 mb-1 fw-semibold",
               "Running — use the tabs above to explore."),
        if (!is.null(counts_str))
          tags$p(class = "text-muted small mb-0", counts_str)
      )
    } else {
      div(
        class = "mt-4 text-center",
        actionButton("start_titan", "EXPLORE TARGETS",
                     icon  = icon("play-circle"),
                     class = "btn-lg btn-primary px-5 py-2 fw-bold",
                     style = "opacity:0.38; pointer-events:none;")
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
    ms        <- ms_data()
    # Base-R rename: avoids !! / tidy-eval entirely, works on any dplyr version
    names(ms)[names(ms) == input$pep_col] <- "matched_peptide"
    score_col <- intersect(PSM_QUALITY_COLS, names(ms))[1L]
    grp <- group_by(ms, matched_peptide)
    # Resolve branch outside dplyr so slice_max receives a concrete expression
    if (!is.na(score_col)) {
      grp %>% slice_max(order_by = .data[[score_col]], n = 1L, with_ties = FALSE) %>% ungroup()
    } else {
      grp %>% slice(1L) %>% ungroup()
    }
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
    gpd <- gene_prioritised_data()
    tq  <- if ("GTEX_tissues_q3_gt1" %in% names(gpd)) gpd$GTEX_tissues_q3_gt1 else NA_character_
    df <- gpd %>%
      mutate(
        spec_category = case_when(
          is.na(GTEX_tumor_only) | is.na(GTEX_tumor_enriched) ~ "Unavailable",
          GTEX_tumor_only %in% TRUE                           ~ "Tumor-only",
          GTEX_tumor_enriched %in% TRUE                       ~ "Tumor-enriched",
          TRUE                                                ~ "Non-specific"
        ),
        score_html   = vapply(priority_score, score_bar_html, character(1)),
        spec_html    = mapply(spec_badge_html, GTEX_tumor_only, GTEX_tumor_enriched),
        biotype_html        = vapply(orf_biotype_single, biotype_badge_html, character(1)),
        off_tissue_label    = mapply(off_tissue_risk, GTEX_tumor_only, tq),
        off_tissue_html     = vapply(off_tissue_label, off_tissue_risk_html, character(1)),
        transl_html  = vapply(target_translation_pct_samples,
                              function(x) pct_bar_html(x, "#28646E"), character(1)),
        expr_html    = vapply(target_expression_pct_samples,
                              function(x) pct_bar_html(x, "#7EB8BF"), character(1))
      )
    sel_spec <- input$prio_spec_filter %||% c("Tumor-only", "Tumor-enriched", "Non-specific")
    sel_risk <- input$prio_risk_filter %||% c("Safe", "Acceptable", "Borderline", "Critical", "Unavailable")
    df <- filter(df,
                 spec_category %in% sel_spec,
                 off_tissue_label %in% sel_risk)
    # DT column layout (0-based, after dropping .orf_id + orf_ids before DT):
    # Sel(0) Gene(1) ORF-biotype(2) Peptides(3) ORF-id(4) Location(5)
    # Specificity(6) Off-tissue risk(7) Score(8) Transl.%(9) Transl.PPM(10) Expr.%(11) Expr.TPM(12)
    # GTEx(13) TCGA T%(14) TCGA T TPM(15) TCGA N%(16) TCGA N TPM(17)
    # RC prim%(18) RC prim PPM(19) RC CL%(20) RC CL PPM(21)
    # .biotype_sort(22) .spec_sort(23) .off_tissue_sort(24) .score_sort(25) .transl_sort(26) .expr_sort(27) .child_rows(28)
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
      Specificity        = spec_html,
      `Off-tissue risk`  = off_tissue_html,
      Score          = score_html,
      `Transl. %`    = transl_html,
      `Transl. PPM`  = round(target_translation_median_PPM,  2),
      `Expr. %`      = expr_html,
      `Expr. TPM`    = round(target_expression_median_TPM,   2),
      `TCGA T%`      = round(TCGA_tumor_pct_samples,         1),
      `TCGA T TPM`   = round(TCGA_tumor_median_TPM,          2),
      `TCGA N%`      = round(TCGA_normal_pct_samples,        1),
      `TCGA N TPM`   = round(TCGA_normal_median_TPM,         2),
      `RC prim %`    = round(ribocrypt_primary_pct_samples,  1),
      `RC prim PPM`  = round(ribocrypt_primary_median_PPM,   2),
      `RC CL %`      = round(`ribocrypt_cell-line_pct_samples`, 1),
      `RC CL PPM`    = round(`ribocrypt_cell-line_median_PPM`,  2),
      .biotype_sort    = orf_biotype_single,
      .spec_sort       = spec_category,
      .off_tissue_sort = unname(.RISK_SORT[off_tissue_label]),
      .score_sort      = round(priority_score, 2),
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
    # Location(5) Specificity(6) Off-tissue risk(7) Score(8) Transl.%(9) PPM(10) Expr.%(11) TPM(12)
    # TCGAT%(13) TCGATPM(14) TCGAN%(15) TCANPM(16) RCprim%(17) RCprimPPM(18)
    # RCCL%(19) RCCLPPM(20) .biotype_sort(21) .spec_sort(22) .off_tissue_sort(23)
    # .score_sort(24) .transl_sort(25) .expr_sort(26) .child_rows(27)
    n_vis <- ncol(df) - 7L   # 7 hidden cols: 6 sort + .child_rows
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
          list(visible   = FALSE, targets = seq(n_vis, n_vis + 6L)),
          list(orderData = n_vis,      targets = 2L),   # ORF-biotype
          list(orderData = n_vis + 1L, targets = 6L),   # Specificity
          list(orderData = n_vis + 2L, targets = 7L),   # Off-tissue risk
          list(orderData = n_vis + 3L, targets = 8L),   # Score
          list(orderData = n_vis + 4L, targets = 9L),   # Transl. %
          list(orderData = n_vis + 5L, targets = 11L),  # Expr. %
          list(orderable = FALSE, targets = 0L)         # Sel not sortable
        ),
        lengthMenu = list(c(10, 20, 50), c("10", "20", "50"))
      )
    )
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
                   n_orfs_label)
          ),
          tags$p(class = "prio-detail-pep text-muted small mb-0", n_pep_label)
        ),
        div(class = "d-flex flex-column align-items-end gap-1 flex-shrink-0",
          actionButton("close_prio_detail", label = NULL, icon = icon("xmark"),
                       class = "btn-sm btn-outline-primary",
                       title = "Close detail panel"),
          div(class = "text-end mt-1",
            HTML(paste0(
              spec_badge_html(row$GTEX_tumor_only, row$GTEX_tumor_enriched), " ",
              off_tissue_risk_html(off_tissue_risk(row$GTEX_tumor_only, row$GTEX_tissues_q3_gt1))
            )),
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

  # When the user navigates away from the Prioritisation tab, close the detail
  # panel. The close button lives inside a renderUI: if the tab is hidden while
  # the panel is open, Shiny's resume-cycle Shiny.bindAll() pass can clear the
  # button's input binding without re-registering it (because renderUI didn't
  # re-execute), leaving the button stuck and non-functional on return.
  # Resetting prio_row_id here guarantees the panel is empty when the tab is
  # re-shown, so there is never a stale binding to worry about.
  observeEvent(input$main_nav, {
    if (!isTRUE(input$main_nav == "Prioritisation")) prio_row_id(NULL)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

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

  # ── Expression modal (from xreact / BLAST hit row clicks) ─────────────────
  output$modal_expr_plot <- renderPlotly({
    mg <- modal_gene_rv()
    req(!is.null(mg))
    gid       <- mg$gid
    log_scale <- isTRUE((input$modal_expr_scale %||% "log") == "log")
    y_label   <- if (log_scale) "log(TPM+1)" else "TPM"
    y_ref     <- if (log_scale) log(2) else 1

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

    # Plot 1: Target tumor
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
      p1 <- no_data_plot("Target tumor", "No RNA-seq matrix") %>%
        layout(yaxis = list(title = y_label))
    }

    # Plot 2: GTEx
    gtex_mat  <- tryCatch(gtex_tpm_rv(), error = function(e) NULL)
    gtex_meta <- tryCatch(gtex_meta_rv(), error = function(e) NULL)
    if (!is.null(gtex_mat) && isTRUE(gid %in% rownames(gtex_mat))) {
      gtex_raw   <- as.numeric(gtex_mat[gid, ])
      tissue_raw <- gtex_meta$tissue_type[match(colnames(gtex_mat), gtex_meta$sample_id)]
      tissue_raw[is.na(tissue_raw)] <- "Unknown"
      tissue     <- gsub("_", " ", tissue_raw)
      gtex_y     <- apply_scale(gtex_raw)
      med_by_tis <- tapply(gtex_y, tissue, median, na.rm = TRUE)
      sorted_tis <- names(sort(med_by_tis, decreasing = TRUE))
      tis_col_map <- setNames(
        sapply(unique(tissue_raw), function(t) {
          idx <- match(t, gtex_colors_subtissue$Tissue)
          if (!is.na(idx)) return(gtex_colors_subtissue$ColorHex[idx])
          gtex_keys <- names(gtex_colors)
          hits_c <- gtex_keys[startsWith(t, gtex_keys)]
          if (length(hits_c)) gtex_colors[[hits_c[which.max(nchar(hits_c))]]] else "#BBBBBB"
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
      p2 <- no_data_plot("GTEx", "GTEx matrix not in app data")
    }

    # Plot 3: TCGA
    tcga_mat  <- tryCatch(tcga_tpm_rv(), error = function(e) NULL)
    tcga_meta <- tryCatch(tcga_meta_rv(), error = function(e) NULL)
    if (!is.null(tcga_mat) && isTRUE(gid %in% rownames(tcga_mat))) {
      tcga_raw    <- as.numeric(tcga_mat[gid, ])
      grp         <- tcga_meta$group[match(colnames(tcga_mat), tcga_meta$sample_id)]
      grp[is.na(grp)] <- "Unknown"
      tcga_y      <- apply_scale(tcga_raw)
      unique_grps <- unique(grp)
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
      grp_factor  <- factor(grp, levels = grp_order)
      df_c <- data.frame(g = grp_factor, y = tcga_y)
      p3 <- plot_ly()
      for (lvl in grp_order) {
        d <- df_c[as.character(df_c$g) == lvl, , drop = FALSE]
        if (nrow(d) == 0L) next
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
      p3 <- no_data_plot("TCGA", "TCGA matrix not in app data")
    }

    out <- subplot(p1, p2, p3, nrows = 1, shareY = TRUE, titleX = TRUE,
                   widths = c(0.12, 0.53, 0.35)) %>%
      layout(
        paper_bgcolor = "white", plot_bgcolor = "white",
        margin = list(l = 5, r = 5, t = 10, b = 5),
        font   = list(family = "Inter", size = 11)
      ) %>%
      config(
        toImageButtonOptions = list(format = "svg", filename = "expression_modal"),
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

  # ── Cross-reactivity — Ensembl 114 pep (self), auto-trigger, immediate ───────
  # Results stored gene-level: one row per (peptide, gene) after collapsing isoforms.
  # Columns: Peptide, Gene_sym, ENSG, Mismatches, N_isoforms

  observeEvent(input$detail_orf_id, {
    req(input$detail_orf_id)
    oid <- input$detail_orf_id
    modal_gene_rv(NULL)

    # xreact
    if (is.null(xreact_cache_rv()[[oid]])) {
      md   <- tryCatch(matched_data(), error = function(e) NULL)
      peps <- if (!is.null(md)) unique(md$matched_peptide[md$orf_id == oid]) else character(0)

      store_xr <- function(result) {
        cache <- xreact_cache_rv(); cache[[oid]] <- result; xreact_cache_rv(cache)
      }

      if (length(peps) == 0L) {
        store_xr(data.frame())
      } else if (is.null(ensembl_pep_index) || length(ensembl_pep_index$seqs) == 0L) {
        store_xr(data.frame(Error = "Ensembl 114 pep index not loaded — run scripts/01_prep_ensembl_pep.sbatch first."))
      } else {
        ref_set <- Biostrings::AAStringSet(ensembl_pep_index$seqs)
        ref_md5 <- names(ensembl_pep_index$seqs)   # md5 hashes as names
        hits <- do.call(rbind, Filter(Negate(is.null), lapply(peps, function(pep) {
          pep_aa <- tryCatch(Biostrings::AAString(pep), error = function(e) NULL)
          if (is.null(pep_aa)) return(NULL)
          m0 <- tryCatch(Biostrings::vmatchPattern(pep_aa, ref_set, max.mismatch = 0L, fixed = TRUE),
                         error = function(e) NULL)
          m1 <- tryCatch(Biostrings::vmatchPattern(pep_aa, ref_set, max.mismatch = 1L, fixed = TRUE),
                         error = function(e) NULL)
          exact_md5 <- if (!is.null(m0)) ref_md5[which(lengths(m0) > 0L)] else character(0)
          near_md5  <- if (!is.null(m1)) setdiff(ref_md5[which(lengths(m1) > 0L)], exact_md5) else character(0)
          if (!length(exact_md5) && !length(near_md5)) return(NULL)
          collapse_to_genes <- function(md5s, mm) {
            if (!length(md5s)) return(NULL)
            rows <- do.call(rbind, lapply(md5s, function(md5) {
              ensg_vec <- ensembl_pep_index$md5_to_ensg[[md5]]
              sym_vec  <- ensembl_pep_index$md5_to_sym[[md5]]
              n_ensp   <- length(ensembl_pep_index$md5_to_ensp[[md5]])
              data.frame(Peptide = pep, Gene_sym = sym_vec %||% "unknown",
                         ENSG = ensg_vec %||% NA_character_,
                         Mismatches = mm, N_isoforms = n_ensp,
                         stringsAsFactors = FALSE)
            }))
            rows[!duplicated(rows$ENSG), , drop = FALSE]
          }
          rbind(collapse_to_genes(exact_md5, 0L), collapse_to_genes(near_md5, 1L))
        })))
        # Final dedup across peptides: keep worst (lowest) Mismatches per gene
        if (!is.null(hits) && nrow(hits) > 0L) {
          hits <- hits[order(hits$ENSG, hits$Mismatches), ]
          hits <- hits[!duplicated(hits$ENSG), ]
        }
        store_xr(if (is.null(hits)) data.frame() else hits)
      }
    }

    # allergen
    if (is.null(allergen_cache_rv()[[oid]])) {
      md   <- tryCatch(matched_data(), error = function(e) NULL)
      peps <- if (!is.null(md)) unique(md$matched_peptide[md$orf_id == oid]) else character(0)

      store_al <- function(result) {
        cache <- allergen_cache_rv(); cache[[oid]] <- result; allergen_cache_rv(cache)
      }

      if (length(peps) == 0L) {
        store_al(data.frame())
      } else if (is.null(allergen_index) || length(allergen_index$seqs) == 0L) {
        store_al(data.frame(Error = "Allergen index not loaded — run scripts/02_prep_allergen.sbatch first."))
      } else {
        al_set   <- Biostrings::AAStringSet(allergen_index$seqs)
        al_names <- names(allergen_index$seqs)  # "ENTRY|ACC"
        al_hits <- do.call(rbind, Filter(Negate(is.null), lapply(peps, function(pep) {
          pep_aa <- tryCatch(Biostrings::AAString(pep), error = function(e) NULL)
          if (is.null(pep_aa)) return(NULL)
          m0 <- tryCatch(Biostrings::vmatchPattern(pep_aa, al_set, max.mismatch = 0L, fixed = TRUE),
                         error = function(e) NULL)
          m1 <- tryCatch(Biostrings::vmatchPattern(pep_aa, al_set, max.mismatch = 1L, fixed = TRUE),
                         error = function(e) NULL)
          exact_idx <- if (!is.null(m0)) which(lengths(m0) > 0L) else integer(0)
          near_idx  <- if (!is.null(m1)) setdiff(which(lengths(m1) > 0L), exact_idx) else integer(0)
          if (!length(exact_idx) && !length(near_idx)) return(NULL)
          nm <- al_names[c(exact_idx, near_idx)]
          acc      <- sub(".*\\|", "", nm)
          gene_sym <- ifelse(acc %in% names(allergen_index$acc_to_sym),
                             allergen_index$acc_to_sym[acc], sub("\\|.*", "", nm))
          data.frame(Peptide = pep, Hit = nm, Gene_sym = gene_sym,
                     Mismatches = c(rep(0L, length(exact_idx)), rep(1L, length(near_idx))),
                     stringsAsFactors = FALSE)
        })))
        store_al(if (is.null(al_hits)) data.frame() else al_hits)
      }
    }
  }, ignoreInit = FALSE, ignoreNULL = TRUE)

  output$xreact_status_ui <- renderUI({
    req(input$detail_orf_id)
    xr <- xreact_cache_rv()[[input$detail_orf_id]]
    if (is.null(xr))
      return(div(class = "d-flex align-items-center gap-2 text-muted small mb-1",
                 tags$span(class = "spinner-border spinner-border-sm"), " Checking…"))
    if ("Error" %in% names(xr))
      return(tags$p(class = "text-danger small mb-1", icon("circle-xmark"), " ", xr$Error[1L]))
    if (nrow(xr) == 0L)
      return(tags$p(class = "text-success small mb-1",
                    icon("circle-check"), " No canonical matches (exact or 1 mismatch)."))
    n_exact <- sum(xr$Mismatches == 0L, na.rm = TRUE)
    n_near  <- sum(xr$Mismatches == 1L, na.rm = TRUE)
    label <- paste0(
      if (n_exact > 0L) paste0(n_exact, " exact gene(s)"),
      if (n_exact > 0L && n_near > 0L) ", ",
      if (n_near  > 0L) paste0(n_near, " 1-mismatch gene(s)")
    )
    tags$p(class = "fw-semibold small text-warning mb-1",
           icon("triangle-exclamation"), " ", label,
           ". Click a row to view expression.")
  })

  output$xreact_hits_dt <- DT::renderDT({
    req(input$detail_orf_id)
    xr <- xreact_cache_rv()[[input$detail_orf_id]]
    if (is.null(xr) || "Error" %in% names(xr) || nrow(xr) == 0L) return(NULL)
    xrs <- xr[order(xr$Mismatches, xr$Gene_sym), ]
    DT::datatable(
      data.frame(Peptide = xrs$Peptide, Gene = xrs$Gene_sym, ENSG = xrs$ENSG,
                 Mismatches = xrs$Mismatches, Isoforms = xrs$N_isoforms,
                 stringsAsFactors = FALSE),
      selection = "single", rownames = FALSE,
      options   = list(pageLength = 5, dom = "tp", scrollX = TRUE),
      class     = "compact stripe"
    )
  })

  observeEvent(input$xreact_hits_dt_rows_selected, {
    idx <- input$xreact_hits_dt_rows_selected
    req(length(idx) > 0L, input$detail_orf_id)
    xr <- xreact_cache_rv()[[input$detail_orf_id]]
    req(!is.null(xr), nrow(xr) > 0L, !"Error" %in% names(xr))
    xrs      <- xr[order(xr$Mismatches, xr$Gene_sym), ]
    row      <- xrs[idx, ]
    gid      <- row$ENSG
    gene_sym <- row$Gene_sym
    if (is.na(gid) || !nzchar(gid)) {
      showNotification(paste0("No ENSG ID for ", gene_sym, " — cannot load expression."),
                       type = "message", duration = 4)
      return()
    }
    modal_gene_rv(list(gid = gid, sym = gene_sym))
    showModal(modalDialog(
      title = paste("Expression:", gene_sym), size = "xl", easyClose = TRUE,
      footer = modalButton("Close"),
      radioButtons("modal_expr_scale", "Y axis:",
                   choices = c("log(TPM+1)" = "log", "Raw TPM" = "raw"),
                   selected = "log", inline = TRUE),
      plotlyOutput("modal_expr_plot", height = "420px")
    ))
  })

  output$allergen_status_ui <- renderUI({
    req(input$detail_orf_id)
    al <- allergen_cache_rv()[[input$detail_orf_id]]
    if (is.null(al))
      return(div(class = "d-flex align-items-center gap-2 text-muted small mb-1",
                 tags$span(class = "spinner-border spinner-border-sm"), " Checking…"))
    if ("Error" %in% names(al))
      return(tags$p(class = "text-danger small mb-1", icon("circle-xmark"), " ", al$Error[1L]))
    if (nrow(al) == 0L)
      return(tags$p(class = "text-success small mb-1",
                    icon("circle-check"), " No allergen matches (exact or 1 mismatch)."))
    n_exact <- sum(al$Mismatches == 0L, na.rm = TRUE)
    n_near  <- sum(al$Mismatches == 1L, na.rm = TRUE)
    label <- paste0(
      if (n_exact > 0L) paste0(n_exact, " exact"),
      if (n_exact > 0L && n_near > 0L) ", ",
      if (n_near  > 0L) paste0(n_near, " 1-mismatch")
    )
    tags$p(class = "fw-semibold small text-danger mb-1",
           icon("triangle-exclamation"), " ", label, " allergen hit(s).")
  })

  output$allergen_hits_dt <- DT::renderDT({
    req(input$detail_orf_id)
    al <- allergen_cache_rv()[[input$detail_orf_id]]
    if (is.null(al) || "Error" %in% names(al) || nrow(al) == 0L) return(NULL)
    als <- al[order(al$Mismatches, al$Gene_sym), ]
    DT::datatable(
      data.frame(Peptide = als$Peptide, Allergen = als$Hit, Gene = als$Gene_sym,
                 Mismatches = als$Mismatches, stringsAsFactors = FALSE),
      selection = "none", rownames = FALSE,
      options   = list(pageLength = 5, dom = "tp", scrollX = TRUE),
      class     = "compact stripe"
    )
  })

  # ── BLAST homology — Ensembl 114 pep, debounced 750 ms ───────────────────────
  # Threshold: ≥50% identity AND ≥30% alignment coverage (coverage = alnlen/qlen).
  # Result columns: Gene_sym, ENSG, pident, coverage, E_value, Description.

  observeEvent(detail_orf_id_debounced(), {
    oid <- detail_orf_id_debounced()
    req(oid)
    if (!is.null(blast_cache_rv()[[oid]])) return()   # cache hit — skip
    orf_seq <- tryCatch(detail_orf()$protein_seq, error = function(e) NA_character_)
    if (is.na(orf_seq) || !nzchar(orf_seq)) return()

    store_bl <- function(result) {
      cache <- blast_cache_rv(); cache[[oid]] <- result; blast_cache_rv(cache)
    }

    db <- REF_DB_ENSEMBL
    if (!file.exists(paste0(db, ".pdb")) && !file.exists(paste0(db, ".phr"))) {
      store_bl(data.frame(Error = paste0("BLAST database not found: ", db,
                                         " — run scripts/01_prep_ensembl_pep.sbatch first.")))
      return()
    }
    if (!nzchar(Sys.which("blastp"))) {
      store_bl(data.frame(Error = "blastp not found in PATH."))
      return()
    }

    qlen <- nchar(orf_seq)
    hits <- tryCatch({
      bl  <- rBLAST::blast(db = db, type = "blastp")
      qry <- Biostrings::AAStringSet(c(query = orf_seq))
      predict(bl, qry, BLAST_args = "-evalue 0.001 -max_target_seqs 50")
    }, error = function(e) data.frame(Error = e$message))

    if ("Error" %in% names(hits)) { store_bl(hits); return() }
    if (is.null(hits) || nrow(hits) == 0L) { store_bl(data.frame()); return() }

    # Apply thresholds: ≥50% identity AND ≥30% coverage
    pident_num <- as.numeric(hits$pident)
    aln_len    <- as.integer(hits$length)
    coverage   <- round(aln_len / qlen * 100, 1)
    keep       <- !is.na(pident_num) & pident_num >= 50 & coverage >= 30
    hits       <- hits[keep, , drop = FALSE]
    pident_num <- pident_num[keep]
    coverage   <- coverage[keep]

    if (nrow(hits) == 0L) { store_bl(data.frame()); return() }

    # Back-map ENSP (sseqid) → ENSG + gene symbol via index
    ensp       <- hits$sseqid
    md5s       <- if (!is.null(ensembl_pep_index))
                    ensembl_pep_index$ensp_to_md5[ensp] else rep(NA_character_, length(ensp))
    ensg <- vapply(md5s, function(md5) {
      if (is.na(md5) || is.null(ensembl_pep_index$md5_to_ensg[[md5]])) return(NA_character_)
      ensembl_pep_index$md5_to_ensg[[md5]][1L]
    }, character(1))
    gene_sym <- vapply(md5s, function(md5) {
      if (is.na(md5) || is.null(ensembl_pep_index$md5_to_sym[[md5]])) return(NA_character_)
      ensembl_pep_index$md5_to_sym[[md5]][1L]
    }, character(1))

    # Annotate with description from offline gene annotation map
    desc <- if (!is.null(ensembl_gene_annot) && "ensembl_gene_id" %in% names(ensembl_gene_annot)) {
      idx_a <- match(ensg, ensembl_gene_annot$ensembl_gene_id)
      ifelse(!is.na(idx_a), ensembl_gene_annot$description[idx_a], NA_character_)
    } else rep(NA_character_, length(ensg))

    result <- data.frame(
      Gene_sym    = gene_sym,
      ENSG        = ensg,
      Identity    = paste0(round(pident_num, 1), "%"),
      pident      = pident_num,
      Coverage    = paste0(coverage, "%"),
      coverage    = coverage,
      E_value     = formatC(as.numeric(hits$evalue), format = "e", digits = 1),
      Description = ifelse(is.na(desc), "", desc),
      stringsAsFactors = FALSE
    )
    # Deduplicate to one row per gene (keep highest identity hit)
    result <- result[order(result$ENSG, -result$pident), ]
    result <- result[!duplicated(result$ENSG), ]
    store_bl(result)
  }, ignoreInit = FALSE, ignoreNULL = TRUE)

  output$blast_status_ui <- renderUI({
    req(input$detail_orf_id)
    br <- blast_cache_rv()[[input$detail_orf_id]]
    if (is.null(br))
      return(div(class = "d-flex align-items-center gap-2 text-muted small mb-1",
                 tags$span(class = "spinner-border spinner-border-sm"), " Running BLAST…"))
    if ("Error" %in% names(br))
      return(tags$p(class = "text-danger small mb-1", icon("circle-xmark"), " ", br$Error[1L]))
    if (nrow(br) == 0L)
      return(tags$p(class = "text-success small mb-1",
                    icon("circle-check"), " No homologous genes (≥50% identity, ≥30% coverage)."))
    tags$p(class = "fw-semibold small text-warning mb-1", icon("dna"), " ",
           sprintf("%d homologous gene(s). Click a row to view expression.", nrow(br)))
  })

  output$blast_hits_dt <- DT::renderDT({
    req(input$detail_orf_id)
    br <- blast_cache_rv()[[input$detail_orf_id]]
    if (is.null(br) || "Error" %in% names(br) || nrow(br) == 0L) return(NULL)
    DT::datatable(
      data.frame(Gene = br$Gene_sym, ENSG = br$ENSG, Identity = br$Identity,
                 Coverage = br$Coverage, E_value = br$E_value,
                 Description = br$Description, stringsAsFactors = FALSE),
      selection = "single", rownames = FALSE,
      options   = list(pageLength = 5, dom = "tp", scrollX = TRUE),
      class     = "compact stripe"
    )
  })

  observeEvent(input$blast_hits_dt_rows_selected, {
    idx <- input$blast_hits_dt_rows_selected
    req(length(idx) > 0L, input$detail_orf_id)
    br  <- blast_cache_rv()[[input$detail_orf_id]]
    req(!is.null(br), nrow(br) > 0L, !"Error" %in% names(br))
    row      <- br[idx, ]
    gid      <- row$ENSG
    gene_sym <- row$Gene_sym
    if (is.na(gid) || !nzchar(gid %||% "")) {
      showNotification(paste0("No ENSG ID for ", gene_sym, " — cannot load expression."),
                       type = "message", duration = 4)
      return()
    }
    modal_gene_rv(list(gid = gid, sym = gene_sym))
    showModal(modalDialog(
      title = paste("Expression:", gene_sym), size = "xl", easyClose = TRUE,
      footer = modalButton("Close"),
      radioButtons("modal_expr_scale", "Y axis:",
                   choices = c("log(TPM+1)" = "log", "Raw TPM" = "raw"),
                   selected = "log", inline = TRUE),
      plotlyOutput("modal_expr_plot", height = "420px")
    ))
  })

  # ── Report ───────────────────────────────────────────────────────────────────
  output$report_status_ui <- renderUI({
    sel <- prio_selected_rowids()
    n   <- length(sel)
    if (!isTRUE(started_rv()))
      return(div(class = "alert alert-secondary py-2",
                 icon("circle-info"), " Run START on the Data tab first."))
    if (n == 0L)
      return(div(class = "alert alert-warning py-2",
                 icon("triangle-exclamation"),
                 " No candidates selected. Use the checkboxes in the Prioritisation table."))
    div(class = "alert alert-success py-2",
        icon("circle-check"),
        sprintf(" %d candidate%s selected.", n, if (n == 1L) "" else "s"))
  })

  output$dl_report <- downloadHandler(
    filename = function() paste0("titan_report_", format(Sys.time(), "%Y%m%d_%H%M"), ".html"),
    content  = function(file) {
      sel <- prio_selected_rowids()
      shiny::validate(shiny::need(length(sel) > 0, "No candidates selected."),
                      shiny::need(isTRUE(started_rv()), "Run START first."))

      gdata <- gene_prioritised_data()
      rows  <- gdata[gdata$.row_id %in% sel, , drop = FALSE]
      shiny::validate(shiny::need(nrow(rows) > 0, "Selected rows not found."))

      rna_mat   <- tryCatch(rna_tpm_rv(),        error = function(e) NULL)
      rna_meta  <- tryCatch(rna_meta_rv(),        error = function(e) NULL)
      gtex_mat  <- tryCatch(gtex_tpm_rv(),        error = function(e) NULL)
      gtex_meta <- tryCatch(gtex_meta_rv(),       error = function(e) NULL)
      tcga_mat  <- tryCatch(tcga_tpm_rv(),        error = function(e) NULL)
      tcga_meta <- tryCatch(tcga_meta_rv(),       error = function(e) NULL)
      ribo_m    <- tryCatch(ribo_ppm_rv(),        error = function(e) NULL)
      ribo_sm   <- tryCatch(ribo_meta_rv(),       error = function(e) NULL)
      rc_mat    <- tryCatch(ribocrypt_mat_rv(),   error = function(e) NULL)
      rc_meta   <- tryCatch(ribocrypt_smeta_rv(), error = function(e) NULL)
      md        <- tryCatch(matched_data(),       error = function(e) NULL)
      ot        <- tryCatch(orf_table_rv(),       error = function(e) orf_table)
      log_scale <- isTRUE((input$rpt_scale %||% "log") == "log")
      gen_date  <- format(Sys.time(), "%Y-%m-%d %H:%M")

      # Attach protein sequence from orf_table (may not be in gene_prioritised_data).
      # Must use a single vectorised assignment — row-by-row assignment on a new column
      # triggers R recycling: NULL[1]<-x creates length-1 c(x), which when written back
      # as a column fills every row with the same value.
      rows$protein_seq <- ot$protein_seq[match(rows$orf_id, ot$orf_id)]
      if (!is.null(md)) {
        na_mask <- is.na(rows$protein_seq)
        if (any(na_mask))
          rows$protein_seq[na_mask] <- md$protein_seq[match(rows$orf_id[na_mask], md$orf_id)]
      }

      .read_svg_uri <- function(path) {
        if (!file.exists(path)) return("")
        raw_b <- readBin(path, "raw", n = file.info(path)$size)
        paste0("data:image/svg+xml;base64,", jsonlite::base64_enc(raw_b))
      }
      logo_uri    <- .read_svg_uri("www/titan_logo_blue.svg")
      vh_logo_uri <- .read_svg_uri("www/vanHeesch_logo_petrol.svg")

      withProgress(message = "Generating report…", value = 0, {
        pages  <- vector("list", nrow(rows) + 1L)
        preset <- active_preset() %||% "Custom"
        for (i in seq_len(nrow(rows))) {
          setProgress(i / (nrow(rows) + 1), detail = paste("Candidate", i, "of", nrow(rows)))
          row      <- rows[i, , drop = FALSE]
          pep_list <- if (!is.null(md))
            unique(md$matched_peptide[md$orf_id == row$orf_id])
          else character(0)
          pages[[i]] <- .rpt_build_page(
            row         = row,
            pep_list    = pep_list,
            rna_mat     = rna_mat,  rna_meta  = rna_meta,
            gtex_mat    = gtex_mat, gtex_meta = gtex_meta,
            tcga_mat    = tcga_mat, tcga_meta = tcga_meta,
            ribo_m      = ribo_m,   ribo_sm   = ribo_sm,
            rc_mat      = rc_mat,   rc_meta   = rc_meta,
            logo_uri    = logo_uri, vh_logo_uri = vh_logo_uri,
            gen_date    = gen_date, log_scale   = log_scale
          )
        }

        setProgress(1, detail = "Building summary page…")
        w <- current_weights()
        weight_rows <- lapply(WEIGHT_META, function(m)
          list(label = m$label, weight = as.numeric(w[[m$id]])))
        filter_rows <- list(
          list("PPM threshold (ribo-seq)",     as.character(input$ppm_threshold %||% 1)),
          list("Min. samples ≥ PPM threshold", as.character(input$ppm_n_samples %||% "")),
          list("TPM threshold (RNA-seq)",      as.character(input$tpm_threshold %||% 1)),
          list("Min. samples ≥ TPM threshold", as.character(input$tpm_n_samples %||% "")),
          list("Biotypes shown",               paste(input$biotype_filter %||% "all", collapse = ", ")),
          list("Specificity filter",           paste(input$prio_spec_filter %||% "all", collapse = ", "))
        )
        d         <- app_data_rv()
        data_info <- list(
          list("ORF table prepared",  if (!is.null(d$prepared_on)) format(d$prepared_on, "%Y-%m-%d %H:%M") else "—"),
          list("Total ORFs",          formatC(nrow(ot), big.mark = ",")),
          list("Ribo-seq samples",    if (!is.null(ribo_m))  as.character(ncol(ribo_m))  else "—"),
          list("RNA-seq samples",     if (!is.null(rna_mat)) as.character(ncol(rna_mat)) else "—"),
          list("Candidates selected", as.character(nrow(rows))),
          list("Scoring preset",      preset),
          list("Report generated",    gen_date)
        )
        pages[[nrow(rows) + 1L]] <- .rpt_build_summary_page(
          params_list = list(preset = preset, weights = weight_rows, filters = filter_rows),
          data_info   = data_info,
          logo_uri    = logo_uri, vh_logo_uri = vh_logo_uri,
          gen_date    = gen_date
        )
      })

      writeLines(.rpt_wrap_html(pages), file, useBytes = FALSE)
    }
  )

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
