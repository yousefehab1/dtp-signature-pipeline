# ==============================================================================
# dashboard/R/mod_about.R - About, Methods & Citations Module
#
# Presents methodological rationale, mathematical definitions, live configuration
# gating thresholds, signature citations, and session info.
# ==============================================================================

#' About / Methods UI Module
#'
#' @param id Module namespace identifier.
#' @return A Shiny UI tag list.
#' @export
mod_about_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    bslib::layout_column_wrap(
      width = 1/2,

      # Methods Summary Card
      bslib::card(
        bslib::card_header(shiny::strong("Methodology Summary")),
        bslib::card_body(
          shiny::tags$h6("3-Year Landmark Survival Endpoints", class = "fw-bold text-primary"),
          shiny::tags$p(
            "Overall Survival (OS) and Recurrence-Free Survival (RFS) are evaluated as 3-year (36-month) ",
            "landmark endpoints. Kaplan-Meier log-rank survival curves split cohorts at the median score, ",
            "and continuous Cox proportional hazards models report Hazard Ratios scaled per 1 standard deviation (SD)."
          ),
          shiny::tags$h6("Single-Sample GSEA (ssGSEA)", class = "fw-bold text-primary mt-3"),
          shiny::tags$p(
            "Single-sample Gene Set Enrichment Analysis (ssGSEA) is performed via ",
            shiny::tags$code("GSVA"),
            " to quantify coordinated gene expression across positive and negative signature sets for each patient."
          ),
          shiny::tags$h6("Composite Score Definition", class = "fw-bold text-primary mt-3"),
          shiny::tags$p(
            "Composite scores evaluate the difference between positive and negative signature terms: ",
            shiny::tags$br(),
            shiny::tags$code("Score_Composite = Score_Positive - Score_Negative"),
            shiny::tags$br(),
            "This subtraction cancels shared baseline transcriptomic variance and amplifies phenotypic signal."
          ),
          shiny::tags$h6("Molecular Subtyping", class = "fw-bold text-primary mt-3"),
          shiny::tags$p(
            "Consensus Molecular Subtypes (CMS1-4) are classified using ",
            shiny::tags$code("CMScaller"),
            ", and Pathway-Derived Subtypes (PDS1-3) are classified using ",
            shiny::tags$code("PDSclassifier"),
            "."
          )
        )
      ),

      # Live vs Precomputed Tradeoff & Live Gating Thresholds
      bslib::card(
        bslib::card_header(shiny::strong("Architecture & Quality Gate Thresholds")),
        bslib::card_body(
          shiny::tags$h6("Live Re-scoring vs Precomputed Tradeoff", class = "fw-bold text-primary"),
          shiny::tags$div(
            class = "alert alert-info",
            shiny::tags$strong("Interactive Scope: "),
            "The ", shiny::tags$strong("Signature Explorer"), " tab recomputes statistical tests dynamically on-the-fly ",
            "from the precomputed single-sample score vectors of all signatures in the panel. ",
            "In contrast, the ", shiny::tags$strong("Pan-Cancer Survival"), " (ssGSEA across 25 cohorts, n=9,264) and ",
            shiny::tags$strong("Metastasis DE/GSEA"), " arms require whole-transcriptome computation and are presented ",
            "as static, publication-ready precomputed composite figures and tables. ",
            "Adding a genuinely new gene signature requires editing ",
            shiny::tags$code("inst/signatures/panel.csv"),
            " and executing ",
            shiny::tags$code("pipeline/build_results_bundle.R"),
            " once."
          ),
          shiny::tags$h6("Pipeline Quality Gate Thresholds (Loaded Live)", class = "fw-bold text-primary mt-3"),
          DT::dataTableOutput(ns("tbl_cfg"))
        )
      )
    ),

    # Signature Citations Table
    bslib::card(
      bslib::card_header(shiny::strong("Signature Panel Citations & Provenance")),
      DT::dataTableOutput(ns("tbl_citations")),
      bslib::card_footer(
        class = "text-muted small",
        "Source citations extracted from inst/signatures/panel.csv for all active signatures."
      )
    ),

    # Session Info Accordion
    bslib::card(
      bslib::card_header(shiny::strong("Runtime Session Information")),
      shiny::verbatimTextOutput(ns("session_info"))
    )
  )
}

#' About / Methods Server Module
#'
#' @param id Module namespace identifier.
#' @param bundle Precomputed results bundle loaded from `results_bundle.rds`.
#' @export
mod_about_server <- function(id, bundle) {
  shiny::moduleServer(id, function(input, output, session) {

    # Live config table
    output$tbl_cfg <- DT::renderDataTable({
      if (is.null(bundle) || is.null(bundle$cfg)) {
        return(DT::datatable(data.frame(Message = "Configuration not available in bundle"), rownames = FALSE))
      }

      cfg_df <- data.frame(
        Parameter = names(bundle$cfg),
        Value = vapply(bundle$cfg, function(x) paste(as.character(x), collapse = ", "), character(1)),
        stringsAsFactors = FALSE
      )

      param_descriptions <- c(
        score_suffix    = "Suffix appended to signature score columns",
        crc_modifiers   = "Subgroup and modifier variables evaluated in CRC",
        os_cutpoint     = "Overall Survival landmark cutoff (months)",
        rfs_cutpoint    = "Recurrence-Free Survival landmark cutoff (months)",
        min_group_n     = "Minimum patients per discrete group for Wilcoxon test",
        min_km_group_n  = "Minimum patients per arm for Kaplan-Meier log-rank",
        min_events      = "Minimum required survival events for Cox regression",
        min_cox_n       = "Minimum total cohort sample size for Cox models",
        min_epv         = "Minimum events per variable required in multivariable Cox"
      )

      cfg_df$Description <- ifelse(cfg_df$Parameter %in% names(param_descriptions),
                                  param_descriptions[cfg_df$Parameter], "Pipeline parameter")

      DT::datatable(
        cfg_df,
        options = list(
          dom = "t",
          paging = FALSE,
          ordering = FALSE,
          autoWidth = TRUE
        ),
        rownames = FALSE
      )
    })

    # Citations table
    output$tbl_citations <- DT::renderDataTable({
      if (is.null(bundle) || is.null(bundle$panel_tbl)) {
        return(DT::datatable(data.frame(Message = "Panel citations not available"), rownames = FALSE))
      }

      p_tbl <- bundle$panel_tbl
      has_cite <- "Source_Citation" %in% colnames(p_tbl)

      df <- p_tbl %>%
        dplyr::group_by(Signature_ID, Signature_Name) %>%
        dplyr::summarise(
          Genes = dplyr::n(),
          Source_Citation = if (has_cite) dplyr::first(stats::na.omit(Source_Citation)) else "N/A",
          .groups = "drop"
        )

      DT::datatable(
        df,
        options = list(
          pageLength = 10,
          dom = "tip",
          autoWidth = TRUE
        ),
        rownames = FALSE
      )
    })

    # Preformatted session info
    output$session_info <- shiny::renderPrint({
      utils::sessionInfo()
    })
  })
}
