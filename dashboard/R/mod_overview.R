# ==============================================================================
# dashboard/R/mod_overview.R - Overview & Landing Page Shiny Module
#
# Presents pipeline run metadata, panel coverage, data availability status,
# and high-level orientation for each analysis section.
# ==============================================================================

#' Overview UI Module
#'
#' @param id Module namespace identifier.
#' @return A Shiny UI tag list.
#' @export
mod_overview_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::uiOutput(ns("status_alert")),

    # Overview header and metrics
    bslib::layout_column_wrap(
      width = 1/4,
      bslib::value_box(
        title = "Panel Signatures",
        value = shiny::textOutput(ns("vb_sig_count")),
        showcase = NULL,
        shiny::textOutput(ns("vb_sig_sub")),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Bundle Built",
        value = shiny::textOutput(ns("vb_built_at")),
        showcase = NULL,
        shiny::textOutput(ns("vb_built_sub")),
        theme = "primary"
      ),
      bslib::value_box(
        title = "CRC Cohorts",
        value = shiny::textOutput(ns("vb_crc_status")),
        showcase = NULL,
        shiny::textOutput(ns("vb_crc_sub")),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Pan-Cancer & Mets",
        value = shiny::textOutput(ns("vb_pancan_status")),
        showcase = NULL,
        shiny::textOutput(ns("vb_pancan_sub")),
        theme = "primary"
      )
    ),

    # Orientation cards
    bslib::layout_column_wrap(
      width = 1/2,
      bslib::card(
        bslib::card_header(shiny::strong("Interactive Exploration")),
        bslib::card_body(
          shiny::tags$h6("Signature Explorer", class = "fw-bold text-primary"),
          shiny::tags$p(
            "Perform live, on-the-fly re-scoring of any canonical signature or custom ",
            shiny::tags$code("Positive - Negative"),
            " signature pair across whole-cohort and adjuvant-treated CRC cohorts (GSE39582, TCGA-COAD). ",
            "Dynamically computes Wilcoxon rank-sum effect sizes, Kaplan-Meier median-split log-rank curves, ",
            "and continuous univariate Cox proportional hazards metrics."
          ),
          shiny::tags$h6("Comprehensive Data Tables & Excel Export", class = "fw-bold text-primary mt-3"),
          shiny::tags$p(
            "Filter, search, and sort all precomputed analysis tables across CRC survival, subgroup analyses, ",
            "pan-cancer continuous Cox modeling, and metastasis GSEA/ROC assessments. Directly download ",
            "the complete multi-sheet publication Excel report workbook (.xlsx)."
          )
        )
      ),
      bslib::card(
        bslib::card_header(shiny::strong("Precomputed Presentation Layer")),
        bslib::card_body(
          shiny::tags$h6("CRC Survival & Molecular Subtyping", class = "fw-bold text-primary"),
          shiny::tags$p(
            "Static presentation of 7 publication-ready composite figures (Fig 1 outcome associations, ",
            "Fig 2 Kaplan-Meier panels, Fig 3 subgroup forest plots, Fig 3.3A-C molecular subtype characterizations, ",
            "and Fig 5 confounding summaries) with complete false discovery rate (BH-FDR) corrections."
          ),
          shiny::tags$h6("Pan-Cancer Survival & Metastasis DE/GSEA", class = "fw-bold text-primary mt-3"),
          shiny::tags$p(
            "Precomputed continuous-score Cox proportional hazards modeling across 25 solid tumor TCGA cohorts (Table 3.5), ",
            "batch-corrected aggregate survival estimates, chemotherapy/radiation-treated subsets, and paired ",
            "primary-vs-metastasis differential expression and GSEA (GSE50760) with liver purity adjustment."
          )
        )
      )
    ),

    # Panel Summary Table
    bslib::card(
      bslib::card_header(shiny::strong("Active Signature Panel Summary")),
      DT::dataTableOutput(ns("panel_table")),
      bslib::card_footer(
        class = "text-muted small",
        "Configurable signature panel loaded from inst/signatures/panel.csv and evaluated across all pipeline modules."
      )
    )
  )
}

#' Overview Server Module
#'
#' @param id Module namespace identifier.
#' @param bundle Precomputed results bundle loaded from `results_bundle.rds`.
#' @export
mod_overview_server <- function(id, bundle) {
  shiny::moduleServer(id, function(input, output, session) {

    output$status_alert <- shiny::renderUI({
      if (is.null(bundle)) {
        bslib::card(
          class = "bg-danger text-white mb-3",
          bslib::card_body(
            "⚠️ Results bundle missing! Please run `pipeline/build_results_bundle.R` to generate precomputed data."
          )
        )
      } else {
        NULL
      }
    })

    output$vb_sig_count <- shiny::renderText({
      if (is.null(bundle) || is.null(bundle$panel_tbl)) return("0")
      n_sigs <- length(unique(bundle$panel_tbl$Signature_ID))
      sprintf("%d Signatures", n_sigs)
    })
    output$vb_sig_sub <- shiny::renderText({
      if (is.null(bundle) || is.null(bundle$panel_tbl)) return("No panel loaded")
      n_genes <- nrow(bundle$panel_tbl)
      sprintf("%d total gene entries", n_genes)
    })

    output$vb_built_at <- shiny::renderText({
      if (is.null(bundle) || is.null(bundle$built_at)) return("Unknown")
      format(bundle$built_at, "%b %d, %Y")
    })
    output$vb_built_sub <- shiny::renderText({
      if (is.null(bundle) || is.null(bundle$built_at)) return("")
      format(bundle$built_at, "%H:%M:%S %Z")
    })

    output$vb_crc_status <- shiny::renderText({
      if (is.null(bundle)) return("Unavailable")
      has_gse <- !is.null(bundle$gse) && nrow(bundle$gse) > 0
      has_tcga <- !is.null(bundle$tcga) && nrow(bundle$tcga) > 0
      if (has_gse && has_tcga) "Available (2 Cohorts)" else if (has_gse) "GSE39582 Only" else "Unavailable"
    })
    output$vb_crc_sub <- shiny::renderText({
      if (is.null(bundle)) return("")
      n_gse <- if (!is.null(bundle$gse)) nrow(bundle$gse) else 0
      n_tcga <- if (!is.null(bundle$tcga)) nrow(bundle$tcga) else 0
      sprintf("GSE: %d | TCGA: %d samples", n_gse, n_tcga)
    })

    output$vb_pancan_status <- shiny::renderText({
      if (is.null(bundle)) return("Unavailable")
      has_pancan <- !is.null(bundle$pancan_stats_df) || !is.null(bundle$pancan_table_tbl)
      has_mets <- !is.null(bundle$mets_gsea_pn) || !is.null(bundle$mets_pvm_cmp)
      if (has_pancan && has_mets) "Available" else if (has_pancan) "Pan-Cancer Only" else if (has_mets) "Mets Only" else "Not in bundle"
    })
    output$vb_pancan_sub <- shiny::renderText({
      if (is.null(bundle)) return("")
      has_pancan <- !is.null(bundle$pancan_table_tbl)
      has_mets <- !is.null(bundle$mets_gsea_pn)
      sprintf("Pancan: %s | Mets: %s", ifelse(has_pancan, "Yes", "No"), ifelse(has_mets, "Yes", "No"))
    })

    output$panel_table <- DT::renderDataTable({
      if (is.null(bundle) || is.null(bundle$panel_tbl)) {
        return(DT::datatable(data.frame(Message = "No signature panel available"), rownames = FALSE))
      }

      df <- bundle$panel_tbl %>%
        dplyr::group_by(Signature_ID, Signature_Name) %>%
        dplyr::summarise(
          Gene_Count = dplyr::n(),
          Source_Citation = if ("Source_Citation" %in% colnames(bundle$panel_tbl)) dplyr::first(stats::na.omit(Source_Citation)) else "N/A",
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
  })
}
