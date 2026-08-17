# ==============================================================================
# dashboard/R/mod_tables.R - Comprehensive Data Tables & Excel Export Module
#
# Organizes all precomputed data tables into sub-tabs and provides direct
# download of the multi-sheet Excel report workbook.
# ==============================================================================

#' Tables UI Module
#'
#' @param id Module namespace identifier.
#' @return A Shiny UI tag list.
#' @export
mod_tables_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::uiOutput(ns("tables_alert")),

    # Download Banner Card
    bslib::card(
      bslib::card_body(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        shiny::tags$div(
          shiny::tags$h5("Analysis Results Tables & Reports", class = "mb-1 fw-bold"),
          shiny::tags$p(
            "Explore precomputed statistics across all analysis arms or download the consolidated 18-sheet Excel workbook.",
            class = "text-muted small mb-0"
          )
        ),
        shiny::uiOutput(ns("download_btn_ui"))
      )
    ),

    # Navset of all tables
    bslib::navset_card_tab(
      title = "Precomputed Data Tables",
      bslib::nav_panel(
        title = "Signature Panel",
        DT::dataTableOutput(ns("tbl_panel")),
        shiny::tags$div(class = "text-muted small mt-2", "Canonical signature panel definitions and gene lists.")
      ),
      bslib::nav_panel(
        title = "Composite Definitions",
        DT::dataTableOutput(ns("tbl_composites")),
        shiny::tags$div(class = "text-muted small mt-2", "Composite signature definitions and positive/negative term weights.")
      ),
      bslib::nav_panel(
        title = "CRC Univariable Survival",
        DT::dataTableOutput(ns("tbl_crc_stats")),
        shiny::tags$div(class = "text-muted small mt-2", "CRC univariate survival statistics (Wilcoxon, KM, Cox) with BH-FDR correction.")
      ),
      bslib::nav_panel(
        title = "Pan-Cancer Survival (Table 3.5)",
        DT::dataTableOutput(ns("tbl_pancan_surv")),
        shiny::tags$div(class = "text-muted small mt-2", "Continuous Cox proportional hazards models across 25 solid tumor TCGA cohorts.")
      ),
      bslib::nav_panel(
        title = "Pan-Cancer Treated",
        DT::dataTableOutput(ns("tbl_pancan_treated")),
        shiny::tags$div(class = "text-muted small mt-2", "Survival models restricted to chemotherapy- or radiation-treated patients.")
      ),
      bslib::nav_panel(
        title = "Mets GSEA (Normal vs Primary)",
        DT::dataTableOutput(ns("tbl_mets_pn")),
        shiny::tags$div(class = "text-muted small mt-2", "GSEA comparing Primary CRC against Normal mucosa (GSE50760).")
      ),
      bslib::nav_panel(
        title = "Mets GSEA (PvM Comparison)",
        DT::dataTableOutput(ns("tbl_mets_pvm")),
        shiny::tags$div(class = "text-muted small mt-2", "Primary vs Metastasis GSEA NES before and after liver score adjustment.")
      ),
      bslib::nav_panel(
        title = "Mets Liver Purity",
        DT::dataTableOutput(ns("tbl_mets_liver")),
        shiny::tags$div(class = "text-muted small mt-2", "Per-sample liver purity ssGSEA scores across tissue types.")
      ),
      bslib::nav_panel(
        title = "Mets ROC AUC",
        DT::dataTableOutput(ns("tbl_mets_roc")),
        shiny::tags$div(class = "text-muted small mt-2", "Single-sample paired ROC AUC distinguishing Primary CRC from Liver Metastasis.")
      ),
      bslib::nav_panel(
        title = "Mets Single-Sample Wilcoxon",
        DT::dataTableOutput(ns("tbl_mets_wilcox")),
        shiny::tags$div(class = "text-muted small mt-2", "Paired Wilcoxon test statistics and BH-adjusted p-values for ssGSEA scores.")
      )
    )
  )
}

#' Tables Server Module
#'
#' @param id Module namespace identifier.
#' @param bundle Precomputed results bundle loaded from `results_bundle.rds`.
#' @export
mod_tables_server <- function(id, bundle) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$tables_alert <- shiny::renderUI({
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

    # Locate Excel workbook path safely
    resolved_excel <- shiny::reactive({
      if (is.null(bundle)) return(NULL)
      excel_cand <- bundle$excel_path
      candidates <- c(
        excel_cand,
        if (!is.null(excel_cand)) file.path("..", excel_cand) else NULL,
        "output/DTP_Results.xlsx",
        "../output/DTP_Results.xlsx"
      )
      for (cand in candidates) {
        if (!is.null(cand) && file.exists(cand)) {
          return(cand)
        }
      }
      NULL
    })

    # Render Download Button UI
    output$download_btn_ui <- shiny::renderUI({
      excel_file <- resolved_excel()
      if (!is.null(excel_file) && file.exists(excel_file)) {
        shiny::downloadButton(
          ns("download_excel"),
          label = "Download Excel Report (.xlsx)",
          class = "btn btn-success"
        )
      } else {
        shiny::tags$button(
          "Excel Report Not Available",
          class = "btn btn-secondary disabled",
          disabled = "disabled",
          title = "Run pipeline/build_results_bundle.R to build the Excel report."
        )
      }
    })

    # Download handler
    output$download_excel <- shiny::downloadHandler(
      filename = function() {
        paste0("DTP_Results_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
      },
      content = function(file) {
        excel_file <- resolved_excel()
        if (!is.null(excel_file) && file.exists(excel_file)) {
          file.copy(excel_file, file)
        }
      }
    )

    render_dt <- function(data_obj, empty_msg = "Table not available in bundle") {
      if (is.null(data_obj) || nrow(as.data.frame(data_obj)) == 0) {
        return(DT::datatable(data.frame(Message = empty_msg), rownames = FALSE))
      }
      df <- as.data.frame(data_obj)
      num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
      for (col in num_cols) {
        if (grepl("P|FDR|p_adj|pvalue|qvalue", col, ignore.case = TRUE)) {
          df[[col]] <- signif(df[[col]], 3)
        } else {
          df[[col]] <- round(df[[col]], 3)
        }
      }
      DT::datatable(
        df,
        options = list(
          pageLength = 15,
          dom = "ftip",
          autoWidth = TRUE
        ),
        filter = "top",
        rownames = FALSE
      )
    }

    # Render all sub-tab tables
    output$tbl_panel          <- DT::renderDataTable(render_dt(bundle$panel_tbl, "Signature panel not available"))
    output$tbl_composites     <- DT::renderDataTable(render_dt(bundle$composite_defs, "Composite definitions not available"))
    output$tbl_crc_stats      <- DT::renderDataTable(render_dt(bundle$stats_df, "CRC survival stats not available"))
    output$tbl_pancan_surv    <- DT::renderDataTable(render_dt(if (!is.null(bundle$pancan_table_tbl)) bundle$pancan_table_tbl else bundle$pancan_stats_df, "Pan-cancer survival table not available"))
    output$tbl_pancan_treated <- DT::renderDataTable(render_dt(bundle$treated_stats_df, "Pan-cancer treated stats not available"))
    output$tbl_mets_pn        <- DT::renderDataTable(render_dt(bundle$mets_gsea_pn, "Mets GSEA (Normal vs Primary) table not available"))
    output$tbl_mets_pvm       <- DT::renderDataTable(render_dt(bundle$mets_pvm_cmp, "Mets GSEA (PvM comparison) table not available"))
    output$tbl_mets_liver     <- DT::renderDataTable(render_dt(bundle$mets_liver_df, "Mets Liver purity table not available"))
    output$tbl_mets_roc       <- DT::renderDataTable(render_dt(bundle$mets_roc_summary, "Mets ROC AUC table not available"))
    output$tbl_mets_wilcox    <- DT::renderDataTable(render_dt(bundle$mets_ssgsea_pvalues, "Mets ssGSEA Wilcoxon table not available"))
  })
}
