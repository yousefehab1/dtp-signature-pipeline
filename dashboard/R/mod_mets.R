# ==============================================================================
# dashboard/R/mod_mets.R - Metastasis DE & GSEA Presentation Module
#
# Displays metastasis composite figures (GSEA and Liver Purity) and
# associated summary statistics tables.
# ==============================================================================

#' Metastasis DE/GSEA UI Module
#'
#' @param id Module namespace identifier.
#' @return A Shiny UI tag list.
#' @export
mod_mets_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::uiOutput(ns("mets_alert")),

    bslib::navset_card_tab(
      title = "Metastasis Composite Figures",
      bslib::nav_panel(
        title = "Fig: Metastasis GSEA Composite",
        shiny::imageOutput(ns("fig_mets_gsea"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig: Liver Purity Composite",
        shiny::imageOutput(ns("fig_mets_liver"), height = "auto")
      )
    ),

    bslib::navset_card_tab(
      title = "Metastasis Analysis Tables",
      bslib::nav_panel(
        title = "Normal vs Primary GSEA",
        DT::dataTableOutput(ns("tbl_gsea_pn")),
        shiny::tags$div(class = "text-muted small mt-2", "Gene Set Enrichment Analysis comparing Primary CRC vs Normal Mucosa (GSE50760).")
      ),
      bslib::nav_panel(
        title = "PvM GSEA (Adjusted vs Unadjusted)",
        DT::dataTableOutput(ns("tbl_pvm_cmp")),
        shiny::tags$div(class = "text-muted small mt-2", "Comparison of GSEA Normalized Enrichment Scores (NES) and BH-FDR before and after liver score adjustment.")
      ),
      bslib::nav_panel(
        title = "Liver Purity Scores",
        DT::dataTableOutput(ns("tbl_liver_df")),
        shiny::tags$div(class = "text-muted small mt-2", "Per-sample liver-specific gene ssGSEA scores across Normal, Primary, and Metastasis tissues.")
      ),
      bslib::nav_panel(
        title = "Single-Sample ROC AUC",
        DT::dataTableOutput(ns("tbl_roc_summary")),
        shiny::tags$div(class = "text-muted small mt-2", "Paired single-sample ROC AUC distinguishing Primary Tumor from Liver Metastasis.")
      )
    )
  )
}

#' Metastasis DE/GSEA Server Module
#'
#' @param id Module namespace identifier.
#' @param bundle Precomputed results bundle loaded from `results_bundle.rds`.
#' @export
mod_mets_server <- function(id, bundle) {
  shiny::moduleServer(id, function(input, output, session) {

    # Missing data guard
    output$mets_alert <- shiny::renderUI({
      if (is.null(bundle)) {
        bslib::card(
          class = "bg-danger text-white mb-3",
          bslib::card_body(
            "⚠️ Results bundle missing! Please run `pipeline/build_results_bundle.R` to generate precomputed data."
          )
        )
      } else if (is.null(bundle$mets_gsea_pn) && is.null(bundle$mets_pvm_cmp)) {
        bslib::card(
          class = "bg-warning text-dark mb-3",
          bslib::card_body(
            "⚠️ Metastasis DE/GSEA data not found in this bundle. Please run `pipeline/build_results_bundle.R` to compute metastasis results."
          )
        )
      } else {
        NULL
      }
    })

    render_fig <- function(fig_name, alt_text) {
      shiny::renderImage({
        fig_dir <- if (!is.null(bundle$fig_dir)) bundle$fig_dir else "output/figures"
        cand_dirs <- c(
          fig_dir,
          file.path("..", fig_dir),
          "output/figures",
          "../output/figures"
        )
        found_path <- NULL
        for (d in cand_dirs) {
          p <- file.path(d, paste0(fig_name, ".png"))
          if (file.exists(p)) {
            found_path <- p
            break
          }
        }

        if (is.null(found_path) || !file.exists(found_path)) {
          tmp_ph <- file.path(tempdir(), paste0("placeholder_", fig_name, ".png"))
          if (!file.exists(tmp_ph)) {
            grDevices::png(tmp_ph, width = 800, height = 400, res = 100)
            graphics::par(mar = c(1, 1, 1, 1), bg = "#f8f9fa")
            graphics::plot.new()
            graphics::text(0.5, 0.55, paste0("Figure '", fig_name, ".png' not found on disk."), cex = 1.2, font = 2, col = "#6c757d")
            graphics::text(0.5, 0.40, "Run pipeline/build_results_bundle.R to generate figure files.", cex = 0.95, col = "#adb5bd")
            grDevices::dev.off()
          }
          found_path <- tmp_ph
        }

        list(
          src = found_path,
          contentType = "image/png",
          width = "100%",
          alt = alt_text
        )
      }, deleteFile = FALSE)
    }

    # 2 Metastasis figures
    output$fig_mets_gsea  <- render_fig("Fig_Mets_GSEA_composite", "Metastasis GSEA Composite")
    output$fig_mets_liver <- render_fig("Fig_Mets_LiverPurity_composite", "Liver Purity Composite")

    # Helper for formatting data tables
    render_table <- function(df_data, empty_msg) {
      if (is.null(df_data) || nrow(as.data.frame(df_data)) == 0) {
        return(DT::datatable(data.frame(Message = empty_msg), rownames = FALSE))
      }
      df <- as.data.frame(df_data)
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
          pageLength = 10,
          dom = "ftip",
          autoWidth = TRUE
        ),
        filter = "top",
        rownames = FALSE
      )
    }

    output$tbl_gsea_pn     <- DT::renderDataTable(render_table(bundle$mets_gsea_pn, "Normal vs Primary GSEA table not available"))
    output$tbl_pvm_cmp     <- DT::renderDataTable(render_table(bundle$mets_pvm_cmp, "PvM GSEA comparison table not available"))
    output$tbl_liver_df    <- DT::renderDataTable(render_table(bundle$mets_liver_df, "Liver purity score table not available"))
    output$tbl_roc_summary <- DT::renderDataTable(render_table(bundle$mets_roc_summary, "ROC AUC summary table not available"))
  })
}
