# ==============================================================================
# dashboard/R/mod_pancan.R - Pan-Cancer Survival Presentation Module
#
# Displays 6 pan-cancer composite figures, Table 3.5 summary, and treated
# patient survival statistics.
# ==============================================================================

#' Pan-Cancer Survival UI Module
#'
#' @param id Module namespace identifier.
#' @return A Shiny UI tag list.
#' @export
mod_pancan_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::uiOutput(ns("pancan_alert")),

    bslib::navset_card_tab(
      title = "Pan-Cancer Composite Figures",
      bslib::nav_panel(
        title = "Fig 6: Pan-Cancer Composite",
        shiny::imageOutput(ns("fig6_comp"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 6A: Up Forest",
        shiny::imageOutput(ns("fig6a"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 6B: Composite Forest",
        shiny::imageOutput(ns("fig6b"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 6C: Down Forest",
        shiny::imageOutput(ns("fig6c"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig Aggregate Forest",
        shiny::imageOutput(ns("fig_agg"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 7: Pan-Cancer Violins",
        shiny::imageOutput(ns("fig7"), height = "auto")
      )
    ),

    bslib::layout_column_wrap(
      width = 1,
      bslib::card(
        bslib::card_header(shiny::strong("Table 3.5: Pan-Cancer Cox Survival Analysis (All Patients)")),
        DT::dataTableOutput(ns("pancan_table")),
        bslib::card_footer(
          class = "text-muted small",
          "Continuous univariate Cox proportional hazards models across 25 solid tumor TCGA cohorts and batch-corrected aggregate meta-estimates."
        )
      ),
      bslib::card(
        bslib::card_header(shiny::strong("Pan-Cancer Survival in Chemotherapy/Radiation-Treated Patients")),
        DT::dataTableOutput(ns("treated_table")),
        bslib::card_footer(
          class = "text-muted small",
          "Univariable Cox survival models evaluated exclusively in patients documented to have received pharmaceutical or radiation therapy."
        )
      )
    )
  )
}

#' Pan-Cancer Survival Server Module
#'
#' @param id Module namespace identifier.
#' @param bundle Precomputed results bundle loaded from `results_bundle.rds`.
#' @export
mod_pancan_server <- function(id, bundle) {
  shiny::moduleServer(id, function(input, output, session) {

    # Missing data guard
    output$pancan_alert <- shiny::renderUI({
      if (is.null(bundle)) {
        bslib::card(
          class = "bg-danger text-white mb-3",
          bslib::card_body(
            "⚠️ Results bundle missing! Please run `pipeline/build_results_bundle.R` to generate precomputed data."
          )
        )
      } else if (is.null(bundle$pancan_table_tbl) && is.null(bundle$pancan_stats_df)) {
        bslib::card(
          class = "bg-warning text-dark mb-3",
          bslib::card_body(
            "⚠️ Pan-cancer survival data not found in this bundle. Please run `pipeline/build_results_bundle.R` to compute pan-cancer results."
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

    # 6 Pan-cancer figures
    output$fig6_comp <- render_fig("Fig6_pancan_survival_composite", "Fig 6: Pan-Cancer Survival Composite")
    output$fig6a     <- render_fig("Fig6A_pancan_forest_Up", "Fig 6A: Up Forest")
    output$fig6b     <- render_fig("Fig6B_pancan_forest_Composite", "Fig 6B: Composite Forest")
    output$fig6c     <- render_fig("Fig6C_pancan_forest_Down", "Fig 6C: Down Forest")
    output$fig_agg   <- render_fig("Fig_Pancan_aggregate_forest", "Fig Aggregate Forest")
    output$fig7      <- render_fig("Fig7_pancan_violins", "Fig 7: Pan-Cancer Violins")

    # Table 3.5 DataTable
    output$pancan_table <- DT::renderDataTable({
      tbl_data <- if (!is.null(bundle$pancan_table_tbl)) bundle$pancan_table_tbl else bundle$pancan_stats_df
      if (is.null(tbl_data)) {
        return(DT::datatable(data.frame(Message = "Pan-cancer survival table not available in bundle"), rownames = FALSE))
      }

      df <- as.data.frame(tbl_data)
      num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
      for (col in num_cols) {
        if (grepl("P|FDR", col, ignore.case = TRUE)) {
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
    })

    # Treated Table DataTable
    output$treated_table <- DT::renderDataTable({
      if (is.null(bundle) || is.null(bundle$treated_stats_df)) {
        return(DT::datatable(data.frame(Message = "Pan-cancer treated statistics not available in bundle"), rownames = FALSE))
      }

      df <- as.data.frame(bundle$treated_stats_df)
      num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
      for (col in num_cols) {
        if (grepl("P|FDR", col, ignore.case = TRUE)) {
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
    })
  })
}
