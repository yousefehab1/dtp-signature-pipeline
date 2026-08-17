# ==============================================================================
# dashboard/R/mod_crc.R - CRC Survival Presentation Module
#
# Displays 7 composite figures and the CRC per-test statistics table.
# ==============================================================================

#' CRC Survival UI Module
#'
#' @param id Module namespace identifier.
#' @return A Shiny UI tag list.
#' @export
mod_crc_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::uiOutput(ns("crc_alert")),

    bslib::navset_card_tab(
      title = "CRC Composite Figures",
      bslib::nav_panel(
        title = "Fig 1: Outcome Association",
        shiny::imageOutput(ns("fig1"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 2: Kaplan-Meier Curves",
        shiny::imageOutput(ns("fig2"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 3: Subgroup Forest Plots",
        shiny::imageOutput(ns("fig3"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 3.3A: Hazard Ratio Matrix",
        shiny::imageOutput(ns("fig3_3a"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 3.3B: Subtype Score Variation",
        shiny::imageOutput(ns("fig3_3b"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 3.3C: Subtype Violins",
        shiny::imageOutput(ns("fig3_3c"), height = "auto")
      ),
      bslib::nav_panel(
        title = "Fig 5: Confounding Summary",
        shiny::imageOutput(ns("fig5"), height = "auto")
      )
    ),

    bslib::card(
      bslib::card_header(shiny::strong("CRC Univariable Survival Statistics Table")),
      DT::dataTableOutput(ns("crc_stats_table")),
      bslib::card_footer(
        class = "text-muted small",
        "Precomputed univariate survival statistics across cohorts (GSE39582, TCGA-COAD) and tests (Wilcoxon, KM, Cox) with Benjamini-Hochberg FDR correction."
      )
    )
  )
}

#' CRC Survival Server Module
#'
#' @param id Module namespace identifier.
#' @param bundle Precomputed results bundle loaded from `results_bundle.rds`.
#' @export
mod_crc_server <- function(id, bundle) {
  shiny::moduleServer(id, function(input, output, session) {

    # Guard for missing bundle or CRC data
    output$crc_alert <- shiny::renderUI({
      if (is.null(bundle)) {
        bslib::card(
          class = "bg-danger text-white mb-3",
          bslib::card_body(
            "⚠️ Results bundle missing! Please run `pipeline/build_results_bundle.R` to generate precomputed data."
          )
        )
      } else if (is.null(bundle$stats_df)) {
        bslib::card(
          class = "bg-warning text-dark mb-3",
          bslib::card_body(
            "⚠️ CRC survival statistics not found in this bundle. Please re-run `pipeline/build_results_bundle.R`."
          )
        )
      } else {
        NULL
      }
    })

    # Helper function to render on-disk PNGs or generate a graceful placeholder
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

    # Render 7 CRC figures
    output$fig1    <- render_fig("Fig1_outcome_composite", "Fig 1: 3-Year Outcome Association")
    output$fig2    <- render_fig("Fig2_km_composite", "Fig 2: Kaplan-Meier Survival Curves")
    output$fig3    <- render_fig("Fig3_subgroup_composite", "Fig 3: Subgroup Forest Plots")
    output$fig3_3a <- render_fig("Fig3_3A_subgroup_hr_matrix", "Fig 3.3A: Subgroup HR Matrix")
    output$fig3_3b <- render_fig("Fig3_3B_score_across_subtype", "Fig 3.3B: Subtype Score Variation")
    output$fig3_3c <- render_fig("Fig3_3C_subtype_violins", "Fig 3.3C: Subtype Violins")
    output$fig5    <- render_fig("Fig5_confounding_summary", "Fig 5: Confounding Summary")

    # Render DataTable
    output$crc_stats_table <- DT::renderDataTable({
      if (is.null(bundle) || is.null(bundle$stats_df)) {
        return(DT::datatable(data.frame(Message = "No CRC survival statistics available"), rownames = FALSE))
      }

      df <- as.data.frame(bundle$stats_df)

      # Round numeric columns if present
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
  })
}
