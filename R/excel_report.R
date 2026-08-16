# ==============================================================================
# R/excel_report.R - Excel Report Workbook Generation (openxlsx)
# ==============================================================================

#' Build the comprehensive DTP analysis Excel workbook
#'
#' Generates an 8-sheet Excel workbook with run metadata, table of contents,
#' canonical signature definitions, and all CRC survival and subtyping statistics tables.
#' Includes bold frozen headers, autofilters, auto-sized columns, and conditional
#' formatting highlighting FDR < 0.05.
#'
#' @param cfg Config from [dtp_config()].
#' @param out_path Path to write the `.xlsx` file.
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param crc Named list of CRC survival analysis results containing `stats_df`,
#'   `adj_df`, `int_df`, `subtype_stats`, `subtype_survival_tbl`, and `subtype_pairwise_tbl`.
#' @export
build_excel_workbook <- function(cfg, out_path, panel_tbl, composite_defs, crc = list()) {
  wb <- openxlsx::createWorkbook()

  header_style    <- openxlsx::createStyle(textDecoration = "bold")
  title_style     <- openxlsx::createStyle(fontSize = 14, textDecoration = "bold")
  section_style   <- openxlsx::createStyle(fontSize = 11, textDecoration = "bold")
  highlight_style <- openxlsx::createStyle(bgFill = "#FFC7CE", fontColour = "#9C0006")

  # ----------------------------------------------------------------------------
  # Sheet 1: README
  # ----------------------------------------------------------------------------
  openxlsx::addWorksheet(wb, "README")

  pkg_version <- tryCatch(
    as.character(utils::packageVersion("dtpsig")),
    error = function(e) "0.1.0"
  )

  readme_meta <- data.frame(
    Property = c("Generated At", "Gene ID Type", "dtpsig Version", "Score Suffix"),
    Value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      as.character(if (!is.null(cfg$id_type)) cfg$id_type else "ensembl"),
      pkg_version,
      as.character(if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA")
    ),
    stringsAsFactors = FALSE
  )

  readme_toc <- data.frame(
    Sheet_Name = c(
      "Signature_Panel",
      "Composite_Definitions",
      "CRC_Univariable_Survival",
      "CRC_Subgroup_Survival",
      "CRC_Adjusted_Cox",
      "CRC_Interaction_Cox",
      "CRC_Subtype_Characterization"
    ),
    Description = c(
      "Canonical gene signature definitions and member genes (panel.csv)",
      "Composite score definitions and weights (composite_defs.csv)",
      "Univariable survival analysis statistics across CRC cohorts (Wilcoxon, KM, Cox)",
      "Univariable survival analysis within clinical and molecular subgroups (Table 3.3)",
      "Multivariable Cox proportional hazards models adjusted for clinical/molecular covariates",
      "Cox interaction tests for effect modification across subgroup modifiers",
      "Kruskal-Wallis tests characterizing signature score variation across molecular subtypes"
    ),
    stringsAsFactors = FALSE
  )

  openxlsx::writeData(wb, "README", "DTP Signature Pipeline - Analysis Report", startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, "README", style = title_style, rows = 1, cols = 1)

  openxlsx::writeData(wb, "README", "Run Metadata", startRow = 3, startCol = 1)
  openxlsx::addStyle(wb, "README", style = section_style, rows = 3, cols = 1)
  openxlsx::writeData(wb, "README", readme_meta, startRow = 4, startCol = 1)
  openxlsx::addStyle(wb, "README", style = header_style, rows = 4, cols = 1:2)

  openxlsx::writeData(wb, "README", "Table of Contents", startRow = 10, startCol = 1)
  openxlsx::addStyle(wb, "README", style = section_style, rows = 10, cols = 1)
  openxlsx::writeData(wb, "README", readme_toc, startRow = 11, startCol = 1)
  openxlsx::addStyle(wb, "README", style = header_style, rows = 11, cols = 1:2)

  openxlsx::setColWidths(wb, "README", cols = 1:2, widths = "auto")

  # ----------------------------------------------------------------------------
  # Helper to format and add data sheets
  # ----------------------------------------------------------------------------
  add_data_sheet <- function(sheet_name, data_tbl) {
    df <- if (is.null(data_tbl)) data.frame() else as.data.frame(data_tbl)
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, df, startRow = 1, startCol = 1)

    if (ncol(df) > 0) {
      openxlsx::addStyle(
        wb, sheet_name, style = header_style,
        rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE
      )
      openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
      openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")

      if (nrow(df) > 0) {
        openxlsx::addFilter(wb, sheet_name, row = 1, cols = seq_len(ncol(df)))

        fdr_col_idx <- which(colnames(df) == "FDR_P")
        if (length(fdr_col_idx) == 1) {
          col_letter <- openxlsx::int2col(fdr_col_idx)
          rule_str <- sprintf("$%s2 < 0.05", col_letter)
          openxlsx::conditionalFormatting(
            wb, sheet_name,
            cols = seq_len(ncol(df)),
            rows = 2:(nrow(df) + 1),
            rule = rule_str,
            style = highlight_style
          )
        }
      }
    }
  }

  # ----------------------------------------------------------------------------
  # Sheets 2-8 in exact requested order
  # ----------------------------------------------------------------------------
  add_data_sheet("Signature_Panel", panel_tbl)
  add_data_sheet("Composite_Definitions", composite_defs)
  add_data_sheet("CRC_Univariable_Survival", crc$stats_df)
  add_data_sheet("CRC_Subgroup_Survival", crc$subtype_survival_tbl)
  add_data_sheet("CRC_Adjusted_Cox", crc$adj_df)
  add_data_sheet("CRC_Interaction_Cox", crc$int_df)
  add_data_sheet("CRC_Subtype_Characterization", crc$subtype_stats)

  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  message("  wrote Excel workbook -> ", out_path)
  invisible(out_path)
}
