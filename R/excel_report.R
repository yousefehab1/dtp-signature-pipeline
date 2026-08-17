# ==============================================================================
# R/excel_report.R - Excel Report Workbook Generation (openxlsx)
# ==============================================================================

#' Build the comprehensive DTP analysis Excel workbook
#'
#' Generates an Excel workbook with run metadata, table of contents,
#' canonical signature definitions, CRC survival/subtyping statistics tables,
#' and pan-cancer survival tables.
#' Includes bold frozen headers, autofilters, auto-sized columns, and conditional
#' formatting highlighting FDR < 0.05.
#'
#' @param cfg Config from [dtp_config()].
#' @param out_path Path to write the `.xlsx` file.
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param crc Named list of CRC survival analysis results containing `stats_df`,
#'   `adj_df`, `int_df`, `subtype_stats`, `subtype_survival_tbl`, and `subtype_pairwise_tbl`.
#' @param pancan Named list of pan-cancer survival analysis results containing
#'   `pancan_table_tbl` (from [build_pancan_composites()]) and `treated_stats_df`
#'   (from [run_pancan_treated()]).
#' @param mets Named list of metastasis DE/GSEA analysis results containing
#'   `gsea_pn`, `gsea_pm_adj`, `pvm_cmp`, `liver_df`, `roc_summary`,
#'   `scores_df`, and `ssgsea_pvalues` (from [run_mets_de()]).
#' @export
build_excel_workbook <- function(cfg, out_path, panel_tbl, composite_defs,
                                 crc = list(), pancan = list(), mets = list()) {
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

  if (!is.null(pancan$pancan_table_tbl) || !is.null(pancan$treated_stats_df)) {
    pancan_toc <- data.frame(
      Sheet_Name = c(
        "PanCancer_Survival",
        "PanCancer_Treated"
      ),
      Description = c(
        "Per-cohort and aggregate continuous-score Cox survival analysis across TCGA cohorts (Table 3.5)",
        "Univariable survival analysis within treated patient subsets across TCGA cohorts"
      ),
      stringsAsFactors = FALSE
    )
    readme_toc <- rbind(readme_toc, pancan_toc)
  }

  has_mets <- !is.null(mets$gsea_pn) || !is.null(mets$gsea_pm_adj) || !is.null(mets$pvm_cmp) ||
              !is.null(mets$liver_df) || !is.null(mets$roc_summary) || !is.null(mets$scores_df) ||
              !is.null(mets$ssgsea_pvalues)
  if (has_mets) {
    mets_toc <- data.frame(
      Sheet_Name = c(
        "Mets_GSEA_NormalVsPrimary",
        "Mets_GSEA_PvM_Adjusted",
        "Mets_GSEA_Adj_vs_Unadj",
        "Mets_LiverPurity_Scores",
        "Mets_ROC_AUC",
        "Mets_SingleSample_Scores",
        "Mets_SingleSample_Wilcoxon"
      ),
      Description = c(
        "Gene Set Enrichment Analysis summary for Primary CRC vs Normal mucosa (GSE50760)",
        "Gene Set Enrichment Analysis summary for Primary CRC vs Liver Metastasis purity-adjusted for liver score (GSE50760)",
        "Comparison of GSEA NES and FDR between liver-unadjusted and purity-adjusted Primary vs Metastasis models",
        "Per-sample liver-specific gene ssGSEA validation scores across Normal, Primary, and Metastasis tissues",
        "Single-sample paired ROC AUC distinguishing Primary CRC from Liver Metastasis",
        "Single-sample paired ssGSEA scores on liver-corrected expression matrix",
        "Paired Wilcoxon test statistics and BH-adjusted p-values for single-sample ssGSEA scores"
      ),
      stringsAsFactors = FALSE
    )
    readme_toc <- rbind(readme_toc, mets_toc)
  }

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
  # Sheets in exact requested order
  # ----------------------------------------------------------------------------
  add_data_sheet("Signature_Panel", panel_tbl)
  add_data_sheet("Composite_Definitions", composite_defs)
  add_data_sheet("CRC_Univariable_Survival", crc$stats_df)
  add_data_sheet("CRC_Subgroup_Survival", crc$subtype_survival_tbl)
  add_data_sheet("CRC_Adjusted_Cox", crc$adj_df)
  add_data_sheet("CRC_Interaction_Cox", crc$int_df)
  add_data_sheet("CRC_Subtype_Characterization", crc$subtype_stats)
  if (!is.null(pancan$pancan_table_tbl) || !is.null(pancan$treated_stats_df)) {
    add_data_sheet("PanCancer_Survival", pancan$pancan_table_tbl)
    add_data_sheet("PanCancer_Treated", pancan$treated_stats_df)
  }
  if (has_mets) {
    add_data_sheet("Mets_GSEA_NormalVsPrimary", mets$gsea_pn)
    add_data_sheet("Mets_GSEA_PvM_Adjusted", mets$gsea_pm_adj)
    add_data_sheet("Mets_GSEA_Adj_vs_Unadj", mets$pvm_cmp)
    add_data_sheet("Mets_LiverPurity_Scores", mets$liver_df)
    add_data_sheet("Mets_ROC_AUC", mets$roc_summary)
    add_data_sheet("Mets_SingleSample_Scores", mets$scores_df)
    add_data_sheet("Mets_SingleSample_Wilcoxon", mets$ssgsea_pvalues)
  }

  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  message("  wrote Excel workbook -> ", out_path)
  invisible(out_path)
}
