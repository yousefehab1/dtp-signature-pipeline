#!/usr/bin/env Rscript
# ==============================================================================
# pipeline/build_results_bundle.R - End-to-end master orchestrator for the
# dtpsig pipeline: runs CRC survival, Pan-Cancer survival & treated, Metastasis
# DE/GSEA, generates all 15 composite figures, builds the comprehensive Excel
# report, and saves the consolidated results bundle for the Shiny dashboard.
# ==============================================================================

suppressPackageStartupMessages({
  library(pkgload)
  library(openxlsx)
  library(readr)
  library(dplyr)
})

cat("======================================================================\n")
cat("Building Full DTP Pipeline Results Bundle & Artifacts\n")
cat("======================================================================\n\n")

# 1. Load package and configuration
pkgload::load_all(".")
cfg <- dtp_config()
panel_tbl <- load_signature_panel(cfg)
composite_defs <- load_composite_defs(cfg)

out_root   <- "output"
fig_dir    <- file.path(out_root, "figures")
excel_path <- file.path(out_root, "DTP_Results.xlsx")
bundle_path <- file.path(out_root, "results_bundle.rds")

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

total_t0 <- Sys.time()

# 2. Run CRC Survival Analysis (GSE39582 & TCGA-COAD)
cat("\n--- Step 1: Running run_crc_survival() ---\n")
t0 <- Sys.time()
crc_result <- run_crc_survival(cfg, panel_tbl, out_root = out_root)
t1 <- Sys.time()
cat(sprintf("  -> run_crc_survival completed in %.2f seconds.\n", as.numeric(difftime(t1, t0, units = "secs"))))

# 3. Run Pan-Cancer Survival Analysis (25 TCGA Cohorts + Treated Subset)
cat("\n--- Step 2: Running run_pancan_survival() + run_pancan_treated() ---\n")
t0 <- Sys.time()
pancan_result <- run_pancan_survival(cfg, panel_tbl)
treated_result <- run_pancan_treated(cfg, pancan_result)
t1 <- Sys.time()
cat(sprintf("  -> Pan-cancer survival & treated analysis completed in %.2f seconds.\n", as.numeric(difftime(t1, t0, units = "secs"))))

# 4. Run Metastasis DE & GSEA Analysis (GSE50760)
cat("\n--- Step 3: Running run_mets_de() ---\n")
t0 <- Sys.time()
mets_result <- run_mets_de(cfg, panel_tbl)
t1 <- Sys.time()
cat(sprintf("  -> run_mets_de completed in %.2f seconds.\n", as.numeric(difftime(t1, t0, units = "secs"))))

# 5. Build and Save All 15 Composite Figures (7 CRC + 6 Pan-Cancer + 2 Mets)
cat("\n--- Step 4: Building All Composite Figures -> ", fig_dir, " ---\n")
t0 <- Sys.time()
crc_comp    <- build_crc_composites(crc_result, panel_tbl, composite_defs, cfg, out_dir = fig_dir)
pancan_comp <- build_pancan_composites(pancan_result, cfg, out_dir = fig_dir)
mets_comp   <- build_mets_composites(mets_result, cfg, out_dir = fig_dir)
t1 <- Sys.time()
cat(sprintf("  -> All composite figures generated in %.2f seconds.\n", as.numeric(difftime(t1, t0, units = "secs"))))

# Check figure counts
expected_figs <- c(
  # CRC (7)
  "Fig1_outcome_composite", "Fig2_km_composite", "Fig3_subgroup_composite",
  "Fig3_3A_subgroup_hr_matrix", "Fig3_3B_score_across_subtype",
  "Fig3_3C_subtype_violins", "Fig5_confounding_summary",
  # Pan-Cancer (6)
  "Fig6A_pancan_forest_Up", "Fig6B_pancan_forest_Composite", "Fig6C_pancan_forest_Down",
  "Fig6_pancan_survival_composite", "Fig_Pancan_aggregate_forest", "Fig7_pancan_violins",
  # Metastasis (2)
  "Fig_Mets_GSEA_composite", "Fig_Mets_LiverPurity_composite"
)

fig_files <- list.files(fig_dir, full.names = FALSE)
cat(sprintf("Found %d total files in %s\n", length(fig_files), fig_dir))
all_figs_ok <- TRUE
for (f in expected_figs) {
  has_pdf <- paste0(f, ".pdf") %in% fig_files
  has_png <- paste0(f, ".png") %in% fig_files
  if (!has_pdf || !has_png) {
    all_figs_ok <- FALSE
    cat(sprintf("  MISSING: %s (PDF=%s, PNG=%s)\n", f, has_pdf, has_png))
  }
}
cat(sprintf("All 15 composite PDF+PNG pairs verified: %s (30 figure files expected)\n", ifelse(all_figs_ok, "YES", "NO")))

# 6. Build Consolidated Multi-Sheet Excel Workbook
cat("\n--- Step 5: Building Comprehensive Excel Report Workbook -> ", excel_path, " ---\n")
crc_for_excel <- list(
  stats_df             = crc_result$stats_df,
  subtype_survival_tbl = crc_comp$subtype_survival_tbl,
  adj_df               = crc_result$adj_df,
  int_df               = crc_result$int_df,
  subtype_stats        = crc_result$subtype_stats,
  subtype_pairwise_tbl = crc_comp$subtype_pairwise_tbl
)

pancan_for_excel <- list(
  pancan_table_tbl = pancan_comp$pancan_table_tbl,
  treated_stats_df = treated_result$stats_df
)

mets_for_excel <- list(
  gsea_pn        = mets_result$gsea_pn,
  gsea_pm_adj    = mets_result$gsea_pm_adj,
  pvm_cmp        = mets_result$pvm_cmp,
  liver_df       = mets_result$liver_df,
  roc_summary    = mets_result$roc_summary,
  scores_df      = mets_result$scores_df,
  ssgsea_pvalues = mets_result$ssgsea_pvalues
)

build_excel_workbook(
  cfg            = cfg,
  out_path       = excel_path,
  panel_tbl      = panel_tbl,
  composite_defs = composite_defs,
  crc            = crc_for_excel,
  pancan         = pancan_for_excel,
  mets           = mets_for_excel
)

if (file.exists(excel_path)) {
  sheet_names <- openxlsx::getSheetNames(excel_path)
  cat(sprintf("Excel report built successfully with %d sheets:\n", length(sheet_names)))
  for (i in seq_along(sheet_names)) {
    cat(sprintf("  [%2d] %s\n", i, sheet_names[i]))
  }
}

# 7. Save Consolidated Results Bundle
cat("\n--- Step 6: Saving Precomputed Results Bundle -> ", bundle_path, " ---\n")
save_results_bundle(
  cfg              = cfg,
  crc_result       = crc_result,
  panel_tbl        = panel_tbl,
  composite_defs   = composite_defs,
  path             = bundle_path,
  pancan_result    = pancan_result,
  treated_result   = treated_result,
  mets_result      = mets_result,
  fig_dir          = fig_dir,
  excel_path       = excel_path,
  pancan_table_tbl = pancan_comp$pancan_table_tbl
)

file_bytes <- file.size(bundle_path)
file_size_kb <- file_bytes / 1024
file_size_mb <- file_size_kb / 1024
cat(sprintf("Results bundle saved successfully: %d bytes (%.2f KB / %.2f MB)\n\n", file_bytes, file_size_kb, file_size_mb))

bundle <- load_results_bundle(bundle_path)
cat("--- Bundle Top-Level Fields ---\n")
print(names(bundle))

total_t1 <- Sys.time()
cat(sprintf("\nTotal pipeline execution time: %.2f seconds (%.2f minutes).\n",
            as.numeric(difftime(total_t1, total_t0, units = "secs")),
            as.numeric(difftime(total_t1, total_t0, units = "mins"))))

cat("\n======================================================================\n")
cat("DTP Master Results Bundle & Artifacts Build Complete!\n")
cat("======================================================================\n")
