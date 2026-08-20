#!/usr/bin/env Rscript
# ==============================================================================
# pipeline/validate_phase7b_mets_composites.R - End-to-end validation of
# Metastasis presentation layer (2 composite figures, 7 Excel sheets,
# and warm cache verification).
# ==============================================================================

suppressPackageStartupMessages({
  library(devtools)
  library(openxlsx)
  library(readr)
  library(dplyr)
})

cat("======================================================================\n")
cat("Phase 7b Metastasis Composites Validation Pipeline\n")
cat("======================================================================\n\n")

# 1. Load package and configurations
devtools::load_all(".")
cfg <- dtp_config()
panel_tbl <- load_signature_panel(cfg)
composite_defs <- load_composite_defs(cfg)

out_root <- "output/phase7b_validate"
fig_dir  <- file.path(out_root, "figures")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Run Metastasis DE & GSEA Analysis against warm cache
cat("\n--- Step 2: Running run_mets_de() ---\n")
t0 <- Sys.time()
cache_logs <- character()
mets_result <- withCallingHandlers(
  run_mets_de(cfg, panel_tbl),
  message = function(m) {
    txt <- conditionMessage(m)
    if (grepl("^\\[cache\\]", txt)) cache_logs <<- c(cache_logs, txt)
  }
)
t1 <- Sys.time()
cat(sprintf("run_mets_de completed in %.2f seconds.\n",
            as.numeric(difftime(t1, t0, units = "secs"))))

# 3. Cache Verification (Zero fresh GEO downloads)
cat("\n--- Step 3: Cache Verification (Zero Fresh GEO Downloads) ---\n")
for (cl in unique(cache_logs)) cat(" ", cl, "\n")
fresh_builds <- grep("\\(build\\)", cache_logs, value = TRUE)
if (length(fresh_builds) == 0) {
  cat("  [PASS] Zero fresh GEO downloads triggered (all loaded from cache).\n")
} else {
  cat("  [!] Fresh cache builds detected:\n")
  for (fb in fresh_builds) cat("   ", fb, "\n")
}

# 4. Build and Save 2 Metastasis Composites
cat("\n--- Step 4: Building Metastasis Composites ---\n")
comp_res <- build_mets_composites(mets_result, cfg, out_dir = fig_dir)

expected_figs <- c(
  "Fig_Mets_GSEA_composite",
  "Fig_Mets_LiverPurity_composite"
)

fig_files <- list.files(fig_dir, full.names = FALSE)
cat(sprintf("Found %d files in %s\n", length(fig_files), fig_dir))

all_figs_ok <- TRUE
for (f in expected_figs) {
  has_pdf <- paste0(f, ".pdf") %in% fig_files
  has_png <- paste0(f, ".png") %in% fig_files
  cat(sprintf("  %-32s: PDF=%s, PNG=%s\n", f,
              ifelse(has_pdf, "OK", "MISSING"), ifelse(has_png, "OK", "MISSING")))
  if (!has_pdf || !has_png) all_figs_ok <- FALSE
}
cat(sprintf("All figure PDF+PNG pairs present: %s (4 figure files total)\n",
            ifelse(all_figs_ok, "YES", "NO")))

# 5. Build Excel Report Workbook
cat("\n--- Step 5: Building Excel Workbook with Metastasis Sheets ---\n")
excel_path <- file.path(out_root, "DTP_Mets_Results.xlsx")
mets_for_excel <- list(
  gsea_pn        = mets_result$gsea_pn,
  gsea_pm_adj    = mets_result$gsea_pm_adj,
  pvm_cmp        = mets_result$pvm_cmp,
  liver_df       = mets_result$liver_df,
  roc_summary    = mets_result$roc_summary,
  scores_df      = mets_result$scores_df,
  ssgsea_pvalues = mets_result$ssgsea_pvalues
)
build_excel_workbook(cfg, excel_path, panel_tbl, composite_defs, crc = list(), pancan = list(), mets = mets_for_excel)

expected_sheets <- c(
  "Mets_GSEA_NormalVsPrimary",
  "Mets_GSEA_PvM_Adjusted",
  "Mets_GSEA_Adj_vs_Unadj",
  "Mets_LiverPurity_Scores",
  "Mets_ROC_AUC",
  "Mets_SingleSample_Scores",
  "Mets_SingleSample_Wilcoxon"
)

if (file.exists(excel_path)) {
  sheet_names <- openxlsx::getSheetNames(excel_path)
  cat("Excel workbook created successfully. Sheets:\n")
  for (i in seq_along(sheet_names)) {
    s_data <- openxlsx::read.xlsx(excel_path, sheet = sheet_names[i])
    cat(sprintf("  [%2d] %-30s (%d rows x %d cols)\n",
                i, sheet_names[i], nrow(s_data), ncol(s_data)))
  }
  all_sheets_ok <- all(expected_sheets %in% sheet_names)
  cat(sprintf("All 7 metastasis sheets present: %s\n", ifelse(all_sheets_ok, "YES", "NO")))
} else {
  stop("Failed to create Excel workbook at ", excel_path)
}

# 6. Spot-print rows from each metastasis sheet
cat("\n--- Step 6: Spot-Printing Sample Rows From Each Metastasis Sheet ---\n")
for (s in expected_sheets) {
  s_data <- openxlsx::read.xlsx(excel_path, sheet = s)
  cat(sprintf("\n=== Sheet: %s (first 3 rows) ===\n", s))
  print(utils::head(s_data, 3))
}

cat("\n======================================================================\n")
cat("Phase 7b Metastasis Composites Validation Complete!\n")
cat("======================================================================\n")
