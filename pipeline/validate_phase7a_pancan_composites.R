#!/usr/bin/env Rscript
# ==============================================================================
# pipeline/validate_phase7a_pancan_composites.R - End-to-end validation of
# Pan-Cancer presentation layer (6 composite figures, Table 3.5, Excel report,
# reference diff, and live captions).
# ==============================================================================

suppressPackageStartupMessages({
  library(devtools)
  library(openxlsx)
  library(readr)
  library(dplyr)
})

cat("======================================================================\n")
cat("Phase 7a Pan-Cancer Composites Validation Pipeline\n")
cat("======================================================================\n\n")

# 1. Load package and configurations
devtools::load_all(".")
cfg <- dtp_config()
panel_tbl <- load_signature_panel(cfg)
composite_defs <- load_composite_defs(cfg)

out_root <- "output/phase7a_validate"
fig_dir  <- file.path(out_root, "figures")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Run Pan-Cancer Survival Analysis (All Patients + Treated)
cat("\n--- Step 2: Running run_pancan_survival() + run_pancan_treated() ---\n")
t0 <- Sys.time()
pancan_result <- run_pancan_survival(cfg, panel_tbl)
treated_result <- run_pancan_treated(cfg, pancan_result)
t1 <- Sys.time()
cat(sprintf("Pan-cancer survival and treated analysis completed in %.2f seconds.\n",
            as.numeric(difftime(t1, t0, units = "secs"))))

# 3. Build and Save 6 Pan-Cancer Composites
cat("\n--- Step 3: Building Pan-Cancer Composites ---\n")
comp_res <- build_pancan_composites(pancan_result, cfg, out_dir = fig_dir)

expected_figs <- c(
  "Fig6A_pancan_forest_Up",
  "Fig6B_pancan_forest_Composite",
  "Fig6C_pancan_forest_Down",
  "Fig6_pancan_survival_composite",
  "Fig_Pancan_aggregate_forest",
  "Fig7_pancan_violins"
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
cat(sprintf("All 6 PDF+PNG pairs present: %s (12 figure files total)\n",
            ifelse(all_figs_ok, "YES", "NO")))

# 4. Build Excel Report Workbook
cat("\n--- Step 4: Building Excel Workbook with PanCancer Sheets ---\n")
excel_path <- file.path(out_root, "DTP_PanCancer_Results.xlsx")
pancan_for_excel <- list(
  pancan_table_tbl = comp_res$pancan_table_tbl,
  treated_stats_df = treated_result$stats_df
)
build_excel_workbook(cfg, excel_path, panel_tbl, composite_defs, crc = list(), pancan = pancan_for_excel)

if (file.exists(excel_path)) {
  sheet_names <- openxlsx::getSheetNames(excel_path)
  cat("Excel workbook created successfully. Sheets:\n")
  for (i in seq_along(sheet_names)) {
    s_data <- openxlsx::read.xlsx(excel_path, sheet = sheet_names[i])
    cat(sprintf("  [%2d] %-30s (%d rows x %d cols)\n",
                i, sheet_names[i], nrow(s_data), ncol(s_data)))
  }
  has_pancan_surv <- "PanCancer_Survival" %in% sheet_names
  has_pancan_trt  <- "PanCancer_Treated" %in% sheet_names
  cat(sprintf("PanCancer_Survival present: %s\n", ifelse(has_pancan_surv, "YES", "NO")))
  cat(sprintf("PanCancer_Treated present:  %s\n", ifelse(has_pancan_trt, "YES", "NO")))
} else {
  stop("Failed to create Excel workbook at ", excel_path)
}

# 5. Persist Table 3.5 to CSV
cat("\n--- Step 5: Table 3.5 Dimensions & Output ---\n")
tbl_csv <- file.path(out_root, "Table_3_5_pancan_survival.csv")
readr::write_csv(comp_res$pancan_table_tbl, tbl_csv)
d_tbl <- dim(comp_res$pancan_table_tbl)
cat(sprintf("pancan_table_tbl: %d rows x %d cols -> %s\n", d_tbl[1], d_tbl[2], tbl_csv))

# 6. Numerical comparison of Table 3.5 against Thesis reference
cat("\n--- Step 6: Numerical Comparison Against Reference (Table 3.5) ---\n")
ref_tbl_file <- "/Users/elabd/Master's Work/Thesis/CRC_DTP_20260815_1705/pancan_survival/composites/Table_3_5_pancan_survival.csv"

if (file.exists(ref_tbl_file)) {
  ref_tbl <- readr::read_csv(ref_tbl_file, show_col_types = FALSE)
  new_tbl <- comp_res$pancan_table_tbl
  cat(sprintf("Table 3.5 comparison:\n  Ref rows: %d, New rows: %d\n", nrow(ref_tbl), nrow(new_tbl)))

  key_cols <- c("Project", "Score", "Endpoint")
  merged_tbl <- dplyr::inner_join(new_tbl, ref_tbl, by = key_cols, suffix = c(".new", ".ref"))
  cat(sprintf("  Matched rows on (%s): %d / %d\n",
              paste(key_cols, collapse = ", "), nrow(merged_tbl), nrow(ref_tbl)))

  num_cols <- c("HR_perSD", "HR_lower", "HR_upper", "C_index", "Cox_FDR_P")
  for (col in num_cols) {
    c_new <- merged_tbl[[paste0(col, ".new")]]
    c_ref <- merged_tbl[[paste0(col, ".ref")]]
    ok <- !is.na(c_new) & !is.na(c_ref)
    if (sum(ok) > 0) {
      rel_diff <- abs(c_new[ok] - c_ref[ok]) / pmax(abs(c_ref[ok]), 1e-6)
      max_rd <- max(rel_diff)
      mean_rd <- mean(rel_diff)
      cat(sprintf("    %-12s: max rel diff = %.6f, mean rel diff = %.6f (n=%d evaluated)\n",
                  col, max_rd, mean_rd, sum(ok)))
    }
  }
} else {
  cat("  Reference Table_3_5_pancan_survival.csv not found at: ", ref_tbl_file, "\n")
}

# 7. Print Caption Files
cat("\n--- Step 7: Figure Captions ---\n")
cap6_file <- file.path(fig_dir, "Fig6_pancan_survival_caption.txt")
cap7_file <- file.path(fig_dir, "Fig7_pancan_violins_caption.txt")

if (file.exists(cap6_file)) {
  cat("\n=================== Figure 6 Caption ===================\n")
  cat(readr::read_file(cap6_file))
  cat("========================================================\n")
}

if (file.exists(cap7_file)) {
  cat("\n=================== Figure 7 Caption ===================\n")
  cat(readr::read_file(cap7_file))
  cat("========================================================\n")
}

cat("\n======================================================================\n")
cat("Phase 7a Pan-Cancer Composites Validation Complete!\n")
cat("======================================================================\n")
