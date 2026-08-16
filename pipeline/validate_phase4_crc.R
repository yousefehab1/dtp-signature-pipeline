#!/usr/bin/env Rscript
# ==============================================================================
# pipeline/validate_phase4_crc.R - End-to-end validation of Phase 4 CRC presentation
# layer (7 composite figures, Table 3.3/3.3C persistence, Excel report, reference diff).
# ==============================================================================

suppressPackageStartupMessages({
  library(devtools)
  library(openxlsx)
  library(readr)
  library(dplyr)
})

cat("======================================================================\n")
cat("Phase 4 CRC Validation Pipeline\n")
cat("======================================================================\n\n")

# 1. Load package and configurations
devtools::load_all(".")
cfg <- dtp_config()
panel_tbl <- load_signature_panel(cfg)
composite_defs <- load_composite_defs(cfg)

out_root <- "output/phase4_validate"
fig_dir  <- file.path(out_root, "figures")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Run CRC Survival Analysis (Phase 3)
cat("\n--- Step 2: Running run_crc_survival() ---\n")
t0 <- Sys.time()
crc_result <- run_crc_survival(cfg, panel_tbl, out_root = out_root)
t1 <- Sys.time()
cat(sprintf("run_crc_survival completed in %.2f seconds.\n", as.numeric(difftime(t1, t0, units = "secs"))))

# 3. Build and Save 7 CRC Composites
cat("\n--- Step 3: Building 7 CRC Composites ---\n")
comp_res <- build_crc_composites(crc_result, panel_tbl, composite_defs, cfg, out_dir = fig_dir)

expected_figs <- c(
  "Fig1_outcome_composite",
  "Fig2_km_composite",
  "Fig3_subgroup_composite",
  "Fig3_3A_subgroup_hr_matrix",
  "Fig3_3B_score_across_subtype",
  "Fig3_3C_subtype_violins",
  "Fig5_confounding_summary"
)

fig_files <- list.files(fig_dir, full.names = FALSE)
cat(sprintf("Found %d files in %s\n", length(fig_files), fig_dir))

all_figs_ok <- TRUE
for (f in expected_figs) {
  has_pdf <- paste0(f, ".pdf") %in% fig_files
  has_png <- paste0(f, ".png") %in% fig_files
  cat(sprintf("  %s: PDF=%s, PNG=%s\n", f, ifelse(has_pdf, "OK", "MISSING"), ifelse(has_png, "OK", "MISSING")))
  if (!has_pdf || !has_png) all_figs_ok <- FALSE
}
cat(sprintf("All 7 PDF+PNG pairs present: %s (14 files total)\n", ifelse(all_figs_ok, "YES", "NO")))

# 4. Build Excel Report Workbook
cat("\n--- Step 4: Building Excel Workbook ---\n")
excel_path <- file.path(out_root, "DTP_Results.xlsx")
crc_for_excel <- list(
  stats_df             = crc_result$stats_df,
  subtype_survival_tbl = comp_res$subtype_survival_tbl,
  adj_df               = crc_result$adj_df,
  int_df               = crc_result$int_df,
  subtype_stats        = crc_result$subtype_stats,
  subtype_pairwise_tbl = comp_res$subtype_pairwise_tbl
)
build_excel_workbook(cfg, excel_path, panel_tbl, composite_defs, crc = crc_for_excel)

if (file.exists(excel_path)) {
  sheet_names <- openxlsx::getSheetNames(excel_path)
  cat("Excel workbook created successfully. Sheets:\n")
  for (i in seq_along(sheet_names)) {
    cat(sprintf("  [%d] %s\n", i, sheet_names[i]))
  }
} else {
  stop("Failed to create Excel workbook at ", excel_path)
}

# 5. Persist extracted tables to CSV
cat("\n--- Step 5: Persisting extracted tables to CSV ---\n")
sub_surv_csv <- file.path(out_root, "subtype_survival.csv")
sub_pair_csv <- file.path(out_root, "subtype_pairwise.csv")
readr::write_csv(comp_res$subtype_survival_tbl, sub_surv_csv)
readr::write_csv(comp_res$subtype_pairwise_tbl, sub_pair_csv)
cat("Wrote:\n  ", sub_surv_csv, "\n  ", sub_pair_csv, "\n")

# 6. Table Dimensions Summary
cat("\n--- Step 6: Table Dimensions Summary ---\n")
tables <- list(
  stats_df             = crc_result$stats_df,
  adj_df               = crc_result$adj_df,
  int_df               = crc_result$int_df,
  subtype_stats        = crc_result$subtype_stats,
  subtype_survival_tbl = comp_res$subtype_survival_tbl,
  subtype_pairwise_tbl = comp_res$subtype_pairwise_tbl
)

for (nm in names(tables)) {
  d <- dim(tables[[nm]])
  cat(sprintf("  %-22s: %4d rows x %2d cols\n", nm, d[1], d[2]))
}

# 7. Numerical comparison against Thesis reference outputs
cat("\n--- Step 7: Numerical Comparison Against Reference ---\n")
ref_dir <- "/Users/elabd/Master's Work/Thesis"
ref_surv_file <- file.path(ref_dir, "CRC_DTP_20260815_1705", "crc_survival", "composites", "Table_3_3_subtype_survival.csv")
ref_pair_file <- file.path(ref_dir, "tables", "AppendixTableA7_subtype_pairwise.csv")

if (file.exists(ref_surv_file)) {
  ref_surv <- readr::read_csv(ref_surv_file, show_col_types = FALSE)
  new_surv <- comp_res$subtype_survival_tbl
  cat(sprintf("Table 3.3 (Subgroup Survival) comparison:\n  Ref rows: %d, New rows: %d\n", nrow(ref_surv), nrow(new_surv)))
  
  # Check matching key columns
  key_cols <- c("Dataset", "Subgroup", "Score", "Endpoint")
  merged_surv <- dplyr::inner_join(new_surv, ref_surv, by = key_cols, suffix = c(".new", ".ref"))
  cat(sprintf("  Matched rows on (%s): %d\n", paste(key_cols, collapse = ", "), nrow(merged_surv)))
  
  # Compare numeric columns: HR_perSD, HR_lower, HR_upper, Effect_r, Cox_FDR_P, C_index
  num_cols <- c("HR_perSD", "HR_lower", "HR_upper", "Effect_r", "Cox_FDR_P", "C_index")
  for (col in num_cols) {
    c_new <- merged_surv[[paste0(col, ".new")]]
    c_ref <- merged_surv[[paste0(col, ".ref")]]
    ok <- !is.na(c_new) & !is.na(c_ref)
    if (sum(ok) > 0) {
      rel_diff <- abs(c_new[ok] - c_ref[ok]) / pmax(abs(c_ref[ok]), 1e-6)
      max_rd <- max(rel_diff)
      mean_rd <- mean(rel_diff)
      cat(sprintf("    %-12s: max rel diff = %.4f, mean rel diff = %.4f (n=%d evaluated)\n",
                  col, max_rd, mean_rd, sum(ok)))
    }
  }
} else {
  cat("  Reference Table_3_3_subtype_survival.csv not found at: ", ref_surv_file, "\n")
}

if (file.exists(ref_pair_file)) {
  ref_pair <- readr::read_csv(ref_pair_file, show_col_types = FALSE)
  new_pair <- comp_res$subtype_pairwise_tbl
  cat(sprintf("\nTable 3.3C (Subtype Pairwise) comparison:\n  Ref rows: %d, New rows: %d\n", nrow(ref_pair), nrow(new_pair)))
  
  key_cols_p <- c("Dataset", "Subtype_Axis", "Score", "Group_A", "Group_B")
  merged_pair <- dplyr::inner_join(new_pair, ref_pair, by = key_cols_p, suffix = c(".new", ".ref"))
  cat(sprintf("  Matched rows on (%s): %d\n", paste(key_cols_p, collapse = ", "), nrow(merged_pair)))
  
  num_cols_p <- c("Median_A", "Median_B", "Effect_r", "Raw_P", "FDR_P")
  for (col in num_cols_p) {
    c_new <- merged_pair[[paste0(col, ".new")]]
    c_ref <- merged_pair[[paste0(col, ".ref")]]
    ok <- !is.na(c_new) & !is.na(c_ref)
    if (sum(ok) > 0) {
      rel_diff <- abs(c_new[ok] - c_ref[ok]) / pmax(abs(c_ref[ok]), 1e-6)
      max_rd <- max(rel_diff)
      mean_rd <- mean(rel_diff)
      cat(sprintf("    %-12s: max rel diff = %.4f, mean rel diff = %.4f (n=%d evaluated)\n",
                  col, max_rd, mean_rd, sum(ok)))
    }
  }
} else {
  cat("  Reference AppendixTableA7_subtype_pairwise.csv not found at: ", ref_pair_file, "\n")
}

cat("\n======================================================================\n")
cat("Phase 4 CRC Validation Complete!\n")
cat("======================================================================\n")
