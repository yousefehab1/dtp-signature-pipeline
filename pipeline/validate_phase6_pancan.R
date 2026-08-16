#!/usr/bin/env Rscript
# ==============================================================================
# pipeline/validate_phase6_pancan.R - End-to-end validation of Phase 6 Pan-Cancer
# survival analysis module (all-patients and treated-only) against Thesis reference.
# ==============================================================================

suppressPackageStartupMessages({
  library(devtools)
  library(readr)
  library(dplyr)
})

cat("======================================================================\n")
cat("Phase 6 Pan-Cancer Validation Pipeline\n")
cat("======================================================================\n\n")

# 1. Load package and configuration
devtools::load_all(".")
cfg <- dtp_config()
panel_tbl <- load_signature_panel(cfg)

# 2. Run Pan-Cancer Survival Analysis (All Patients)
cat("\n--- Step 2: Running run_pancan_survival() ---\n")
t0 <- Sys.time()
pancan_result <- run_pancan_survival(cfg, panel_tbl)
t1 <- Sys.time()
pancan_duration <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("run_pancan_survival completed in %.2f seconds.\n", pancan_duration))

# 3. Run Pan-Cancer Survival Analysis (Treated Patients Only)
cat("\n--- Step 3: Running run_pancan_treated() ---\n")
t2 <- Sys.time()
treated_result <- run_pancan_treated(cfg, pancan_result)
t3 <- Sys.time()
treated_duration <- as.numeric(difftime(t3, t2, units = "secs"))
cat(sprintf("run_pancan_treated completed in %.2f seconds.\n", treated_duration))

# 4. Print dimensions of all returned tables
cat("\n--- Step 4: Table Dimensions Summary ---\n")
cat("pancan_result tables:\n")
for (nm in names(pancan_result)) {
  if (is.data.frame(pancan_result[[nm]])) {
    d <- dim(pancan_result[[nm]])
    cat(sprintf("  %-22s: %5d rows x %2d cols\n", nm, d[1], d[2]))
  }
}
cat("\ntreated_result tables:\n")
for (nm in names(treated_result)) {
  if (is.data.frame(treated_result[[nm]])) {
    d <- dim(treated_result[[nm]])
    cat(sprintf("  %-22s: %5d rows x %2d cols\n", nm, d[1], d[2]))
  }
}

# 5. Numerical validation of pancan_result$stats_df against reference
cat("\n--- Step 5: Numerical Validation of Pan-Cancer Survival (All Patients) ---\n")
ref_pancan_file <- "/Users/elabd/Master's Work/Thesis/CRC_DTP_20260815_1705/pancan_survival/FDR_Stats_Summary.csv"

if (file.exists(ref_pancan_file)) {
  ref_pancan <- readr::read_csv(ref_pancan_file, show_col_types = FALSE)
  
  # Target score columns present in the reference run
  target_scores <- c(
    "Up_ssGSEA", "Down_ssGSEA", "Composite_ssGSEA",
    "Up_ssGSEA_Corrected", "Down_ssGSEA_Corrected", "Composite_ssGSEA_Corrected"
  )
  new_pancan_subset <- pancan_result$stats_df %>% dplyr::filter(Score %in% target_scores)
  
  cat(sprintf("Reference rows: %d, New total rows: %d, New shared-score rows: %d\n",
              nrow(ref_pancan), nrow(pancan_result$stats_df), nrow(new_pancan_subset)))
  
  key_cols <- c("Project", "Test", "Metric", "Score")
  merged_pancan <- dplyr::inner_join(new_pancan_subset, ref_pancan, by = key_cols, suffix = c(".new", ".ref"))
  cat(sprintf("Matched rows on (%s): %d / %d reference rows\n",
              paste(key_cols, collapse = ", "), nrow(merged_pancan), nrow(ref_pancan)))
  
  num_cols <- c("Raw_P", "Effect_r", "HR", "HR_lower", "HR_upper", "C_index", "FDR_P")
  for (col in num_cols) {
    c_new <- merged_pancan[[paste0(col, ".new")]]
    c_ref <- merged_pancan[[paste0(col, ".ref")]]
    ok <- !is.na(c_new) & !is.na(c_ref)
    if (sum(ok) > 0) {
      rel_diff <- abs(c_new[ok] - c_ref[ok]) / pmax(abs(c_ref[ok]), 1e-6)
      max_rd <- max(rel_diff)
      mean_rd <- mean(rel_diff)
      cat(sprintf("  %-12s: max rel diff = %.6f, mean rel diff = %.6f (n=%d evaluated)\n",
                  col, max_rd, mean_rd, sum(ok)))
    }
  }
} else {
  cat("  Reference FDR_Stats_Summary.csv not found at: ", ref_pancan_file, "\n")
}

# 6. Numerical validation of treated_result$stats_df against reference
cat("\n--- Step 6: Numerical Validation of Pan-Cancer Survival (Treated Only) ---\n")
ref_treated_file <- "/Users/elabd/Master's Work/Thesis/CRC_DTP_20260815_1705/pancan_survival/treated/FDR_Stats_Summary_Treated.csv"

if (file.exists(ref_treated_file)) {
  ref_treated <- readr::read_csv(ref_treated_file, show_col_types = FALSE)
  
  target_scores <- c(
    "Up_ssGSEA", "Down_ssGSEA", "Composite_ssGSEA",
    "Up_ssGSEA_Corrected", "Down_ssGSEA_Corrected", "Composite_ssGSEA_Corrected"
  )
  new_treated_subset <- treated_result$stats_df %>% dplyr::filter(Score %in% target_scores)
  
  cat(sprintf("Reference rows: %d, New total rows: %d, New shared-score rows: %d\n",
              nrow(ref_treated), nrow(treated_result$stats_df), nrow(new_treated_subset)))
  
  key_cols <- c("Project", "Test", "Metric", "Score")
  merged_treated <- dplyr::inner_join(new_treated_subset, ref_treated, by = key_cols, suffix = c(".new", ".ref"))
  cat(sprintf("Matched rows on (%s): %d / %d reference rows\n",
              paste(key_cols, collapse = ", "), nrow(merged_treated), nrow(ref_treated)))
  
  num_cols <- c("Raw_P", "Effect_r", "HR", "HR_lower", "HR_upper", "C_index", "FDR_P")
  for (col in num_cols) {
    c_new <- merged_treated[[paste0(col, ".new")]]
    c_ref <- merged_treated[[paste0(col, ".ref")]]
    ok <- !is.na(c_new) & !is.na(c_ref)
    if (sum(ok) > 0) {
      rel_diff <- abs(c_new[ok] - c_ref[ok]) / pmax(abs(c_ref[ok]), 1e-6)
      max_rd <- max(rel_diff)
      mean_rd <- mean(rel_diff)
      cat(sprintf("  %-12s: max rel diff = %.6f, mean rel diff = %.6f (n=%d evaluated)\n",
                  col, max_rd, mean_rd, sum(ok)))
    }
  }
} else {
  cat("  Reference FDR_Stats_Summary_Treated.csv not found at: ", ref_treated_file, "\n")
}

cat("\n======================================================================\n")
cat("Phase 6 Pan-Cancer Validation Complete!\n")
cat("======================================================================\n")
