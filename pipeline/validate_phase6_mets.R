#!/usr/bin/env Rscript
# ==============================================================================
# pipeline/validate_phase6_mets.R - End-to-end validation of Phase 6 Metastasis
# Differential Expression and GSEA analysis module against Thesis reference.
# ==============================================================================

suppressPackageStartupMessages({
  library(devtools)
  library(readr)
  library(dplyr)
  library(tibble)
})

cat("======================================================================\n")
cat("Phase 6 Metastasis Differential Expression Validation Pipeline\n")
cat("======================================================================\n\n")

# 1. Load package and configuration
devtools::load_all(".")
cfg <- dtp_config()
panel_tbl <- load_signature_panel(cfg)

# 2. Run Metastasis DE & GSEA Analysis
cat("\n--- Step 2: Running run_mets_de() ---\n")
t0 <- Sys.time()
mets_result <- run_mets_de(cfg, panel_tbl)
t1 <- Sys.time()
duration <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("run_mets_de completed in %.2f seconds.\n", duration))

# 3. Print dimensions of all returned structures
cat("\n--- Step 3: Table Dimensions & Structure Summary ---\n")
for (nm in names(mets_result)) {
  val <- mets_result[[nm]]
  if (is.data.frame(val)) {
    d <- dim(val)
    cat(sprintf("  %-18s: %5d rows x %2d cols\n", nm, d[1], d[2]))
  } else if (is.numeric(val)) {
    cat(sprintf("  %-18s: numeric vector (length %d)\n", nm, length(val)))
  } else if (is.list(val)) {
    cat(sprintf("  %-18s: list of %d elements\n", nm, length(val)))
  } else {
    cat(sprintf("  %-18s: class %s\n", nm, paste(class(val), collapse = "/")))
  }
}

ref_dir <- "/Users/elabd/Master's Work/Thesis/CRC_DTP_20260815_1705/mets_de/results"

compare_tables <- function(new_df, ref_path, key_cols, num_cols, label, map_csc_to_cbc = FALSE) {
  cat(sprintf("\n=== Validation: %s ===\n", label))
  if (!file.exists(ref_path)) {
    cat(sprintf("  [!] Reference file not found: %s\n", ref_path))
    return(invisible(NULL))
  }
  ref_df <- readr::read_csv(ref_path, show_col_types = FALSE)
  if (map_csc_to_cbc) {
    if ("Signature" %in% colnames(ref_df)) {
      ref_df$Signature[ref_df$Signature == "CSC"] <- "CBC"
    }
    if ("ID" %in% colnames(ref_df)) {
      ref_df$ID[ref_df$ID == "CSC"] <- "CBC"
    }
    if ("CSC_ssGSEA" %in% colnames(ref_df)) {
      colnames(ref_df)[colnames(ref_df) == "CSC_ssGSEA"] <- "CBC_ssGSEA"
    }
  }
  cat(sprintf("  Reference rows: %d, New rows: %d\n", nrow(ref_df), nrow(new_df)))

  merged <- dplyr::inner_join(new_df, ref_df, by = key_cols, suffix = c(".new", ".ref"))
  cat(sprintf("  Matched rows on (%s): %d / %d reference rows\n",
              paste(key_cols, collapse = ", "), nrow(merged), nrow(ref_df)))

  for (col in num_cols) {
    new_col_name <- paste0(col, ".new")
    ref_col_name <- paste0(col, ".ref")
    if (!new_col_name %in% colnames(merged) || !ref_col_name %in% colnames(merged)) {
      cat(sprintf("  %-18s: [!] Column %s not found in merged data\n", col, col))
      next
    }
    c_new <- merged[[new_col_name]]
    c_ref <- merged[[ref_col_name]]
    ok <- !is.na(c_new) & !is.na(c_ref)
    if (sum(ok) > 0) {
      denom <- pmax(abs(c_ref[ok]), 1e-6)
      rel_diff <- abs(c_new[ok] - c_ref[ok]) / denom
      max_rd <- max(rel_diff)
      mean_rd <- mean(rel_diff)
      cat(sprintf("  %-18s: max rel diff = %.6f, mean rel diff = %.6f (n=%d)\n",
                  col, max_rd, mean_rd, sum(ok)))
    } else {
      cat(sprintf("  %-18s: [!] No non-NA values to compare\n", col))
    }
  }
}

# 4. Compare liver validation scores
compare_tables(
  new_df   = mets_result$liver_df,
  ref_path = file.path(ref_dir, "PurityAdjusted/liver_validation_scores_3way.csv"),
  key_cols = "GSM",
  num_cols = "Liver",
  label    = "Liver Validation Scores (3-Way)"
)

# 5. Compare DE unadjusted vs purity-adjusted
compare_tables(
  new_df   = mets_result$comparison_de,
  ref_path = file.path(ref_dir, "PurityAdjusted/DE_unadjusted_vs_purityAdjusted.csv"),
  key_cols = "gene",
  num_cols = c("logFC_unadj", "t_unadj", "pval_unadj", "padj_unadj",
               "logFC_adj", "t_adj", "pval_adj", "padj_adj"),
  label    = "DE Unadjusted vs Purity-Adjusted"
)

# 6. Compare single-sample paired signature scores
score_cols <- intersect(
  c("Up_ssGSEA", "Down_ssGSEA", "Fetal_ssGSEA", "revSC_ssGSEA",
    "RSC_ssGSEA", "IBD_ssGSEA", "CBC_ssGSEA", "MYC_ssGSEA", "Composite_ssGSEA"),
  colnames(mets_result$scores_df)
)
compare_tables(
  new_df          = mets_result$scores_df,
  ref_path        = file.path(ref_dir, "SingleSampleDiagnostics/signature_scores_paired.csv"),
  key_cols        = "GSM",
  num_cols        = score_cols,
  label           = "Single-Sample Paired ssGSEA Scores (CSC mapped to canonical CBC)",
  map_csc_to_cbc  = TRUE
)

# 7. Compare ROC AUC summary
compare_tables(
  new_df          = mets_result$roc_summary,
  ref_path        = file.path(ref_dir, "SingleSampleDiagnostics/ROC_AUC_summary.csv"),
  key_cols        = "Signature",
  num_cols        = "AUC",
  label           = "Single-Sample ROC AUC Summary (CSC mapped to canonical CBC)",
  map_csc_to_cbc  = TRUE
)

# 8. Compare ssGSEA paired p-values
compare_tables(
  new_df          = mets_result$ssgsea_pvalues,
  ref_path        = file.path(ref_dir, "SingleSampleDiagnostics/ssGSEA_paired_pvalues_BHadjusted.csv"),
  key_cols        = "Signature",
  num_cols        = c("p_raw", "p_adj_BH"),
  label           = "ssGSEA Paired Wilcoxon P-values (BH adjusted, CSC mapped to CBC)",
  map_csc_to_cbc  = TRUE
)

# 9. Compare GSEA Primary vs Metastasis (unadjusted)
compare_tables(
  new_df          = mets_result$gsea_pm,
  ref_path        = file.path(ref_dir, "GSEA_PrimaryVsMetastasis_Summary.csv"),
  key_cols        = "ID",
  num_cols        = c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue"),
  label           = "GSEA Primary vs Metastasis (Unadjusted, CSC mapped to CBC)",
  map_csc_to_cbc  = TRUE
)

# 10. Compare GSEA Primary vs Metastasis (purity-adjusted)
compare_tables(
  new_df          = mets_result$gsea_pm_adj,
  ref_path        = file.path(ref_dir, "GSEA_PrimaryVsMetastasis_purityAdjusted_Summary.csv"),
  key_cols        = "ID",
  num_cols        = c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue"),
  label           = "GSEA Primary vs Metastasis (Purity-Adjusted, CSC mapped to CBC)",
  map_csc_to_cbc  = TRUE
)

# 11. Compare GSEA Primary vs Normal (sign-aligned comparison against legacy NormalVsPrimary file)
cat("\n=== Validation: GSEA Primary vs Normal (Sign-Aligned) ===\n")
cat("  Note: Old disk CSV 'GSEA_NormalVsPrimary_Summary.csv' used Normal - Primary ranking.\n")
cat("  Per old code lines 335-341, Primary vs Normal negates the coefficient so positive = Up in Primary.\n")
ref_pn <- readr::read_csv(file.path(ref_dir, "GSEA_NormalVsPrimary_Summary.csv"), show_col_types = FALSE)
ref_pn$ID[ref_pn$ID == "CSC"] <- "CBC"
# Sign-flip the reference ES and NES to align polarity (positive = Up in Primary)
ref_pn_aligned <- ref_pn %>%
  dplyr::mutate(
    enrichmentScore = -enrichmentScore,
    NES             = -NES
  )

merged_pn <- dplyr::inner_join(mets_result$gsea_pn, ref_pn_aligned, by = "ID", suffix = c(".new", ".ref"))
cat(sprintf("  Matched rows on (ID): %d / %d reference rows\n", nrow(merged_pn), nrow(ref_pn)))
pn_cols <- c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue")
for (col in pn_cols) {
  c_new <- merged_pn[[paste0(col, ".new")]]
  c_ref <- merged_pn[[paste0(col, ".ref")]]
  ok <- !is.na(c_new) & !is.na(c_ref)
  denom <- pmax(abs(c_ref[ok]), 1e-6)
  rel_diff <- abs(c_new[ok] - c_ref[ok]) / denom
  cat(sprintf("  %-18s: max rel diff = %.6f, mean rel diff = %.6f (n=%d)\n",
              col, max(rel_diff), mean(rel_diff), sum(ok)))
}

cat("\n======================================================================\n")
cat("Phase 6 Metastasis Validation Complete!\n")
cat("======================================================================\n")
