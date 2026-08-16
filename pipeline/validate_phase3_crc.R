#!/usr/bin/env Rscript
# Validation script for Phase 3: CRC Survival Analysis Module

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

cfg <- dtp_config()
init_run(cfg)

out_dir <- "output/phase3_validate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_tbl <- load_signature_panel(cfg)
composite_defs <- load_composite_defs(cfg)

message("Running CRC survival analysis pipeline...")
result <- run_crc_survival(cfg, panel_tbl, out_root = out_dir)

# Write output CSVs
write.csv(csv_safe(result$stats_df), file.path(out_dir, "stats_df.csv"), row.names = FALSE)
write.csv(csv_safe(result$adj_df), file.path(out_dir, "adj_df.csv"), row.names = FALSE)
write.csv(csv_safe(result$int_df), file.path(out_dir, "int_df.csv"), row.names = FALSE)
if (!is.null(result$level_df)) {
  write.csv(csv_safe(result$level_df), file.path(out_dir, "level_df.csv"), row.names = FALSE)
}
write.csv(csv_safe(result$subtype_stats), file.path(out_dir, "subtype_stats.csv"), row.names = FALSE)

cat("\n=== Validation Results for Phase 3 (CRC Survival) ===\n\n")
cat("dim(stats_df):", dim(result$stats_df), "\n")
cat("dim(adj_df):", dim(result$adj_df), "\n")
cat("dim(int_df):", dim(result$int_df), "\n")
cat("dim(level_df):", dim(result$level_df), "\n")
cat("dim(subtype_stats):", dim(result$subtype_stats), "\n\n")
cat("head(stats_df, 5):\n")
print(head(as.data.frame(result$stats_df), 5))
