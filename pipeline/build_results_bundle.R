#!/usr/bin/env Rscript
# ==============================================================================
# pipeline/build_results_bundle.R - Builds and saves the precomputed results
# bundle for the interactive dashboard.
# ==============================================================================

suppressPackageStartupMessages({
  library(devtools)
})

cat("======================================================================\n")
cat("Building DTP Results Bundle (CRC Data)\n")
cat("======================================================================\n\n")

# 1. Load package and configuration
devtools::load_all(".")
cfg <- dtp_config()
panel_tbl <- load_signature_panel(cfg)
composite_defs <- load_composite_defs(cfg)

# 2. Run CRC survival analysis (uses warm cache)
cat("\n--- Running run_crc_survival() ---\n")
t0 <- Sys.time()
crc_result <- run_crc_survival(cfg, panel_tbl, out_root = "output")
t1 <- Sys.time()
cat(sprintf("run_crc_survival completed in %.2f seconds.\n", as.numeric(difftime(t1, t0, units = "secs"))))

# 3. Save results bundle
cat("\n--- Saving Results Bundle ---\n")
out_path <- "output/results_bundle.rds"
save_results_bundle(cfg, crc_result, panel_tbl, composite_defs, path = out_path)

# 4. Verify and print diagnostics
file_bytes <- file.size(out_path)
file_size_kb <- file_bytes / 1024
file_size_mb <- file_size_kb / 1024

cat(sprintf("Bundle saved to: %s\n", out_path))
cat(sprintf("File size: %d bytes (%.2f KB / %.2f MB)\n\n", file_bytes, file_size_kb, file_size_mb))

bundle <- load_results_bundle(out_path)
cat("--- Bundle Structure (str) ---\n")
utils::str(bundle, max.level = 2)

cat("\n======================================================================\n")
cat("DTP Results Bundle Build Complete!\n")
cat("======================================================================\n")
