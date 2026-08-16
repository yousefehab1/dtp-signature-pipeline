# ==============================================================================
# dashboard/global.R - Global setup for the dtpsig Shiny dashboard.
#
# Supports launching via `shiny::runApp("dashboard")` from the repository root
# or launching from within the `dashboard/` directory.
# ==============================================================================

# 1. Package loading
if (requireNamespace("devtools", quietly = TRUE)) {
  if (file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else if (file.exists("../DESCRIPTION")) {
    devtools::load_all("..", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(dtpsig)
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(survival)
  library(dplyr)
})

# 2. Source modules if not already sourced
mod_files <- c(
  "dashboard/R/mod_signature_explorer.R",
  "R/mod_signature_explorer.R"
)
for (mf in mod_files) {
  if (file.exists(mf)) {
    source(mf, local = FALSE)
    break
  }
}

# 3. Locate and load results bundle
bundle_paths <- c(
  "output/results_bundle.rds",
  "../output/results_bundle.rds"
)
bundle_path <- NULL
for (bp in bundle_paths) {
  if (file.exists(bp)) {
    bundle_path <- bp
    break
  }
}

panel_paths <- c(
  "inst/signatures/panel.csv",
  "../inst/signatures/panel.csv"
)
panel_path <- NULL
for (pp in panel_paths) {
  if (file.exists(pp)) {
    panel_path <- pp
    break
  }
}

bundle <- tryCatch({
  if (!is.null(bundle_path)) {
    load_results_bundle(bundle_path)
  } else {
    warning("No results bundle found. Please run pipeline/build_results_bundle.R first.")
    NULL
  }
}, error = function(e) {
  warning("Failed to load results bundle: ", e$message)
  NULL
})

is_bundle_stale <- if (!is.null(bundle) && !is.null(panel_path)) {
  check_bundle_staleness(bundle, panel_path)
} else {
  FALSE
}
