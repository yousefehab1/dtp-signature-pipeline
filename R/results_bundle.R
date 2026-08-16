# ==============================================================================
# R/results_bundle.R - Results bundle creation, loading, staleness checking,
# and pure live-score calculation for the interactive dashboard.
# ==============================================================================

#' Compute live score for a custom signature pair
#'
#' Evaluates `positive - negative` (or just `positive` if `negative_id` is
#' `"None"`, `NULL`, or empty) on a clinical data frame containing precomputed
#' ssGSEA score columns.
#'
#' @param clin Clinical data frame containing score columns.
#' @param positive_id Required Signature_ID for positive term.
#' @param negative_id Optional Signature_ID for negative term, or `"None"`.
#' @param suffix Score column suffix (defaults to `"_ssGSEA"`).
#' @return Numeric vector of live scores.
#' @export
.compute_live_score <- function(clin, positive_id, negative_id = "None", suffix = "_ssGSEA") {
  if (missing(positive_id) || is.null(positive_id) || is.na(positive_id) || trimws(positive_id) == "") {
    stop("positive_id must be specified.")
  }

  pos_col <- paste0(positive_id, suffix)
  if (!pos_col %in% colnames(clin)) {
    stop("Positive signature '", positive_id, "' (column '", pos_col, "') not found in clinical table.")
  }
  pos_vec <- suppressWarnings(as.numeric(clin[[pos_col]]))

  has_neg <- !is.null(negative_id) && !is.na(negative_id) &&
    trimws(negative_id) != "" && negative_id != "None"

  if (!has_neg) {
    return(pos_vec)
  }

  neg_col <- paste0(negative_id, suffix)
  if (!neg_col %in% colnames(clin)) {
    stop("Negative signature '", negative_id, "' (column '", neg_col, "') not found in clinical table.")
  }
  neg_vec <- suppressWarnings(as.numeric(clin[[neg_col]]))

  pos_vec - neg_vec
}

#' Trim clinical data frame for results bundle
#'
#' Selects only the sample ID, precomputed score columns, endpoint variables,
#' chemotherapy indicator (if present), and harmonised clinical modifiers.
#'
#' @param df Clinical data frame.
#' @param cohort "GSE39582" or "TCGA-COAD".
#' @param cfg Config from [dtp_config()].
#' @return Trimmed tibble.
#' @noRd
.trim_cohort_clinical <- function(df, cohort, cfg) {
  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  harm_df <- harmonize_crc_modifiers(df, cohort)

  # Sample identifier
  id_col <- intersect(c("Sample_ID", "ID", "Patient_ID"), colnames(harm_df))[1]
  if (is.na(id_col)) {
    id_col <- character(0)
  }

  score_cols <- grep(paste0(suffix, "$"), colnames(harm_df), value = TRUE)
  endpoint_cols <- c("Surv_3yr", "Recurrence_3yr", "OS3Y_delay", "OS3Y_event", "RFS3Y_delay", "RFS3Y_event")
  extra_cols <- intersect(c("Chemo_adj"), colnames(harm_df))
  mod_cols <- intersect(cfg$crc_modifiers, colnames(harm_df))

  keep_cols <- unique(c(id_col, score_cols, endpoint_cols, extra_cols, mod_cols))
  avail_cols <- intersect(keep_cols, colnames(harm_df))

  tibble::as_tibble(harm_df[, avail_cols, drop = FALSE])
}

#' Save precomputed results bundle
#'
#' Trims the CRC clinical data frames and bundles them alongside signature
#' annotations, precomputed survival statistics, configuration thresholds, and
#' panel modification timestamp for dashboard consumption.
#'
#' @param cfg Config from [dtp_config()].
#' @param crc_result Result list returned by [run_crc_survival()].
#' @param panel_tbl Canonical signature panel table from [load_signature_panel()].
#' @param composite_defs Composite definitions table from [load_composite_defs()].
#' @param path Destination RDS path; defaults to `"output/results_bundle.rds"`.
#' @return The saved path, invisibly.
#' @export
save_results_bundle <- function(cfg, crc_result, panel_tbl, composite_defs,
                                path = "output/results_bundle.rds") {
  gse_trimmed  <- .trim_cohort_clinical(crc_result$gse_clinical, "GSE39582", cfg)
  tcga_trimmed <- .trim_cohort_clinical(crc_result$tcga_clinical, "TCGA-COAD", cfg)

  cfg_keys <- c("score_suffix", "crc_modifiers", "os_cutpoint", "rfs_cutpoint",
                "min_group_n", "min_km_group_n", "min_events", "min_cox_n", "min_epv")
  cfg_subset <- cfg[intersect(cfg_keys, names(cfg))]

  panel_file <- cfg$signature_panel_path
  panel_mtime <- if (!is.null(panel_file) && file.exists(panel_file)) {
    file.mtime(panel_file)
  } else {
    NA
  }

  bundle <- list(
    gse            = gse_trimmed,
    tcga           = tcga_trimmed,
    panel_tbl      = panel_tbl,
    composite_defs = composite_defs,
    stats_df       = crc_result$stats_df,
    cfg            = cfg_subset,
    built_at       = Sys.time(),
    panel_mtime    = panel_mtime
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bundle, path)
  invisible(path)
}

#' Load precomputed results bundle
#'
#' @param path Path to the RDS file; defaults to `"output/results_bundle.rds"`.
#' @return The loaded bundle list.
#' @export
load_results_bundle <- function(path = "output/results_bundle.rds") {
  if (!file.exists(path)) {
    stop("Results bundle not found at '", path, "'. Run pipeline/build_results_bundle.R first.")
  }
  readRDS(path)
}

#' Check if results bundle is stale relative to signature panel
#'
#' Compares the bundle's recorded `panel_mtime` to the current on-disk modification
#' time of `panel_path`.
#'
#' @param bundle Results bundle list loaded via [load_results_bundle()].
#' @param panel_path Path to the signature panel CSV file; defaults to `"inst/signatures/panel.csv"`.
#' @return Logical `TRUE` if panel file exists and is newer than `bundle$panel_mtime`, `FALSE` otherwise.
#' @export
check_bundle_staleness <- function(bundle, panel_path = "inst/signatures/panel.csv") {
  if (is.null(bundle$panel_mtime) || is.na(bundle$panel_mtime) || !file.exists(panel_path)) {
    return(FALSE)
  }
  curr_mtime <- file.mtime(panel_path)
  isTRUE(curr_mtime > bundle$panel_mtime)
}
