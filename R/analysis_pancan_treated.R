# Pan-cancer survival analysis in treated patients only (pharmaceutical or radiation therapy).
# Reuses the already-scored in-memory master clinical data frame from run_pancan_survival().

# Per-cohort TCGA treatment flags. Cached per project so the GDC clinical pull
# happens once and is incremental. We cache the RAW pharmaceutical / radiation
# yes-no flags rather than a derived label, so the "treated" definition can be
# changed downstream without re-pulling from GDC.
#
# IMPORTANT BUSINESS RULE: A network FAILURE is NOT cached (it would otherwise pin
# the cohort at zero treated forever) — only a real response is, including a cohort
# that genuinely lacks one or both treatment fields. We use manual caching here
# instead of cache_rds() because cache_rds() cannot distinguish "don't cache this"
# from a normal return value.
.pancan_treatment_flags <- function(cfg, projects) {
  has_value <- function(col, val) vapply(col, function(x) val %in% x, logical(1))

  # A cohort can be missing either field independently; absent => FALSE both ways,
  # which lands the patient in "Unknown/Not Reported" downstream.
  flag <- function(cl, field, val) {
    if (field %in% names(cl)) has_value(cl[[field]], val) else rep(FALSE, nrow(cl))
  }

  empty <- data.frame(Patient_ID = character(), Project_ID = character(),
                      pharma_yes = logical(), pharma_no = logical(),
                      rad_yes = logical(), rad_no = logical(), stringsAsFactors = FALSE)
  PHARMA <- "treatments_pharmaceutical_treatment_or_therapy"
  RAD    <- "treatments_radiation_treatment_or_therapy"

  cfile <- function(proj) {
    fname <- paste0(gsub("[^A-Za-z0-9_]+", "_", paste0("treat_flags_", proj)), ".rds")
    file.path(cfg$cache_dir, fname)
  }

  dplyr::bind_rows(lapply(projects, function(proj) {
    f <- cfile(proj)
    if (file.exists(f)) {
      message("[cache] treat_flags_", proj, " (load)")
      return(readRDS(f))
    }

    # Warm-through from legacy cache if present
    if (!is.null(cfg$legacy_cache_dir)) {
      fname    <- paste0(gsub("[^A-Za-z0-9_]+", "_", paste0("treat_flags_", proj)), ".rds")
      f_legacy <- file.path(cfg$legacy_cache_dir, fname)
      if (file.exists(f_legacy)) {
        message("[cache] treat_flags_", proj, " (warm from legacy cache)")
        dir.create(cfg$cache_dir, recursive = TRUE, showWarnings = FALSE)
        file.copy(f_legacy, f)
        return(readRDS(f))
      }
    }

    cl <- tryCatch(suppressWarnings(TCGAbiolinks::GDCquery_clinic(project = proj, type = "clinical")),
                   error = function(e) e)
    if (inherits(cl, "error")) {
      message("  [!] clinical pull failed for ", proj, ": ", conditionMessage(cl),
              " — NOT cached, will retry next run")
      return(empty)
    }

    if (!any(c(PHARMA, RAD) %in% names(cl))) {
      message("  [!] ", proj, ": no pharmaceutical or radiation field — all Unknown")
    }

    df <- data.frame(
      Patient_ID = trimws(toupper(gsub("\\.", "-", cl$submitter_id))),
      Project_ID = proj,
      pharma_yes = flag(cl, PHARMA, "yes"), pharma_no = flag(cl, PHARMA, "no"),
      rad_yes    = flag(cl, RAD,    "yes"), rad_no    = flag(cl, RAD,    "no"),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::distinct(Patient_ID, .keep_all = TRUE)

    dir.create(cfg$cache_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(df, f) # Cache successful pulls only (incl. legit no-field cohorts)
    df
  }))
}

# Run per-cohort and batch-corrected aggregate statistics on treated patients with FDR.
.pancan_treated_stats_block <- function(mt_all, counts, score_cols, cfg) {
  rows <- list()

  # ---- Phase 4a: per-cohort stats (uncorrected) on EVERY treated cohort ----
  for (proj in sort(unique(mt_all$Project_ID))) {
    cd <- mt_all %>% dplyr::filter(Project_ID == proj)
    for (sc in score_cols) {
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "OS",  "PanCancer", proj, sc, "Surv_3yr", cfg)
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "RFS", "PanCancer", proj, sc, "Recurrence_3yr", cfg)
    }
  }

  # ---- Phase 4b: batch-corrected aggregate over well-populated treated cohorts ----
  agg_min <- cfg$min_cox_n
  keep    <- counts$Project_ID[counts$N_treated >= agg_min]
  mt_agg  <- mt_all %>% dplyr::filter(Project_ID %in% keep)
  message(sprintf("  aggregate: %d of %d treated cohorts have >= %d treated patients",
                  length(keep), dplyr::n_distinct(mt_all$Project_ID), agg_min))

  if (length(keep) >= 2) {
    sm        <- t(as.matrix(mt_agg[, score_cols, drop = FALSE]))
    corrected <- limma::removeBatchEffect(sm, batch = mt_agg$Project_ID)
    corr_df   <- as.data.frame(t(corrected))
    corr_cols <- paste0(score_cols, "_Corrected")
    colnames(corr_df) <- corr_cols

    # Update or append corrected score columns for the aggregate cohort
    mt_agg <- mt_agg %>%
      dplyr::select(-dplyr::any_of(corr_cols)) %>%
      dplyr::bind_cols(corr_df)

    for (sc in corr_cols) {
      rows[[length(rows) + 1]] <- run_all_stats(mt_agg, sc, "OS",  "PanCancer", "Pan-Cancer", sc, "Surv_3yr", cfg)
      rows[[length(rows) + 1]] <- run_all_stats(mt_agg, sc, "RFS", "PanCancer", "Pan-Cancer", sc, "Recurrence_3yr", cfg)
    }
  } else {
    message("  [!] < 2 cohorts meet the treated floor — skipping batch-corrected aggregate.")
  }

  # ---- FDR (Family x Test grouping) ----
  stats_df <- do.call(rbind, rows) %>%
    dplyr::mutate(Family = ifelse(Project == "Pan-Cancer", "PanCancer", "PerCohort"))
  apply_fdr(stats_df, cfg, family = "pancan")
}

#' Run pan-cancer survival analysis in treated patients
#'
#' Filters the scored pan-cancer master data to patients who received pharmaceutical
#' or radiation therapy, evaluating per-cohort and batch-corrected aggregate survival statistics.
#'
#' @param cfg Config from [dtp_config()].
#' @param pancan_result Result list returned by [run_pancan_survival()].
#' @param out_root Output directory path (reserved for figure generation in future phases).
#' @return A named list containing:
#'   \item{stats_df}{FDR-corrected survival statistics for treated cohorts and aggregate.}
#'   \item{master_treated}{Clinical data frame restricted to treated patients.}
#'   \item{treatment_counts}{Per-cohort attrition summary of treated, untreated, and unknown patients.}
#' @export
run_pancan_treated <- function(cfg, pancan_result, out_root = "output") {
  force(out_root)
  stopifnot(!is.null(pancan_result$master), !is.null(pancan_result$stats_df))

  message("\n======== PAN-CANCER SURVIVAL — TREATED (pharmaceutical or radiation) ========")

  master     <- pancan_result$master
  score_cols <- grep(paste0(cfg$score_suffix, "$"), colnames(master), value = TRUE)

  # ---- Treatment status (pharmaceutical OR radiation), cached per cohort ----
  treat <- .pancan_treatment_flags(cfg, sort(unique(master$Project_ID)))
  master <- master %>%
    dplyr::left_join(
      dplyr::select(treat, Patient_ID, pharma_yes, pharma_no, rad_yes, rad_no),
      by = "Patient_ID"
    ) %>%
    dplyr::mutate(Treatment_Status = dplyr::case_when(
      pharma_yes | rad_yes ~ "Treated",
      pharma_no  & rad_no  ~ "Not Treated",
      TRUE                 ~ "Unknown/Not Reported"
    )) %>%
    dplyr::select(-pharma_yes, -pharma_no, -rad_yes, -rad_no)

  # ---- Attrition report ----
  counts <- master %>%
    dplyr::group_by(Project_ID) %>%
    dplyr::summarise(
      N_total     = dplyr::n(),
      N_treated   = sum(Treatment_Status == "Treated"),
      N_untreated = sum(Treatment_Status == "Not Treated"),
      N_unknown   = sum(Treatment_Status == "Unknown/Not Reported"),
      .groups     = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(N_treated))

  message(sprintf("  Treated: %d / %d (%.0f%% pharma or radiation)  |  Not treated: %d  |  Unknown: %d",
                  sum(master$Treatment_Status == "Treated"), nrow(master),
                  100 * mean(master$Treatment_Status == "Treated"),
                  sum(master$Treatment_Status == "Not Treated"),
                  sum(master$Treatment_Status == "Unknown/Not Reported")))

  mt_all <- master %>% dplyr::filter(Treatment_Status == "Treated")
  if (!nrow(mt_all)) {
    stop("run_pancan_treated(): no treated patients found.")
  }

  # ---- Statistics ----
  stats_df <- .pancan_treated_stats_block(mt_all, counts, score_cols, cfg)

  message("Pan-cancer treated section complete.")
  invisible(list(
    stats_df         = stats_df,
    master_treated   = mt_all,
    treatment_counts = counts
  ))
}
