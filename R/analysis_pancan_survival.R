# Pan-cancer survival analysis module across TCGA cohorts.
# Decomposed into small, testable internal helpers (.pancan_*) and orchestrated
# by run_pancan_survival().

# Returns TRUE for errors that may resolve on retry (network, server-side).
# vctrs_error_subscript_oob and similar parsing errors are NOT transient.
.pancan_is_transient_gdc_error <- function(e) {
  cls <- class(e)
  if (any(grepl("vctrs_error|subscript|parse|type|coerce|match", cls, ignore.case = TRUE))) {
    return(FALSE)
  }
  msg <- tolower(paste(conditionMessage(e), collapse = " "))
  any(grepl("timeout|connection|network|curl|ssl|503|502|500|temporarily", msg))
}

# Build a minimal SummarizedExperiment directly from the downloaded STAR-Counts
# TSV files, bypassing GDCprepare entirely. Used as a fallback when
# GDCprepare's clinical-annotation step crashes with vctrs_error_subscript_oob.
# The resulting SE has exactly what prep_tcga_tpm() needs:
#   - assay "tpm_unstrand"
#   - rowData: gene_id, gene_name
#   - colData: sample_type (inferred from the TCGA barcode sample-code field)
.pancan_read_star_tpm_as_se <- function(query, cfg = NULL) {
  manifest <- as.data.frame(TCGAbiolinks::getResults(query))
  proj     <- manifest$project[1]
  gdc_dir <- file.path("GDCdata", proj,
                       "Transcriptome_Profiling",
                       "Gene_Expression_Quantification")

  # Locate one TSV per file-UUID directory
  tsv_paths <- vapply(manifest$id, function(fid) {
    d <- file.path(gdc_dir, fid)
    f <- list.files(d, pattern = "\\.tsv$", full.names = TRUE)
    if (length(f) == 0) stop("No TSV found in: ", d)
    f[1]
  }, character(1))

  # Gene annotation from the first file (identical across all samples)
  gene_ann <- readr::read_tsv(tsv_paths[1], comment = "#", show_col_types = FALSE) %>%
    dplyr::filter(!grepl("^N_", gene_id)) %>%
    dplyr::select(gene_id, gene_name)

  # TPM matrix: genes x samples (vapply -> nrow(gene_ann) x n_samples)
  tpm_mat <- vapply(tsv_paths, function(f) {
    df <- readr::read_tsv(f, comment = "#", show_col_types = FALSE,
                          col_select = c("gene_id", "tpm_unstranded")) %>%
      dplyr::filter(!grepl("^N_", gene_id))
    df$tpm_unstranded[match(gene_ann$gene_id, df$gene_id)]
  }, numeric(nrow(gene_ann)))

  # Column names: prefer sample barcode, fall back to case/patient ID.
  barcode_col <- intersect(c("sample.submitter_id", "cases.submitter_id", "cases"),
                           colnames(manifest))[1]
  if (is.na(barcode_col)) {
    stop("Cannot find a barcode column in GDC manifest. Available: ",
         paste(colnames(manifest), collapse = ", "))
  }
  barcodes <- manifest[[barcode_col]]
  if (anyNA(barcodes) || length(barcodes) != ncol(tpm_mat)) {
    stop(sprintf("Barcode vector length (%d) != matrix columns (%d); NAs: %d",
                 length(barcodes), ncol(tpm_mat), sum(is.na(barcodes))))
  }
  colnames(tpm_mat) <- barcodes

  # sample_type: use manifest column if present, else infer from barcode
  # TCGA barcode field 4 (1-based, 14-15 chars): "01"=Primary, "06"=Metastatic, "11"=Normal
  sample_type <- if ("sample_type" %in% colnames(manifest)) {
    manifest$sample_type
  } else {
    codes <- substr(barcodes, 14, 15)
    dplyr::case_when(
      codes == "01" ~ "Primary Tumor",
      codes == "06" ~ "Metastatic",
      codes == "11" ~ "Solid Tissue Normal",
      TRUE          ~ paste0("Sample_", codes)
    )
  }

  SummarizedExperiment::SummarizedExperiment(
    assays  = list(tpm_unstrand = tpm_mat),
    rowData = S4Vectors::DataFrame(gene_id   = gene_ann$gene_id,
                                   gene_name = gene_ann$gene_name,
                                   row.names = gene_ann$gene_id),
    colData = S4Vectors::DataFrame(sample_type = sample_type,
                                   row.names   = barcodes)
  )
}

# GDCprepare with direct-read fallback.
# If GDCprepare fails with vctrs_error_subscript_oob (TCGAbiolinks clinical
# annotation bug), we fall back to reading the already-downloaded TSV files
# directly. All other errors are re-raised.
.pancan_gdc_prepare_safe <- function(query, cfg = NULL) {
  se <- tryCatch(TCGAbiolinks::GDCprepare(query), error = function(e) e)
  if (!inherits(se, "error")) return(se)

  is_clinical_err <- inherits(se, "vctrs_error_subscript_oob") ||
    any(grepl("vctrs_error_subscript", class(se), ignore.case = TRUE))
  if (!is_clinical_err) stop(se)

  message("   -> clinical annotation failed; reading TSV files directly")
  .pancan_read_star_tpm_as_se(query, cfg = cfg)
}

# Download + prepare one TCGA project with up to 3 retries on transient
# network/server errors. Cached via cache_rds so it only runs on cache miss.
.pancan_gdc_with_retry <- function(cfg, query, proj, retries = 3, wait_sec = 15) {
  cache_rds(cfg, paste0("se_", proj), function() {
    last_err <- NULL
    for (attempt in seq_len(retries)) {
      result <- tryCatch({
        TCGAbiolinks::GDCdownload(query, method = "api", files.per.chunk = 6)
        .pancan_gdc_prepare_safe(query, cfg = cfg)
      }, error = function(e) e)
      if (!inherits(result, "error")) return(result)
      last_err <- result
      msg <- if (nzchar(trimws(last_err$message))) last_err$message else class(last_err)[1]
      # Don't retry non-transient errors (data/parsing failures won't resolve)
      if (!.pancan_is_transient_gdc_error(last_err)) {
        message(sprintf("   [!] non-transient error, skipping retries: %s", msg))
        break
      }
      if (attempt < retries) {
        message(sprintf("   [retry %d/%d] %s - waiting %ds", attempt, retries, msg, wait_sec))
        Sys.sleep(wait_sec)
      }
    }
    stop(last_err)
  })
}

# Build a TPM matrix for proj by chunking valid_patients into groups of
# chunk_size, running GDCprepare separately on each, and cbind-ing results.
# Files are already downloaded; only the in-memory SE construction is chunked.
.pancan_prep_tpm_chunked <- function(cfg, proj, valid_patients, chunk_size = 300L) {
  chunks <- split(valid_patients, ceiling(seq_along(valid_patients) / chunk_size))
  message(sprintf("   -> %d samples split into %d chunks of ~%d (memory optimisation)",
                  length(valid_patients), length(chunks), chunk_size))
  tpm_chunks <- lapply(seq_along(chunks), function(k) {
    bc <- chunks[[k]]
    message(sprintf("   -> chunk %d/%d  (%d samples)", k, length(chunks), length(bc)))
    q <- TCGAbiolinks::GDCquery(project       = proj,
                                data.category = "Transcriptome Profiling",
                                data.type     = "Gene Expression Quantification",
                                workflow.type = "STAR - Counts",
                                barcode       = bc)
    se <- .pancan_gdc_prepare_safe(q, cfg = cfg)
    m  <- prep_tcga_tpm(cfg, se)
    rm(se); gc()
    m
  })
  genes <- rownames(tpm_chunks[[1]])
  do.call(cbind, lapply(tpm_chunks, function(m) m[genes, , drop = FALSE]))
}

# CDR inventory and valid/excluded project partitioning.
.pancan_cohort_inventory <- function(cfg) {
  cdr_clinical <- load_tcga_cdr(cfg)
  cohort_summary <- cdr_clinical %>%
    dplyr::group_by(Project_ID) %>%
    dplyr::summarise(CDR_Patients = dplyr::n(),
                     OS_Available  = sum(!is.na(OS_months)  & !is.na(OS_event)),
                     DFS_Available = sum(!is.na(DFS_months) & !is.na(DFS_event)),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(CDR_Patients))

  valid_projects <- cohort_summary %>%
    dplyr::filter(!Project_ID %in% cfg$exclude_projects, CDR_Patients >= cfg$min_cohort_n) %>%
    dplyr::pull(Project_ID)

  excluded_df <- cohort_summary %>%
    dplyr::mutate(
      Excluded = Project_ID %in% cfg$exclude_projects | CDR_Patients < cfg$min_cohort_n,
      Reason   = dplyr::case_when(
        Project_ID %in% cfg$exclude_projects ~ paste0("Manually excluded (non-solid / DTP N/A): ",
                                                      paste(cfg$exclude_projects, collapse = ", ")),
        CDR_Patients < cfg$min_cohort_n      ~ sprintf("Too few CDR patients (%d < threshold %d)",
                                                        CDR_Patients, cfg$min_cohort_n),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(Excluded) %>%
    dplyr::select(Project_ID, CDR_Patients, Reason)

  if (nrow(excluded_df) > 0) {
    message(sprintf("\n--- Pre-queue exclusions (%d cancer type(s)) ---", nrow(excluded_df)))
    for (i in seq_len(nrow(excluded_df))) {
      message(sprintf("  [EXCLUDED] %-15s  n=%4d  |  %s",
                      excluded_df$Project_ID[i], excluded_df$CDR_Patients[i], excluded_df$Reason[i]))
    }
    message("---")
  }
  message(sprintf("%d projects queued.", length(valid_projects)))

  list(
    cdr_clinical   = cdr_clinical,
    cohort_summary = cohort_summary,
    valid_projects = valid_projects,
    excluded_df    = excluded_df
  )
}

# Single cohort scoring, coverage gating, ssGSEA, composite addition, and CDR join.
.pancan_score_cohort <- function(cfg, proj, cdr_proj, gs_default, gs_rest, composite_defs, panel_list) {
  valid_patients <- unique(cdr_proj$Patient_ID)
  if (length(valid_patients) == 0) {
    message("   -> no CDR patients.")
    return(list(status = "Skipped", n = 0L, reason = "No CDR patients", clinical = NULL))
  }

  # Two-level cache: build TPM matrix once, reload the small RDS on reruns.
  # Large cohorts (> 300 samples) are chunked to keep peak memory bounded.
  mat <- cache_rds(cfg, paste0("tpm_", proj, "_", cfg$id_type), function() {
    if (length(valid_patients) > 300L) {
      .pancan_prep_tpm_chunked(cfg, proj, valid_patients, chunk_size = 300L)
    } else {
      query <- TCGAbiolinks::GDCquery(project       = proj,
                                      data.category = "Transcriptome Profiling",
                                      data.type     = "Gene Expression Quantification",
                                      workflow.type = "STAR - Counts",
                                      barcode       = valid_patients)
      se <- .pancan_gdc_with_retry(cfg, query, proj)
      m  <- prep_tcga_tpm(cfg, se)
      rm(se); gc()
      m
    }
  })

  # Coverage gate: check coverage for each signature in the panel
  covs <- vapply(panel_list, function(genes) gene_set_coverage(mat, genes), numeric(1))
  cov_strs <- paste0(names(covs), "=", sprintf("%.1f%%", covs))
  message(sprintf("   -> coverage: %s", paste(cov_strs, collapse = ", ")))

  # Strict gating: cohort is skipped if ANY panel signature falls below threshold (mirrors legacy Up/Down check)
  low_covs <- covs[covs < cfg$min_gene_cov]
  if (length(low_covs) > 0) {
    message("   -> coverage below threshold.")
    return(list(
      status   = "Skipped",
      n        = NA_integer_,
      reason   = sprintf("Gene coverage too low (%s; threshold=%.0f%%)",
                         paste0(names(low_covs), "=", sprintf("%.1f%%", low_covs), collapse = ", "),
                         cfg$min_gene_cov),
      clinical = NULL
    ))
  }

  # GSVA's ssGSEA normalizes each sample's enrichment scores against the range
  # of ALL gene sets scored in that same call: (ES - min(ES)) / (max(ES) - min(ES)).
  # Scoring the default composite pair (e.g. Up/Down) in an isolated call preserves
  # exact numerical parity with the legacy thesis pipeline, while the remaining panel
  # signatures are scored in a separate isolated call before combining.
  scores_def  <- if (!is.null(gs_default)) run_ssgsea(mat, gs_default, suffix = cfg$score_suffix) else NULL
  scores_rest <- if (!is.null(gs_rest))    run_ssgsea(mat, gs_rest,    suffix = cfg$score_suffix) else NULL

  scores <- if (!is.null(scores_def) && !is.null(scores_rest)) {
    stopifnot(identical(sort(rownames(scores_def)), sort(rownames(scores_rest))))
    cbind(scores_def, scores_rest[rownames(scores_def), , drop = FALSE])
  } else if (!is.null(scores_def)) {
    scores_def
  } else {
    scores_rest
  }

  scores <- add_all_composites(scores, composite_defs, suffix = cfg$score_suffix)
  scores <- scores %>%
    tibble::rownames_to_column("Full_Barcode") %>%
    dplyr::arrange(Full_Barcode) %>%
    dplyr::mutate(Patient_ID = tcga_patient_id(Full_Barcode)) %>%
    dplyr::distinct(Patient_ID, .keep_all = TRUE) %>%
    dplyr::select(-Full_Barcode)

  clinical <- scores %>%
    dplyr::inner_join(
      cdr_proj %>% dplyr::select(Patient_ID, OS_event, OS_months, DFS_event, DFS_months),
      by = "Patient_ID"
    ) %>%
    dplyr::mutate(
      Project_ID = proj,
      OS_delay = OS_months,
      Recurrence_event = DFS_event,
      Recurrence_delay = DFS_months
    ) %>%
    derive_endpoints(
      os_event = "OS_event", os_time = "OS_delay",
      rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay",
      os_cut = cfg$os_cutpoint, rfs_cut = cfg$rfs_cutpoint
    )

  if (nrow(clinical) < cfg$min_cohort_n) {
    message(sprintf("   -> %d matched < MIN_COHORT_N=%d.", nrow(clinical), cfg$min_cohort_n))
    return(list(
      status   = "Skipped",
      n        = nrow(clinical),
      reason   = sprintf("Only %d patients after CDR join (threshold=%d)", nrow(clinical), cfg$min_cohort_n),
      clinical = NULL
    ))
  }

  message(sprintf("   -> processed %d patients.", nrow(clinical)))
  list(status = "Success", n = nrow(clinical), reason = NA_character_, clinical = clinical)
}

# Run scoring across all valid projects in-memory and compile cohort execution summary.
.pancan_run_all_cohorts <- function(cfg, valid_projects, cdr_clinical, gs_default, gs_rest, composite_defs, panel_list) {
  cohort_outcomes <- vector("list", length(valid_projects))
  names(cohort_outcomes) <- valid_projects
  cohort_clinicals <- list()

  for (i in seq_along(valid_projects)) {
    proj <- valid_projects[i]
    message(sprintf("\n[%d/%d] %s", i, length(valid_projects), proj))
    cdr_proj <- cdr_clinical %>% dplyr::filter(Project_ID == proj)

    res <- tryCatch({
      .pancan_score_cohort(cfg, proj, cdr_proj, gs_default, gs_rest, composite_defs, panel_list)
    }, error = function(e) {
      err_msg <- if (nzchar(trimws(e$message))) e$message else class(e)[1]
      message("   [!] failed: ", err_msg)
      list(status = "Failed", n = NA_integer_, reason = err_msg, clinical = NULL)
    })

    cohort_outcomes[[proj]] <- list(status = res$status, n = res$n, reason = res$reason)
    if (!is.null(res$clinical)) {
      cohort_clinicals[[proj]] <- res$clinical
    }
  }

  summary_df <- dplyr::bind_rows(lapply(names(cohort_outcomes), function(p) {
    o <- cohort_outcomes[[p]]
    data.frame(Project = p,
               Status  = o$status,
               N       = if (is.null(o$n))      NA_integer_   else o$n,
               Reason  = if (is.null(o$reason)) NA_character_ else o$reason,
               stringsAsFactors = FALSE)
  }))

  n_ok   <- sum(grepl("^Success", summary_df$Status))
  n_skip <- sum(summary_df$Status == "Skipped")
  n_fail <- sum(summary_df$Status == "Failed")
  message(sprintf(
    "\n========== PAN-CANCER COHORT SUMMARY ==========\n  Succeeded : %d\n  Skipped   : %d\n  Failed    : %d\n  Total queued: %d",
    n_ok, n_skip, n_fail, nrow(summary_df)))

  if (n_skip + n_fail > 0) {
    not_ok <- summary_df[summary_df$Status %in% c("Skipped", "Failed"), , drop = FALSE]
    max_proj <- max(nchar(not_ok$Project))
    message("")
    for (j in seq_len(nrow(not_ok))) {
      r <- not_ok[j, ]
      message(sprintf("  [%-*s]  %-8s  %s", max_proj, r$Project, r$Status, r$Reason))
    }
  }
  message("================================================\n")

  if (length(cohort_clinicals) == 0) {
    stop("No per-cohort data produced successfully.")
  }

  master <- dplyr::bind_rows(cohort_clinicals)
  list(master = master, cohort_run_summary = summary_df)
}

# Apply batch correction (limma::removeBatchEffect) across all score columns.
.pancan_batch_correct <- function(master, score_cols) {
  score_matrix <- t(as.matrix(master[, score_cols, drop = FALSE]))
  corrected    <- limma::removeBatchEffect(score_matrix, batch = master$Project_ID)
  stopifnot(ncol(corrected) == nrow(master))
  corr_df      <- as.data.frame(t(corrected))
  corr_cols    <- paste0(score_cols, "_Corrected")
  colnames(corr_df) <- corr_cols
  master       <- dplyr::bind_cols(master, corr_df)
  message("Batch correction applied (corrected scores: pan-cancer aggregate only).")
  list(master = master, corr_cols = corr_cols)
}

# Run statistical tests (Wilcoxon, KM, Cox) for aggregate and per-cohort scores with FDR.
.pancan_stats_block <- function(master, score_cols, corr_cols, cfg) {
  rows <- list()
  # Pan-cancer aggregate uses CORRECTED scores
  for (sc in corr_cols) {
    rows[[length(rows) + 1]] <- run_all_stats(master, sc, "OS",  "PanCancer", "Pan-Cancer", sc, "Surv_3yr", cfg)
    rows[[length(rows) + 1]] <- run_all_stats(master, sc, "RFS", "PanCancer", "Pan-Cancer", sc, "Recurrence_3yr", cfg)
  }
  # Per-cohort uses UNCORRECTED scores
  for (proj in unique(master$Project_ID)) {
    cd <- master %>% dplyr::filter(Project_ID == proj)
    for (sc in score_cols) {
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "OS",  "PanCancer", proj, sc, "Surv_3yr", cfg)
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "RFS", "PanCancer", proj, sc, "Recurrence_3yr", cfg)
    }
  }

  stats_df <- do.call(rbind, rows) %>%
    dplyr::mutate(Family = ifelse(Project == "Pan-Cancer", "PanCancer", "PerCohort"))
  stats_df <- apply_fdr(stats_df, cfg, family = "pancan")

  cox_summary <- stats_df %>%
    dplyr::filter(Test == "Cox", !is.na(HR)) %>%
    dplyr::select(Project, Metric, Score, N, N_events, HR, HR_lower, HR_upper,
                  C_index, Raw_P, FDR_P, Is_Significant) %>%
    dplyr::arrange(Metric, Score, FDR_P)

  list(stats_df = stats_df, cox_summary = cox_summary)
}

#' Run pan-cancer survival analysis
#'
#' Evaluates the canonical DTP signature panel across TCGA solid tumor cohorts,
#' performing ssGSEA scoring, landmark survival analysis, batch-corrected aggregate
#' evaluation, and Benjamini-Hochberg FDR correction. Scores the default composite
#' signature pair (e.g. Up/Down) in an isolated ssGSEA call to preserve exact numerical
#' parity with legacy thesis results under GSVA within-sample normalization.
#'
#' @param cfg Config from [dtp_config()].
#' @param panel_tbl Canonical signature panel table from [load_signature_panel()].
#' @param out_root Output directory path (reserved for figure generation in future phases).
#' @return A named list containing:
#'   \item{stats_df}{FDR-corrected survival statistics table (Wilcoxon, KM, Cox) across cohorts and aggregate.}
#'   \item{master}{Full scored and joined clinical data frame with batch-corrected columns.}
#'   \item{cohort_summary}{Summary of TCGA-CDR inventory across cancer types.}
#'   \item{cohort_run_summary}{Execution status and patient counts per cohort.}
#'   \item{cox_summary}{Cox proportional hazards summary table sorted by FDR p-value.}
#' @export
run_pancan_survival <- function(cfg, panel_tbl, out_root = "output") {
  force(out_root)
  message("\n================== MODULE: PAN-CANCER SURVIVAL ==================")

  panel_list     <- panel_to_list(panel_tbl, id_type = cfg$id_type)
  composite_defs <- load_composite_defs(cfg)

  # Partition gene sets into default composite pair (e.g. Up/Down) and remaining signatures.
  # This guarantees isolated ssGSEA evaluation for the default pair to preserve
  # exact numerical scale parity with legacy results under GSVA's within-sample normalization.
  def_row <- composite_defs[composite_defs$Is_Default == TRUE, , drop = FALSE]
  default_sigs <- if (nrow(def_row) >= 1) {
    intersect(c(def_row$Positive_Signature[1], def_row$Negative_Signature[1]), names(panel_list))
  } else {
    character(0)
  }
  rest_sigs <- setdiff(names(panel_list), default_sigs)

  gs_default <- if (length(default_sigs) > 0) build_gene_sets(panel_list, which = default_sigs) else NULL
  gs_rest    <- if (length(rest_sigs) > 0)    build_gene_sets(panel_list, which = rest_sigs)    else NULL

  # ---- Phase 1: CDR + project inventory ----
  inv <- .pancan_cohort_inventory(cfg)

  # ---- Phase 2: per-cohort RNA-seq -> ssGSEA -> clinical ----
  run_res <- .pancan_run_all_cohorts(cfg, inv$valid_projects, inv$cdr_clinical,
                                     gs_default, gs_rest, composite_defs, panel_list)
  master  <- run_res$master

  # Identify dynamic score columns (all uncorrected score columns)
  score_cols <- grep(paste0(cfg$score_suffix, "$"), colnames(master), value = TRUE)

  # ---- Phase 3: batch correction across cancer types ----
  bc_res <- .pancan_batch_correct(master, score_cols)
  master <- bc_res$master

  # ---- Phase 4: statistics ----
  stats_res <- .pancan_stats_block(master, score_cols, bc_res$corr_cols, cfg)

  message("Pan-cancer survival module complete.")
  invisible(list(
    stats_df           = stats_res$stats_df,
    master             = master,
    cohort_summary     = inv$cohort_summary,
    cohort_run_summary = run_res$cohort_run_summary,
    cox_summary        = stats_res$cox_summary
  ))
}
