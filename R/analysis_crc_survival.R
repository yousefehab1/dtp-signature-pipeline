# CRC survival analysis module (GSE39582 microarray + TCGA-COAD STAR TPM).
# Decomposed into small, testable internal helpers (.crc_*) and orchestrated
# by run_crc_survival().

.crc_load_gse39582 <- function(cfg) {
  cache_rds(cfg, "GSE39582_eset", function() {
    g <- GEOquery::getGEO("GSE39582", GSEMatrix = TRUE)
    g[[1]]
  })
}

.crc_score_gse <- function(cfg, eset, gs, composite_defs) {
  mat_sym <- prep_microarray_symbols(cfg, Biobase::exprs(eset))
  scores_gse <- run_ssgsea(mat_sym, gs, suffix = cfg$score_suffix)
  scores_gse <- add_all_composites(scores_gse, composite_defs, suffix = cfg$score_suffix)
  scores_gse$Sample_ID <- rownames(scores_gse)
  scores_gse
}

.crc_check_gse39582_consistency <- function(clinical_df) {
  characteristics_cols_used <- c(
    "characteristics_ch1.2",  "characteristics_ch1.3",  "characteristics_ch1.4",
    "characteristics_ch1.9",  "characteristics_ch1.11", "characteristics_ch1.12",
    "characteristics_ch1.13", "characteristics_ch1.14", "characteristics_ch1.15",
    "characteristics_ch1.18", "characteristics_ch1.22", "characteristics_ch1.26"
  )
  consistency_issues <- character(0)
  for (col in characteristics_cols_used) {
    if (!col %in% colnames(clinical_df)) next
    vals   <- clinical_df[[col]]
    vals   <- vals[!is.na(vals)]
    labels <- unique(trimws(sub(":.*$", "", vals)))
    if (length(labels) > 1) {
      consistency_issues <- c(
        consistency_issues,
        paste0(col, " has ", length(labels), " labels: ", paste(labels, collapse = ", "))
      )
    }
  }
  if (length(consistency_issues) > 0) {
    stop("GSE39582 characteristics_ch1.N field order inconsistent:\n",
         paste(consistency_issues, collapse = "\n"))
  }
  invisible(TRUE)
}

.crc_parse_numeric_field <- function(raw, field_name) {
  stripped <- stringr::str_remove(raw, "^.*: ")
  out      <- suppressWarnings(as.numeric(stripped))
  newly_na <- is.na(out) & !is.na(stripped)
  if (any(newly_na)) {
    message(sprintf("  %s: %d value(s) unparseable; e.g. %s", field_name, sum(newly_na),
                    paste(unique(stripped[newly_na])[seq_len(min(5, sum(newly_na)))], collapse = " | ")))
  }
  out
}

.crc_parse_gse39582_clinical <- function(eset, scores_gse, cfg) {
  clinical_df <- as.data.frame(Biobase::pData(eset))
  .crc_check_gse39582_consistency(clinical_df)

  clinical <- clinical_df %>%
    tibble::rownames_to_column(var = "Sample_ID") %>%
    dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
    dplyr::select(
      Sample_ID,
      Age_raw = characteristics_ch1.3,       Sex_raw = characteristics_ch1.2,
      Chemo_adj_raw = characteristics_ch1.9, TNM_stage_raw = characteristics_ch1.4,
      MMR_raw = characteristics_ch1.15,
      RFS_event_raw = characteristics_ch1.11, RFS_delay_raw = characteristics_ch1.12,
      OS_event_raw = characteristics_ch1.13,  OS_delay_raw = characteristics_ch1.14,
      TP53_raw = characteristics_ch1.18,      KRAS_raw = characteristics_ch1.22,
      BRAF_raw = characteristics_ch1.26
    ) %>%
    dplyr::mutate(
      Age       = .crc_parse_numeric_field(Age_raw, "Age"),
      RFS_event = .crc_parse_numeric_field(RFS_event_raw, "RFS_event"),
      RFS_delay = .crc_parse_numeric_field(RFS_delay_raw, "RFS_delay"),
      OS_event  = .crc_parse_numeric_field(OS_event_raw, "OS_event"),
      OS_delay  = .crc_parse_numeric_field(OS_delay_raw, "OS_delay"),
      Sex       = stringr::str_trim(stringr::str_remove(Sex_raw, "^.*: ")),
      Chemo_adj = stringr::str_trim(stringr::str_remove(Chemo_adj_raw, "^.*: ")),
      TNM_stage = stringr::str_trim(stringr::str_remove(TNM_stage_raw, "^.*: ")),
      MMR       = stringr::str_trim(stringr::str_remove(MMR_raw, "^.*: "))
    ) %>%
    dplyr::select(-dplyr::ends_with("_raw")) %>%
    dplyr::left_join(scores_gse, by = "Sample_ID") %>%
    dplyr::mutate(Recurrence_event = RFS_event, Recurrence_delay = RFS_delay) %>%
    derive_endpoints(
      os_event = "OS_event", os_time = "OS_delay",
      rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay",
      os_cut = cfg$os_cutpoint, rfs_cut = cfg$rfs_cutpoint
    )
  clinical
}

.crc_subtype_gse <- function(cfg, eset, clinical) {
  emat_gse_sym <- build_symbol_matrix(cfg, Biobase::exprs(eset), "microarray")
  cms_gse <- cache_rds(cfg, "cms_gse", function() call_cms(emat_gse_sym))
  pds_gse <- cache_rds(cfg, "pds_gse", function() call_pds(emat_gse_sym))
  clinical %>%
    attach_subtypes(cms_gse, clinical_key = "Sample_ID") %>%
    attach_subtypes(pds_gse, clinical_key = "Sample_ID")
}

.crc_build_gse_cohorts <- function(clinical) {
  gse_violin_cohorts <- list(
    "GSE_All_Patients"          = clinical,
    "GSE_Treated"               = clinical %>% dplyr::filter(Chemo_adj == "Y"),
    "GSE_Stage3_4_Treated"      = clinical %>% dplyr::filter(TNM_stage %in% c(3, 4), Chemo_adj == "Y"),
    "GSE_Stage3_Treated"        = clinical %>% dplyr::filter(TNM_stage %in% c(3),    Chemo_adj == "Y"),
    "GSE_Stage3_4_Treated_MSS"  = clinical %>% dplyr::filter(TNM_stage %in% c(3, 4), Chemo_adj == "Y", MMR == "pMMR")
  )
  gse_km_cohorts <- list(
    "GSE_All_Treated"           = clinical %>% dplyr::filter(Chemo_adj == "Y"),
    "GSE_All_MSS"               = clinical %>% dplyr::filter(MMR == "pMMR"),
    "GSE_Stage2_Untreated_MSS"  = clinical %>% dplyr::filter(MMR == "pMMR", TNM_stage == 2, Chemo_adj == "N"),
    "GSE_Stage3_Treated_MSS"    = clinical %>% dplyr::filter(MMR == "pMMR", TNM_stage == 3, Chemo_adj == "Y"),
    "GSE_Stage34_Treated_MSS"   = clinical %>% dplyr::filter(MMR == "pMMR", TNM_stage %in% c(3, 4), Chemo_adj == "Y")
  )
  gse_extra_cohorts <- list(
    "GSE_Stage1"  = clinical %>% dplyr::filter(TNM_stage == 1),
    "GSE_Stage2"  = clinical %>% dplyr::filter(TNM_stage == 2),
    "GSE_Stage3"  = clinical %>% dplyr::filter(TNM_stage == 3),
    "GSE_Stage4"  = clinical %>% dplyr::filter(TNM_stage == 4),
    "GSE_All_MSI" = clinical %>% dplyr::filter(MMR == "dMMR")
  )
  c(gse_violin_cohorts, gse_km_cohorts, gse_extra_cohorts,
    subtype_cohorts(clinical, "CMS", "GSE"),
    subtype_cohorts(clinical, "PDS", "GSE"))
}

.crc_load_tcga_coad <- function(cfg) {
  query <- TCGAbiolinks::GDCquery(
    project = "TCGA-COAD",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type   = "Primary Tumor"
  )
  mat_tcga <- cache_rds(cfg, paste0("TCGA_COAD_tpm_", cfg$id_type), function() {
    se <- cache_rds(cfg, "TCGA_COAD_se", function() {
      TCGAbiolinks::GDCdownload(query)
      TCGAbiolinks::GDCprepare(query)
    })
    m  <- prep_tcga_tpm(cfg, se)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    saveRDS(cd, file.path(cfg$cache_dir, "TCGA_COAD_coldata.rds"))
    m
  })
  coad_coldata <- cache_rds(cfg, "TCGA_COAD_coldata", function() stop("colData cache missing"))
  list(mat_tcga = mat_tcga, coad_coldata = coad_coldata)
}

.crc_score_tcga <- function(cfg, mat_tcga, gs, composite_defs) {
  run_ssgsea(mat_tcga, gs, suffix = cfg$score_suffix) %>%
    add_all_composites(composite_defs, suffix = cfg$score_suffix) %>%
    tibble::rownames_to_column("Full_Barcode") %>%
    dplyr::arrange(Full_Barcode) %>%
    dplyr::mutate(ID = tcga_patient_id(Full_Barcode)) %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::select(-Full_Barcode)
}

.crc_load_tcga_clinical <- function(cfg, coad_coldata) {
  cdr <- load_tcga_cdr(cfg) %>% dplyr::filter(Project_ID == "TCGA-COAD")

  has_value <- function(col, val) vapply(col, function(x) val %in% x, logical(1))
  treat <- cache_rds(cfg, "TCGA_COAD_clinic_treatment", function() {
    TCGAbiolinks::GDCquery_clinic(project = "TCGA-COAD", type = "clinical") %>%
      dplyr::mutate(
        ID = trimws(toupper(gsub("\\.", "-", submitter_id))),
        Treatment_Status = dplyr::case_when(
          has_value(treatments_pharmaceutical_treatment_or_therapy, "yes") |
            has_value(treatments_radiation_treatment_or_therapy, "yes") ~ "Treated",
          has_value(treatments_pharmaceutical_treatment_or_therapy, "no") &
            has_value(treatments_radiation_treatment_or_therapy, "no")  ~ "Not Treated",
          TRUE ~ "Unknown/Not Reported")) %>%
      dplyr::select(ID, Treatment_Status) %>%
      dplyr::distinct(ID, .keep_all = TRUE)
  })

  msi <- coad_coldata %>%
    dplyr::mutate(ID = substr(barcode, 1, 12)) %>%
    dplyr::arrange(barcode) %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::select(ID, dplyr::any_of("paper_MSI_status"))

  list(cdr = cdr, treat = treat, msi = msi)
}

.crc_join_tcga_clinical <- function(scores_tcga, cdr, treat, msi, cfg) {
  scores_tcga %>%
    dplyr::inner_join(
      cdr %>% dplyr::select(ID = Patient_ID, Stage, OS_event, OS_months, DFS_event, DFS_months),
      by = "ID"
    ) %>%
    dplyr::left_join(treat, by = "ID") %>%
    dplyr::left_join(msi,   by = "ID") %>%
    dplyr::mutate(OS_delay = OS_months, Recurrence_event = DFS_event, Recurrence_delay = DFS_months) %>%
    derive_endpoints(
      os_event = "OS_event", os_time = "OS_delay",
      rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay",
      os_cut = cfg$os_cutpoint, rfs_cut = cfg$rfs_cutpoint
    )
}

.crc_subtype_tcga <- function(cfg, clinical_tcga) {
  se_coad <- cache_rds(cfg, "TCGA_COAD_se", function() stop("SE cache missing"))
  emat_tcga_sym <- build_symbol_matrix(cfg, se_coad, "tcga")
  cms_tcga <- cache_rds(cfg, "cms_tcga", function() call_cms(emat_tcga_sym))
  pds_tcga <- cache_rds(cfg, "pds_tcga", function() call_pds(emat_tcga_sym))
  clinical_tcga %>%
    attach_subtypes(cms_tcga, clinical_key = "ID", key_fun = tcga_patient_id) %>%
    attach_subtypes(pds_tcga, clinical_key = "ID", key_fun = tcga_patient_id)
}

.crc_build_tcga_cohorts <- function(clinical_tcga) {
  has_msi <- "paper_MSI_status" %in% colnames(clinical_tcga)
  tcga_cohorts <- list(
    "TCGA_All_Patients"     = clinical_tcga,
    "TCGA_Untreated"        = clinical_tcga %>% dplyr::filter(Treatment_Status == "Not Treated"),
    "TCGA_Treated"          = clinical_tcga %>% dplyr::filter(Treatment_Status == "Treated"),
    "TCGA_Stage1_Untreated" = clinical_tcga %>% dplyr::filter(Stage == "I",  Treatment_Status == "Not Treated"),
    "TCGA_Stage2_Untreated" = clinical_tcga %>% dplyr::filter(Stage == "II", Treatment_Status == "Not Treated"),
    "TCGA_Stage1"           = clinical_tcga %>% dplyr::filter(Stage == "I"),
    "TCGA_Stage2"           = clinical_tcga %>% dplyr::filter(Stage == "II"),
    "TCGA_Stage3"           = clinical_tcga %>% dplyr::filter(Stage == "III"),
    "TCGA_Stage4"           = clinical_tcga %>% dplyr::filter(Stage == "IV")
  )
  if (has_msi) {
    tcga_cohorts[["TCGA_All_MSS"]] <- clinical_tcga %>% dplyr::filter(paper_MSI_status == "MSS")
    tcga_cohorts[["TCGA_All_MSI"]] <- clinical_tcga %>% dplyr::filter(paper_MSI_status %in% c("MSI-H", "MSI-L"))
  }
  c(tcga_cohorts,
    subtype_cohorts(clinical_tcga, "CMS", "TCGA"),
    subtype_cohorts(clinical_tcga, "PDS", "TCGA"))
}

.crc_subtype_characterization <- function(subtype_data, score_cols, cfg) {
  sub_rows <- list()
  for (ds in names(subtype_data)) {
    cd <- subtype_data[[ds]]
    for (ax in c("CMS", "PDS")) {
      if (!ax %in% colnames(cd)) next
      for (sc in score_cols) {
        ks <- get_kruskal_stats(cd, sc, ax, cfg)
        sub_rows[[length(sub_rows) + 1]] <- data.frame(
          Dataset = ds, Subtype_Axis = ax, Score = sc,
          Raw_P = ks$p, Eps2 = ks$eps2, N = ks$n, K = ks$k,
          Is_Testable = !is.na(ks$p), stringsAsFactors = FALSE
        )
      }
    }
  }
  subtype_stats <- if (length(sub_rows) > 0) do.call(rbind, sub_rows) else NULL
  if (!is.null(subtype_stats)) {
    subtype_stats <- subtype_stats %>%
      dplyr::group_by(Dataset, Subtype_Axis) %>%
      dplyr::mutate(FDR_P = stats::p.adjust(Raw_P, method = "BH")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Is_Significant = !is.na(FDR_P) & FDR_P < 0.05)
  }
  subtype_stats
}

.crc_confounding_block <- function(cox_data, score_cols, cfg) {
  core_scores <- intersect(paste0(c("Up", "Down", "Composite"), cfg$score_suffix), score_cols)
  adj_specs   <- list(Clinicopath = c("Stage_bin", "MSI_group"), CMS_adjusted = "CMS", PDS_adjusted = "PDS")

  adj_rows <- list()
  int_rows <- list()
  level_rows <- list()

  for (ds in names(cox_data)) {
    cd <- cox_data[[ds]]
    for (metric in c("OS", "RFS")) {
      mc <- metric_cols(metric)
      tcol <- mc[["t"]]
      ecol <- mc[["e"]]
      for (sc in core_scores) {
        for (mlab in names(adj_specs)) {
          covs <- intersect(adj_specs[[mlab]], colnames(cd))
          a <- if (length(covs) == 0) NULL else get_cox_adjusted(cd, sc, tcol, ecol, covs, cfg)
          adj_rows[[length(adj_rows) + 1]] <- data.frame(
            Dataset = ds, Metric = metric, Score = sc, Model = mlab,
            HR = if (is.null(a)) NA else a$HR,
            HR_lower = if (is.null(a)) NA else a$HR_lower, HR_upper = if (is.null(a)) NA else a$HR_upper,
            Raw_P = if (is.null(a)) NA else a$P, C_index = if (is.null(a)) NA else a$C_index,
            LRT_score_P = if (is.null(a)) NA else a$LRT_P, PH_P = if (is.null(a)) NA else a$PH_P,
            Unadj_HR = if (is.null(a)) NA else a$Unadj_HR,
            Delta_logHR = if (is.null(a)) NA else a$Delta_logHR,
            Delta_logHR_pct = if (is.null(a)) NA else a$Delta_logHR_pct,
            N = if (is.null(a)) NA else a$N, N_events = if (is.null(a)) NA else a$N_events,
            Is_Testable = !is.null(a), stringsAsFactors = FALSE
          )
        }
        for (modf in intersect(cfg$crc_modifiers, colnames(cd))) {
          it <- get_cox_interaction(cd, sc, tcol, ecol, modf, cfg)
          int_rows[[length(int_rows) + 1]] <- data.frame(
            Dataset = ds, Metric = metric, Score = sc, Modifier = modf,
            Raw_P = if (is.null(it)) NA else it$Interaction_P, K = if (is.null(it)) NA else it$K,
            N = if (is.null(it)) NA else it$N, N_events = if (is.null(it)) NA else it$N_events,
            Is_Testable = !is.null(it), stringsAsFactors = FALSE
          )
          if (!is.null(it) && !is.null(it$per_level)) {
            pl <- it$per_level
            pl$Dataset <- ds
            pl$Metric <- metric
            pl$Score <- sc
            pl$Modifier <- modf
            level_rows[[length(level_rows) + 1]] <- pl
          }
        }
      }
    }
  }
  adj_df <- if (length(adj_rows) > 0) do.call(rbind, adj_rows) else NULL
  if (!is.null(adj_df)) {
    adj_df <- adj_df %>%
      dplyr::group_by(Dataset, Metric, Model) %>%
      dplyr::mutate(FDR_P = stats::p.adjust(LRT_score_P, method = "BH")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Is_Significant = !is.na(FDR_P) & FDR_P < 0.05)
  }
  int_df <- if (length(int_rows) > 0) do.call(rbind, int_rows) else NULL
  if (!is.null(int_df)) {
    int_df <- int_df %>%
      dplyr::group_by(Dataset, Metric, Modifier) %>%
      dplyr::mutate(FDR_P = stats::p.adjust(Raw_P, method = "BH")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Is_Significant = !is.na(FDR_P) & FDR_P < 0.05)
  }
  level_df <- if (length(level_rows) > 0) do.call(rbind, level_rows) else NULL

  list(adj_df = adj_df, int_df = int_df, level_df = level_df)
}

#' Run colorectal cancer survival analysis
#'
#' Evaluates the canonical DTP signature panel in GSE39582 (Affymetrix microarray)
#' and TCGA-COAD (STAR TPM) cohorts across multiple sub-cohorts, molecular subtypes
#' (CMS, PDS), and confounder-adjusted models.
#'
#' @param cfg Config from [dtp_config()].
#' @param panel_tbl Canonical signature panel table from [load_signature_panel()].
#' @param out_root Output directory path (reserved for figure generation in future phases).
#' @return A named list containing:
#'   \item{stats_df}{Main survival statistics table (Wilcoxon, KM, Cox) with cohort-level FDR.}
#'   \item{adj_df}{Multivariable adjusted Cox models with model-level FDR.}
#'   \item{int_df}{Cox interaction models with modifier-level FDR.}
#'   \item{level_df}{Stratified Cox statistics per modifier level.}
#'   \item{subtype_stats}{Kruskal-Wallis test statistics across CMS and PDS subtypes.}
#'   \item{gse_clinical}{Annotated clinical data frame for GSE39582.}
#'   \item{tcga_clinical}{Annotated clinical data frame for TCGA-COAD.}
#' @export
run_crc_survival <- function(cfg, panel_tbl, out_root = "output") {
  # out_root is accepted for interface compatibility with future plotting phases
  force(out_root)

  panel_list     <- panel_to_list(panel_tbl, id_type = cfg$id_type)
  gs             <- build_gene_sets(panel_list, which = NULL)
  composite_defs <- load_composite_defs(cfg)
  all_stats      <- list()

  # ---- PART A: GSE39582 (Affymetrix microarray) ----
  eset        <- .crc_load_gse39582(cfg)
  scores_gse  <- .crc_score_gse(cfg, eset, gs, composite_defs)
  score_cols  <- setdiff(colnames(scores_gse), "Sample_ID")
  clinical    <- .crc_parse_gse39582_clinical(eset, scores_gse, cfg)
  clinical    <- .crc_subtype_gse(cfg, eset, clinical)
  gse_cohorts <- .crc_build_gse_cohorts(clinical)
  all_stats[["gse"]] <- run_survival_block(gse_cohorts, "GSE39582", score_cols, cfg)

  # ---- PART B: TCGA-COAD (STAR TPM, primary-only) ----
  tcga_load     <- .crc_load_tcga_coad(cfg)
  mat_tcga      <- tcga_load$mat_tcga
  coad_coldata  <- tcga_load$coad_coldata
  scores_tcga   <- .crc_score_tcga(cfg, mat_tcga, gs, composite_defs)
  tcga_clin_raw <- .crc_load_tcga_clinical(cfg, coad_coldata)
  clinical_tcga <- .crc_join_tcga_clinical(scores_tcga, tcga_clin_raw$cdr,
                                           tcga_clin_raw$treat, tcga_clin_raw$msi, cfg)
  clinical_tcga <- .crc_subtype_tcga(cfg, clinical_tcga)
  tcga_cohorts  <- .crc_build_tcga_cohorts(clinical_tcga)
  all_stats[["tcga"]] <- run_survival_block(tcga_cohorts, "TCGA-COAD", score_cols, cfg)

  # ---- FDR (cohort family for main survival block) ----
  stats_df <- apply_fdr(do.call(rbind, all_stats), cfg, family = "cohort")

  # ---- Score-across-subtype characterisation (Kruskal-Wallis) ----
  subtype_data  <- list("GSE39582" = clinical, "TCGA-COAD" = clinical_tcga)
  subtype_stats <- .crc_subtype_characterization(subtype_data, score_cols, cfg)

  # ---- Confounding + effect-modification (adjusted Cox, interaction Cox) ----
  cox_data <- list(
    "GSE39582"           = harmonize_crc_modifiers(clinical, "GSE39582"),
    "GSE39582 (treated)" = harmonize_crc_modifiers(dplyr::filter(clinical, Chemo_adj == "Y"), "GSE39582"),
    "TCGA-COAD"          = harmonize_crc_modifiers(clinical_tcga, "TCGA-COAD")
  )
  conf_res <- .crc_confounding_block(cox_data, score_cols, cfg)

  invisible(list(
    stats_df      = stats_df,
    adj_df        = conf_res$adj_df,
    int_df        = conf_res$int_df,
    level_df      = conf_res$level_df,
    subtype_stats = subtype_stats,
    gse_clinical  = clinical,
    tcga_clinical = clinical_tcga
  ))
}
