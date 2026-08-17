test_that("save_results_bundle and load_results_bundle round-trip with trimmed columns", {
  cfg <- dtp_config(overrides = list(
    score_suffix = "_ssGSEA",
    crc_modifiers = c("CMS", "PDS", "Stage_bin", "MSI_group")
  ))

  # Synthetic clinical data with extra raw columns that should be trimmed
  gse_clin <- data.frame(
    Sample_ID = c("S1", "S2", "S3", "S4"),
    TNM_stage = c("1", "2", "3", "4"),
    MMR = c("pMMR", "pMMR", "dMMR", "dMMR"),
    CMS = c("CMS1", "CMS2", "CMS3", "CMS4"),
    PDS = c("PDS1", "PDS2", "PDS3", "PDS1"),
    Chemo_adj = c("Y", "N", "Y", "N"),
    Surv_3yr = factor(c("Alive_3yr", "Dead_3yr", "Alive_3yr", "Dead_3yr")),
    Recurrence_3yr = factor(c("Recurrence-Free", "Recurred", "Recurrence-Free", "Recurred")),
    OS3Y_delay = c(36, 12, 36, 20),
    OS3Y_event = c(0, 1, 0, 1),
    RFS3Y_delay = c(36, 10, 36, 18),
    RFS3Y_event = c(0, 1, 0, 1),
    Up_ssGSEA = c(0.5, 1.2, -0.3, 0.8),
    Down_ssGSEA = c(-0.4, 0.1, 0.6, -0.2),
    Raw_unwanted_col1 = c("x", "y", "z", "w"),
    Raw_unwanted_col2 = 1:4,
    stringsAsFactors = FALSE
  )

  tcga_clin <- data.frame(
    ID = c("T1", "T2", "T3", "T4"),
    Stage = c("I", "II", "III", "IV"),
    paper_MSI_status = c("MSS", "MSS", "MSI-H", "MSI-L"),
    CMS = c("CMS2", "CMS1", "CMS4", "CMS3"),
    PDS = c("PDS2", "PDS1", "PDS3", "PDS2"),
    Surv_3yr = factor(c("Alive_3yr", "Alive_3yr", "Dead_3yr", "Dead_3yr")),
    Recurrence_3yr = factor(c("Recurrence-Free", "Recurrence-Free", "Recurred", "Recurred")),
    OS3Y_delay = c(36, 36, 15, 22),
    OS3Y_event = c(0, 0, 1, 1),
    RFS3Y_delay = c(36, 36, 14, 20),
    RFS3Y_event = c(0, 0, 1, 1),
    Up_ssGSEA = c(0.1, -0.2, 1.5, 0.9),
    Down_ssGSEA = c(0.3, 0.5, -0.6, -0.1),
    Raw_unwanted_col = c(10, 20, 30, 40),
    stringsAsFactors = FALSE
  )

  panel_tbl <- data.frame(
    Signature_ID = c("Up", "Down"),
    Signature_Name = c("DTP Up", "DTP Down"),
    Gene_ID = c("G1", "G2"),
    stringsAsFactors = FALSE
  )

  composite_defs <- data.frame(
    Composite_ID = "Composite",
    Display_Name = "DTP Composite",
    Positive_Signature = "Up",
    Negative_Signature = "Down",
    Is_Default = TRUE,
    stringsAsFactors = FALSE
  )

  stats_df <- data.frame(
    Dataset = "GSE39582", Project = "GSE_All_Patients", Test = "Wilcoxon",
    Metric = "OS", Score = "Up_ssGSEA", Raw_P = 0.05, FDR_P = 0.05,
    stringsAsFactors = FALSE
  )

  crc_result <- list(
    gse_clinical = gse_clin,
    tcga_clinical = tcga_clin,
    stats_df = stats_df
  )

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file), add = TRUE)

  saved_path <- save_results_bundle(cfg, crc_result, panel_tbl, composite_defs, path = tmp_file)
  expect_equal(saved_path, tmp_file)
  expect_true(file.exists(tmp_file))

  bundle <- load_results_bundle(tmp_file)

  # Check top-level names including new extension fields
  expected_names <- c(
    "gse", "tcga", "panel_tbl", "composite_defs", "stats_df", "cfg", "built_at", "panel_mtime",
    "pancan_stats_df", "pancan_table_tbl", "treated_stats_df", "mets_gsea_pn", "mets_gsea_pm_adj",
    "mets_pvm_cmp", "mets_liver_df", "mets_roc_summary", "mets_ssgsea_pvalues", "fig_dir", "excel_path"
  )
  expect_setequal(names(bundle), expected_names)

  # Check that optional fields default gracefully to NULL / default paths
  expect_null(bundle$pancan_stats_df)
  expect_null(bundle$pancan_table_tbl)
  expect_null(bundle$treated_stats_df)
  expect_null(bundle$mets_gsea_pn)
  expect_null(bundle$mets_gsea_pm_adj)
  expect_null(bundle$mets_pvm_cmp)
  expect_null(bundle$mets_liver_df)
  expect_null(bundle$mets_roc_summary)
  expect_null(bundle$mets_ssgsea_pvalues)
  expect_null(bundle$excel_path)
  expect_equal(bundle$fig_dir, "output/figures")

  # Check trimmed GSE columns
  expected_gse_cols <- c("Sample_ID", "Up_ssGSEA", "Down_ssGSEA", "Surv_3yr", "Recurrence_3yr",
                         "OS3Y_delay", "OS3Y_event", "RFS3Y_delay", "RFS3Y_event", "Chemo_adj",
                         "CMS", "PDS", "Stage_bin", "MSI_group")
  expect_setequal(colnames(bundle$gse), expected_gse_cols)
  expect_false("Raw_unwanted_col1" %in% colnames(bundle$gse))
  expect_false("Raw_unwanted_col2" %in% colnames(bundle$gse))
  expect_false("TNM_stage" %in% colnames(bundle$gse))
  expect_false("MMR" %in% colnames(bundle$gse))

  # Check trimmed TCGA columns
  expected_tcga_cols <- c("ID", "Up_ssGSEA", "Down_ssGSEA", "Surv_3yr", "Recurrence_3yr",
                          "OS3Y_delay", "OS3Y_event", "RFS3Y_delay", "RFS3Y_event",
                          "CMS", "PDS", "Stage_bin", "MSI_group")
  expect_setequal(colnames(bundle$tcga), expected_tcga_cols)
  expect_false("Raw_unwanted_col" %in% colnames(bundle$tcga))
  expect_false("Stage" %in% colnames(bundle$tcga))
  expect_false("paper_MSI_status" %in% colnames(bundle$tcga))

  # Missing bundle load error
  expect_error(load_results_bundle(tempfile(fileext = ".rds")), "Results bundle not found")
})

test_that("save_results_bundle extended call captures pan-cancer, mets, and artifact paths", {
  cfg <- dtp_config()
  panel_tbl <- data.frame(Signature_ID = "Up", Signature_Name = "DTP Up", Gene_ID = "G1", stringsAsFactors = FALSE)
  composite_defs <- data.frame(Composite_ID = "Comp", Display_Name = "DTP", Positive_Signature = "Up",
                               Negative_Signature = "None", Is_Default = TRUE, stringsAsFactors = FALSE)
  crc_result <- list(
    gse_clinical = data.frame(
      Sample_ID = "S1", TNM_stage = "1", MMR = "pMMR",
      Up_ssGSEA = 1.0, Surv_3yr = factor("Alive_3yr"), stringsAsFactors = FALSE
    ),
    tcga_clinical = data.frame(
      ID = "T1", Stage = "I", paper_MSI_status = "MSS",
      Up_ssGSEA = 0.5, Surv_3yr = factor("Alive_3yr"), stringsAsFactors = FALSE
    ),
    stats_df = data.frame(Test = "Wilcoxon", stringsAsFactors = FALSE)
  )

  pancan_res <- list(stats_df = data.frame(Project = "TCGA-BRCA", HR = 1.5))
  treated_res <- list(stats_df = data.frame(Project = "TCGA-BRCA", HR = 1.8))
  mets_res <- list(
    gsea_pn = data.frame(ID = "Sig1", NES = 1.5),
    gsea_pm_adj = data.frame(ID = "Sig1", NES = 1.2),
    pvm_cmp = data.frame(ID = "Sig1", NES_adj = 1.2),
    liver_df = data.frame(Sample = "S1", Score = 0.8),
    roc_summary = data.frame(ID = "Sig1", AUC = 0.85),
    ssgsea_pvalues = data.frame(Signature = "Sig1", p_adj_BH = 0.01)
  )
  pancan_tbl <- data.frame(Project = "TCGA-BRCA", Table = 3.5)

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file), add = TRUE)

  save_results_bundle(
    cfg = cfg,
    crc_result = crc_result,
    panel_tbl = panel_tbl,
    composite_defs = composite_defs,
    path = tmp_file,
    pancan_result = pancan_res,
    treated_result = treated_res,
    mets_result = mets_res,
    fig_dir = "custom/figures",
    excel_path = "custom/report.xlsx",
    pancan_table_tbl = pancan_tbl
  )

  bundle <- load_results_bundle(tmp_file)
  expect_equal(bundle$pancan_stats_df$Project, "TCGA-BRCA")
  expect_equal(bundle$pancan_table_tbl$Table, 3.5)
  expect_equal(bundle$treated_stats_df$HR, 1.8)
  expect_equal(bundle$mets_gsea_pn$NES, 1.5)
  expect_equal(bundle$mets_gsea_pm_adj$NES, 1.2)
  expect_equal(bundle$mets_pvm_cmp$NES_adj, 1.2)
  expect_equal(bundle$mets_liver_df$Score, 0.8)
  expect_equal(bundle$mets_roc_summary$AUC, 0.85)
  expect_equal(bundle$mets_ssgsea_pvalues$p_adj_BH, 0.01)
  expect_equal(bundle$fig_dir, "custom/figures")
  expect_equal(bundle$excel_path, "custom/report.xlsx")
})

test_that("check_bundle_staleness correctly compares modification timestamps", {
  tmp_panel <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_panel), add = TRUE)
  writeLines("test", tmp_panel)

  panel_mtime <- file.mtime(tmp_panel)

  # Bundle built BEFORE panel was modified (stale)
  bundle_old <- list(panel_mtime = panel_mtime - 10)
  expect_true(check_bundle_staleness(bundle_old, tmp_panel))

  # Bundle built AFTER panel was modified (fresh)
  bundle_fresh <- list(panel_mtime = panel_mtime + 10)
  expect_false(check_bundle_staleness(bundle_fresh, tmp_panel))

  # Non-existent panel file returns FALSE
  expect_false(check_bundle_staleness(bundle_old, tempfile(fileext = ".csv")))

  # Null/NA mtime returns FALSE
  expect_false(check_bundle_staleness(list(panel_mtime = NULL), tmp_panel))
  expect_false(check_bundle_staleness(list(panel_mtime = NA), tmp_panel))
})

test_that(".compute_live_score computes single and composite scores and handles errors", {
  df <- data.frame(
    SigA_ssGSEA = c(1.0, 2.5, -0.5, 0.0),
    SigB_ssGSEA = c(0.5, 1.0, 0.5, -1.0)
  )

  # Positive only (None)
  score_pos <- .compute_live_score(df, positive_id = "SigA", negative_id = "None")
  expect_equal(score_pos, df$SigA_ssGSEA)

  score_pos_null <- .compute_live_score(df, positive_id = "SigA", negative_id = NULL)
  expect_equal(score_pos_null, df$SigA_ssGSEA)

  # Difference: SigA - SigB
  score_diff <- .compute_live_score(df, positive_id = "SigA", negative_id = "SigB")
  expect_equal(score_diff, df$SigA_ssGSEA - df$SigB_ssGSEA)
  expect_equal(score_diff, c(0.5, 1.5, -1.0, 1.0))

  # Missing positive_id
  expect_error(.compute_live_score(df, positive_id = "UnknownSig"), "not found in clinical table")

  # Missing negative_id
  expect_error(.compute_live_score(df, positive_id = "SigA", negative_id = "UnknownSig"), "not found in clinical table")

  # Empty positive_id
  expect_error(.compute_live_score(df, positive_id = ""), "positive_id must be specified")
})
