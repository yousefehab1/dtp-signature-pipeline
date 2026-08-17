test_that(".gsea_running_data and plot_gsea_enrichment operate correctly on gseaResult", {
  set.seed(42)
  genes <- paste0("G", 1:50)
  gene_list <- sort(stats::setNames(rnorm(50), genes), decreasing = TRUE)
  t2g <- data.frame(
    Term = c(rep("SetA", 10), rep("SetB", 10)),
    Gene = c(paste0("G", 1:10), paste0("G", 41:50)),
    stringsAsFactors = FALSE
  )
  gsea_res <- clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = t2g,
    minGSSize = 3,
    maxGSSize = 500,
    pvalueCutoff = 1,
    BPPARAM = BiocParallel::SerialParam()
  )

  # 1. Internal running data helper
  df_run <- .gsea_running_data(gsea_res, "SetA")
  expect_s3_class(df_run, "data.frame")
  expect_true(all(c("x", "runningScore", "position", "geneList") %in% colnames(df_run)))

  # Missing ID produces informative error
  expect_error(.gsea_running_data(gsea_res, "NonExistentSet"), "is not present in the GSEA result")

  # 2. Single gene set enrichment plot
  p_single <- plot_gsea_enrichment(gsea_res, "SetA", title = "Single Set")
  expect_true(inherits(p_single, "patchwork") || inherits(p_single, "ggplot"))

  # 3. Multi-signature overlay
  p_multi <- plot_gsea_enrichment(gsea_res, c("SetA", "SetB"), title = "Multi Set")
  expect_true(inherits(p_multi, "patchwork") || inherits(p_multi, "ggplot"))
})

test_that(".mets_nes_heatmap produces ggplot from synthetic GSEA data frame", {
  df_synth <- data.frame(
    ID = c("Sig1", "Sig2", "Sig3"),
    NES = c(1.85, -1.20, 0.45),
    pvalue = c(0.002, 0.04, 0.60),
    p.adjust = c(0.008, 0.08, 0.60),
    stringsAsFactors = FALSE
  )

  p <- .mets_nes_heatmap(
    gsea_df = df_synth,
    contrast_lab = "P / M",
    pos_label = "Metastasis",
    sig_order = c("Sig1", "Sig2", "Sig3"),
    lim = 2.0
  )

  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$Signature), c("Sig3", "Sig2", "Sig1"))
  expect_true(any(grepl("1.85", p$data$lab)))
  expect_true(any(grepl("q=0.008", p$data$lab)))
})

test_that("build_mets_gsea_composite handles shared signature order and missing signatures gracefully", {
  set.seed(42)
  genes <- paste0("G", 1:60)
  gene_list_pn <- sort(stats::setNames(rnorm(60, mean = 0.5), genes), decreasing = TRUE)
  gene_list_pm <- sort(stats::setNames(rnorm(60, mean = -0.5), genes), decreasing = TRUE)

  t2g_pn <- data.frame(
    Term = c(rep("SigShared", 10), rep("SigPNOnly", 10)),
    Gene = c(paste0("G", 1:10), paste0("G", 11:20)),
    stringsAsFactors = FALSE
  )
  t2g_pm <- data.frame(
    Term = c(rep("SigShared", 10), rep("SigPMOnly", 10)),
    Gene = c(paste0("G", 1:10), paste0("G", 41:50)),
    stringsAsFactors = FALSE
  )

  gsea_pn_obj <- clusterProfiler::GSEA(
    geneList = gene_list_pn, TERM2GENE = t2g_pn,
    minGSSize = 3, maxGSSize = 500, pvalueCutoff = 1,
    BPPARAM = BiocParallel::SerialParam()
  )
  gsea_pm_obj <- clusterProfiler::GSEA(
    geneList = gene_list_pm, TERM2GENE = t2g_pm,
    minGSSize = 3, maxGSSize = 500, pvalueCutoff = 1,
    BPPARAM = BiocParallel::SerialParam()
  )

  mets_res <- list(
    gsea_pn = as.data.frame(gsea_pn_obj),
    gsea_pm_adj = as.data.frame(gsea_pm_obj),
    gsea_pn_obj = gsea_pn_obj,
    gsea_pm_adj_obj = gsea_pm_obj
  )

  # Check that message is emitted when signatures are absent from one contrast,
  # but the composite still builds successfully without crashing.
  expect_message(
    p_comp <- build_mets_gsea_composite(mets_res),
    "\\[mets GSEA\\] signature\\(s\\) absent from one contrast"
  )
  expect_true(inherits(p_comp, "patchwork") || inherits(p_comp, "ggplot"))
})

test_that("build_mets_liver_composite creates 2-panel composite figure", {
  liver_df <- data.frame(
    GSM = paste0("GSM", 1:18),
    Patient = rep(paste0("P", 1:6), each = 3),
    Tissue = rep(c("Normal", "Primary", "Metastasis"), 6),
    Liver = c(rnorm(6, -0.2, 0.1), rnorm(6, -0.1, 0.1), rnorm(6, 0.5, 0.1)),
    stringsAsFactors = FALSE
  )

  pvm_cmp <- data.frame(
    ID = c("Up", "Down", "CBC"),
    NES_unadjusted = c(1.5, -1.2, 0.8),
    FDR_unadjusted = c(0.01, 0.04, 0.15),
    NES_adjusted = c(1.8, -1.4, 0.5),
    FDR_adjusted = c(0.005, 0.02, 0.30),
    dNES = c(0.3, -0.2, -0.3),
    Sig_unadjusted = c(TRUE, TRUE, FALSE),
    Sig_adjusted = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  mets_res <- list(liver_df = liver_df, pvm_cmp = pvm_cmp)

  p_liver <- build_mets_liver_composite(mets_res)
  expect_true(inherits(p_liver, "patchwork") || inherits(p_liver, "ggplot"))
})

test_that("build_mets_composites saves figures and isolates errors via tryCatch", {
  set.seed(42)
  genes <- paste0("G", 1:50)
  gene_list <- sort(stats::setNames(rnorm(50), genes), decreasing = TRUE)
  t2g <- data.frame(Term = rep("SigA", 10), Gene = paste0("G", 1:10), stringsAsFactors = FALSE)
  gsea_obj <- clusterProfiler::GSEA(
    geneList = gene_list, TERM2GENE = t2g,
    minGSSize = 3, maxGSSize = 500, pvalueCutoff = 1,
    BPPARAM = BiocParallel::SerialParam()
  )

  liver_df <- data.frame(
    GSM = paste0("GSM", 1:6),
    Patient = rep(c("P1", "P2"), each = 3),
    Tissue = rep(c("Normal", "Primary", "Metastasis"), 2),
    Liver = rnorm(6),
    stringsAsFactors = FALSE
  )
  pvm_cmp <- data.frame(
    ID = "SigA",
    NES_unadjusted = 1.2,
    FDR_unadjusted = 0.04,
    NES_adjusted = 1.5,
    FDR_adjusted = 0.01,
    stringsAsFactors = FALSE
  )

  mets_res_valid <- list(
    gsea_pn = as.data.frame(gsea_obj),
    gsea_pm_adj = as.data.frame(gsea_obj),
    gsea_pn_obj = gsea_obj,
    gsea_pm_adj_obj = gsea_obj,
    liver_df = liver_df,
    pvm_cmp = pvm_cmp
  )

  tmp_out <- tempfile("mets_fig_test_")
  dir.create(tmp_out, recursive = TRUE)
  on.exit(unlink(tmp_out, recursive = TRUE), add = TRUE)

  res_valid <- build_mets_composites(mets_res_valid, out_dir = tmp_out)
  expect_type(res_valid$status, "list")
  expect_true(all(unlist(res_valid$status)))

  expect_true(file.exists(file.path(tmp_out, "Fig_Mets_GSEA_composite.pdf")))
  expect_true(file.exists(file.path(tmp_out, "Fig_Mets_GSEA_composite.png")))
  expect_true(file.exists(file.path(tmp_out, "Fig_Mets_LiverPurity_composite.pdf")))
  expect_true(file.exists(file.path(tmp_out, "Fig_Mets_LiverPurity_composite.png")))

  # Error isolation: deliberately break liver inputs
  mets_res_broken <- mets_res_valid
  mets_res_broken$liver_df <- NULL

  tmp_broken <- tempfile("mets_fig_broken_")
  dir.create(tmp_broken, recursive = TRUE)
  on.exit(unlink(tmp_broken, recursive = TRUE), add = TRUE)

  suppressWarnings({
    res_broken <- build_mets_composites(mets_res_broken, out_dir = tmp_broken)
  })
  expect_true(res_broken$status[["GSEA composite"]])
  expect_false(res_broken$status[["Liver purity composite"]])
  expect_true(file.exists(file.path(tmp_broken, "Fig_Mets_GSEA_composite.pdf")))
})

test_that("build_excel_workbook creates workbook with all 7 metastasis sheets", {
  cfg <- dtp_config()
  tmp_xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp_xlsx), add = TRUE)

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

  mets_for_excel <- list(
    gsea_pn = data.frame(ID = "Up", NES = 1.5, p.adjust = 0.01),
    gsea_pm_adj = data.frame(ID = "Up", NES = 1.8, p.adjust = 0.005),
    pvm_cmp = data.frame(ID = "Up", NES_unadj = 1.5, NES_adj = 1.8),
    liver_df = data.frame(GSM = "GSM1", Tissue = "Normal", Liver = 0.1),
    roc_summary = data.frame(Signature = "Up", AUC = 0.85),
    scores_df = data.frame(GSM = "GSM1", Up_ssGSEA = 0.5),
    ssgsea_pvalues = data.frame(Signature = "Up", p_raw = 0.01, p_adj_BH = 0.02)
  )

  build_excel_workbook(
    cfg, tmp_xlsx, panel_tbl, composite_defs,
    crc = list(), pancan = list(), mets = mets_for_excel
  )

  expect_true(file.exists(tmp_xlsx))
  sheets <- openxlsx::getSheetNames(tmp_xlsx)
  expected_mets_sheets <- c(
    "Mets_GSEA_NormalVsPrimary",
    "Mets_GSEA_PvM_Adjusted",
    "Mets_GSEA_Adj_vs_Unadj",
    "Mets_LiverPurity_Scores",
    "Mets_ROC_AUC",
    "Mets_SingleSample_Scores",
    "Mets_SingleSample_Wilcoxon"
  )
  for (s in expected_mets_sheets) {
    expect_true(s %in% sheets)
  }
})
