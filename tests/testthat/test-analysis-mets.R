test_that(".mets_paired_wilcox calculates paired Wilcoxon p-values correctly", {
  set.seed(123)
  df <- data.frame(
    Patient = paste0("P", rep(1:10, 2)),
    Tissue  = rep(c("Primary", "Metastasis"), each = 10),
    Score   = c(1:10 + rnorm(10, 0, 0.05), (1:10) + 2 + rnorm(10, 0, 0.05)),
    stringsAsFactors = FALSE
  )
  pval <- .mets_paired_wilcox(df, "Score", "Tissue", "Patient", "Primary", "Metastasis")
  expect_true(is.numeric(pval))
  expect_true(pval < 0.01)
})

test_that(".mets_compute_pca returns tidy coordinates and variance explained", {
  set.seed(42)
  n_genes <- 20
  n_samples <- 6
  mat <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(mat) <- paste0("G", 1:n_genes)
  colnames(mat) <- paste0("S", 1:n_samples)

  cd <- data.frame(
    GSM     = paste0("S", 1:n_samples),
    Tissue  = rep(c("Normal", "Primary", "Metastasis"), 2),
    Patient = rep(c("P1", "P2"), each = 3),
    stringsAsFactors = FALSE
  )

  pca_res <- .mets_compute_pca(mat, cd)
  expect_s3_class(pca_res, "tbl_df")
  expect_equal(nrow(pca_res), n_samples)
  expect_named(pca_res, c("GSM", "Tissue", "Patient", "PC1", "PC2", "PC1_var_pct", "PC2_var_pct"))
  expect_true(all(pca_res$PC1_var_pct >= 0 & pca_res$PC1_var_pct <= 100))
})

test_that("paired limma model and liver purity status classification operate correctly", {
  set.seed(42)
  patients <- paste0("P", 1:6)
  tissues <- c("Normal", "Primary", "Metastasis")
  cd <- expand.grid(Tissue = tissues, Patient = patients, stringsAsFactors = FALSE)
  cd$GSM <- paste0("GSM_", seq_len(nrow(cd)))
  rownames(cd) <- cd$GSM
  cd$Tissue <- factor(cd$Tissue, levels = c("Normal", "Primary", "Metastasis"))
  cd$Patient <- factor(cd$Patient)
  cd$Tissue_for_DE <- stats::relevel(cd$Tissue, ref = "Primary")

  # Synthetic expression matrix: 10 genes x 18 samples
  genes <- paste0("Gene_", 1:10)
  mat <- matrix(rnorm(10 * 18, mean = 5, sd = 1), nrow = 10, ncol = 18,
                dimnames = list(genes, cd$GSM))

  # Add artificial signal
  # Gene_1: robustly up in Metastasis
  mat["Gene_1", cd$Tissue == "Metastasis"] <- mat["Gene_1", cd$Tissue == "Metastasis"] + 4
  # Gene_2: up in Normal
  mat["Gene_2", cd$Tissue == "Normal"] <- mat["Gene_2", cd$Tissue == "Normal"] + 4

  design <- stats::model.matrix(~ Patient + Tissue_for_DE, data = cd)
  fit <- limma::eBayes(limma::lmFit(mat, design))
  res <- limma::topTable(fit, coef = "Tissue_for_DEMetastasis", number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(stat = t, padj = adj.P.Val, pvalue = P.Value)

  expect_named(res, c("gene", "logFC", "AveExpr", "stat", "pvalue", "padj", "B"))
  expect_equal(nrow(res), 10)

  # Check 4-way Status classification logic
  test_comp <- data.frame(
    gene       = paste0("G", 1:4),
    padj_unadj = c(0.01, 0.10, 0.01, 0.50),
    padj_adj   = c(0.10, 0.01, 0.02, 0.60),
    Is_Liver_Gene = c(TRUE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Status = dplyr::case_when(
        padj_unadj < 0.05 & padj_adj >= 0.05 ~ "Lost significance after adjustment",
        padj_unadj >= 0.05 & padj_adj < 0.05 ~ "Gained significance after adjustment",
        padj_unadj < 0.05 & padj_adj < 0.05  ~ "Robust to adjustment",
        TRUE                                  ~ "Not significant either way"
      )
    )

  expect_equal(test_comp$Status, c(
    "Lost significance after adjustment",
    "Gained significance after adjustment",
    "Robust to adjustment",
    "Not significant either way"
  ))
})

test_that("single-sample ROC AUC and Wilcoxon handle presence and absence of Composite", {
  set.seed(42)
  patients <- paste0("P", 1:10)
  cd2 <- expand.grid(Tissue = c("Primary", "Metastasis"), Patient = patients, stringsAsFactors = FALSE)
  cd2$GSM <- paste0("GSM_", seq_len(nrow(cd2)))

  # Scores with Up, Down, and Composite
  scores_with_comp <- data.frame(
    GSM             = cd2$GSM,
    Up_ssGSEA       = ifelse(cd2$Tissue == "Metastasis", rnorm(20, 0.6, 0.1), rnorm(20, 0.2, 0.1)),
    Down_ssGSEA     = rnorm(20, 0.4, 0.1),
    Composite_ssGSEA= ifelse(cd2$Tissue == "Metastasis", rnorm(20, 0.2, 0.1), rnorm(20, -0.2, 0.1)),
    Tissue          = factor(cd2$Tissue, levels = c("Primary", "Metastasis")),
    Patient         = factor(cd2$Patient),
    stringsAsFactors = FALSE
  )

  # Run ROC & Wilcoxon for with-composite
  sig_names <- c("Up", "Down", "Composite")
  roc_summary <- list()
  for (gs in sig_names) {
    col <- paste0(gs, "_ssGSEA")
    v <- scores_with_comp[[col]]
    mp <- mean(v[scores_with_comp$Tissue == "Primary"])
    mm <- mean(v[scores_with_comp$Tissue == "Metastasis"])
    pos <- if (mm >= mp) "Metastasis" else "Primary"
    neg <- if (mm >= mp) "Primary" else "Metastasis"
    pred <- ROCR::prediction(v, scores_with_comp$Tissue, label.ordering = c(neg, pos))
    auc  <- ROCR::performance(pred, "auc")@y.values[[1]]
    roc_summary[[gs]] <- data.frame(Signature = gs, AUC = auc, Predicts_For = pos, stringsAsFactors = FALSE)
  }
  roc_df <- dplyr::bind_rows(roc_summary)
  expect_equal(nrow(roc_df), 3)
  expect_true("Composite" %in% roc_df$Signature)
  expect_true(roc_df$AUC[roc_df$Signature == "Up"] > 0.8)

  # Scores without Composite
  scores_no_comp <- scores_with_comp %>% dplyr::select(-Composite_ssGSEA)
  sig_names_no_comp <- c("Up", "Down")
  roc_summary_no <- list()
  for (gs in sig_names_no_comp) {
    col <- paste0(gs, "_ssGSEA")
    v <- scores_no_comp[[col]]
    mp <- mean(v[scores_no_comp$Tissue == "Primary"])
    mm <- mean(v[scores_no_comp$Tissue == "Metastasis"])
    pos <- if (mm >= mp) "Metastasis" else "Primary"
    neg <- if (mm >= mp) "Primary" else "Metastasis"
    pred <- ROCR::prediction(v, scores_no_comp$Tissue, label.ordering = c(neg, pos))
    auc  <- ROCR::performance(pred, "auc")@y.values[[1]]
    roc_summary_no[[gs]] <- data.frame(Signature = gs, AUC = auc, Predicts_For = pos, stringsAsFactors = FALSE)
  }
  roc_df_no <- dplyr::bind_rows(roc_summary_no)
  expect_equal(nrow(roc_df_no), 2)
  expect_false("Composite" %in% roc_df_no$Signature)
})
