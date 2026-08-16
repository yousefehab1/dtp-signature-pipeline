# Metastasis differential-expression and GSEA analysis module (GSE50760 RNA-seq).
# Decomposed into testable internal helpers (.mets_*) and orchestrated by run_mets_de().
#
# Analytical core preserved verbatim from thesis:
#   - GSE50760 18 patients x 3 tissues (Normal, Primary, Metastasis)
#   - FPKM duplicate symbol summing (additive Cufflinks sub-features)
#   - >=10% sample detection expression filter before GSEA ranking
#   - 3-way paired limma model (~ Patient + Tissue_for_DE, ref = "Primary")
#   - Liver-purity negative-control validation and covariate adjustment (~ Patient + Tissue_for_DE + LiverScore)
#   - clusterProfiler GSEA on canonical signature panel
#   - 3-way PCA (global + core signature)
#   - Single-sample paired ROC AUC and Wilcoxon diagnostics on liver-corrected matrix

.mets_load_metadata <- function(cfg) {
  meta <- cache_rds(cfg, "GSE50760_meta", function() GEOquery::getGEO("GSE50760", GSEMatrix = TRUE))
  pdata <- Biobase::pData(meta[[1]])
  titles <- as.character(pdata$title)
  gsm_ids <- rownames(pdata)

  tissue <- dplyr::case_when(
    grepl("normal",  titles, ignore.case = TRUE) ~ "Normal",
    grepl("primary", titles, ignore.case = TRUE) ~ "Primary",
    grepl("metasta", titles, ignore.case = TRUE) ~ "Metastasis",
    TRUE                                         ~ NA_character_
  )
  patient_code <- sub(".*?(AMC_[0-9]+).*", "\\1", titles)
  patient_code[!grepl("^AMC_[0-9]+$", patient_code)] <- NA_character_
  sample_meta <- data.frame(
    GSM     = gsm_ids,
    Tissue  = tissue,
    Patient = patient_code,
    stringsAsFactors = FALSE
  )
  if (anyNA(sample_meta$Tissue) || anyNA(sample_meta$Patient)) {
    stop("Could not parse Tissue/Patient from titles.")
  }
  if (!all(table(sample_meta$Patient, sample_meta$Tissue) == 1)) {
    stop("Not a clean 18x3 design.")
  }
  message(sprintf("Parsed %d samples / %d patients / %d tissues.",
                  nrow(sample_meta), dplyr::n_distinct(sample_meta$Patient), dplyr::n_distinct(sample_meta$Tissue)))
  sample_meta
}

.mets_load_fpkm <- function(cfg, sample_meta) {
  supp <- cache_rds(cfg, "GSE50760_supp_dir", function() {
    sf  <- GEOquery::getGEOSuppFiles("GSE50760", baseDir = cfg$cache_dir)
    ex  <- file.path(cfg$cache_dir, "GSE50760_extracted")
    dir.create(ex, showWarnings = FALSE, recursive = TRUE)
    utils::untar(rownames(sf)[1], exdir = ex)
    ex
  })
  supp_dir <- supp
  if (!dir.exists(supp_dir)) {
    if (!is.null(cfg$legacy_cache_dir) && dir.exists(file.path(cfg$legacy_cache_dir, basename(supp_dir)))) {
      supp_dir <- file.path(cfg$legacy_cache_dir, basename(supp_dir))
    } else if (dir.exists(file.path(cfg$cache_dir, basename(supp_dir)))) {
      supp_dir <- file.path(cfg$cache_dir, basename(supp_dir))
    }
  }

  all_files <- list.files(supp_dir, full.names = TRUE)
  files_to_read <- character()
  for (id in sample_meta$GSM) {
    m <- all_files[grepl(paste0(id, "([._-]|$)"), all_files)]
    if (length(m)) files_to_read[id] <- m[1]
  }
  if (length(setdiff(sample_meta$GSM, names(files_to_read))) > 0) {
    stop("Missing FPKM file(s): ", paste(setdiff(sample_meta$GSM, names(files_to_read)), collapse = ", "))
  }

  count_list <- list()
  for (id in names(files_to_read)) {
    td <- utils::read.delim(files_to_read[id], header = TRUE, stringsAsFactors = FALSE)[, c(1, 2)]
    colnames(td) <- c("gene_id", id)
    td$gene_id <- as.character(td$gene_id)
    td[[id]] <- suppressWarnings(as.numeric(td[[id]]))
    # Sum duplicate symbols: additive Cufflinks sub-features, not redundant measurements.
    count_list[[id]] <- td %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(!!id := sum(.data[[id]], na.rm = TRUE), .groups = "drop") %>%
      as.data.frame()
  }
  fpkm_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
  fpkm <- as.matrix(fpkm_df[, setdiff(colnames(fpkm_df), "gene_id")])
  rownames(fpkm) <- fpkm_df$gene_id
  fpkm <- fpkm[, sample_meta$GSM]
  if (anyNA(fpkm)) fpkm <- fpkm[rowSums(is.na(fpkm)) == 0, ]
  log2_matrix <- log2(fpkm + 1)

  # Harmonise row IDs to cfg$id_type (Ensembl or Symbol)
  log2_matrix <- harmonise_matrix_ids(cfg, log2_matrix)

  # Drop undetected / invariant genes before any modelling.
  # An all-zero (or near-constant) FPKM row yields limma t = 0, and in the GSEA
  # ranked list all such genes are ordered purely by random tie-breaking —
  # a large contiguous block of noise-ranked genes that dilutes enrichment-score
  # normalisation. Require detectable expression (FPKM > 0, log2 > 0) in >= 10% of samples.
  min_det <- max(3L, floor(0.10 * ncol(log2_matrix)))
  detected <- rowSums(log2_matrix > 0) >= min_det
  message(sprintf("  Expression filter: %d/%d genes retained (detected in >= %d samples).",
                  sum(detected), length(detected), min_det))
  log2_matrix[detected, , drop = FALSE]
}

.mets_paired_wilcox <- function(df, value_col, group_col, id_col, a, b) {
  w <- df %>%
    dplyr::filter(.data[[group_col]] %in% c(a, b)) %>%
    dplyr::select(dplyr::all_of(c(id_col, group_col, value_col))) %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(group_col), values_from = dplyr::all_of(value_col))
  w <- w[stats::complete.cases(w[, c(a, b)]), ]
  stats::wilcox.test(w[[a]], w[[b]], paired = TRUE)$p.value
}

.mets_compute_pca <- function(mat, cd) {
  mv <- mat[apply(mat, 1, stats::var) > 0, , drop = FALSE]
  if (nrow(mv) < 2) return(NULL)
  pca <- stats::prcomp(t(mv), scale. = TRUE)
  pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
  tibble::tibble(
    GSM         = cd$GSM,
    Tissue      = cd$Tissue,
    Patient     = cd$Patient,
    PC1         = pca$x[, 1],
    PC2         = pca$x[, 2],
    PC1_var_pct = pv[1],
    PC2_var_pct = pv[2]
  )
}

.mets_run_gsea <- function(de_res, custom_t2g) {
  rnk <- de_res$stat
  names(rnk) <- de_res$gene
  # Tie-breaker jitter: preserves legacy ranking behavior on tied statistics
  rnk <- sort(rnk + stats::runif(length(rnk), 1e-9, 1e-8), decreasing = TRUE)
  clusterProfiler::GSEA(
    geneList      = rnk,
    TERM2GENE     = custom_t2g,
    pvalueCutoff  = 1,
    eps           = 0,
    minGSSize     = 5,
    maxGSSize     = 500,
    BPPARAM       = BiocParallel::SerialParam()
  )
}

#' Run metastasis differential-expression and GSEA analysis
#'
#' Evaluates paired Normal/Primary/Metastasis samples in GSE50760 (18 patients x 3 tissues)
#' using limma paired models, liver-purity covariate adjustment, GSEA enrichment,
#' PCA decomposition, and single-sample paired ROC/Wilcoxon diagnostics.
#'
#' @param cfg Config from [dtp_config()].
#' @param panel_tbl Canonical signature panel table from [load_signature_panel()].
#' @param out_root Output directory path (reserved for figure generation in future phases).
#' @param core_signature Name of signature to use for core PCA; defaults to `"Up"`.
#' @return A named list containing:
#'   \item{sample_meta}{GSE50760 sample metadata (GSM, Tissue, Patient).}
#'   \item{res}{Primary vs Metastasis unadjusted limma topTable results.}
#'   \item{comparison_de}{DE comparison table between unadjusted and purity-adjusted models.}
#'   \item{liver_df}{Per-sample liver ssGSEA validation scores with metadata.}
#'   \item{liver_pvalues}{Paired Wilcoxon p-values across tissue pairs for liver score.}
#'   \item{pvl}{Alias for `liver_pvalues`.}
#'   \item{gsea_pm}{GSEA summary table for Primary vs Metastasis (unadjusted).}
#'   \item{gsea_pm_adj}{GSEA summary table for Primary vs Metastasis (purity-adjusted).}
#'   \item{gsea_pn}{GSEA summary table for Primary vs Normal.}
#'   \item{gsea_pm_obj}{Raw `gseaResult` object for Primary vs Metastasis (unadjusted).}
#'   \item{gsea_pm_adj_obj}{Raw `gseaResult` object for Primary vs Metastasis (purity-adjusted).}
#'   \item{gsea_pn_obj}{Raw `gseaResult` object for Primary vs Normal.}
#'   \item{pvm_cmp}{Adjusted vs unadjusted GSEA NES/FDR comparison table.}
#'   \item{pca_global}{PCA coordinates and % variance explained for global expression matrix.}
#'   \item{pca_core}{PCA coordinates and % variance explained for core signature genes.}
#'   \item{scores_df}{Paired single-sample ssGSEA scores on liver-corrected matrix.}
#'   \item{roc_summary}{Single-sample ROC AUC and prediction direction per signature.}
#'   \item{ssgsea_pvalues}{Paired Wilcoxon raw and BH-adjusted p-values for ssGSEA scores.}
#' @export
run_mets_de <- function(cfg, panel_tbl, out_root = "output", core_signature = "Up") {
  # out_root is accepted for interface compatibility with future plotting phases
  force(out_root)
  set.seed(cfg$global_seed)

  panel_list     <- panel_to_list(panel_tbl, id_type = cfg$id_type)
  composite_defs <- load_composite_defs(cfg)

  # ---- Part 1: Download & Parse Metadata ----
  sample_meta <- .mets_load_metadata(cfg)

  # ---- Part 2: Build FPKM Matrix ----
  log2_matrix <- .mets_load_fpkm(cfg, sample_meta)

  # ---- Part 3: Paired DE (limma ~ Patient + Tissue) ----
  cd <- sample_meta
  rownames(cd) <- cd$GSM
  cd <- cd[colnames(log2_matrix), ]
  cd$Tissue        <- factor(cd$Tissue, levels = c("Normal", "Primary", "Metastasis"))
  cd$Patient       <- factor(cd$Patient)
  cd$Tissue_for_DE <- stats::relevel(cd$Tissue, ref = "Primary")

  design  <- stats::model.matrix(~ Patient + Tissue_for_DE, data = cd)
  fit     <- limma::eBayes(limma::lmFit(log2_matrix, design))
  coef_nm <- "Tissue_for_DEMetastasis"

  res <- limma::topTable(fit, coef = coef_nm, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(stat = t, padj = adj.P.Val, pvalue = P.Value) %>%
    dplyr::filter(!is.na(stat))
  rownames(res) <- res$gene

  keep2 <- cd$Tissue %in% c("Primary", "Metastasis")
  log2_2way <- log2_matrix[, keep2]
  cd2 <- droplevels(cd[keep2, ])

  # ---- Part 4: Liver-purity Covariate & Negative-Control Validation ----
  liver_genes_sym <- cache_rds(cfg, "msigdbr_hsiao_liver_specific_genes", function() {
    msigdbr::msigdbr(species = "Homo sapiens") %>%
      dplyr::filter(gs_name == "HSIAO_LIVER_SPECIFIC_GENES") %>%
      dplyr::pull(gene_symbol) %>%
      unique()
  })
  liver_genes <- if (cfg$id_type == "ensembl") {
    convert_gene_ids(cfg, liver_genes_sym, from = "SYMBOL", to = "ENSEMBL")
  } else {
    liver_genes_sym
  }
  liver_present <- intersect(liver_genes, rownames(log2_matrix))
  liver_gs <- GSEABase::GeneSetCollection(GSEABase::GeneSet(liver_present, setName = "Liver"))
  liver_scores <- GSVA::gsva(GSVA::ssgseaParam(exprData = log2_matrix, geneSets = liver_gs))
  liver_df <- as.data.frame(t(liver_scores)) %>%
    tibble::rownames_to_column("GSM") %>%
    dplyr::left_join(cd %>% dplyr::select(GSM, Tissue, Patient), by = "GSM") %>%
    dplyr::arrange(Patient, Tissue)

  LiverScore <- as.numeric(liver_scores["Liver", colnames(log2_2way)])
  cd2$LiverScore <- LiverScore

  comp <- list(c("Normal", "Primary"), c("Primary", "Metastasis"), c("Normal", "Metastasis"))
  pvl  <- vapply(comp, function(x) {
    .mets_paired_wilcox(liver_df, "Liver", "Tissue", "Patient", x[1], x[2])
  }, numeric(1))
  names(pvl) <- c("Normal_vs_Primary", "Primary_vs_Metastasis", "Normal_vs_Metastasis")
  message(sprintf("Liver validation — N/P p=%.3g | P/M p=%.3g | N/M p=%.3g", pvl[1], pvl[2], pvl[3]))

  # ---- Part 5: Purity Sensitivity Model + Corrected Matrix ----
  design_adj <- stats::model.matrix(~ Patient + Tissue_for_DE + LiverScore, data = cd2)
  res_adj    <- limma::topTable(limma::eBayes(limma::lmFit(log2_2way, design_adj)), coef = coef_nm,
                                number = Inf, sort.by = "none")

  comparison_de <- res %>%
    dplyr::select(gene, logFC, stat, pvalue, padj) %>%
    dplyr::rename(logFC_unadj = logFC, t_unadj = stat, pval_unadj = pvalue, padj_unadj = padj) %>%
    dplyr::left_join(
      res_adj %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::select(gene, logFC, t, P.Value, adj.P.Val) %>%
        dplyr::rename(logFC_adj = logFC, t_adj = t, pval_adj = P.Value, padj_adj = adj.P.Val),
      by = "gene"
    ) %>%
    dplyr::mutate(
      Is_Liver_Gene = gene %in% liver_genes,
      Status = dplyr::case_when(
        padj_unadj < 0.05 & padj_adj >= 0.05 ~ "Lost significance after adjustment",
        padj_unadj >= 0.05 & padj_adj < 0.05 ~ "Gained significance after adjustment",
        padj_unadj < 0.05 & padj_adj < 0.05  ~ "Robust to adjustment",
        TRUE                                  ~ "Not significant either way"
      )
    ) %>%
    dplyr::arrange(padj_adj)

  preserve <- stats::model.matrix(~ Patient + Tissue, data = cd2)
  log2_2way_corrected <- limma::removeBatchEffect(log2_2way, covariates = LiverScore, design = preserve)

  # ---- Part 6: GSEA on the Canonical Signature Panel ----
  custom_t2g <- dplyr::bind_rows(lapply(
    names(panel_list),
    function(nm) data.frame(Term = nm, Gene = panel_list[[nm]], stringsAsFactors = FALSE)
  ))

  # Comparison 1: Primary vs Metastasis (liver-unadjusted)
  message("  Running GSEA: Primary vs Metastasis (liver-unadjusted)")
  gsea_pm <- .mets_run_gsea(res, custom_t2g)

  # Comparison 1b: Primary vs Metastasis (purity-adjusted, PRIMARY figure)
  message("  Running GSEA: Primary vs Metastasis (purity-adjusted for liver score)")
  res_adj_ranked <- res_adj %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(stat = t, padj = adj.P.Val, pvalue = P.Value) %>%
    dplyr::filter(!is.na(stat))
  gsea_pm_adj <- .mets_run_gsea(res_adj_ranked, custom_t2g)

  # Comparison 2: Primary vs Normal
  # Extracted from the 3-way fit; coefficient Tissue_for_DENormal gives Normal - Primary,
  # so negate stat/logFC so positive = up in Primary.
  message("  Running GSEA: Primary vs Normal")
  res_pn <- limma::topTable(fit, coef = "Tissue_for_DENormal", number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(stat = t, padj = adj.P.Val, pvalue = P.Value) %>%
    dplyr::filter(!is.na(stat)) %>%
    dplyr::mutate(stat = -stat, logFC = -logFC)
  gsea_pn <- .mets_run_gsea(res_pn, custom_t2g)

  # Comparison table: P/M unadjusted vs adjusted
  .pm_u <- as.data.frame(gsea_pm)     %>% dplyr::select(ID, NES, p.adjust)
  .pm_a <- as.data.frame(gsea_pm_adj) %>% dplyr::select(ID, NES, p.adjust)
  pvm_cmp <- dplyr::full_join(
    dplyr::rename(.pm_u, NES_unadjusted = NES, FDR_unadjusted = p.adjust),
    dplyr::rename(.pm_a, NES_adjusted   = NES, FDR_adjusted   = p.adjust),
    by = "ID"
  ) %>%
    dplyr::mutate(
      dNES           = NES_adjusted - NES_unadjusted,
      Sig_unadjusted = FDR_unadjusted < 0.05,
      Sig_adjusted   = FDR_adjusted   < 0.05
    ) %>%
    dplyr::arrange(FDR_adjusted)

  # ---- Part 7: PCA (3-way) ----
  pca_global <- .mets_compute_pca(log2_matrix, cd)
  pca_core <- NULL
  if (core_signature %in% names(panel_list)) {
    cg <- intersect(panel_list[[core_signature]], rownames(log2_matrix))
    if (length(cg) >= 2) {
      pca_core <- .mets_compute_pca(log2_matrix[cg, , drop = FALSE], cd)
    }
  }

  # ---- Part 8: Single-Sample Diagnostics (ROC + Paired Wilcoxon) ----
  # Scores all panel signatures together in one ssGSEA call (no Up/Down isolation)
  gs_collection <- build_gene_sets(panel_list, which = NULL)
  scores <- run_ssgsea(log2_2way_corrected, gs_collection, suffix = cfg$score_suffix)
  scores <- add_all_composites(scores, composite_defs, suffix = cfg$score_suffix)

  scores_df <- scores %>%
    tibble::rownames_to_column("GSM") %>%
    dplyr::left_join(cd2 %>% dplyr::select(GSM, Tissue, Patient), by = "GSM") %>%
    dplyr::arrange(Patient, Tissue)
  scores_df$Tissue <- factor(scores_df$Tissue, levels = c("Primary", "Metastasis"))

  composite_names <- if (!is.null(composite_defs) && nrow(composite_defs) > 0) {
    intersect(composite_defs$Composite_ID, sub(paste0(cfg$score_suffix, "$"), "", colnames(scores_df)))
  } else {
    character(0)
  }
  sig_names <- unique(c(names(panel_list), composite_names))

  roc_summary <- list()
  for (gs in sig_names) {
    col <- paste0(gs, cfg$score_suffix)
    if (!col %in% colnames(scores_df)) next
    v <- scores_df[[col]]
    mp <- mean(v[scores_df$Tissue == "Primary"])
    mm <- mean(v[scores_df$Tissue == "Metastasis"])
    if (mm >= mp) {
      pos <- "Metastasis"
      neg <- "Primary"
    } else {
      pos <- "Primary"
      neg <- "Metastasis"
    }
    pred <- ROCR::prediction(v, scores_df$Tissue, label.ordering = c(neg, pos))
    auc  <- ROCR::performance(pred, "auc")@y.values[[1]]
    roc_summary[[gs]] <- data.frame(
      Signature    = gs,
      AUC          = auc,
      Predicts_For = pos,
      stringsAsFactors = FALSE
    )
  }
  roc_df <- dplyr::bind_rows(roc_summary) %>% dplyr::arrange(dplyr::desc(AUC))

  raw_p <- vapply(sig_names, function(gs) {
    col <- paste0(gs, cfg$score_suffix)
    if (!col %in% colnames(scores_df)) return(NA_real_)
    .mets_paired_wilcox(scores_df, col, "Tissue", "Patient", "Primary", "Metastasis")
  }, numeric(1))
  adj_p <- stats::p.adjust(raw_p, method = "BH")
  ssgsea_pvalues <- data.frame(
    Signature = sig_names,
    p_raw     = raw_p,
    p_adj_BH  = adj_p,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(p_adj_BH)

  message("Mets DE module complete.")
  invisible(list(
    sample_meta     = sample_meta,
    res             = res,
    comparison_de   = comparison_de,
    liver_df        = liver_df,
    liver_pvalues   = pvl,
    pvl             = pvl,
    gsea_pm         = as.data.frame(gsea_pm),
    gsea_pm_adj     = as.data.frame(gsea_pm_adj),
    gsea_pn         = as.data.frame(gsea_pn),
    gsea_pm_obj     = gsea_pm,
    gsea_pm_adj_obj = gsea_pm_adj,
    gsea_pn_obj     = gsea_pn,
    pvm_cmp         = pvm_cmp,
    pca_global      = pca_global,
    pca_core        = pca_core,
    scores_df       = scores_df,
    roc_summary     = roc_df,
    ssgsea_pvalues  = ssgsea_pvalues
  ))
}
