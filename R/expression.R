# Per-platform normalisation, converging on a log2-scale, gene-indexed matrix
# ready for ssGSEA. ssGSEA is a within-sample rank statistic, so log2 vs.
# linear scale is rank-invariant -- log2 is used everywhere purely so every
# platform shares one convention.

#' Prepare an Affymetrix microarray matrix for scoring
#'
#' Maps probes to the target gene-ID type (`cfg$id_type`), then collapses
#' duplicate genes by averaging in log space.
#'
#' @param cfg Config from [dtp_config()].
#' @param exp_matrix Probe x sample expression matrix.
#' @param anno_db An AnnotationDbi platform annotation object; defaults to
#'   `hgu133plus2.db::hgu133plus2.db`.
#' @export
prep_microarray_symbols <- function(cfg, exp_matrix, anno_db = NULL) {
  if (is.null(anno_db)) anno_db <- hgu133plus2.db::hgu133plus2.db
  target_col <- if (cfg$id_type == "ensembl") "ENSEMBL" else "SYMBOL"
  id_label   <- cfg$id_type

  # multiVals="first" is deterministic within a fixed AnnotationDbi/platform-db
  # version (pinned via renv.lock) but is an accepted, order-dependent
  # simplification for probes mapping to more than one gene -- changing the
  # selection strategy would change scientific results, which is out of scope
  # for this rewrite.
  probe_ids <- AnnotationDbi::mapIds(anno_db, keys = rownames(exp_matrix),
                                      column = target_col, keytype = "PROBEID",
                                      multiVals = "first")
  probe_all <- AnnotationDbi::mapIds(anno_db, keys = rownames(exp_matrix),
                                      column = target_col, keytype = "PROBEID",
                                      multiVals = "CharacterList")
  n_multi <- sum(lengths(probe_all) > 1)
  message(sprintf("  probe annotation: %d/%d (%.1f%%) multi-map; multiVals='first'.",
                   n_multi, length(probe_all), 100 * n_multi / length(probe_all)))

  mat <- exp_matrix
  rownames(mat) <- unname(probe_ids)
  n_unmapped <- sum(is.na(rownames(mat)))
  message(sprintf("  probe->%s: %d/%d (%.1f%%) unmapped, dropped.",
                   id_label, n_unmapped, nrow(mat), 100 * n_unmapped / nrow(mat)))
  mat <- mat[!is.na(rownames(mat)), , drop = FALSE]
  if (cfg$id_type == "ensembl") rownames(mat) <- strip_ensembl_version(rownames(mat))
  message(sprintf("  collapsing duplicates: %d probes -> %d %ss.",
                   nrow(mat), length(unique(rownames(mat))), id_label))
  limma::avereps(mat)
}

#' Prepare a TCGA STAR-Counts SummarizedExperiment for scoring
#'
#' log2(TPM+1)-transforms, collapses duplicate genes in log space, and
#' filters to primary-tumour samples. Used identically for CRC and
#' pan-cancer cohorts.
#'
#' @param cfg Config from [dtp_config()].
#' @param se A `SummarizedExperiment` with a `tpm_unstrand` assay.
#' @export
prep_tcga_tpm <- function(cfg, se) {
  if (!"tpm_unstrand" %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay 'tpm_unstrand' not found. Available: ",
         paste(SummarizedExperiment::assayNames(se), collapse = ", "))
  }

  rd  <- SummarizedExperiment::rowData(se)
  mat <- SummarizedExperiment::assay(se, "tpm_unstrand")

  if (cfg$id_type == "ensembl") {
    if (!"gene_id" %in% colnames(rd)) {
      stop("rowData has no 'gene_id' column; cannot use id_type='ensembl' for TCGA.")
    }
    rownames(mat) <- strip_ensembl_version(rd$gene_id)
  } else {
    rownames(mat) <- rd$gene_name
  }

  mat <- mat[!is.na(rownames(mat)), , drop = FALSE]
  mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]

  mat <- log2(mat + 1)
  mat <- limma::avereps(mat)

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  keep <- cd$sample_type %in% cfg$tumour_codes
  if (sum(keep) == 0) stop("No samples of type(s): ", paste(cfg$tumour_codes, collapse = ", "))
  mat[, keep, drop = FALSE]
}

#' Percentage of a gene set present in an expression matrix
#'
#' @param mat An expression matrix with gene IDs as row names.
#' @param genes Character vector of gene IDs.
#' @export
gene_set_coverage <- function(mat, genes) mean(genes %in% rownames(mat)) * 100
