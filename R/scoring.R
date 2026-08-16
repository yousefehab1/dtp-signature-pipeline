# ssGSEA scoring and the (now data-driven) Composite score.

#' Score a log2-scale expression matrix with ssGSEA
#'
#' @param mat A log2-scale, gene-ID-indexed expression matrix (genes x samples).
#' @param gene_sets A `GSEABase::GeneSetCollection` from [build_gene_sets()].
#' @param suffix Column-name suffix; result columns are `"<Signature>"` + `suffix`.
#' @return A data frame (samples x signatures).
#' @export
run_ssgsea <- function(mat, gene_sets, suffix = "_ssGSEA") {
  param <- GSVA::ssgseaParam(exprData = as.matrix(mat), geneSets = gene_sets)
  res <- GSVA::gsva(param)
  df <- as.data.frame(t(res))
  colnames(df) <- paste0(colnames(df), suffix)
  df
}

#' Add one composite score (positive signature minus negative signature)
#'
#' Composite is computed post-hoc from already-scored columns -- an
#' arithmetic operation on scores, not a gene set of its own -- so which
#' signatures form a composite is data ([load_composite_defs()]), not a
#' hardcoded default.
#'
#' @param scores_df Data frame of `<Signature>_ssGSEA` columns, as returned
#'   by [run_ssgsea()].
#' @param positive,negative Signature_IDs to subtract.
#' @param name Name of the resulting composite (unsuffixed).
#' @param suffix Score column suffix, matching the one used in `scores_df`.
#' @export
add_composite <- function(scores_df, positive, negative, name = "Composite", suffix = "_ssGSEA") {
  pc <- paste0(positive, suffix)
  nc <- paste0(negative, suffix)
  if (all(c(pc, nc) %in% colnames(scores_df))) {
    scores_df[[paste0(name, suffix)]] <- scores_df[[pc]] - scores_df[[nc]]
  }
  scores_df
}

#' Add every composite defined in a composite-definitions table
#'
#' @param scores_df Data frame of `<Signature>_ssGSEA` columns.
#' @param composite_defs Table from [load_composite_defs()].
#' @param suffix Score column suffix.
#' @export
add_all_composites <- function(scores_df, composite_defs, suffix = "_ssGSEA") {
  for (i in seq_len(nrow(composite_defs))) {
    scores_df <- add_composite(scores_df,
                                positive = composite_defs$Positive_Signature[i],
                                negative = composite_defs$Negative_Signature[i],
                                name = composite_defs$Composite_ID[i],
                                suffix = suffix)
  }
  scores_df
}

#' Extract the 12-character TCGA patient ID from a full barcode
#'
#' @param barcodes Character vector of TCGA barcodes/aliquot IDs.
#' @export
tcga_patient_id <- function(barcodes) {
  stringr::str_extract(gsub("\\.", "-", toupper(barcodes)), "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
}
