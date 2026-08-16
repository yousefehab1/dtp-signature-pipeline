# Molecular subtyping (CMS and PDS) for colorectal cancer cohorts.
# Classifiers require gene symbols, so matrix preparation forces
# id_type = "symbol" regardless of the run configuration.

#' Build a symbol-keyed expression matrix for subtyping
#'
#' CMScaller and PDSclassifier require gene symbols regardless of `cfg$id_type`.
#' This helper forces `id_type = "symbol"` internally and dispatches to the
#' appropriate platform normalisation function.
#'
#' @param cfg Config from [dtp_config()].
#' @param x Input object: a probe x sample matrix for microarray or a
#'   `SummarizedExperiment` for TCGA.
#' @param source Platform source: `"microarray"` or `"tcga"`.
#' @return A log2-scale expression matrix with gene symbols as row names.
#' @export
build_symbol_matrix <- function(cfg, x, source = c("microarray", "tcga")) {
  source <- match.arg(source)
  cfg_symbol <- utils::modifyList(cfg, list(id_type = "symbol"))
  if (source == "microarray") {
    prep_microarray_symbols(cfg_symbol, x)
  } else {
    prep_tcga_tpm(cfg_symbol, x)
  }
}

#' Classify CRC samples into Consensus Molecular Subtypes (CMS)
#'
#' @param emat_symbol_log2 A log2-scale expression matrix with gene symbols as row names.
#' @return A data frame with columns `Sample_ID` and `CMS` (CMS1..CMS4 or NA).
#' @export
call_cms <- function(emat_symbol_log2) {
  message(sprintf("  CMS: classifying %d samples x %d symbols (RNAseq=FALSE).",
                  ncol(emat_symbol_log2), nrow(emat_symbol_log2)))
  res <- CMScaller::CMScaller(emat = emat_symbol_log2, rowNames = "symbol",
                              RNAseq = FALSE, doPlot = FALSE)
  data.frame(Sample_ID = rownames(res),
             CMS = as.character(res$prediction),
             stringsAsFactors = FALSE)
}

#' Classify CRC samples into Pathway-Derived Subtypes (PDS)
#'
#' Runs PDSclassifier if installed; degrades gracefully returning `NULL` if absent.
#'
#' @param emat_symbol_log2 A log2-scale expression matrix with gene symbols as row names.
#' @param species Species name passed to `PDSclassifier::PDSpredict` (default `"human"`).
#' @param threshold Posterior probability threshold for assigning PDS1/2/3 (default 0.6).
#' @return A data frame with columns `Sample_ID` and `PDS`, or `NULL` if PDSclassifier is unavailable.
#' @export
call_pds <- function(emat_symbol_log2, species = "human", threshold = 0.6) {
  if (!requireNamespace("PDSclassifier", quietly = TRUE)) {
    message("  PDS: PDSclassifier not installed; skipping ",
            "(remotes::install_github('sidmall/PDSclassifier@c89a19c')).")
    return(NULL)
  }
  test_df <- data.frame(Gene = rownames(emat_symbol_log2),
                        as.data.frame(emat_symbol_log2, check.names = FALSE),
                        check.names = FALSE, stringsAsFactors = FALSE)
  res <- PDSclassifier::PDSpredict(test_df, species = species, threshold = threshold)
  pds_col <- intersect(c("PDS_call", "prediction", "PDS"), colnames(res))[1]
  id_col  <- intersect(c("Sample_ID", "Sample", "SampleID", "sample"), colnames(res))[1]
  if (is.na(pds_col)) {
    message("  PDS: unexpected PDSpredict output columns (",
            paste(colnames(res), collapse = ", "), "); skipping.")
    return(NULL)
  }
  ids <- if (is.na(id_col)) rownames(res) else res[[id_col]]
  data.frame(Sample_ID = as.character(ids),
             PDS = as.character(res[[pds_col]]),
             stringsAsFactors = FALSE)
}

#' Attach molecular subtype calls to a clinical data frame
#'
#' Maps `Sample_ID` to `clinical_key` via `key_fun` (e.g. barcode->patient collapse)
#' and performs deterministic deduplication before left-joining with clinical annotations.
#'
#' @param clinical A clinical data frame.
#' @param subtype_df Subtype table with `Sample_ID` and a subtype column (from [call_cms()] or [call_pds()]).
#' @param clinical_key Character name of the key column in `clinical`.
#' @param key_fun Function mapping `Sample_ID` values to `clinical_key` values (default: identity).
#' @return Clinical data frame with attached subtype annotations.
#' @export
attach_subtypes <- function(clinical, subtype_df, clinical_key = "Sample_ID", key_fun = identity) {
  if (is.null(subtype_df) || nrow(subtype_df) == 0) return(clinical)
  sd <- subtype_df
  sd[[clinical_key]] <- key_fun(sd$Sample_ID)
  sd <- sd %>%
    dplyr::arrange(Sample_ID) %>%
    dplyr::distinct(.data[[clinical_key]], .keep_all = TRUE)
  if (clinical_key != "Sample_ID") sd$Sample_ID <- NULL
  dplyr::left_join(clinical, sd, by = clinical_key)
}

#' Split a clinical data frame into per-subtype cohorts
#'
#' @param clinical Clinical data frame with a subtype column.
#' @param subtype_col Name of the subtype column (e.g. `"CMS"` or `"PDS"`).
#' @param prefix Cohort name prefix (e.g. `"GSE"` or `"TCGA"`).
#' @return Named list of data frames, one per non-missing subtype level, named `"<prefix>_<level>"`.
#' @export
subtype_cohorts <- function(clinical, subtype_col, prefix) {
  if (!subtype_col %in% colnames(clinical)) return(list())
  lvls <- sort(unique(stats::na.omit(clinical[[subtype_col]])))
  lvls <- lvls[lvls != ""]
  out <- lapply(lvls, function(lv)
    clinical[!is.na(clinical[[subtype_col]]) & clinical[[subtype_col]] == lv, , drop = FALSE])
  names(out) <- paste0(prefix, "_", lvls)
  out
}
