# Signature panel loading. The panel is a long-format table (one row per gene
# per signature, inst/signatures/panel.csv) rather than the old project's wide,
# blank-padded CSV -- adding or removing a signature is editing that one file,
# no code change or column-list update required.

#' Load the canonical signature panel
#'
#' @param cfg Config from [dtp_config()].
#' @param path Panel CSV path; defaults to `cfg$signature_panel_path`.
#' @return A tibble with columns `Signature_ID`, `Signature_Name`, `Gene_ID`,
#'   `Gene_Symbol`, `Direction`, `Source_Citation`, `Source_Gene_Count`, `Notes`.
#' @export
load_signature_panel <- function(cfg, path = cfg$signature_panel_path) {
  stopifnot(file.exists(path))
  panel <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)

  required_cols <- c("Signature_ID", "Gene_ID")
  missing_cols <- setdiff(required_cols, names(panel))
  if (length(missing_cols)) {
    stop("Signature panel is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  if (any(is.na(panel$Signature_ID) | trimws(panel$Signature_ID) == "")) {
    stop("Signature panel has row(s) with a blank Signature_ID.")
  }
  if (anyDuplicated(panel[c("Signature_ID", "Gene_ID")])) {
    stop("Signature panel has duplicate (Signature_ID, Gene_ID) rows.")
  }

  n_genes <- table(panel$Signature_ID)
  small <- names(n_genes)[n_genes < 5]
  if (length(small)) {
    warning("Signature(s) with fewer than 5 genes: ", paste(small, collapse = ", "), call. = FALSE)
  }

  message(sprintf("Signature panel: %d sets loaded (%s).",
                   length(n_genes), paste(names(n_genes), collapse = ", ")))
  panel
}

#' Convert the long-format panel into a named list of gene vectors
#'
#' @param panel_tbl Panel table from [load_signature_panel()].
#' @param id_type "ensembl" or "symbol" -- selects `Gene_ID` (always Ensembl,
#'   the panel's canonical storage namespace) or `Gene_Symbol`.
#' @export
panel_to_list <- function(panel_tbl, id_type = c("ensembl", "symbol")) {
  id_type <- match.arg(id_type)
  gene_col <- if (id_type == "ensembl") "Gene_ID" else "Gene_Symbol"
  split(panel_tbl[[gene_col]], panel_tbl$Signature_ID) |>
    lapply(function(g) unique(g[!is.na(g) & trimws(g) != ""]))
}

#' Build a GSEABase::GeneSetCollection from selected signatures
#'
#' @param panel_list Named list of gene vectors from [panel_to_list()].
#' @param which Character vector of `Signature_ID`s to include; `NULL` = all.
#' @export
build_gene_sets <- function(panel_list, which = NULL) {
  if (!is.null(which)) {
    miss <- setdiff(which, names(panel_list))
    if (length(miss)) stop("Signatures not found in panel: ", paste(miss, collapse = ", "))
    panel_list <- panel_list[which]
  }
  GSEABase::GeneSetCollection(
    lapply(names(panel_list), function(nm) GSEABase::GeneSet(panel_list[[nm]], setName = nm))
  )
}

#' Warn if any gene appears in both signatures of a composite
#'
#' A gene in both the positive and negative signature of a composite
#' contributes +1 and -1 simultaneously, which is very likely unintended.
#'
#' @param panel_list Named list of gene vectors from [panel_to_list()].
#' @param positive,negative Signature_IDs to compare.
#' @export
check_updown_overlap <- function(panel_list, positive = "Up", negative = "Down") {
  if (all(c(positive, negative) %in% names(panel_list))) {
    ov <- intersect(panel_list[[positive]], panel_list[[negative]])
    if (length(ov)) {
      warning(length(ov), " gene(s) in both ", positive, " and ", negative, ": ",
              paste(utils::head(ov, 10), collapse = ", "), call. = FALSE)
    }
  }
}

#' Load the composite-signature definitions
#'
#' @param cfg Config from [dtp_config()].
#' @param path Composite-definitions CSV path; defaults to `cfg$composite_defs_path`.
#' @return A tibble with columns `Composite_ID`, `Display_Name`,
#'   `Positive_Signature`, `Negative_Signature`, `Is_Default`.
#' @export
load_composite_defs <- function(cfg, path = cfg$composite_defs_path) {
  stopifnot(file.exists(path))
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

#' Hash the signature panel's content
#'
#' Used as a cache-key component so a change to the panel (add/remove a
#' signature or gene) never silently reuses a cached score/ID-mapping built
#' from a different panel.
#'
#' @param panel_tbl Panel table from [load_signature_panel()].
#' @export
panel_hash <- function(panel_tbl) {
  ordered <- panel_tbl[order(panel_tbl$Signature_ID, panel_tbl$Gene_ID),
                        c("Signature_ID", "Gene_ID")]
  substr(digest::digest(ordered), 1, 10)
}
