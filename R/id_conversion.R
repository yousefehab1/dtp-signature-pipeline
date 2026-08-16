# Gene-ID type detection and conversion, controlled by cfg$id_type ("symbol"
# or "ensembl"). Conversions use org.Hs.eg.db (no internet required); the
# full organism-wide map is cached on disk via cache_rds() so it only builds
# once (and warm-hits the legacy project's cache immediately, since it does
# not depend on the signature panel or on cfg$id_type -- only on which
# `from`/`to` keytypes are requested, which are already encoded in the key).

#' Heuristically detect a gene-ID vector's type
#'
#' @param genes Character vector of gene IDs.
#' @return `"ensembl"` if more than 80% of non-NA entries start with `ENSG`,
#'   else `"symbol"`.
#' @export
detect_id_type <- function(genes) {
  g <- genes[!is.na(genes) & nzchar(genes)]
  if (length(g) == 0) return("symbol")
  if (mean(grepl("^ENSG", g)) > 0.8) "ensembl" else "symbol"
}

#' Strip an Ensembl gene-ID version suffix
#'
#' `ENSG00000123456.3` -> `ENSG00000123456`.
#'
#' @param ids Character vector of Ensembl IDs.
#' @export
strip_ensembl_version <- function(ids) sub("\\.\\d+$", "", ids)

.get_full_id_map <- function(cfg, from, to) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db is required for ID conversion. ",
         "Install with: BiocManager::install('org.Hs.eg.db')")
  }
  key <- paste0("full_id_map_", tolower(from), "_to_", tolower(to))
  # Organism-wide, independent of the signature panel and of cfg$id_type (the
  # direction is already encoded in `key`) -- no vary_on, so this always
  # warm-hits the legacy project's cache instead of rebuilding.
  cache_rds(cfg, key, function() {
    all_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = from)
    map <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                             keys = all_keys, column = to, keytype = from,
                             multiVals = "first")
    )
    map[!is.na(map)]
  })
}

#' Convert gene IDs via a cached organism-wide map
#'
#' @param cfg Config from [dtp_config()].
#' @param genes Character vector of gene IDs to convert.
#' @param from,to AnnotationDbi keytype strings, e.g. `"SYMBOL"`/`"ENSEMBL"`.
#' @export
convert_gene_ids <- function(cfg, genes, from = "SYMBOL", to = "ENSEMBL") {
  map <- .get_full_id_map(cfg, from, to)
  genes_clean <- unique(strip_ensembl_version(genes))
  result <- map[genes_clean]
  n_lost <- sum(is.na(result))
  if (n_lost > 0) {
    message(sprintf("  ID conversion %s->%s: %d/%d genes unmapped (dropped).",
                     from, to, n_lost, length(genes_clean)))
  }
  unique(unname(result[!is.na(result)]))
}

#' Harmonise an expression matrix's row names to `cfg$id_type`
#'
#' Converts row names to the target ID type if they are not already, with an
#' ALIAS-based fallback (symbol-keyed matrices only) for HGNC symbols that
#' have since been renamed (e.g. GSE50760, 2014: `SDPR`->`CAVIN2`), accepting
#' only alias resolutions that are unambiguous.
#'
#' @param cfg Config from [dtp_config()].
#' @param mat A matrix with gene IDs as row names.
#' @param current_type Row-name ID type; auto-detected if `NULL`.
#' @export
harmonise_matrix_ids <- function(cfg, mat, current_type = NULL) {
  if (is.null(current_type)) current_type <- detect_id_type(rownames(mat))
  target <- cfg$id_type
  if (current_type == target) {
    if (target == "ensembl") rownames(mat) <- strip_ensembl_version(rownames(mat))
    return(mat)
  }
  from_col <- if (target == "ensembl") "SYMBOL" else "ENSEMBL"
  to_col   <- if (target == "ensembl") "ENSEMBL" else "SYMBOL"

  orig_ids <- strip_ensembl_version(rownames(mat))
  new_ids <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                           keys = orig_ids, column = to_col, keytype = from_col,
                           multiVals = "first")
  )

  if (from_col == "SYMBOL" && anyNA(new_ids)) {
    miss <- which(is.na(new_ids))
    resid <- orig_ids[miss]
    al <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                             keys = unique(resid), column = to_col,
                             keytype = "ALIAS", multiVals = "list"))
    al1 <- vapply(al, function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 1L) x else NA_character_
    }, character(1))
    filled <- unname(al1[resid])
    new_ids[miss] <- filled
    n_rec <- sum(!is.na(filled))
    if (n_rec > 0) {
      message(sprintf("  harmonise_matrix_ids: +%d rows recovered via ALIAS (renamed symbols).", n_rec))
    }
  }

  keep <- !is.na(new_ids)
  n_lost <- sum(!keep)
  if (n_lost > 0) {
    message(sprintf("  harmonise_matrix_ids %s->%s: %d/%d rows unmapped (dropped).",
                     current_type, target, n_lost, nrow(mat)))
  }
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- unname(new_ids[keep])
  limma::avereps(mat)
}
