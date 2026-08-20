# Run scaffolding, disk caching, reproducibility.

#' Start a run
#'
#' Sets the global seed and creates a timestamped output directory under
#' `out_dir` so a run never silently overwrites a previous one.
#'
#' @param cfg Config from [dtp_config()].
#' @param prefix Output directory name prefix.
#' @param out_dir Parent directory for run output.
#' @export
init_run <- function(cfg, prefix = "run", out_dir = "output") {
  set.seed(cfg$global_seed)
  dir.create(cfg$cache_dir, showWarnings = FALSE, recursive = TRUE)
  root <- file.path(out_dir, paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M")))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  message("Run output root: ", root)
  root
}

#' Finish a run
#'
#' Writes `sessionInfo()` to the run directory for exact package-version
#' reproducibility.
#'
#' @param root Run output directory, as returned by [init_run()].
#' @export
finalize_run <- function(root) {
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(root, "sessionInfo.txt"))
  message("sessionInfo saved to ", file.path(root, "sessionInfo.txt"))
}

#' Disk-backed cache
#'
#' `fun` only runs on a cache miss. Entries whose content depends on config
#' (e.g. `id_type` or the signature panel) must pass `vary_on` -- this hashes
#' those values into the filename so a config change can never silently load
#' a stale, incompatible cache entry.
#'
#' @param cfg Config from [dtp_config()].
#' @param key Cache key; sanitized to a safe filename.
#' @param fun Zero-argument function producing the value to cache.
#' @param vary_on Optional list of config-derived values the cached value
#'   depends on (e.g. `list(id_type = cfg$id_type)`). Included in the cache
#'   filename via a short hash. Omit for raw downloads that are independent
#'   of `cfg`.
#' @export
cache_rds <- function(cfg, key, fun, vary_on = NULL) {
  safe_key <- gsub("[^A-Za-z0-9_]+", "_", key)
  fname <- if (is.null(vary_on)) {
    paste0(safe_key, ".rds")
  } else {
    paste0(safe_key, "__", substr(digest::digest(vary_on), 1, 10), ".rds")
  }

  f_new <- file.path(cfg$cache_dir, fname)
  if (file.exists(f_new)) {
    message("[cache] ", key, " (load)")
    return(readRDS(f_new))
  }

  message("[cache] ", key, " (build)")
  v <- fun()
  dir.create(cfg$cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(v, f_new)
  v
}

#' Make a data frame safe to write as CSV
#'
#' Flattens list columns to pipe-separated strings and strips embedded
#' newlines from character columns.
#'
#' @param df A data frame.
#' @export
csv_safe <- function(df) {
  df %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.list),
                                 ~ sapply(., \(x) paste(unlist(x), collapse = " | ")))) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ gsub("\n|\r", " ", .)))
}
