# Run configuration. Every analysis/plotting function in this package takes a
# `cfg` list (built here) instead of reading bare global constants -- this is
# what lets cache keys depend on config (id_type, panel version) and lets the
# FDR grouping be a named, non-overridable enum instead of an ad-hoc `by=`
# argument scattered across call sites (see apply_fdr() in R/stats.R).

#' Build the run configuration
#'
#' Reads `pipeline/pipeline_config.yml`, merges it with package defaults, and
#' resolves relative paths against `root`. Pass `overrides` to change values
#' for a single run without editing the YAML (e.g. in tests).
#'
#' @param config_path Path to the YAML config file.
#' @param root Directory relative paths in the YAML are resolved against.
#' @param overrides Named list of values to override after loading the YAML.
#' @export
dtp_config <- function(config_path = "pipeline/pipeline_config.yml",
                        root = ".",
                        overrides = list()) {
  defaults <- list(
    id_type = "ensembl",
    score_suffix = "_ssGSEA",
    os_cutpoint = 36,
    rfs_cutpoint = 36,
    min_group_n = 10,
    min_km_group_n = 5,
    min_events = 5,
    min_cox_n = 10,
    min_epv = 10,
    min_cohort_n = 100,
    min_gene_cov = 70,
    tumour_codes = c("Primary Tumor"),
    exclude_projects = c("TCGA-LAML", "TCGA-DLBC"),
    fdr_by_cohort = c("Dataset", "Project", "Test", "Metric"),
    fdr_by_pancan = c("Family", "Test"),
    cdr_file = "data/TCGA-CDR.csv",
    cdr_rfs_type = "PFI",
    signature_panel_path = "inst/signatures/panel.csv",
    composite_defs_path = "inst/signatures/composite_defs.csv",
    crc_modifiers = c("CMS", "PDS", "Stage_bin", "MSI_group"),
    cache_dir = "cache",
    global_seed = 42
  )

  from_yaml <- if (!is.null(config_path) && file.exists(config_path)) {
    yaml::read_yaml(config_path)
  } else {
    list()
  }

  cfg <- utils::modifyList(defaults, from_yaml)
  cfg <- utils::modifyList(cfg, overrides)

  # yaml::read_yaml() returns sequences as lists; the rest of the package
  # expects plain character vectors for these.
  vector_fields <- c("tumour_codes", "exclude_projects", "fdr_by_cohort",
                      "fdr_by_pancan", "crc_modifiers")
  for (f in vector_fields) cfg[[f]] <- as.character(unlist(cfg[[f]]))

  path_fields <- c("cdr_file", "signature_panel_path", "composite_defs_path",
                    "cache_dir")
  for (f in path_fields) {
    if (!is.null(cfg[[f]]) && !is_absolute_path(cfg[[f]])) {
      cfg[[f]] <- file.path(root, cfg[[f]])
    }
  }

  cfg$panel_hash <- NULL  # filled in by load_signature_panel() once the panel is read
  cfg
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[\\\\/]|~)", path)
}
