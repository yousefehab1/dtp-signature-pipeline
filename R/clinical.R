# TCGA Clinical Data Resource loading and the shared 3-year landmark endpoint
# derivation used by both the CRC and pan-cancer survival modules.

#' Load the TCGA Clinical Data Resource (Liu et al. 2018)
#'
#' Returns one row per patient with overall survival and the configured
#' recurrence endpoint (`cfg$cdr_rfs_type`: OS/DSS/DFI/PFI).
#'
#' @param cfg Config from [dtp_config()].
#' @export
load_tcga_cdr <- function(cfg) {
  path <- cfg$cdr_file
  rfs_type <- cfg$cdr_rfs_type
  stopifnot(file.exists(path))
  if (!rfs_type %in% c("OS", "DSS", "DFI", "PFI")) {
    stop("cfg$cdr_rfs_type must be one of OS/DSS/DFI/PFI")
  }

  rfs_time_col <- paste0(rfs_type, ".time")
  cdr <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                          na.strings = c("", "NA", "#N/A",
                                         "[Not Available]", "[Not Applicable]", "[Not Evaluated]"))
  colnames(cdr) <- trimws(gsub("^﻿", "", colnames(cdr)))

  required <- c("bcr_patient_barcode", "type", "OS", "OS.time", rfs_type, rfs_time_col)
  miss <- setdiff(required, colnames(cdr))
  if (length(miss)) stop("CDR missing columns: ", paste(miss, collapse = ", "))

  has_stage <- "ajcc_pathologic_tumor_stage" %in% colnames(cdr)

  neg_os <- sum(cdr$OS.time < 0, na.rm = TRUE)
  neg_rfs <- sum(cdr[[rfs_time_col]] < 0, na.rm = TRUE)
  if (neg_os > 0) warning(neg_os, " patients with negative OS.time", call. = FALSE)
  if (neg_rfs > 0) warning(neg_rfs, " patients with negative ", rfs_time_col, call. = FALSE)
  message(sprintf("CDR: %d patients, %d cancer types. Recurrence endpoint = %s (%d non-null).",
                   nrow(cdr), dplyr::n_distinct(cdr$type), rfs_type, sum(!is.na(cdr[[rfs_type]]))))

  cdr %>%
    dplyr::mutate(
      Patient_ID = toupper(trimws(bcr_patient_barcode)),
      Project_ID = paste0("TCGA-", toupper(trimws(type))),
      Stage = if (has_stage) {
        gsub("Stage |[ABC]$| ", "", ajcc_pathologic_tumor_stage)
      } else {
        NA_character_
      },
      OS_event = as.integer(OS),
      OS_months = as.numeric(OS.time) / 30.44,
      DFS_event = as.integer(.data[[rfs_type]]),
      DFS_months = as.numeric(.data[[rfs_time_col]]) / 30.44
    ) %>%
    dplyr::filter(!is.na(Patient_ID), Patient_ID != "") %>%
    dplyr::distinct(Patient_ID, .keep_all = TRUE)
}

#' Derive censoring-aware 3-year landmark endpoints
#'
#' Adds `Surv_3yr`/`Recurrence_3yr` factors (for Wilcoxon/violin grouping)
#' and `OS3Y_event`/`OS3Y_delay`/`RFS3Y_event`/`RFS3Y_delay` numeric columns
#' (for KM/Cox), landmarked at `os_cut`/`rfs_cut` months.
#'
#' @param df A data frame with numeric event/time columns.
#' @param os_event,os_time,rfs_event,rfs_time Column names in `df`.
#' @param os_cut,rfs_cut Landmark cutpoints in months (e.g. `cfg$os_cutpoint`).
#' @export
derive_endpoints <- function(df,
                              os_event = "OS_event", os_time = "OS_delay",
                              rfs_event = "Recurrence_event", rfs_time = "Recurrence_delay",
                              os_cut, rfs_cut) {
  df %>%
    dplyr::mutate(
      .oe = as.numeric(.data[[os_event]]), .ot = as.numeric(.data[[os_time]]),
      .re = as.numeric(.data[[rfs_event]]), .rt = as.numeric(.data[[rfs_time]]),

      Surv_3yr = factor(dplyr::case_when(
        !is.na(.ot) & !is.na(.oe) & .ot <= os_cut & .oe == 1 ~ "Dead_3yr",
        !is.na(.ot) & !is.na(.oe) ~ "Alive_3yr",
        TRUE ~ NA_character_),
        levels = c("Dead_3yr", "Alive_3yr")),

      Recurrence_3yr = factor(dplyr::case_when(
        .re == 1 & .rt <= rfs_cut ~ "Recurred",
        .re == 1 & .rt > rfs_cut ~ "Recurrence-Free",
        .re == 0 & .rt >= rfs_cut ~ "Recurrence-Free",
        TRUE ~ NA_character_),
        levels = c("Recurred", "Recurrence-Free")),

      OS3Y_event = dplyr::case_when(is.na(.oe) | is.na(.ot) ~ NA_real_,
                                     .ot > os_cut | .oe == 0 ~ 0, TRUE ~ 1),
      OS3Y_delay = ifelse(is.na(.ot), NA_real_, ifelse(.ot >= os_cut, os_cut, .ot)),

      RFS3Y_event = dplyr::case_when(is.na(.re) | is.na(.rt) ~ NA_real_,
                                      .rt > rfs_cut | .re == 0 ~ 0, TRUE ~ 1),
      RFS3Y_delay = ifelse(is.na(.rt), NA_real_, ifelse(.rt >= rfs_cut, rfs_cut, .rt))
    ) %>%
    dplyr::select(-.oe, -.ot, -.re, -.rt)
}
