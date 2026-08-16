# Statistical tests and survival analysis. Every function that gates on sample
# size or event count takes an explicit `cfg` argument from `dtp_config()`
# rather than reading bare globals.

#' Get survival time and event column names for a metric
#'
#' @param metric "OS" or "RFS".
#' @return Named character vector with elements `t` and `e`.
#' @export
metric_cols <- function(metric) {
  if (metric == "OS") {
    c(t = "OS3Y_delay", e = "OS3Y_event")
  } else {
    c(t = "RFS3Y_delay", e = "RFS3Y_event")
  }
}

#' Calculate Wilcoxon rank-sum test statistics and effect size
#'
#' @param df Data frame containing score and group columns.
#' @param score_col Character name of the numeric score column.
#' @param group_col Character name of the binary group column.
#' @param cfg Config from [dtp_config()].
#' @return List with `p` (p-value), `r` (rank biserial correlation), and `n` (sample size).
#' @export
get_wilcox_stats <- function(df, score_col, group_col, cfg) {
  d <- df %>% dplyr::filter(!is.na(.data[[group_col]]) & .data[[group_col]] != "")
  if (!is.numeric(d[[score_col]])) return(list(p = NA, r = NA, n = NA))
  tbl <- table(droplevels(factor(d[[group_col]])))
  if (length(tbl) < 2 || any(tbl < cfg$min_group_n)) return(list(p = NA, r = NA, n = NA))
  lvls <- names(tbl)
  g1 <- d[[score_col]][d[[group_col]] == lvls[1]]
  g2 <- d[[score_col]][d[[group_col]] == lvls[2]]
  wt <- suppressWarnings(stats::wilcox.test(g1, g2, exact = FALSE))
  r  <- as.numeric(1 - 2 * wt$statistic / (length(g1) * length(g2)))
  list(p = wt$p.value, r = r, n = length(g1) + length(g2))
}

#' Calculate Kruskal-Wallis test statistics and epsilon-squared effect size
#'
#' @param df Data frame containing score and group columns.
#' @param score_col Character name of the numeric score column.
#' @param group_col Character name of the multi-level group column.
#' @param cfg Config from [dtp_config()].
#' @return List with `p`, `eps2` (epsilon-squared), `n`, and `k` (number of groups).
#' @export
get_kruskal_stats <- function(df, score_col, group_col, cfg) {
  d <- df %>% dplyr::filter(!is.na(.data[[group_col]]) & .data[[group_col]] != "" &
                            !is.na(.data[[score_col]]))
  if (!is.numeric(d[[score_col]])) return(list(p = NA, eps2 = NA, n = NA, k = NA))
  tbl  <- table(droplevels(factor(d[[group_col]])))
  keep <- names(tbl)[tbl >= cfg$min_group_n]
  d    <- d[d[[group_col]] %in% keep, , drop = FALSE]
  if (length(keep) < 2) return(list(p = NA, eps2 = NA, n = nrow(d), k = length(keep)))
  kt <- tryCatch(stats::kruskal.test(d[[score_col]], factor(d[[group_col]])),
                 error = function(e) NULL)
  if (is.null(kt)) return(list(p = NA, eps2 = NA, n = nrow(d), k = length(keep)))
  n    <- nrow(d)
  eps2 <- as.numeric(kt$statistic) / ((n^2 - 1) / (n + 1))
  list(p = kt$p.value, eps2 = eps2, n = n, k = length(keep))
}

#' Calculate Kaplan-Meier log-rank test p-value (median split)
#'
#' @param df Data frame containing score and survival time/event columns.
#' @param score_col Character name of the numeric score column.
#' @param t_col Character name of the follow-up time column.
#' @param e_col Character name of the binary event status column (1=event, 0=censored).
#' @param cfg Config from [dtp_config()].
#' @return Numeric log-rank p-value, or NA if gated/untestable.
#' @export
get_km_pval <- function(df, score_col, t_col, e_col, cfg) {
  d <- df[!is.na(df[[t_col]]) & !is.na(df[[e_col]]) & !is.na(df[[score_col]]), ]
  if (nrow(d) == 0) return(NA)
  d$Group <- factor(ifelse(d[[score_col]] > stats::median(d[[score_col]], na.rm = TRUE), "High", "Low"),
                    levels = c("Low", "High"))
  if (length(table(d$Group)) < 2 || any(table(d$Group) < cfg$min_km_group_n) ||
      sum(d[[e_col]]) < cfg$min_events) return(NA)
  tryCatch({
    sd <- survival::survdiff(stats::as.formula(paste0("Surv(", t_col, ",", e_col, ") ~ Group")), data = d)
    1 - stats::pchisq(sd$chisq, length(sd$n) - 1)
  }, error = function(e) NA)
}

#' Calculate univariate Cox proportional hazards statistics
#'
#' @param df Data frame containing score and survival time/event columns.
#' @param score_col Character name of the numeric score column.
#' @param t_col Character name of the follow-up time column.
#' @param e_col Character name of the binary event status column (1=event, 0=censored).
#' @param cfg Config from [dtp_config()].
#' @return List with `HR`, `HR_lower`, `HR_upper`, `P`, `C_index`, `N`, `N_events`, or NULL if gated/failed.
#' @export
get_cox_stats <- function(df, score_col, t_col, e_col, cfg) {
  d <- df[!is.na(df[[t_col]]) & !is.na(df[[e_col]]) & !is.na(df[[score_col]]), ]
  if (nrow(d) < cfg$min_cox_n || sum(d[[e_col]]) < cfg$min_events) return(NULL)
  tryCatch({
    fit <- survival::coxph(stats::as.formula(paste0("Surv(", t_col, ",", e_col, ") ~ ", score_col)), data = d)
    s   <- summary(fit)
    list(HR = exp(stats::coef(fit)[[1]]),
         HR_lower = s$conf.int[1, "lower .95"], HR_upper = s$conf.int[1, "upper .95"],
         P = s$coefficients[1, "Pr(>|z|)"], C_index = unname(s$concordance["C"]),
         N = nrow(d), N_events = sum(d[[e_col]]))
  }, error = function(e) NULL)
}

#' Harmonise CRC clinical modifiers
#'
#' Harmonises the four CRC modifiers into clean, unified factor columns so one
#' code path serves both cohorts.
#'
#' @param df Clinical data frame.
#' @param cohort "GSE39582" or "TCGA-COAD".
#' @return Data frame with harmonised modifier columns (`Stage_bin`, `MSI_group`, `CMS`, `PDS`).
#' @export
harmonize_crc_modifiers <- function(df, cohort) {
  stage_raw <- if (cohort == "GSE39582") suppressWarnings(as.numeric(df$TNM_stage))
               else df$Stage
  df$Stage_bin <- if (cohort == "GSE39582")
    factor(dplyr::case_when(stage_raw %in% 0:2 ~ "Early",
                            stage_raw %in% 3:4 ~ "Late", TRUE ~ NA_character_),
           levels = c("Early", "Late"))
  else
    factor(dplyr::case_when(df$Stage %in% c("I", "II")   ~ "Early",
                            df$Stage %in% c("III", "IV")  ~ "Late", TRUE ~ NA_character_),
           levels = c("Early", "Late"))

  msi_raw <- if (cohort == "GSE39582") df$MMR else df[["paper_MSI_status"]]
  # MSI-L is grouped with MSS, not MSI: MSI-L lacks the hypermutator/immune phenotype that
  # defines the distinct MSI-high (dMMR) entity and behaves like MSS (Rantanen 2023 pooled
  # MSS/MSI-low vs MSI-high). This also makes GSE (dMMR/pMMR only) and TCGA consistent.
  df$MSI_group <- factor(dplyr::case_when(
    msi_raw %in% c("pMMR", "MSS", "MSI-L") ~ "MSS",
    msi_raw %in% c("dMMR", "MSI-H")        ~ "MSI",
    TRUE                                   ~ NA_character_), levels = c("MSS", "MSI"))

  if ("CMS" %in% colnames(df)) {
    cms <- df$CMS; cms[!cms %in% paste0("CMS", 1:4)] <- NA
    df$CMS <- stats::relevel(factor(cms), ref = "CMS2")
  }
  if ("PDS" %in% colnames(df)) {
    pds <- df$PDS; pds[!pds %in% paste0("PDS", 1:3)] <- NA
    df$PDS <- droplevels(factor(pds))
    if (nlevels(df$PDS) > 0) {
      ref <- names(which.max(table(df$PDS)))
      df$PDS <- stats::relevel(df$PDS, ref = ref)
    }
  }
  df
}

# Second-row p-value from an anova.coxph LRT table, robust to the p-column name. Internal.
.anova_lrt_p <- function(reduced, full) {
  av   <- stats::anova(reduced, full)
  pcol <- grep("Chi", colnames(av), value = TRUE)
  if (length(pcol) == 0) return(NA_real_)
  av[[utils::tail(pcol, 1)]][2]
}

# Events-per-variable gate. Internal.
.cox_epv_ok <- function(n, n_events, n_params, cfg) {
  n >= cfg$min_cox_n && n_events >= cfg$min_events && n_events >= cfg$min_epv * n_params
}

# Drop rows missing any modelled column, then drop unused factor levels. Internal.
.cox_complete <- function(df, cols) {
  d <- df[stats::complete.cases(df[, cols, drop = FALSE]), , drop = FALSE]
  droplevels(d)
}

#' Calculate multivariable adjusted Cox proportional hazards statistics
#'
#' Fits a multivariable Cox model for `score_col` adjusted for `covars`.
#' Returns the score-term HR/CI/p plus a model LRT (testing added signal
#' beyond covariates) and proportional hazards test p-value (`cox.zph`).
#'
#' @param df Data frame containing score, time, event, and covariate columns.
#' @param score_col Character name of the numeric score column.
#' @param t_col Character name of the follow-up time column.
#' @param e_col Character name of the binary event status column.
#' @param covars Character vector of covariate column names.
#' @param cfg Config from [dtp_config()].
#' @return List of adjusted Cox statistics, or NULL if inadmissible under EPV gate.
#' @export
get_cox_adjusted <- function(df, score_col, t_col, e_col, covars, cfg) {
  cols <- c(score_col, t_col, e_col, covars)
  d <- .cox_complete(df, cols)
  covars <- covars[vapply(covars, function(c) length(unique(d[[c]])) > 1, logical(1))]
  n_params <- length(covars) + 1L +
    sum(vapply(covars, function(c) if (is.factor(d[[c]])) nlevels(d[[c]]) - 2L else 0L, integer(1)))
  n_ev <- sum(d[[e_col]])
  if (!.cox_epv_ok(nrow(d), n_ev, n_params, cfg)) return(NULL)
  rhs  <- paste(c(score_col, covars), collapse = " + ")
  base <- paste(covars, collapse = " + ")
  tryCatch({
    full <- survival::coxph(stats::as.formula(sprintf("Surv(%s,%s) ~ %s", t_col, e_col, rhs)),
                            data = d, ties = "efron")
    s    <- summary(full)
    lrt  <- if (length(covars) == 0) unname(s$logtest[["pvalue"]]) else {
      red <- survival::coxph(stats::as.formula(sprintf("Surv(%s,%s) ~ %s", t_col, e_col, base)), data = d)
      .anova_lrt_p(red, full)
    }
    ph <- tryCatch(survival::cox.zph(full)$table[score_col, "p"], error = function(e) NA_real_)
    adj_hr <- exp(stats::coef(full)[[score_col]])
    unadj  <- survival::coxph(stats::as.formula(sprintf("Surv(%s,%s) ~ %s", t_col, e_col, score_col)), data = d)
    unadj_hr <- exp(stats::coef(unadj)[[1]])
    delta_loghr <- log(adj_hr) - log(unadj_hr)
    delta_pct   <- 100 * delta_loghr / max(abs(log(unadj_hr)), 0.05)
    list(HR = adj_hr,
         HR_lower = s$conf.int[score_col, "lower .95"],
         HR_upper = s$conf.int[score_col, "upper .95"],
         P = s$coefficients[score_col, "Pr(>|z|)"], C_index = unname(s$concordance["C"]),
         LRT_P = lrt, PH_P = ph, Unadj_HR = unadj_hr,
         Delta_logHR = delta_loghr, Delta_logHR_pct = delta_pct,
         N = nrow(d), N_events = n_ev)
  }, error = function(e) NULL)
}

#' Calculate Cox proportional hazards interaction statistics
#'
#' Tests effect modification (`score * modifier` vs `score + modifier` LRT)
#' and fits stratified models for each level of `modifier`.
#'
#' @param df Data frame containing score, time, event, and modifier columns.
#' @param score_col Character name of the numeric score column.
#' @param t_col Character name of the follow-up time column.
#' @param e_col Character name of the binary event status column.
#' @param modifier Character name of the modifier factor column.
#' @param cfg Config from [dtp_config()].
#' @return List with `Interaction_P`, `K`, `N`, `N_events`, and `per_level` data frame, or NULL.
#' @export
get_cox_interaction <- function(df, score_col, t_col, e_col, modifier, cfg) {
  cols <- c(score_col, t_col, e_col, modifier)
  d <- .cox_complete(df, cols)
  if (!is.factor(d[[modifier]])) d[[modifier]] <- factor(d[[modifier]])
  d <- droplevels(d)
  k <- nlevels(d[[modifier]])
  if (k < 2) return(NULL)
  n_params <- 2L * k - 1L
  n_ev <- sum(d[[e_col]])
  if (!.cox_epv_ok(nrow(d), n_ev, n_params, cfg)) return(NULL)
  tryCatch({
    red  <- survival::coxph(stats::as.formula(sprintf("Surv(%s,%s) ~ %s + %s", t_col, e_col, score_col, modifier)), data = d)
    full <- survival::coxph(stats::as.formula(sprintf("Surv(%s,%s) ~ %s * %s", t_col, e_col, score_col, modifier)), data = d)
    lrt  <- .anova_lrt_p(red, full)
    per_level <- lapply(levels(d[[modifier]]), function(lv) {
      cs <- get_cox_stats(d[d[[modifier]] == lv, , drop = FALSE], score_col, t_col, e_col, cfg)
      if (is.null(cs)) return(NULL)
      data.frame(Level = lv, HR = cs$HR, HR_lower = cs$HR_lower, HR_upper = cs$HR_upper,
                 N = cs$N, N_events = cs$N_events, stringsAsFactors = FALSE)
    })
    list(Interaction_P = lrt, K = k, N = nrow(d), N_events = n_ev,
         per_level = do.call(rbind, per_level))
  }, error = function(e) NULL)
}

# One tidy stats row. Internal.
.stat_row <- function(dataset, proj, test, metric, score,
                      raw_p = NA, effect_r = NA, hr = NA, hr_l = NA, hr_u = NA,
                      c_index = NA, n = NA, n_events = NA, is_testable = !is.na(raw_p)) {
  data.frame(Dataset = dataset, Project = proj, Test = test, Metric = metric, Score = score,
             Raw_P = raw_p, Effect_r = effect_r, HR = hr, HR_lower = hr_l, HR_upper = hr_u,
             C_index = c_index, N = n, N_events = n_events, Is_Testable = is_testable,
             stringsAsFactors = FALSE)
}

#' Run Wilcoxon, KM, and Cox tests for a single score and survival metric
#'
#' @param df Data frame containing clinical and score data.
#' @param score_col Column name of the score to test.
#' @param metric "OS" or "RFS".
#' @param dataset Name of dataset (e.g. "TCGA", "GSE39582").
#' @param proj Project or cohort name.
#' @param score_name Display name of the score.
#' @param group_col Column name for Wilcoxon 3-year binary grouping.
#' @param cfg Config from [dtp_config()].
#' @return A data frame with up to 3 tidy result rows (Wilcoxon, KM, Cox).
#' @export
run_all_stats <- function(df, score_col, metric, dataset, proj, score_name, group_col, cfg) {
  mc <- metric_cols(metric)
  ws <- get_wilcox_stats(df, score_col, group_col, cfg)
  km <- get_km_pval(df, score_col, mc[["t"]], mc[["e"]], cfg)
  cs <- get_cox_stats(df, score_col, mc[["t"]], mc[["e"]], cfg)
  rows <- list(
    .stat_row(dataset, proj, "Wilcoxon", metric, score_name,
              raw_p = ws$p, effect_r = ws$r, n = ws$n, is_testable = !is.na(ws$p)),
    .stat_row(dataset, proj, "KM", metric, score_name, raw_p = km, is_testable = !is.na(km)),
    if (is.null(cs))
      .stat_row(dataset, proj, "Cox", metric, score_name, is_testable = FALSE)
    else
      .stat_row(dataset, proj, "Cox", metric, score_name,
                raw_p = cs$P, hr = cs$HR, hr_l = cs$HR_lower, hr_u = cs$HR_upper,
                c_index = cs$C_index, n = cs$N, n_events = cs$N_events)
  )
  do.call(rbind, rows)
}

#' Run survival statistics across cohorts and score columns
#'
#' @param cohorts Named list of data frames (one per cohort).
#' @param dataset Dataset identifier.
#' @param score_cols Character vector of score column names.
#' @param cfg Config from [dtp_config()].
#' @param group_os Column name for 3-year overall survival binary grouping.
#' @param group_rfs Column name for 3-year recurrence-free survival binary grouping.
#' @return A tidy data frame with raw (un-FDR-corrected) statistical test results.
#' @export
run_survival_block <- function(cohorts, dataset, score_cols, cfg,
                               group_os = "Surv_3yr", group_rfs = "Recurrence_3yr") {
  rows <- list()
  for (cn in names(cohorts)) {
    cd <- cohorts[[cn]]
    for (sc in score_cols) {
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "OS",  dataset, cn, sc, group_os, cfg)
      rows[[length(rows) + 1]] <- run_all_stats(cd, sc, "RFS", dataset, cn, sc, group_rfs, cfg)
    }
  }
  do.call(rbind, rows)
}

#' Apply Benjamini-Hochberg FDR correction within a named grouping
#'
#' The grouping is a closed enum, not a free-form `by=` vector: `family = "cohort"` groups by
#' `cfg$fdr_by_cohort` (used by the single-cohort CRC arm), `family = "pancan"` groups by
#' `cfg$fdr_by_pancan` (used by the cross-cohort pan-cancer arm). The two groupings are deliberately
#' different -- they answer different statistical questions -- but making the choice an enum backed by
#' cfg instead of an ad-hoc argument at each call site means no call site can silently pass a third,
#' inconsistent grouping.
#'
#' @param stats_df A stats table from [run_all_stats()]/[run_survival_block()] or equivalent.
#' @param cfg Config from [dtp_config()].
#' @param family "cohort" or "pancan" -- selects `cfg$fdr_by_cohort`/`cfg$fdr_by_pancan`.
#' @param p_col P-value column BH acts on. Most families use "Raw_P"; the adjusted-Cox (confounding)
#'   family should pass `p_col = "LRT_score_P"` so its FDR is on the likelihood-ratio test of score vs.
#'   the covariate-only model, matching the interaction family's use of an LRT rather than a Wald p.
#' @export
apply_fdr <- function(stats_df, cfg, family = c("cohort", "pancan"), p_col = "Raw_P") {
  family <- match.arg(family)
  by <- if (family == "cohort") cfg$fdr_by_cohort else cfg$fdr_by_pancan
  stats_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::mutate(FDR_P = stats::p.adjust(.data[[p_col]], method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Is_Significant = !is.na(FDR_P) & FDR_P < 0.05)
}
