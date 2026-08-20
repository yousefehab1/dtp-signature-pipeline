# ==============================================================================
# R/pancan_composites.R - Pan-Cancer Survival Presentation Layer (Fig 6, 7, App)
#
# Generates publication-ready figures and summary tables directly from in-memory
# analysis tibbles (no CSV reads inside builders).
# ==============================================================================

FOREST_XLIM <- c(0.1, 10)

.PANCAN_DIR_PAL <- c("Higher hazard (HR>1)" = "#E64B35", "Lower hazard (HR<1)" = "#3B6EA5")
.PANCAN_CORR_LAB <- c(
  Up_ssGSEA_Corrected        = "DTP Up",
  Down_ssGSEA_Corrected      = "DTP Down",
  Composite_ssGSEA_Corrected = "DTP Composite"
)
.PANCAN_HIGHLIGHT <- "TCGA-COAD"
.PANCAN_ENDPT_LAB <- c(OS = "Overall survival", RFS = "Recurrence")
.PANCAN_METRIC_EVT <- c(OS = "OS3Y_event", RFS = "RFS3Y_event")
.PANCAN_METRIC_TIM <- c(OS = "OS3Y_delay", RFS = "RFS3Y_delay")

# ---- helpers -----------------------------------------------------------------

.pancan_num <- function(x) suppressWarnings(as.numeric(x))

.fmt_p <- function(prefix, p) {
  if (length(p) != 1 || is.na(p)) {
    NA_character_
  } else {
    paste0(prefix, if (p < 0.001) "<.001" else sub("^0", "", sprintf("=%.3f", p)))
  }
}

# Within-cohort SD of one score over the same subset get_cox_stats() models:
# non-NA score AND non-NA truncated follow-up time AND event for that metric.
.pancan_sd <- function(master, project, score, metric) {
  ev <- .PANCAN_METRIC_EVT[[metric]]
  tm <- .PANCAN_METRIC_TIM[[metric]]
  d  <- master[master$Project_ID == project, , drop = FALSE]
  if (!all(c(score, tm, ev) %in% names(d))) return(NA_real_)
  v  <- .pancan_num(d[[score]])
  t  <- .pancan_num(d[[tm]])
  e  <- .pancan_num(d[[ev]])
  ok <- !is.na(v) & !is.na(t) & !is.na(e)
  if (sum(ok) < 2L) return(NA_real_)
  .score_sd(v[ok])
}

# Aggregate pooled SD of a score over the non-NA follow-up and event subset.
.pancan_agg_sd <- function(master, score, metric) {
  ev <- .PANCAN_METRIC_EVT[[metric]]
  tm <- .PANCAN_METRIC_TIM[[metric]]
  if (!all(c(score, tm, ev) %in% names(master))) return(NA_real_)
  v  <- .pancan_num(master[[score]])
  t  <- .pancan_num(master[[tm]])
  e  <- .pancan_num(master[[ev]])
  ok <- !is.na(v) & !is.na(t) & !is.na(e)
  if (sum(ok) < 2L) return(NA_real_)
  .score_sd(v[ok])
}

# Add per-SD HR columns (point + both bounds) + direction + off-scale squish.
# Non-convergent fits (non-finite per-SD HR or CI) are NOT squished to the
# FOREST_XLIM boundary - display coords stay NA so geoms drop them (na.rm).
.pancan_rescale <- function(d, sd_vec) {
  hr <- .pancan_num(d$HR)
  lo <- .pancan_num(d$HR_lower)
  hi <- .pancan_num(d$HR_upper)
  ok <- !is.na(hr) & !is.na(sd_vec)
  d$hr_sd <- ifelse(ok, exp(log(hr) * sd_vec), NA_real_)
  d$lo_sd <- ifelse(ok, exp(log(lo) * sd_vec), NA_real_)
  d$hi_sd <- ifelse(ok, exp(log(hi) * sd_vec), NA_real_)
  nonconv <- (!is.na(d$hr_sd) & (!is.finite(d$hr_sd) | d$hr_sd <= 0)) |
             (!is.na(d$lo_sd) & (!is.finite(d$lo_sd) | d$lo_sd <= 0)) |
             (!is.na(d$hi_sd) & (!is.finite(d$hi_sd) | d$hi_sd <= 0))
  if (any(nonconv)) {
    message(sprintf("[pancan forest] %d non-convergent Cox row(s) not squished to axis boundary (non-finite per-SD HR or CI)",
                    sum(nonconv)))
  }
  d$dir <- ifelse(is.na(d$hr_sd) | nonconv, NA_character_,
                  ifelse(d$hr_sd > 1, "Higher hazard (HR>1)", "Lower hazard (HR<1)"))
  d$offscale <- !nonconv & !is.na(d$hr_sd) &
                (d$hr_sd < FOREST_XLIM[1] | d$hr_sd > FOREST_XLIM[2])
  d$hr_d <- ifelse(nonconv, NA_real_, scales::squish(d$hr_sd, FOREST_XLIM))
  d$lo_d <- ifelse(nonconv, NA_real_, scales::squish(d$lo_sd, FOREST_XLIM))
  d$hi_d <- ifelse(nonconv, NA_real_, scales::squish(d$hi_sd, FOREST_XLIM))
  d$nonconv <- nonconv
  d
}

# Cohort row order shared by panels A, B, C: ascending Up OS per-SD HR,
# not-testable OS sorted to the bottom via na.last = FALSE.
.pancan_cohort_order <- function(stats_df, master, order_score = "Up_ssGSEA") {
  d <- stats_df %>%
    dplyr::filter(Family == "PerCohort", Test == "Cox",
                  Score == order_score, Metric == "OS")
  if (!nrow(d)) {
    warning("[pancan composites] no PerCohort/Cox/OS rows for ordering score '",
            order_score, "'; cohort row order falls back to input order and is NOT ",
            "the ascending Up-OS per-SD HR order the figure claims.", call. = FALSE)
    return(unique(stats_df$Project[stats_df$Family == "PerCohort"]))
  }
  sd_vec <- mapply(function(pr, mt) .pancan_sd(master, pr, order_score, mt),
                   d$Project, d$Metric)
  d <- .pancan_rescale(d, sd_vec)
  ord <- d$Project[order(d$hr_sd, na.last = FALSE)]
  unique(c(ord, setdiff(unique(stats_df$Project[stats_df$Family == "PerCohort"]), ord)))
}

# ==============================================================================
# PER-COHORT FOREST (Figure 6 panels A, B, C)
# ==============================================================================

#' Build per-cohort pan-cancer survival forest plot for a single score
#'
#' Generates a univariable continuous-score Cox proportional hazards forest plot
#' across all testable TCGA cohorts (per 1 SD), faceted by survival endpoint
#' (OS and RFS), highlighting TCGA-COAD.
#'
#' @param stats_df FDR stats summary tibble from [run_pancan_survival()].
#' @param master Scored master clinical tibble from [run_pancan_survival()].
#' @param score Name of the score column to plot (e.g. "Up_ssGSEA").
#' @param score_lab Human-readable score label (e.g. "DTP Up").
#' @param order_score Name of the score column used to sort cohort rows (defaults to `score`).
#' @return A ggplot object.
#' @export
build_pancan_forest <- function(stats_df, master, score = "Up_ssGSEA",
                                score_lab = "DTP Up", order_score = score) {
  d <- stats_df %>%
    dplyr::filter(Family == "PerCohort", Test == "Cox", Score == score, Metric %in% c("OS", "RFS")) %>%
    dplyr::mutate(
      Is_Testable = as.logical(Is_Testable),
      Is_Significant = as.logical(Is_Significant),
      FDR_P = .pancan_num(FDR_P),
      Raw_P = .pancan_num(Raw_P),
      N = .pancan_num(N),
      N_events = .pancan_num(N_events)
    )

  if (nrow(d) == 0) {
    warning("[pancan forest] no PerCohort Cox rows found for score '", score, "'.", call. = FALSE)
    return(ggplot2::ggplot() + ggplot2::labs(title = paste(score_lab, "- no data")))
  }

  sd_vec <- mapply(function(pr, mt) .pancan_sd(master, pr, score, mt), d$Project, d$Metric)
  d <- .pancan_rescale(d, sd_vec)
  d$ScoreLab <- factor(score_lab)

  cohort_levels <- .pancan_cohort_order(stats_df, master, order_score)
  d$Project <- factor(d$Project, levels = cohort_levels)
  d$Metric  <- factor(d$Metric, levels = c("OS", "RFS"))

  coad_pos <- match(.PANCAN_HIGHLIGHT, levels(d$Project))
  lab_cols <- ifelse(levels(d$Project) == .PANCAN_HIGHLIGHT, "#B2182B", "grey20")
  lab_face <- ifelse(levels(d$Project) == .PANCAN_HIGHLIGHT, "bold", "plain")

  d$ne_clean <- sprintf("n=%d, e=%d", as.integer(d$N), as.integer(d$N_events))
  d$ne_full  <- ifelse(is.na(d$hr_sd), d$ne_clean,
                       sprintf("%s\nHR=%.2f [%.2f-%.2f]", d$ne_clean, d$hr_sd, d$lo_sd, d$hi_sd))
  d$rm_full  <- mapply(function(rp, fp) {
    paste(stats::na.omit(c(.fmt_p("p", rp),
                           paste0(.fmt_p("q", fp), .sig_stars(fp, "")))), collapse = "  ")
  }, d$Raw_P, d$FDR_P)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = hr_d, y = Project))
  if (!is.na(coad_pos)) {
    p <- p + ggplot2::annotate("rect", xmin = FOREST_XLIM[1], xmax = FOREST_XLIM[2],
                               ymin = coad_pos - 0.5, ymax = coad_pos + 0.5,
                               fill = "#FDE7C9", alpha = 0.7)
  }
  p <- p +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo_d, xmax = hi_d, colour = dir),
                           width = 0.3, linewidth = 0.5, orientation = "y", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = dir, shape = offscale), size = 2.0, na.rm = TRUE) +
    ggplot2::geom_text(data = subset(d, Is_Testable & !is.na(N)), inherit.aes = FALSE,
                       ggplot2::aes(x = FOREST_XLIM[1], y = Project, label = ne_full),
                       hjust = 0, vjust = -0.6, size = 1.9, colour = "grey45", lineheight = 0.82) +
    ggplot2::geom_text(data = subset(d, !Is_Testable), inherit.aes = FALSE,
                       ggplot2::aes(x = 1, y = Project, label = "n/t"),
                       colour = "grey55", size = 2.6, fontface = "italic") +
    ggplot2::geom_text(data = subset(d, Is_Testable & !is.na(FDR_P)), inherit.aes = FALSE,
                       ggplot2::aes(x = FOREST_XLIM[2], y = Project, label = rm_full),
                       hjust = 1, vjust = -0.6, size = 1.7, colour = "grey30") +
    ggplot2::facet_grid(ScoreLab ~ Metric, labeller = ggplot2::labeller(Metric = .PANCAN_ENDPT_LAB)) +
    ggplot2::scale_x_log10(limits = FOREST_XLIM, breaks = c(0.1, 0.3, 1, 3, 10),
                           labels = c("0.1", "0.3", "1", "3", "10")) +
    ggplot2::scale_y_discrete(labels = function(x) ifelse(x == .PANCAN_HIGHLIGHT, paste0(x, " (CRC)"), x)) +
    ggplot2::scale_colour_manual(values = .PANCAN_DIR_PAL, name = "Score effect",
                                 na.value = "grey70", na.translate = FALSE, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21),
                                labels = c(`FALSE` = "in range", `TRUE` = "off-scale (squished)"),
                                name = NULL) +
    ggplot2::labs(
      title = sprintf("%s score vs survival across %d TCGA cohorts", score_lab, nlevels(d$Project)),
      subtitle = paste0(
        "Univariable continuous-score Cox HR per 1 SD within each cohort (95% CI); x truncated to [0.1, 10].\n",
        "Left: n/e + per-SD HR [CI]; right: raw p, FDR q (BH); * FDR<0.05; grey n/t = not testable. TCGA-COAD (CRC) highlighted."
      ),
      x = "Hazard ratio per 1 SD (log scale)", y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7.5),
      panel.grid.major.x = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
      panel.spacing.x = grid::unit(9, "pt")
    )

  p
}

# ==============================================================================
# BATCH-CORRECTED AGGREGATE FOREST
# ==============================================================================

#' Build pan-cancer batch-corrected aggregate Cox forest plot
#'
#' Evaluates batch-corrected continuous Cox proportional hazards models
#' pooled across all solid-tumor cohorts for Up, Down, and Composite scores.
#'
#' @param stats_df FDR stats summary tibble from [run_pancan_survival()].
#' @param master Scored master clinical tibble from [run_pancan_survival()].
#' @param cfg Config from [dtp_config()].
#' @param composite_defs Composite definitions table.
#' @return A ggplot object.
#' @export
build_pancan_aggregate <- function(stats_df, master, cfg = dtp_config(), composite_defs = NULL) {
  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  def_row <- if (!is.null(composite_defs)) composite_defs[composite_defs$Is_Default == TRUE, , drop = FALSE] else NULL
  pos_sig <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Positive_Signature[1] else "Up"
  neg_sig <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Negative_Signature[1] else "Down"
  comp_id <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Composite_ID[1] else "Composite"

  up_raw  <- paste0(pos_sig, suffix)
  dn_raw  <- paste0(neg_sig, suffix)
  cmp_raw <- paste0(comp_id, suffix)

  up_corr  <- paste0(up_raw, "_Corrected")
  dn_corr  <- paste0(dn_raw, "_Corrected")
  cmp_corr <- paste0(cmp_raw, "_Corrected")

  raw_cols  <- c(up_raw, dn_raw, cmp_raw)
  corr_cols <- c(up_corr, dn_corr, cmp_corr)

  if (!all(corr_cols %in% colnames(master))) {
    sm <- t(as.matrix(sapply(master[, raw_cols, drop = FALSE], .pancan_num)))
    corrected <- limma::removeBatchEffect(sm, batch = master$Project_ID)
    corr_df <- as.data.frame(t(corrected))
    colnames(corr_df) <- corr_cols
    master <- dplyr::bind_cols(master, corr_df)
  }

  score_lab_map <- setNames(c("DTP Up", "DTP Down", "DTP Composite"), corr_cols)

  d <- stats_df %>%
    dplyr::filter(Family == "PanCancer", Test == "Cox", Metric %in% c("OS", "RFS"), Score %in% corr_cols) %>%
    dplyr::mutate(
      C_index = .pancan_num(C_index),
      FDR_P = .pancan_num(FDR_P),
      Raw_P = .pancan_num(Raw_P),
      Is_Significant = as.logical(Is_Significant),
      ScoreLab = score_lab_map[Score]
    )

  sd_vec <- mapply(function(sc, mt) .pancan_agg_sd(master, sc, mt), d$Score, d$Metric)
  d <- .pancan_rescale(d, sd_vec)
  d$Metric   <- factor(d$Metric, levels = c("OS", "RFS"))
  d$ScoreLab <- factor(d$ScoreLab, levels = rev(c("DTP Up", "DTP Down", "DTP Composite")))

  # Direction-consistency tally for the Up score (from the per-cohort family).
  up <- stats_df %>%
    dplyr::filter(Family == "PerCohort", Test == "Cox", Score == up_raw) %>%
    dplyr::mutate(Is_Testable = as.logical(Is_Testable),
                  Is_Significant = as.logical(Is_Significant),
                  HR = .pancan_num(HR))
  tal <- function(mt) {
    x <- up %>% dplyr::filter(Metric == mt, Is_Testable, !is.na(HR))
    list(gt = sum(x$HR > 1), n = nrow(x), sig = sum(x$Is_Significant, na.rm = TRUE))
  }
  os <- tal("OS")
  rf <- tal("RFS")
  tally_str <- sprintf("DTP Up direction across cohorts - HR > 1: %d of %d (OS), %d of %d (RFS); FDR-significant (OS): %d",
                       os$gt, os$n, rf$gt, rf$n, os$sig)

  lab_x <- FOREST_XLIM[2]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = hr_d, y = ScoreLab)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo_d, xmax = hi_d, colour = dir),
                           width = 0.22, linewidth = 0.6, orientation = "y", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = dir, shape = offscale), size = 2.8, na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = lab_x,
        label = sprintf("HR %.2f  C=%.2f\n%s  %s%s", hr_sd, C_index,
                        ifelse(is.na(Raw_P), "", ifelse(Raw_P < 0.001, "p<.001", sprintf("p=%.3f", Raw_P))),
                        ifelse(FDR_P < 0.001, "q<.001", sprintf("q=%.3f", FDR_P)),
                        .sig_stars(FDR_P, ""))
      ),
      hjust = 1, vjust = -0.6, size = 2.4, colour = "grey20", lineheight = 0.82
    ) +
    ggplot2::facet_grid(. ~ Metric, labeller = ggplot2::labeller(Metric = .PANCAN_ENDPT_LAB)) +
    ggplot2::scale_x_log10(limits = FOREST_XLIM, breaks = c(0.1, 0.3, 1, 3, 10),
                           labels = c("0.1", "0.3", "1", "3", "10")) +
    ggplot2::scale_colour_manual(values = .PANCAN_DIR_PAL, name = "Score effect", drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
    ggplot2::labs(
      title = "Pan-cancer aggregate (cancer type as batch)",
      subtitle = "Pooled continuous-score Cox HR per 1 SD of the batch-corrected score (95% CI). C = concordance; cells add raw p, FDR q.",
      caption = tally_str,
      x = "Hazard ratio per 1 SD (log scale)", y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "grey25"),
      panel.grid.major.x = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
      panel.spacing.x = grid::unit(9, "pt")
    )

  p
}

# ==============================================================================
# TABLE (Table 3.5 pan-cancer survival)
# ==============================================================================

#' Build Table 3.5 pan-cancer survival summary tibble
#'
#' Rescales univariable Cox proportional hazards estimates to per-1-SD HRs
#' for all per-cohort and pooled batch-corrected aggregate models.
#'
#' @param stats_df FDR stats summary tibble from [run_pancan_survival()].
#' @param master Scored master clinical tibble from [run_pancan_survival()].
#' @param cfg Config from [dtp_config()].
#' @param composite_defs Composite definitions table.
#' @return A tibble with rescaled per-1-SD Cox survival results.
#' @export
build_pancan_table <- function(stats_df, master, cfg = dtp_config(), composite_defs = NULL) {
  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  def_row <- if (!is.null(composite_defs)) composite_defs[composite_defs$Is_Default == TRUE, , drop = FALSE] else NULL
  pos_sig <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Positive_Signature[1] else "Up"
  neg_sig <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Negative_Signature[1] else "Down"
  comp_id <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Composite_ID[1] else "Composite"

  up_col  <- paste0(pos_sig, suffix)
  dn_col  <- paste0(neg_sig, suffix)
  cmp_col <- paste0(comp_id, suffix)

  corr_up  <- paste0(up_col, "_Corrected")
  corr_dn  <- paste0(dn_col, "_Corrected")
  corr_cmp <- paste0(cmp_col, "_Corrected")

  raw_cols  <- c(up_col, dn_col, cmp_col)
  corr_cols <- c(corr_up, corr_dn, corr_cmp)

  if (!all(corr_cols %in% colnames(master))) {
    sm <- t(as.matrix(sapply(master[, raw_cols, drop = FALSE], .pancan_num)))
    corrected <- limma::removeBatchEffect(sm, batch = master$Project_ID)
    corr_df <- as.data.frame(t(corrected))
    colnames(corr_df) <- corr_cols
    master <- dplyr::bind_cols(master, corr_df)
  }

  rd3 <- function(x) round(x, 3)
  rd4 <- function(x) round(x, 4)

  # Per-cohort rows (all three scores, both endpoints), per-SD rescaled.
  per <- stats_df %>%
    dplyr::filter(Family == "PerCohort", Test == "Cox", Metric %in% c("OS", "RFS"), Score %in% raw_cols) %>%
    dplyr::mutate(Is_Testable = as.logical(Is_Testable))
  sd_vec <- mapply(function(pr, sc, mt) .pancan_sd(master, pr, sc, mt),
                   per$Project, per$Score, per$Metric)
  per <- .pancan_rescale(per, sd_vec)

  # Aggregate rows (batch-corrected), per-SD rescaled by the corrected pooled SD.
  agg <- stats_df %>%
    dplyr::filter(Family == "PanCancer", Test == "Cox", Metric %in% c("OS", "RFS"), Score %in% corr_cols) %>%
    dplyr::mutate(Is_Testable = as.logical(Is_Testable))
  agg_sd_vec <- mapply(function(sc, mt) .pancan_agg_sd(master, sc, mt), agg$Score, agg$Metric)
  agg <- .pancan_rescale(agg, agg_sd_vec)

  mk_rows <- function(df, project_label) df %>%
    dplyr::transmute(
      Project = if (is.null(project_label)) Project else project_label,
      Score,
      Endpoint = Metric,
      N = as.integer(.pancan_num(N)),
      N_events = as.integer(.pancan_num(N_events)),
      Is_Testable,
      HR_perSD = rd3(hr_sd),
      HR_lower = rd3(lo_sd),
      HR_upper = rd3(hi_sd),
      C_index = rd3(.pancan_num(C_index)),
      Cox_FDR_P = rd4(.pancan_num(FDR_P)),
      Is_Significant = as.logical(Is_Significant)
    )

  tbl <- dplyr::bind_rows(
    mk_rows(per, NULL) %>%
      dplyr::arrange(Project, factor(Score, raw_cols), factor(Endpoint, c("OS", "RFS"))),
    mk_rows(agg, "Pan-Cancer (batch-corrected)") %>%
      dplyr::arrange(factor(Score, corr_cols), factor(Endpoint, c("OS", "RFS")))
  )
  tibble::as_tibble(tbl)
}

# ==============================================================================
# FIGURE 7 - POOLED PAN-CANCER VIOLINS
# ==============================================================================

#' Build Figure 7: Pooled pan-cancer DTP score violins by 3-year outcome
#'
#' Displays batch-corrected DTP scores (Up, Down, Composite) split by 3-year
#' outcome groups (No event vs Event) across OS and RFS, annotated with
#' Pan-Cancer Wilcoxon rank-biserial r and FDR.
#'
#' @param stats_df FDR stats summary tibble from [run_pancan_survival()].
#' @param master Scored master clinical tibble from [run_pancan_survival()].
#' @param cfg Config from [dtp_config()].
#' @param composite_defs Composite definitions table.
#' @return A list with elements `plot` (ggplot object) and `caption` (character string).
#' @export
build_pancan_violins <- function(stats_df, master, cfg = dtp_config(), composite_defs = NULL) {
  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  def_row <- if (!is.null(composite_defs)) composite_defs[composite_defs$Is_Default == TRUE, , drop = FALSE] else NULL
  pos_sig <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Positive_Signature[1] else "Up"
  neg_sig <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Negative_Signature[1] else "Down"
  comp_id <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Composite_ID[1] else "Composite"

  up_col  <- paste0(pos_sig, suffix)
  dn_col  <- paste0(neg_sig, suffix)
  cmp_col <- paste0(comp_id, suffix)

  corr_up  <- paste0(up_col, "_Corrected")
  corr_dn  <- paste0(dn_col, "_Corrected")
  corr_cmp <- paste0(cmp_col, "_Corrected")

  raw_cols  <- c(up_col, dn_col, cmp_col)
  corr_cols <- c(corr_up, corr_dn, corr_cmp)

  if (!all(corr_cols %in% colnames(master))) {
    sm <- t(as.matrix(sapply(master[, raw_cols, drop = FALSE], .pancan_num)))
    corrected <- limma::removeBatchEffect(sm, batch = master$Project_ID)
    corr_df <- as.data.frame(t(corrected))
    colnames(corr_df) <- corr_cols
    master <- dplyr::bind_cols(master, corr_df)
  }

  SCORES <- setNames(c("DTP Up", "DTP Down", "DTP Composite"), corr_cols)
  SCORE_LEVELS <- c("DTP Up", "DTP Down", "DTP Composite")
  PAL <- c("No event" = "#4DBBD5", "Event" = "#E64B35")

  n_coh <- length(unique(master$Project_ID))

  mk_long <- function(score_col, score_lab) {
    dplyr::bind_rows(
      data.frame(ScoreLab = score_lab, Endpoint = "OS", value = .pancan_num(master[[score_col]]),
                 group = c(Dead_3yr = "Event", Alive_3yr = "No event")[as.character(master$Surv_3yr)],
                 stringsAsFactors = FALSE),
      data.frame(ScoreLab = score_lab, Endpoint = "RFS", value = .pancan_num(master[[score_col]]),
                 group = c(Recurred = "Event", `Recurrence-Free` = "No event")[as.character(master$Recurrence_3yr)],
                 stringsAsFactors = FALSE)
    )
  }
  long <- dplyr::bind_rows(Map(mk_long, names(SCORES), unname(SCORES))) %>%
    dplyr::filter(!is.na(group), !is.na(value)) %>%
    dplyr::mutate(
      ScoreLab = factor(ScoreLab, levels = SCORE_LEVELS),
      Endpoint = factor(Endpoint, levels = c("OS", "RFS")),
      group    = factor(group, levels = c("No event", "Event"))
    )

  w <- stats_df %>%
    dplyr::filter(Family == "PanCancer", Test == "Wilcoxon", Metric %in% c("OS", "RFS"), Score %in% corr_cols) %>%
    dplyr::mutate(
      Effect_r = .pancan_num(Effect_r),
      FDR_P = .pancan_num(FDR_P),
      Raw_P = .pancan_num(Raw_P),
      ScoreLab = factor(SCORES[Score], levels = SCORE_LEVELS),
      Endpoint = factor(Metric, levels = c("OS", "RFS")),
      lab_full = sprintf("r = %+.3f\n%s  FDR %s%s", Effect_r,
                         ifelse(is.na(Raw_P), "", ifelse(Raw_P < 0.001, "p<.001", sprintf("p=%.3f", Raw_P))),
                         ifelse(FDR_P < 0.001, "< 0.001", sprintf("= %.3f", FDR_P)),
                         .sig_stars(FDR_P, ""))
    )

  nlab <- long %>%
    dplyr::group_by(ScoreLab, Endpoint, group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  p <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_violin(alpha = 0.85, trim = FALSE, scale = "width", linewidth = 0.3, colour = "grey30") +
    ggplot2::geom_boxplot(width = 0.16, fill = "white", alpha = 0.9, outlier.shape = NA, linewidth = 0.3) +
    ggplot2::geom_text(data = w, inherit.aes = FALSE, ggplot2::aes(x = 1.5, y = Inf, label = lab_full),
                       vjust = 1.15, size = 2.7, colour = "grey20", lineheight = 0.9) +
    ggplot2::geom_text(data = nlab, inherit.aes = FALSE, ggplot2::aes(x = group, y = -Inf, label = paste0("n=", n)),
                       vjust = -0.6, size = 2.3, colour = "grey40") +
    ggplot2::facet_grid(ScoreLab ~ Endpoint, scales = "free_y",
                        labeller = ggplot2::labeller(Endpoint = .PANCAN_ENDPT_LAB)) +
    ggplot2::scale_fill_manual(values = PAL, name = "3-year outcome",
                               labels = c("No event" = "No event", "Event" = "Event (died / recurred)")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::labs(
      title = "Pooled pan-cancer DTP score by 3-year outcome",
      subtitle = paste0(
        sprintf("Batch-corrected score (cancer type as batch), all %d cohorts pooled. ", n_coh),
        "Wilcoxon rank-biserial r, raw p and BH-FDR q (Pan-Cancer family).\n",
        "Negative r = higher score in the event group. * FDR<0.05, ** <0.01, *** <0.001."
      ),
      caption = "Negative rank-biserial r = higher batch-corrected score in the event (died / recurred) group.",
      x = NULL, y = "Batch-corrected ssGSEA score"
    ) +
    .base_theme +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 8.5, colour = "grey30"),
      panel.spacing = grid::unit(7, "pt"),
      panel.grid.major.y = ggplot2::element_line(colour = "grey93", linewidth = 0.3)
    )

  wc <- stats_df[stats_df$Family == "PanCancer" & stats_df$Test == "Wilcoxon" &
                 stats_df$Metric %in% c("OS", "RFS") & stats_df$Score %in% corr_cols, , drop = FALSE]
  wc_txt <- if (!nrow(wc)) "no Pan-Cancer Wilcoxon rows in this run" else {
    lab <- SCORES[wc$Score]
    lab <- ifelse(is.na(lab), wc$Score, sub("^DTP ", "", lab))
    f   <- .pancan_num(wc$FDR_P)
    paste(sprintf("%s %s r = %+.3f (FDR %s)", lab, wc$Metric, .pancan_num(wc$Effect_r),
                  ifelse(is.na(f), "n/a",
                    ifelse(f < 0.001, "< 0.001",
                      paste0(sprintf("%.3f", f), ifelse(f >= 0.05, ", n.s.", ""))))), collapse = ", ")
  }

  cap7 <- c(
    "Figure 7. Pooled pan-cancer DTP-score association by 3-year outcome (Section 3.5).",
    "",
    "Batch-corrected ssGSEA scores (limma::removeBatchEffect with cancer type as batch, recomputed from the master",
    sprintf("exactly as the aggregate) for all %d solid-tumour TCGA cohorts pooled, shown as violins split by 3-year", n_coh),
    "outcome group -- no event vs event (died for OS; recurred for RFS). Six panels: the DTP Up, Down and Composite",
    "scores (rows) by overall survival and recurrence (columns). Each panel is annotated with the Pan-Cancer",
    "Wilcoxon rank-biserial effect size r and its Benjamini-Hochberg FDR; a NEGATIVE r means a higher score in the",
    paste0("event group. In this run: ", wc_txt, "."),
    "Violin widths are scaled equal so the unequal group sizes (annotated as n) are comparable."
  )

  list(plot = p, caption = paste(cap7, collapse = "\n"))
}

# ==============================================================================
# FIGURE 6 COMPOSITE ASSEMBLY
# ==============================================================================

#' Assemble Figure 6: Pan-Cancer Survival Forest Composite
#'
#' Arranges Up, Composite, and Down per-cohort forest panels in a "AA\nBC" layout
#' using [patchwork::wrap_elements()].
#'
#' @param plot_up ggplot object from [build_pancan_forest()] for Up score.
#' @param plot_composite ggplot object from [build_pancan_forest()] for Composite score.
#' @param plot_down ggplot object from [build_pancan_forest()] for Down score.
#' @return A patchwork object.
#' @export
build_fig6_pancan_survival_composite <- function(plot_up, plot_composite, plot_down) {
  if (is.null(plot_up) || is.null(plot_composite) || is.null(plot_down)) {
    stop("plot_up, plot_composite, and plot_down must all be non-NULL.")
  }
  nl <- ggplot2::theme(legend.position = "none")
  design <- "AA\nBC"
  p <- patchwork::wrap_elements(plot_up) +
    patchwork::wrap_elements(plot_composite + nl) +
    patchwork::wrap_elements(plot_down + nl) +
    patchwork::plot_layout(design = design, heights = c(1, 1)) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))
  p
}

# ==============================================================================
# ORCHESTRATOR
# ==============================================================================

#' Orchestrate generation and saving of all 6 pan-cancer composite figures and tables
#'
#' Evaluates all 6 pan-cancer figure builders using in-memory results from
#' [run_pancan_survival()], saves PDF and PNG copies into `out_dir`, writes
#' caption sidecars, and returns the rescaled Table 3.5 tibble.
#'
#' @param pancan_result Result list from [run_pancan_survival()].
#' @param cfg Config from [dtp_config()].
#' @param out_dir Output directory path for figures.
#' @return A named list containing:
#'   \item{pancan_table_tbl}{Table 3.5 pan-cancer survival summary tibble.}
#'   \item{fig6_caption}{Figure 6 caption character string.}
#'   \item{fig7_caption}{Figure 7 caption character string.}
#'   \item{status}{Named logical vector of execution statuses.}
#' @export
build_pancan_composites <- function(pancan_result, cfg = dtp_config(),
                                    out_dir = file.path("output", "figures")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("Building pan-cancer composite figures -> ", out_dir)

  stats_df <- pancan_result$stats_df
  master   <- pancan_result$master
  stopifnot(!is.null(stats_df), !is.null(master))

  composite_defs <- tryCatch(
    load_composite_defs(cfg),
    error = function(e) {
      cpath <- system.file("signatures/composite_defs.csv", package = "dtpsig")
      if (nzchar(cpath) && file.exists(cpath)) {
        readr::read_csv(cpath, show_col_types = FALSE, progress = FALSE)
      } else if (file.exists("inst/signatures/composite_defs.csv")) {
        readr::read_csv("inst/signatures/composite_defs.csv", show_col_types = FALSE, progress = FALSE)
      } else if (file.exists("../../inst/signatures/composite_defs.csv")) {
        readr::read_csv("../../inst/signatures/composite_defs.csv", show_col_types = FALSE, progress = FALSE)
      } else {
        NULL
      }
    }
  )
  suffix  <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  def_row <- if (!is.null(composite_defs)) composite_defs[composite_defs$Is_Default == TRUE, , drop = FALSE] else NULL
  pos_sig <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Positive_Signature[1] else "Up"
  neg_sig <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Negative_Signature[1] else "Down"
  comp_id <- if (!is.null(def_row) && nrow(def_row) >= 1) def_row$Composite_ID[1] else "Composite"

  up_col  <- paste0(pos_sig, suffix)
  dn_col  <- paste0(neg_sig, suffix)
  cmp_col <- paste0(comp_id, suffix)
  order_score <- up_col

  panels  <- new.env(parent = emptyenv())
  results <- list()
  status  <- list()

  # Panel A: Up forest
  status[["A (Up forest)"]] <- tryCatch({
    pA <- build_pancan_forest(stats_df, master, score = up_col, score_lab = "DTP Up", order_score = order_score)
    panels$A <- pA
    save_fig(pA, "Fig6A_pancan_forest_Up", out_dir, 9.0, 8.0)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[pancan composites] group A FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Panel B: Composite forest
  status[["B (Composite forest)"]] <- tryCatch({
    pB <- build_pancan_forest(stats_df, master, score = cmp_col, score_lab = "DTP Composite", order_score = order_score)
    panels$B <- pB
    save_fig(pB, "Fig6B_pancan_forest_Composite", out_dir, 9.0, 8.0)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[pancan composites] group B FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Panel C: Down forest
  status[["C (Down forest)"]] <- tryCatch({
    pC <- build_pancan_forest(stats_df, master, score = dn_col, score_lab = "DTP Down", order_score = order_score)
    panels$C <- pC
    save_fig(pC, "Fig6C_pancan_forest_Down", out_dir, 9.0, 8.0)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[pancan composites] group C FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Fig 6 Composite
  status[["composite"]] <- tryCatch({
    if (is.null(panels$A) || is.null(panels$B) || is.null(panels$C)) {
      stop("panel A/B/C did not all build; composite skipped")
    }
    p_comp <- build_fig6_pancan_survival_composite(panels$A, panels$B, panels$C)
    save_fig(p_comp, "Fig6_pancan_survival_composite", out_dir, 15.0, 16.5)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[pancan composites] group composite FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Aggregate forest
  status[["aggregate (appendix)"]] <- tryCatch({
    p_agg <- build_pancan_aggregate(stats_df, master, cfg = cfg, composite_defs = composite_defs)
    save_fig(p_agg, "Fig_Pancan_aggregate_forest", out_dir, 8.5, 3.2)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[pancan composites] group aggregate FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Violins (Fig 7)
  status[["violins (Fig7)"]] <- tryCatch({
    res_v <- build_pancan_violins(stats_df, master, cfg = cfg, composite_defs = composite_defs)
    results$fig7_caption <- res_v$caption
    save_fig(res_v$plot, "Fig7_pancan_violins", out_dir, 8.5, 9.0)
    writeLines(res_v$caption, file.path(out_dir, "Fig7_pancan_violins_caption.txt"))
    TRUE
  }, error = function(e) {
    msg <- sprintf("[pancan composites] group violins FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Table 3.5
  status[["table"]] <- tryCatch({
    tbl <- build_pancan_table(stats_df, master, cfg = cfg, composite_defs = composite_defs)
    results$pancan_table_tbl <- tbl
    TRUE
  }, error = function(e) {
    msg <- sprintf("[pancan composites] group table FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Run-specific numbers for Fig 6 caption
  .tv <- function(x) { b <- suppressWarnings(as.logical(x)); !is.na(b) & b }
  .cohort_endpoints <- function(sc, want_sig = TRUE) {
    r <- stats_df[stats_df$Family == "PerCohort" & stats_df$Test == "Cox" & stats_df$Score == sc &
                  stats_df$Metric %in% c("OS", "RFS"), , drop = FALSE]
    r <- if (want_sig) r[.tv(r$Is_Significant), , drop = FALSE]
         else          r[!.tv(r$Is_Testable), , drop = FALSE]
    if (!nrow(r)) return(character(0))
    r <- r[order(r$Project, r$Metric), , drop = FALSE]
    paste(r$Project, r$Metric)
  }
  .listify <- function(x, none = "none") if (!length(x)) none else paste(x, collapse = ", ")
  nt_up   <- .cohort_endpoints(up_col, want_sig = FALSE)
  sig_up  <- .cohort_endpoints(up_col)
  sig_cmp <- .cohort_endpoints(cmp_col)
  sig_dn  <- .cohort_endpoints(dn_col)
  n_coh   <- length(unique(master$Project_ID))

  cap6 <- c(
    "Figure 6. Per-cohort DTP-score survival association across the pan-cancer TCGA cohort (Section 3.5).",
    "Per-cohort forests, one per DTP score, on a shared cohort order: DTP Up on top, DTP Composite and DTP Down below.",
    "",
    "(A) DTP Up, (B) DTP Composite, (C) DTP Down. Each panel is a forest of the univariable continuous-score Cox",
    sprintf("hazard ratio per 1 standard deviation of that score within each of %d solid-tumour TCGA cohorts, with 95%%", n_coh),
    "confidence interval on a log scale (truncated to [0.1, 10]), faceted by endpoint: overall survival (OS) and",
    "recurrence (RFS, derived from the TCGA-CDR disease-free interval). Red = HR > 1 (higher score, higher",
    "hazard), blue = HR < 1; an asterisk marks FDR < 0.05; a grey \"n/t\" marks a cohort that failed the Cox gate",
    sprintf("(n >= 10 and >= 5 events; in this run: %s). All three panels share the SAME row order -- ascending DTP Up",
            .listify(nt_up, "no cohort failed the gate")),
    "overall-survival per-SD HR -- so the tissues line up; TCGA-COAD (the colorectal tissue) is highlighted in each.",
    sprintf("Up is FDR-significant in %d cohort-endpoint(s) (%s);", length(sig_up), .listify(sig_up)),
    sprintf("Composite in %d (%s);", length(sig_cmp), .listify(sig_cmp)),
    sprintf("Down in %d (%s).", length(sig_dn), .listify(sig_dn)),
    "",
    "Hazard ratios are per 1 SD; the raw Cox HRs are per unit of the (log2) ssGSEA score and are not shown. FDR is",
    "Benjamini-Hochberg within each family (per-cohort vs. batch-corrected aggregate) and test, pooling OS and RFS.",
    "The pooled association is shown as violins in Figure 7; the batch-corrected aggregate Cox forest is an appendix",
    "figure (Fig_Pancan_aggregate_forest). Full per-cohort and aggregate statistics are in Appendix Table A4",
    "(Table_3_5_pancan_survival.csv).")

  results$fig6_caption <- paste(cap6, collapse = "\n")
  writeLines(cap6, file.path(out_dir, "Fig6_pancan_survival_caption.txt"))

  n_ok  <- sum(unlist(status))
  n_tot <- length(status)
  summ  <- sprintf("Pan-cancer composites: %d/%d groups OK%s", n_ok, n_tot,
                   if (n_ok < n_tot) paste0(" -- FAILED: ", paste(names(status)[!unlist(status)], collapse = ", ")) else "")
  message(summ)
  if (n_ok < n_tot) warning(summ, call. = FALSE)
  results$status <- status
  message("Pan-cancer composite figures done -> ", out_dir)

  results
}
