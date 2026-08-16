# ==============================================================================
# R/composites/crc_composites.R - CRC Survival Presentation Layer (7 Composites)
#
# Generates publication-ready figures and summary tables directly from in-memory
# analysis tibbles (no CSV reads inside builders).
# ==============================================================================

#' Build Figure 1: 3-year outcome composite (effect-size matrix + violins)
#'
#' @param stats_df Survival statistics summary tibble from [run_crc_survival()].
#' @param gse_clinical Annotated clinical tibble for GSE39582.
#' @param tcga_clinical Annotated clinical tibble for TCGA-COAD.
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param cfg Config from [dtp_config()].
#' @return A patchwork object combining the effect-size matrix and violin grid.
#' @export
build_fig1_outcome_composite <- function(stats_df, gse_clinical, tcga_clinical,
                                         panel_tbl, composite_defs, cfg) {
  PAL <- c(Poor = "#E64B35", Good = "#4DBBD5")

  CONTRASTS <- tibble::tribble(
    ~key,          ~dataset,             ~cohort,             ~endpoint, ~stats_ds,   ~stats_cohort,        ~clin,   ~group_col,       ~poor,       ~good,
    "M_whole_OS",  "Marisa (GSE39582)", "Whole cohort",     "OS",      "GSE39582",  "GSE_All_Patients",  "gse",   "Surv_3yr",       "Dead_3yr", "Alive_3yr",
    "M_whole_RFS", "Marisa (GSE39582)", "Whole cohort",     "RFS",     "GSE39582",  "GSE_All_Patients",  "gse",   "Recurrence_3yr", "Recurred", "Recurrence-Free",
    "M_trt_OS",    "Marisa (GSE39582)", "Adjuvant-treated", "OS",      "GSE39582",  "GSE_Treated",       "gse",   "Surv_3yr",       "Dead_3yr", "Alive_3yr",
    "M_trt_RFS",   "Marisa (GSE39582)", "Adjuvant-treated", "RFS",     "GSE39582",  "GSE_Treated",       "gse",   "Recurrence_3yr", "Recurred", "Recurrence-Free",
    "T_whole_OS",  "TCGA-COAD",         "Whole cohort",     "OS",      "TCGA-COAD", "TCGA_All_Patients", "tcga",  "Surv_3yr",       "Dead_3yr", "Alive_3yr",
    "T_whole_RFS", "TCGA-COAD",         "Whole cohort",     "RFS",     "TCGA-COAD", "TCGA_All_Patients", "tcga",  "Recurrence_3yr", "Recurred", "Recurrence-Free"
  )
  CONTRASTS$key     <- factor(CONTRASTS$key, levels = CONTRASTS$key)
  CONTRASTS$dataset <- factor(CONTRASTS$dataset, levels = unique(CONTRASTS$dataset))
  CONTRASTS$col_lab <- paste0(CONTRASTS$cohort, "\n(", ENDPOINT_LABEL[CONTRASTS$endpoint], ")")
  col_lab_map <- setNames(CONTRASTS$col_lab, CONTRASTS$key)

  clin <- list(gse = gse_clinical, tcga = tcga_clinical)

  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  score_cols <- unique(c(
    grep(paste0(suffix, "$"), names(gse_clinical), value = TRUE),
    grep(paste0(suffix, "$"), names(tcga_clinical), value = TRUE)
  ))
  MOD <- .crc_modules(score_cols, panel_tbl, composite_defs, cfg)

  cohort_subset <- function(row) {
    d <- clin[[row$clin]]
    if (row$cohort == "Adjuvant-treated") d <- d[which(d$Chemo_adj == "Y"), , drop = FALSE]
    d
  }

  stat_lookup <- function(row, score) {
    s <- stats_df[stats_df$Dataset == row$stats_ds & stats_df$Project == row$stats_cohort &
                  stats_df$Test == "Wilcoxon" & stats_df$Metric == row$endpoint &
                  stats_df$Score == score, ]
    if (!nrow(s)) return(list(fdr = NA, r = NA, raw = NA, n = NA, sig = FALSE))
    list(fdr = s$FDR_P[1], r = s$Effect_r[1], raw = s$Raw_P[1], n = s$N[1],
         sig = isTRUE(as.logical(s$Is_Significant[1])))
  }

  long_rows <- list(); cell_rows <- list()
  for (i in seq_len(nrow(CONTRASTS))) {
    row <- CONTRASTS[i, ]
    d   <- cohort_subset(row)
    g   <- d[[row$group_col]]
    keep <- g %in% c(row$poor, row$good)
    d <- d[keep, ]; g <- g[keep]
    grp <- factor(ifelse(g == row$poor, "Poor", "Good"), levels = c("Poor", "Good"))
    for (sc in MOD$cols) {
      if (!sc %in% names(d)) next
      v  <- suppressWarnings(as.numeric(d[[sc]]))
      ok <- !is.na(v)
      long_rows[[length(long_rows) + 1]] <- tibble::tibble(
        key = row$key, module = MOD$labels[[sc]], score_col = sc, group = grp[ok], value = v[ok]
      )
      st <- stat_lookup(row, sc)
      cell_rows[[length(cell_rows) + 1]] <- tibble::tibble(
        key = row$key, module = MOD$labels[[sc]], score_col = sc,
        fdr = st$fdr, raw = st$raw, signed_r = -st$r, sig = st$sig, stars = .sig_stars(st$fdr, "ns")
      )
    }
  }
  long  <- dplyr::bind_rows(long_rows) %>%
    dplyr::mutate(module = factor(module, levels = MOD$levels), key = factor(key, levels = levels(CONTRASTS$key)))
  cells <- dplyr::bind_rows(cell_rows) %>%
    dplyr::mutate(module = factor(module, levels = MOD$levels), key = factor(key, levels = levels(CONTRASTS$key))) %>%
    dplyr::left_join(dplyr::select(CONTRASTS, key, dataset), by = "key")

  # Panel A: Effect matrix
  cells$lab_clean <- ifelse(cells$sig, cells$stars, "")
  lim <- max(abs(cells$signed_r), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1.0

  p_matrix <- ggplot2::ggplot(cells, ggplot2::aes(x = key, y = forcats::fct_rev(module), fill = signed_r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = lab_clean), size = 4.2, fontface = "bold",
      colour = "grey15", vjust = 0.5, lineheight = 0.82
    ) +
    ggplot2::facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
      limits = c(-lim, lim), name = "Rank-biserial r",
      guide = ggplot2::guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)
    ) +
    ggplot2::scale_x_discrete(labels = col_lab_map) +
    ggplot2::labs(
      title = "DTP-module scores vs 3-year outcome across cohorts",
      subtitle = "Colour = effect size (red: higher in poor-outcome group).  * marks FDR-significant Wilcoxon (BH, p<0.05).",
      x = NULL, y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, lineheight = 0.85),
      panel.spacing.x = grid::unit(6, "pt")
    )

  # Panel B: Violins
  curated_levels <- MOD$levels
  vl <- long  %>% dplyr::filter(module %in% curated_levels) %>% dplyr::mutate(module = factor(module, levels = curated_levels))
  vc <- cells %>% dplyr::filter(module %in% curated_levels) %>% dplyr::mutate(module = factor(module, levels = curated_levels))
  yr <- vl %>%
    dplyr::group_by(module, key) %>%
    dplyr::summarise(ymin = min(value), ymax = max(value), .groups = "drop") %>%
    dplyr::mutate(rng = ymax - ymin)
  vc <- vc %>% dplyr::left_join(yr, by = c("module", "key"))
  vc$lab_clean_v <- ifelse(vc$stars == "", "ns", vc$stars)

  nlab <- vl %>%
    dplyr::group_by(module, key, group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(dplyr::select(yr, module, key, ymin, rng), by = c("module", "key"))

  build_violin_block <- function(keys, title, show_y) {
    dl <- dplyr::filter(vl, key %in% keys)
    dc <- dplyr::filter(vc, key %in% keys)
    dn <- dplyr::filter(nlab, key %in% keys)

    ggplot2::ggplot(dl, ggplot2::aes(x = group, y = value, fill = group)) +
      ggplot2::geom_violin(alpha = 0.85, trim = FALSE, linewidth = 0.3, colour = "grey30") +
      ggplot2::geom_boxplot(width = 0.16, fill = "white", alpha = 0.9, outlier.shape = NA, linewidth = 0.3) +
      ggplot2::geom_text(
        data = dc, inherit.aes = FALSE,
        ggplot2::aes(x = 1.5, y = ymax + 0.12 * rng, label = lab_clean_v),
        fontface = "bold", size = 3.2, lineheight = 0.82
      ) +
      ggplot2::geom_text(
        data = dn, inherit.aes = FALSE,
        ggplot2::aes(x = group, y = ymin - 0.10 * rng, label = paste0("n=", n)),
        size = 2.4, colour = "grey35"
      ) +
      ggplot2::facet_grid(
        module ~ key, scales = "free_y", switch = "y",
        labeller = ggplot2::labeller(key = col_lab_map, module = scales::label_wrap(14))
      ) +
      ggplot2::scale_fill_manual(
        values = PAL, name = "3-year outcome",
        labels = c(Poor = "Deceased / recurred", Good = "Alive / recurrence-free")
      ) +
      ggplot2::scale_x_discrete(labels = c(Poor = "Poor", Good = "Good")) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.16)) +
      ggplot2::labs(title = title, x = NULL, y = if (show_y) "ssGSEA score" else NULL) +
      .base_theme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
        strip.placement = "outside",
        strip.text.y.left = if (show_y) ggplot2::element_text(angle = 0, face = "bold", size = 8.5) else ggplot2::element_blank(),
        strip.text.x = ggplot2::element_text(size = 8, lineheight = 0.85),
        axis.text.x = ggplot2::element_text(size = 8),
        panel.spacing = grid::unit(4, "pt")
      )
  }

  block_m <- build_violin_block(c("M_whole_OS", "M_whole_RFS", "M_trt_OS", "M_trt_RFS"), "Marisa (GSE39582)", TRUE)
  block_t <- build_violin_block(c("T_whole_OS", "T_whole_RFS"), "TCGA-COAD", FALSE)

  p_violins <- block_m + block_t +
    patchwork::plot_layout(widths = c(4, 2), guides = "collect") +
    patchwork::plot_annotation(
      title = "Signature-score distributions by 3-year outcome",
      subtitle = "Full cohort, true n shown; y-axis free within each dataset (ssGSEA scale differs GSE vs TCGA).  * = FDR-significant; ns = not.",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 10, colour = "grey30", hjust = 0)
      )
    ) &
    ggplot2::theme(legend.position = "bottom")

  n_rows <- length(curated_levels)
  viol_h <- 1.5 * n_rows + 1.4

  patchwork::wrap_elements(p_matrix) / patchwork::wrap_elements(p_violins) +
    patchwork::plot_layout(heights = c(4.6, viol_h - 0.4)) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))
}

#' Build Figure 2: Kaplan-Meier survival composite (Cox HR matrix + KM curves)
#'
#' @param stats_df Survival statistics summary tibble from [run_crc_survival()].
#' @param gse_clinical Annotated clinical tibble for GSE39582.
#' @param tcga_clinical Annotated clinical tibble for TCGA-COAD.
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param cfg Config from [dtp_config()].
#' @return A patchwork object combining the Cox HR/SD matrix and KM curve grid.
#' @export
build_fig2_km_composite <- function(stats_df, gse_clinical, tcga_clinical,
                                    panel_tbl, composite_defs, cfg) {
  PAL <- c(Low = "#4DBBD5", High = "#E64B35")
  X_CAP <- if (!is.null(cfg$os_cutpoint)) cfg$os_cutpoint else 36

  CONTRASTS <- tibble::tribble(
    ~key,          ~dataset,             ~cohort,             ~metric, ~stats_ds,   ~stats_cohort,        ~clin,
    "M_whole_OS",  "Marisa (GSE39582)", "Whole cohort",      "OS",    "GSE39582",  "GSE_All_Patients",  "gse",
    "M_whole_RFS", "Marisa (GSE39582)", "Whole cohort",      "RFS",   "GSE39582",  "GSE_All_Patients",  "gse",
    "M_trt_OS",    "Marisa (GSE39582)", "Adjuvant-treated",  "OS",    "GSE39582",  "GSE_Treated",       "gse",
    "M_trt_RFS",   "Marisa (GSE39582)", "Adjuvant-treated",  "RFS",   "GSE39582",  "GSE_Treated",       "gse",
    "T_whole_OS",  "TCGA-COAD",         "Whole cohort",      "OS",    "TCGA-COAD", "TCGA_All_Patients", "tcga",
    "T_whole_RFS", "TCGA-COAD",         "Whole cohort",      "RFS",   "TCGA-COAD", "TCGA_All_Patients", "tcga"
  )
  CONTRASTS$key     <- factor(CONTRASTS$key, levels = CONTRASTS$key)
  CONTRASTS$dataset <- factor(CONTRASTS$dataset, levels = unique(CONTRASTS$dataset))
  CONTRASTS$col_lab <- paste0(CONTRASTS$cohort, "\n(", ENDPOINT_LABEL[CONTRASTS$metric], ")")
  col_lab_map <- setNames(CONTRASTS$col_lab, CONTRASTS$key)

  clin <- list(gse = gse_clinical, tcga = tcga_clinical)

  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  score_cols <- unique(c(
    grep(paste0(suffix, "$"), names(gse_clinical), value = TRUE),
    grep(paste0(suffix, "$"), names(tcga_clinical), value = TRUE)
  ))
  MOD <- .crc_modules(score_cols, panel_tbl, composite_defs, cfg)

  cohort_subset <- function(row) {
    d <- clin[[row$clin]]
    if (row$cohort == "Adjuvant-treated") d <- d[which(d$Chemo_adj == "Y"), , drop = FALSE]
    d
  }

  cox_stat <- function(row, score) {
    s <- stats_df[stats_df$Dataset == row$stats_ds & stats_df$Project == row$stats_cohort &
                  stats_df$Test == "Cox" & stats_df$Metric == row$metric & stats_df$Score == score, ]
    if (!nrow(s) || !isTRUE(as.logical(s$Is_Testable[1]))) {
      return(list(hr_unit = NA_real_, fdr = NA_real_, raw = NA_real_, sig = FALSE, testable = FALSE))
    }
    list(hr_unit = as.numeric(s$HR[1]), fdr = as.numeric(s$FDR_P[1]),
         raw = as.numeric(s$Raw_P[1]), sig = isTRUE(as.logical(s$Is_Significant[1])), testable = TRUE)
  }

  logrank_stat <- function(row, score) {
    s <- stats_df[stats_df$Dataset == row$stats_ds & stats_df$Project == row$stats_cohort &
                  stats_df$Test == "KM" & stats_df$Metric == row$metric & stats_df$Score == score, ]
    if (!nrow(s)) return(list(fdr = NA_real_, raw = NA_real_, sig = FALSE))
    list(fdr = as.numeric(s$FDR_P[1]), raw = as.numeric(s$Raw_P[1]),
         sig = isTRUE(as.logical(s$Is_Significant[1])))
  }

  curve_rows <- list(); cell_rows <- list()
  for (i in seq_len(nrow(CONTRASTS))) {
    row <- CONTRASTS[i, ]
    d0  <- cohort_subset(row)
    mc  <- metric_cols(row$metric)
    tcol <- mc[["t"]]; ecol <- mc[["e"]]

    for (sc in MOD$cols) {
      if (!all(c(sc, tcol, ecol) %in% names(d0))) next
      d <- d0[!is.na(d0[[sc]]) & !is.na(d0[[tcol]]) & !is.na(d0[[ecol]]), ]
      score <- as.numeric(d[[sc]])
      d <- data.frame(
        time = pmin(as.numeric(d[[tcol]]), X_CAP),
        event = as.numeric(d[[ecol]]),
        Group = factor(ifelse(score > stats::median(score), "High", "Low"), levels = c("Low", "High"))
      )
      tb <- table(d$Group)
      min_grp <- if (!is.null(cfg$min_km_group_n)) cfg$min_km_group_n else 5
      min_ev  <- if (!is.null(cfg$min_events)) cfg$min_events else 5
      drawable <- length(tb) == 2 && all(tb >= min_grp) && sum(d$event) >= min_ev

      cx <- cox_stat(row, sc)
      lr <- logrank_stat(row, sc)
      hr_sd <- if (cx$testable && is.finite(cx$hr_unit)) exp(log(cx$hr_unit) * .score_sd(score)) else NA_real_

      cell_rows[[length(cell_rows) + 1]] <- tibble::tibble(
        key = row$key, module = MOD$labels[[sc]], score_col = sc, testable = cx$testable,
        hr = hr_sd, log2hr = log2(hr_sd), fdr = cx$fdr, raw = cx$raw, sig = cx$sig,
        stars = .sig_stars(cx$fdr, "ns"), lr_fdr = lr$fdr, lr_raw = lr$raw, lr_sig = lr$sig
      )

      if (!drawable) next
      sf <- survival::survfit(survival::Surv(time, event) ~ Group, data = d)
      strata <- rep(sub(".*=", "", names(sf$strata)), sf$strata)
      cd <- rbind(
        data.frame(time = 0, surv = 1, ncens = 0, Group = levels(d$Group)),
        data.frame(time = sf$time, surv = sf$surv, ncens = sf$n.censor, Group = strata)
      )
      cd$key <- row$key
      cd$module <- MOD$labels[[sc]]
      cd$small <- nrow(d) < 500
      curve_rows[[length(curve_rows) + 1]] <- cd
    }
  }

  cells <- dplyr::bind_rows(cell_rows) %>%
    dplyr::mutate(module = factor(module, levels = MOD$levels), key = factor(key, levels = levels(CONTRASTS$key))) %>%
    dplyr::left_join(dplyr::select(CONTRASTS, key, dataset), by = "key")

  curves <- dplyr::bind_rows(curve_rows) %>%
    dplyr::mutate(
      module = factor(module, levels = MOD$levels),
      key = factor(key, levels = levels(CONTRASTS$key)),
      Group = factor(Group, levels = c("Low", "High"))
    )

  # Panel A: Cox HR matrix
  cells$lab_clean <- ifelse(cells$sig, cells$stars, "")
  lim <- min(max(abs(cells$log2hr[is.finite(cells$log2hr)]), na.rm = TRUE), 2)
  if (!is.finite(lim) || lim == 0) lim <- 1.0

  p_matrix <- ggplot2::ggplot(
    cells,
    ggplot2::aes(x = key, y = forcats::fct_rev(module), fill = pmax(pmin(log2hr, lim), -lim))
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = lab_clean), size = 4.2, fontface = "bold",
      colour = "grey15", vjust = 0.5, lineheight = 0.82
    ) +
    ggplot2::facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
      limits = c(-lim, lim), name = expression(log[2]~"HR (per 1 SD of score)"),
      guide = ggplot2::guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1), na.value = "grey85"
    ) +
    ggplot2::scale_x_discrete(labels = col_lab_map) +
    ggplot2::labs(
      title = "Survival association of DTP-module scores across cohorts",
      subtitle = "Univariable continuous-score Cox PH, follow-up truncated at 36 mo.  Colour = log2 HR per 1 SD (red = higher score worse).  * = FDR<0.05 (BH).",
      x = NULL, y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, lineheight = 0.85),
      panel.spacing.x = grid::unit(6, "pt")
    )

  # Panel B: KM curves
  km_levels <- MOD$levels
  cv <- curves %>% dplyr::filter(module %in% km_levels) %>% dplyr::mutate(module = factor(module, levels = km_levels))
  drawn <- dplyr::distinct(cv, key, module)
  cc <- cells %>%
    dplyr::filter(module %in% km_levels) %>%
    dplyr::semi_join(drawn, by = c("key", "module")) %>%
    dplyr::mutate(
      module = factor(module, levels = km_levels),
      lab_clean = ifelse(is.na(lr_fdr), "log-rank n/a", paste0("log-rank ", .sig_stars(lr_fdr, "n.s.")))
    )

  build_km_block <- function(keys, title, show_y) {
    dc <- dplyr::filter(cv, key %in% keys)
    da <- dplyr::filter(cc, key %in% keys)

    ggplot2::ggplot(dc, ggplot2::aes(x = time, y = surv, colour = Group)) +
      ggplot2::geom_step(linewidth = 0.6) +
      ggplot2::geom_point(data = subset(dc, ncens > 0 & small), shape = 3, size = 1.1, show.legend = FALSE) +
      ggplot2::geom_text(
        data = da, inherit.aes = FALSE, ggplot2::aes(x = 0.5, y = 0.06, label = lab_clean),
        hjust = 0, vjust = 0, size = 2.5, colour = "grey25", lineheight = 0.9
      ) +
      ggplot2::facet_grid(
        module ~ key, switch = "y",
        labeller = ggplot2::labeller(key = col_lab_map, module = scales::label_wrap(14))
      ) +
      ggplot2::scale_colour_manual(
        values = PAL, name = "ssGSEA score",
        labels = c(Low = "Low (< median)", High = "High (> median)")
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, X_CAP), breaks = c(0, 12, 24, 36),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1), labels = scales::percent_format(accuracy = 1),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        title = title, x = "Months from diagnosis (follow-up truncated at 36)",
        y = if (show_y) "Survival probability" else NULL
      ) +
      .base_theme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
        strip.placement = "outside",
        strip.text.y.left = if (show_y) ggplot2::element_text(angle = 0, face = "bold", size = 8.5) else ggplot2::element_blank(),
        strip.text.x = ggplot2::element_text(size = 8, lineheight = 0.85),
        panel.spacing.x = grid::unit(11, "pt"), panel.spacing.y = grid::unit(5, "pt"),
        panel.grid.major.y = ggplot2::element_line(colour = "grey93", linewidth = 0.3)
      )
  }

  block_m <- build_km_block(c("M_whole_OS", "M_whole_RFS", "M_trt_OS", "M_trt_RFS"), "Marisa (GSE39582)", TRUE)
  block_t <- build_km_block(c("T_whole_OS", "T_whole_RFS"), "TCGA-COAD", FALSE)

  p_km <- block_m + block_t +
    patchwork::plot_layout(widths = c(4, 2), guides = "collect") +
    patchwork::plot_annotation(
      title = "Kaplan-Meier survival by signature score (median High/Low split)",
      subtitle = "Median High/Low split; log-rank BH-FDR significance shown as stars (* <0.05, ** <0.01, *** <0.001; n.s. = not significant). Continuous-Cox HR effect sizes are in the matrix panel. Time runs from diagnosis with follow-up truncated at 36 mo (NOT a 3-yr landmark).",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 10, colour = "grey30", hjust = 0)
      )
    ) &
    ggplot2::theme(legend.position = "bottom")

  n_rows <- length(km_levels)
  km_h <- 2.0 * n_rows + 1.4

  patchwork::wrap_elements(p_matrix) / patchwork::wrap_elements(p_km) +
    patchwork::plot_layout(heights = c(4.6, km_h - 0.4)) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))
}

#' Build Figure 3: Subgroup composite (interaction matrix + per-score forest)
#'
#' @param int_df Interaction statistics tibble from [run_crc_survival()].
#' @param level_df Per-level subgroup HR tibble from [run_crc_survival()].
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param cfg Config from [dtp_config()].
#' @param primary_score Primary score column to plot in forest panel (e.g. "Up_ssGSEA").
#' @param gse_clinical Optional annotated clinical tibble for GSE39582 (used for exact within-subgroup score SD).
#' @param tcga_clinical Optional annotated clinical tibble for TCGA-COAD (used for exact within-subgroup score SD).
#' @return A patchwork object combining the interaction matrix and subgroup forest.
#' @export
build_fig3_subgroup_composite <- function(int_df, level_df, panel_tbl, composite_defs, cfg,
                                          primary_score = "Up_ssGSEA",
                                          gse_clinical = NULL, tcga_clinical = NULL) {
  FOREST_XLIM <- c(0.1, 10)
  MOD_LABEL  <- c(CMS = "CMS subtype", PDS = "PDS subtype", Stage_bin = "Stage", MSI_group = "MSI status")
  MOD_LEVELS <- names(MOD_LABEL)
  LEVEL_ORDER <- c("CMS1", "CMS2", "CMS3", "CMS4", "PDS1", "PDS2", "PDS3", "Early", "Late", "MSS", "MSI")

  COLS <- tibble::tibble(
    ckey    = c("GSE_OS", "GSE_RFS", "GSEt_OS", "GSEt_RFS", "TCGA_OS", "TCGA_RFS"),
    Dataset = c("GSE39582", "GSE39582", "GSE39582 (treated)", "GSE39582 (treated)", "TCGA-COAD", "TCGA-COAD"),
    Metric  = c("OS", "RFS", "OS", "RFS", "OS", "RFS")
  )
  COLS$clab <- paste0(DATASET_SHORT[COLS$Dataset], "\n(", ENDPOINT_LABEL[COLS$Metric], ")")
  col_lab_map <- setNames(COLS$clab, COLS$ckey)
  DS_PREFIX <- c("GSE39582" = "GSE", "GSE39582 (treated)" = "GSEt", "TCGA-COAD" = "TCGA")
  DIR_PAL <- c("Higher hazard (HR>1)" = "#E64B35", "Lower hazard (HR<1)" = "#3B6EA5")

  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  core_score_ids <- c("Up", "Composite", "Down")
  core_scores <- intersect(paste0(core_score_ids, suffix), unique(int_df$Score))
  if (length(core_scores) == 0) core_scores <- unique(int_df$Score)
  MOD <- .crc_modules(core_scores, panel_tbl, composite_defs, cfg)

  # Within-subgroup SD calculation
  if (!is.null(gse_clinical) && !is.null(tcga_clinical)) {
    clin <- list(
      GSE39582 = harmonize_crc_modifiers(gse_clinical, "GSE39582"),
      "TCGA-COAD" = harmonize_crc_modifiers(tcga_clinical, "TCGA-COAD")
    )
    clin[["GSE39582 (treated)"]] <- clin$GSE39582[which(clin$GSE39582$Chemo_adj == "Y"), , drop = FALSE]
    metric_evt <- c(OS = "OS3Y_event", RFS = "RFS3Y_event")
    metric_tim <- c(OS = "OS3Y_delay", RFS = "RFS3Y_delay")

    sd_of <- function(ds, metric, score, modifier, level) {
      d <- clin[[ds]]; ev <- metric_evt[[metric]]; tm <- metric_tim[[metric]]
      if (is.null(d) || !all(c(score, ev, tm, modifier) %in% names(d))) return(1.0)
      d <- d[which(as.character(d[[modifier]]) == level), , drop = FALSE]
      v <- suppressWarnings(as.numeric(d[[score]]))
      ok <- !is.na(v) & !is.na(d[[tm]]) & !is.na(d[[ev]])
      sd_val <- .score_sd(v[ok])
      if (is.na(sd_val) || sd_val <= 0) 1.0 else sd_val
    }
    sd_sc <- mapply(sd_of, level_df$Dataset, level_df$Metric, level_df$Score, level_df$Modifier, level_df$Level)
  } else {
    sd_sc <- rep(1.0, nrow(level_df))
  }

  subg <- level_df %>%
    dplyr::mutate(
      ckey  = paste0(unname(DS_PREFIX[Dataset]), "_", Metric),
      sd_sc = sd_sc,
      hr_sd = exp(log(HR) * sd_sc),
      lo_sd = exp(log(HR_lower) * sd_sc),
      hi_sd = exp(log(HR_upper) * sd_sc),
      Modifier = factor(Modifier, levels = MOD_LEVELS),
      Level    = factor(Level, levels = LEVEL_ORDER),
      ckey     = factor(ckey, levels = COLS$ckey),
      nonconv  = (!is.na(hr_sd) & !is.finite(hr_sd)) |
                 (!is.na(lo_sd) & !is.finite(lo_sd)) |
                 (!is.na(hi_sd) & !is.finite(hi_sd)),
      offscale = !nonconv & !is.na(hr_sd) & (hr_sd < FOREST_XLIM[1] | hr_sd > FOREST_XLIM[2]),
      hr_d     = ifelse(nonconv, NA_real_, scales::squish(hr_sd, FOREST_XLIM)),
      lo_d     = ifelse(nonconv, NA_real_, scales::squish(lo_sd, FOREST_XLIM)),
      hi_d     = ifelse(nonconv, NA_real_, scales::squish(hi_sd, FOREST_XLIM)),
      ne_clean = sprintf("n=%d, e=%d", as.integer(N), as.integer(N_events))
    )

  inter <- int_df %>%
    dplyr::mutate(
      ckey = paste0(unname(DS_PREFIX[Dataset]), "_", Metric),
      Modifier = factor(Modifier, levels = MOD_LEVELS),
      ckey = factor(ckey, levels = COLS$ckey),
      int_lab = ifelse(is.na(FDR_P), "interaction n/a", paste0("interaction ", .sig_stars(FDR_P, "n.s."))),
      int_sig = !is.na(Is_Significant) & as.logical(Is_Significant)
    )

  # Panel A: Interaction Matrix
  mat <- inter %>%
    dplyr::filter(Score %in% MOD$cols) %>%
    dplyr::mutate(
      Score    = factor(MOD$labels[Score], levels = MOD$levels),
      Modifier = forcats::fct_rev(factor(MOD_LABEL[as.character(Modifier)], levels = unname(MOD_LABEL))),
      mlog     = pmin(-log10(FDR_P), 4),
      lab_clean = .sig_stars(FDR_P)
    )

  p_matrix <- ggplot2::ggplot(mat, ggplot2::aes(x = ckey, y = Modifier, fill = mlog)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = lab_clean), fontface = "bold", size = 4.2, vjust = 0.5, lineheight = 0.82) +
    ggplot2::facet_grid(. ~ Score) +
    ggplot2::scale_fill_gradient(
      low = "white", high = "#762A83", limits = c(0, 4),
      name = expression(-log[10]~"interaction FDR"),
      guide = ggplot2::guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)
    ) +
    ggplot2::scale_x_discrete(labels = col_lab_map) +
    ggplot2::labs(
      title = "Effect modification of the DTP signature (score-by-subgroup interaction)",
      subtitle = "Likelihood-ratio test of score*modifier vs score+modifier (Cox, follow-up truncated at 36 mo).  * marks FDR-significant (BH, p<0.05).",
      x = NULL, y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, lineheight = 0.85),
      panel.spacing.x = grid::unit(6, "pt")
    )

  # Panel B: Subgroup Forest (primary score)
  ds <- subg  %>% dplyr::filter(Score == primary_score, !nonconv)
  di <- inter %>% dplyr::filter(Score == primary_score)
  ds$dir <- ifelse(ds$hr_sd > 1, "Higher hazard (HR>1)", "Lower hazard (HR<1)")

  primary_label <- if (primary_score %in% names(MOD$labels)) MOD$labels[[primary_score]] else primary_score

  p_forest <- ggplot2::ggplot(ds, ggplot2::aes(x = hr_d, y = forcats::fct_rev(Level))) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo_d, xmax = hi_d, colour = dir), height = 0.28, linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(colour = dir, shape = offscale), size = 2.1) +
    ggplot2::geom_text(
      inherit.aes = FALSE,
      ggplot2::aes(x = FOREST_XLIM[1], y = forcats::fct_rev(Level), label = ne_clean),
      hjust = 0, vjust = -0.5, size = 1.95, colour = "grey45", lineheight = 0.82
    ) +
    ggplot2::geom_text(
      data = di, inherit.aes = FALSE,
      ggplot2::aes(x = FOREST_XLIM[1], y = Inf, label = int_lab),
      hjust = 0, vjust = 1.25, size = 2.5,
      colour = ifelse(di$int_sig, "#B2182B", "grey35"), lineheight = 0.82
    ) +
    ggplot2::scale_y_discrete() +
    ggplot2::facet_grid(
      Modifier ~ ckey, scales = "free_y", space = "free_y", switch = "y",
      labeller = ggplot2::labeller(ckey = col_lab_map, Modifier = function(x) MOD_LABEL[x])
    ) +
    ggplot2::scale_x_log10(
      limits = FOREST_XLIM, breaks = c(0.1, 0.3, 1, 3, 10),
      labels = c("0.1", "0.3", "1", "3", "10")
    ) +
    ggplot2::scale_colour_manual(values = DIR_PAL, name = "Score effect", drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 21),
      labels = c(`FALSE` = "in range", `TRUE` = "off-scale (squished)"), name = NULL
    ) +
    ggplot2::labs(
      title = paste0("Per-subgroup hazard ratio of ", primary_label, " score"),
      subtitle = "Continuous-score Cox HR per 1 SD within each subgroup (95% CI); x truncated to [0.1, 10]. n / e = subgroup size / events in that cohort.",
      x = "Hazard ratio per 1 SD (log scale)", y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 8.5),
      strip.text.x = ggplot2::element_text(size = 8, lineheight = 0.85),
      axis.text.y = ggplot2::element_text(size = 7.5),
      panel.spacing = grid::unit(4, "pt"),
      panel.grid.major.x = ggplot2::element_line(colour = "grey93", linewidth = 0.3)
    )

  patchwork::wrap_elements(p_matrix) / patchwork::wrap_elements(p_forest) +
    patchwork::plot_layout(heights = c(3.4, 8.2)) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))
}

#' Build Figure 3_3A: Univariable subgroup Cox HR matrix and Table 3.3
#'
#' @param stats_df Survival statistics summary tibble from [run_crc_survival()].
#' @param subtype_stats Kruskal-Wallis stats tibble from [run_crc_survival()].
#' @param gse_clinical Annotated clinical tibble for GSE39582.
#' @param tcga_clinical Annotated clinical tibble for TCGA-COAD.
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param cfg Config from [dtp_config()].
#' @return Named list with `plot` (ggplot object) and `table` (tibble for Table 3.3).
#' @export
build_fig3_3a_subgroup_hr_matrix <- function(stats_df, subtype_stats,
                                             gse_clinical, tcga_clinical,
                                             panel_tbl, composite_defs, cfg) {
  clin <- list(GSE39582 = gse_clinical, "TCGA-COAD" = tcga_clinical)

  sd_of <- function(ds, project, score, metric) {
    d <- clin[[ds]]; mc <- metric_cols(metric)
    if (!all(c(score, mc[["t"]], mc[["e"]]) %in% names(d))) return(NA_real_)
    mem <- .subgroup_members(d, project, ds)
    d <- d[mem, , drop = FALSE]
    v <- suppressWarnings(as.numeric(d[[score]]))
    ok <- !is.na(v) & !is.na(d[[mc[["t"]]]]) & !is.na(d[[mc[["e"]]]])
    .score_sd(v[ok])
  }

  cell <- function(ds, project, test, metric, score) {
    s <- stats_df[stats_df$Dataset == ds & stats_df$Project == project & stats_df$Test == test &
                  stats_df$Metric == metric & stats_df$Score == score, ]
    if (!nrow(s)) return(NULL)
    s[1, ]
  }

  proj <- unique(stats_df[, c("Dataset", "Project")])
  proj$Category <- vapply(proj$Project, .subgroup_category, character(1))
  proj <- proj[!is.na(proj$Category), , drop = FALSE]
  if (!nrow(proj)) stop("[build_fig3_3a] no subgroup cohorts found in stats_df")
  proj$Subgroup <- vapply(proj$Project, .subgroup_label, character(1))
  proj$dtag     <- ifelse(proj$Dataset == "GSE39582", "GSE", "TCGA")
  proj$ylab     <- paste0(unname(DATASET_SHORT[proj$Dataset]), " \u00B7 ", proj$Subgroup)
  CAT_LEVELS    <- intersect(c("CMS", "PDS", "Stage", "MSI"), unique(proj$Category))

  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  core_scores <- intersect(paste0(c("Up", "Down", "Composite"), suffix), unique(stats_df$Score))
  if (length(core_scores) == 0) core_scores <- unique(stats_df$Score)
  metrics <- c("OS", "RFS")

  rows <- list()
  for (i in seq_len(nrow(proj))) {
    p <- proj[i, ]
    for (sc in core_scores) for (mt in metrics) {
      cx <- cell(p$Dataset, p$Project, "Cox",      mt, sc)
      wl <- cell(p$Dataset, p$Project, "Wilcoxon", mt, sc)
      km <- cell(p$Dataset, p$Project, "KM",       mt, sc)
      if (is.null(cx)) next
      testable <- isTRUE(as.logical(cx$Is_Testable))
      sd_sc <- if (testable) sd_of(p$Dataset, p$Project, sc, mt) else NA_real_
      to_sd <- function(x) {
        x <- suppressWarnings(as.numeric(x))
        if (!testable || is.na(x) || is.na(sd_sc)) NA_real_ else exp(log(x) * sd_sc)
      }
      hr_sd <- to_sd(cx$HR)
      lo_sd <- to_sd(cx$HR_lower)
      hi_sd <- to_sd(cx$HR_upper)

      rows[[length(rows) + 1]] <- tibble::tibble(
        Dataset = p$Dataset, Category = p$Category, Subgroup = p$Subgroup, ylab = p$ylab,
        Score = sc, Endpoint = mt, N = as.integer(cx$N), N_events = as.integer(cx$N_events),
        Is_Testable = testable,
        Effect_r = suppressWarnings(as.numeric(if (!is.null(wl$Effect_r)) wl$Effect_r else NA)),
        Wilcoxon_FDR_P = suppressWarnings(as.numeric(if (!is.null(wl$FDR_P)) wl$FDR_P else NA)),
        HR_perSD = hr_sd, HR_lower_perSD = lo_sd, HR_upper_perSD = hi_sd,
        log2HR = ifelse(is.na(hr_sd), NA_real_, log2(hr_sd)),
        Cox_Raw_P = suppressWarnings(as.numeric(cx$Raw_P)),
        Cox_FDR_P = suppressWarnings(as.numeric(cx$FDR_P)),
        Cox_sig = isTRUE(as.logical(cx$Is_Significant)),
        C_index = suppressWarnings(as.numeric(cx$C_index)),
        KM_logrank_FDR_P = suppressWarnings(as.numeric(if (!is.null(km$FDR_P)) km$FDR_P else NA))
      )
    }
  }
  cells <- dplyr::bind_rows(rows)
  cells$Category <- factor(cells$Category, levels = CAT_LEVELS)
  cells$ylab     <- factor(cells$ylab, levels = rev(unique(proj$ylab[order(proj$Category, proj$dtag, proj$Subgroup)])))

  # Panel A Plot
  A_SCORES <- intersect(paste0(c("Up", "Composite"), suffix), core_scores)
  if (length(A_SCORES) == 0) A_SCORES <- core_scores[1:min(2, length(core_scores))]

  ck_lab <- c(
    setNames("DTP Up\n(3-yr OS)", paste0(paste0("Up", suffix), "|OS")),
    setNames("DTP Up\n(3-yr RFS)", paste0(paste0("Up", suffix), "|RFS")),
    setNames("DTP Comp.\n(3-yr OS)", paste0(paste0("Composite", suffix), "|OS")),
    setNames("DTP Comp.\n(3-yr RFS)", paste0(paste0("Composite", suffix), "|RFS"))
  )

  pa <- cells %>%
    dplyr::filter(Score %in% A_SCORES) %>%
    dplyr::mutate(
      ckey = factor(paste0(Score, "|", Endpoint), levels = names(ck_lab)),
      tile_clean = ifelse(!Is_Testable, "n/t", ifelse(Cox_sig, .sig_stars(Cox_FDR_P, ""), ""))
    )

  lim <- min(max(abs(pa$log2HR[is.finite(pa$log2HR)]), na.rm = TRUE), 2)
  if (!is.finite(lim) || lim == 0) lim <- 1.0

  min_ev  <- if (!is.null(cfg$min_events)) cfg$min_events else 5
  min_cox <- if (!is.null(cfg$min_cox_n)) cfg$min_cox_n else 10

  p_mat <- ggplot2::ggplot(pa, ggplot2::aes(x = ckey, y = ylab, fill = pmax(pmin(log2HR, lim), -lim))) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = tile_clean),
      colour = ifelse(pa$tile_clean == "n/t", "grey45", "grey10"),
      fontface = "bold", size = 3.1, vjust = 0.5, lineheight = 0.82
    ) +
    ggplot2::facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradient2(
      low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
      limits = c(-lim, lim), name = expression(log[2]~"HR (per 1 SD of score)"),
      guide = ggplot2::guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1), na.value = "grey85"
    ) +
    ggplot2::scale_x_discrete(labels = ck_lab) +
    ggplot2::labs(
      title = "DTP score vs 3-year outcome within molecular and clinical subgroups",
      subtitle = paste0(
        "Univariable Cox HR per 1 SD (follow-up truncated at 36 mo). * = FDR<0.05.\n",
        "Grey \"n/t\" = not testable (< ", min_ev, " events / ", min_cox,
        " n), distinct from tested-but-null.\n",
        "Most subgroups are underpowered: a non-significant or not-testable cell is NOT evidence of absence."
      ),
      x = NULL, y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, lineheight = 0.85),
      axis.text.y = ggplot2::element_text(size = 8),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 9),
      panel.spacing.y = grid::unit(4, "pt")
    )

  # Table 3.3
  rd2 <- function(x) round(x, 2)
  rd3 <- function(x) round(x, 3)

  tbl <- cells %>%
    dplyr::transmute(
      Dataset, Category = as.character(Category), Subgroup, Score, Endpoint,
      N, N_events, Is_Testable,
      Effect_r = rd2(Effect_r), Wilcoxon_FDR_P = rd3(Wilcoxon_FDR_P),
      HR_perSD = rd2(HR_perSD), HR_lower = rd2(HR_lower_perSD), HR_upper = rd2(HR_upper_perSD),
      Cox_FDR_P = rd3(Cox_FDR_P), C_index = rd2(C_index),
      KM_logrank_FDR_P = rd3(KM_logrank_FDR_P)
    ) %>%
    dplyr::arrange(
      factor(Category, levels = CAT_LEVELS), Subgroup,
      factor(Score, levels = core_scores), factor(Endpoint, levels = metrics)
    )

  list(plot = p_mat, table = tbl)
}

#' Build Figure 3_3B: Score across molecular subtypes Kruskal-Wallis heatmap
#'
#' @param subtype_stats Kruskal-Wallis stats tibble from [run_crc_survival()].
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param cfg Config from [dtp_config()].
#' @return A ggplot object.
#' @export
build_fig3_3b_score_across_subtype <- function(subtype_stats, panel_tbl, composite_defs, cfg) {
  .subtype_eps2_heatmap(subtype_stats, panel_tbl, composite_defs, cfg)
}

#' Build Figure 3_3C: Subtype pairwise violin composite and Table 3.3C
#'
#' @param gse_clinical Annotated clinical tibble for GSE39582.
#' @param tcga_clinical Annotated clinical tibble for TCGA-COAD.
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param cfg Config from [dtp_config()].
#' @return Named list with `plot` (patchwork composite) and `table` (pairwise stats tibble).
#' @export
build_fig3_3c_subtype_violins <- function(gse_clinical, tcga_clinical,
                                          panel_tbl, composite_defs, cfg) {
  clin_list <- list(GSE39582 = gse_clinical, "TCGA-COAD" = tcga_clinical)
  axes <- names(.SUBTYPE_LEVELS)

  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  core_scores <- intersect(
    paste0(c("Up", "Down", "Composite"), suffix),
    unique(c(names(gse_clinical), names(tcga_clinical)))
  )
  if (length(core_scores) == 0) {
    core_scores <- grep(paste0(suffix, "$"), names(gse_clinical), value = TRUE)
  }
  MOD <- .crc_modules(core_scores, panel_tbl, composite_defs, cfg)

  # Pairwise table
  pw <- do.call(rbind, lapply(names(clin_list), function(ds) {
    do.call(rbind, lapply(axes, function(ax) {
      do.call(rbind, lapply(core_scores, function(sc) {
        .subtype_pair_stats(clin_list[[ds]], ds, ax, sc, cfg = cfg)
      }))
    }))
  }))

  pw <- pw %>%
    dplyr::group_by(Dataset, Subtype_Axis, Score) %>%
    dplyr::mutate(FDR_P = stats::p.adjust(Raw_P, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      FDR_P_global = stats::p.adjust(Raw_P, method = "BH"),
      Is_Significant = !is.na(FDR_P) & FDR_P < 0.05
    )

  out_tbl <- pw %>%
    dplyr::transmute(
      Dataset, Subtype_Axis, Score, Group_A, Group_B, N_A, N_B,
      Median_A = round(Median_A, 4), Median_B = round(Median_B, 4),
      Effect_r = round(Effect_r, 3), Raw_P, FDR_P, FDR_P_global,
      Is_Testable, Is_Significant
    )

  base_ylim <- function(dataset, score) {
    v <- suppressWarnings(as.numeric(clin_list[[dataset]][[score]]))
    v <- v[!is.na(v)]
    rng <- range(v)
    pad <- diff(rng) * 0.04
    c(rng[1] - pad, rng[2] + pad)
  }

  make_cell <- function(score, dataset, axis, is_top_row, is_left_col) {
    lv <- .SUBTYPE_LEVELS[[axis]]
    cl <- clin_list[[dataset]]
    d  <- cl[cl[[axis]] %in% lv, c(axis, score), drop = FALSE]
    d[[score]] <- suppressWarnings(as.numeric(d[[score]]))
    d  <- d[!is.na(d[[score]]), , drop = FALSE]
    d[[axis]] <- factor(d[[axis]], levels = lv)
    tb   <- table(d[[axis]])
    lbls <- setNames(paste0(names(tb), "\n(n=", as.integer(tb), ")"), names(tb))

    bylim <- base_ylim(dataset, score)
    sig   <- pw[pw$Dataset == dataset & pw$Subtype_Axis == axis & pw$Score == score &
                pw$Is_Significant, , drop = FALSE]
    br    <- .violin_brackets(sig, lv, bylim[2], diff(bylim))

    score_title <- if (score %in% names(MOD$labels)) MOD$labels[[score]] else score

    ggplot2::ggplot(d, ggplot2::aes(x = .data[[axis]], y = .data[[score]], fill = .data[[axis]])) +
      ggplot2::geom_violin(alpha = 0.8, trim = FALSE, scale = "width", linewidth = 0.3, colour = "grey30") +
      ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, linewidth = 0.3) +
      br$layers +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::scale_x_discrete(labels = lbls) +
      ggplot2::coord_cartesian(ylim = c(bylim[1], br$ylim_top)) +
      ggplot2::labs(
        title = if (is_top_row) axis else NULL,
        x = NULL,
        y = if (is_left_col) paste0(score_title, " ssGSEA") else NULL
      ) +
      .base_theme +
      ggplot2::theme(
        legend.position = "none",
        plot.title  = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
        axis.text.x = ggplot2::element_text(size = 7.5, lineheight = 0.85)
      )
  }

  .panel_header <- function(label) {
    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = label, fontface = "bold", size = 4.6, hjust = 0.5) +
      ggplot2::theme(plot.margin = ggplot2::margin(t = 2, b = 4))
  }

  make_panel <- function(dataset) {
    cells <- list()
    for (sc in core_scores) {
      is_top <- sc == core_scores[[1]]
      for (i in seq_along(axes)) {
        cells[[length(cells) + 1]] <- make_cell(sc, dataset, axes[[i]], is_top, is_left_col = (i == 1))
      }
    }
    grid <- patchwork::wrap_plots(cells, ncol = length(axes), byrow = TRUE)
    .panel_header(DATASET_FULL[[dataset]]) / grid + patchwork::plot_layout(heights = c(0.045, 1))
  }

  CAP <- paste0(
    "(A) Marisa (GSE39582); (B) TCGA-COAD. CMS and PDS subtype axes within each panel. CMS-unclassified and PDS \"Mixed\" are excluded.\n",
    "Pairwise Mann-Whitney tests with rank-biserial r, BH-corrected within each cohort x axis x score;\n",
    "only FDR<0.05 pairs are bracketed. Full pairwise table: Table_3_3C_subtype_pairwise.csv."
  )

  panel_A <- patchwork::wrap_elements(make_panel("GSE39582"))
  panel_B <- patchwork::wrap_elements(make_panel("TCGA-COAD"))

  comp <- patchwork::wrap_plots(panel_A, panel_B, ncol = 2) +
    patchwork::plot_annotation(
      tag_levels = "A",
      caption = CAP,
      theme = ggplot2::theme(
        plot.title   = ggplot2::element_text(face = "bold", size = 15, hjust = 0),
        plot.caption = ggplot2::element_text(size = 8, colour = "grey30", hjust = 0)
      )
    ) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 15))

  list(plot = comp, table = out_tbl)
}

#' Build Figure 5: Confounder adjustment summary heatmap and caption
#'
#' @param adj_df Adjusted Cox statistics summary tibble from [run_crc_survival()].
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param cfg Config from [dtp_config()].
#' @param primary_score Primary score column to summarize (e.g. "Up_ssGSEA").
#' @return Named list with `plot` (ggplot object) and `caption` (character string).
#' @export
build_fig5_confounding_summary <- function(adj_df, panel_tbl, composite_defs, cfg,
                                          primary_score = "Up_ssGSEA") {
  MODEL_LEVELS  <- c("Clinicopath", "CMS_adjusted", "PDS_adjusted")
  MODEL_LABEL   <- c(Clinicopath = "Clinicopathology", CMS_adjusted = "+CMS", PDS_adjusted = "+PDS")
  DS_LEVELS     <- c("GSE39582", "GSE39582 (treated)", "TCGA-COAD")
  METRIC_LEVELS <- c("OS", "RFS")

  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  score_lab <- sub(paste0(suffix, "$"), "", primary_score)

  num <- function(x) suppressWarnings(as.numeric(x))

  d <- adj_df %>%
    dplyr::filter(
      Score == primary_score, Model %in% MODEL_LEVELS,
      Dataset %in% DS_LEVELS, Metric %in% METRIC_LEVELS
    ) %>%
    dplyr::mutate(
      Model   = factor(Model, levels = MODEL_LEVELS),
      Dataset = factor(Dataset, levels = DS_LEVELS),
      Metric  = factor(Metric, levels = METRIC_LEVELS),
      FDR_P   = num(FDR_P),
      LRT_score_P = num(LRT_score_P),
      Delta_logHR_pct = num(Delta_logHR_pct),
      Is_Testable = as.logical(Is_Testable)
    )
  if (!nrow(d)) stop("[build_fig5] no rows for score '", primary_score, "' in adj_df")

  pa <- d %>%
    dplyr::mutate(
      lab_clean = ifelse(!Is_Testable, "n/t", paste0(sprintf("%.3f", FDR_P), .sig_stars(FDR_P)))
    )

  lim <- ceiling(max(abs(d$Delta_logHR_pct), na.rm = TRUE) / 10) * 10
  if (!is.finite(lim) || lim == 0) lim <- 50
  dark_at <- 0.65 * lim

  pa <- pa %>%
    dplyr::mutate(
      lab_col = ifelse(!Is_Testable, "grey45",
                ifelse(!is.na(Delta_logHR_pct) & abs(Delta_logHR_pct) > dark_at, "white", "grey10"))
    )

  p <- ggplot2::ggplot(pa, ggplot2::aes(x = Model, y = forcats::fct_rev(Metric), fill = Delta_logHR_pct)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = lab_clean), colour = pa$lab_col, fontface = "bold",
      size = 2.9, lineheight = 0.8
    ) +
    ggplot2::facet_grid(Dataset ~ ., switch = "y", labeller = ggplot2::labeller(Dataset = DATASET_FULL)) +
    ggplot2::scale_fill_gradient2(
      low = "#3B6EA5", mid = "white", high = "#E64B35", midpoint = 0,
      limits = c(-lim, lim), na.value = "grey85",
      name = "Change in score log-HR after adjustment (%)  [-] attenuated  <->  [+] strengthened",
      guide = ggplot2::guide_colourbar(barwidth = 12, barheight = 0.5, title.position = "top", title.hjust = 0.5)
    ) +
    ggplot2::scale_x_discrete(labels = MODEL_LABEL) +
    ggplot2::scale_y_discrete(labels = ENDPOINT_LABEL) +
    ggplot2::labs(
      title = "Independence of the DTP score after adjustment",
      subtitle = paste0(
        "Tile colour = % change in the DTP ", score_lab,
        " score log-HR after adjustment (blue = attenuated toward null, red = strengthened).\n",
        "Number = adjusted-model LRT BH-FDR (score vs covariate-only model); * marks FDR < 0.05.\n",
        "Per 1 SD, follow-up truncated at 36 mo.  |shift| > 10% ~ meaningful confounding (Greenland).  Exploratory."
      ),
      x = NULL, y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 8.5),
      plot.title.position = "plot",
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 9),
      panel.spacing.y = grid::unit(4, "pt")
    )

  cap <- c(
    "Figure 5. Independence of the DTP score after adjustment (Section 3.4; exploratory).",
    "",
    paste0("Adjusted-Cox confounding heatmap. Each cell is the DTP ", score_lab,
           " score term from a Cox model (follow-up truncated at 36 months) of overall"),
    "survival (OS) or recurrence-free survival (RFS), adjusted for clinicopathology (stage + microsatellite",
    "status), CMS subtype, or PDS subtype, in GSE39582 (whole cohort and adjuvant-treated subset) and TCGA-COAD.",
    "",
    "Tile colour is the percent change in the score's log-hazard ratio (the Cox coefficient) after adjustment,",
    "relative to the unadjusted model: blue = attenuated toward the null (the adjustor absorbs part of the",
    "signal), red = strengthened, white = unchanged. Following Greenland's rule of thumb a shift beyond +/-10%",
    "is treated as materially confounded; cells within +/-10% therefore render near-white.",
    "",
    "The cell number is the Benjamini-Hochberg FDR of that score term's likelihood-ratio test (the adjusted model",
    "with vs. without the score; residual significance after adjustment) -- the same nested-model test the",
    "effect-modification analysis uses, and more robust than the Wald test for the large per-SD hazard ratios",
    "here. An asterisk marks FDR < 0.05. The FDR family is the three DTP scores (Up / Composite / Down) within",
    "each cohort x endpoint x model. A grey \"n/t\" cell denotes a model the events-per-variable gate (>=10 events",
    "per parameter, n>=10, >=5 events) never fitted (not modelled, as opposed to fitted-but-null); all cells here",
    "were testable.",
    "",
    "Exploratory analysis: adjustment removes only measured confounders, so residual confounding remains",
    "possible; multiple-testing correction reduces but does not eliminate false positives."
  )

  list(plot = p, caption = paste(cap, collapse = "\n"))
}

#' Orchestrate generation and saving of all 7 CRC composite figures
#'
#' Evaluates all 7 composite figure builders using in-memory results from
#' [run_crc_survival()], saves PDF and PNG copies into `out_dir`, and returns
#' derived subgroup tables and caption text.
#'
#' @param crc_result Result list from [run_crc_survival()].
#' @param panel_tbl Canonical signature panel table.
#' @param composite_defs Composite definitions table.
#' @param cfg Config from [dtp_config()].
#' @param out_dir Output directory path for figures.
#' @return A named list containing:
#'   \item{subtype_survival_tbl}{Table 3.3 subgroup univariable survival tibble.}
#'   \item{subtype_pairwise_tbl}{Table 3.3C subtype pairwise statistics tibble.}
#'   \item{fig5_caption}{Figure 5 caption character string.}
#' @export
build_crc_composites <- function(crc_result, panel_tbl, composite_defs, cfg,
                                 out_dir = file.path("output", "figures")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("Building CRC composite figures -> ", out_dir)

  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  primary_score <- paste0("Up", suffix)

  results <- list()
  status  <- list()

  # Group 1: Outcomes
  status[["1 (outcomes)"]] <- tryCatch({
    p1 <- build_fig1_outcome_composite(
      crc_result$stats_df, crc_result$gse_clinical, crc_result$tcga_clinical,
      panel_tbl, composite_defs, cfg
    )
    score_cols <- unique(c(
      grep(paste0(suffix, "$"), names(crc_result$gse_clinical), value = TRUE),
      grep(paste0(suffix, "$"), names(crc_result$tcga_clinical), value = TRUE)
    ))
    n_rows <- length(unique(score_cols))
    h1 <- 4.6 + (1.5 * n_rows + 1.4) + 0.6
    save_fig(p1, "Fig1_outcome_composite", out_dir, 11.5, h1)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[composite figures] group 1 FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Group 2: KM
  status[["2 (KM)"]] <- tryCatch({
    p2 <- build_fig2_km_composite(
      crc_result$stats_df, crc_result$gse_clinical, crc_result$tcga_clinical,
      panel_tbl, composite_defs, cfg
    )
    score_cols <- unique(c(
      grep(paste0(suffix, "$"), names(crc_result$gse_clinical), value = TRUE),
      grep(paste0(suffix, "$"), names(crc_result$tcga_clinical), value = TRUE)
    ))
    n_rows <- length(unique(score_cols))
    h2 <- 4.6 + (2.0 * n_rows + 1.4) + 0.6
    save_fig(p2, "Fig2_km_composite", out_dir, 11.5, h2)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[composite figures] group 2 FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Group 3: Subgroups
  status[["3 (subgroups)"]] <- tryCatch({
    p3 <- build_fig3_subgroup_composite(
      crc_result$int_df, crc_result$level_df, panel_tbl, composite_defs, cfg,
      primary_score = primary_score,
      gse_clinical = crc_result$gse_clinical, tcga_clinical = crc_result$tcga_clinical
    )
    save_fig(p3, "Fig3_subgroup_composite", out_dir, 14.0, 12.2)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[composite figures] group 3 FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Group 4: Subtype survival matrix (Fig 3_3A)
  status[["4 (subtype survival)"]] <- tryCatch({
    res4 <- build_fig3_3a_subgroup_hr_matrix(
      crc_result$stats_df, crc_result$subtype_stats,
      crc_result$gse_clinical, crc_result$tcga_clinical,
      panel_tbl, composite_defs, cfg
    )
    results$subtype_survival_tbl <- res4$table
    n_sub <- length(unique(res4$table$Subgroup))
    h4 <- max(4.5, 0.32 * n_sub + 2.6)
    save_fig(res4$plot, "Fig3_3A_subgroup_hr_matrix", out_dir, 8.8, h4)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[composite figures] group 4 FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Group 5: Score across subtypes (Fig 3_3B)
  status[["4b (subtype score heatmap)"]] <- tryCatch({
    p5 <- build_fig3_3b_score_across_subtype(crc_result$subtype_stats, panel_tbl, composite_defs, cfg)
    save_fig(p5, "Fig3_3B_score_across_subtype", out_dir, 7.5, 3.2)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[composite figures] group 4b FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Group 6: Subtype violins (Fig 3_3C)
  status[["4c (subtype violins)"]] <- tryCatch({
    res6 <- build_fig3_3c_subtype_violins(
      crc_result$gse_clinical, crc_result$tcga_clinical,
      panel_tbl, composite_defs, cfg
    )
    results$subtype_pairwise_tbl <- res6$table
    save_fig(res6$plot, "Fig3_3C_subtype_violins", out_dir, 13.0, 15.5)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[composite figures] group 4c FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Group 7: Confounding summary (Fig 5)
  status[["5 (confounding)"]] <- tryCatch({
    res7 <- build_fig5_confounding_summary(
      crc_result$adj_df, panel_tbl, composite_defs, cfg,
      primary_score = primary_score
    )
    results$fig5_caption <- res7$caption
    save_fig(res7$plot, "Fig5_confounding_summary", out_dir, 9.5, 5.6)
    writeLines(res7$caption, file.path(out_dir, "Fig5_confounding_summary_caption.txt"))
    TRUE
  }, error = function(e) {
    msg <- sprintf("[composite figures] group 5 FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  n_ok <- sum(unlist(status))
  n_tot <- length(status)
  summ <- sprintf("Composite figures: %d/%d groups OK%s", n_ok, n_tot,
    if (n_ok < n_tot) paste0(" -- FAILED: ", paste(names(status)[!unlist(status)], collapse = ", ")) else "")
  message(summ)
  if (n_ok < n_tot) warning(summ, call. = FALSE)

  results
}
