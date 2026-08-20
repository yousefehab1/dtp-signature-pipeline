# ==============================================================================
# R/mets_composites.R - Metastasis GSEA and Liver Purity Composite Figures
#
# Generates publication-ready figures directly from in-memory metastasis DE/GSEA
# analysis objects returned by run_mets_de() (pure presentation layer, no disk
# reads).
# ==============================================================================

.fmt_p <- function(prefix, p) {
  if (length(p) != 1 || is.na(p)) {
    NA_character_
  } else {
    paste0(prefix, if (p < 0.001) "<.001" else sub("^0", "", sprintf("=%.3f", p)))
  }
}

# ---------------------------------------------------------------------------
# .mets_nes_heatmap() - one-column NES heatmap (fill = NES, label = NES + raw p +
# FDR stars/q) over all signatures for a single GSEA contrast. Signatures are
# ordered top->bottom by `sig_order`; fill limits `lim` are shared across the
# two contrasts so the colours are directly comparable. The NES sign convention
# lives on the x-axis title so it remains clear in the figure.
# ---------------------------------------------------------------------------
.mets_nes_heatmap <- function(gsea_df, contrast_lab, pos_label, sig_order, lim) {
  d <- gsea_df %>% dplyr::filter(.data$ID %in% sig_order)
  d$Signature <- factor(d$ID, levels = rev(sig_order))
  d$Contrast  <- contrast_lab

  pval_vec <- if ("pvalue" %in% names(d)) d$pvalue else rep(NA_real_, nrow(d))
  fdr_vec  <- if ("p.adjust" %in% names(d)) d$p.adjust else rep(NA_real_, nrow(d))

  d$lab <- mapply(function(nes, pval, fdr) {
    parts <- c(
      if (!is.na(nes)) sprintf("%.2f", nes) else NA_character_,
      if (!is.na(pval)) .fmt_p("p", pval) else NA_character_,
      if (!is.na(fdr)) paste0(.fmt_p("q", fdr), .sig_stars(fdr, "")) else NA_character_
    )
    parts <- parts[!is.na(parts) & nzchar(parts)]
    paste(parts, collapse = "\n")
  }, d$NES, pval_vec, fdr_vec)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$Contrast, y = .data$Signature, fill = .data$NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$lab), size = 2.3,
                       fontface = "bold", colour = "grey10", lineheight = 0.82) +
    ggplot2::scale_fill_gradient2(
      low = "#3B6EA5", mid = "white", high = "#E64B35",
      midpoint = 0, limits = c(-lim, lim), name = "NES"
    ) +
    ggplot2::labs(
      title = "NES + FDR",
      x = sprintf("+NES -> up in %s", pos_label),
      y = NULL
    ) +
    pub_theme +
    ggplot2::theme(
      legend.position = "right",
      panel.grid   = ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 9, face = "bold"),
      plot.title   = ggplot2::element_text(size = 11, face = "bold", hjust = 0)
    )
}

# ---------------------------------------------------------------------------
# build_mets_gsea_composite()
# ---------------------------------------------------------------------------

#' Build metastasis GSEA 4-panel composite figure
#'
#' Arranges Primary-vs-Normal and Primary-vs-Metastasis (purity-adjusted) GSEA
#' running enrichment score plots alongside corresponding NES/FDR heatmaps.
#'
#' The whole figure is anchored on the primary tumour, matching the thesis
#' narrative (primary -> metastasis): panels A/B read Primary vs Normal (positive
#' NES = up in Primary) and panels C/D read Primary vs Metastasis. Panels C/D read
#' the liver-score-adjusted ranking: the metastasis arm carries heavy hepatic
#' signal (liver-validation figure), so the metastasis-direction enrichment must
#' be read off the adjusted fit. Panels A/B (Primary-vs-Normal) have no liver
#' tissue and stay unadjusted.
#'
#' @param mets_result Result list from [run_mets_de()].
#' @param cfg Config from [dtp_config()].
#' @return A patchwork object.
#' @export
build_mets_gsea_composite <- function(mets_result, cfg = dtp_config()) {
  gsea_pn     <- mets_result$gsea_pn
  gsea_pm     <- mets_result$gsea_pm_adj
  gsea_pn_obj <- mets_result$gsea_pn_obj
  gsea_pm_obj <- mets_result$gsea_pm_adj_obj

  if (is.null(gsea_pn) || is.null(gsea_pm) || is.null(gsea_pn_obj) || is.null(gsea_pm_obj)) {
    stop("build_mets_gsea_composite(): missing required GSEA results in mets_result.")
  }

  df_pn <- as.data.frame(gsea_pn)
  df_pm <- as.data.frame(gsea_pm)

  if (!nrow(df_pn) || !nrow(df_pm)) {
    stop("build_mets_gsea_composite(): empty GSEA result(s).")
  }

  # Shared signature order (top->bottom) = descending Primary-vs-Metastasis NES,
  # the headline contrast; any signature absent there is appended.
  ord       <- df_pm$ID[order(-df_pm$NES)]
  sig_order <- c(ord, setdiff(union(df_pn$ID, df_pm$ID), ord))
  glim      <- max(abs(c(df_pn$NES, df_pm$NES)), na.rm = TRUE)
  if (is.na(glim) || glim <= 0) glim <- 1

  # sig_order is the UNION of the two results, so it can name a signature that one
  # object lacks (a set can drop out of one contrast and not the other). Each
  # enrichment plot may only be asked for IDs its own object holds -
  # .gsea_running_data() stop()s otherwise, which would take the whole composite down.
  # The shared order is preserved; only the missing rows are omitted, and flagged.
  pn_ids <- intersect(sig_order, df_pn$ID)
  pm_ids <- intersect(sig_order, df_pm$ID)
  if (length(pn_ids) < length(sig_order) || length(pm_ids) < length(sig_order)) {
    message("  [mets GSEA] signature(s) absent from one contrast, omitted from that panel: ",
            paste(setdiff(sig_order, intersect(pn_ids, pm_ids)), collapse = ", "))
  }

  A <- plot_gsea_enrichment(
    gsea_pn_obj, gene_set_ids = pn_ids,
    title = "Primary vs Normal: all signatures",
    label_left = "Up in Primary", label_right = "Up in Normal"
  )
  b <- .mets_nes_heatmap(df_pn, "P / N", pos_label = "Primary", sig_order = sig_order, lim = glim)

  C <- plot_gsea_enrichment(
    gsea_pm_obj, gene_set_ids = pm_ids,
    title = "Primary vs Metastasis, purity-adjusted for liver score: all signatures",
    label_left = "Up in Metastasis", label_right = "Up in Primary"
  )
  d <- .mets_nes_heatmap(df_pm, "P / M (adj)", pos_label = "Metastasis, purity-adjusted",
                         sig_order = sig_order, lim = glim)

  p_comp <- patchwork::wrap_plots(
    patchwork::wrap_elements(A), patchwork::wrap_elements(b),
    patchwork::wrap_elements(C), patchwork::wrap_elements(d),
    ncol = 2, widths = c(3, 1)
  ) +
    patchwork::plot_annotation(
      tag_levels = "A",
      title = "Mets GSEA: primary tumour vs normal mucosa and vs liver metastasis (metastasis arm purity-adjusted for liver score)",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15))
    ) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))

  p_comp
}

# ---------------------------------------------------------------------------
# build_mets_liver_composite()
# ---------------------------------------------------------------------------

#' Build metastasis liver purity composite figure
#'
#' Generates a 2-panel composite figure showing:
#' (A) Per-patient liver ssGSEA score trajectory across Normal, Primary, and Metastasis tissues.
#' (B) Primary-vs-Metastasis GSEA NES comparison before vs after liver-purity adjustment.
#'
#' @param mets_result Result list from [run_mets_de()].
#' @return A patchwork object.
#' @export
build_mets_liver_composite <- function(mets_result) {
  liver_df <- mets_result$liver_df
  cmp      <- mets_result$pvm_cmp

  if (is.null(liver_df) || is.null(cmp)) {
    stop("build_mets_liver_composite(): missing liver_df or pvm_cmp in mets_result.")
  }

  # -- Panel A: per-patient liver-score trajectory (the purity problem) --------
  liver_df <- as.data.frame(liver_df)
  liver_df$Tissue <- factor(liver_df$Tissue, levels = c("Normal", "Primary", "Metastasis"))
  pA <- ggplot2::ggplot(liver_df, ggplot2::aes(x = .data$Tissue, y = .data$Liver, group = .data$Patient)) +
    ggplot2::geom_line(colour = "gray60", alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Tissue), size = 3) +
    ggplot2::scale_colour_manual(values = pub_palette) +
    ggplot2::labs(
      title = "Liver signal across Normal -> Primary -> Metastasis",
      subtitle = "HSIAO_LIVER_SPECIFIC_GENES ssGSEA per sample; lines join each patient's three tissues.",
      x = NULL, y = "Liver score"
    ) +
    .base_theme +
    ggplot2::theme(legend.position = "none")

  # -- Panel B: unadjusted vs purity-adjusted P/M NES (what adjustment changes) -
  cmp <- as.data.frame(cmp)
  hm <- dplyr::bind_rows(
    dplyr::transmute(cmp, ID = .data$ID, Model = "Unadjusted", NES = .data$NES_unadjusted, FDR = .data$FDR_unadjusted),
    dplyr::transmute(cmp, ID = .data$ID, Model = "Adjusted",   NES = .data$NES_adjusted,   FDR = .data$FDR_adjusted)
  ) %>%
    dplyr::mutate(
      Model = factor(.data$Model, levels = c("Unadjusted", "Adjusted")),
      Signature = factor(.data$ID, levels = rev(cmp$ID[order(cmp$NES_adjusted)])),
      lab = mapply(function(nes, fdr) {
        parts <- c(
          if (!is.na(nes)) sprintf("%.2f", nes) else NA_character_,
          if (!is.na(fdr)) paste0(.fmt_p("q", fdr), .sig_stars(fdr, "")) else NA_character_
        )
        parts <- parts[!is.na(parts) & nzchar(parts)]
        paste(parts, collapse = "\n")
      }, .data$NES, .data$FDR)
    )

  lim <- max(abs(hm$NES), na.rm = TRUE)
  if (is.na(lim) || lim <= 0) lim <- 1

  pB <- ggplot2::ggplot(hm, ggplot2::aes(x = .data$Model, y = .data$Signature, fill = .data$NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$lab), size = 2.5, fontface = "bold",
                       colour = "grey10", lineheight = 0.82) +
    ggplot2::scale_fill_gradient2(
      low = "#3B6EA5", mid = "white", high = "#E64B35",
      midpoint = 0, limits = c(-lim, lim), name = "NES"
    ) +
    ggplot2::labs(
      title = "Primary vs Metastasis NES: unadjusted vs purity-adjusted",
      subtitle = "NES with BH-FDR q per cell; * <0.05, ** <0.01, *** <0.001.",
      x = "Primary vs Metastasis ranking", y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "right"
    )

  CAP <- paste0(
    "(A) HSIAO_LIVER_SPECIFIC_GENES ssGSEA score per sample; lines join each patient's three tissues - metastases carry heavy hepatic signal.\n",
    "(B) Primary-vs-Metastasis GSEA NES before vs after adjusting the ranking for the liver score; FDR stars * <0.05, ** <0.01, *** <0.001 (BH)."
  )

  p_comp <- (patchwork::wrap_elements(pA) + patchwork::wrap_elements(pB)) +
    patchwork::plot_layout(widths = c(1.25, 1)) +
    patchwork::plot_annotation(
      tag_levels = "A",
      title   = "Liver purity in GSE50760 metastases and its effect on Primary-vs-Metastasis GSEA",
      caption = CAP,
      theme = ggplot2::theme(
        plot.title   = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
        plot.caption = ggplot2::element_text(hjust = 0, size = 8.5, colour = "grey30")
      )
    ) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16))

  p_comp
}

# ---------------------------------------------------------------------------
# build_mets_composites()
# ---------------------------------------------------------------------------

#' Orchestrate generation and saving of metastasis composite figures
#'
#' Evaluates metastasis figure builders using in-memory results from
#' [run_mets_de()], saves PDF and PNG copies into `out_dir`, and returns
#' an execution status list.
#'
#' @param mets_result Result list from [run_mets_de()].
#' @param cfg Config from [dtp_config()].
#' @param out_dir Output directory path for figures.
#' @return A named list containing:
#'   \item{status}{Named logical list of execution statuses.}
#' @export
build_mets_composites <- function(mets_result, cfg = dtp_config(),
                                  out_dir = file.path("output", "figures")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("Building metastasis composite figures -> ", out_dir)

  status <- list()

  # Figure 1: GSEA composite
  status[["GSEA composite"]] <- tryCatch({
    p_gsea <- build_mets_gsea_composite(mets_result, cfg = cfg)
    df_pm <- as.data.frame(mets_result$gsea_pm_adj)
    df_pn <- as.data.frame(mets_result$gsea_pn)
    n_sig <- length(union(df_pn$ID, df_pm$ID))
    h <- max(14, 2 * (n_sig * 0.55 + 4))
    save_fig(p_gsea, "Fig_Mets_GSEA_composite", out_dir, 16.0, h)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[mets composites] GSEA composite FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  # Figure 2: Liver purity composite
  status[["Liver purity composite"]] <- tryCatch({
    p_liver <- build_mets_liver_composite(mets_result)
    save_fig(p_liver, "Fig_Mets_LiverPurity_composite", out_dir, 12.5, 6.2)
    TRUE
  }, error = function(e) {
    msg <- sprintf("[mets composites] Liver purity composite FAILED: %s", conditionMessage(e))
    message("  ", msg); warning(msg, call. = FALSE)
    FALSE
  })

  n_ok  <- sum(unlist(status))
  n_tot <- length(status)
  summ  <- sprintf("Metastasis composites: %d/%d figures OK%s", n_ok, n_tot,
                   if (n_ok < n_tot) paste0(" -- FAILED: ", paste(names(status)[!unlist(status)], collapse = ", ")) else "")
  message(summ)
  if (n_ok < n_tot) warning(summ, call. = FALSE)
  message("Metastasis composite figures done -> ", out_dir)

  list(status = status)
}
