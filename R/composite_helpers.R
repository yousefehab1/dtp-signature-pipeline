# ==============================================================================
# R/composites/composite_helpers.R - Shared helpers for CRC composite figures:
# dynamic signature labeling/ordering, cohort membership classification,
# subtype pairing, and significance brackets.
# ==============================================================================

# Derive dynamic signature score columns, labels, and display order. Internal.
.crc_modules <- function(score_cols, panel_tbl, composite_defs, cfg = list()) {
  suffix <- if (!is.null(cfg$score_suffix)) cfg$score_suffix else "_ssGSEA"
  score_cols <- unique(score_cols)

  # Composite mappings
  comp_map <- list()
  comp_order <- character(0)
  default_comp_order <- character(0)
  if (!is.null(composite_defs) && nrow(composite_defs) > 0) {
    for (i in seq_len(nrow(composite_defs))) {
      cid <- composite_defs$Composite_ID[i]
      col_name <- paste0(cid, suffix)
      dname <- composite_defs$Display_Name[i]
      comp_map[[col_name]] <- dname
      is_def <- isTRUE(as.logical(composite_defs$Is_Default[i]))
      if (is_def) {
        default_comp_order <- c(default_comp_order, col_name)
      } else {
        comp_order <- c(comp_order, col_name)
      }
    }
  }

  # Base signature panel mappings (first appearance in panel_tbl)
  panel_map <- list()
  panel_order <- character(0)
  if (!is.null(panel_tbl) && nrow(panel_tbl) > 0) {
    sig_ids <- unique(panel_tbl$Signature_ID)
    for (sid in sig_ids) {
      col_name <- paste0(sid, suffix)
      sname <- panel_tbl$Signature_Name[match(sid, panel_tbl$Signature_ID)]
      panel_map[[col_name]] <- sname
      panel_order <- c(panel_order, col_name)
    }
  }

  # Order:
  # 1. Default composites present in score_cols
  # 2. Panel signatures present in score_cols (order of panel_tbl appearance)
  # 3. Non-default composites present in score_cols
  # 4. Unmatched novel score columns
  def_comp_present <- intersect(default_comp_order, score_cols)
  panel_present <- intersect(panel_order, score_cols)
  other_comp_present <- intersect(comp_order, score_cols)
  known_cols <- c(def_comp_present, panel_present, other_comp_present)
  unmatched_cols <- setdiff(score_cols, known_cols)

  cols <- c(known_cols, unmatched_cols)

  # Label lookup
  labels <- setNames(character(length(cols)), cols)
  for (col in cols) {
    if (col %in% names(comp_map)) {
      labels[[col]] <- comp_map[[col]]
    } else if (col %in% names(panel_map)) {
      labels[[col]] <- panel_map[[col]]
    } else {
      labels[[col]] <- col
    }
  }

  list(cols = cols, labels = labels, levels = unname(labels))
}

# Classify a project name into a subgroup category (CMS / PDS / Stage / MSI). Internal.
.subgroup_category <- function(project) {
  if (grepl("_CMS[1-4]$", project))                                return("CMS")
  if (grepl("_PDS[1-3]$", project) || grepl("_Mixed$", project))   return("PDS")
  if (grepl("_Stage[1-4]$", project))                              return("Stage")
  if (grepl("_All_MSS$|_All_MSI$", project))                       return("MSI")
  NA_character_
}

# Standardized human-readable label for a subgroup cohort project. Internal.
.subgroup_label <- function(project) {
  s <- gsub("_", " ", sub("^(GSE|TCGA)_", "", project))
  for (r in list(c("Stage1", "Stage I"), c("Stage2", "Stage II"),
                 c("Stage3", "Stage III"), c("Stage4", "Stage IV"))) {
    s <- sub(r[[1]], r[[2]], s, fixed = TRUE)
  }
  s
}

# Reconstruct subgroup membership boolean mask from clinical frame. Internal.
.subgroup_members <- function(clin, project, dataset) {
  na0 <- function(x) { x[is.na(x)] <- FALSE; x }
  if (grepl("_CMS[1-4]$", project)) return(na0(clin$CMS == sub(".*_", "", project)))
  if (grepl("_PDS[1-3]$", project)) return(na0(clin$PDS == sub(".*_", "", project)))
  if (grepl("_Mixed$", project))    return(na0(clin$PDS == "Mixed"))

  m <- if (dataset == "GSE39582") {
    switch(project,
      "GSE_All_MSS" = clin$MMR == "pMMR",
      "GSE_All_MSI" = clin$MMR == "dMMR",
      "GSE_Stage1"  = clin$TNM_stage %in% c("1", 1),
      "GSE_Stage2"  = clin$TNM_stage %in% c("2", 2),
      "GSE_Stage3"  = clin$TNM_stage %in% c("3", 3),
      "GSE_Stage4"  = clin$TNM_stage %in% c("4", 4),
      NULL
    )
  } else {
    switch(project,
      "TCGA_All_MSS" = clin$paper_MSI_status == "MSS",
      "TCGA_All_MSI" = clin$paper_MSI_status %in% c("MSI-H", "MSI-L"),
      "TCGA_Stage1"  = clin$Stage == "I",
      "TCGA_Stage2"  = clin$Stage == "II",
      "TCGA_Stage3"  = clin$Stage == "III",
      "TCGA_Stage4"  = clin$Stage == "IV",
      NULL
    )
  }
  if (is.null(m)) {
    stop("[subgroup_members] no membership rule for kept subgroup '", project, "'.")
  }
  na0(m)
}

# Within-subset score SD for rescaling continuous Cox unit HRs to per-1-SD HRs. Internal.
.score_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

# Subtype levels dictionary for CMS and PDS axes. Internal.
.SUBTYPE_LEVELS <- list(
  CMS = c("CMS1", "CMS2", "CMS3", "CMS4"),
  PDS = c("PDS1", "PDS2", "PDS3")
)

# Pairwise Mann-Whitney + rank-biserial r for one (Dataset, Axis, Score). Internal.
.subtype_pair_stats <- function(clin, dataset, axis, score, cfg = dtp_config()) {
  lv <- .SUBTYPE_LEVELS[[axis]]
  if (!all(c(axis, score) %in% colnames(clin))) return(NULL)
  d  <- clin[clin[[axis]] %in% lv, c(axis, score), drop = FALSE]
  d[[score]] <- suppressWarnings(as.numeric(d[[score]]))
  d  <- d[!is.na(d[[score]]), , drop = FALSE]
  rows <- lapply(utils::combn(lv, 2, simplify = FALSE), function(pr) {
    dd <- d[d[[axis]] %in% pr, , drop = FALSE]
    dd[[axis]] <- factor(dd[[axis]], levels = pr)
    a  <- dd[[score]][dd[[axis]] == pr[1]]
    b  <- dd[[score]][dd[[axis]] == pr[2]]
    ws <- get_wilcox_stats(dd, score, axis, cfg = cfg)
    data.frame(
      Dataset = dataset, Subtype_Axis = axis, Score = score,
      Group_A = pr[1], Group_B = pr[2], N_A = length(a), N_B = length(b),
      Median_A = if (length(a)) stats::median(a) else NA_real_,
      Median_B = if (length(b)) stats::median(b) else NA_real_,
      Effect_r = ws$r, Raw_P = ws$p, Is_Testable = !is.na(ws$p),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Nested significance brackets (clean-only, asterisk labels) for subtype violins. Internal.
.violin_brackets <- function(sig_df, lv, y0, span) {
  if (!nrow(sig_df)) return(list(layers = list(), ylim_top = y0 + span * 0.06))
  sig_df$xa <- match(sig_df$Group_A, lv)
  sig_df$xb <- match(sig_df$Group_B, lv)
  sig_df <- sig_df[order(abs(sig_df$xb - sig_df$xa)), ]
  step   <- span * 0.11
  tick_h <- step * 0.22
  layers <- list()
  for (i in seq_len(nrow(sig_df))) {
    r <- sig_df[i, ]
    y <- y0 + i * step
    lab <- .sig_stars(r$FDR_P, "")
    layers[[length(layers) + 1]] <- ggplot2::annotate(
      "segment",
      x = r$xa, xend = r$xb, y = y, yend = y, linewidth = 0.4, colour = "grey20"
    )
    layers[[length(layers) + 1]] <- ggplot2::annotate(
      "segment",
      x = c(r$xa, r$xb), xend = c(r$xa, r$xb), y = y, yend = y - tick_h,
      linewidth = 0.4, colour = "grey20"
    )
    layers[[length(layers) + 1]] <- ggplot2::annotate(
      "text",
      x = (r$xa + r$xb) / 2, y = y, label = lab, vjust = -0.15,
      size = 3, fontface = "bold", colour = "grey10", lineheight = 0.85
    )
  }
  list(layers = layers, ylim_top = y0 + nrow(sig_df) * step + step * 0.9)
}

# Score-across-subtype eps^2 Kruskal-Wallis heatmap (clean-only). Internal.
.subtype_eps2_heatmap <- function(subtype_stats, panel_tbl, composite_defs, cfg) {
  core_scores <- intersect(
    paste0(c("Up", "Down", "Composite"), cfg$score_suffix),
    unique(subtype_stats$Score)
  )
  if (length(core_scores) == 0) core_scores <- unique(subtype_stats$Score)
  MOD <- .crc_modules(core_scores, panel_tbl, composite_defs, cfg)

  pb <- subtype_stats %>%
    dplyr::filter(Score %in% MOD$cols, Subtype_Axis %in% c("CMS", "PDS")) %>%
    dplyr::mutate(
      Score = factor(MOD$labels[Score], levels = rev(MOD$levels)),
      Axis  = factor(paste0(Subtype_Axis, " subtype"), levels = c("CMS subtype", "PDS subtype")),
      eps2  = suppressWarnings(as.numeric(Eps2)),
      star  = .sig_stars(FDR_P, ""),
      lab_clean = ifelse(star == "", sprintf("%.02f", eps2),
                         paste0(sprintf("%.02f", eps2), " ", star))
    )

  ggplot2::ggplot(pb, ggplot2::aes(x = Dataset, y = Score, fill = eps2)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = lab_clean), size = 3.0, fontface = "bold",
      colour = "grey10", lineheight = 0.82
    ) +
    ggplot2::facet_grid(. ~ Axis) +
    ggplot2::scale_fill_gradient(
      low = "white", high = "#762A83", limits = c(0, NA),
      name = expression(epsilon^2~"(effect size)"),
      guide = ggplot2::guide_colourbar(barwidth = 9, barheight = 0.5, title.vjust = 1)
    ) +
    ggplot2::labs(
      title = "Does the DTP score itself differ across subtypes? (Kruskal-Wallis)",
      subtitle = "Epsilon-squared effect size across subtype levels; * = FDR<0.05 (BH).",
      x = NULL, y = NULL
    ) +
    .base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 9),
      panel.spacing.x = grid::unit(6, "pt")
    )
}
