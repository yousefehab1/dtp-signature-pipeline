# ==============================================================================
# R/plotting.R - Core plotting themes, palettes, dataset/endpoint labels,
# and figure output utilities.
# ==============================================================================

#' Publication theme for ggplot2
#'
#' Standardized classic theme with bold titles and bottom legend.
#' @export
pub_theme <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = 12),
    axis.title    = ggplot2::element_text(face = "bold"),
    axis.text     = ggplot2::element_text(colour = "black"),
    legend.position = "bottom",
    legend.title  = ggplot2::element_text(face = "bold"),
    strip.text    = ggplot2::element_text(face = "bold")
  )

#' Standardized color palette for outcomes, groups, and sample types
#' @export
pub_palette <- c(
  "Dead_3yr"        = "#E64B35", "Alive_3yr"       = "#4DBBD5",
  "Recurred"        = "#E64B35", "Recurrence-Free" = "#4DBBD5",
  "High"            = "#E64B35", "Low"             = "#4DBBD5",
  "Normal"          = "#55A868", "Primary"         = "#4C72B0", "Metastasis" = "#C44E52"
)

#' Canonical full dataset display labels
#' @export
DATASET_FULL <- c(
  "GSE39582"           = "Marisa (GSE39582)",
  "GSE39582 (treated)" = "Marisa treated",
  "TCGA-COAD"          = "TCGA-COAD"
)

#' Canonical short dataset display labels
#' @export
DATASET_SHORT <- c(
  "GSE39582"           = "Marisa",
  "GSE39582 (treated)" = "Marisa treated",
  "TCGA-COAD"          = "TCGA"
)

#' Canonical survival endpoint display labels
#' @export
ENDPOINT_LABEL <- c(
  OS  = "3-yr OS",
  RFS = "3-yr RFS"
)

# Significance stars helper: returns asterisks for p < 0.05 / 0.01 / 0.001. Internal.
.sig_stars <- function(p, ns = "") {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ns))))
}

# Shared composite base theme. Internal.
.base_theme <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey30", hjust = 0),
    strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
    strip.text       = ggplot2::element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    plot.tag         = ggplot2::element_text(face = "bold", size = 16)
  )

#' Save a figure as both vector PDF and 300-dpi PNG
#'
#' Writes `<name>.pdf` via [grDevices::pdf()] and `<name>.png` via [ragg::agg_png()]
#' if installed (falling back to [grDevices::png()]).
#'
#' @param plot A ggplot or patchwork object to save.
#' @param name Base file name (without extension).
#' @param out_dir Output directory.
#' @param w Figure width in inches.
#' @param h Figure height in inches.
#' @export
save_fig <- function(plot, name, out_dir, w, h) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))

  ggplot2::ggsave(
    pdf_path, plot, width = w, height = h,
    device = grDevices::pdf, limitsize = FALSE
  )

  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png
  ggplot2::ggsave(
    png_path, plot, width = w, height = h,
    dpi = 300, bg = "white", limitsize = FALSE, device = png_dev
  )

  message("  wrote ", name, " (", w, "x", h, " in)")
}
