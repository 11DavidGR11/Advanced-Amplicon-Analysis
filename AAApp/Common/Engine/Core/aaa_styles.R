TRIPLE_A_VISUAL_IDENTITY <- list(
  primary = "#5B3F7A",
  secondary = "#8C6BB1",
  accent = "#2A9D8F",
  positive = "#B23A48",
  negative = "#2878B5",
  warning = "#D89B2B",
  neutral = "#B7B7B7",
  dark = "#2E2E38",
  light = "#F6F3F8",
  background = "#FFFFFF",
  grid = "#E4DDE8",
  text = "#2E2E38",
  font_family = "sans"
)

aaa_visual_identity <- function() {
  TRIPLE_A_VISUAL_IDENTITY
}

aaa_treatment_palette <- function(n) {
  identity <- aaa_visual_identity()

  base_palette <- c(
    identity$primary,
    identity$accent,
    identity$secondary,
    identity$warning,
    identity$negative,
    identity$positive,
    "#4C956C",
    "#D67AB1",
    "#6C8EBF",
    "#A67C52"
  )

  if (n <= length(base_palette)) {
    return(base_palette[seq_len(n)])
  }

  grDevices::colorRampPalette(
    base_palette
  )(n)
}

# =============================================================================
# Shared visual system for Advanced_Amplicon_Analysis
# Inspired by the consistent publication-oriented style used in
# Metabolomic_Analysis.Rmd.
# =============================================================================

aaa_colours <- function() {
  c(
    increase = "#F2A51A",
    decrease = "#26A7E8",
    neutral = "#B8B8B8",
    dark = "#303030",
    accent = "#6A4C93"
  )
}

aaa_palette <- function(n) {
  grDevices::hcl.colors(n, palette = "Dark 3")
}

aaa_heatmap_palette <- function(n = 100) {
  grDevices::colorRampPalette(
    c("#FFFFFF", "#FFF1B8", "#F2A51A", "#7A3E9D")
  )(n)
}

aaa_theme <- function(base_size = 12,
                      legend_position = "bottom") {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = base_size + 2,
        margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, colour = "grey35",
        margin = ggplot2::margin(b = 10)
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.box.margin = ggplot2::margin(t = 4, r = 4, b = 4, l = 4),
      plot.margin = ggplot2::margin(t = 10, r = 12, b = 10, l = 10),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey92", linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        fill = "grey95", colour = NA
      ),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

aaa_legend_rows <- function(n_items, max_per_row = 4L) {
  if (is.null(n_items) || is.na(n_items) || n_items < 1) {
    return(1L)
  }
  as.integer(ceiling(n_items / max(1L, max_per_row)))
}

aaa_legend_columns <- function(n_items, max_per_row = 4L) {
  if (is.null(n_items) || is.na(n_items) || n_items < 1) {
    return(1L)
  }
  as.integer(min(n_items, max(1L, max_per_row)))
}

aaa_wrap_legend_labels <- function(labels, width = 34L) {
  vapply(
    as.character(labels),
    function(label) paste(strwrap(label, width = width), collapse = "\n"),
    character(1)
  )
}

aaa_legend_theme <- function(
  n_items,
  position = "bottom",
  max_per_row = 4L,
  text_size = 9
) {
  rows <- aaa_legend_rows(n_items, max_per_row)

  ggplot2::theme(
    legend.position = position,
    legend.box = if (position == "bottom") "vertical" else "horizontal",
    legend.box.just = "center",
    legend.text = ggplot2::element_text(
      size = text_size,
      lineheight = 0.95
    ),
    legend.key.height = grid::unit(0.45, "cm"),
    legend.key.width = grid::unit(0.55, "cm"),
    plot.margin = ggplot2::margin(
      t = 10,
      r = 14,
      b = 10 + 8 * rows,
      l = 10
    )
  )
}

aaa_plot_height_with_legend <- function(
  base_height = 6,
  n_items = 0L,
  max_per_row = 4L,
  row_height = 0.42
) {
  base_height +
    row_height * aaa_legend_rows(n_items, max_per_row)
}

aaa_save_plot <- function(plot, filename,
                          width = 8, height = 6,
                          dpi = 300) {
  # ggplot2 >= 3.5/4.0 no longer creates the destination directory implicitly
  # and errors if it is missing; create.dir = TRUE restores that behaviour so a
  # plot whose analysis folder has not been pre-created still saves.
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE,
    create.dir = TRUE
  )
  invisible(filename)
}

#' Measure the rendered width of the widest string, in inches, at a given
#' font size. Used to size heatmaps so long row labels (e.g. full taxonomic
#' lineages) never get clipped by a fixed image width.
aaa_measure_text_width_inches <- function(labels, fontsize) {
  # Callers routinely pass factor labels (e.g. a ggplot Taxon column ordered
  # via factor() for plotting order); nzchar() rejects factors outright with
  # "requires a character vector", so coerce before filtering.
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (!length(labels)) {
    return(0)
  }
  grDevices::png(tempfile(fileext = ".png"), width = 4, height = 4, units = "in", res = 96)
  on.exit(grDevices::dev.off())
  graphics::par(ps = fontsize)
  max(graphics::strwidth(labels, units = "inches", cex = 1, font = 1))
}

#' Width (inches) a coord_flip() categorical-axis plot needs so long labels
#' (taxon names, full taxonomic lineages) are not clipped by a fixed image
#' width; mirrors aaa_plot_heatmap()'s dynamic row-label sizing below.
aaa_flipped_axis_plot_width <- function(labels, base_width = 5, fontsize = 10) {
  base_width + aaa_measure_text_width_inches(labels, fontsize)
}

aaa_plot_heatmap <- function(
  matrix, filename, title,
  labels = NULL,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  width = NULL,
  height = NULL,
  fontsize = 11
) {
  matrix <- as.matrix(matrix)

  if (is.null(height)) {
    height <- max(3.5, 0.42 * nrow(matrix) + 1.5)
  }

  if (is.null(width)) {
    # A fixed image width clips long row labels (full taxonomic lineages
    # routinely exceed 100 characters) instead of the plot growing to fit
    # them, so measure the actual rendered label width and size around it.
    row_label_inches <- aaa_measure_text_width_inches(rownames(matrix), fontsize)
    body_inches <- max(3, 0.55 * ncol(matrix))
    width <- body_inches + row_label_inches + 1.6
  }

  finite_values <- as.numeric(matrix[is.finite(matrix)])

  if (length(finite_values) == 0L) {
    plot_matrix <- matrix
    plot_matrix[, ] <- 0
    heatmap_breaks <- seq(-0.5, 0.5, length.out = 101L)
  } else {
    plot_matrix <- matrix
    value_range <- range(finite_values, na.rm = TRUE)

    if (!all(is.finite(value_range)) || diff(value_range) == 0) {
      centre <- if (all(is.finite(value_range))) value_range[[1L]] else 0
      padding <- max(abs(centre) * 0.01, 0.5)
      value_range <- c(centre - padding, centre + padding)
    }

    heatmap_breaks <- seq(
      from = value_range[[1L]],
      to = value_range[[2L]],
      length.out = 101L
    )
  }

  heatmap_breaks <- unique(as.numeric(heatmap_breaks))
  if (length(heatmap_breaks) < 2L) {
    heatmap_breaks <- c(-0.5, 0.5)
  }

  pheatmap::pheatmap(
    plot_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    main = title,
    angle_col = 315,
    # pheatmap's own display_numbers handling does a raw if(x) internally
    # and errors on NULL ("argument is of length zero"); this function's
    # own labels = NULL default must not hit that.
    display_numbers = if (is.null(labels)) FALSE else labels,
    number_color = "black",
    fontsize_number = 8,
    fontsize = fontsize,
    fontsize_col = fontsize,
    border_color = "white",
    color = aaa_heatmap_palette(length(heatmap_breaks) - 1L),
    breaks = heatmap_breaks,
    filename = filename,
    width = width,
    height = height,
    show_rownames = TRUE
  )

  invisible(filename)
}
