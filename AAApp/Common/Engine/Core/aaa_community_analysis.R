# =============================================================================
# Triple_A community-structure analyses
# =============================================================================

aaa_ordination_plot <- function(
  scores,
  x,
  y,
  x_label,
  y_label,
  title,
  show_labels = FALSE,
  draw_ellipses = TRUE,
  subtitle_extra = NULL,
  fixed_aspect = FALSE,
  legend_title = "Treatment"
) {
  plot_scores <- as.data.frame(scores, stringsAsFactors = FALSE)

  if (!all(c(x, y) %in% names(plot_scores))) {
    stop("Ordination score columns were not found: ", x, ", ", y)
  }

  # Preserve the analytical coordinates in result tables, but separate samples
  # that would otherwise be drawn at exactly the same position. The offset is
  # deterministic and only affects the visual representation.
  x_values <- suppressWarnings(as.numeric(plot_scores[[x]]))
  y_values <- suppressWarnings(as.numeric(plot_scores[[y]]))
  x_span <- diff(range(x_values, na.rm = TRUE))
  y_span <- diff(range(y_values, na.rm = TRUE))
  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1

  coordinate_key <- ifelse(
    is.finite(x_values) & is.finite(y_values),
    paste(formatC(x_values, digits = 12, format = "fg"),
      formatC(y_values, digits = 12, format = "fg"),
      sep = "|"
    ),
    paste0("missing|", seq_along(x_values))
  )
  overlap_size <- ave(seq_along(coordinate_key), coordinate_key, FUN = length)
  overlap_index <- ave(seq_along(coordinate_key), coordinate_key, FUN = seq_along)
  overlap_groups <- unique(coordinate_key[overlap_size > 1])
  overlapping_samples <- sum(overlap_size > 1)

  angle <- 2 * pi * (overlap_index - 1) / pmax(overlap_size, 1)
  plot_scores$.plot_x <- x_values
  plot_scores$.plot_y <- y_values
  duplicated_position <- overlap_size > 1 & is.finite(x_values) & is.finite(y_values)
  plot_scores$.plot_x[duplicated_position] <-
    x_values[duplicated_position] + 0.014 * x_span * cos(angle[duplicated_position])
  plot_scores$.plot_y[duplicated_position] <-
    y_values[duplicated_position] + 0.014 * y_span * sin(angle[duplicated_position])

  subtitle <- paste0(
    "Samples displayed: ",
    sum(stats::complete.cases(plot_scores[, c(x, y), drop = FALSE])),
    " / ", nrow(plot_scores)
  )
  if (length(overlap_groups) > 0) {
    subtitle <- paste0(
      subtitle,
      " · ", overlapping_samples,
      " samples in overlapping positions were separated visually"
    )
  }
  if (!is.null(subtitle_extra) && length(subtitle_extra) > 0L &&
    any(nzchar(trimws(as.character(subtitle_extra))))) {
    subtitle <- paste(
      c(subtitle, as.character(subtitle_extra)[nzchar(trimws(as.character(subtitle_extra)))]),
      collapse = " · "
    )
  }

  # Colour represents the analytical group. Point shape is deliberately fixed:
  # mapping both colour and shape to a factor with many levels produces invalid
  # or blank symbols after ggplot2 exhausts its discrete shape palette.
  plot <- ggplot2::ggplot(
    plot_scores,
    ggplot2::aes(
      x = .plot_x,
      y = .plot_y,
      colour = Treatment
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      colour = "grey80"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      colour = "grey80"
    ) +
    ggplot2::geom_point(
      shape = 16,
      size = 3.2,
      alpha = 0.9
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      colour = legend_title
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.10)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.10)) +
    aaa_theme() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 18, 10, 10),
      legend.position = "bottom"
    )

  ellipse_data <- plot_scores[
    stats::complete.cases(
      plot_scores[, c(".plot_x", ".plot_y", "Treatment"), drop = FALSE]
    ), ,
    drop = FALSE
  ]

  eligible_groups <- vapply(
    split(ellipse_data, ellipse_data$Treatment),
    function(group_data) {
      coordinates <- unique(
        group_data[, c(".plot_x", ".plot_y"), drop = FALSE]
      )
      if (nrow(coordinates) < 3) {
        return(FALSE)
      }
      covariance <- tryCatch(stats::cov(coordinates), error = function(e) NULL)
      !is.null(covariance) && all(is.finite(covariance)) &&
        abs(det(covariance)) > sqrt(.Machine$double.eps)
    },
    logical(1)
  )

  ellipse_groups <- names(eligible_groups[eligible_groups])

  if (isTRUE(draw_ellipses) && length(ellipse_groups) > 0) {
    plot <- plot +
      ggplot2::stat_ellipse(
        data = ellipse_data[
          ellipse_data$Treatment %in% ellipse_groups, ,
          drop = FALSE
        ],
        ggplot2::aes(x = .plot_x, y = .plot_y, group = Treatment),
        linewidth = 0.65,
        alpha = 0.7,
        show.legend = FALSE,
        na.rm = TRUE
      )
  }

  if (isTRUE(fixed_aspect)) {
    plot <- plot + ggplot2::coord_fixed(clip = "off")
  }

  if (isTRUE(show_labels)) {
    plot <- plot +
      ggrepel::geom_text_repel(
        ggplot2::aes(x = .plot_x, y = .plot_y, label = Sample_column),
        size = 3,
        max.overlaps = Inf,
        seed = 1,
        show.legend = FALSE
      )
  }

  plot
}

# =============================================================================
# Optional graphical summaries derived from the community-structure tables.
# These add publication-ready figures for results that were previously exposed
# only as tables (pairwise PERMANOVA, beta dispersion, ANOSIM, PERMANOVA
# variance). They introduce no new statistical method: every value plotted is
# taken directly from the tables already computed in aaa_community_analysis().
# =============================================================================

# Uniform "analysis ran, but this figure has no data" panel so a catalogued
# figure output always exists (its methodology card and result-tree node stay
# consistent) even when, e.g., there are fewer than three groups.
aaa_community_placeholder_plot <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, size = 4.2) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    aaa_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

aaa_pairwise_permanova_plot <- function(pairwise_table, group_levels,
                                        title = "Pairwise PERMANOVA") {
  ok <- is.data.frame(pairwise_table) && nrow(pairwise_table) > 0L &&
    all(c("Group_1", "Group_2", "R2", "P_adjusted_BH") %in% names(pairwise_table)) &&
    any(is.finite(pairwise_table$P_adjusted_BH))
  if (!ok) {
    return(aaa_community_placeholder_plot(
      title,
      "Pairwise PERMANOVA requires at least three treatment groups."
    ))
  }

  levels_present <- if (length(group_levels)) {
    group_levels
  } else {
    sort(unique(c(pairwise_table$Group_1, pairwise_table$Group_2)))
  }

  mirrored <- rbind(
    data.frame(
      Row = pairwise_table$Group_1, Col = pairwise_table$Group_2,
      R2 = pairwise_table$R2, P = pairwise_table$P_adjusted_BH,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Row = pairwise_table$Group_2, Col = pairwise_table$Group_1,
      R2 = pairwise_table$R2, P = pairwise_table$P_adjusted_BH,
      stringsAsFactors = FALSE
    )
  )
  mirrored$Row <- factor(mirrored$Row, levels = levels_present)
  mirrored$Col <- factor(mirrored$Col, levels = rev(levels_present))
  mirrored$Label <- ifelse(
    is.finite(mirrored$P),
    sprintf("p=%.3f\nR²=%.2f", mirrored$P, mirrored$R2),
    ""
  )

  ident <- aaa_visual_identity()
  ggplot2::ggplot(mirrored, ggplot2::aes(x = Row, y = Col, fill = P)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = Label),
      size = 3, colour = ident$dark,
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_gradient(
      low = ident$positive, high = ident$light,
      limits = c(0, 1), na.value = "grey90", name = "Adjusted P"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "BH-adjusted permutation P-values (darker = stronger separation)",
      x = NULL, y = NULL
    ) +
    ggplot2::coord_fixed() +
    aaa_theme()
}

aaa_beta_dispersion_plot <- function(dispersion_df, title = "Beta dispersion") {
  ok <- is.data.frame(dispersion_df) && nrow(dispersion_df) > 0L &&
    all(c("Treatment", "Distance") %in% names(dispersion_df)) &&
    any(is.finite(dispersion_df$Distance))
  if (!ok) {
    return(aaa_community_placeholder_plot(
      title,
      "Beta-dispersion distances require at least two treatment groups."
    ))
  }

  ggplot2::ggplot(
    dispersion_df,
    ggplot2::aes(x = Treatment, y = Distance, colour = Treatment)
  ) +
    ggplot2::geom_boxplot(ggplot2::aes(group = Treatment),
      outlier.shape = NA,
      alpha = 0.16, show.legend = FALSE
    ) +
    ggplot2::geom_jitter(
      width = 0.12, size = 2.4, alpha = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      fun = mean, geom = "point", shape = 18, size = 3.4,
      colour = aaa_visual_identity()$dark, show.legend = FALSE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Distance to group centroid (homogeneity-of-dispersion diagnostic)",
      x = NULL, y = "Distance to centroid"
    ) +
    aaa_theme()
}

aaa_anosim_plot <- function(anosim_table, title = "ANOSIM") {
  ok <- is.data.frame(anosim_table) && nrow(anosim_table) > 0L &&
    "R" %in% names(anosim_table) && is.finite(anosim_table$R[1])
  if (!ok) {
    return(aaa_community_placeholder_plot(
      title,
      "ANOSIM requires at least two treatment groups."
    ))
  }

  R <- anosim_table$R[1]
  P <- if ("P_value" %in% names(anosim_table)) anosim_table$P_value[1] else NA_real_
  perms <- if ("Permutations" %in% names(anosim_table)) anosim_table$Permutations[1] else NA_integer_
  ident <- aaa_visual_identity()
  marker <- if (R >= 0.75) ident$positive else if (R >= 0.25) ident$warning else ident$negative

  subtitle <- if (is.finite(P)) {
    sprintf(
      "R = %.3f, P = %.3f%s", R, P,
      if (is.finite(perms)) sprintf(" (%d permutations)", as.integer(perms)) else ""
    )
  } else {
    sprintf("R = %.3f", R)
  }

  ggplot2::ggplot() +
    ggplot2::annotate("rect",
      xmin = -1, xmax = 1, ymin = 0.35, ymax = 0.65,
      fill = "grey92", colour = "grey60"
    ) +
    ggplot2::annotate("segment",
      x = 0, xend = 0, y = 0.28, yend = 0.72,
      linetype = "dashed", colour = ident$dark
    ) +
    ggplot2::annotate("point", x = R, y = 0.5, size = 9, colour = marker) +
    ggplot2::annotate("text",
      x = R, y = 0.5, label = sprintf("%.2f", R),
      colour = "white", fontface = "bold", size = 3.4
    ) +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "ANOSIM R statistic (0 = no separation, 1 = full separation)", y = NULL
    ) +
    aaa_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

aaa_permanova_variance_plot <- function(permanova_table,
                                        title = "PERMANOVA explained variance") {
  ok <- is.data.frame(permanova_table) && "R2" %in% names(permanova_table) &&
    "Term" %in% names(permanova_table) && any(is.finite(permanova_table$R2))
  if (!ok) {
    return(aaa_community_placeholder_plot(
      title,
      "PERMANOVA variance requires at least two treatment groups."
    ))
  }

  df <- permanova_table[is.finite(permanova_table$R2) &
    !permanova_table$Term %in% c("Total"), , drop = FALSE]
  if (!nrow(df)) {
    return(aaa_community_placeholder_plot(title, "No PERMANOVA variance terms available."))
  }
  df$Term <- factor(df$Term, levels = df$Term[order(df$R2)])

  ident <- aaa_visual_identity()
  ggplot2::ggplot(df, ggplot2::aes(x = Term, y = R2)) +
    ggplot2::geom_col(width = 0.62, fill = ident$primary) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", 100 * R2)),
      hjust = -0.12, size = 3.4
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      title = title,
      subtitle = "Proportion of community variation explained (R²) by model term",
      x = NULL, y = "Variance explained (R²)"
    ) +
    aaa_theme()
}

# Dendrogram drawing -----------------------------------------------------------
# The segments are derived from the hclust object itself rather than through an
# extra plotting package, so a single figure does not add a dependency the
# distribution would otherwise not need.
aaa_dendrogram_segments <- function(hc) {
  leaf_x <- numeric(length(hc$order))
  leaf_x[hc$order] <- seq_along(hc$order)
  node_x <- numeric(nrow(hc$merge))
  segments <- vector("list", nrow(hc$merge) * 3L)
  k <- 0L

  for (i in seq_len(nrow(hc$merge))) {
    # A negative entry in hc$merge is a leaf (drawn at height 0); a positive one
    # is an earlier merge, drawn at the height that merge occurred.
    corner <- lapply(hc$merge[i, ], function(m) {
      if (m < 0L) c(x = leaf_x[-m], y = 0) else c(x = node_x[m], y = hc$height[m])
    })
    height <- hc$height[i]
    node_x[i] <- mean(c(corner[[1L]][["x"]], corner[[2L]][["x"]]))

    for (side in corner) {
      k <- k + 1L
      segments[[k]] <- data.frame(
        x = side[["x"]], y = side[["y"]],
        xend = side[["x"]], yend = height
      )
    }
    k <- k + 1L
    segments[[k]] <- data.frame(
      x = corner[[1L]][["x"]], y = height,
      xend = corner[[2L]][["x"]], yend = height
    )
  }

  do.call(rbind, segments[seq_len(k)])
}

aaa_dendrogram_plot <- function(hc, sample_metadata,
                                distance_method = "bray",
                                linkage = "average",
                                title = "Hierarchical clustering") {
  if (!inherits(hc, "hclust") || length(hc$order) < 3L) {
    return(aaa_community_placeholder_plot(
      title,
      "Hierarchical clustering requires at least three samples."
    ))
  }

  segments <- aaa_dendrogram_segments(hc)
  leaves <- data.frame(
    Sample_column = hc$labels[hc$order],
    x = seq_along(hc$order),
    stringsAsFactors = FALSE
  )

  treatment <- rep(NA_character_, nrow(leaves))
  if (is.data.frame(sample_metadata) && "Treatment" %in% names(sample_metadata)) {
    matched <- match(leaves$Sample_column, rownames(sample_metadata))
    treatment <- as.character(sample_metadata$Treatment)[matched]
  }
  leaves$Treatment <- ifelse(is.na(treatment), "Unassigned", treatment)

  groups <- unique(leaves$Treatment)
  colours <- stats::setNames(aaa_treatment_palette(length(groups)), groups)
  label_offset <- 0.04 * max(hc$height, na.rm = TRUE)

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segments,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = aaa_visual_identity()$dark, linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = leaves,
      ggplot2::aes(x = x, y = 0, colour = Treatment),
      size = 2.6
    ) +
    ggplot2::geom_text(
      data = leaves,
      ggplot2::aes(x = x, y = -label_offset, label = Sample_column, colour = Treatment),
      angle = 90, hjust = 1, size = 2.8, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.28, 0.05))) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "%s distance, %s linkage",
        if (identical(distance_method, "jaccard")) "Jaccard" else "Bray-Curtis",
        linkage
      ),
      x = NULL, y = "Dissimilarity at which samples join"
    ) +
    aaa_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )
}

aaa_community_analysis <- function(
  dataset,
  abundance_type = c("proportion", "percentage", "counts"),
  transformation = c("hellinger", "relative", "log1p"),
  distance_method = c("bray", "jaccard"),
  permutations = 999, significance_alpha = 0.05,
  nmds_trymax = 50, show_sample_labels = FALSE,
  project_dir, analysis_name = "Community_structure"
) {
  abundance_type <- match.arg(abundance_type)
  transformation <- match.arg(transformation)
  distance_method <- match.arg(distance_method)

  if (!is.numeric(permutations) ||
    length(permutations) != 1 ||
    is.na(permutations) ||
    permutations < 99) {
    stop("'permutations' must be a number of at least 99.")
  }

  if (!is.numeric(significance_alpha) ||
    length(significance_alpha) != 1 ||
    is.na(significance_alpha) ||
    significance_alpha <= 0 ||
    significance_alpha >= 1) {
    stop("'significance_alpha' must be between 0 and 1.")
  }

  if (!is.numeric(nmds_trymax) ||
    length(nmds_trymax) != 1 ||
    is.na(nmds_trymax) ||
    nmds_trymax < 1) {
    stop("'nmds_trymax' must be a positive number.")
  }

  if (!is.logical(show_sample_labels) ||
    length(show_sample_labels) != 1 ||
    is.na(show_sample_labels)) {
    stop("'show_sample_labels' must be TRUE or FALSE.")
  }

  aaa_check_packages(analyses = "community_structure")

  # Shares the sample x taxon preparation pipeline (naming, zero/negative
  # clipping, relative-abundance transform, near-constant taxa removal) with
  # aaa_multivariate_input() rather than re-deriving it here.
  input <- aaa_multivariate_input(
    dataset, abundance_type, transformation, project_dir, analysis_name
  )

  if (nrow(input$relative) < 3) {
    stop(
      "At least three samples are required for two-dimensional ",
      "community ordination."
    )
  }

  transformed <- input$matrix
  relative_matrix <- input$relative
  sample_metadata <- input$metadata

  pca <- stats::prcomp(
    transformed,
    center = TRUE,
    scale. = FALSE
  )

  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  pca_scores <- data.frame(
    Sample_column = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      sample_metadata,
      by = "Sample_column"
    )

  pca_variance <- data.frame(
    Component = paste0("PC", seq_along(explained)),
    Variance_explained = explained,
    Cumulative_variance = cumsum(explained),
    stringsAsFactors = FALSE
  )

  pca_loadings <- data.frame(
    Taxon = rownames(pca$rotation),
    PC1 = pca$rotation[, 1],
    PC2 = if (ncol(pca$rotation) >= 2L) pca$rotation[, 2] else 0,
    stringsAsFactors = FALSE
  )
  pca_loadings$Contribution_PC1 <- 100 * pca_loadings$PC1^2 / sum(pca_loadings$PC1^2)
  pca_loadings$Contribution_PC2 <- 100 * pca_loadings$PC2^2 / max(sum(pca_loadings$PC2^2), .Machine$double.eps)
  pca_loadings$Cos2 <- pca_loadings$PC1^2 + pca_loadings$PC2^2
  pca_loadings <- pca_loadings[order(pca_loadings$Cos2, decreasing = TRUE), , drop = FALSE]
  pca_summary <- data.frame(
    Metric = c("Samples", "Taxa", "PC1 variance (%)", "PC2 variance (%)", "Cumulative PC1-PC2 (%)"),
    Value = c(nrow(transformed), ncol(transformed), explained[1], explained[2], sum(explained[1:2])),
    stringsAsFactors = FALSE
  )

  pca_plot <- aaa_ordination_plot(
    scores = pca_scores,
    x = "PC1",
    y = "PC2",
    x_label = sprintf("PC1 (%.1f%%)", explained[1]),
    y_label = sprintf("PC2 (%.1f%%)", explained[2]),
    title = "Principal component analysis",
    show_labels = show_sample_labels
  )

  scree_data <- utils::head(
    pca_variance,
    min(10, nrow(pca_variance))
  )

  pca_scree_plot <- ggplot2::ggplot(
    scree_data,
    ggplot2::aes(
      x = stats::reorder(Component, seq_along(Component)),
      y = Variance_explained,
      group = 1
    )
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::labs(
      title = "PCA explained variance",
      x = "Principal component",
      y = "Variance explained (%)"
    ) +
    aaa_theme()

  pca_loadings_top <- utils::head(pca_loadings, 20)
  pca_loading_plot <- ggplot2::ggplot(
    pca_loadings_top,
    ggplot2::aes(x = stats::reorder(Taxon, Cos2), y = Cos2)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "PCA taxon contributions", x = NULL, y = "Cos2 (PC1 + PC2)") +
    aaa_theme()

  beta_distance <- vegan::vegdist(
    relative_matrix,
    method = distance_method,
    binary = identical(distance_method, "jaccard")
  )

  if (!any(as.numeric(beta_distance) > 0)) {
    stop(
      "All samples have zero beta-diversity distance; ",
      "ordination cannot be calculated."
    )
  }

  # Hierarchical clustering of the same distance matrix the ordinations use, so
  # the dendrogram and the PCoA cannot disagree about how far apart samples are.
  # Average linkage (UPGMA) is the conventional choice for ecological
  # dissimilarities: it does not assume Euclidean geometry, unlike Ward.
  clustering_linkage <- "average"
  hierarchical_clustering <- stats::hclust(beta_distance, method = clustering_linkage)

  # Cut at the number of declared treatments so the table answers the question a
  # reader actually has: does an unsupervised grouping recover the design?
  declared_groups <- unique(stats::na.omit(as.character(sample_metadata$Treatment)))
  cut_k <- max(2L, min(length(declared_groups), length(hierarchical_clustering$order) - 1L))
  cluster_membership <- stats::cutree(hierarchical_clustering, k = cut_k)

  clustering_table <- data.frame(
    Sample_column = hierarchical_clustering$labels[hierarchical_clustering$order],
    stringsAsFactors = FALSE
  )
  clustering_table$Treatment <- as.character(
    sample_metadata$Treatment
  )[match(clustering_table$Sample_column, rownames(sample_metadata))]
  clustering_table$Cluster <- unname(
    cluster_membership[clustering_table$Sample_column]
  )
  clustering_table$Dendrogram_order <- seq_len(nrow(clustering_table))
  clustering_table$Linkage <- clustering_linkage
  clustering_table$Distance <- distance_method
  clustering_table$Clusters_requested <- cut_k

  pcoa <- stats::cmdscale(
    beta_distance,
    k = 2,
    eig = TRUE,
    add = TRUE
  )

  positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
  pcoa_explained <- 100 * positive_eigenvalues /
    sum(positive_eigenvalues)

  if (length(pcoa_explained) < 2) {
    pcoa_explained <- c(pcoa_explained, NA_real_)
  }

  pcoa_scores <- data.frame(
    Sample_column = rownames(pcoa$points),
    Axis1 = pcoa$points[, 1],
    Axis2 = pcoa$points[, 2],
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      sample_metadata,
      by = "Sample_column"
    )

  pcoa_plot <- aaa_ordination_plot(
    scores = pcoa_scores,
    x = "Axis1",
    y = "Axis2",
    x_label = sprintf(
      "PCoA1 (%.1f%%)",
      pcoa_explained[1]
    ),
    y_label = sprintf(
      "PCoA2 (%.1f%%)",
      pcoa_explained[2]
    ),
    title = paste0(
      "Principal coordinates analysis (",
      tools::toTitleCase(distance_method),
      ")"
    ),
    show_labels = show_sample_labels
  )

  set.seed(1)

  nmds_warnings <- character()

  nmds <- withCallingHandlers(
    vegan::metaMDS(
      relative_matrix,
      distance = distance_method,
      binary = identical(
        distance_method,
        "jaccard"
      ),
      k = 2,
      trymax = as.integer(nmds_trymax),
      autotransform = FALSE,
      trace = FALSE
    ),
    warning = function(warning_condition) {
      nmds_warnings <<- c(
        nmds_warnings,
        conditionMessage(warning_condition)
      )

      invokeRestart("muffleWarning")
    }
  )

  nmds_points <- vegan::scores(
    nmds,
    display = "sites"
  )

  nmds_scores <- data.frame(
    Sample_column = rownames(nmds_points),
    NMDS1 = nmds_points[, 1],
    NMDS2 = nmds_points[, 2],
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      sample_metadata,
      by = "Sample_column"
    )

  nmds_plot <- aaa_ordination_plot(
    scores = nmds_scores,
    x = "NMDS1",
    y = "NMDS2",
    x_label = "NMDS1",
    y_label = "NMDS2",
    title = sprintf(
      "NMDS (%s; stress = %.3f)",
      tools::toTitleCase(distance_method),
      nmds$stress
    ),
    show_labels = show_sample_labels
  )

  if (length(unique(sample_metadata$Treatment)) >= 2) {
    permanova <- vegan::adonis2(
      beta_distance ~ Treatment,
      data = sample_metadata,
      permutations = as.integer(permutations)
    )

    permanova_table <- data.frame(
      Term = rownames(permanova),
      permanova,
      row.names = NULL,
      check.names = FALSE
    )

    p_column <- intersect(
      c("Pr(>F)", "Pr..F."),
      names(permanova_table)
    )

    if (length(p_column) > 0) {
      permanova_table$Significance_alpha <- significance_alpha
      permanova_table$Significant <- ifelse(
        is.na(permanova_table[[p_column[1]]]),
        NA_character_,
        ifelse(
          permanova_table[[p_column[1]]] <= significance_alpha,
          "Yes",
          "No"
        )
      )
    }
  } else {
    permanova_table <- data.frame(
      Term = "Not performed",
      Reason = "PERMANOVA requires at least two treatments.",
      stringsAsFactors = FALSE
    )
  }


  # Community-comparison extensions. These tests reuse the same
  # validated distance matrix and grouping factor to avoid redundant work.
  anosim_table <- data.frame()
  beta_dispersion_table <- data.frame()
  # Per-sample distances to the group centroid, retained so the beta-dispersion
  # boxplot can be drawn from the same betadisper() fit that produces the test.
  beta_dispersion_df <- data.frame()
  pairwise_permanova_table <- data.frame()

  valid_groups <- !is.na(sample_metadata$Treatment)
  group_factor <- droplevels(sample_metadata$Treatment[valid_groups])
  grouped_distance <- stats::as.dist(as.matrix(beta_distance)[valid_groups, valid_groups, drop = FALSE])

  if (length(unique(group_factor)) >= 2L) {
    anosim_result <- tryCatch(
      vegan::anosim(grouped_distance, grouping = group_factor, permutations = as.integer(permutations)),
      error = function(e) e
    )
    anosim_table <- if (inherits(anosim_result, "error")) {
      data.frame(R = NA_real_, P_value = NA_real_, Permutations = as.integer(permutations), Reason = conditionMessage(anosim_result), stringsAsFactors = FALSE)
    } else {
      data.frame(R = unname(anosim_result$statistic), P_value = anosim_result$signif, Permutations = as.integer(permutations), Reason = NA_character_, stringsAsFactors = FALSE)
    }

    dispersion_result <- tryCatch(vegan::betadisper(grouped_distance, group_factor), error = function(e) e)
    if (inherits(dispersion_result, "error")) {
      beta_dispersion_table <- data.frame(F = NA_real_, P_value = NA_real_, Permutations = as.integer(permutations), Reason = conditionMessage(dispersion_result), stringsAsFactors = FALSE)
    } else {
      dispersion_test <- tryCatch(vegan::permutest(dispersion_result, permutations = as.integer(permutations)), error = function(e) e)
      if (inherits(dispersion_test, "error")) {
        beta_dispersion_table <- data.frame(F = NA_real_, P_value = NA_real_, Permutations = as.integer(permutations), Reason = conditionMessage(dispersion_test), stringsAsFactors = FALSE)
      } else {
        beta_dispersion_table <- data.frame(F = unname(dispersion_test$tab[1, "F"]), P_value = unname(dispersion_test$tab[1, "Pr(>F)"]), Permutations = as.integer(permutations), Reason = NA_character_, stringsAsFactors = FALSE)
      }
      beta_dispersion_df <- data.frame(
        Sample_column = names(dispersion_result$distances),
        Treatment = as.character(dispersion_result$group),
        Distance = as.numeric(dispersion_result$distances),
        stringsAsFactors = FALSE
      )
    }

    levels_present <- levels(group_factor)
    if (length(levels_present) >= 3L) {
      comparisons <- utils::combn(levels_present, 2L, simplify = FALSE)
      pair_rows <- lapply(comparisons, function(pair) {
        keep <- group_factor %in% pair
        pair_group <- droplevels(group_factor[keep])
        pair_dist <- stats::as.dist(as.matrix(grouped_distance)[keep, keep, drop = FALSE])
        fit <- tryCatch(vegan::adonis2(pair_dist ~ pair_group, permutations = as.integer(permutations)), error = function(e) e)
        if (inherits(fit, "error")) {
          return(data.frame(Group_1 = pair[1], Group_2 = pair[2], F = NA_real_, R2 = NA_real_, P_value = NA_real_, Reason = conditionMessage(fit), stringsAsFactors = FALSE))
        }
        data.frame(Group_1 = pair[1], Group_2 = pair[2], F = unname(fit$F[1]), R2 = unname(fit$R2[1]), P_value = unname(fit$`Pr(>F)`[1]), Reason = NA_character_, stringsAsFactors = FALSE)
      })
      pairwise_permanova_table <- do.call(rbind, pair_rows)
      pairwise_permanova_table$P_adjusted_BH <- stats::p.adjust(pairwise_permanova_table$P_value, method = "BH")
    }
  }

  alpha_table <- data.frame(
    Sample_column = rownames(relative_matrix),
    Observed_taxa = rowSums(relative_matrix > 0),
    Shannon = vegan::diversity(
      relative_matrix,
      index = "shannon"
    ),
    Simpson = vegan::diversity(
      relative_matrix,
      index = "simpson"
    ),
    Inverse_Simpson = vegan::diversity(
      relative_matrix,
      index = "invsimpson"
    ),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      sample_metadata,
      by = "Sample_column"
    )

  alpha_long <- alpha_table |>
    tidyr::pivot_longer(
      cols = c(
        Observed_taxa,
        Shannon,
        Simpson,
        Inverse_Simpson
      ),
      names_to = "Metric",
      values_to = "Value"
    )

  alpha_plot <- ggplot2::ggplot(
    alpha_long,
    ggplot2::aes(
      x = Treatment,
      y = Value,
      colour = Treatment
    )
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = Treatment),
      outlier.shape = NA,
      alpha = 0.16,
      show.legend = FALSE
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      size = 2.5,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ~Metric,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::labs(
      title = "Alpha diversity",
      x = NULL,
      y = NULL
    ) +
    aaa_theme()

  beta_matrix <- as.matrix(beta_distance)
  beta_table <- data.frame(
    Sample_column = rownames(beta_matrix),
    beta_matrix,
    check.names = FALSE,
    row.names = NULL
  )

  project <- aaa_create_project_structure(
    project_dir,
    analysis_name
  )

  # Optional graphical summaries of the community-comparison tables. Values are
  # taken directly from the tables computed above; no new statistics are run.
  pairwise_permanova_plot <- aaa_pairwise_permanova_plot(
    pairwise_permanova_table, levels(group_factor)
  )
  beta_dispersion_plot <- aaa_beta_dispersion_plot(beta_dispersion_df)
  anosim_plot <- aaa_anosim_plot(anosim_table)
  permanova_variance_plot <- aaa_permanova_variance_plot(permanova_table)
  dendrogram_plot <- aaa_dendrogram_plot(
    hierarchical_clustering,
    sample_metadata,
    distance_method = distance_method,
    linkage = clustering_linkage
  )

  files <- c(
    dendrogram = file.path(project$analysis, "Hierarchical_clustering.png"),
    pca = file.path(project$analysis, "PCA.png"),
    pca_scree = file.path(project$analysis, "PCA_scree.png"),
    pca_loadings = file.path(project$analysis, "PCA_loadings.png"),
    pcoa = file.path(project$analysis, "PCoA.png"),
    nmds = file.path(project$analysis, "NMDS.png"),
    alpha_diversity = file.path(
      project$analysis,
      "Alpha_diversity.png"
    ),
    pairwise_permanova = file.path(project$analysis, "Pairwise_PERMANOVA_heatmap.png"),
    beta_dispersion = file.path(project$analysis, "Beta_dispersion_boxplot.png"),
    anosim = file.path(project$analysis, "ANOSIM_R.png"),
    permanova_variance = file.path(project$analysis, "PERMANOVA_variance.png"),
    summary = file.path(
      project$analysis,
      "Community_structure_summary.xlsx"
    )
  )

  aaa_save_plot(pca_plot, files[["pca"]], width = 8, height = 6)
  aaa_save_plot(
    pca_scree_plot,
    files[["pca_scree"]],
    width = 8,
    height = 5.5
  )
  aaa_save_plot(
    pca_loading_plot, files[["pca_loadings"]],
    width = max(9, aaa_flipped_axis_plot_width(pca_loadings_top$Taxon)), height = 7
  )
  aaa_save_plot(pcoa_plot, files[["pcoa"]], width = 8, height = 6)
  aaa_save_plot(nmds_plot, files[["nmds"]], width = 8, height = 6)
  aaa_save_plot(
    alpha_plot,
    files[["alpha_diversity"]],
    width = 10,
    height = 7
  )
  aaa_save_plot(
    pairwise_permanova_plot, files[["pairwise_permanova"]],
    width = max(7, 1.1 * length(levels(group_factor)) + 3.5),
    height = max(6, 1.0 * length(levels(group_factor)) + 3)
  )
  aaa_save_plot(beta_dispersion_plot, files[["beta_dispersion"]], width = 8, height = 6)
  aaa_save_plot(anosim_plot, files[["anosim"]], width = 8, height = 5)
  # Leaf labels are rotated 90 degrees, so they consume height, not width. Width
  # has to grow with the number of samples instead, or the leaves overlap.
  aaa_save_plot(
    dendrogram_plot,
    files[["dendrogram"]],
    width = max(8, min(24, 0.38 * nrow(clustering_table))),
    height = 6.5
  )
  aaa_save_plot(permanova_variance_plot, files[["permanova_variance"]], width = 8, height = 5.5)

  openxlsx::write.xlsx(
    list(
      Sample_metadata = sample_metadata,
      Summary = pca_summary,
      PCA_scores = pca_scores,
      PCA_variance = pca_variance,
      PCA_loadings = pca_loadings,
      PCoA_scores = pcoa_scores,
      NMDS_scores = nmds_scores,
      PERMANOVA = permanova_table,
      Pairwise_PERMANOVA = pairwise_permanova_table,
      ANOSIM = anosim_table,
      Beta_dispersion = beta_dispersion_table,
      Alpha_diversity = alpha_table,
      Beta_distance = beta_table
    ),
    files[["summary"]],
    overwrite = TRUE
  )

  aaa_result(
    tables = list(
      summary = pca_summary,
      sample_metadata = sample_metadata,
      pca_scores = pca_scores,
      pca_variance = pca_variance,
      pca_loadings = pca_loadings,
      pcoa_scores = pcoa_scores,
      nmds_scores = nmds_scores,
      nmds_diagnostics = data.frame(
        Stress = nmds$stress,
        Warning = if (length(nmds_warnings) == 0) {
          NA_character_
        } else {
          paste(
            unique(nmds_warnings),
            collapse = " | "
          )
        },
        stringsAsFactors = FALSE
      ),
      permanova = permanova_table,
      pairwise_permanova = pairwise_permanova_table,
      anosim = anosim_table,
      beta_dispersion = beta_dispersion_table,
      alpha_diversity = alpha_table,
      beta_distance = beta_table,
      hierarchical_clustering = clustering_table
    ),
    plots = list(
      dendrogram = dendrogram_plot,
      pca = pca_plot,
      pca_scree = pca_scree_plot,
      pca_loadings = pca_loading_plot,
      pcoa = pcoa_plot,
      nmds = nmds_plot,
      alpha_diversity = alpha_plot,
      pairwise_permanova = pairwise_permanova_plot,
      beta_dispersion = beta_dispersion_plot,
      anosim = anosim_plot,
      permanova_variance = permanova_variance_plot
    ),
    files = files,
    output_dir = project$analysis,
    metadata = list(
      transformation = transformation,
      distance_method = distance_method,
      permutations = as.integer(permutations),
      significance_alpha = significance_alpha,
      nmds_stress = nmds$stress,
      clustering_linkage = clustering_linkage,
      clusters_requested = cut_k,
      n_samples = nrow(relative_matrix),
      n_taxa = ncol(relative_matrix)
    )
  )
}
