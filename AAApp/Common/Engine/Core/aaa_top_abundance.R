#' Analyze the most abundant taxa
#'
#' Generates a heatmap, lollipop plot, stacked composition plot and
#' replicate-level distribution plot using the shared Advanced_Amplicon_Analysis style.
#'
#' @return A Triple_A_result object.
aaa_top_abundance <- function(
  dataset, top_n = 20,
  abundance_type = c("proportion", "percentage", "counts"),
  project_dir, analysis_name = "aaa_top_abundance",
  graph_main = "Most abundant microorganisms"
) {
  abundance_type <- match.arg(abundance_type)
  aaa_check_packages(analyses = "top_abundance")

  if (!is.numeric(top_n) || length(top_n) != 1 ||
    is.na(top_n) || top_n < 1 || top_n %% 1 != 0) {
    stop("'top_n' must be a positive integer.")
  }

  prepared <- aaa_prepare_amplicon_data(
    dataset = dataset, abundance_type = abundance_type,
    project_dir = project_dir, analysis_name = analysis_name
  )

  ranked <- prepared$wide
  ranked$Mean_abundance <- rowMeans(
    ranked[prepared$sample_columns],
    na.rm = TRUE
  )
  ranked <- ranked |>
    dplyr::arrange(dplyr::desc(Mean_abundance)) |>
    dplyr::slice_head(n = min(top_n, nrow(ranked))) |>
    dplyr::mutate(
      # Taxonomy/Genus are coerced to character before nzchar(): if either
      # column is read in as a factor (which can happen depending on the
      # source file/import path), nzchar() on a raw factor fails with
      # "'nzchar()' requires a character vector" instead of just working.
      Taxon = make.unique(
        ifelse(
          !is.na(Taxonomy) & nzchar(as.character(Taxonomy)),
          as.character(Taxonomy), as.character(Genus)
        )
      )
    )

  matrix <- ranked |>
    dplyr::select(Taxon, dplyr::all_of(prepared$sample_columns)) |>
    tibble::column_to_rownames("Taxon") |>
    as.data.frame()

  if (isTRUE(prepared$has_replicates)) {
    summary <- aaa_replicate_summary(
      matrix, prepared$sample_map, prepared$samples_name
    )
    heatmap_data <- summary$mean
    heatmap_labels <- summary$labels
  } else {
    heatmap_data <- matrix
    names(heatmap_data) <- prepared$samples_name
    heatmap_labels <- matrix(
      sprintf("%.1f%%", as.matrix(heatmap_data)),
      nrow = nrow(heatmap_data), dimnames = dimnames(heatmap_data)
    )
  }

  long <- ranked |>
    dplyr::select(Taxon, dplyr::all_of(prepared$sample_columns)) |>
    tidyr::pivot_longer(
      cols = -Taxon,
      names_to = "Sample_column",
      values_to = "Abundance"
    ) |>
    dplyr::left_join(prepared$sample_map, by = "Sample_column")

  taxon_levels <- rev(ranked$Taxon)
  long$Taxon <- factor(long$Taxon, levels = taxon_levels)
  long$Treatment <- factor(long$Treatment, levels = prepared$samples_name)

  heatmap_file <- file.path(
    prepared$project$analysis,
    "Top_taxa_heatmap.png"
  )
  aaa_plot_heatmap(
    heatmap_data,
    filename = heatmap_file,
    title = graph_main,
    labels = heatmap_labels,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )

  overall <- long |>
    dplyr::group_by(Taxon) |>
    dplyr::summarise(
      Mean = mean(Abundance, na.rm = TRUE),
      SD = stats::sd(Abundance, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(Mean) |>
    dplyr::mutate(Taxon = factor(Taxon, levels = Taxon))

  accent <- aaa_colours()[["accent"]]

  lollipop_plot <- ggplot2::ggplot(
    overall,
    ggplot2::aes(x = Taxon, y = Mean)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = Taxon, y = 0, yend = Mean),
      colour = "grey75",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      colour = accent,
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(Mean - SD, 0),
        ymax = Mean + SD
      ),
      width = 0.15,
      colour = "grey35"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = graph_main,
      subtitle = "Global mean abundance ± standard deviation",
      x = NULL,
      y = "Mean relative abundance (%)"
    ) +
    aaa_theme()

  composition <- long |>
    dplyr::group_by(Treatment, Taxon) |>
    dplyr::summarise(
      Abundance = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    )

  taxon_colours <- stats::setNames(
    aaa_palette(length(unique(composition$Taxon))),
    unique(composition$Taxon)
  )

  stacked_plot <- ggplot2::ggplot(
    composition,
    ggplot2::aes(
      x = Treatment,
      y = Abundance,
      fill = Taxon
    )
  ) +
    ggplot2::geom_col(
      width = 0.75,
      colour = "white",
      linewidth = 0.15
    ) +
    ggplot2::scale_fill_manual(
      values = taxon_colours,
      labels = aaa_wrap_legend_labels(
        names(taxon_colours),
        width = 32
      ),
      guide = ggplot2::guide_legend(
        ncol = aaa_legend_columns(
          length(taxon_colours),
          max_per_row = 4
        ),
        byrow = TRUE,
        title.position = "top"
      )
    ) +
    ggplot2::labs(
      title = paste0(graph_main, ": community composition"),
      subtitle = "Only the selected top taxa are represented",
      x = NULL,
      y = "Mean relative abundance (%)",
      fill = "Taxon"
    ) +
    aaa_theme() +
    aaa_legend_theme(
      n_items = length(taxon_colours),
      position = "bottom",
      max_per_row = 4,
      text_size = 8.5
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
    )

  distribution_plot <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = Taxon,
      y = Abundance,
      colour = Treatment
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      colour = "grey45",
      fill = "white",
      linewidth = 0.45
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      alpha = 0.8,
      size = 1.8
    ) +
    ggplot2::scale_colour_manual(
      values = aaa_palette(length(prepared$samples_name)),
      labels = aaa_wrap_legend_labels(
        prepared$samples_name,
        width = 28
      ),
      guide = ggplot2::guide_legend(
        ncol = aaa_legend_columns(
          length(prepared$samples_name),
          max_per_row = 4
        ),
        byrow = TRUE
      )
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0(graph_main, ": replicate distribution"),
      x = NULL,
      y = "Relative abundance (%)",
      colour = "Treatment"
    ) +
    aaa_theme() +
    aaa_legend_theme(
      n_items = length(prepared$samples_name),
      position = "bottom",
      max_per_row = 4
    )

  lollipop_file <- file.path(
    prepared$project$analysis,
    "Top_taxa_lollipop.png"
  )
  stacked_file <- file.path(
    prepared$project$analysis,
    "Top_taxa_stacked_composition.png"
  )
  distribution_file <- file.path(
    prepared$project$analysis,
    "Top_taxa_distribution.png"
  )

  aaa_save_plot(
    lollipop_plot,
    lollipop_file,
    width = max(8, aaa_flipped_axis_plot_width(overall$Taxon)),
    height = max(5, 0.28 * nrow(ranked) + 2)
  )
  aaa_save_plot(
    stacked_plot,
    stacked_file,
    width = max(
      9,
      1.15 * length(prepared$samples_name) + 5
    ),
    height = aaa_plot_height_with_legend(
      base_height = 5.5,
      n_items = length(taxon_colours),
      max_per_row = 4,
      row_height = 0.55
    )
  )
  aaa_save_plot(
    distribution_plot,
    distribution_file,
    width = max(9, aaa_flipped_axis_plot_width(overall$Taxon)),
    height = max(
      5.5,
      0.30 * nrow(ranked) +
        aaa_plot_height_with_legend(
          base_height = 1.5,
          n_items = length(prepared$samples_name),
          max_per_row = 4,
          row_height = 0.4
        )
    )
  )

  summary_file <- file.path(
    prepared$project$analysis,
    "aaa_top_abundance_summary.xlsx"
  )
  openxlsx::write.xlsx(
    list(
      Top_taxa = ranked,
      Heatmap_values = tibble::rownames_to_column(
        heatmap_data,
        "Taxon"
      ),
      Long_format = long
    ),
    summary_file,
    overwrite = TRUE
  )

  aaa_result(
    tables = list(
      top_taxa = ranked,
      heatmap_values = heatmap_data,
      long = long
    ),
    plots = list(
      heatmap = heatmap_file,
      lollipop = lollipop_plot,
      stacked_composition = stacked_plot,
      distribution = distribution_plot
    ),
    files = c(
      summary = summary_file,
      heatmap = heatmap_file,
      lollipop = lollipop_file,
      stacked_composition = stacked_file,
      distribution = distribution_file
    ),
    output_dir = prepared$project$analysis,
    metadata = list(
      top_n = top_n,
      abundance_type = abundance_type
    )
  )
}
