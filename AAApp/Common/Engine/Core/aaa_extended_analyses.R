# Shared by the volcano and MA plots below: both colour points by the same
# colour_by/point_size choice and only differ in title/axis labels.
aaa_apply_differential_colour_scale <- function(plot, colour_by, point_size, group_a, group_b) {
  plot <- plot +
    ggplot2::labs(
      colour = switch(colour_by,
        significance = "Classification",
        log2FC = "log2FC",
        abundance = "Mean abundance (%)"
      )
    )

  if (identical(point_size, "mean_abundance")) {
    plot <- plot + ggplot2::labs(size = "Mean abundance (%)")
  }

  if (colour_by == "significance") {
    plot <- plot +
      ggplot2::scale_colour_manual(
        values = stats::setNames(
          c(
            aaa_visual_identity()$positive,
            aaa_visual_identity()$negative,
            aaa_visual_identity()$neutral
          ),
          c(
            paste0("Higher in ", group_b),
            paste0("Higher in ", group_a),
            "Not significant"
          )
        ),
        drop = FALSE
      )
  } else if (colour_by == "log2FC") {
    plot <- plot +
      ggplot2::scale_colour_gradient2(
        low = aaa_visual_identity()$negative,
        mid = aaa_visual_identity()$neutral,
        high = aaa_visual_identity()$positive,
        midpoint = 0
      )
  } else {
    plot <- plot +
      ggplot2::scale_colour_gradient(
        low = aaa_visual_identity()$light,
        high = aaa_visual_identity()$primary
      )
  }

  plot
}

#' Pairwise differential-abundance analysis
#'
#' Performs all pairwise treatment comparisons and exports one volcano plot
#' per comparison. This is an exploratory analysis of relative abundances.
#'
#' @return A Triple_A_result object.
aaa_differential_abundance <- function(
  dataset,
  abundance_type = c("proportion", "percentage", "counts"),
  method = c("wilcox", "t_test"), paired = FALSE,
  pseudocount = 1e-06, min_prevalence = 0,
  min_mean_abundance = 0, alpha = 0.05,
  log2fc_threshold = 1, top_n_table = 25,
  max_labels = 10, label_only_significant = TRUE,
  colour_by = c("significance", "log2FC", "abundance"),
  point_size = c("constant", "mean_abundance"),
  x_limit = NULL, dynamic_ylim = TRUE, comparisons = NULL,
  project_dir,
  analysis_name = "aaa_differential_abundance"
) {
  abundance_type <- match.arg(abundance_type)
  method <- match.arg(method)
  colour_by <- match.arg(colour_by)
  point_size <- match.arg(point_size)
  aaa_check_packages(analyses = "differential_abundance")

  if (!is.logical(paired) || length(paired) != 1 || is.na(paired)) {
    stop("'paired' must be TRUE or FALSE.")
  }
  if (!is.numeric(pseudocount) || length(pseudocount) != 1 ||
    is.na(pseudocount) || pseudocount <= 0) {
    stop("'pseudocount' must be positive.")
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 ||
    is.na(min_prevalence) || min_prevalence < 0 || min_prevalence > 1) {
    stop("'min_prevalence' must be between 0 and 1.")
  }
  if (!is.numeric(min_mean_abundance) ||
    length(min_mean_abundance) != 1 ||
    is.na(min_mean_abundance) || min_mean_abundance < 0) {
    stop("'min_mean_abundance' must be non-negative.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 ||
    is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1.")
  }
  if (!is.numeric(log2fc_threshold) ||
    length(log2fc_threshold) != 1 ||
    is.na(log2fc_threshold) || log2fc_threshold < 0) {
    stop("'log2fc_threshold' must be non-negative.")
  }
  if (!is.numeric(top_n_table) || length(top_n_table) != 1 ||
    is.na(top_n_table) || top_n_table < 1) {
    stop("'top_n_table' must be a positive integer.")
  }
  if (!is.numeric(max_labels) || length(max_labels) != 1 ||
    is.na(max_labels) || max_labels < 0) {
    stop("'max_labels' must be non-negative.")
  }
  if (!is.logical(label_only_significant) ||
    length(label_only_significant) != 1 ||
    is.na(label_only_significant)) {
    stop("'label_only_significant' must be TRUE or FALSE.")
  }
  if (!is.null(x_limit) &&
    (!is.numeric(x_limit) || length(x_limit) != 1 ||
      is.na(x_limit) || x_limit <= 0)) {
    stop("'x_limit' must be NULL or a positive number.")
  }
  if (!is.logical(dynamic_ylim) || length(dynamic_ylim) != 1 ||
    is.na(dynamic_ylim)) {
    stop("'dynamic_ylim' must be TRUE or FALSE.")
  }

  prepared <- aaa_prepare_amplicon_data(
    dataset = dataset, abundance_type = abundance_type,
    project_dir = project_dir, analysis_name = analysis_name
  )
  if (length(prepared$samples_name) < 2) stop("At least two treatments are required.")
  if (prepared$n_replicates < 2L) {
    stop(
      "Differential abundance requires biological replication: at least two samples per treatment are required.",
      call. = FALSE
    )
  }

  available_comparisons <- utils::combn(prepared$samples_name, 2, simplify = FALSE)
  if (is.null(comparisons) || length(comparisons) == 0L) {
    comparisons <- available_comparisons
  } else {
    requested <- lapply(as.character(comparisons), function(x) strsplit(x, "\\|\\|", perl = TRUE)[[1]])
    requested <- requested[vapply(requested, length, integer(1)) == 2L]
    comparisons <- Filter(function(pair) all(pair %in% prepared$samples_name), requested)
    if (!length(comparisons)) stop("None of the selected pairwise comparisons matches the current treatment groups.")
  }

  all_results <- list()
  top_results <- list()
  volcano_plots <- list()
  ma_plots <- list()
  qq_plots <- list()
  volcano_files <- character()
  ma_files <- character()
  qq_files <- character()

  for (pair in comparisons) {
    a <- pair[1]
    b <- pair[2]
    comparison <- paste0(b, "_vs_", a)

    data <- prepared$long |>
      dplyr::filter(Treatment %in% c(a, b))

    split_taxa <- split(
      data,
      interaction(data$Taxonomy, data$Genus,
        data$Tax_level,
        drop = TRUE
      )
    )

    result <- dplyr::bind_rows(lapply(split_taxa, function(df) {
      va <- df$Abundance[df$Treatment == a]
      vb <- df$Abundance[df$Treatment == b]

      mean_a <- mean(va, na.rm = TRUE)
      mean_b <- mean(vb, na.rm = TRUE)
      mean_abundance <- mean(c(va, vb), na.rm = TRUE)
      prevalence_a <- mean(va > 0, na.rm = TRUE)
      prevalence_b <- mean(vb > 0, na.rm = TRUE)

      tested <- max(prevalence_a, prevalence_b) >= min_prevalence &&
        max(mean_a, mean_b) >= min_mean_abundance

      p_value <- NA_real_

      if (tested) {
        test <- tryCatch(
          if (method == "wilcox") {
            stats::wilcox.test(
              vb, va,
              paired = paired, exact = FALSE
            )
          } else {
            stats::t.test(vb, va, paired = paired)
          },
          error = function(e) NULL
        )
        if (!is.null(test)) p_value <- test$p.value
      }

      # A plain list is bound just as well as a one-row tibble by
      # dplyr::bind_rows() below, but avoids tibble::tibble()'s per-call
      # quasiquotation/glue overhead. Profiling showed tibble::tibble()
      # alone accounted for ~48% of this function's total runtime when
      # called once per taxon per pairwise comparison.
      list(
        Taxonomy = df$Taxonomy[1],
        Genus = df$Genus[1],
        Tax_level = df$Tax_level[1],
        Treatment_A = a,
        Treatment_B = b,
        Mean_A = mean_a,
        Mean_B = mean_b,
        Mean_abundance = mean_abundance,
        Difference = mean_b - mean_a,
        log2FC = log2(
          (mean_b + pseudocount) /
            (mean_a + pseudocount)
        ),
        Prevalence_A = prevalence_a,
        Prevalence_B = prevalence_b,
        P_value = p_value,
        Tested = tested
      )
    }))

    result <- result |>
      dplyr::mutate(
        Adjusted_P = stats::p.adjust(P_value, method = "BH"),
        Significance = dplyr::case_when(
          !is.na(Adjusted_P) & Adjusted_P < alpha &
            log2FC >= log2fc_threshold ~ paste0("Higher in ", b),
          !is.na(Adjusted_P) & Adjusted_P < alpha &
            log2FC <= -log2fc_threshold ~ paste0("Higher in ", a),
          TRUE ~ "Not significant"
        ),
        Comparison = comparison,
        Volcano_Y = -log10(
          pmax(Adjusted_P, .Machine$double.xmin)
        ),
        MA_X = log10(Mean_abundance + pseudocount)
      ) |>
      dplyr::arrange(Adjusted_P, dplyr::desc(abs(log2FC)))

    representative <- result |>
      dplyr::filter(!is.na(Adjusted_P)) |>
      dplyr::arrange(
        dplyr::desc(Significance != "Not significant"),
        Adjusted_P,
        dplyr::desc(abs(log2FC)),
        dplyr::desc(Mean_abundance)
      ) |>
      dplyr::slice_head(n = top_n_table)

    label_data <- result |>
      dplyr::filter(!is.na(Adjusted_P))

    if (label_only_significant) {
      label_data <- label_data |>
        dplyr::filter(Significance != "Not significant")
    }

    label_data <- label_data |>
      dplyr::arrange(
        Adjusted_P,
        dplyr::desc(abs(log2FC)),
        dplyr::desc(Mean_abundance)
      ) |>
      dplyr::slice_head(n = max_labels)

    plot_data <- result |>
      dplyr::filter(
        Tested, !is.na(Adjusted_P),
        is.finite(log2FC), is.finite(Volcano_Y),
        is.finite(MA_X)
      )

    point_mapping <- switch(paste(colour_by, point_size, sep = "_"),
      significance_constant =
        ggplot2::aes(colour = Significance),
      significance_mean_abundance =
        ggplot2::aes(
          colour = Significance,
          size = Mean_abundance
        ),
      log2FC_constant =
        ggplot2::aes(colour = log2FC),
      log2FC_mean_abundance =
        ggplot2::aes(
          colour = log2FC,
          size = Mean_abundance
        ),
      abundance_constant =
        ggplot2::aes(colour = Mean_abundance),
      abundance_mean_abundance =
        ggplot2::aes(
          colour = Mean_abundance,
          size = Mean_abundance
        )
    )

    point_layer <- if (
      identical(point_size, "constant")
    ) {
      ggplot2::geom_point(
        mapping = point_mapping,
        alpha = 0.75,
        size = 2.2
      )
    } else {
      ggplot2::geom_point(
        mapping = point_mapping,
        alpha = 0.75
      )
    }

    volcano <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = log2FC, y = Volcano_Y)
    ) +
      point_layer +
      ggplot2::geom_vline(
        xintercept = c(-log2fc_threshold, log2fc_threshold),
        linetype = "dashed",
        colour = aaa_visual_identity()$dark,
        linewidth = 0.55
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(alpha),
        linetype = "dashed",
        colour = aaa_visual_identity()$dark,
        linewidth = 0.55
      ) +
      ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(
          x = log2FC, y = Volcano_Y, label = Taxonomy
        ),
        inherit.aes = FALSE,
        size = 3,
        max.overlaps = max_labels,
        box.padding = 0.35,
        point.padding = 0.25,
        min.segment.length = 0,
        segment.colour = aaa_visual_identity()$secondary,
        show.legend = FALSE
      ) +
      ggplot2::labs(
        title = paste0(
          "Differential abundance: ", b, " vs ", a
        ),
        subtitle = paste0(
          "Positive log2FC indicates higher abundance in ", b
        ),
        x = paste0(
          "log2 fold change (", b, " / ", a, ")"
        ),
        y = expression(-log[10]("adjusted P value"))
      ) +
      aaa_theme() +
      aaa_legend_theme(
        n_items = if (colour_by == "significance") 3 else 1,
        position = "bottom",
        max_per_row = 3
      )

    volcano <- aaa_apply_differential_colour_scale(volcano, colour_by, point_size, a, b)

    volcano_x <- if (is.null(x_limit)) NULL else c(-x_limit, x_limit)
    volcano_y <- NULL

    if (dynamic_ylim && nrow(plot_data) > 0) {
      ymax <- max(plot_data$Volcano_Y, na.rm = TRUE)
      volcano_y <- c(
        0,
        max(ymax * 1.15, -log10(alpha) * 1.05)
      )
    }

    if (!is.null(volcano_x) || !is.null(volcano_y)) {
      volcano <- volcano +
        ggplot2::coord_cartesian(
          xlim = volcano_x,
          ylim = volcano_y
        )
    }

    ma_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = MA_X, y = log2FC)
    ) +
      point_layer +
      ggplot2::geom_hline(
        yintercept = c(
          -log2fc_threshold, 0, log2fc_threshold
        ),
        linetype = c("dashed", "solid", "dashed")
      ) +
      ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(
          x = MA_X, y = log2FC, label = Taxonomy
        ),
        inherit.aes = FALSE,
        size = 3,
        max.overlaps = max_labels,
        box.padding = 0.35,
        point.padding = 0.25,
        min.segment.length = 0,
        segment.colour = aaa_visual_identity()$secondary,
        show.legend = FALSE
      ) +
      ggplot2::labs(
        title = paste0("MA plot: ", b, " vs ", a),
        subtitle = paste0(
          "Positive log2FC indicates higher abundance in ", b
        ),
        x = expression(
          log[10]("mean relative abundance (%)")
        ),
        y = paste0(
          "log2 fold change (", b, " / ", a, ")"
        )
      ) +
      aaa_theme() +
      aaa_legend_theme(
        n_items = if (colour_by == "significance") 3 else 1,
        position = "bottom",
        max_per_row = 3
      )

    ma_plot <- aaa_apply_differential_colour_scale(ma_plot, colour_by, point_size, a, b)

    if (!is.null(x_limit)) {
      ma_plot <- ma_plot +
        ggplot2::coord_cartesian(
          ylim = c(-x_limit, x_limit)
        )
    }

    qq_data <- result[is.finite(result$P_value) & !is.na(result$P_value), , drop = FALSE]
    qq_data <- qq_data[order(qq_data$P_value), , drop = FALSE]
    if (nrow(qq_data) > 0L) {
      qq_data$Expected <- -log10(stats::ppoints(nrow(qq_data)))
      qq_data$Observed <- -log10(pmax(qq_data$P_value, .Machine$double.xmin))
      qq_plot <- ggplot2::ggplot(qq_data, ggplot2::aes(Expected, Observed)) +
        ggplot2::geom_point(alpha = 0.65) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
        ggplot2::labs(title = paste0("QQ plot — ", comparison), x = "Expected -log10(P)", y = "Observed -log10(P)") +
        aaa_theme()
    } else {
      qq_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No finite P-values available") +
        aaa_theme()
    }

    volcano_file <- file.path(
      prepared$project$analysis,
      paste0("Volcano_", aaa_safe_name(comparison), ".png")
    )

    ma_file <- file.path(
      prepared$project$analysis,
      paste0("MA_", aaa_safe_name(comparison), ".png")
    )
    qq_file <- file.path(
      prepared$project$analysis,
      paste0("QQ_", aaa_safe_name(comparison), ".png")
    )

    aaa_save_plot(
      volcano, volcano_file,
      width = 9,
      height = 6.8
    )

    aaa_save_plot(
      ma_plot, ma_file,
      width = 9,
      height = 6.8
    )
    aaa_save_plot(qq_plot, qq_file, width = 7, height = 6)

    all_results[[comparison]] <- result
    top_results[[comparison]] <- representative
    volcano_plots[[comparison]] <- volcano
    ma_plots[[comparison]] <- ma_plot
    qq_plots[[comparison]] <- qq_plot
    volcano_files[comparison] <- volcano_file
    ma_files[comparison] <- ma_file
    qq_files[comparison] <- qq_file
  }

  combined <- dplyr::bind_rows(all_results)
  representatives <- dplyr::bind_rows(top_results)
  differential_summary <- data.frame(
    Metric = c("Comparisons", "Taxa tested", "Significant results", "Higher in treatment B", "Higher in treatment A", "Adjusted P threshold", "Absolute log2FC threshold"),
    Value = c(
      length(comparisons), sum(combined$Tested, na.rm = TRUE),
      sum(combined$Significance != "Not significant", na.rm = TRUE),
      sum(grepl("^Higher in", combined$Significance) & combined$log2FC > 0, na.rm = TRUE),
      sum(grepl("^Higher in", combined$Significance) & combined$log2FC < 0, na.rm = TRUE),
      alpha, log2fc_threshold
    ),
    stringsAsFactors = FALSE
  )

  workbook <- c(
    list(
      Summary = differential_summary,
      Combined_results = combined,
      Representative_taxa = representatives
    ),
    stats::setNames(
      all_results,
      paste0("All_", seq_along(all_results))
    ),
    stats::setNames(
      top_results,
      paste0("Top_", seq_along(top_results))
    )
  )

  summary_file <- file.path(
    prepared$project$analysis,
    "aaa_differential_abundance_summary.xlsx"
  )

  openxlsx::write.xlsx(
    workbook,
    summary_file,
    overwrite = TRUE
  )

  aaa_result(
    tables = list(
      comparisons = all_results,
      representative = top_results,
      summary = differential_summary,
      combined = combined
    ),
    plots = list(
      volcano = volcano_plots,
      ma = ma_plots,
      qq = qq_plots
    ),
    files = c(
      summary = summary_file,
      stats::setNames(
        volcano_files,
        paste0("volcano_", names(volcano_files))
      ),
      stats::setNames(
        ma_files,
        paste0("ma_", names(ma_files))
      ),
      stats::setNames(
        qq_files,
        paste0("qq_", names(qq_files))
      )
    ),
    output_dir = prepared$project$analysis,
    metadata = list(
      method = method,
      paired = paired,
      alpha = alpha,
      log2fc_threshold = log2fc_threshold,
      max_labels = max_labels,
      label_only_significant = label_only_significant,
      colour_by = colour_by,
      point_size = point_size,
      x_limit = x_limit,
      dynamic_ylim = dynamic_ylim
    )
  )
}

#' Summarize selected functional pathways
#'
#' @param pathways Named list. Each entry contains `results` (data frame or
#'   aaa_functional_potential Excel file) and `include` (class names or a logical
#'   selector function).
#' @return A Triple_A_result object.
aaa_functional_abundance <- function(
  dataset, pathways,
  abundance_type = c("proportion", "percentage", "counts"),
  top_taxa_per_pathway = 10, project_dir,
  analysis_name = "Potential_metabolomic_pathways_abundance",
  graph_main = "Potential metabolomic pathways abundance"
) {
  abundance_type <- match.arg(abundance_type)
  aaa_check_packages(analyses = "functional_abundance")

  if (!is.list(pathways) || length(pathways) == 0 ||
    is.null(names(pathways)) || any(!nzchar(names(pathways)))) {
    stop("'pathways' must be a non-empty named list.")
  }

  prepared <- aaa_prepare_amplicon_data(
    dataset = dataset, abundance_type = abundance_type,
    project_dir = project_dir, analysis_name = analysis_name
  )

  read_results <- function(x) {
    if (inherits(x, "Triple_A_result")) {
      return(x$tables$taxa)
    }
    if (is.data.frame(x)) {
      return(x)
    }
    if (is.character(x) && length(x) == 1 && file.exists(x)) {
      extension <- tolower(tools::file_ext(x))
      if (!extension %in% c("xls", "xlsx")) {
        stop("Excel result paths must use .xls or .xlsx.")
      }
      sheets <- aaa_list_workbook_sheets(x)
      sheet <- if ("Taxon_results" %in% sheets) "Taxon_results" else sheets[[1L]]
      return(aaa_import_table(x, sheet = sheet))
    }
    stop("Each pathway 'results' must be a result object, data frame or Excel path.")
  }

  select_classes <- function(potential, include) {
    if (is.function(include)) {
      selected <- include(potential)
      if (!is.logical(selected) || length(selected) != length(potential)) {
        stop("A selector function must return one logical value per row.")
      }
      selected[is.na(selected)] <- FALSE
      selected
    } else if (is.character(include) && length(include) > 0) {
      !is.na(potential) & potential %in% include
    } else {
      stop("'include' must be class names or a selector function.")
    }
  }

  pathway_tables <- list()
  contribution_tables <- list()
  definitions <- list()

  for (name in names(pathways)) {
    cfg <- pathways[[name]]
    if (is.null(cfg$results) || is.null(cfg$include)) {
      stop("Pathway '", name, "' must contain 'results' and 'include'.")
    }
    functional <- read_results(cfg$results)
    required <- c("Taxonomy", "Genus", "Tax_level", "Potential")
    missing <- setdiff(required, names(functional))
    if (length(missing) > 0) {
      stop(
        "Pathway '", name, "' is missing: ",
        paste(missing, collapse = ", ")
      )
    }

    selected <- functional[
      select_classes(functional$Potential, cfg$include),
      required,
      drop = FALSE
    ] |>
      dplyr::distinct(Taxonomy, Genus, Tax_level, .keep_all = TRUE)

    contributions <- prepared$wide |>
      dplyr::inner_join(
        selected,
        by = c("Taxonomy", "Genus", "Tax_level")
      ) |>
      dplyr::mutate(Pathway = name)

    values <- if (nrow(contributions) == 0) {
      row <- tibble::tibble(Pathway = name)
      for (column in prepared$sample_columns) row[[column]] <- 0
      row
    } else {
      contributions |>
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(prepared$sample_columns),
            \(x) sum(x, na.rm = TRUE)
          )
        ) |>
        dplyr::mutate(Pathway = name, .before = 1)
    }

    contribution_long <- contributions |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(prepared$sample_columns),
        names_to = "Sample_column", values_to = "Abundance"
      ) |>
      dplyr::left_join(prepared$sample_map, by = "Sample_column") |>
      dplyr::group_by(
        Pathway, Taxonomy, Genus, Tax_level, Potential, Treatment
      ) |>
      dplyr::summarise(
        Mean_abundance = mean(Abundance, na.rm = TRUE),
        SD = stats::sd(Abundance, na.rm = TRUE),
        .groups = "drop"
      )

    pathway_tables[[name]] <- values
    contribution_tables[[name]] <- contribution_long
    definitions[[name]] <- tibble::tibble(
      Pathway = name,
      Included_classes = if (is.function(cfg$include)) {
        "Custom selector function"
      } else {
        paste(cfg$include, collapse = " | ")
      },
      Selected_taxa = nrow(selected)
    )
  }

  pathway_by_sample <- dplyr::bind_rows(pathway_tables)
  matrix <- pathway_by_sample |>
    tibble::column_to_rownames("Pathway") |>
    as.data.frame()

  if (isTRUE(prepared$has_replicates)) {
    summary <- aaa_replicate_summary(
      matrix, prepared$sample_map, prepared$samples_name
    )
    heatmap_data <- summary$mean
    heatmap_labels <- summary$labels
    sd_data <- summary$sd
  } else {
    heatmap_data <- matrix
    names(heatmap_data) <- prepared$samples_name
    sd_data <- heatmap_data
    sd_data[, ] <- NA_real_
    heatmap_labels <- matrix(
      sprintf("%.1f%%", as.matrix(heatmap_data)),
      nrow = nrow(heatmap_data), dimnames = dimnames(heatmap_data)
    )
  }

  treatment_order <- colnames(heatmap_data)

  summary_long <- tibble::rownames_to_column(heatmap_data, "Pathway") |>
    tidyr::pivot_longer(-Pathway,
      names_to = "Treatment",
      values_to = "Mean_abundance"
    ) |>
    dplyr::left_join(
      tibble::rownames_to_column(sd_data, "Pathway") |>
        tidyr::pivot_longer(-Pathway,
          names_to = "Treatment",
          values_to = "SD"
        ),
      by = c("Pathway", "Treatment")
    ) |>
    dplyr::mutate(
      Treatment = factor(
        Treatment,
        levels = treatment_order,
        ordered = TRUE
      )
    )

  heatmap_file <- file.path(
    prepared$project$analysis, "Potential_metabolomic_pathways_abundance_heatmap.png"
  )
  aaa_plot_heatmap(
    heatmap_data,
    filename = heatmap_file,
    title = graph_main,
    labels = heatmap_labels,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    height = max(3.5, 0.5 * nrow(heatmap_data) + 1.5)
  )

  bar_plot <- ggplot2::ggplot(
    summary_long,
    ggplot2::aes(x = Treatment, y = Mean_abundance, fill = Pathway)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8), width = 0.75
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(Mean_abundance - SD, 0),
        ymax = Mean_abundance + SD
      ),
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.2, na.rm = TRUE
    ) +
    ggplot2::labs(
      title = graph_main,
      subtitle = "Potential pathways may overlap taxonomically",
      x = NULL, y = "Mean relative abundance (%)", fill = "Pathway"
    ) +
    ggplot2::scale_fill_manual(
      values = aaa_treatment_palette(length(unique(summary_long$Pathway))),
      labels = aaa_wrap_legend_labels(
        unique(summary_long$Pathway),
        width = 30
      ),
      guide = ggplot2::guide_legend(
        ncol = aaa_legend_columns(
          length(unique(summary_long$Pathway)),
          max_per_row = 4
        ),
        byrow = TRUE,
        title.position = "top"
      )
    ) +
    aaa_theme() +
    aaa_legend_theme(
      n_items = length(unique(summary_long$Pathway)),
      position = "bottom",
      max_per_row = 4,
      text_size = 8.5
    )

  bar_file <- file.path(
    prepared$project$analysis, "Potential_metabolomic_pathways_abundance_barplot.png"
  )
  aaa_save_plot(
    bar_plot, bar_file,
    width = max(
      9,
      1.1 * length(unique(summary_long$Treatment)) + 5
    ),
    height = aaa_plot_height_with_legend(
      base_height = 5.5,
      n_items = length(unique(summary_long$Pathway)),
      max_per_row = 4,
      row_height = 0.5
    )
  )

  contributions <- dplyr::bind_rows(contribution_tables)
  top_contributors <- contributions |>
    dplyr::group_by(Pathway, Taxonomy) |>
    dplyr::summarise(
      Mean_abundance = mean(Mean_abundance, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(Pathway) |>
    dplyr::slice_max(
      Mean_abundance,
      n = top_taxa_per_pathway, with_ties = FALSE
    ) |>
    dplyr::ungroup()

  # A pathway definition can be valid even when no taxa from the input table
  # match it. In that case, keep the analysis successful and generate an
  # informative placeholder instead of asking ggplot2 to facet an empty table.
  has_contributors <- nrow(top_contributors) > 0L &&
    "Pathway" %in% names(top_contributors) &&
    any(!is.na(top_contributors$Pathway) & nzchar(as.character(top_contributors$Pathway)))

  if (has_contributors) {
    contributor_plot <- ggplot2::ggplot(
      top_contributors,
      ggplot2::aes(
        x = stats::reorder(Taxonomy, Mean_abundance),
        y = Mean_abundance,
        fill = Pathway
      )
    ) +
      ggplot2::geom_col(
        width = 0.72,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = aaa_treatment_palette(length(unique(top_contributors$Pathway)))
      ) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~Pathway, scales = "free", ncol = 1) +
      ggplot2::labs(
        title = "Main taxonomic contributors to potential metabolomic pathways",
        x = NULL, y = "Mean relative abundance (%)"
      ) +
      aaa_theme() +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 10, face = "bold"), axis.text.y = ggplot2::element_text(size = 9), plot.title = ggplot2::element_text(size = 14))
  } else {
    contributor_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0, y = 0,
        label = paste(
          "No matching taxonomic contributors were detected",
          "for the selected potential metabolomic pathways."
        ),
        size = 4.2
      ) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::labs(
        title = "Main taxonomic contributors to potential metabolomic pathways",
        subtitle = "The analysis completed successfully; no contributor bars are available.",
        x = NULL, y = NULL
      ) +
      aaa_theme() +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }

  contributor_file <- file.path(
    prepared$project$analysis, "Potential_metabolic_pathways_top_contributors.png"
  )
  # After coord_flip(), Taxonomy values (routinely full lineage strings well
  # past 100 characters) become the left-margin axis labels; a fixed width
  # clips them instead of the plot growing to fit, same failure mode already
  # fixed for pheatmap row labels via aaa_measure_text_width_inches().
  aaa_save_plot(
    contributor_plot, contributor_file,
    width = if (has_contributors) {
      max(13, aaa_flipped_axis_plot_width(unique(top_contributors$Taxonomy)))
    } else {
      13
    },
    height = max(7, 2.2 * length(unique(top_contributors$Pathway)))
  )

  dot_plot <- ggplot2::ggplot(
    summary_long,
    ggplot2::aes(
      x = Treatment,
      y = Pathway,
      size = Mean_abundance,
      colour = Mean_abundance
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_colour_gradient(
      low = aaa_visual_identity()$light,
      high = aaa_visual_identity()$primary
    ) +
    ggplot2::scale_size_continuous(range = c(2, 12)) +
    ggplot2::labs(
      title = paste0(graph_main, ": dot plot"),
      x = NULL,
      y = NULL,
      colour = "Mean abundance (%)",
      size = "Mean abundance (%)"
    ) +
    aaa_theme() +
    aaa_legend_theme(
      n_items = 2,
      position = "bottom",
      max_per_row = 2
    )

  dot_file <- file.path(
    prepared$project$analysis,
    "Potential_metabolomic_pathways_abundance_dotplot.png"
  )
  aaa_save_plot(
    dot_plot, dot_file,
    width = max(
      9,
      1.05 * length(unique(summary_long$Treatment)) + 5
    ),
    height = max(5.5, 0.60 * length(pathways) + 2.5)
  )

  definitions <- dplyr::bind_rows(definitions)
  summary_file <- file.path(
    prepared$project$analysis, "Potential_metabolomic_pathways_abundance_summary.xlsx"
  )
  openxlsx::write.xlsx(
    list(
      Pathway_definitions = definitions,
      Pathway_by_sample = pathway_by_sample,
      Pathway_summary = summary_long,
      Taxon_contributions = contributions,
      Top_contributors = top_contributors
    ),
    summary_file,
    overwrite = TRUE
  )

  aaa_result(
    tables = list(
      definitions = definitions,
      pathway_by_sample = pathway_by_sample,
      summary = summary_long,
      contributions = contributions,
      top_contributors = top_contributors
    ),
    plots = list(
      heatmap = heatmap_file,
      barplot = bar_plot,
      dotplot = dot_plot,
      top_contributors = contributor_plot
    ),
    files = c(
      summary = summary_file,
      heatmap = heatmap_file,
      barplot = bar_file,
      dotplot = dot_file,
      top_contributors = contributor_file
    ),
    output_dir = prepared$project$analysis,
    metadata = list(pathways = names(pathways))
  )
}
