# =============================================================================
# Comparison of inferred functional profiles between experimental groups.
#
# Both analyses here consume results the engine already produces rather than
# recomputing anything: differential functions tests the pathway-by-sample
# matrix built by aaa_functional_abundance(), and enrichment crosses the
# per-taxon functional calls from aaa_functional_potential() with the
# differential-abundance results.
#
# Neither reports on measured activity. They describe inferred genomic
# potential, so a "differential function" means the potential is differentially
# represented, not that the pathway is more active.
# =============================================================================

# Differential functions --------------------------------------------------------
# Wilcoxon rank-sum per pathway for every pair of treatments, Benjamini-Hochberg
# corrected within each comparison. The non-parametric test matches the one used
# for pairwise differential abundance: pathway abundances inherit the skew of the
# taxon abundances they are summed from, so normality cannot be assumed.
aaa_differential_functions <- function(
  pathway_by_sample,
  sample_design,
  alpha = 0.05,
  project_dir,
  analysis_name = "Differential_functions"
) {
  if (!is.data.frame(pathway_by_sample) || !nrow(pathway_by_sample)) {
    stop("'pathway_by_sample' must be the pathway-by-sample table of the functional abundance analysis.")
  }
  if (!"Pathway" %in% names(pathway_by_sample)) {
    stop("'pathway_by_sample' must contain a 'Pathway' column.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number strictly between 0 and 1.")
  }

  design <- as.data.frame(sample_design, stringsAsFactors = FALSE)
  required <- c("Sample_column", "Treatment")
  if (!all(required %in% names(design))) {
    stop("'sample_design' must contain the columns: ", paste(required, collapse = ", "), ".")
  }

  sample_columns <- intersect(as.character(design$Sample_column), names(pathway_by_sample))
  if (length(sample_columns) < 2L) {
    stop("The functional table and the sample design share fewer than two samples.")
  }
  design <- design[match(sample_columns, design$Sample_column), , drop = FALSE]
  treatments <- unique(as.character(design$Treatment))
  if (length(treatments) < 2L) {
    stop("Differential functions requires at least two treatment groups.")
  }

  values <- as.matrix(pathway_by_sample[, sample_columns, drop = FALSE])
  rownames(values) <- as.character(pathway_by_sample$Pathway)
  storage.mode(values) <- "double"

  project <- aaa_create_project_structure(project_dir, analysis_name)

  pairs <- utils::combn(treatments, 2L, simplify = FALSE)
  per_comparison <- lapply(pairs, function(pair) {
    a <- pair[[1L]]
    b <- pair[[2L]]
    idx_a <- which(design$Treatment == a)
    idx_b <- which(design$Treatment == b)

    rows <- lapply(seq_len(nrow(values)), function(i) {
      x <- values[i, idx_a]
      y <- values[i, idx_b]
      mean_a <- mean(x, na.rm = TRUE)
      mean_b <- mean(y, na.rm = TRUE)
      # A pathway that is constant across every sample carries no information and
      # would make the test emit a warning rather than a usable P value.
      testable <- length(idx_a) >= 2L && length(idx_b) >= 2L &&
        stats::var(c(x, y), na.rm = TRUE) > 0
      p <- if (testable) {
        suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
      } else {
        NA_real_
      }
      data.frame(
        Pathway = rownames(values)[i],
        Treatment_A = a, Treatment_B = b,
        Mean_A = mean_a, Mean_B = mean_b,
        Difference = mean_b - mean_a,
        log2FC = log2((mean_b + 1e-6) / (mean_a + 1e-6)),
        P_value = p, Tested = testable,
        stringsAsFactors = FALSE
      )
    })

    out <- do.call(rbind, rows)
    out$Adjusted_P <- stats::p.adjust(out$P_value, method = "BH")
    out$Significance <- ifelse(
      !is.na(out$Adjusted_P) & out$Adjusted_P < alpha,
      ifelse(out$log2FC > 0, paste0("Higher in ", b), paste0("Higher in ", a)),
      "Not significant"
    )
    out$Comparison <- paste(a, "vs", b)
    out
  })

  results <- do.call(rbind, per_comparison)
  results <- results[order(results$Adjusted_P, na.last = TRUE), , drop = FALSE]
  rownames(results) <- NULL

  significant <- results[results$Significance != "Not significant", , drop = FALSE]
  plot_data <- if (nrow(significant)) significant else utils::head(results, 20L)

  figure <- if (!nrow(plot_data) || all(is.na(plot_data$Adjusted_P))) {
    aaa_community_placeholder_plot(
      "Differential functions",
      "No pathway could be tested with the current design."
    )
  } else {
    plot_data$Label <- paste0(plot_data$Pathway, " (", plot_data$Comparison, ")")
    plot_data$Label <- factor(plot_data$Label, levels = rev(unique(plot_data$Label)))
    ident <- aaa_visual_identity()
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = log2FC, y = Label, fill = Significance != "Not significant")
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_vline(xintercept = 0, colour = ident$dark, linewidth = 0.4) +
      ggplot2::scale_fill_manual(
        values = c(`TRUE` = ident$positive, `FALSE` = "grey75"),
        labels = c(`TRUE` = "Significant", `FALSE` = "Not significant"),
        name = NULL
      ) +
      ggplot2::labs(
        title = "Differential functions",
        subtitle = sprintf(
          "Wilcoxon rank-sum per pathway, Benjamini-Hochberg adjusted (alpha = %s)", alpha
        ),
        x = "log2 fold change", y = NULL
      ) +
      aaa_theme()
  }

  files <- c(
    differential_functions = file.path(project$analysis, "Differential_functions.png"),
    summary = file.path(project$analysis, "Differential_functions_summary.xlsx")
  )
  aaa_save_plot(
    figure, files[["differential_functions"]],
    width = aaa_flipped_axis_plot_width(as.character(plot_data$Pathway), base_width = 7),
    height = max(4, min(14, 0.32 * max(1L, nrow(plot_data)) + 2))
  )
  openxlsx::write.xlsx(list(Differential_functions = results), files[["summary"]], overwrite = TRUE)

  aaa_result(

    tables = list(comparisons = results, significant = significant),
    plots = list(differential_functions = figure),
    files = files,
    output_dir = project$analysis,
    metadata = list(
      alpha = alpha,
      n_pathways = nrow(values),
      n_comparisons = length(pairs),
      n_significant = nrow(significant),
      test = "Wilcoxon rank-sum",
      correction = "BH"
    )
  )
}


# Functional enrichment ---------------------------------------------------------
# Fisher's exact test asking whether the taxa carrying a given function are
# over-represented among the differentially abundant taxa. The background is
# every taxon that received a functional call, not every taxon in the dataset:
# a taxon with no reference genome could never have been called positive, so
# including it would inflate the apparent enrichment.
aaa_functional_enrichment <- function(
  functional_potential,
  differential_abundance,
  alpha = 0.05,
  project_dir,
  analysis_name = "Functional_enrichment"
) {
  if (!is.list(functional_potential) || !length(functional_potential)) {
    stop("'functional_potential' must be the named list of functional-potential results.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number strictly between 0 and 1.")
  }

  comparisons <- differential_abundance$tables$comparisons
  if (!is.data.frame(comparisons) || !nrow(comparisons)) {
    stop("Functional enrichment requires the comparison table of a differential-abundance result.")
  }

  significant_taxa <- unique(as.character(
    comparisons$Taxonomy[comparisons$Significance != "Not significant"]
  ))

  rows <- lapply(names(functional_potential), function(id) {
    taxa <- functional_potential[[id]]$tables$taxa
    if (!is.data.frame(taxa) || !nrow(taxa) || !"Potential" %in% names(taxa)) return(NULL)

    # A call is positive when it is neither a negative nor an insufficient-evidence
    # verdict; those two are recorded as text rather than as NA.
    call <- tolower(trimws(as.character(taxa$Potential)))
    negative <- is.na(call) | !nzchar(call) |
      grepl("^no |^not |insufficient|negative|unassigned|unknown", call)
    positive_taxa <- unique(as.character(taxa$Taxonomy[!negative]))
    background <- unique(as.character(taxa$Taxonomy))

    sig <- intersect(significant_taxa, background)
    a <- length(intersect(positive_taxa, sig))          # positive and significant
    b <- length(setdiff(positive_taxa, sig))            # positive, not significant
    c <- length(setdiff(sig, positive_taxa))            # significant, not positive
    d <- length(background) - a - b - c                 # neither

    if (any(c(a, b, c, d) < 0)) return(NULL)
    test <- suppressWarnings(tryCatch(
      stats::fisher.test(matrix(c(a, b, c, d), nrow = 2L)),
      error = function(e) NULL
    ))

    data.frame(
      Function = id,
      Taxa_with_function = length(positive_taxa),
      Taxa_in_background = length(background),
      Significant_taxa = length(sig),
      Overlap = a,
      Expected = if (length(background) > 0L) {
        length(positive_taxa) * length(sig) / length(background)
      } else NA_real_,
      Odds_ratio = if (is.null(test)) NA_real_ else unname(test$estimate),
      P_value = if (is.null(test)) NA_real_ else test$p.value,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(results) || !nrow(results)) {
    stop("No functional category could be tested for enrichment.")
  }
  results$Adjusted_P <- stats::p.adjust(results$P_value, method = "BH")
  results$Enriched <- !is.na(results$Adjusted_P) & results$Adjusted_P < alpha &
    !is.na(results$Odds_ratio) & results$Odds_ratio > 1
  results <- results[order(results$Adjusted_P, na.last = TRUE), , drop = FALSE]
  rownames(results) <- NULL

  project <- aaa_create_project_structure(project_dir, analysis_name)

  plot_data <- results[!is.na(results$P_value), , drop = FALSE]
  figure <- if (!nrow(plot_data)) {
    aaa_community_placeholder_plot(
      "Functional enrichment",
      "No functional category could be tested."
    )
  } else {
    ident <- aaa_visual_identity()
    plot_data$Function <- factor(
      plot_data$Function, levels = rev(unique(plot_data$Function))
    )
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = -log10(pmax(Adjusted_P, .Machine$double.xmin)),
        y = Function, size = Overlap, colour = Enriched
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(
        xintercept = -log10(alpha), linetype = "dashed", colour = ident$dark
      ) +
      ggplot2::scale_colour_manual(
        values = c(`TRUE` = ident$positive, `FALSE` = "grey65"),
        labels = c(`TRUE` = "Enriched", `FALSE` = "Not enriched"), name = NULL
      ) +
      ggplot2::scale_size_continuous(name = "Overlapping taxa") +
      ggplot2::labs(
        title = "Functional enrichment",
        subtitle = sprintf(
          "Fisher's exact test per function, Benjamini-Hochberg adjusted (alpha = %s)", alpha
        ),
        x = "-log10 adjusted P", y = NULL
      ) +
      aaa_theme()
  }

  files <- c(
    functional_enrichment = file.path(project$analysis, "Functional_enrichment.png"),
    summary = file.path(project$analysis, "Functional_enrichment_summary.xlsx")
  )
  aaa_save_plot(
    figure, files[["functional_enrichment"]],
    width = aaa_flipped_axis_plot_width(as.character(results$Function), base_width = 7),
    height = max(4, min(14, 0.35 * nrow(results) + 2))
  )
  openxlsx::write.xlsx(list(Functional_enrichment = results), files[["summary"]], overwrite = TRUE)

  aaa_result(

    tables = list(enrichment = results),
    plots = list(functional_enrichment = figure),
    files = files,
    output_dir = project$analysis,
    metadata = list(
      alpha = alpha,
      n_functions = nrow(results),
      n_enriched = sum(results$Enriched, na.rm = TRUE),
      significant_taxa = length(significant_taxa),
      test = "Fisher's exact test",
      correction = "BH"
    )
  )
}
