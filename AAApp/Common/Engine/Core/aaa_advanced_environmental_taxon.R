# =============================================================================
# Triple_A advanced environmental and taxon-association analyses
# =============================================================================

aaa_analysis_metadata_frame <- function(dataset, variables = NULL) {
  aaa_validate_dataset(dataset)
  if (is.null(dataset$metadata)) stop("No metadata table is attached to this dataset.")
  roles <- dataset$metadata_roles
  id_col <- roles$Column[roles$Role == "identifier"]
  if (length(id_col) != 1L) stop("Exactly one metadata identifier is required.")
  if (is.null(variables)) variables <- setdiff(names(dataset$metadata), id_col)
  missing <- setdiff(variables, names(dataset$metadata))
  if (length(missing)) stop("Metadata variables not found: ", paste(missing, collapse = ", "))
  out <- dataset$metadata[, c(id_col, variables), drop = FALSE]
  names(out)[1] <- "Sample_column"
  out
}

aaa_align_analysis_inputs <- function(dataset, abundance_type, transformation,
                                      variables, project_dir, analysis_name) {
  inp <- aaa_multivariate_input(dataset, abundance_type, transformation, project_dir, analysis_name)
  meta <- aaa_analysis_metadata_frame(dataset, variables)
  meta <- meta[match(rownames(inp$matrix), meta$Sample_column), , drop = FALSE]
  if (anyNA(meta$Sample_column)) stop("Some abundance-table samples are missing from metadata.")
  rownames(meta) <- meta$Sample_column
  list(X = inp$matrix, raw = inp$relative, metadata = meta, input = inp)
}

aaa_clean_model_frame <- function(metadata, variables) {
  d <- as.data.frame(metadata[, variables, drop = FALSE], stringsAsFactors = FALSE)
  for (nm in names(d)) {
    x <- d[[nm]]
    num <- suppressWarnings(as.numeric(x))
    if (sum(is.finite(num)) >= max(3L, ceiling(.8 * length(x)))) d[[nm]] <- num else d[[nm]] <- factor(x)
  }
  keep <- vapply(d, function(x) {
    y <- x[!is.na(x)]
    length(unique(y)) > 1L
  }, logical(1))
  d <- d[, keep, drop = FALSE]
  if (!ncol(d)) stop("No informative metadata variables remain after validation.")
  d
}

aaa_save_standard_result <- function(project_dir, analysis_name, tables, plots = list(), metadata = list()) {
  project <- aaa_create_project_structure(project_dir, analysis_name)
  files <- c(summary = file.path(project$analysis, paste0(analysis_name, "_summary.xlsx")))
  openxlsx::write.xlsx(tables, files[["summary"]], overwrite = TRUE)
  if (length(plots)) {
    for (nm in names(plots)) {
      f <- file.path(project$analysis, paste0(analysis_name, "_", nm, ".png"))
      aaa_save_plot(plots[[nm]], f, 8, 6)
      files[nm] <- f
    }
  }
  aaa_result(tables = tables, plots = plots, files = files, output_dir = project$analysis, metadata = metadata)
}

# =============================================================================
# Optional graphical summaries for the taxon-association engines. Both draw
# only from the result tables the engines already produce (no new statistics):
# ANCOM-BC2 gets a volcano plot; MaAsLin2 gets a coefficient forest plot. The
# figures appear automatically in the result tree because they are stored in
# the result's files list (these analyses are not part of the figure-selection
# catalogue, so they always accompany their table).
# =============================================================================

aaa_ancombc2_volcano_plot <- function(res, alpha = 0.05,
                                      title = "ANCOM-BC2 differential abundance") {
  lfc_cols <- setdiff(grep("^lfc_", names(res), value = TRUE), "lfc_(Intercept)")
  if (!length(lfc_cols)) {
    return(aaa_community_placeholder_plot(
      title, "No non-intercept ANCOM-BC2 contrasts are available for a volcano plot."
    ))
  }
  taxon <- if ("taxon" %in% names(res)) {
    as.character(res$taxon)
  } else if ("Taxon" %in% names(res)) {
    as.character(res$Taxon)
  } else {
    as.character(seq_len(nrow(res)))
  }

  long <- do.call(rbind, lapply(lfc_cols, function(lc) {
    term <- sub("^lfc_", "", lc)
    qc <- paste0("q_", term)
    dc <- paste0("diff_", term)
    q <- if (qc %in% names(res)) suppressWarnings(as.numeric(res[[qc]])) else NA_real_
    diff <- if (dc %in% names(res)) as.logical(res[[dc]]) else (is.finite(q) & q <= alpha)
    data.frame(
      Taxon = taxon, Term = term,
      lfc = suppressWarnings(as.numeric(res[[lc]])), q = q,
      Significant = ifelse(is.na(diff), FALSE, diff), stringsAsFactors = FALSE
    )
  }))
  long <- long[is.finite(long$lfc) & is.finite(long$q), , drop = FALSE]
  if (!nrow(long)) {
    return(aaa_community_placeholder_plot(
      title, "ANCOM-BC2 produced no finite effect sizes to plot."
    ))
  }
  long$Volcano_Y <- -log10(pmax(long$q, .Machine$double.xmin))

  ident <- aaa_visual_identity()
  labels <- long[long$Significant, , drop = FALSE]
  labels <- labels[order(labels$q, -abs(labels$lfc)), , drop = FALSE]
  labels <- utils::head(labels, 12)

  plot <- ggplot2::ggplot(long, ggplot2::aes(x = lfc, y = Volcano_Y, colour = Significant)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey70") +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = ident$dark) +
    ggplot2::geom_point(alpha = 0.8, size = 2.2) +
    ggrepel::geom_text_repel(
      data = labels, ggplot2::aes(label = Taxon), size = 3,
      max.overlaps = 12, box.padding = 0.35, min.segment.length = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = ident$positive, `FALSE` = ident$neutral),
      labels = c(`TRUE` = paste0("q ≤ ", alpha), `FALSE` = "Not significant"),
      name = NULL, drop = FALSE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Bias-corrected log fold change vs BH-adjusted significance",
      x = "Log fold change", y = expression(-log[10]("adjusted P value"))
    ) +
    aaa_theme()

  if (length(unique(long$Term)) > 1L) {
    plot <- plot + ggplot2::facet_wrap(~Term, scales = "free")
  }
  plot
}

aaa_maaslin_forest_plot <- function(res, alpha = 0.05, top_n = 25L,
                                    title = "MaAsLin2 associations") {
  if (!is.data.frame(res) || !all(c("coef", "stderr", "qval") %in% names(res)) ||
    !nrow(res)) {
    return(aaa_community_placeholder_plot(
      title, "MaAsLin2 produced no coefficients to plot."
    ))
  }
  df <- res[is.finite(res$coef) & is.finite(res$stderr), , drop = FALSE]
  if (!nrow(df)) {
    return(aaa_community_placeholder_plot(
      title, "MaAsLin2 produced no finite coefficients to plot."
    ))
  }
  df$Feature <- if ("feature" %in% names(df)) {
    as.character(df$feature)
  } else {
    as.character(seq_len(nrow(df)))
  }
  df$Term <- if (all(c("metadata", "value") %in% names(df))) {
    paste0(as.character(df$metadata), ": ", as.character(df$value))
  } else if ("metadata" %in% names(df)) {
    as.character(df$metadata)
  } else {
    "association"
  }
  df$Row <- make.unique(paste0(df$Feature, "  (", df$Term, ")"))
  df$Significant <- is.finite(df$qval) & df$qval <= alpha

  df <- df[order(df$qval, -abs(df$coef)), , drop = FALSE]
  df <- utils::head(df, top_n)

  ident <- aaa_visual_identity()
  ggplot2::ggplot(df, ggplot2::aes(
    x = coef, y = stats::reorder(Row, coef),
    colour = Significant
  )) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = coef - 1.96 * stderr, xmax = coef + 1.96 * stderr)
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = ident$positive, `FALSE` = ident$neutral),
      labels = c(`TRUE` = paste0("q ≤ ", alpha), `FALSE` = "Not significant"),
      name = NULL, drop = FALSE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Model coefficients with approximate 95% confidence intervals",
      x = "Coefficient (effect size)", y = NULL
    ) +
    aaa_theme()
}

aaa_envfit_analysis <- function(dataset, abundance_type, environmental_variables,
                                transformation = "hellinger", permutations = 999,
                                project_dir, analysis_name = "envfit") {
  aaa_check_packages(analyses = "envfit")
  a <- aaa_align_analysis_inputs(dataset, abundance_type, transformation, environmental_variables, project_dir, analysis_name)
  env <- aaa_clean_model_frame(a$metadata, environmental_variables)
  complete <- complete.cases(env)
  X <- a$X[complete, , drop = FALSE]
  env <- env[complete, , drop = FALSE]
  if (nrow(X) < 4L) stop("envfit requires at least four complete samples.")
  ord <- vegan::rda(X)
  fit <- vegan::envfit(ord, env, permutations = as.integer(permutations), na.rm = TRUE)
  scores <- as.data.frame(vegan::scores(ord, display = "sites", choices = 1:2))
  names(scores)[1:2] <- c("Axis1", "Axis2")
  scores$Sample_column <- rownames(scores)
  scores$Treatment <- a$input$metadata[rownames(scores), "Treatment"]
  vectors <- if (!is.null(fit$vectors)) data.frame(Variable = rownames(fit$vectors$arrows), fit$vectors$arrows, R2 = fit$vectors$r, P_value = fit$vectors$pvals, row.names = NULL, check.names = FALSE) else data.frame()
  if (nrow(vectors) && ncol(vectors) >= 3L) names(vectors)[2:3] <- c("Axis1", "Axis2")
  # vegan::envfit() names factors$centroids rows "<Variable><Level>" (no
  # separator) and factors$pvals by variable only (one P-value per factor,
  # not per level); stripping up to "=" was a no-op against that format, so
  # this always produced NA. factors$var.id maps each centroid row to its
  # parent variable and is the correct key into pvals.
  factors <- if (!is.null(fit$factors)) data.frame(Level = rownames(fit$factors$centroids), fit$factors$centroids, P_value = unname(fit$factors$pvals[fit$factors$var.id]), row.names = NULL, check.names = FALSE) else data.frame()
  plot <- aaa_ordination_plot(scores, "Axis1", "Axis2", "PC1", "PC2", "Environmental fit (envfit)", show_labels = FALSE)
  if (nrow(vectors)) plot <- plot + ggplot2::geom_segment(data = vectors, ggplot2::aes(x = 0, y = 0, xend = Axis1, yend = Axis2), inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.18, "cm"))) + ggplot2::geom_text(data = vectors, ggplot2::aes(x = Axis1, y = Axis2, label = Variable), inherit.aes = FALSE)
  aaa_save_standard_result(project_dir, analysis_name, list(Vectors = vectors, Factor_centroids = factors, Sites = scores), list(envfit = plot), list(permutations = permutations))
}

aaa_constrained_analysis <- function(dataset, abundance_type, environmental_variables,
                                     conditioning_variables = character(), transformation = "hellinger",
                                     distance = NULL, permutations = 999, project_dir,
                                     analysis_name = "Constrained_ordination") {
  aaa_check_packages(analyses = if (is.null(distance)) "partial_rda" else "dbrda")
  vars <- unique(c(environmental_variables, conditioning_variables))
  a <- aaa_align_analysis_inputs(dataset, abundance_type, transformation, vars, project_dir, analysis_name)
  d <- aaa_clean_model_frame(a$metadata, vars)
  complete <- complete.cases(d)
  X <- a$X[complete, , drop = FALSE]
  d <- d[complete, , drop = FALSE]
  env <- intersect(environmental_variables, names(d))
  cond <- intersect(conditioning_variables, names(d))
  if (!length(env)) stop("Select at least one informative explanatory variable.")
  rhs <- paste(sprintf("`%s`", env), collapse = " + ")
  if (length(cond)) rhs <- paste0(rhs, " + Condition(", paste(sprintf("`%s`", cond), collapse = " + "), ")")
  form <- stats::as.formula(paste("X ~", rhs))
  model <- if (is.null(distance)) vegan::rda(form, data = d) else vegan::capscale(form, data = d, distance = distance, add = TRUE)
  anova_global <- as.data.frame(vegan::anova.cca(model, permutations = permutations))
  anova_global$Term <- rownames(anova_global)
  rownames(anova_global) <- NULL
  anova_terms <- as.data.frame(vegan::anova.cca(model, by = "term", permutations = permutations))
  anova_terms$Term <- rownames(anova_terms)
  rownames(anova_terms) <- NULL
  scr <- as.data.frame(vegan::scores(model, display = "sites", choices = 1:2))
  if (ncol(scr) < 2) scr$Axis2 <- 0
  names(scr)[1:2] <- c("Axis1", "Axis2")
  scr$Sample_column <- rownames(scr)
  scr$Treatment <- a$input$metadata[rownames(scr), "Treatment"]
  eig <- vegan::eigenvals(model)
  variance <- data.frame(Axis = names(eig), Eigenvalue = as.numeric(eig), Percent = 100 * as.numeric(eig) / sum(abs(eig)))
  ttl <- if (is.null(distance)) if (length(cond)) "Partial redundancy analysis" else "Redundancy analysis" else if (length(cond)) "Partial distance-based RDA" else "Distance-based RDA"
  plot <- aaa_ordination_plot(scr, "Axis1", "Axis2", names(scr)[1], names(scr)[2], ttl, show_labels = FALSE)
  aaa_save_standard_result(project_dir, analysis_name, list(Global_test = anova_global, Terms = anova_terms, Sites = scr, Variance = variance, Model_variables = data.frame(Explanatory = env, Conditioning = paste(cond, collapse = ", "))), list(ordination = plot), list(distance = distance, permutations = permutations))
}

aaa_variance_partitioning <- function(dataset, abundance_type, environmental_variables,
                                      experimental_factors, transformation = "hellinger",
                                      project_dir, analysis_name = "Variance_partitioning") {
  aaa_check_packages(analyses = "variance_partitioning")
  if (!length(environmental_variables) || !length(experimental_factors)) stop("Variance partitioning requires at least one environmental variable and one experimental factor.")
  vars <- unique(c(environmental_variables, experimental_factors))
  a <- aaa_align_analysis_inputs(dataset, abundance_type, transformation, vars, project_dir, analysis_name)
  d <- aaa_clean_model_frame(a$metadata, vars)
  complete <- complete.cases(d)
  X <- a$X[complete, , drop = FALSE]
  d <- d[complete, , drop = FALSE]
  e <- intersect(environmental_variables, names(d))
  f <- intersect(experimental_factors, names(d))
  if (!length(e) || !length(f)) stop("Both predictor sets must remain informative.")
  env_matrix <- stats::model.matrix(stats::as.formula(paste("~", paste(sprintf("`%s`", e), collapse = "+"))), d)[, -1, drop = FALSE]
  factor_matrix <- stats::model.matrix(stats::as.formula(paste("~", paste(sprintf("`%s`", f), collapse = "+"))), d)[, -1, drop = FALSE]
  vp <- vegan::varpart(X, env_matrix, factor_matrix)
  fractions <- as.data.frame(vp$part$indfract)
  fractions$Fraction <- rownames(fractions)
  rownames(fractions) <- NULL
  plot <- ggplot2::ggplot(fractions, ggplot2::aes(x = Fraction, y = Adj.R.squared)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Variance partitioning", y = "Adjusted R-squared", x = NULL) +
    aaa_theme()
  aaa_save_standard_result(project_dir, analysis_name, list(Fractions = fractions, Environmental_variables = data.frame(Variable = e), Experimental_factors = data.frame(Variable = f)), list(fractions = plot))
}

aaa_taxon_association_input <- function(dataset, abundance_type, variables, project_dir, analysis_name) {
  a <- aaa_align_analysis_inputs(dataset, abundance_type, "relative", variables, project_dir, analysis_name)
  # a$input$matrix is aaa_multivariate_input()'s ordination-oriented matrix:
  # it drops near-constant/low-variance taxa (a filter meant for PCA/RDA
  # numerics, not for differential testing) before ANCOM-BC2/MaAsLin2 ever
  # see them. a$raw is the same relative-abundance table without that
  # variance filter (only genuinely all-zero taxa are dropped), so
  # proportion/percentage inputs keep the full taxon set here.
  # Relative inputs carry no sequencing depth, so counts per million are the
  # closest honest approximation available.
  counts <- round(a$raw * 1e6)
  if (identical(abundance_type, "counts")) {
    prepared <- a$input$prepared
    wide <- prepared$wide
    # prepared$wide has already been rescaled to percentages by
    # aaa_prepare_amplicon_data(), so rounding it directly produced a table
    # whose sample totals were ~100 instead of the real library size: a taxon
    # at 0.4% collapsed to 0 and ANCOM-BC2's bias correction and structural-zero
    # detection, both of which model sampling depth, lost their input. Undo the
    # percentage rescaling with the library sizes captured before it happened.
    percentages <- as.matrix(wide[, prepared$sample_columns, drop = FALSE])
    sizes <- prepared$library_sizes
    if (!is.null(sizes) && length(sizes) == ncol(percentages) &&
      all(is.finite(sizes)) && all(sizes > 0)) {
      counts <- round(t(sweep(percentages, 2L, sizes / 100, "*")))
    } else {
      # Depth could not be recovered; fall back to the relative-input path
      # rather than silently handing over rounded percentages.
      counts <- round(t(sweep(
        percentages, 2L, pmax(colSums(percentages), .Machine$double.eps), "/"
      ) * 1e6))
    }
    # prepared$wide carries no rownames; recover taxon labels the same way
    # aaa_multivariate_input() does, or ANCOM-BC2/MaAsLin2 receive an
    # unlabeled feature table.
    # Coerce to character before nzchar(): a Taxonomy column read in as a
    # factor would otherwise break nzchar() with "requires a character vector".
    taxonomy_chr <- as.character(wide$Taxonomy)
    colnames(counts) <- make.unique(ifelse(!is.na(taxonomy_chr) & nzchar(taxonomy_chr), taxonomy_chr, paste0("Taxon_", seq_len(nrow(wide)))))
  }
  counts[counts < 0 | !is.finite(counts)] <- 0
  storage.mode(counts) <- "integer"
  meta <- a$metadata
  rownames(meta) <- meta$Sample_column
  meta$Sample_column <- NULL
  list(counts = counts, relative = a$raw, metadata = meta, input = a$input)
}

aaa_ancombc2_analysis <- function(dataset, abundance_type, variables, group_variable,
                                  min_prevalence = .2, alpha = .05, project_dir,
                                  analysis_name = "ANCOM_BC2") {
  aaa_check_packages(analyses = "ancombc2")
  z <- aaa_taxon_association_input(dataset, abundance_type, variables, project_dir, analysis_name)
  keep <- colMeans(z$counts > 0) >= min_prevalence
  X <- z$counts[, keep, drop = FALSE]
  if (ncol(X) < 2) stop("ANCOM-BC2 retained fewer than two taxa.")
  md <- z$metadata[rownames(X), , drop = FALSE]
  formula_vars <- intersect(variables, names(md))
  if (!length(formula_vars)) stop("No valid ANCOM-BC2 model variables.")
  # ancombc2() matches fix_formula variables against colData by literal name;
  # backtick-quoting (needed for base-R formulas) makes its own metadata
  # lookup fail, so build the formula unquoted like ANCOMBC's own examples.
  fix_formula <- paste(formula_vars, collapse = " + ")
  # ancombc2() no longer accepts separate data/meta_data matrices; current
  # ANCOMBC releases require a (Tree)SummarizedExperiment bundling both.
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(assays = list(counts = t(X)), colData = S4Vectors::DataFrame(md))
  out <- ANCOMBC::ancombc2(data = tse, assay_name = "counts", fix_formula = fix_formula, group = group_variable, p_adj_method = "BH", prv_cut = min_prevalence, lib_cut = 0, struc_zero = TRUE, neg_lb = TRUE, alpha = alpha, global = TRUE, pairwise = TRUE, dunnet = FALSE, trend = FALSE, iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE), em_control = list(tol = 1e-5, max_iter = 100), lme_control = NULL, mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), trend_control = NULL)
  res <- as.data.frame(out$res)
  res$Taxon <- rownames(res)
  rownames(res) <- NULL
  sig_cols <- grep("^diff_", names(res), value = TRUE)
  significant <- if (length(sig_cols)) res[rowSums(res[, sig_cols, drop = FALSE], na.rm = TRUE) > 0, , drop = FALSE] else res[0, , drop = FALSE]
  volcano <- tryCatch(aaa_ancombc2_volcano_plot(res, alpha), error = function(e) aaa_community_placeholder_plot("ANCOM-BC2 differential abundance", paste("Volcano plot unavailable:", conditionMessage(e))))
  aaa_save_standard_result(project_dir, analysis_name, list(Results = res, Significant = significant, Model = data.frame(Fixed_formula = fix_formula, Group = group_variable)), list(volcano = volcano))
}

aaa_maaslin_analysis <- function(dataset, abundance_type, variables, random_effects = character(),
                                 min_prevalence = .2, alpha = .05, project_dir,
                                 analysis_name = "MaAsLin2") {
  aaa_check_packages(analyses = "maaslin")
  z <- aaa_taxon_association_input(dataset, abundance_type, variables, project_dir, analysis_name)
  keep <- colMeans(z$relative > 0) >= min_prevalence
  data <- as.data.frame(z$relative[, keep, drop = FALSE])
  md <- z$metadata[rownames(data), , drop = FALSE]
  fixed <- setdiff(intersect(variables, names(md)), random_effects)
  random <- intersect(random_effects, names(md))
  if (!length(fixed)) stop("MaAsLin2 requires at least one fixed effect.")
  project <- aaa_create_project_structure(project_dir, analysis_name)
  outdir <- file.path(project$analysis, "MaAsLin2_output")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  # MaAsLin2 errors out on any fixed effect with more than two levels unless a
  # reference is supplied; default to the first level (alphabetically) so
  # analyses with 3+ groups do not fail outright.
  ref_terms <- vapply(fixed, function(v) {
    levels <- unique(stats::na.omit(md[[v]]))
    if (is.numeric(md[[v]]) || length(levels) <= 2L) {
      return(NA_character_)
    }
    paste0(v, ",", sort(as.character(levels))[1])
  }, character(1))
  ref_terms <- ref_terms[!is.na(ref_terms)]
  reference <- if (length(ref_terms)) paste(ref_terms, collapse = ";") else NULL
  fit <- Maaslin2::Maaslin2(input_data = data, input_metadata = md, output = outdir, fixed_effects = fixed, random_effects = random, reference = reference, normalization = "NONE", transform = "LOG", analysis_method = "LM", correction = "BH", standardize = TRUE, min_prevalence = min_prevalence, max_significance = alpha, plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1)
  res <- as.data.frame(fit$results)
  sig <- res[is.finite(res$qval) & res$qval <= alpha, , drop = FALSE]
  files <- c(summary = file.path(project$analysis, "MaAsLin2_summary.xlsx"))
  openxlsx::write.xlsx(list(Results = res, Significant = sig, Fixed_effects = data.frame(Variable = fixed), Random_effects = data.frame(Variable = random)), files[[1]], overwrite = TRUE)
  forest <- tryCatch(aaa_maaslin_forest_plot(res, alpha), error = function(e) aaa_community_placeholder_plot("MaAsLin2 associations", paste("Forest plot unavailable:", conditionMessage(e))))
  forest_file <- file.path(project$analysis, "MaAsLin2_forest.png")
  aaa_save_plot(forest, forest_file, width = 9, height = max(5, 0.35 * min(nrow(res), 25L) + 2.5))
  files["forest"] <- forest_file
  aaa_result(tables = list(results = res, significant = sig), plots = list(forest = forest), files = files, output_dir = project$analysis, metadata = list(fixed_effects = fixed, random_effects = random))
}
