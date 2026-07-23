# =============================================================================
# Triple_A supervised and environmentally constrained multivariate analyses
# =============================================================================

aaa_multivariate_input <- function(
  dataset, abundance_type, transformation, project_dir, analysis_name
) {
  prepared <- aaa_prepare_amplicon_data(
    dataset = dataset, abundance_type = abundance_type,
    project_dir = project_dir, analysis_name = analysis_name,
    filter_genus = FALSE
  )
  matrix <- t(as.matrix(prepared$wide[, prepared$sample_columns, drop = FALSE]))
  # Coerce to character before nzchar(): this feeds every multivariate
  # analysis (RDA, PLS-DA, sPLS-DA, community structure...) through this one
  # helper, so a Taxonomy column read in as a factor would otherwise break
  # nzchar() here with "'nzchar()' requires a character vector".
  taxonomy_chr <- as.character(prepared$wide$Taxonomy)
  colnames(matrix) <- make.unique(ifelse(
    !is.na(taxonomy_chr) & nzchar(taxonomy_chr),
    taxonomy_chr, paste0("Taxon_", seq_len(nrow(prepared$wide)))
  ))
  matrix[!is.finite(matrix)] <- 0
  matrix[matrix < 0] <- 0
  matrix <- matrix[, colSums(matrix) > 0, drop = FALSE]
  totals <- rowSums(matrix)
  if (any(totals <= 0)) stop("All samples must have positive total abundance.")
  relative <- sweep(matrix, 1, totals, "/")
  transformed <- switch(transformation,
    hellinger = sqrt(relative),
    relative = relative,
    log1p = log1p(relative * 100)
  )
  variable_variance <- vapply(
    seq_len(ncol(transformed)),
    function(index) {
      stats::var(
        transformed[, index],
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  keep_variable <- is.finite(variable_variance) &
    variable_variance >
      sqrt(.Machine$double.eps)

  transformed <- transformed[
    ,
    keep_variable,
    drop = FALSE
  ]

  if (ncol(transformed) < 2L) {
    stop(
      paste(
        "At least two informative taxa are required after",
        "removing constant or non-finite variables."
      )
    )
  }
  metadata <- as.data.frame(prepared$sample_map, stringsAsFactors = FALSE)
  rownames(metadata) <- metadata$Sample_column
  metadata <- metadata[rownames(transformed), , drop = FALSE]
  metadata$Treatment <- factor(metadata$Treatment, levels = prepared$samples_name)
  list(matrix = transformed, relative = relative, metadata = metadata, prepared = prepared)
}

aaa_prepare_plsda_matrix <- function(X, scale_predictors = TRUE) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  if (nrow(X) < 3L) {
    stop("PLS-DA requires at least three samples.")
  }

  if (ncol(X) < 2L) {
    stop("PLS-DA requires at least two informative taxa.")
  }

  # Treat every non-finite value consistently before imputation.
  X[!is.finite(X)] <- NA_real_

  all_missing <- vapply(
    seq_len(ncol(X)),
    function(index) all(is.na(X[, index])),
    logical(1)
  )
  X <- X[, !all_missing, drop = FALSE]

  if (ncol(X) < 2L) {
    stop("PLS-DA has fewer than two finite taxa after preprocessing.")
  }

  # Deterministic median imputation for isolated missing/non-finite values.
  for (index in seq_len(ncol(X))) {
    missing <- is.na(X[, index])
    if (any(missing)) {
      replacement <- stats::median(X[, index], na.rm = TRUE)
      if (!is.finite(replacement)) {
        stop(
          sprintf(
            "PLS-DA could not impute predictor '%s'.",
            colnames(X)[index] %||% paste0("column_", index)
          )
        )
      }
      X[missing, index] <- replacement
    }
  }

  # Remove constant and numerically near-constant predictors.  The relative
  # tolerance prevents internal scaling from producing Inf/NaN values.
  centers <- colMeans(X)
  scales <- apply(X, 2L, stats::sd)
  magnitudes <- pmax(1, apply(abs(X), 2L, max))
  scale_tolerance <- sqrt(.Machine$double.eps) * magnitudes
  keep <- is.finite(scales) & scales > scale_tolerance
  X <- X[, keep, drop = FALSE]
  centers <- centers[keep]
  scales <- scales[keep]

  if (ncol(X) < 2L) {
    stop(
      paste(
        "PLS-DA has fewer than two informative taxa after",
        "removing constant or numerically unstable variables."
      )
    )
  }

  if (isTRUE(scale_predictors)) {
    X <- sweep(X, 2L, centers, FUN = "-")
    X <- sweep(X, 2L, scales, FUN = "/")
  }

  storage.mode(X) <- "double"
  if (any(!is.finite(X))) {
    stop(
      "PLS-DA input still contains missing or infinite values after preprocessing."
    )
  }

  attr(X, "plsda_preprocessing") <- list(
    imputation = "column_median",
    centered = isTRUE(scale_predictors),
    scaled = isTRUE(scale_predictors),
    retained_predictors = ncol(X)
  )
  X
}


aaa_plsda_component_matrix <- function(x, component, n_samples, n_responses,
                                       label = "PLS-DA response") {
  component <- as.integer(component)
  arr <- as.array(x)
  dims <- dim(arr)

  if (is.null(dims)) {
    if (length(arr) != n_samples * n_responses) {
      stop(label, " has an unexpected length.", call. = FALSE)
    }
    return(matrix(as.numeric(arr), nrow = n_samples, ncol = n_responses))
  }

  if (length(dims) == 2L) {
    if (identical(as.integer(dims), c(as.integer(n_samples), as.integer(n_responses)))) {
      return(matrix(as.numeric(arr),
        nrow = n_samples, ncol = n_responses,
        dimnames = dimnames(arr)
      ))
    }
    if (identical(as.integer(dims), c(as.integer(n_responses), as.integer(n_samples)))) {
      return(t(matrix(as.numeric(arr),
        nrow = n_responses, ncol = n_samples,
        dimnames = dimnames(arr)
      )))
    }
    if (n_responses == 1L && dims[[1L]] == n_samples && dims[[2L]] >= component) {
      return(matrix(arr[, component], nrow = n_samples, ncol = 1L))
    }
  }

  if (length(dims) >= 3L && dims[[1L]] == n_samples && dims[[2L]] == n_responses) {
    if (dims[[3L]] < component) {
      stop(label, " does not contain component ", component, ".", call. = FALSE)
    }
    index <- lapply(dims, seq_len)
    index[[1L]] <- seq_len(n_samples)
    index[[2L]] <- seq_len(n_responses)
    index[[3L]] <- component
    if (length(index) > 3L) {
      for (j in 4:length(index)) index[[j]] <- 1L
    }
    value <- do.call(`[`, c(list(arr), index, list(drop = TRUE)))
    return(matrix(as.numeric(value),
      nrow = n_samples, ncol = n_responses,
      dimnames = list(dimnames(arr)[[1L]], dimnames(arr)[[2L]])
    ))
  }

  stop(
    label, " has incompatible dimensions: ",
    paste(dims, collapse = " x "),
    "; expected ", n_samples, " x ", n_responses,
    " (optionally followed by a component dimension).",
    call. = FALSE
  )
}


aaa_plsda_analysis <- function(
  dataset, abundance_type = c("proportion", "percentage", "counts"),
  transformation = c("hellinger", "relative", "log1p"),
  n_components = 2, cv_folds = 5, seed = 123,
  permutation_repetitions = 99,
  show_sample_labels = FALSE, project_dir,
  analysis_name = "PLS_DA"
) {
  transformation <- match.arg(transformation)
  aaa_check_packages(analyses = "plsda")
  input <- aaa_multivariate_input(dataset, abundance_type, transformation, project_dir, analysis_name)
  X <- aaa_prepare_plsda_matrix(input$matrix)
  group <- droplevels(input$metadata$Treatment)
  if (any(is.na(group))) stop("PLS-DA treatment labels contain missing values.")
  if (nlevels(group) < 2) stop("PLS-DA requires at least two treatment groups.")
  if (any(table(group) < 2)) stop("PLS-DA requires at least two samples per group.")

  folds <- max(2L, min(as.integer(cv_folds), nrow(X), min(table(group))))
  # PLS-DA is not limited to (number of classes - 1) components: that bound
  # belongs to LDA. Capping it there collapsed every two-group analysis to a
  # single component, which left the score plot as a flat line. The real limits
  # are the rank of X and the largest number of components the cross-validation
  # segments can support (pls rejects ncomp > n - largest segment - 1).
  max_cv_components <- nrow(X) - ceiling(nrow(X) / folds) - 1L
  n_components <- max(1L, min(
    as.integer(n_components), nrow(X) - 1L, ncol(X), max_cv_components
  ))
  Y <- stats::model.matrix(~ group - 1)
  colnames(Y) <- levels(group)
  set.seed(as.integer(seed))
  model <- tryCatch(
    pls::plsr(Y ~ X,
      ncomp = n_components, validation = "CV",
      segments = folds, scale = FALSE, center = TRUE
    ),
    error = function(e) stop("PLS-DA model fitting failed after data validation: ", conditionMessage(e), call. = FALSE)
  )

  score_matrix <- pls::scores(model)
  scores <- data.frame(
    Sample_column = rownames(X), Component1 = score_matrix[, 1],
    Component2 = if (n_components >= 2) score_matrix[, 2] else 0,
    Treatment = group, stringsAsFactors = FALSE
  )
  predictions <- aaa_plsda_component_matrix(
    model$validation$pred,
    component = n_components,
    n_samples = nrow(Y),
    n_responses = ncol(Y),
    label = "PLS-DA cross-validation predictions"
  )
  rownames(predictions) <- rownames(Y)
  colnames(predictions) <- colnames(Y)
  predicted_class <- colnames(Y)[max.col(predictions, ties.method = "first")]
  performance <- data.frame(
    Sample_column = rownames(X), Observed = as.character(group),
    Predicted = predicted_class, Correct = predicted_class == as.character(group),
    stringsAsFactors = FALSE
  )
  accuracy <- mean(performance$Correct)

  confusion <- as.data.frame.matrix(table(Observed = performance$Observed, Predicted = performance$Predicted))
  confusion <- data.frame(Observed = rownames(confusion), confusion, row.names = NULL, check.names = FALSE)
  classes <- levels(group)
  class_metrics <- do.call(rbind, lapply(classes, function(cl) {
    obs <- performance$Observed == cl
    pred <- performance$Predicted == cl
    tp <- sum(obs & pred)
    fn <- sum(obs & !pred)
    fp <- sum(!obs & pred)
    tn <- sum(!obs & !pred)
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    f1 <- if (is.finite(precision + sensitivity) && (precision + sensitivity) > 0) 2 * precision * sensitivity / (precision + sensitivity) else NA_real_
    data.frame(
      Class = cl, Support = sum(obs), Sensitivity = sensitivity, Specificity = specificity,
      Precision = precision, F1 = f1, stringsAsFactors = FALSE
    )
  }))
  balanced_accuracy <- mean(class_metrics$Sensitivity, na.rm = TRUE)

  loadings_matrix <- as.matrix(pls::loadings(model))
  loadings <- data.frame(
    Taxon = rownames(loadings_matrix), Component1 = loadings_matrix[, 1],
    Component2 = if (n_components >= 2) loadings_matrix[, 2] else 0, stringsAsFactors = FALSE
  )
  loadings$Importance <- sqrt(loadings$Component1^2 + loadings$Component2^2)

  # Standard VIP approximation from X weights and Y variance explained.
  W <- as.matrix(model$loading.weights)[, seq_len(n_components), drop = FALSE]
  Tscore <- as.matrix(pls::scores(model))[, seq_len(n_components), drop = FALSE]
  Q <- as.matrix(model$Yloadings)[, seq_len(n_components), drop = FALSE]
  ssy <- vapply(seq_len(n_components), function(h) sum(Tscore[, h]^2) * sum(Q[, h]^2), numeric(1))
  vip <- sqrt(ncol(X) * rowSums(sweep(W^2, 2, ssy, "*")) / max(sum(ssy), .Machine$double.eps))
  vip_table <- data.frame(Taxon = rownames(W), VIP = vip, stringsAsFactors = FALSE)
  vip_table <- vip_table[order(vip_table$VIP, decreasing = TRUE), , drop = FALSE]
  loadings <- loadings[order(loadings$Importance, decreasing = TRUE), , drop = FALSE]

  # Model fit and predictive summary.
  # stats::fitted() on a pls model ignores ncomp and returns the full
  # sample x response x component array, so the component to extract has to be
  # requested explicitly. Taking component 1 made R2Y describe a one-component
  # model while Q2 described the n-component one.
  fitted_y <- aaa_plsda_component_matrix(
    stats::fitted(model),
    component = n_components,
    n_samples = nrow(Y),
    n_responses = ncol(Y),
    label = "PLS-DA fitted responses"
  )
  rownames(fitted_y) <- rownames(Y)
  colnames(fitted_y) <- colnames(Y)
  # Both statistics are referred to the same total sum of squares. colMeans()
  # has to be expanded row-wise: subtracting the bare vector recycles it down
  # the columns and only happens to be correct when every class has the same
  # number of replicates.
  total_sum_of_squares <- max(
    sum((Y - matrix(colMeans(Y), nrow(Y), ncol(Y), byrow = TRUE))^2),
    .Machine$double.eps
  )
  r2y <- 1 - sum((Y - fitted_y)^2) / total_sum_of_squares
  press <- sum((Y - predictions)^2)
  q2 <- 1 - press / total_sum_of_squares

  permutation_repetitions <- max(0L, as.integer(permutation_repetitions))
  permutation_accuracy <- numeric(permutation_repetitions)
  if (permutation_repetitions > 0L) {
    set.seed(as.integer(seed) + 1L)
    for (i in seq_len(permutation_repetitions)) {
      gp <- sample(group)
      yp <- stats::model.matrix(~ gp - 1)
      colnames(yp) <- levels(group)
      pm <- try(pls::plsr(yp ~ X, ncomp = n_components, validation = "CV", segments = folds, scale = FALSE, center = TRUE), silent = TRUE)
      if (inherits(pm, "try-error")) {
        permutation_accuracy[i] <- NA_real_
        next
      }
      pp <- aaa_plsda_component_matrix(
        pm$validation$pred,
        component = n_components,
        n_samples = nrow(yp),
        n_responses = ncol(yp),
        label = "PLS-DA permutation predictions"
      )
      colnames(pp) <- colnames(yp)
      permutation_accuracy[i] <- mean(
        colnames(yp)[max.col(pp, ties.method = "first")] == as.character(gp)
      )
    }
  }
  permutation_p <- if (any(is.finite(permutation_accuracy))) {
    (1 + sum(permutation_accuracy >= accuracy, na.rm = TRUE)) / (1 + sum(is.finite(permutation_accuracy)))
  } else {
    NA_real_
  }
  permutation_table <- data.frame(Iteration = seq_along(permutation_accuracy), Accuracy = permutation_accuracy)

  summary_table <- data.frame(
    Metric = c(
      "Samples", "Classes", "Components", "CV folds", "Cross-validated accuracy",
      "Balanced accuracy", "R2Y", "Q2", "Permutation repetitions", "Permutation p-value"
    ),
    Value = c(
      nrow(X), nlevels(group), n_components, folds, accuracy, balanced_accuracy, r2y, q2,
      permutation_repetitions, permutation_p
    ), stringsAsFactors = FALSE
  )

  plot <- aaa_ordination_plot(scores, "Component1", "Component2", "PLS component 1", "PLS component 2",
    sprintf("PLS-DA (cross-validated accuracy = %.1f%%)", 100 * accuracy),
    show_labels = show_sample_labels
  )
  confusion_long <- tidyr::pivot_longer(confusion, -Observed, names_to = "Predicted", values_to = "Count")
  confusion_plot <- ggplot2::ggplot(confusion_long, ggplot2::aes(Predicted, Observed, fill = Count)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = Count), fontface = "bold") +
    ggplot2::labs(title = "PLS-DA cross-validation confusion matrix", x = "Predicted", y = "Observed") +
    aaa_theme()
  vip_top <- utils::head(vip_table, 20)
  vip_plot <- ggplot2::ggplot(vip_top, ggplot2::aes(x = stats::reorder(Taxon, VIP), y = VIP)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::labs(title = "PLS-DA variable importance", x = NULL, y = "VIP") +
    aaa_theme()

  project <- aaa_create_project_structure(project_dir, analysis_name)
  files <- c(
    plot = file.path(project$analysis, "PLS_DA.png"),
    confusion = file.path(project$analysis, "PLS_DA_confusion_matrix.png"),
    vip = file.path(project$analysis, "PLS_DA_VIP.png"),
    summary = file.path(project$analysis, "PLS_DA_summary.xlsx")
  )
  aaa_save_plot(plot, files[["plot"]], 8, 6)
  aaa_save_plot(confusion_plot, files[["confusion"]], 7, 6)
  aaa_save_plot(vip_plot, files[["vip"]], max(8, aaa_flipped_axis_plot_width(vip_top$Taxon)), 7)
  openxlsx::write.xlsx(list(
    Summary = summary_table, Scores = scores, Cross_validation = performance,
    Confusion_matrix = confusion, Classification_metrics = class_metrics, VIP_scores = vip_table,
    Taxon_loadings = loadings, Permutation_test = permutation_table
  ), files[["summary"]], overwrite = TRUE)
  aaa_result(
    tables = list(
      summary = summary_table, scores = scores, performance = performance,
      confusion_matrix = confusion, classification_metrics = class_metrics, vip_scores = vip_table,
      loadings = loadings, permutation_test = permutation_table
    ),
    plots = list(plsda = plot, confusion_matrix = confusion_plot, vip = vip_plot), files = files,
    output_dir = project$analysis, metadata = list(
      accuracy = accuracy, balanced_accuracy = balanced_accuracy,
      r2y = r2y, q2 = q2, permutation_p = permutation_p
    )
  )
}

aaa_splsda_analysis <- function(
  dataset, abundance_type = c("proportion", "percentage", "counts"),
  transformation = c("hellinger", "relative", "log1p"),
  n_components = 2, cv_folds = 5, repeats = 10,
  keepx_candidates = c(5, 10, 20), tune = TRUE, seed = 123,
  show_sample_labels = FALSE, project_dir, analysis_name = "sPLS_DA"
) {
  transformation <- match.arg(transformation)
  aaa_check_packages(analyses = "splsda")
  input <- aaa_multivariate_input(dataset, abundance_type, transformation, project_dir, analysis_name)
  X <- aaa_prepare_plsda_matrix(input$matrix)
  group <- droplevels(input$metadata$Treatment)
  if (any(is.na(group))) stop("sPLS-DA treatment labels contain missing values.")
  if (nlevels(group) < 2L) stop("sPLS-DA requires at least two treatment groups.")
  if (any(table(group) < 2L)) stop("sPLS-DA requires at least two samples per group.")

  # Same reasoning as in aaa_plsda_analysis(): the (classes - 1) bound is an LDA
  # limit, not a PLS one, and it flattened every two-group score plot.
  n_components <- max(1L, min(as.integer(n_components), nrow(X) - 1L, ncol(X)))
  folds <- max(2L, min(as.integer(cv_folds), min(table(group))))
  repeats <- max(1L, as.integer(repeats))
  keepx_candidates <- sort(unique(as.integer(keepx_candidates)))
  keepx_candidates <- keepx_candidates[is.finite(keepx_candidates) & keepx_candidates >= 1L & keepx_candidates <= ncol(X)]
  if (!length(keepx_candidates)) keepx_candidates <- unique(pmin(ncol(X), c(5L, 10L, 20L)))
  keepX <- rep(min(keepx_candidates), n_components)

  set.seed(as.integer(seed))
  tuning <- NULL
  if (isTRUE(tune) && length(keepx_candidates) > 1L) {
    tuning <- mixOmics::tune.splsda(
      X = X, Y = group, ncomp = n_components,
      test.keepX = keepx_candidates,
      validation = "Mfold", folds = folds, nrepeat = repeats,
      dist = "centroids.dist", measure = "BER",
      progressBar = FALSE
    )
    selected <- tuning$choice.keepX
    if (!is.null(selected) && length(selected)) keepX[seq_len(min(length(selected), n_components))] <- as.integer(selected[seq_len(min(length(selected), n_components))])
  }

  model <- mixOmics::splsda(X = X, Y = group, ncomp = n_components, keepX = keepX)
  set.seed(as.integer(seed) + 1L)
  performance_object <- mixOmics::perf(
    model,
    validation = "Mfold", folds = folds, nrepeat = repeats,
    dist = "centroids.dist", progressBar = FALSE
  )
  # predict() is a base R generic (stats::predict); mixOmics registers an S3
  # method for it but does not re-export the generic itself under its own
  # namespace, so mixOmics::predict(...) fails with "not an exported object"
  # even though dispatch to mixOmics's method works fine through the generic.
  prediction <- stats::predict(model, X)
  predicted <- prediction$class$centroids.dist[, n_components]
  performance <- data.frame(
    Sample_column = rownames(X), Observed = as.character(group),
    Predicted = as.character(predicted), Correct = as.character(predicted) == as.character(group),
    stringsAsFactors = FALSE
  )
  confusion <- as.data.frame.matrix(table(Observed = performance$Observed, Predicted = performance$Predicted))
  confusion <- data.frame(Observed = rownames(confusion), confusion, row.names = NULL, check.names = FALSE)
  accuracy <- mean(performance$Correct)
  class_recall <- vapply(levels(group), function(cl) mean(performance$Predicted[performance$Observed == cl] == cl), numeric(1))
  balanced_accuracy <- mean(class_recall, na.rm = TRUE)

  variates <- model$variates$X
  scores <- data.frame(
    Sample_column = rownames(X), Component1 = variates[, 1],
    Component2 = if (n_components >= 2L) variates[, 2] else 0,
    Treatment = group, stringsAsFactors = FALSE
  )
  selected_rows <- do.call(rbind, lapply(seq_len(n_components), function(comp) {
    vals <- mixOmics::selectVar(model, comp = comp)$value
    if (is.null(vals) || !nrow(vals)) {
      return(NULL)
    }
    data.frame(Taxon = rownames(vals), Component = comp, Loading = as.numeric(vals[, 1]), stringsAsFactors = FALSE)
  }))
  if (is.null(selected_rows)) selected_rows <- data.frame(Taxon = character(), Component = integer(), Loading = numeric())
  selected_rows$Absolute_loading <- abs(selected_rows$Loading)
  selected_rows <- selected_rows[order(selected_rows$Component, -selected_rows$Absolute_loading), , drop = FALSE]

  stability <- aggregate(Component ~ Taxon, data = selected_rows, FUN = length)
  names(stability)[2] <- "Components_selected"
  stability$Selection_frequency <- stability$Components_selected / n_components
  stability <- stability[order(-stability$Selection_frequency, stability$Taxon), , drop = FALSE]

  ber <- tryCatch(performance_object$error.rate$BER$centroids.dist[n_components, "mean"], error = function(e) NA_real_)
  summary_table <- data.frame(
    Metric = c("Samples", "Classes", "Components", "CV folds", "CV repeats", "Selected taxa per component", "Fitted accuracy (descriptive only)", "Balanced fitted accuracy (descriptive only)", "Cross-validated BER"),
    Value = c(nrow(X), nlevels(group), n_components, folds, repeats, paste(keepX, collapse = ", "), accuracy, balanced_accuracy, ber),
    stringsAsFactors = FALSE
  )
  score_plot <- aaa_ordination_plot(scores, "Component1", "Component2", "sPLS component 1", "sPLS component 2",
    sprintf("sPLS-DA (selected signature; CV BER = %s)", ifelse(is.finite(ber), sprintf("%.3f", ber), "NA")),
    show_labels = show_sample_labels
  )
  confusion_long <- tidyr::pivot_longer(confusion, -Observed, names_to = "Predicted", values_to = "Count")
  confusion_plot <- ggplot2::ggplot(confusion_long, ggplot2::aes(Predicted, Observed, fill = Count)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = Count), fontface = "bold") +
    ggplot2::labs(title = "sPLS-DA fitted confusion matrix", x = "Predicted", y = "Observed") +
    aaa_theme()
  feature_top <- utils::head(selected_rows, 30)
  feature_plot <- ggplot2::ggplot(feature_top, ggplot2::aes(x = stats::reorder(Taxon, Absolute_loading), y = Absolute_loading, fill = factor(Component))) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "sPLS-DA selected microbial signature", x = NULL, y = "Absolute loading", fill = "Component") +
    aaa_theme()

  project <- aaa_create_project_structure(project_dir, analysis_name)
  files <- c(
    plot = file.path(project$analysis, "sPLS_DA.png"), confusion = file.path(project$analysis, "sPLS_DA_confusion_matrix.png"),
    selected = file.path(project$analysis, "sPLS_DA_selected_features.png"), summary = file.path(project$analysis, "sPLS_DA_summary.xlsx")
  )
  aaa_save_plot(score_plot, files[["plot"]], 8, 6)
  aaa_save_plot(confusion_plot, files[["confusion"]], 7, 6)
  aaa_save_plot(feature_plot, files[["selected"]], max(8, aaa_flipped_axis_plot_width(feature_top$Taxon)), 8)
  openxlsx::write.xlsx(list(
    Summary = summary_table, Scores = scores, Predictions = performance, Confusion_matrix = confusion,
    Selected_features = selected_rows, Component_recurrence = stability
  ), files[["summary"]], overwrite = TRUE)
  aaa_result(
    tables = list(
      summary = summary_table, scores = scores, performance = performance, confusion_matrix = confusion,
      selected_features = selected_rows, component_recurrence = stability
    ), plots = list(splsda = score_plot, confusion_matrix = confusion_plot, selected_features = feature_plot),
    files = files, output_dir = project$analysis, metadata = list(accuracy = accuracy, balanced_accuracy = balanced_accuracy, ber = ber, keepX = keepX)
  )
}

aaa_rda_analysis <- function(
  dataset, abundance_type = c("proportion", "percentage", "counts"),
  environmental_variables,
  transformation = c("hellinger", "relative", "log1p"),
  permutations = 999, significance_alpha = 0.05,
  show_sample_labels = FALSE, project_dir,
  analysis_name = "RDA"
) {
  transformation <- match.arg(transformation)
  aaa_check_packages(analyses = "rda")
  input <- aaa_multivariate_input(
    dataset, abundance_type, transformation, project_dir, analysis_name
  )
  environmental <- aaa_dataset_environment(dataset, environmental_variables)

  if (is.null(environmental_variables) || length(environmental_variables) == 0) {
    stop("Select at least one environmental variable for RDA.")
  }

  missing_variables <- setdiff(environmental_variables, names(environmental))
  if (length(missing_variables)) {
    stop(
      "Environmental variables not found: ",
      paste(missing_variables, collapse = ", ")
    )
  }

  env <- environmental[
    ,
    c("Sample_column", environmental_variables),
    drop = FALSE
  ]
  env <- env[
    match(rownames(input$matrix), env$Sample_column), ,
    drop = FALSE
  ]

  if (anyNA(env$Sample_column)) {
    stop("Some abundance-table samples are missing from environmental metadata.")
  }

  for (variable in environmental_variables) {
    env[[variable]] <- suppressWarnings(as.numeric(env[[variable]]))
    if (all(is.na(env[[variable]]))) {
      stop("RDA variable is not numeric: ", variable)
    }
  }

  complete <- stats::complete.cases(
    env[, environmental_variables, drop = FALSE]
  )
  X <- input$matrix[complete, , drop = FALSE]
  env_complete <- env[complete, , drop = FALSE]

  if (nrow(X) < 4L) {
    stop("RDA requires at least four complete samples.")
  }

  raw_predictors <- as.data.frame(
    env_complete[, environmental_variables, drop = FALSE],
    stringsAsFactors = FALSE
  )

  predictor_sd <- vapply(
    raw_predictors,
    function(x) stats::sd(x, na.rm = TRUE),
    numeric(1)
  )
  informative <- is.finite(predictor_sd) & predictor_sd > 0
  candidate_variables <- names(raw_predictors)[informative]
  removed_constant <- names(raw_predictors)[!informative]

  if (length(candidate_variables) == 0L) {
    stop("No informative environmental variables remain after removing constant variables.")
  }

  # Keep at least one residual degree of freedom. Variables are considered in
  # the order selected by the user; linearly dependent variables are skipped.
  max_predictors <- max(1L, nrow(X) - 2L)
  used_variables <- character()
  removed_collinear <- character()
  removed_sample_limit <- character()

  for (variable in candidate_variables) {
    if (length(used_variables) >= max_predictors) {
      removed_sample_limit <- c(removed_sample_limit, variable)
      next
    }

    trial_variables <- c(used_variables, variable)
    trial_data <- raw_predictors[, trial_variables, drop = FALSE]
    design <- stats::model.matrix(~., data = trial_data)

    if (qr(design)$rank == ncol(design)) {
      used_variables <- trial_variables
    } else {
      removed_collinear <- c(removed_collinear, variable)
    }
  }

  if (length(used_variables) == 0L) {
    stop(
      "No independent environmental variables remain after RDA predictor filtering."
    )
  }

  predictors <- as.data.frame(
    scale(raw_predictors[, used_variables, drop = FALSE]),
    stringsAsFactors = FALSE
  )
  rownames(predictors) <- rownames(X)

  removed_variables <- unique(c(
    removed_constant,
    removed_collinear,
    removed_sample_limit
  ))

  if (length(removed_variables)) {
    warning(
      "RDA automatically excluded environmental variables: ",
      paste(removed_variables, collapse = ", "),
      ". See Predictor_selection in the RDA summary workbook for details.",
      call. = FALSE
    )
  }

  selection_table <- data.frame(
    Variable = environmental_variables,
    Status = ifelse(
      environmental_variables %in% used_variables,
      "Used",
      "Excluded"
    ),
    Reason = vapply(
      environmental_variables,
      function(variable) {
        if (variable %in% used_variables) {
          return("Included in RDA model")
        }
        if (variable %in% removed_constant) {
          return("Constant or zero variance")
        }
        if (variable %in% removed_collinear) {
          return("Linearly dependent on variables already included")
        }
        if (variable %in% removed_sample_limit) {
          return("Excluded to preserve residual degrees of freedom")
        }
        "Excluded"
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  model <- vegan::rda(X ~ ., data = predictors)

  raw_site_scores <- as.data.frame(
    vegan::scores(model, display = "sites", choices = 1:2)
  )
  site_axis_names <- names(raw_site_scores)
  if (ncol(raw_site_scores) < 2L) {
    raw_site_scores$Axis2 <- 0
    site_axis_names <- c(site_axis_names, "Axis2")
  }
  site_scores <- raw_site_scores[, 1:2, drop = FALSE]
  names(site_scores) <- c("Axis1", "Axis2")
  site_scores$Sample_column <- rownames(X)
  site_scores <- dplyr::left_join(
    site_scores,
    input$metadata,
    by = "Sample_column"
  )

  raw_biplot_scores <- as.data.frame(
    vegan::scores(model, display = "bp", choices = 1:2)
  )
  if (ncol(raw_biplot_scores) < 2L) raw_biplot_scores$Axis2 <- 0
  biplot_scores <- raw_biplot_scores[, 1:2, drop = FALSE]
  names(biplot_scores) <- c("Axis1", "Axis2")
  biplot_scores$Variable <- rownames(raw_biplot_scores)

  eigenvalues <- model$CCA$eig
  explained <- if (length(eigenvalues) && sum(eigenvalues) > 0) {
    100 * eigenvalues / sum(eigenvalues)
  } else {
    numeric()
  }

  total_inertia <- model$tot.chi
  unconstrained_eigenvalues <- model$CA$eig %||% numeric()
  axis_label <- function(axis_name, axis_number) {
    if (grepl("^RDA", axis_name)) {
      axis_index <- suppressWarnings(as.integer(sub("^RDA", "", axis_name)))
      percent <- if (is.finite(axis_index) && axis_index <= length(explained)) explained[axis_index] else NA_real_
      if (is.finite(percent)) {
        return(sprintf("%s (%.1f%% constrained)", axis_name, percent))
      }
      return(paste0(axis_name, " (constrained)"))
    }
    if (grepl("^PC", axis_name)) {
      axis_index <- suppressWarnings(as.integer(sub("^PC", "", axis_name)))
      percent <- if (is.finite(axis_index) && axis_index <= length(unconstrained_eigenvalues) &&
        is.finite(total_inertia) && total_inertia > 0) {
        100 * unconstrained_eigenvalues[axis_index] / total_inertia
      } else {
        NA_real_
      }
      if (is.finite(percent)) {
        return(sprintf("%s (%.1f%% total; unconstrained)", axis_name, percent))
      }
      return(paste0(axis_name, " (unconstrained)"))
    }
    paste0("Ordination axis ", axis_number)
  }
  x_axis_label <- axis_label(site_axis_names[1], 1L)
  y_axis_label <- axis_label(site_axis_names[2], 2L)

  explained_table <- data.frame(
    Axis = paste0("RDA", seq_along(explained)),
    Percent = explained,
    Cumulative = cumsum(explained)
  )

  safe_anova <- function(...) {
    tryCatch(
      {
        result <- as.data.frame(vegan::anova.cca(...))
        data.frame(
          Term = rownames(result),
          result,
          row.names = NULL,
          check.names = FALSE
        )
      },
      error = function(e) {
        data.frame(
          Term = "Not available",
          Message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }

  permutation_overall <- safe_anova(
    model,
    permutations = as.integer(permutations)
  )
  permutation_terms <- safe_anova(
    model,
    by = "term",
    permutations = as.integer(permutations)
  )
  permutation_axes <- safe_anova(
    model,
    by = "axis",
    permutations = as.integer(permutations)
  )

  vif_values <- tryCatch(
    vegan::vif.cca(model),
    error = function(e) {
      stats::setNames(rep(NA_real_, length(used_variables)), used_variables)
    }
  )
  vif_table <- data.frame(
    Variable = names(vif_values),
    VIF = as.numeric(vif_values),
    stringsAsFactors = FALSE
  )

  adj <- tryCatch(
    vegan::RsquareAdj(model),
    error = function(e) list(r.squared = NA_real_, adj.r.squared = NA_real_)
  )

  pcol <- intersect(c("Pr(>F)", "Pr..F."), names(permutation_terms))
  significant_terms <- if (length(pcol)) {
    sum(
      permutation_terms[[pcol[1]]] <= significance_alpha,
      na.rm = TRUE
    )
  } else {
    NA_integer_
  }

  # Share of the TOTAL inertia captured by the constrained component. The
  # per-axis percentages in `explained` are relative to the constrained part
  # only, so summing them always yields 100 and says nothing about how much of
  # the community variation the environment actually explains.
  constrained_percent <- if (is.finite(model$tot.chi) && model$tot.chi > 0 &&
    !is.null(model$CCA$tot.chi)) {
    100 * model$CCA$tot.chi / model$tot.chi
  } else {
    NA_real_
  }

  summary_table <- data.frame(
    Metric = c(
      "Samples",
      "Environmental variables requested",
      "Environmental variables used",
      "Environmental variables excluded",
      "Residual degrees of freedom",
      "Constrained variance (% of total inertia)",
      "R2",
      "Adjusted R2",
      "Significant terms",
      "Maximum VIF",
      "Permutations"
    ),
    Value = c(
      nrow(X),
      length(environmental_variables),
      length(used_variables),
      length(removed_variables),
      nrow(X) - ncol(stats::model.matrix(~., data = predictors)),
      constrained_percent,
      adj$r.squared,
      adj$adj.r.squared,
      significant_terms,
      if (any(is.finite(vif_table$VIF))) {
        max(vif_table$VIF, na.rm = TRUE)
      } else {
        NA_real_
      },
      as.integer(permutations)
    ),
    stringsAsFactors = FALSE
  )

  # Scale environmental vectors to fit inside the site-score panel. The
  # analytical scores exported below remain unchanged; only this plotting copy
  # is rescaled. Axis-wise limits avoid a single extreme vector dominating the
  # figure and provide the same practical result as vegan biplot arrow scaling
  # without depending on a base-graphics device.
  biplot_plot_scores <- biplot_scores
  site_abs <- vapply(
    site_scores[, c("Axis1", "Axis2"), drop = FALSE],
    function(values) max(abs(values), na.rm = TRUE),
    numeric(1)
  )
  vector_abs <- vapply(
    biplot_plot_scores[, c("Axis1", "Axis2"), drop = FALSE],
    function(values) max(abs(values), na.rm = TRUE),
    numeric(1)
  )
  valid_axes <- is.finite(site_abs) & site_abs > 0 &
    is.finite(vector_abs) & vector_abs > 0
  arrow_multiplier <- if (any(valid_axes)) {
    min(0.72 * site_abs[valid_axes] / vector_abs[valid_axes])
  } else {
    1
  }
  if (!is.finite(arrow_multiplier) || arrow_multiplier <= 0) {
    arrow_multiplier <- 1
  }
  biplot_plot_scores$Axis1 <- biplot_plot_scores$Axis1 * arrow_multiplier
  biplot_plot_scores$Axis2 <- biplot_plot_scores$Axis2 * arrow_multiplier

  overall_p_column <- intersect(c("Pr(>F)", "Pr..F."), names(permutation_overall))
  overall_p <- if (length(overall_p_column) > 0L) {
    suppressWarnings(as.numeric(permutation_overall[[overall_p_column[[1L]]]][[1L]]))
  } else {
    NA_real_
  }
  rda_annotations <- c(
    if (is.finite(constrained_percent)) {
      sprintf("Constrained variance: %.1f%%", constrained_percent)
    },
    if (is.finite(overall_p)) {
      paste0("Permutation p = ", format.pval(overall_p, digits = 3, eps = 0.001))
    },
    paste0("Predictors retained: ", length(used_variables))
  )

  plot <- aaa_ordination_plot(
    site_scores,
    "Axis1",
    "Axis2",
    x_axis_label,
    y_axis_label,
    "Redundancy Analysis (RDA)",
    show_labels = show_sample_labels,
    subtitle_extra = rda_annotations,
    fixed_aspect = TRUE,
    legend_title = "Treatment"
  ) +
    ggplot2::geom_segment(
      data = biplot_plot_scores,
      ggplot2::aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(.18, "cm")),
      colour = "grey30",
      linewidth = .7
    ) +
    ggrepel::geom_text_repel(
      data = biplot_plot_scores,
      ggplot2::aes(x = Axis1, y = Axis2, label = Variable),
      inherit.aes = FALSE,
      colour = "grey20",
      fontface = "bold",
      size = 3.4,
      seed = 1,
      max.overlaps = Inf,
      box.padding = 0.45,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.colour = "grey55"
    )

  variance_plot <- ggplot2::ggplot(
    utils::head(explained_table, 10),
    ggplot2::aes(
      x = stats::reorder(Axis, seq_along(Axis)),
      y = Percent,
      group = 1
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = "RDA constrained variance",
      x = "Axis",
      y = "Variance explained (%)"
    ) +
    aaa_theme()

  project <- aaa_create_project_structure(project_dir, analysis_name)
  files <- c(
    plot = file.path(project$analysis, "RDA.png"),
    variance = file.path(project$analysis, "RDA_variance.png"),
    summary = file.path(project$analysis, "RDA_summary.xlsx")
  )

  aaa_save_plot(plot, files[["plot"]], 8, 6)
  aaa_save_plot(variance_plot, files[["variance"]], 8, 5.5)

  openxlsx::write.xlsx(
    list(
      Summary = summary_table,
      Predictor_selection = selection_table,
      Environmental_data = env_complete,
      Scaled_predictors_used = predictors,
      Site_scores = site_scores,
      Environmental_vectors = biplot_scores,
      Overall_test = permutation_overall,
      Term_tests = permutation_terms,
      Axis_tests = permutation_axes,
      Explained_variance = explained_table,
      VIF = vif_table
    ),
    files[["summary"]],
    overwrite = TRUE
  )

  aaa_result(
    tables = list(
      summary = summary_table,
      predictor_selection = selection_table,
      environmental_data = env_complete,
      scaled_predictors_used = predictors,
      scores = site_scores,
      vectors = biplot_scores,
      overall_test = permutation_overall,
      term_tests = permutation_terms,
      axis_tests = permutation_axes,
      explained_variance = explained_table,
      vif = vif_table
    ),
    plots = list(rda = plot, variance = variance_plot),
    files = files,
    output_dir = project$analysis,
    metadata = list(
      r2 = adj$r.squared,
      adjusted_r2 = adj$adj.r.squared,
      requested_variables = environmental_variables,
      used_variables = used_variables,
      excluded_variables = removed_variables
    )
  )
}
