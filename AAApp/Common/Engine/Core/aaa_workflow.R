# =============================================================================
# Triple_A internal workflow
# =============================================================================

aaa_remove_output <- function(result, key) {
  if (is.null(result)) {
    return(result)
  }

  path <- NULL

  if (!is.null(result$files) &&
    key %in% names(result$files)) {
    path <- result$files[[key]]
  }

  if (!is.null(path) && length(path) > 0) {
    path <- as.character(path)
    path <- path[!is.na(path) & nzchar(path)]
    existing <- unname(path[vapply(
      path,
      file.exists,
      logical(1)
    )])
    if (length(existing) > 0) {
      unlink(existing, force = TRUE)
    }
  }

  if (!is.null(result$files)) {
    result$files <- result$files[
      names(result$files) != key
    ]
  }

  if (!is.null(result$plots) &&
    key %in% names(result$plots)) {
    result$plots[[key]] <- NULL
  }

  result
}

aaa_filter_selected_outputs <- function(
  results,
  selected_outputs
) {
  if (!is.null(results$functional_potential) &&
    !"functional_heatmap" %in% selected_outputs) {
    for (id in names(results$functional_potential)) {
      results$functional_potential[[id]] <-
        aaa_remove_output(
          results$functional_potential[[id]],
          "heatmap"
        )
    }
  }

  if (!is.null(results$top_abundance)) {
    mapping <- c(
      top_heatmap = "heatmap",
      top_lollipop = "lollipop",
      top_composition = "stacked_composition",
      top_distribution = "distribution"
    )

    for (selection in names(mapping)) {
      if (!selection %in% selected_outputs) {
        results$top_abundance <- aaa_remove_output(
          results$top_abundance,
          mapping[[selection]]
        )
      }
    }
  }

  if (!is.null(results$differential_abundance)) {
    if (!"volcano_plot" %in% selected_outputs) {
      volcano_keys <- grep(
        "^volcano_",
        names(results$differential_abundance$files),
        value = TRUE
      )
      for (key in volcano_keys) {
        results$differential_abundance <-
          aaa_remove_output(
            results$differential_abundance,
            key
          )
      }
      results$differential_abundance$plots$volcano <- NULL
    }

    if (!"ma_plot" %in% selected_outputs) {
      ma_keys <- grep(
        "^ma_",
        names(results$differential_abundance$files),
        value = TRUE
      )
      for (key in ma_keys) {
        results$differential_abundance <-
          aaa_remove_output(
            results$differential_abundance,
            key
          )
      }
      results$differential_abundance$plots$ma <- NULL
    }

    if (!"qq_plot" %in% selected_outputs) {
      qq_keys <- grep(
        "^qq_",
        names(results$differential_abundance$files),
        value = TRUE
      )
      for (key in qq_keys) {
        results$differential_abundance <-
          aaa_remove_output(
            results$differential_abundance,
            key
          )
      }
      results$differential_abundance$plots$qq <- NULL
    }
  }

  if (!is.null(results$functional_abundance)) {
    mapping <- c(
      functional_abundance_heatmap = "heatmap",
      functional_abundance_barplot = "barplot",
      functional_abundance_dotplot = "dotplot",
      functional_contributors = "top_contributors"
    )

    for (selection in names(mapping)) {
      if (!selection %in% selected_outputs) {
        results$functional_abundance <-
          aaa_remove_output(
            results$functional_abundance,
            mapping[[selection]]
          )
      }
    }
  }


  if (!is.null(results$community_structure)) {
    mapping <- c(
      pca_plot = "pca",
      pca_scree = "pca_scree",
      pcoa_plot = "pcoa",
      nmds_plot = "nmds",
      alpha_diversity_plot = "alpha_diversity",
      pairwise_permanova_heatmap = "pairwise_permanova",
      beta_dispersion_boxplot = "beta_dispersion",
      anosim_plot = "anosim",
      permanova_variance_plot = "permanova_variance",
      dendrogram_plot = "dendrogram"
    )

    for (selection in names(mapping)) {
      if (!selection %in% selected_outputs) {
        results$community_structure <- aaa_remove_output(
          results$community_structure,
          mapping[[selection]]
        )
      }
    }
  }
  if (!is.null(results$plsda)) {
    if (!"plsda_plot" %in% selected_outputs) {
      results$plsda <- aaa_remove_output(results$plsda, "plot")
      results$plsda$plots$plsda <- NULL
    }
    if (!"plsda_confusion_plot" %in% selected_outputs) {
      results$plsda <- aaa_remove_output(results$plsda, "confusion")
      results$plsda$plots$confusion_matrix <- NULL
    }
    if (!"plsda_vip_plot" %in% selected_outputs) {
      results$plsda <- aaa_remove_output(results$plsda, "vip")
      results$plsda$plots$vip <- NULL
    }
  }

  if (!is.null(results$splsda)) {
    if (!"splsda_plot" %in% selected_outputs) {
      results$splsda <- aaa_remove_output(results$splsda, "plot")
      results$splsda$plots$splsda <- NULL
    }
    if (!"splsda_confusion_plot" %in% selected_outputs) {
      results$splsda <- aaa_remove_output(results$splsda, "confusion")
      results$splsda$plots$confusion_matrix <- NULL
    }
    if (!"splsda_selected_features" %in% selected_outputs) {
      results$splsda <- aaa_remove_output(results$splsda, "selected")
      results$splsda$plots$selected_features <- NULL
    }
  }

  if (!is.null(results$rda)) {
    if (!"rda_plot" %in% selected_outputs) {
      results$rda <- aaa_remove_output(results$rda, "plot")
      results$rda$plots$rda <- NULL
    }
    if (!"rda_variance_plot" %in% selected_outputs) {
      results$rda <- aaa_remove_output(results$rda, "variance")
      results$rda$plots$variance <- NULL
    }
  }
  results
}

aaa_run_workflow <- function(
  dataset,
  abundance_type,
  output_dir,
  analyses,
  outputs,
  functional_functions = NULL,
  top_abundance = list(), differential_abundance = list(),
  functional_abundance = list(), community_structure = list(),
  supervised_multivariate = list(), environmental = list(),
  progress_verbosity = c("standard", "detailed", "developer"),
  verbose = TRUE, progress_callback = NULL
) {
  aaa_validate_selections(analyses, outputs)
  progress_verbosity <- match.arg(progress_verbosity)

  # Persist the complete, explicit workflow configuration.  Earlier builds
  # attempted to save an undeclared global object named `config`, which made
  # every plugin execution fail after the scientific analysis had completed.
  config <- list(
    abundance_type = abundance_type,
    analyses = analyses,
    outputs = outputs,
    functional_functions = functional_functions,
    top_abundance = top_abundance,
    differential_abundance = differential_abundance,
    functional_abundance = functional_abundance,
    community_structure = community_structure,
    supervised_multivariate = supervised_multivariate,
    environmental = environmental,
    progress_verbosity = progress_verbosity,
    verbose = verbose
  )

  aaa_validate_dataset(dataset)
  # Same declared order and same empty-group handling the analyses use, so the
  # reproducibility record matches what was actually run.
  samples_name <- aaa_treatment_levels(dataset$sample_design$Treatment)
  sample_schema <- "Triple_A_dataset_v2"
  replicate_counts <- aaa_treatment_counts(dataset$sample_design$Treatment)
  replicates <- if (length(unique(as.integer(replicate_counts))) == 1L) as.integer(replicate_counts[[1]]) else NA_integer_
  aaa_check_packages(analyses = analyses, input_files = character())

  workflow_started <- Sys.time()
  run_id <- paste0(
    format(
      workflow_started,
      "%Y%m%d_%H%M%S"
    ),
    "_",
    sprintf(
      "%03d",
      as.integer(
        (as.numeric(workflow_started) %% 1) *
          1000
      )
    )
  )

  results_root <- normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(
    results_root,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_dir <- file.path(
    results_root,
    "Runs",
    paste0(
      "Run_",
      run_id
    )
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  paths <- aaa_create_project_structure(
    project_dir = output_dir,
    analysis_name = NULL,
    results_root = results_root
  )

  logs_dir <- file.path(output_dir, "Logs")

  dir.create(
    logs_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  log_file <- file.path(
    logs_dir,
    paste0(
      "Triple_A_",
      run_id,
      "_",
      Sys.getpid(),
      ".log"
    )
  )

  aaa_append_log(
    log_file,
    "initialization",
    paste0(
      "Triple_A workflow started."
    )
  )

  aaa_write_run_status(
    output_dir,
    "running",
    "Workflow started.",
    list(run_id = run_id, analyses = analyses)
  )

  report_progress <- function(stage, detail, completed, total) {
    aaa_append_log(
      path = log_file,
      stage = stage,
      detail = detail
    )

    if (is.function(progress_callback)) {
      progress_callback(
        stage = stage,
        detail = detail,
        completed = completed,
        total = total
      )
    }

    invisible(NULL)
  }

  functional_steps <- if (
    any(c(
      "functional_potential",
      "functional_abundance"
    ) %in% analyses)
  ) {
    max(1L, length(functional_functions))
  } else {
    0L
  }

  total_steps <- functional_steps +
    as.integer("top_abundance" %in% analyses) +
    as.integer("differential_abundance" %in% analyses) +
    as.integer("functional_abundance" %in% analyses) +
    as.integer("community_structure" %in% analyses) +
    as.integer("plsda" %in% analyses) +
    as.integer("splsda" %in% analyses) +
    as.integer("rda" %in% analyses) +
    2L

  completed_steps <- 0L

  advance_progress <- function(stage, detail) {
    completed_steps <<- completed_steps + 1L
    report_progress(
      stage,
      detail,
      completed_steps,
      total_steps
    )
  }

  # Isolate each independent analysis. A failure in one of them (for example
  # ANCOM-BC2 declining to estimate a model with collinear predictors) is
  # logged and reported, and the remaining analyses, the reproducibility
  # metadata and the final report are still produced, instead of the single
  # error aborting the whole workflow. `expr` is passed unevaluated and only
  # forced inside tryCatch(), so the error is caught here.
  run_analysis_step <- function(label, expr) {
    tryCatch(
      expr,
      error = function(e) {
        detail <- conditionMessage(e)
        aaa_append_log(
          log_file,
          paste0(label, ".error"),
          paste0("Analysis skipped after an error: ", detail)
        )
        report_progress(
          stage = paste0(label, ".error"),
          detail = paste0(
            "The '", label,
            "' analysis could not be completed and was skipped: ", detail
          ),
          completed = completed_steps,
          total = total_steps
        )
        NULL
      }
    )
  }

  cache_status <- aaa_cache_status(
    results_root
  )

  report_progress(
    "initialization",
    paste0(
      "Preparing the selected analyses. Shared cache: ",
      cache_status$Value[cache_status$Metric == "Cache status"],
      "; ",
      cache_status$Value[cache_status$Metric == "Reusable cached objects"],
      " reusable object(s)."
    ),
    0L,
    total_steps
  )

  if (!exists(
    "aaa_get_defaults",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      "Triple_A internal configuration function ",
      "'aaa_get_defaults()' was not loaded."
    )
  }

  defaults <- aaa_get_defaults()

  top_cfg <- utils::modifyList(
    defaults$top_abundance,
    top_abundance
  )

  diff_cfg <- utils::modifyList(
    defaults$differential_abundance,
    differential_abundance
  )

  functional_cfg <- utils::modifyList(
    defaults$functional_abundance,
    functional_abundance
  )

  community_cfg <- utils::modifyList(
    defaults$community_structure,
    community_structure
  )

  supervised_cfg <- utils::modifyList(
    defaults$supervised_multivariate,
    supervised_multivariate
  )

  results <- list(
    functional_potential = NULL,
    top_abundance = NULL,
    differential_abundance = NULL,
    functional_abundance = NULL,
    community_structure = NULL,
    plsda = NULL,
    rda = NULL
  )

  need_functional_results <- any(
    c(
      "functional_potential",
      "functional_abundance"
    ) %in% analyses
  )

  if (need_functional_results) {
    if (is.null(functional_functions) ||
      length(functional_functions) == 0) {
      stop(
        "Select at least one biological function for ",
        "functional analyses."
      )
    }

    definitions <- aaa_registry_functional_definitions(
      functional_functions
    )

    results$functional_potential <- list()

    for (id in names(definitions)) {
      cfg <- definitions[[id]]

      report_progress(
        "functional_potential",
        paste0(
          "Analysing ",
          id,
          ". NCBI and GFF retrieval may take several minutes."
        ),
        completed_steps,
        total_steps
      )

      results$functional_potential[[id]] <-
        aaa_functional_potential(
          dataset = dataset,
          genes = cfg$genes,
          graph_main = cfg$graph_main,
          graph_note = cfg$graph_note %||% NULL,
          classification_function =
            cfg$classification_function,
          abundance_type = abundance_type,
          project_dir = output_dir,
          analysis_name = cfg$analysis_name,
          biological_function_id =
            cfg$registry_id,
          verbose = verbose,
          progress_verbosity = progress_verbosity,
          progress_callback = function(stage,
                                       detail,
                                       completed,
                                       total) {
            local_fraction <- max(
              0,
              min(
                1,
                completed / max(1, total)
              )
            )

            is_monitor_event <- identical(
              stage,
              "functional_monitor"
            )

            report_progress(
              stage = if (is_monitor_event) {
                stage
              } else {
                paste0(
                  "functional_potential.",
                  stage
                )
              },
              detail = if (is_monitor_event) {
                detail
              } else {
                paste0(
                  id,
                  ": ",
                  detail
                )
              },
              completed =
                completed_steps +
                  local_fraction,
              total = total_steps
            )
          }
        )

      advance_progress(
        "functional_potential",
        paste0("Completed functional analysis: ", id)
      )
    }
  }

  if ("top_abundance" %in% analyses) {
    report_progress(
      "top_abundance",
      "Calculating the most abundant taxa.",
      completed_steps,
      total_steps
    )

    report_progress(
      "top_abundance.input",
      "Reading the abundance table and ranking taxa.",
      completed_steps + 0.10,
      total_steps
    )

    results$top_abundance <- run_analysis_step("top_abundance", aaa_top_abundance(
      dataset = dataset,
      top_n = top_cfg$top_n,
      abundance_type = abundance_type,
      project_dir = output_dir,
      analysis_name = "Top_abundance",
      graph_main = "Most abundant microorganisms"
    ))

    report_progress(
      "top_abundance.figures",
      "Saving the selected heatmap, composition and distribution figures.",
      completed_steps + 0.90,
      total_steps
    )

    advance_progress(
      "top_abundance",
      "Completed top-abundance analysis."
    )
  }

  if ("differential_abundance" %in% analyses) {
    report_progress(
      "differential_abundance",
      "Running pairwise differential-abundance comparisons.",
      completed_steps,
      total_steps
    )

    report_progress(
      "differential_abundance.preprocessing",
      "Filtering taxa by prevalence and mean abundance.",
      completed_steps + 0.10,
      total_steps
    )

    report_progress(
      "differential_abundance.comparisons",
      "Preparing all pairwise treatment comparisons.",
      completed_steps + 0.25,
      total_steps
    )

    results$differential_abundance <-
      run_analysis_step("differential_abundance", aaa_differential_abundance(
        dataset = dataset,
        abundance_type = abundance_type,
        method = diff_cfg$method,
        paired = diff_cfg$paired,
        pseudocount = diff_cfg$pseudocount,
        min_prevalence = diff_cfg$min_prevalence,
        min_mean_abundance =
          diff_cfg$min_mean_abundance,
        alpha = diff_cfg$alpha,
        log2fc_threshold =
          diff_cfg$log2fc_threshold,
        top_n_table = diff_cfg$top_n_table,
        max_labels = diff_cfg$max_labels,
        label_only_significant =
          diff_cfg$label_only_significant,
        colour_by = diff_cfg$colour_by,
        point_size = diff_cfg$point_size,
        x_limit = diff_cfg$x_limit,
        dynamic_ylim = diff_cfg$dynamic_ylim,
        comparisons = diff_cfg$comparisons %||% NULL,
        project_dir = output_dir,
        analysis_name = "Differential_abundance"
      ))

    report_progress(
      "differential_abundance.outputs",
      "Saving statistical tables, volcano plots and MA plots.",
      completed_steps + 0.90,
      total_steps
    )

    advance_progress(
      "differential_abundance",
      "Completed differential-abundance analysis."
    )
  }

  if ("functional_abundance" %in% analyses) {
    report_progress(
      "functional_abundance",
      "Preparing potential metabolomic pathway definitions.",
      completed_steps,
      total_steps
    )

    report_progress(
      "functional_abundance.definitions",
      "Building selectors from the curated biological-function registry.",
      completed_steps + 0.15,
      total_steps
    )

    pathway_definitions <-
      aaa_registry_pathway_definitions(
        functional_functions,
        functional_results =
          results$functional_potential,
        use_analysis_reference = FALSE
      )

    report_progress(
      "functional_abundance.summary",
      "Calculating potential metabolomic pathways abundance and principal taxonomic contributors.",
      completed_steps + 0.40,
      total_steps
    )

    results$functional_abundance <-
      run_analysis_step("functional_abundance", aaa_functional_abundance(
        dataset = dataset,
        pathways = pathway_definitions,
        abundance_type = abundance_type,
        top_taxa_per_pathway =
          functional_cfg$top_taxa_per_pathway,
        project_dir = output_dir,
        analysis_name = "Potential_metabolomic_pathways_abundance",
        graph_main = "Potential metabolomic pathways abundance"
      ))

    report_progress(
      "functional_abundance.outputs",
      "Saving potential metabolomic pathways abundance plots and contributor tables.",
      completed_steps + 0.90,
      total_steps
    )

    advance_progress(
      "functional_abundance",
      "Completed potential metabolomic pathways abundance analysis."
    )
  }

  if ("community_structure" %in% analyses) {
    report_progress(
      "community_structure",
      "Preparing community ordination and diversity analyses.",
      completed_steps,
      total_steps
    )

    report_progress(
      "community_structure.ordination",
      "Calculating PCA, PCoA and NMDS coordinates.",
      completed_steps + 0.20,
      total_steps
    )

    results$community_structure <- run_analysis_step("community_structure", aaa_community_analysis(
      dataset = dataset,
      abundance_type = abundance_type,
      transformation = community_cfg$transformation,
      distance_method = community_cfg$distance_method,
      permutations = community_cfg$permutations,
      significance_alpha = community_cfg$significance_alpha,
      nmds_trymax = community_cfg$nmds_trymax,
      show_sample_labels = community_cfg$show_sample_labels,
      project_dir = output_dir,
      analysis_name = "Community_structure"
    ))

    report_progress(
      "community_structure.outputs",
      "Saving ordination plots, diversity tables and PERMANOVA results.",
      completed_steps + 0.90,
      total_steps
    )

    advance_progress(
      "community_structure",
      "Completed community-structure analysis."
    )
  }


  if ("plsda" %in% analyses) {
    report_progress("plsda", "Preparing supervised PLS-DA.", completed_steps, total_steps)
    results$plsda <- run_analysis_step("plsda", aaa_plsda_analysis(
      dataset = dataset, abundance_type = abundance_type,
      transformation = supervised_cfg$transformation,
      n_components = supervised_cfg$plsda_components,
      cv_folds = supervised_cfg$plsda_cv_folds,
      seed = supervised_cfg$plsda_seed,
      show_sample_labels = supervised_cfg$show_sample_labels,
      project_dir = output_dir, analysis_name = "PLS_DA"
    ))
    advance_progress("plsda", "Completed PLS-DA and cross-validation.")
  }

  if ("splsda" %in% analyses) {
    report_progress("splsda", "Tuning sparse PLS-DA and validating the selected microbial signature.", completed_steps, total_steps)
    results$splsda <- run_analysis_step("splsda", aaa_splsda_analysis(
      dataset = dataset, abundance_type = abundance_type,
      transformation = supervised_cfg$transformation,
      n_components = supervised_cfg$splsda_components,
      cv_folds = supervised_cfg$splsda_cv_folds,
      repeats = supervised_cfg$splsda_repeats,
      keepx_candidates = supervised_cfg$splsda_keepx,
      tune = supervised_cfg$splsda_tune,
      seed = supervised_cfg$splsda_seed,
      show_sample_labels = supervised_cfg$show_sample_labels,
      project_dir = output_dir, analysis_name = "sPLS_DA"
    ))
    advance_progress("splsda", "Completed sPLS-DA tuning, validation and feature selection.")
  }

  if ("rda" %in% analyses) {
    if (is.null(dataset$metadata)) {
      stop("RDA was selected but no metadata table is attached to the dataset.")
    }
    report_progress("rda", "Relating community composition to environmental variables.", completed_steps, total_steps)
    results$rda <- run_analysis_step("rda", aaa_rda_analysis(
      dataset = dataset, abundance_type = abundance_type,
      environmental_variables = environmental$variables,
      transformation = supervised_cfg$transformation,
      permutations = supervised_cfg$rda_permutations,
      significance_alpha = supervised_cfg$rda_alpha,
      show_sample_labels = supervised_cfg$show_sample_labels,
      project_dir = output_dir, analysis_name = "RDA"
    ))
    advance_progress("rda", "Completed RDA and permutation tests.")
  }


  if ("envfit" %in% analyses) results$envfit <- run_analysis_step("envfit", aaa_envfit_analysis(dataset, abundance_type, environmental$variables, supervised_cfg$transformation, supervised_cfg$rda_permutations, output_dir))
  if ("partial_rda" %in% analyses) results$partial_rda <- run_analysis_step("partial_rda", aaa_constrained_analysis(dataset, abundance_type, environmental$variables, environmental$experimental_factors, supervised_cfg$transformation, NULL, supervised_cfg$rda_permutations, output_dir, "Partial_RDA"))
  if ("dbrda" %in% analyses) results$dbrda <- run_analysis_step("dbrda", aaa_constrained_analysis(dataset, abundance_type, environmental$variables, character(), supervised_cfg$transformation, environmental$distance %||% "bray", supervised_cfg$rda_permutations, output_dir, "dbRDA"))
  if ("partial_dbrda" %in% analyses) results$partial_dbrda <- run_analysis_step("partial_dbrda", aaa_constrained_analysis(dataset, abundance_type, environmental$variables, environmental$experimental_factors, supervised_cfg$transformation, environmental$distance %||% "bray", supervised_cfg$rda_permutations, output_dir, "Partial_dbRDA"))
  if ("variance_partitioning" %in% analyses) results$variance_partitioning <- run_analysis_step("variance_partitioning", aaa_variance_partitioning(dataset, abundance_type, environmental$variables, environmental$experimental_factors, supervised_cfg$transformation, output_dir))
  if ("ancombc2" %in% analyses) {
    results$ancombc2 <- run_analysis_step("ancombc2", {
      group <- (environmental$experimental_factors %||% character())[1]
      if (is.na(group) || !nzchar(group)) stop("ANCOM-BC2 requires an experimental factor.")
      aaa_ancombc2_analysis(dataset, abundance_type, unique(c(environmental$experimental_factors, environmental$variables)), group, diff_cfg$min_prevalence, diff_cfg$alpha, output_dir)
    })
  }
  if ("maaslin" %in% analyses) results$maaslin <- run_analysis_step("maaslin", aaa_maaslin_analysis(dataset, abundance_type, unique(c(environmental$experimental_factors, environmental$variables)), character(), diff_cfg$min_prevalence, diff_cfg$alpha, output_dir))

  # Both of these consume results computed above rather than the dataset, so they
  # are skipped with an explanatory message when their input was not produced,
  # instead of failing the whole run.
  if ("differential_functions" %in% analyses) {
    if (is.null(results$functional_abundance)) {
      aaa_append_log(
        log_file, "differential_functions",
        "Skipped: differential functions requires the potential metabolomic pathways abundance analysis."
      )
    } else {
      results$differential_functions <- aaa_differential_functions(
        pathway_by_sample = results$functional_abundance$tables$pathway_by_sample,
        sample_design = dataset$sample_design,
        alpha = diff_cfg$alpha,
        project_dir = output_dir
      )
    }
  }

  if ("functional_enrichment" %in% analyses) {
    if (is.null(results$functional_potential) || is.null(results$differential_abundance)) {
      aaa_append_log(
        log_file, "functional_enrichment",
        "Skipped: functional enrichment requires both the functional potential and the differential abundance analyses."
      )
    } else {
      results$functional_enrichment <- aaa_functional_enrichment(
        functional_potential = results$functional_potential,
        differential_abundance = results$differential_abundance,
        alpha = diff_cfg$alpha,
        project_dir = output_dir
      )
    }
  }

  report_progress(
    "outputs",
    "Filtering and indexing selected outputs.",
    completed_steps,
    total_steps
  )

  results <- aaa_filter_selected_outputs(
    results,
    outputs
  )

  # Functional inference is a prerequisite for functional abundance.
  # When it was not selected as an independent analysis, remove its
  # intermediate analytical outputs after the functional summary has been
  # generated. Reference data remain available under Cache/.
  if (
    !"functional_potential" %in% analyses &&
      !is.null(results$functional_potential)
  ) {
    for (id in names(results$functional_potential)) {
      directory <-
        results$functional_potential[[id]]$output_dir

      if (!is.null(directory) &&
        dir.exists(directory)) {
        unlink(
          directory,
          recursive = TRUE,
          force = TRUE
        )
      }
    }

    results$functional_potential <- NULL
  }

  workflow_finished <- Sys.time()
  elapsed_seconds <- as.numeric(
    difftime(
      workflow_finished,
      workflow_started,
      units = "secs"
    )
  )

  manifest_file <- file.path(
    paths$metadata,
    "Triple_A_manifest.xlsx"
  )

  session_file <- file.path(
    paths$metadata,
    "Session_information.txt"
  )

  aaa_write_session_information(session_file)

  package_names <- aaa_required_packages()

  package_versions <- data.frame(
    Package = package_names,
    Version = vapply(
      package_names,
      function(package) {
        if (!requireNamespace(package, quietly = TRUE)) {
          return(NA_character_)
        }

        as.character(
          utils::packageVersion(package)
        )
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  input_information <- tryCatch(
    data.frame(size = NA_real_, mtime = NA),
    error = function(e) NULL
  )

  settings_table <- data.frame(
    Parameter = c(
      "software_title",
      "software_author",
      "software_institution",
      "software_license",
      "software_citation",
      "workflow_started",
      "workflow_finished",
      "elapsed_seconds",
      "input_file",
      "input_file_size_bytes",
      "input_file_modified",
      "sample_identifier",
      "treatments",
      "replicates",
      "abundance_type",
      "output_dir",
      "cache_directory",
      "r_version",
      "operating_system",
      "run_id",
      "run_directory",
      "selected_analyses",
      "selected_outputs",
      "random_seed"
    ),
    Value = c(
      aaa_project_metadata()$title,
      paste(
        aaa_project_metadata()$authors$name,
        collapse = " | "
      ),
      paste(
        unique(
          aaa_project_metadata()$authors$institution
        ),
        collapse = " | "
      ),
      aaa_project_metadata()$license,
      aaa_project_metadata()$citation,
      format(workflow_started),
      format(workflow_finished),
      format(round(elapsed_seconds, 3)),
      as.character(dataset$source$file %||% NA_character_),
      if (is.null(input_information)) {
        NA_character_
      } else {
        as.character(input_information$size[1])
      },
      if (is.null(input_information)) {
        NA_character_
      } else {
        as.character(input_information$mtime[1])
      },
      sample_schema,
      paste(samples_name, collapse = " | "),
      replicates,
      abundance_type,
      normalizePath(
        output_dir,
        winslash = "/",
        mustWork = TRUE
      ),
      paths$cache,
      R.version.string,
      paste(
        Sys.info()[c("sysname", "release", "machine")],
        collapse = " | "
      ),
      run_id,
      normalizePath(output_dir, winslash = "/", mustWork = TRUE),
      paste(analyses, collapse = " | "),
      paste(outputs, collapse = " | "),
      as.character(supervised_cfg$plsda_seed)
    ),
    stringsAsFactors = FALSE
  )

  # Each analysis stores its result under a slot named exactly like its id, so a
  # NULL slot means the isolated analysis was skipped after an error (logged
  # above by run_analysis_step()). Report that honestly instead of always
  # claiming completion.
  analysis_table <- data.frame(
    Analysis = analyses,
    Status = vapply(
      analyses,
      function(a) if (is.null(results[[a]])) "Skipped (see log)" else "Completed",
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  output_table <- data.frame(
    Output = outputs,
    stringsAsFactors = FALSE
  )

  function_table <- data.frame(
    Biological_function =
      functional_functions %||% character(),
    stringsAsFactors = FALSE
  )

  methodology_table <- aaa_methodology_table(
    analyses = analyses,
    outputs = outputs,
    parameters = list(
      differential_abundance = diff_cfg,
      community_structure = community_cfg,
      top_abundance = top_cfg,
      functional_abundance = functional_cfg,
      supervised_multivariate = supervised_cfg
    )
  )

  methodology_files <- aaa_write_methodology_files(
    methodology = methodology_table,
    metadata_dir = paths$metadata
  )

  openxlsx::write.xlsx(
    list(
      Run_metadata = settings_table,
      Analyses = analysis_table,
      Selected_outputs = output_table,
      Analysis_methods = methodology_table,
      Biological_functions = function_table,
      Authors = aaa_project_metadata()$authors,
      Package_versions = package_versions
    ),
    manifest_file,
    overwrite = TRUE
  )

  run_cache_status <- aaa_result_cache_status(results)
  run_cache_counts <- aaa_result_cache_counts(results)
  snapshot_file <- file.path(paths$metadata, "Run_snapshot.rds")
  saveRDS(list(
    dataset = dataset,
    config = config,
    results = results,
    has_biological_replicates = any(aaa_treatment_counts(dataset$sample_design$Treatment) > 1L),
    saved_at = Sys.time()
  ), snapshot_file)
  metadata_json <- file.path(paths$metadata, "Run_metadata.json")
  jsonlite::write_json(
    list(
      run_id = run_id, status = "completed",
      cache_status = run_cache_status,
      cache_hits = run_cache_counts$hits, cache_misses = run_cache_counts$misses,
      started = format(workflow_started),
      finished = format(workflow_finished), elapsed_seconds = elapsed_seconds,
      input_file = dataset$source$file %||% NA_character_,
      sample_schema = sample_schema, treatments = samples_name, analyses = analyses,
      outputs = outputs, functions = functional_functions %||% character(),
      run_directory = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
      r_version = R.version.string, operating_system = Sys.info()
    ), metadata_json,
    pretty = TRUE, auto_unbox = TRUE
  )

  figure_files <- unique(as.character(unlist(
    aaa_collect_files(results),
    recursive = TRUE,
    use.names = FALSE
  )))
  figure_files <- figure_files[
    !is.na(figure_files) & nzchar(figure_files) &
      vapply(figure_files, file.exists, logical(1)) &
      tolower(tools::file_ext(figure_files)) %in% c("png", "jpg", "jpeg")
  ]
  report_files <- aaa_write_analysis_report(
    report_dir = paths$metadata,
    run_id = run_id,
    elapsed_seconds = elapsed_seconds,
    analyses = analyses,
    methodology_table = methodology_table,
    settings_table = settings_table,
    function_table = function_table,
    package_versions = package_versions,
    figure_files = figure_files
  )
  report_file <- report_files$html

  result_manifest_file <- aaa_write_result_manifest(output_dir)

  aaa_write_run_status(
    output_dir,
    "completed",
    "Workflow completed successfully.",
    list(
      run_id = run_id, elapsed_seconds = elapsed_seconds, cache_status = run_cache_status,
      cache_hits = run_cache_counts$hits, cache_misses = run_cache_counts$misses
    )
  )

  latest_file <- file.path(results_root, "latest_run")
  writeLines(normalizePath(output_dir, winslash = "/", mustWork = TRUE), latest_file)

  results$metadata <- list(
    manifest = manifest_file,
    result_manifest = result_manifest_file,
    session_information = session_file,
    log = log_file,
    workflow_started = workflow_started,
    workflow_finished = workflow_finished,
    elapsed_seconds = elapsed_seconds,
    methods = methodology_table,
    methods_workbook = methodology_files$workbook,
    methods_text = methodology_files$text,
    metadata_json = metadata_json,
    report = report_file,
    report_pdf = report_files$pdf,
    run_id = run_id,
    run_directory = output_dir
  )

  results$manifest_file <- manifest_file
  results$session_information_file <- session_file
  results$log_file <- log_file
  results$methods_file <- methodology_files$workbook
  results$methods_text_file <- methodology_files$text

  advance_progress(
    "completed",
    paste0(
      "Workflow completed. Metadata: ",
      paths$metadata
    )
  )

  class(results) <- c("Triple_A_analysis", "list")
  results
}

# =============================================================================
# Automatic analysis report
#
# Produces a self-contained HTML report (figures embedded as base64 data URIs,
# so it can be moved or emailed as a single file) that also serves as the print
# source for an automatic PDF rendered through headless Chrome (chromote). The
# report bundles the exact run parameters and the full methodological summary,
# covering the "automatic report", "parameter export" and "methodological
# summary" requests. PDF generation is strictly best-effort: any failure leaves
# the HTML report in place and never interrupts the analysis.
# =============================================================================

aaa_html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

aaa_html_table <- function(df) {
  if (is.null(df) || !nrow(df)) {
    return("")
  }
  header <- paste0("<th>", aaa_html_escape(names(df)), "</th>", collapse = "")
  body <- vapply(seq_len(nrow(df)), function(i) {
    cells <- paste0("<td>", aaa_html_escape(unlist(df[i, , drop = TRUE])), "</td>", collapse = "")
    paste0("<tr>", cells, "</tr>")
  }, character(1))
  paste0(
    "<table><thead><tr>", header, "</tr></thead><tbody>",
    paste(body, collapse = ""), "</tbody></table>"
  )
}

aaa_embed_image_data_uri <- function(path) {
  ext <- tolower(tools::file_ext(path))
  mime <- if (ext == "png") "image/png" else if (ext %in% c("jpg", "jpeg")) "image/jpeg" else "application/octet-stream"
  paste0("data:", mime, ";base64,", base64enc::base64encode(path))
}

aaa_try_render_report_pdf <- function(html_file) {
  if (!isTRUE(getOption("triple_a_render_pdf", TRUE))) {
    return(NA_character_)
  }
  if (!requireNamespace("chromote", quietly = TRUE) ||
    !requireNamespace("base64enc", quietly = TRUE)) {
    return(NA_character_)
  }
  pdf_file <- sub("\\.html?$", ".pdf", html_file)
  ok <- tryCatch(
    {
      session <- chromote::ChromoteSession$new()
      on.exit(try(session$close(), silent = TRUE), add = TRUE)
      session$Page$navigate(paste0("file:///", normalizePath(html_file, winslash = "/", mustWork = TRUE)))
      Sys.sleep(0.5)
      result <- session$Page$printToPDF(printBackground = TRUE, preferCSSPageSize = TRUE)
      writeBin(base64enc::base64decode(result$data), pdf_file)
      file.exists(pdf_file)
    },
    error = function(e) FALSE
  )
  if (isTRUE(ok)) pdf_file else NA_character_
}

aaa_write_analysis_report <- function(report_dir, run_id, elapsed_seconds, analyses,
                                      methodology_table, settings_table,
                                      function_table, package_versions, figure_files) {
  esc <- aaa_html_escape

  css <- paste0(
    "@page{size:A4;margin:16mm}",
    "body{font-family:Arial,Helvetica,sans-serif;color:#2E2E38;font-size:12px;line-height:1.45;max-width:1100px;margin:auto;padding:24px}",
    "h1{color:#5B3F7A;margin-bottom:2px}h2{color:#5B3F7A;border-bottom:2px solid #E4DDE8;padding-bottom:4px;margin-top:26px}",
    "h3{color:#5B3F7A;margin-bottom:4px}",
    ".meta{background:#F6F3F8;padding:14px 16px;border-radius:8px}",
    "table{border-collapse:collapse;width:100%;margin:8px 0;font-size:11px}",
    "th,td{border:1px solid #ddd;padding:5px 7px;text-align:left;vertical-align:top}",
    "th{background:#F6F3F8;color:#5B3F7A}",
    ".method{page-break-inside:avoid;break-inside:avoid;margin-bottom:12px;border-left:3px solid #E4DDE8;padding-left:10px}",
    "figure{margin:10px 0;page-break-inside:avoid;break-inside:avoid}figcaption{color:#5B3F7A;font-weight:bold;margin-bottom:4px}",
    "img{max-width:100%;border:1px solid #ddd}code{background:#F6F3F8;padding:1px 4px;border-radius:3px}",
    ".hint{background:#fff7e6;border:1px solid #f0d9a8;padding:8px 12px;border-radius:6px;margin:10px 0}",
    "@media print{.hint{display:none}}"
  )

  method_html <- if (!is.null(methodology_table) && nrow(methodology_table)) {
    vapply(seq_len(nrow(methodology_table)), function(i) {
      r <- methodology_table[i, , drop = FALSE]
      paste0(
        "<div class='method'><h3>", esc(r$Analysis), " — ", esc(r$Output), "</h3>",
        "<p><b>Objective:</b> ", esc(r$Objective), "</p>",
        "<p><b>Method:</b> ", esc(r$Method), "</p>",
        "<p><b>Statistical test:</b> ", esc(r$Statistical_test), "</p>",
        "<p><b>Significance criterion:</b> ", esc(r$Significance_criterion), "</p>",
        "<p><b>Interpretation:</b> ", esc(r$Interpretation), "</p>",
        "<p><b>Run parameters:</b> <code>", esc(r$Run_parameters), "</code></p></div>"
      )
    }, character(1))
  } else {
    character()
  }

  figure_html <- if (length(figure_files)) {
    vapply(figure_files, function(f) {
      uri <- tryCatch(aaa_embed_image_data_uri(f), error = function(e) NA_character_)
      if (is.na(uri)) {
        return("")
      }
      paste0("<figure><figcaption>", esc(basename(f)), "</figcaption><img src='", uri, "'></figure>")
    }, character(1))
  } else {
    character()
  }

  html <- c(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<title>Triple_A analysis report</title>",
    paste0("<style>", css, "</style></head><body>"),
    "<h1>Advanced Amplicon Analysis — Triple_A</h1>",
    "<div class='hint'>Tip: use your browser's Print (Ctrl/Cmd+P) and choose “Save as PDF” if the automatic PDF was not generated.</div>",
    paste0(
      "<div class='meta'><b>Run:</b> ", esc(run_id),
      "<br><b>Elapsed:</b> ", round(elapsed_seconds, 1), " s",
      "<br><b>Analyses:</b> ", esc(paste(analyses, collapse = ", ")), "</div>"
    ),
    "<h2>Run parameters</h2>",
    aaa_html_table(settings_table),
    if (!is.null(function_table) && nrow(function_table)) {
      c(
        "<h2>Biological functions evaluated</h2>",
        "<p>Each functional call is reported in the Taxon_results table with its diagnostic-module completeness and a confidence tier (High / Medium / Low / Insufficient evidence).</p>",
        aaa_html_table(function_table)
      )
    } else {
      character()
    },
    "<h2>Methods and statistical criteria</h2>",
    method_html,
    "<h2>Generated figures</h2>",
    figure_html,
    "<h2>Reproducibility</h2>",
    aaa_html_table(package_versions),
    "</body></html>"
  )

  html_file <- file.path(report_dir, "Triple_A_report.html")
  writeLines(html, html_file, useBytes = TRUE)
  pdf_file <- aaa_try_render_report_pdf(html_file)
  list(html = html_file, pdf = pdf_file)
}
