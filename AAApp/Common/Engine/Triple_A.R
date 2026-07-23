# =============================================================================
# Triple_A public API
#
# Public functions:
#   triple_a_load()
#   triple_a_config()
#   triple_a_run()
#   triple_a_run_config()
#   triple_a_list_functions()
#   triple_a_list_analyses()
#   triple_a_list_outputs()
#   triple_a_list_methods()
# =============================================================================

.triple_a_find_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "Run_Triple_A.R")) &&
      file.exists(file.path(current, "AAApp", "Common", "Engine", "Triple_A.R")) &&
      dir.exists(file.path(current, "AAApp", "Common", "Engine", "Core"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("The Triple_A distribution root could not be located.", call. = FALSE)
    }
    current <- parent
  }
}

.triple_a_resolve_root <- function(
  input_file = NULL,
  start = getwd()
) {
  loaded_root <- get0(
    "TRIPLE_A_ROOT",
    envir = globalenv(),
    inherits = TRUE,
    ifnotfound = NULL
  )

  if (!is.null(loaded_root) &&
    dir.exists(loaded_root)) {
    return(normalizePath(
      loaded_root,
      winslash = "/",
      mustWork = TRUE
    ))
  }

  candidate_starts <- character()

  if (!is.null(input_file) &&
    length(input_file) == 1L &&
    !is.na(input_file) &&
    nzchar(input_file) &&
    file.exists(input_file)) {
    candidate_starts <- c(
      candidate_starts,
      dirname(normalizePath(
        input_file,
        winslash = "/",
        mustWork = TRUE
      ))
    )
  }

  if (!is.null(start) &&
    length(start) == 1L &&
    !is.na(start) &&
    nzchar(start) &&
    dir.exists(start)) {
    candidate_starts <- c(
      candidate_starts,
      normalizePath(
        start,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  candidate_starts <- unique(
    candidate_starts
  )

  for (candidate in candidate_starts) {
    resolved <- tryCatch(
      .triple_a_find_root(candidate),
      error = function(e) NULL
    )

    if (!is.null(resolved)) {
      return(resolved)
    }
  }

  stop(
    "The Advanced_Amplicon_Analysis root could not be located. ",
    "Load Triple_A with triple_a_load(project_root = ...) or provide ",
    "an input file stored inside the project."
  )
}


triple_a_load <- function(
  project_root = NULL,
  install_missing = FALSE,
  envir = globalenv(),
  verbose = TRUE
) {
  if (is.null(project_root)) {
    project_root <- .triple_a_find_root()
  }

  project_root <- normalizePath(
    project_root,
    winslash = "/",
    mustWork = TRUE
  )

  core <- file.path(
    project_root,
    "AAApp", "Common", "Engine", "Core"
  )

  extensions <- file.path(
    project_root,
    "AAApp", "Common", "Engine", "Extensions"
  )

  files <- c(
    file.path(core, "aaa_globals.R"),
    file.path(core, "aaa_settings.R"),
    file.path(core, "aaa_dependencies.R"),
    file.path(core, "aaa_importer.R"),
    file.path(core, "aaa_data_model.R"),
    file.path(core, "aaa_utils.R"),
    file.path(core, "aaa_cache_database.R"),
    file.path(core, "aaa_projects.R"),
    file.path(core, "aaa_analysis_cache.R"),
    file.path(core, "aaa_environment_diagnostics.R"),
    file.path(core, "aaa_validation.R"),
    file.path(core, "aaa_results.R"),
    file.path(core, "aaa_background_runner.R"),
    file.path(core, "aaa_metadata.R"),
    file.path(core, "aaa_analysis_registry.R"),
    file.path(core, "aaa_plugin_system.R"),
    file.path(core, "aaa_styles.R"),
    file.path(extensions, "gene_dictionary.R"),
    file.path(core, "aaa_biological_registry_validation.R"),
    file.path(
      extensions,
      "classification_functions.R"
    ),
    file.path(
      extensions,
      "biological_functions_registry.R"
    ),
    file.path(extensions, "custom_biological_functions.R"),
    file.path(core, "aaa_functional_potential.R"),
    file.path(core, "aaa_top_abundance.R"),
    file.path(core, "aaa_extended_analyses.R"),
    file.path(core, "aaa_community_analysis.R"),
    file.path(core, "aaa_supervised_analysis.R"),
    file.path(core, "aaa_advanced_environmental_taxon.R"),
    # Consumes results produced by the functional and differential engines
    # above, so it must load after them and before the workflow that calls it.
    file.path(core, "aaa_functional_comparison.R"),
    file.path(core, "aaa_workflow.R")
  )

  missing <- files[!file.exists(files)]

  if (length(missing) > 0) {
    stop(
      "Missing Triple_A files:\n",
      paste(missing, collapse = "\n")
    )
  }

  for (file in files) {
    if (verbose) {
      message("Triple_A loading: ", basename(file))
    }
    sys.source(file, envir = envir)
  }

  if (isTRUE(install_missing)) {
    aaa_check_packages(
      install_missing = TRUE,
      include_core = TRUE
    )
  }

  required_internal_functions <- c(
    "aaa_get_defaults",
    "aaa_analysis_method_registry",
    "aaa_methodology_table",
    "aaa_create_project_structure",
    "aaa_reference_cache_connect",
    "aaa_load_reference_cache",
    "aaa_community_analysis",
    "aaa_plsda_analysis",
    "aaa_rda_analysis",
    "aaa_run_workflow"
  )

  invalid_internal_functions <- required_internal_functions[
    !vapply(
      required_internal_functions,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
    )
  ]

  if (length(invalid_internal_functions) > 0) {
    stop(
      "Triple_A failed to load internal functions: ",
      paste(invalid_internal_functions, collapse = ", ")
    )
  }

  aaa_load_plugins(project_root)

  assign(
    "TRIPLE_A_ROOT",
    project_root,
    envir = envir
  )

  invisible(project_root)
}

triple_a_config <- function(
  dataset,
  abundance_type = c("proportion", "percentage", "counts"),
  output_dir = NULL,
  progress_verbosity = c("standard", "detailed", "developer"),
  analyses = NULL, outputs = NULL, functional_functions = NULL,
  top_abundance = list(), differential_abundance = list(),
  functional_abundance = list(), community_structure = list(),
  supervised_multivariate = list(), environmental = list()
) {
  aaa_validate_dataset(dataset)
  abundance_type <- match.arg(abundance_type)
  progress_verbosity <- match.arg(progress_verbosity)
  defaults <- aaa_get_defaults()
  if (is.null(output_dir)) {
    root <- .triple_a_resolve_root(start = getwd())
    output_dir <- file.path(root, "Results")
  }
  if (is.null(analyses)) analyses <- defaults$analyses
  if (is.null(outputs)) {
    catalogue <- aaa_output_catalogue()
    outputs <- catalogue$ID[catalogue$Module %in% analyses]
  }
  aaa_validate_selections(analyses, outputs)
  structure(list(
    project = utils::modifyList(list(
      name = TRIPLE_A_NAME, short_name = TRIPLE_A_SHORT_NAME
    ), aaa_project_metadata()),
    dataset = dataset,
    input = list(abundance_type = abundance_type),
    output = list(directory = output_dir),
    execution = list(progress_verbosity = progress_verbosity),
    selection = list(
      analyses = analyses, outputs = outputs,
      functional_functions = functional_functions
    ),
    parameters = list(
      top_abundance = utils::modifyList(defaults$top_abundance, top_abundance),
      differential_abundance = utils::modifyList(defaults$differential_abundance, differential_abundance),
      functional_abundance = utils::modifyList(defaults$functional_abundance, functional_abundance),
      community_structure = utils::modifyList(defaults$community_structure, community_structure),
      supervised_multivariate = utils::modifyList(defaults$supervised_multivariate, supervised_multivariate),
      environmental = environmental
    )
  ), class = c("Triple_A_config", "list"))
}

triple_a_validate <- function(config) {
  if (!inherits(config, "Triple_A_config")) {
    stop(
      "'config' must be created with triple_a_config()."
    )
  }

  aaa_validate_dataset(config$dataset)
  treatments <- unique(as.character(config$dataset$sample_design$Treatment))
  if (length(treatments) == 0L || any(!nzchar(treatments))) {
    stop("At least one treatment name is required in the canonical sample design.")
  }

  aaa_validate_selections(
    config$selection$analyses,
    config$selection$outputs
  )

  invisible(TRUE)
}

triple_a_run_config <- function(
  config,
  verbose = TRUE,
  progress_callback = NULL
) {
  triple_a_validate(config)

  aaa_run_workflow(
    dataset = config$dataset,
    abundance_type = config$input$abundance_type,
    output_dir =
      config$output$directory,
    analyses =
      config$selection$analyses,
    outputs =
      config$selection$outputs,
    functional_functions =
      config$selection$functional_functions,
    top_abundance =
      config$parameters$top_abundance,
    differential_abundance =
      config$parameters$differential_abundance,
    functional_abundance =
      config$parameters$functional_abundance,
    community_structure =
      config$parameters$community_structure,
    supervised_multivariate =
      config$parameters$supervised_multivariate,
    environmental =
      config$parameters$environmental,
    progress_verbosity =
      config$execution$progress_verbosity %||% "standard",
    verbose = verbose,
    progress_callback = progress_callback
  )
}

triple_a_run <- function(
  dataset, abundance_type = "proportion", output_dir = NULL,
  progress_verbosity = "standard", analyses = NULL, outputs = NULL,
  functional_functions = NULL, top_abundance = list(),
  differential_abundance = list(), functional_abundance = list(),
  community_structure = list(), supervised_multivariate = list(),
  environmental = list(), verbose = TRUE, progress_callback = NULL
) {
  config <- triple_a_config(
    dataset = dataset, abundance_type = abundance_type, output_dir = output_dir,
    progress_verbosity = progress_verbosity, analyses = analyses, outputs = outputs,
    functional_functions = functional_functions, top_abundance = top_abundance,
    differential_abundance = differential_abundance,
    functional_abundance = functional_abundance,
    community_structure = community_structure,
    supervised_multivariate = supervised_multivariate, environmental = environmental
  )
  triple_a_run_config(config, verbose = verbose, progress_callback = progress_callback)
}

triple_a_dependency_status <- function(
  analyses = character(),
  input_files = character()
) {
  aaa_dependency_status(
    analyses = analyses,
    input_files = input_files
  )
}

triple_a_install_dependencies <- function(
  analyses = character(),
  input_files = character()
) {
  aaa_check_packages(
    analyses = analyses,
    input_files = input_files,
    install_missing = TRUE,
    include_core = TRUE
  )
}


triple_a_list_functions <- function() {
  aaa_registry_function_catalogue()
}

triple_a_list_analyses <- function() {
  aaa_analysis_catalogue()
}

triple_a_list_outputs <- function() {
  aaa_output_catalogue()
}

triple_a_list_methods <- function(
  analyses = aaa_analysis_catalogue()$ID,
  outputs = aaa_output_catalogue()$ID,
  parameters = list()
) {
  aaa_methodology_table(
    analyses = analyses,
    outputs = outputs,
    parameters = parameters
  )
}

print.Triple_A_config <- function(x, ...) {
  cat("Triple_A configuration\n")
  cat("----------------------\n")
  cat("Dataset schema: ", x$dataset$schema_version, "\n", sep = "")
  cat(
    "Treatments: ",
    paste(unique(x$dataset$sample_design$Treatment), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "Analyses: ",
    paste(x$selection$analyses, collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "Outputs: ",
    paste(x$selection$outputs, collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "Results: ",
    x$output$directory,
    "\n",
    sep = ""
  )
  invisible(x)
}


triple_a_list_plugins <- function() aaa_list_plugins()
triple_a_create_project <- function(...) aaa_create_project(...)
triple_a_open_project <- function(...) aaa_open_project(...)
