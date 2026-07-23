# Triple_A plugin system -----------------------------------------------
# Plugins are declarative wrappers around the existing scientific engine.
.aaa_plugin_registry <- new.env(parent = emptyenv())

aaa_plugin_validate_definition <- function(plugin) {
  required <- c(
    "id", "name", "category", "description",
    "requires_environmental", "parameters", "outputs", "runner"
  )
  missing <- setdiff(required, names(plugin))
  if (length(missing)) stop("Invalid plugin; missing fields: ", paste(missing, collapse = ", "))
  if (!grepl("^[a-z][a-z0-9_]*$", plugin$id)) stop("Invalid plugin id: ", plugin$id)
  if (!is.function(plugin$runner)) stop("Plugin runner must be a function: ", plugin$id)
  if (!is.null(plugin$validator) && !is.function(plugin$validator)) stop("Plugin validator must be a function: ", plugin$id)
  if (!is.null(plugin$test_configuration) && !is.function(plugin$test_configuration)) {
    stop("Plugin test_configuration must be a function: ", plugin$id)
  }
  if (!is.null(plugin$required_parameters) && !is.character(plugin$required_parameters)) {
    stop("Plugin required_parameters must be a character vector: ", plugin$id)
  }
  TRUE
}

aaa_register_plugin <- function(plugin, overwrite = FALSE) {
  aaa_plugin_validate_definition(plugin)
  if (exists(plugin$id, envir = .aaa_plugin_registry, inherits = FALSE) && !overwrite) {
    stop("Plugin already registered: ", plugin$id)
  }
  assign(plugin$id, plugin, envir = .aaa_plugin_registry)
  invisible(plugin)
}

aaa_get_plugin <- function(id) {
  if (!exists(id, envir = .aaa_plugin_registry, inherits = FALSE)) stop("Unknown plugin: ", id)
  get(id, envir = .aaa_plugin_registry, inherits = FALSE)
}

aaa_list_plugins <- function() {
  ids <- sort(ls(.aaa_plugin_registry, all.names = TRUE))
  if (!length(ids)) {
    return(data.frame())
  }
  do.call(rbind, lapply(ids, function(id) {
    p <- aaa_get_plugin(id)
    data.frame(
      ID = p$id, Name = p$name, Category = p$category,
      Requires_environmental = isTRUE(p$requires_environmental),
      Description = p$description, stringsAsFactors = FALSE
    )
  }))
}

aaa_validate_plugin <- function(id, context = list()) {
  p <- aaa_get_plugin(id)
  base <- data.frame(
    Check = paste0("Plugin: ", p$name), Status = "Valid",
    Message = paste0("Plugin '", p$id, "' is registered."), Blocking = FALSE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(p$requires_environmental) && is.null(context$environmental_file)) {
    base <- data.frame(
      Check = paste0("Plugin: ", p$name), Status = "Error",
      Message = "Environmental metadata are required.", Blocking = TRUE,
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(p$validator)) {
    extra <- p$validator(context)
    if (!is.null(extra)) base <- rbind(base, extra)
  }
  base
}

# `%||%` is defined once in aaa_globals.R.


aaa_default_test_configuration <- function(plugin) {
  if (is.character(plugin) && length(plugin) == 1L) {
    plugin <- aaa_get_plugin(plugin)
  }
  if (!is.list(plugin) || is.null(plugin$id)) {
    stop("A valid plugin definition is required.")
  }

  list(
    parameters = plugin$parameters %||% list(),
    workflow_arguments = list()
  )
}


aaa_plugin_test_configuration <- function(id, dataset, context = list()) {
  p <- aaa_get_plugin(id)
  default_cfg <- aaa_default_test_configuration(p)

  cfg <- if (is.function(p$test_configuration)) {
    p$test_configuration(dataset = dataset, context = context)
  } else {
    NULL
  }

  # Simple plugins may return NULL and inherit the safe default structure.
  if (is.null(cfg)) cfg <- default_cfg
  if (!is.list(cfg)) {
    stop("Plugin test configuration must return a list: ", id)
  }

  parameters <- cfg$parameters
  workflow_arguments <- cfg$workflow_arguments

  if (is.null(parameters)) parameters <- list()
  if (is.null(workflow_arguments)) workflow_arguments <- list()

  if (!is.list(parameters)) {
    stop("Plugin test configuration field 'parameters' must be a list: ", id)
  }
  if (!is.list(workflow_arguments)) {
    stop("Plugin test configuration field 'workflow_arguments' must be a list: ", id)
  }

  cfg$parameters <- utils::modifyList(default_cfg$parameters, parameters)
  cfg$workflow_arguments <- workflow_arguments
  cfg
}

aaa_run_plugin <- function(id, context, parameters = list()) {
  p <- aaa_get_plugin(id)
  if (!is.list(context)) stop("Plugin context must be a list: ", id)
  if (!is.list(parameters)) stop("Plugin parameters must be a list: ", id)

  merged <- utils::modifyList(p$parameters %||% list(), parameters)

  # Keep one explicit configuration contract for plugins and QA callers.
  # No runner should depend on a free/global variable named `config`.
  context$parameters <- merged
  context$params <- merged
  context$config <- merged
  context$workflow_arguments <- context$workflow_arguments %||% list()

  validation <- aaa_validate_plugin(id, context)
  if (any(validation$Blocking)) {
    stop(
      "Plugin validation failed for ", id, ": ",
      paste(validation$Message[validation$Blocking], collapse = "; ")
    )
  }

  missing_required <- setdiff(p$required_parameters %||% character(), names(merged))
  if (length(missing_required)) {
    stop("Plugin parameters are missing: ", paste(missing_required, collapse = ", "))
  }

  p$runner(context = context, parameters = merged)
}

aaa_load_plugins <- function(project_root = get0("TRIPLE_A_ROOT", envir = globalenv(), inherits = TRUE)) {
  plugin_dir <- file.path(project_root, "Plugins")
  files <- sort(list.files(plugin_dir, pattern = "plugin\\.R$", recursive = TRUE, full.names = TRUE))
  for (f in files) sys.source(f, envir = globalenv())
  invisible(aaa_list_plugins())
}
