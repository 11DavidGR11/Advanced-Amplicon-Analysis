# Auto-discovered Triple_A plugin: rda
aaa_register_plugin(list(
  id = "rda",
  name = "RDA with environmental variables",
  category = "Environmental",
  description = "Canonical-data plugin for the Triple_A scientific backend.",
  requires_environmental = TRUE,
  workflow_analysis_id = "rda",
  parameters = list(),
  required_parameters = character(),
  test_configuration = function(dataset, context = list()) {
    roles <- dataset$metadata_roles
    vars <- roles$Column[roles$Role == "environmental_variable"]
    numeric_vars <- vars[vapply(dataset$metadata[vars], is.numeric, logical(1))]
    if (!length(numeric_vars)) stop("No numeric environmental variables are available for RDA testing.")
    # Keep the fixture well over the sample-to-predictor minimum.
    list(parameters = list(), workflow_arguments = list(environmental = list(variables = head(numeric_vars, 2L))))
  },
  outputs = tryCatch(aaa_analysis_method_registry()[["rda"]]$outputs, error = function(e) list()),
  validator = function(context, ...) {
    aaa_validate_dataset(context$dataset)
    NULL
  },
  runner = function(context, parameters = list(), ...) {
    args <- utils::modifyList(context$workflow_arguments, list(dataset = context$dataset, analyses = "rda"))
    do.call(aaa_run_workflow, args)
  }
), overwrite = TRUE)
