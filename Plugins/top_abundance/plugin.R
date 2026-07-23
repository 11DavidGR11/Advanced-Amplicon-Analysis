# Auto-discovered Triple_A plugin: top_abundance
aaa_register_plugin(list(
  id = "top_abundance",
  name = "Top abundance",
  category = "Composition",
  description = "Canonical-data plugin for the Triple_A scientific backend.",
  requires_environmental = FALSE,
  workflow_analysis_id = "top_abundance",
  parameters = list(),
  required_parameters = character(),
  test_configuration = function(dataset, context = list()) list(parameters = list(), workflow_arguments = list()),
  outputs = tryCatch(aaa_analysis_method_registry()[["top_abundance"]]$outputs, error = function(e) list()),
  validator = function(context, ...) {
    aaa_validate_dataset(context$dataset)
    NULL
  },
  runner = function(context, parameters = list(), ...) {
    args <- utils::modifyList(context$workflow_arguments, list(dataset = context$dataset, analyses = "top_abundance"))
    do.call(aaa_run_workflow, args)
  }
), overwrite = TRUE)
