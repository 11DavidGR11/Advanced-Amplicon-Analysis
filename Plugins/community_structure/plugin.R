# Auto-discovered Triple_A plugin: community_structure
aaa_register_plugin(list(
  id = "community_structure",
  name = "Community structure and diversity",
  category = "Community",
  description = "Canonical-data plugin for the Triple_A scientific backend.",
  requires_environmental = FALSE,
  workflow_analysis_id = "community_structure",
  parameters = list(),
  required_parameters = character(),
  test_configuration = function(dataset, context = list()) list(parameters = list(), workflow_arguments = list()),
  outputs = tryCatch(aaa_analysis_method_registry()[["community_structure"]]$outputs, error = function(e) list()),
  validator = function(context, ...) {
    aaa_validate_dataset(context$dataset)
    NULL
  },
  runner = function(context, parameters = list(), ...) {
    args <- utils::modifyList(context$workflow_arguments, list(dataset = context$dataset, analyses = "community_structure"))
    do.call(aaa_run_workflow, args)
  }
), overwrite = TRUE)
