# Auto-discovered Triple_A plugin: potential_metabolomic_pathways_abundance
aaa_register_plugin(list(
  id = "potential_metabolomic_pathways_abundance",
  name = "Potential metabolomic pathways abundance",
  category = "Functional",
  description = "Canonical-data plugin for the Triple_A scientific backend.",
  requires_environmental = FALSE,
  workflow_analysis_id = "functional_abundance",
  parameters = list(),
  required_parameters = character(),
  test_configuration = function(dataset, context = list()) {
    ids <- names(aaa_registry())
    if (!length(ids)) stop("No biological functions are registered for plugin testing.")
    list(parameters = list(), workflow_arguments = list(functional_functions = ids[[1]]))
  },
  outputs = tryCatch(aaa_analysis_method_registry()[["functional_abundance"]]$outputs, error = function(e) list()),
  validator = function(context, ...) {
    aaa_validate_dataset(context$dataset)
    NULL
  },
  runner = function(context, parameters = list(), ...) {
    args <- utils::modifyList(context$workflow_arguments, list(dataset = context$dataset, analyses = "functional_abundance"))
    do.call(aaa_run_workflow, args)
  }
), overwrite = TRUE)
