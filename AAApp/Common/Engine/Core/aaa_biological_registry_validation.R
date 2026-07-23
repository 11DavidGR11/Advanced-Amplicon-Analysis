# =============================================================================
# Triple_A biological-function registry validation and scoring helpers
# =============================================================================

aaa_validate_biological_registry <- function(registry = aaa_registry(), strict = FALSE) {
  required <- c(
    "display_name", "category", "description", "gene_roles", "genes",
    "classifier", "graph_main", "analysis_name", "pathway_selector"
  )
  issues <- character()
  if (!is.list(registry) || !length(registry)) issues <- c(issues, "Registry is empty")
  if (anyDuplicated(names(registry))) issues <- c(issues, "Registry IDs are duplicated")
  for (id in names(registry)) {
    x <- registry[[id]]
    missing <- setdiff(required, names(x))
    if (length(missing)) issues <- c(issues, sprintf("%s: missing %s", id, paste(missing, collapse = ", ")))
    if (!is.function(x$classifier)) issues <- c(issues, sprintf("%s: classifier is not a function", id))
    if (!length(x$genes)) issues <- c(issues, sprintf("%s: no genes", id))
    role_genes <- unique(unlist(x$gene_roles, use.names = FALSE))
    if (!all(x$genes %in% role_genes)) issues <- c(issues, sprintf("%s: genes and gene_roles differ", id))
  }
  out <- list(
    valid = !length(issues), issues = issues, n_functions = length(registry),
    n_genes = length(unique(unlist(lapply(registry, `[[`, "genes"), use.names = FALSE)))
  )
  if (strict && !out$valid) stop(paste(issues, collapse = "\n"), call. = FALSE)
  out
}
