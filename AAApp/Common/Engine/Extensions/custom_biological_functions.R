# Declarative custom biological-function loader.

aaa_custom_function_directory <- function(project_root = getOption("triple_a_root", getwd())) {
  file.path(project_root, "Resources", "FunctionalDB", "CustomFunctions")
}

aaa_custom_character <- function(x, field, allow_empty = TRUE) {
  if (is.null(x) || length(x) == 0L) {
    if (allow_empty) {
      return(character())
    }
    stop("The field '", field, "' cannot be empty.", call. = FALSE)
  }
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  if (!is.atomic(x)) stop("The field '", field, "' must contain text values.", call. = FALSE)
  values <- trimws(as.character(x))
  values <- unique(values[nzchar(values) & !is.na(values)])
  if (!allow_empty && !length(values)) stop("The field '", field, "' cannot be empty.", call. = FALSE)
  values
}

aaa_custom_scalar <- function(x, field, allow_empty = FALSE) {
  values <- aaa_custom_character(x, field, allow_empty = allow_empty)
  if (allow_empty && !length(values)) {
    return("")
  }
  if (length(values) != 1L) stop("The field '", field, "' must contain one value.", call. = FALSE)
  values[[1L]]
}

aaa_tokenize_boolean_rule <- function(expression) {
  expression <- trimws(as.character(if (is.null(expression)) "" else expression))
  if (!nzchar(expression)) stop("A complex rule requires a Boolean expression.", call. = FALSE)
  expression <- gsub("([()])", " \\1 ", expression, perl = TRUE)
  tokens <- strsplit(trimws(gsub("[[:space:]]+", " ", expression)), " ", fixed = TRUE)[[1L]]
  tokens[nzchar(tokens)]
}

aaa_parse_boolean_rule <- function(expression, allowed_genes) {
  tokens <- aaa_tokenize_boolean_rule(expression)
  position <- 1L
  peek <- function() if (position <= length(tokens)) tokens[[position]] else NA_character_
  consume <- function() {
    value <- peek()
    position <<- position + 1L
    value
  }

  parse_primary <- NULL
  parse_not <- NULL
  parse_and <- NULL
  parse_or <- NULL

  parse_primary <- function() {
    token <- peek()
    if (is.na(token)) stop("Unexpected end of complex rule.", call. = FALSE)
    if (identical(token, "(")) {
      consume()
      node <- parse_or()
      if (!identical(peek(), ")")) stop("Missing closing parenthesis in complex rule.", call. = FALSE)
      consume()
      return(node)
    }
    if (toupper(token) %in% c("AND", "OR", "NOT") || identical(token, ")")) {
      stop("Unexpected token '", token, "' in complex rule.", call. = FALSE)
    }
    consume()
    if (!(token %in% allowed_genes)) {
      stop("The complex rule references undeclared gene '", token, "'.", call. = FALSE)
    }
    list(type = "gene", gene = token)
  }

  parse_not <- function() {
    if (identical(toupper(peek()), "NOT")) {
      consume()
      return(list(type = "not", child = parse_not()))
    }
    parse_primary()
  }

  parse_and <- function() {
    node <- parse_not()
    while (identical(toupper(peek()), "AND")) {
      consume()
      node <- list(type = "and", left = node, right = parse_not())
    }
    node
  }

  parse_or <- function() {
    node <- parse_and()
    while (identical(toupper(peek()), "OR")) {
      consume()
      node <- list(type = "or", left = node, right = parse_and())
    }
    node
  }

  tree <- parse_or()
  if (position <= length(tokens)) stop("Unexpected token '", tokens[[position]], "' in complex rule.", call. = FALSE)
  tree
}

aaa_evaluate_boolean_rule <- function(tree, values) {
  switch(tree$type,
    gene = isTRUE(values[[tree$gene]]),
    not = !aaa_evaluate_boolean_rule(tree$child, values),
    and = aaa_evaluate_boolean_rule(tree$left, values) && aaa_evaluate_boolean_rule(tree$right, values),
    or = aaa_evaluate_boolean_rule(tree$left, values) || aaa_evaluate_boolean_rule(tree$right, values),
    stop("Unsupported Boolean-rule node.", call. = FALSE)
  )
}

aaa_normalize_custom_function_definition <- function(x) {
  if (!is.list(x) || is.data.frame(x)) stop("A custom biological function must be a named JSON object.", call. = FALSE)
  required <- c("id", "display_name", "category", "description", "diagnostic", "supporting", "accessory", "rule")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing custom-function fields: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.list(x$rule) || is.data.frame(x$rule)) stop("The field 'rule' must be a JSON object.", call. = FALSE)

  minimum <- suppressWarnings(as.integer(unlist(x$rule$minimum, recursive = TRUE, use.names = FALSE)[1L]))
  if (!length(minimum) || is.na(minimum)) minimum <- 1L
  expression <- aaa_custom_scalar(x$rule$expression, "rule.expression", allow_empty = TRUE)

  list(
    id = aaa_custom_scalar(x$id, "id"),
    display_name = aaa_custom_scalar(x$display_name, "display_name"),
    category = aaa_custom_scalar(x$category, "category"),
    description = aaa_custom_scalar(x$description, "description"),
    diagnostic = aaa_custom_character(x$diagnostic, "diagnostic"),
    supporting = aaa_custom_character(x$supporting, "supporting"),
    accessory = aaa_custom_character(x$accessory, "accessory"),
    rule = list(type = aaa_custom_scalar(x$rule$type, "rule.type"), minimum = minimum, expression = expression),
    evidence_note = aaa_custom_character(x$evidence_note, "evidence_note"),
    references = aaa_custom_character(x$references, "references")
  )
}

aaa_validate_custom_function_definition <- function(x, dictionary = gene_aliases, require_dictionary = FALSE) {
  x <- aaa_normalize_custom_function_definition(x)
  if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", x$id)) {
    stop("The function ID must start with a letter and contain only letters, numbers and underscores.", call. = FALSE)
  }
  genes <- unique(c(x$diagnostic, x$supporting, x$accessory))
  if (!length(genes)) stop("At least one marker gene is required.", call. = FALSE)
  duplicated_roles <- genes[table(c(x$diagnostic, x$supporting, x$accessory))[genes] > 1L]
  if (length(duplicated_roles)) stop("Genes cannot occur in more than one role: ", paste(duplicated_roles, collapse = ", "), call. = FALSE)

  allowed_rules <- c("all_diagnostic", "minimum_diagnostic", "minimum_total", "expression")
  if (!(x$rule$type %in% allowed_rules)) stop("Unsupported custom-function rule.", call. = FALSE)
  if (x$rule$type != "expression" && (!is.finite(x$rule$minimum) || x$rule$minimum < 1L)) {
    stop("The rule minimum must be a positive integer.", call. = FALSE)
  }
  if (x$rule$type %in% c("all_diagnostic", "minimum_diagnostic") && !length(x$diagnostic)) {
    stop("Diagnostic rules require at least one diagnostic gene.", call. = FALSE)
  }
  if (x$rule$type == "minimum_diagnostic" && x$rule$minimum > length(x$diagnostic)) {
    stop("The diagnostic threshold cannot exceed the number of diagnostic genes.", call. = FALSE)
  }
  if (x$rule$type == "minimum_total" && x$rule$minimum > length(genes)) {
    stop("The total threshold cannot exceed the number of marker genes.", call. = FALSE)
  }
  if (x$rule$type == "expression") {
    tree <- aaa_parse_boolean_rule(x$rule$expression, genes)
    absent <- as.list(stats::setNames(rep(FALSE, length(genes)), genes))
    if (isTRUE(aaa_evaluate_boolean_rule(tree, absent))) {
      stop("A complex rule cannot be positive when every declared marker is absent.", call. = FALSE)
    }
  }

  if (isTRUE(require_dictionary)) {
    unknown <- setdiff(genes, names(dictionary))
    if (length(unknown)) stop("Genes absent from the Gene Dictionary: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  invisible(x)
}

aaa_custom_marker_classifier <- function(definition) {
  definition <- aaa_validate_custom_function_definition(definition)
  genes <- unique(c(definition$diagnostic, definition$supporting, definition$accessory))
  expression_tree <- if (definition$rule$type == "expression") aaa_parse_boolean_rule(definition$rule$expression, genes) else NULL
  force(definition)
  force(expression_tree)

  function(x) {
    marker_values <- x[genes]
    if (!length(marker_values) || all(is.na(marker_values))) {
      return(paste("Unknown", definition$display_name, "potential"))
    }
    marker_values[is.na(marker_values)] <- FALSE
    marker_values <- as.list(stats::setNames(as.logical(marker_values), names(marker_values)))
    diagnostic_values <- unlist(marker_values[definition$diagnostic], use.names = FALSE)
    total_present <- sum(unlist(marker_values, use.names = FALSE) %in% TRUE)
    diagnostic_present <- sum(diagnostic_values %in% TRUE)
    positive <- switch(definition$rule$type,
      all_diagnostic = diagnostic_present == length(definition$diagnostic),
      minimum_diagnostic = diagnostic_present >= definition$rule$minimum,
      minimum_total = total_present >= definition$rule$minimum,
      expression = aaa_evaluate_boolean_rule(expression_tree, marker_values),
      FALSE
    )
    if (isTRUE(positive)) {
      return(paste(definition$display_name, "potential"))
    }
    if (total_present > 0L) {
      return(paste("Partial", definition$display_name, "evidence"))
    }
    paste("No detected", definition$display_name, "potential")
  }
}

aaa_custom_definition_to_registry <- function(x) {
  x <- aaa_validate_custom_function_definition(x)
  roles <- aaa_gene_roles(x$diagnostic, x$supporting, x$accessory)
  list(
    display_name = x$display_name,
    category = "Custom functions",
    biological_category = x$category,
    description = x$description,
    gene_roles = roles,
    genes = aaa_registry_genes(roles),
    classifier = aaa_custom_marker_classifier(x),
    graph_main = paste(x$display_name, "potential"),
    analysis_name = x$id,
    pathway_selector = function(value) !grepl("^(No detected|Unknown)", value, ignore.case = TRUE),
    evidence_note = if (length(x$evidence_note)) x$evidence_note[[1L]] else "User-defined declarative marker rule. Review biological specificity before interpretation.",
    custom_definition = TRUE,
    references = x$references,
    custom_rule = x$rule
  )
}

aaa_read_custom_function_definition <- function(file, dictionary = gene_aliases) {
  raw <- jsonlite::read_json(file, simplifyVector = FALSE)
  aaa_validate_custom_function_definition(raw, dictionary = dictionary, require_dictionary = TRUE)
}

aaa_load_custom_biological_functions <- function(project_root = getOption("triple_a_root", getwd()), registry = biological_function_registry, dictionary = gene_aliases) {
  directory <- aaa_custom_function_directory(project_root)
  if (!dir.exists(directory)) {
    return(registry)
  }
  files <- sort(list.files(directory, pattern = "\\.json$", full.names = TRUE, ignore.case = TRUE))
  for (file in files) {
    definition <- tryCatch(aaa_read_custom_function_definition(file, dictionary = dictionary), error = function(e) {
      warning("Skipping invalid custom biological function '", basename(file), "': ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(definition)) next
    if (definition$id %in% names(registry)) {
      warning("Skipping duplicate biological-function ID '", definition$id, "' from ", basename(file), ".", call. = FALSE)
      next
    }
    registry[[definition$id]] <- aaa_custom_definition_to_registry(definition)
  }
  registry
}

biological_function_registry <- aaa_load_custom_biological_functions(
  project_root = getOption("triple_a_root", getwd()), registry = biological_function_registry
)
