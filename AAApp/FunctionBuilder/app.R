library(shiny)
library(bslib)

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (basename(app_dir) != "FunctionBuilder") app_dir <- normalizePath(file.path(app_dir, "AAApp", "FunctionBuilder"), winslash = "/", mustWork = FALSE)
source(file.path(dirname(app_dir), "Common", "paths.R"), local = TRUE)
project_root <- normalizePath(aaa_distribution_root(app_dir), winslash = "/", mustWork = TRUE)
options(triple_a_root = project_root)
source(file.path(project_root, "AAApp", "Common", "help.R"), local = TRUE)
source(file.path(project_root, "AAApp", "Common", "Engine", "Triple_A.R"), local = TRUE)
triple_a_load(project_root, install_missing = FALSE, verbose = FALSE)

definition_dir <- aaa_custom_function_directory(project_root)
dir.create(definition_dir, recursive = TRUE, showWarnings = FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x

parse_entries <- function(x) {
  values <- unlist(strsplit(x %||% "", "[,;\\n]+"), use.names = FALSE)
  unique(trimws(values[nzchar(trimws(values))]))
}

parse_alias_updates <- function(x) {
  lines <- trimws(unlist(strsplit(x %||% "", "\\n", fixed = FALSE), use.names = FALSE))
  lines <- lines[nzchar(lines)]
  output <- list()
  for (line in lines) {
    parts <- strsplit(line, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop("Each alias line must use: gene = alias 1; alias 2", call. = FALSE)
    gene <- trimws(parts[[1L]])
    aliases <- parse_entries(parts[[2L]])
    if (!nzchar(gene) || !length(aliases)) stop("Each alias line requires a gene and at least one synonym.", call. = FALSE)
    output[[gene]] <- unique(c(output[[gene]], aliases))
  }
  aaa_normalize_gene_alias_updates(output)
}

write_definition_safely <- function(definition, destination, dictionary = gene_aliases) {
  temporary <- tempfile(pattern = ".custom_function_", tmpdir = dirname(destination), fileext = ".json")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  jsonlite::write_json(definition, temporary, pretty = TRUE, auto_unbox = TRUE, null = "null")
  aaa_read_custom_function_definition(temporary, dictionary = dictionary)
  if (file.exists(destination) && !file.copy(destination, paste0(destination, ".bak"), overwrite = TRUE)) {
    stop("The previous definition could not be backed up.", call. = FALSE)
  }
  if (!file.rename(temporary, destination) && !file.copy(temporary, destination, overwrite = TRUE)) {
    stop("The validated definition could not be installed.", call. = FALSE)
  }
  invisible(destination)
}

ui <- page_sidebar(
  title = "Biological Function Builder (Triple A)",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#5B3F7A"),
  sidebar = sidebar(
    width = 440,
    aaa_help_button(), tags$hr(),
    textInput("id", "Stable function ID", placeholder = "nitrite_respiration_custom"),
    textInput("display_name", "Display name", placeholder = "Nitrite respiration"),
    textInput("category", "Biological category", placeholder = "Nitrogen metabolism"),
    helpText("Custom functions are displayed together under 'Custom functions'; this field preserves their scientific category."),
    textAreaInput("description", "Scientific description", rows = 3),
    textAreaInput("diagnostic", "Diagnostic genes", rows = 3, placeholder = "nirK, nirS"),
    textAreaInput("supporting", "Supporting genes", rows = 2, placeholder = "norB"),
    textAreaInput("accessory", "Accessory genes", rows = 2, placeholder = "nosZ"),
    textAreaInput("aliases", "New synonyms for the Gene Dictionary", rows = 4, placeholder = "newGene = accepted product name; alternative annotation"),
    selectInput("rule_type", "Decision rule", choices = c(
      "All diagnostic genes" = "all_diagnostic",
      "Minimum diagnostic genes" = "minimum_diagnostic",
      "Minimum genes across all roles" = "minimum_total",
      "Complex Boolean expression" = "expression"
    )),
    conditionalPanel("input.rule_type != 'expression'", numericInput("minimum", "Minimum markers", value = 1, min = 1, step = 1)),
    conditionalPanel("input.rule_type == 'expression'", textAreaInput("expression", "Boolean rule", rows = 3, placeholder = "(nirK OR nirS) AND norB"), helpText("Use declared gene names with AND, OR, NOT and parentheses.")),
    textAreaInput("evidence_note", "Interpretation and limitations", rows = 3),
    textAreaInput("references", "References (one per line)", rows = 3),
    actionButton("validate", "Validate definition", class = "btn-outline-primary w-100"),
    actionButton("save", "Save and activate", class = "btn-primary w-100 mt-2")
  ),
  navset_card_tab(
    nav_panel("Definition", card(card_header("Validation and preview"), uiOutput("status"), verbatimTextOutput("preview"))),
    nav_panel("Gene Dictionary", card(card_header("Dictionary status"), tableOutput("gene_status"))),
    nav_panel("Installed custom functions", card(tableOutput("installed")))
  )
)

server <- function(input, output, session) {
  aaa_register_context_help(input, output, session, project_root, "FunctionBuilder")
  status <- reactiveVal(tags$div(class = "alert alert-light", "Complete the definition and validate it before saving."))

  alias_updates <- reactive(parse_alias_updates(input$aliases))
  effective_dictionary <- reactive(aaa_merge_gene_aliases(gene_aliases, alias_updates()))
  definition <- reactive(aaa_normalize_custom_function_definition(list(
    id = trimws(input$id), display_name = trimws(input$display_name), category = trimws(input$category),
    description = trimws(input$description), diagnostic = parse_entries(input$diagnostic),
    supporting = parse_entries(input$supporting), accessory = parse_entries(input$accessory),
    rule = list(type = input$rule_type, minimum = as.integer(input$minimum %||% 1L), expression = trimws(input$expression %||% "")),
    evidence_note = trimws(input$evidence_note), references = parse_entries(input$references)
  )))

  validate_all <- function() {
    updates <- alias_updates()
    dictionary <- effective_dictionary()
    value <- aaa_validate_custom_function_definition(definition(), dictionary = dictionary, require_dictionary = TRUE)
    list(definition = value, aliases = updates, dictionary = dictionary)
  }

  output$status <- renderUI(status())
  output$preview <- renderText(tryCatch(jsonlite::toJSON(definition(), pretty = TRUE, auto_unbox = TRUE), error = function(e) paste("Definition incomplete:", conditionMessage(e))))
  output$gene_status <- renderTable({
    genes <- unique(c(parse_entries(input$diagnostic), parse_entries(input$supporting), parse_entries(input$accessory)))
    updates <- tryCatch(alias_updates(), error = function(e) list())
    data.frame(
      Gene = genes,
      Dictionary_status = ifelse(genes %in% names(gene_aliases), "Existing", ifelse(genes %in% names(updates), "New entry", "Missing")),
      Alias_count = vapply(genes, function(g) length(unique(c(gene_aliases[[g]], updates[[g]]))), integer(1)),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE)

  installed_data <- reactiveVal(data.frame())
  refresh_installed <- function() {
    files <- sort(list.files(definition_dir, pattern = "\\.json$", full.names = TRUE, ignore.case = TRUE))
    rows <- lapply(files, function(file) {
      value <- tryCatch(aaa_read_custom_function_definition(file), error = function(e) NULL)
      if (is.null(value)) return(NULL)
      data.frame(ID = value$id, Name = value$display_name, Interface_group = "Custom functions", Biological_category = value$category, Rule = value$rule$type, File = basename(file), stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    installed_data(if (length(rows)) do.call(rbind, rows) else data.frame())
  }
  refresh_installed()
  output$installed <- renderTable(installed_data(), striped = TRUE)

  observeEvent(input$validate, ignoreInit = TRUE, {
    tryCatch({
      result <- validate_all()
      new_genes <- setdiff(unique(c(result$definition$diagnostic, result$definition$supporting, result$definition$accessory)), names(gene_aliases))
      status(tags$div(class = "alert alert-success", paste0("Definition is valid. ", length(new_genes), " new dictionary entr", if (length(new_genes) == 1L) "y" else "ies", " will be saved.")))
    }, error = function(e) status(tags$div(class = "alert alert-danger", conditionMessage(e))))
  })

  observeEvent(input$save, ignoreInit = TRUE, {
    tryCatch({
      result <- validate_all()
      existing <- biological_function_registry[[result$definition$id]]
      if (!is.null(existing) && !isTRUE(existing$custom_definition)) stop("This ID belongs to a built-in biological function.", call. = FALSE)
      destination <- file.path(definition_dir, paste0(result$definition$id, ".json"))
      write_definition_safely(result$definition, destination, dictionary = result$dictionary)
      if (length(result$aliases)) aaa_save_custom_gene_aliases(result$aliases, project_root)
      refresh_installed()
      status(tags$div(class = "alert alert-success", "Function and dictionary updates installed. Restart Biological Analysis to load the updated registry."))
    }, error = function(e) status(tags$div(class = "alert alert-danger", conditionMessage(e))))
  })
}

shinyApp(ui, server)
