library(shiny)
library(bslib)

app_dir <- normalizePath(getOption("triple_a_multiamplicon_dir", getwd()), winslash = "/", mustWork = TRUE)
source(file.path(dirname(app_dir), "Common", "paths.R"), local = TRUE)
root <- normalizePath(getOption("triple_a_root", aaa_distribution_root(app_dir)), winslash = "/", mustWork = TRUE)
source(file.path(root, "AAApp", "Common", "Engine", "Core", "aaa_importer.R"), local = TRUE)
source(file.path(root, "AAApp", "Common", "help.R"), local = TRUE)
source(file.path(app_dir, "multiamplicon_core.R"), local = TRUE)

ui <- page_sidebar(
  title = "Amplicon Integrator (Triple A)",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#5B3F7A"),
  sidebar = sidebar(
    aaa_help_button(),
    tags$hr(),
    fileInput(
      "files", "Amplicon count tables",
      multiple = TRUE,
      accept = c(".csv", ".tsv", ".txt", ".xls", ".xlsx")
    ),
    tags$p(class = "text-muted small",
      "Select two or more tables with the same column names. Column order is aligned automatically after pressing 'Load and validate'."
    ),
    actionButton("load_files", "Load and validate", class = "btn-outline-primary w-100"),
    tags$hr(),
    uiOutput("count_column_ui"),
    checkboxInput("aggregate_duplicates", "Sum rows with identical identifier/taxonomy fields", TRUE),
    actionButton("integrate", "Integrate tables", class = "btn-primary w-100"),
    tags$hr(),
    uiOutput("download_ui")
  ),
  card(
    card_header("Validation and preview"),
    uiOutput("status"),
    tableOutput("preview")
  ),
  card(
    card_header("Integration rules"),
    tags$ul(
      tags$li("Only non-negative integer counts are accepted in selected count columns."),
      tags$li("Missing count cells are converted to zero; textual or fractional values are rejected."),
      tags$li("All files must contain the same column names; differing column order is aligned automatically."),
      tags$li("Duplicate aggregation uses the complete set of non-count columns as the key."),
      tags$li("The exported TSV uses the column order of the first input table.")
    )
  )
)

server <- function(input, output, session) {
  aaa_register_context_help(input, output, session, root, "MultiAmplicon")
  imported <- reactiveVal(NULL)
  integrated <- reactiveVal(NULL)
  message_state <- reactiveVal(list(type = "light", text = "Select at least two count tables."))

  observeEvent(input$files, {
    imported(NULL)
    integrated(NULL)
    files <- input$files
    if (is.null(files) || nrow(files) == 0L) {
      message_state(list(type = "light", text = "Select at least two count tables."))
    } else if (nrow(files) < 2L) {
      message_state(list(
        type = "warning",
        text = "One file selected. Add at least one more file, then press 'Load and validate'."
      ))
    } else {
      total_mb <- sum(files$size, na.rm = TRUE) / 1024^2
      message_state(list(
        type = "info",
        text = sprintf("%d files selected (%.1f MB). Press 'Load and validate' to read them.", nrow(files), total_mb)
      ))
    }
  }, ignoreInit = TRUE)

  observeEvent(input$load_files, {
    files <- input$files
    if (is.null(files) || nrow(files) < 2L) {
      imported(NULL)
      integrated(NULL)
      message_state(list(type = "warning", text = "Select at least two amplicon count tables."))
      return()
    }

    imported(NULL)
    integrated(NULL)
    result <- tryCatch(
      withProgress(message = "Loading amplicon tables", value = 0, {
        tables <- vector("list", nrow(files))
        names(tables) <- files$name
        for (i in seq_len(nrow(files))) {
          incProgress(1 / (nrow(files) + 1), detail = paste("Reading", files$name[[i]]))
          tables[[i]] <- aaa_import_table(files$datapath[[i]], original_name = files$name[[i]])
        }
        incProgress(1 / (nrow(files) + 1), detail = "Checking headers")
        multiamplicon_validate_headers(tables, files$name)
        tables
      }),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      message_state(list(type = "danger", text = conditionMessage(result)))
      return()
    }

    imported(result)
    message_state(list(
      type = "success",
      text = sprintf("%d compatible tables loaded. Confirm the sample-count columns.", length(result))
    ))
  }, ignoreInit = TRUE)

  output$count_column_ui <- renderUI({
    tables <- imported()
    if (is.null(tables)) {
      return(tags$div(class = "text-muted small", "Count-column selection appears after successful loading."))
    }
    headers <- names(tables[[1L]])
    candidates <- multiamplicon_guess_count_columns(tables, sample_rows = 250L)
    selectizeInput(
      "count_columns", "Sample count columns",
      choices = headers, selected = candidates, multiple = TRUE,
      options = list(plugins = list("remove_button"), maxOptions = length(headers))
    )
  })

  observeEvent(input$integrate, {
    tables <- imported()
    if (is.null(tables)) {
      message_state(list(type = "warning", text = "Load and validate the files before integration."))
      return()
    }
    result <- tryCatch(
      withProgress(message = "Integrating count tables", value = 0.2, {
        value <- multiamplicon_validate_and_combine(
          tables = tables,
          count_columns = input$count_columns,
          file_names = names(tables),
          aggregate_duplicates = isTRUE(input$aggregate_duplicates)
        )
        incProgress(0.8, detail = "Preparing preview")
        value
      }),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      integrated(NULL)
      message_state(list(type = "danger", text = conditionMessage(result)))
      return()
    }
    integrated(result)
    message_state(list(
      type = "success",
      text = sprintf(
        "Integration completed: %d input rows produced %d output rows and %d count columns.",
        sum(vapply(tables, nrow, integer(1))), nrow(result), length(input$count_columns)
      )
    ))
  }, ignoreInit = TRUE)

  output$download_ui <- renderUI({
    if (is.null(integrated())) {
      return(tags$button(
        type = "button",
        class = "btn btn-outline-secondary w-100",
        disabled = "disabled",
        "Integrate tables before downloading TSV"
      ))
    }
    downloadButton(
      "download_tsv",
      "Download integrated TSV",
      class = "btn-outline-primary w-100"
    )
  })

  output$status <- renderUI({
    state <- message_state()
    tags$div(class = paste("alert", paste0("alert-", state$type)), state$text)
  })

  output$preview <- renderTable({
    value <- integrated()
    if (is.null(value)) {
      tables <- imported()
      if (is.null(tables)) return(NULL)
      return(utils::head(tables[[1L]], 12L))
    }
    utils::head(value, 20L)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$download_tsv <- downloadHandler(
    filename = function() {
      paste0("Triple_A_MultiAmplicon_", format(Sys.Date(), "%Y%m%d"), ".tsv")
    },
    contentType = "text/tab-separated-values; charset=UTF-8",
    content = function(file) {
      value <- integrated()
      if (is.null(value)) {
        stop("No integrated table is available for download.", call. = FALSE)
      }

      utils::write.table(
        value,
        file = file,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        na = "",
        qmethod = "double",
        fileEncoding = "UTF-8"
      )

      if (!file.exists(file) || is.na(file.info(file)$size) || file.info(file)$size <= 0L) {
        stop("The TSV export could not be created.", call. = FALSE)
      }
    }
  )
  outputOptions(output, "download_tsv", suspendWhenHidden = FALSE)

}

shinyApp(ui, server)
