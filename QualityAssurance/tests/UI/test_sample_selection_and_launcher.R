qa_register_test(
  "UI_009", "regression", "critical",
  "Distribution paths, sample selection and responsive launcher controls are wired",
  function() {
    paths <- c(
      file.path(QA_ROOT, "AAApp", "Biological", "app.R"),
      file.path(QA_ROOT, "AAApp", "FASTQ", "app.R"),
      file.path(QA_ROOT, "AAApp", "CacheManager", "app.R"),
      file.path(QA_ROOT, "AAApp", "FunctionBuilder", "app.R"),
      file.path(QA_ROOT, "AAApp", "Launcher", "app.R"),
      file.path(QA_ROOT, "Run_Triple_A.R")
    )
    qa_expect_files(paths)

    server <- paste(readLines(file.path(QA_ROOT, "AAApp", "Biological", "modules", "20_server.R"), warn = FALSE), collapse = "\n")
    launcher <- paste(readLines(file.path(QA_ROOT, "AAApp", "Launcher", "app.R"), warn = FALSE), collapse = "\n")

    required_server_tokens <- c(
      "output$selected_samples_ui",
      "input$sample_identifier",
      "ignore.case = TRUE",
      "fixed = TRUE",
      "keep_original_sample_names",
      "return(sample_columns())"
    )
    for (token in required_server_tokens) {
      qa_expect_true(grepl(token, server, fixed = TRUE), paste("Missing sample-selection contract:", token))
    }

    required_launcher_tokens <- c(
      "page_fluid(",
      "launcher_card <- function",
      'icon_name = "circle"',
      'button_label = "Open FASTQ Pipeline"',
      'icon_name = "dna"',
      "launcher-grid",
      "@media (max-width:760px)",
      "input$open_analysis",
      "input$open_fastq",
      "input$open_cache_manager",
      "input$open_function_builder",
      'AAApp", "Biological',
      'AAApp", "FASTQ',
      'AAApp", "CacheManager',
      'AAApp", "FunctionBuilder',
      "suppressWarnings(tryCatch("
    )
    for (token in required_launcher_tokens) {
      qa_expect_true(grepl(token, launcher, fixed = TRUE), paste("Missing launcher contract:", token))
    }
    qa_expect_false(grepl("page_fillable(", launcher, fixed = TRUE), "Launcher still uses a fill layout that can collapse cards.")
    qa_expect_false(grepl("layout_columns(", launcher, fixed = TRUE), "Launcher still nests fill-oriented column layouts.")
    TRUE
  }
)
