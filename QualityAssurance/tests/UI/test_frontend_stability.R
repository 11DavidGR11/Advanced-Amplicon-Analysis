qa_register_test(
  "UI_005",
  "regression",
  "critical",
  "Frontend result navigation is bounded, on-demand and session-safe",
  function() {
    app_path <- file.path(QA_ROOT, "AAApp", "Biological", "app.R")
    text <- qa_read_app_source(QA_ROOT)

    qa_expect_true(
      grepl("result_inventory <- reactive", text, fixed = TRUE),
      "The result inventory is not cached by active run."
    )
    qa_expect_true(
      grepl("select_result_file <- function", text, fixed = TRUE),
      "Result selection is not centralized."
    )
    qa_expect_true(
      grepl("too large for a safe in-app preview", text, fixed = TRUE),
      "Large result tables are not protected from blocking previews."
    )
    qa_expect_true(
      grepl("session$onSessionEnded", text, fixed = TRUE),
      "Background processes and resource paths are not cleaned when the session ends."
    )
    qa_expect_true(
      !grepl("output$result_files_ui <- renderUI", text, fixed = TRUE),
      "The obsolete global figure gallery is still registered."
    )
    qa_expect_true(
      grepl("Results are rendered exclusively on demand", text, fixed = TRUE),
      "On-demand result rendering contract is missing."
    )

    TRUE
  }
)
