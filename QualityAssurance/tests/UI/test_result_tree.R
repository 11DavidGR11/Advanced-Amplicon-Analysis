qa_register_test(
  "UI_004",
  "regression",
  "high",
  "Result tree is the only selector and displays one result at a time",
  function() {
    app_path <- file.path(QA_ROOT, "AAApp", "Biological", "app.R")
    text <- qa_read_app_source(QA_ROOT)

    qa_expect_true(grepl('card_header\\("Selected result"\\)', text),
                   "Single selected-result viewer is missing")
    qa_expect_true(grepl('"Summary", "Figures", "Tables", "Reports", "Technical files"', text, fixed = TRUE),
                   "Result categories are not declared in the tree")
    qa_expect_true(!grepl('nav_panel\\("Figures"', text),
                   "Global Figures panel still displays all figures")
    qa_expect_true(!grepl('nav_panel\\(\\s*"Tables"', text, perl = TRUE),
                   "Global Tables panel still displays all tables")
    qa_expect_true(grepl("result_tree_selection", text, fixed = TRUE),
                   "Tree selection event is missing")
    TRUE
  }
)
