qa_register_test(
  "UI_008", "regression", "high",
  "Biological-function interface provides searchable thematic blocks and bulk controls",
  function() {
    ui <- paste(readLines(file.path(QA_ROOT,"AAApp","Biological","modules","10_ui.R"),warn=FALSE),collapse="\n")
    server <- paste(readLines(file.path(QA_ROOT,"AAApp","Biological","modules","20_server.R"),warn=FALSE),collapse="\n")
    qa_expect_true(all(vapply(c("function_search","functions_all","functions_none","function_preset","function_group_accordion"),grepl,logical(1),x=ui,fixed=TRUE)),"function controls missing")
    qa_expect_true(all(vapply(c("update_function_groups","function_selection_summary"),grepl,logical(1),x=server,fixed=TRUE)),"function handlers missing")
    TRUE
  }
)
