qa_register_test(
  "PLUGIN_005", "regression", "critical",
  "Plugin execution uses an explicit configuration contract and workflow snapshots persist it",
  function() {
    plugin_system <- file.path(QA_ROOT, "AAApp", "Common", "Engine", "Core", "aaa_plugin_system.R")
    workflow <- file.path(QA_ROOT, "AAApp", "Common", "Engine", "Core", "aaa_workflow.R")
    qa_expect_files(c(plugin_system, workflow))

    plugin_src <- paste(readLines(plugin_system, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    workflow_src <- paste(readLines(workflow, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

    qa_expect_true(
      grepl("context\\$config[[:space:]]*<-[[:space:]]*merged", plugin_src),
      "aaa_run_plugin() does not provide context$config explicitly"
    )
    qa_expect_true(
      grepl("config[[:space:]]*<-[[:space:]]*list\\(", workflow_src),
      "aaa_run_workflow() does not construct an explicit snapshot configuration"
    )
    qa_expect_true(
      grepl("config[[:space:]]*=[[:space:]]*config", workflow_src),
      "Run_snapshot.rds does not persist the explicit workflow configuration"
    )
    TRUE
  }
)
