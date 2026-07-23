qa_register_test(
  "PLUGIN_003", "release", "high",
  "Every plugin provides a valid self-contained test configuration",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Triple_A.R"))
    triple_a_load(QA_ROOT, verbose = FALSE)

    dataset <- qa_plugin_fixture_dataset()
    ids <- qa_plugin_ids()
    qa_expect_true(length(ids) > 0L, "No plugins were discovered for PLUGIN_003.")

    failures <- character()

    for (id in ids) {
      result <- tryCatch(
        {
          cfg <- aaa_plugin_test_configuration(
            id,
            dataset,
            context = list(testing = TRUE)
          )

          if (!is.list(cfg)) {
            stop("configuration is not a list")
          }
          if (!is.list(cfg$parameters)) {
            stop("parameters is not a list")
          }
          if (!is.list(cfg$workflow_arguments)) {
            stop("workflow_arguments is not a list")
          }
          TRUE
        },
        error = function(e) {
          failures <<- c(
            failures,
            paste0(id, ": ", conditionMessage(e))
          )
          FALSE
        }
      )

      invisible(result)
    }

    qa_expect_true(
      length(failures) == 0L,
      paste("Invalid plugin test configurations:", paste(failures, collapse = " | "))
    )

    TRUE
  }
)
