source(file.path(QA_ROOT, "QualityAssurance", "framework", "plugin_validation.R"))

qa_register_test(
  "PLUGIN_002",
  "release",
  "high",
  "Every discovered plugin validates and delegates to the canonical workflow",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Triple_A.R"))
    triple_a_load(QA_ROOT, verbose = FALSE)

    dataset <- qa_plugin_fixture_dataset()
    ids <- qa_plugin_ids()
    qa_expect_true(length(ids) > 0L, "No plugins were discovered")

    original_workflow <- get("aaa_run_workflow", envir = globalenv(), inherits = TRUE)
    calls <- list()
    assign(
      "aaa_run_workflow",
      function(...) {
        args <- list(...)
        key <- as.character(args$analyses %||% "unknown")[[1L]]
        calls[[key]] <<- args
        structure(
          list(status = "PASS", analyses = args$analyses, output_dir = args$output_dir),
          class = "Triple_A_QA_plugin_result"
        )
      },
      envir = globalenv()
    )
    on.exit(assign("aaa_run_workflow", original_workflow, envir = globalenv()), add = TRUE)

    failures <- character()
    metadata_path <- file.path(
      QA_ROOT, "QualityAssurance", "fixtures", "valid", "environmental_metadata.csv"
    )

    for (id in ids) {
      plugin <- qa_get_plugin_definition(id)
      cfg <- tryCatch(
        aaa_plugin_test_configuration(id, dataset, context = list(testing = TRUE)),
        error = function(e) e
      )
      if (inherits(cfg, "error")) {
        failures <- c(failures, paste0(id, ": test configuration: ", conditionMessage(cfg)))
        next
      }

      context <- list(
        dataset = dataset,
        parameters = cfg$parameters %||% list(),
        workflow_arguments = utils::modifyList(
          list(
            abundance_type = "counts",
            output_dir = file.path(tempdir(), "Triple_A_plugin_contract", id),
            outputs = names(plugin$outputs %||% list()),
            functional_functions = NULL,
            top_abundance = list(),
            differential_abundance = list(),
            functional_abundance = list(),
            community_structure = list(),
            supervised_multivariate = list(),
            environmental = list(),
            progress_verbosity = "standard",
            verbose = FALSE,
            progress_callback = NULL
          ),
          cfg$workflow_arguments %||% list()
        ),
        environmental_file = metadata_path,
        testing = TRUE
      )

      validation <- tryCatch(aaa_validate_plugin(id, context), error = function(e) e)
      if (inherits(validation, "error")) {
        failures <- c(failures, paste0(id, ": validation error: ", conditionMessage(validation)))
        next
      }
      if (is.data.frame(validation) && "Blocking" %in% names(validation) && any(validation$Blocking %in% TRUE)) {
        failures <- c(
          failures,
          paste0(id, ": blocking validation: ", paste(validation$Message[validation$Blocking %in% TRUE], collapse = "; "))
        )
        next
      }

      result <- tryCatch(
        aaa_run_plugin(id = id, context = context, parameters = context$parameters),
        error = function(e) e
      )
      if (inherits(result, "error")) {
        failures <- c(failures, paste0(id, ": execution error: ", conditionMessage(result)))
        next
      }
      expected_analysis <- plugin$workflow_analysis_id %||% id
      if (!identical(as.character(result$analyses), as.character(expected_analysis))) {
        failures <- c(
          failures,
          paste0(id, ": delegated analysis was '", paste(result$analyses, collapse = ","),
                 "' instead of '", expected_analysis, "'")
        )
      }
    }

    qa_expect_true(length(failures) == 0L, paste(failures, collapse = " | "))
    TRUE
  }
)
