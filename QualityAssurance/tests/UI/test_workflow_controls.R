qa_ui_103_sources <- function() {
  ui_path <- file.path(QA_ROOT, "AAApp", "Biological", "modules", "10_ui.R")
  server_path <- file.path(QA_ROOT, "AAApp", "Biological", "modules", "20_server.R")
  qa_expect_files(c(ui_path, server_path))
  list(
    ui = paste(readLines(ui_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
    server = paste(readLines(server_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  )
}

qa_register_test(
  "UI_103A", "regression", "critical",
  "Unavailable analyses use real HTML locks for selectors and parameter panels",
  function() {
    src <- qa_ui_103_sources()

    for (id in c(
      "top_abundance_controls",
      "diversity_controls",
      "community_structure_controls",
      "community_comparison_controls",
      "functional_potential_controls",
      "functional_abundance_controls",
      "differential_abundance_controls",
      "plsda_controls",
      "splsda_controls",
      "rda_controls"
    )) {
      qa_expect_true(
        grepl(
          paste0('tags\\$fieldset\\([\\s\\S]*?id[[:space:]]*=[[:space:]]*["\']', id, '["\']'),
          src$ui,
          perl = TRUE
        ),
        paste0("Analysis controls are not wrapped in a lockable fieldset: ", id)
      )
    }

    qa_expect_true(
      grepl('root\\.disabled[[:space:]]*=[[:space:]]*%s', src$server),
      "The analysis fieldset disabled property is not updated"
    )
    qa_expect_true(
      grepl('shinyjs::toggleState\\([[:space:]]*input_id', src$server),
      "The main analysis selector is not enabled and disabled centrally"
    )
    qa_expect_true(
      grepl('updateCheckboxInput\\([\\s\\S]*?input_id[\\s\\S]*?value[[:space:]]*=[[:space:]]*FALSE', src$server, perl = TRUE),
      "Unavailable analyses are not automatically deselected"
    )
    qa_expect_true(
      grepl('session\\$onFlushed\\(', src$server),
      "Analysis locks are not reapplied after dynamic Shiny UI rendering"
    )
    qa_expect_true(
      grepl('aaa_check_plsda_eligibility\\(', src$server),
      "PLS-DA availability does not use the centralized engine eligibility rule"
    )

    # Every analysis module must be gated by an engine-level eligibility rule,
    # not only the four that originally had one.
    for (check in c(
      "aaa_check_top_abundance_eligibility",
      "aaa_check_community_structure_eligibility",
      "aaa_check_differential_abundance_eligibility",
      "aaa_check_functional_potential_eligibility",
      "aaa_check_taxon_association_eligibility",
      "aaa_check_constrained_ordination_eligibility",
      "aaa_check_variance_partitioning_eligibility",
      "aaa_check_rda_eligibility"
    )) {
      qa_expect_true(
        grepl(paste0(check, "\\("), src$server),
        paste0("Analysis availability does not use the engine eligibility rule: ", check)
      )
    }

    # A locked analysis has to explain how to unlock itself.
    qa_expect_true(
      grepl('How to unlock', src$server, fixed = TRUE),
      "Locked analyses do not render unlocking instructions"
    )
    qa_expect_true(
      grepl('status\\$guidance', src$server),
      "The availability notice does not render the eligibility guidance steps"
    )

    TRUE
  }
)

qa_register_test(
  "UI_103B", "regression", "critical",
  "Pairwise differential comparisons are selectable",
  function() {
    src <- qa_ui_103_sources()
    qa_expect_true(
      grepl('checkboxGroupInput\\([[:space:]]*["\']differential_comparisons["\']', src$server),
      "Selectable pairwise differential comparisons are not wired"
    )
    qa_expect_true(
      grepl('comparisons[[:space:]]*=[[:space:]]*input\\$differential_comparisons', src$server),
      "Selected comparisons are not passed to the differential-analysis configuration"
    )
    TRUE
  }
)

qa_register_test(
  "UI_103C", "regression", "critical",
  "Analysis history supports opening and reloading saved projects",
  function() {
    src <- qa_ui_103_sources()
    qa_expect_true(
      grepl('actionButton\\(["\']reload_history_run["\']', src$ui),
      "The history reload button is missing"
    )
    qa_expect_true(
      grepl('observeEvent\\([[:space:]]*input\\$reload_history_run', src$server),
      "History project reload action is not wired"
    )
    qa_expect_true(
      grepl('Run_snapshot\\.rds', src$server),
      "History reload does not use the saved run snapshot"
    )
    TRUE
  }
)

qa_register_test(
  "UI_103D", "regression", "critical",
  "Cache status uses the centralized project cache",
  function() {
    src <- qa_ui_103_sources()
    qa_expect_true(
      grepl('aaa_cache_status\\([[:space:]]*project_root[[:space:]]*\\)', src$server),
      "The cache panel does not use the centralized project cache"
    )
    TRUE
  }
)

qa_test(
  id = "UI_103E",
  level = "regression",
  severity = "high",
  description = "All analyses start unselected and require explicit user activation",
  {
    ui_file <- file.path(
      QA_PROJECT_ROOT,
      "AAApp", "Biological", "modules", "10_ui.R"
    )
    qa_expect_true(file.exists(ui_file), "Biological analysis UI file not found.")
    ui_src <- paste(readLines(ui_file, warn = FALSE), collapse = "\n")

    analysis_inputs <- c(
      "use_functional_potential",
      "use_functional_abundance",
      "use_top_abundance",
      "use_differential_abundance",
      "use_community_structure",
      "use_plsda",
      "use_rda"
    )

    # analysis_option() is a shared UI helper that always renders
    # checkboxInput(id, label, FALSE); selectors built through it never
    # appear as a literal checkboxInput(..., FALSE) call in the source.
    for (input_id in analysis_inputs) {
      direct_pattern <- paste0(
        "checkboxInput\\([[:space:]]*[\"']", input_id,
        "[\"'][^\\)]*,[[:space:]]*FALSE[[:space:]]*\\)"
      )
      helper_pattern <- paste0(
        "analysis_option\\([[:space:]]*[\"']", input_id, "[\"']"
      )
      qa_expect_true(
        grepl(direct_pattern, ui_src, perl = TRUE) ||
          grepl(helper_pattern, ui_src, perl = TRUE),
        paste0("Analysis selector '", input_id, "' is not unselected at startup.")
      )
    }

    qa_expect_true(
      grepl(
        "analysis_option[[:space:]]*<-[[:space:]]*function[^\\n]*\\{[\\s\\S]*?checkboxInput\\([^\\)]*,[[:space:]]*FALSE[[:space:]]*\\)",
        ui_src,
        perl = TRUE
      ),
      "The analysis_option() helper no longer defaults its checkbox to FALSE."
    )

    qa_expect_true(
      grepl(
        'checkboxGroupInput\\([[:space:]]*[\"\']analyses[\"\'][^\\)]*selected[[:space:]]*=[[:space:]]*character\\(0\\)',
        ui_src,
        perl = TRUE
      ),
      "The hidden legacy analysis selector does not start empty."
    )

    TRUE
  }
)
