# Integration tests that drive the real Shiny reactives with shiny::testServer.
# The other UI tests grep the source files, which cannot tell whether the
# experimental design a user configures actually reaches the engine intact.

qa_design_packages <- c("shiny", "bslib", "DT", "shinyjs")

qa_biological_server <- function() {
  app_dir <- file.path(QA_ROOT, "AAApp", "Biological")
  # 10_ui.R builds the layout at source time, so the UI packages have to be
  # attached exactly as the app entry point attaches them.
  for (package in qa_design_packages) {
    suppressWarnings(suppressMessages(
      library(package, character.only = TRUE)
    ))
  }
  env <- new.env(parent = globalenv())
  env$project_root <- QA_ROOT
  # 20_server.R closes over globals defined by the app entry point and by
  # 00_shared.R, so everything is sourced into one environment.
  suppressWarnings(suppressMessages({
    sys.source(file.path(QA_ROOT, "AAApp", "Common", "help.R"), envir = env)
    for (file in c("00_shared.R", "10_ui.R", "20_server.R")) {
      sys.source(file.path(app_dir, "modules", file), envir = env)
    }
  }))
  env$server
}

qa_design_fixture <- function() {
  abundance <- utils::read.csv(
    file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv"),
    check.names = FALSE
  )
  samples <- grep("^Sample_", names(abundance), value = TRUE)
  metadata <- data.frame(
    # Deliberately a different separator: identifier matching normalises it.
    Sample = gsub("_", ".", samples),
    Grupo = c("Ctl", "Trt", "Ctl", "Trt", "Ctl", "Trt",
              "Ctl", "Trt", "Trt", "Trt", "Ctl", "Trt"),
    pH = seq(6, 8.2, length.out = length(samples)),
    stringsAsFactors = FALSE
  )
  roles <- data.frame(
    Column = names(metadata),
    Detected_type = c("character", "character", "numeric"),
    Suggested_role = c("identifier", "experimental_factor", "environmental_variable"),
    Selected_role = c("identifier", "experimental_factor", "environmental_variable"),
    stringsAsFactors = FALSE
  )
  list(abundance = abundance, samples = samples, metadata = metadata, roles = roles)
}

qa_register_test(
  "UI_105", "regression", "critical",
  "The configured experimental design reaches the engine for both design sources",
  function() {
    available <- vapply(
      qa_design_packages,
      function(package) requireNamespace(package, quietly = TRUE),
      logical(1)
    )
    if (!all(available)) {
      # The interface packages are optional for headless engine testing.
      return(TRUE)
    }
    fixture <- qa_design_fixture()
    server <- qa_biological_server()
    results <- new.env(parent = emptyenv())

    shiny::testServer(server, {
      state$data <- fixture$abundance
      state$environmental_data <- fixture$metadata
      state$metadata_roles <- fixture$roles
      session$setInputs(
        taxonomy_layout = "single", taxonomy_column = "Taxonomy",
        sample_selection_mode = "all_numeric", identifier_column_mode = "none",
        abundance_type = "counts", environmental_sample_id = "Sample",
        metadata_role_1 = "identifier", metadata_role_2 = "experimental_factor",
        metadata_role_3 = "environmental_variable"
      )

      # --- Consecutive-block mode must keep behaving exactly as before.
      session$setInputs(
        design_source = "blocks", replicates = "quadruplicate",
        treatments = "A, B, C"
      )
      # input$treatments is debounced; advance the virtual clock so the
      # debounced reactive fires before the design is read.
      session$elapse(1000)
      session$flushReact()
      results$block_design <- sample_design_table()
      results$block_sizes <- design_group_sizes()

      # --- Metadata mode: interleaved and unbalanced.
      session$setInputs(design_source = "metadata", design_group_column = "Grupo")
      session$flushReact()
      results$meta_status <- sample_design_eligibility()
      results$meta_design <- sample_design_table()
      results$meta_sizes <- design_group_sizes()
      results$meta_names <- treatment_names()
      results$meta_da <- differential_abundance_eligibility()$available
      results$meta_plsda <- plsda_eligibility()$available
      results$meta_validation <- validation_report()

      # --- The replicate keyword is hidden in metadata mode and keeps its last
      # value. A stale "none" must not drop the supervised analyses from a
      # design that is in fact replicated.
      session$setInputs(replicates = "none", use_plsda = TRUE, use_splsda = TRUE)
      session$flushReact()
      results$stale_replicates <- selected_analyses()

      # --- A column with one distinct value per sample is not a grouping
      # column. It must be rejected with an explanation, rather than producing
      # single-sample groups and a puzzling lack-of-replicates message later.
      session$setInputs(design_group_column = "pH")
      session$flushReact()
      results$per_sample_status <- sample_design_eligibility()

      # --- An invalid grouping column must block, with guidance, and the block
      # must propagate to the analyses that depend on the design.
      session$setInputs(design_group_column = "Sample")
      session$flushReact()
      results$bad_status <- sample_design_eligibility()
      results$bad_da <- differential_abundance_eligibility()
    })

    qa_expect_true(
      identical(as.integer(results$block_sizes), c(4L, 4L, 4L)),
      "Consecutive-block mode no longer produces three balanced groups."
    )
    qa_expect_true(
      identical(head(results$block_design$Treatment, 4L), rep("A", 4L)),
      "Consecutive-block mode no longer assigns the first block to the first treatment."
    )

    qa_expect_true(
      isTRUE(results$meta_status$available),
      "A valid metadata-driven design was rejected."
    )
    qa_expect_true(
      identical(sort(names(results$meta_sizes)), c("Ctl", "Trt")) &&
        identical(sort(as.integer(results$meta_sizes)), c(5L, 7L)),
      "The metadata design did not produce the declared unbalanced groups."
    )
    qa_expect_true(
      identical(as.character(head(results$meta_design$Treatment, 4L)),
                c("Ctl", "Trt", "Ctl", "Trt")),
      "An interleaved metadata design was flattened into consecutive blocks."
    )
    # Group order is declared in the metadata, so it must not follow the order
    # in which the samples happen to appear in the abundance table.
    qa_expect_true(
      identical(levels(results$meta_design$Treatment), c("Ctl", "Trt")),
      paste0(
        "Group order was not taken from the metadata: ",
        paste(levels(results$meta_design$Treatment), collapse = ", ")
      )
    )
    qa_expect_true(
      identical(as.character(results$meta_names), c("Ctl", "Trt")),
      "Treatment names were not taken from the metadata column."
    )
    qa_expect_true(
      isTRUE(results$meta_da) && isTRUE(results$meta_plsda),
      "An unbalanced design blocked analyses that support unequal group sizes."
    )
    replicate_row <- results$meta_validation[results$meta_validation$Check == "Replicates", , drop = FALSE]
    qa_expect_true(
      nrow(replicate_row) == 1L && identical(replicate_row$Status, "Valid"),
      "Pre-flight validation does not accept an unbalanced design."
    )

    qa_expect_true(
      all(c("plsda", "splsda") %in% results$stale_replicates),
      paste0(
        "A stale replicate keyword removed the supervised analyses from a ",
        "replicated metadata design; selected: ",
        paste(results$stale_replicates, collapse = ", ")
      )
    )

    qa_expect_true(
      !isTRUE(results$per_sample_status$available) &&
        grepl("every sample", results$per_sample_status$reason, fixed = TRUE),
      paste0(
        "A grouping column with one distinct value per sample was accepted; reason given: ",
        results$per_sample_status$reason %||% "<none>"
      )
    )

    qa_expect_true(
      !isTRUE(results$bad_status$available) && length(results$bad_status$guidance) > 0L,
      "Selecting the identifier column as the grouping column was not blocked with guidance."
    )
    qa_expect_true(
      !isTRUE(results$bad_da$available),
      "An unusable experimental design did not block the analyses that depend on it."
    )
    TRUE
  }
)
