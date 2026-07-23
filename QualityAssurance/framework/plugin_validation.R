# Automatic plugin execution validation for Triple_A.
`%||%` <- function(x, y) if (is.null(x)) y else x
# This file belongs only to QualityAssurance; production code must never source it.

qa_plugin_fixture_dataset <- function() {
  abundance_path <- file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv")
  metadata_path <- file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "environmental_metadata.csv")

  abundance <- utils::read.csv(abundance_path, check.names = FALSE, stringsAsFactors = FALSE)
  metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

  sample_columns <- grep("^Sample_", names(abundance), value = TRUE)
  qa_expect_true(length(sample_columns) > 0L, "No synthetic sample columns were detected")

  parsed <- regexec("^Sample_([^0-9]+)([0-9]+)$", sample_columns)
  pieces <- regmatches(sample_columns, parsed)
  valid_names <- lengths(pieces) == 3L
  qa_expect_true(all(valid_names), "Synthetic sample names do not follow Sample_<Treatment><Replicate>")

  sample_design <- data.frame(
    Sample_column = sample_columns,
    Treatment = vapply(pieces, `[[`, character(1), 2L),
    Replicate = as.integer(vapply(pieces, `[[`, character(1), 3L)),
    stringsAsFactors = FALSE
  )

  # Match metadata sample identifier naming flexibly.
  sample_id_candidates <- c("SampleID", "Sample_ID", "Sample", "sample_id")
  sample_id <- sample_id_candidates[sample_id_candidates %in% names(metadata)][1]
  if (is.na(sample_id)) sample_id <- names(metadata)[1]

  detected_types <- vapply(metadata, function(x) {
    if (is.numeric(x)) "numeric" else if (is.logical(x)) "logical" else "categorical"
  }, character(1))

  roles <- data.frame(
    Column = names(metadata),
    Detected_type = unname(detected_types),
    Suggested_role = "ignore",
    Role = "ignore",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  roles$Suggested_role[roles$Column == sample_id] <- "identifier"
  roles$Suggested_role[tolower(roles$Column) %in% c("treatment", "group", "block", "batch")] <- "experimental_factor"
  roles$Suggested_role[roles$Detected_type == "numeric"] <- "environmental_variable"
  roles$Suggested_role[roles$Column == sample_id] <- "identifier"
  roles$Role <- roles$Suggested_role

  qa_expect_true(
    sum(roles$Role == "identifier") == 1L,
    "Synthetic metadata role mapping must contain exactly one identifier"
  )

  if (exists("aaa_new_dataset", mode = "function")) {
    args <- list(
      abundance = abundance,
      sample_design = sample_design,
      metadata = metadata,
      metadata_roles = roles,
      source = list(
        abundance_file = abundance_path,
        metadata_file = metadata_path,
        testing = TRUE
      )
    )
    fml <- names(formals(aaa_new_dataset))
    args <- args[names(args) %in% fml]
    dataset <- do.call(aaa_new_dataset, args)
  } else {
    dataset <- structure(
      list(
        abundance = abundance,
        sample_design = sample_design,
        metadata = metadata,
        metadata_roles = roles,
        source = list(abundance_file = abundance_path, metadata_file = metadata_path),
        schema_version = "2.0"
      ),
      class = c("Triple_A_dataset", "list")
    )
  }

  if (exists("aaa_validate_dataset", mode = "function")) {
    validation <- aaa_validate_dataset(dataset)
    if (is.logical(validation) && length(validation) == 1L) {
      qa_expect_true(validation, "Canonical synthetic dataset did not validate")
    }
  }

  dataset
}

qa_plugin_ids <- function() {
  plugins <- triple_a_list_plugins()
  if (is.data.frame(plugins)) {
    id_col <- intersect(c("ID", "id", "plugin_id", "Plugin"), names(plugins))[1]
    qa_expect_true(!is.na(id_col), "Plugin list has no identifiable ID column")
    return(as.character(plugins[[id_col]]))
  }
  if (is.character(plugins)) return(plugins)
  if (is.list(plugins) && !is.null(names(plugins))) return(names(plugins))
  stop("Unsupported result returned by triple_a_list_plugins()", call. = FALSE)
}

qa_get_plugin_definition <- function(id) {
  getters <- c("aaa_get_plugin", "triple_a_get_plugin", "aaa_plugin_get")
  for (getter in getters) {
    if (exists(getter, mode = "function")) {
      value <- tryCatch(do.call(getter, list(id)), error = function(e) NULL)
      if (!is.null(value)) return(value)
    }
  }
  NULL
}

qa_plugin_defaults <- function(plugin) {
  if (is.null(plugin) || !is.list(plugin)) return(list())
  params <- plugin$parameters %||% plugin$params %||% list()
  if (!is.list(params)) return(list())

  out <- list()
  for (nm in names(params)) {
    item <- params[[nm]]
    if (is.list(item)) {
      value <- item$default %||% item$value %||% item$initial
      if (!is.null(value)) out[[nm]] <- value
    } else if (length(item) == 1L) {
      out[[nm]] <- item
    }
  }
  out
}

qa_validate_one_plugin <- function(id, context) {
  if (!exists("aaa_validate_plugin", mode = "function")) return(TRUE)
  aaa_validate_plugin(id = id, context = context)
}

qa_run_one_plugin <- function(id, context) {
  qa_expect_true(exists("aaa_run_plugin", mode = "function"), "aaa_run_plugin() is unavailable")
  aaa_run_plugin(
    id = id,
    context = context,
    parameters = context$parameters %||% list()
  )
}

qa_result_has_content <- function(result, output_dir) {
  generated <- if (dir.exists(output_dir)) {
    list.files(output_dir, recursive = TRUE, all.files = FALSE, no.. = TRUE)
  } else character()
  !is.null(result) || length(generated) > 0L
}

qa_cache_evidence <- function(result, max_nodes = 250L, max_depth = 6L) {
  fields <- c("cached", "cache_hit", "from_cache", "reused")
  statuses <- c("cache_hit", "hit", "reused")

  status_is_hit <- function(value) {
    if (is.null(value) || length(value) == 0L) return(FALSE)
    any(tolower(as.character(value)) %in% statuses, na.rm = TRUE)
  }

  inspect_one <- function(x) {
    if (is.null(x)) return(FALSE)

    status <- attr(x, "triple_a_cache_status", exact = TRUE)
    if (status_is_hit(status)) return(TRUE)

    # Atomic objects cannot contain named cache fields. This also prevents
    # `$` access warnings/errors on character vectors and other scalars.
    if (!is.list(x)) return(FALSE)

    object_names <- names(x)
    if (is.null(object_names)) object_names <- character()

    for (field in intersect(fields, object_names)) {
      if (isTRUE(x[[field]])) return(TRUE)
    }

    if ("cache_status" %in% object_names && status_is_hit(x[["cache_status"]])) {
      return(TRUE)
    }

    if ("metadata" %in% object_names) {
      metadata <- x[["metadata"]]
      if (is.list(metadata)) {
        metadata_names <- names(metadata)
        if (!is.null(metadata_names) &&
            "cache_status" %in% metadata_names &&
            status_is_hit(metadata[["cache_status"]])) {
          return(TRUE)
        }
      }
    }

    FALSE
  }

  # Use a bounded iterative traversal. Plugin results can contain large data
  # frames, model objects, environments, or self-references; recursively
  # walking every element can freeze an interactive R session.
  queue <- list(list(value = result, depth = 0L))
  visited <- 0L

  while (length(queue) > 0L && visited < max_nodes) {
    current <- queue[[1L]]
    queue <- queue[-1L]
    visited <- visited + 1L

    x <- current$value
    depth <- current$depth

    if (inspect_one(x)) return(TRUE)

    if (depth >= max_depth || !is.list(x)) next

    # Data frames/tibbles can contain thousands of columns/cells and are not
    # expected to hide cache metadata below their top-level attributes/fields.
    if (is.data.frame(x)) next

    children <- unname(x)
    if (length(children) > 0L) {
      remaining <- max_nodes - visited - length(queue)
      if (remaining <= 0L) break
      children <- children[seq_len(min(length(children), remaining))]
      queue <- c(
        queue,
        lapply(children, function(child) list(value = child, depth = depth + 1L))
      )
    }
  }

  FALSE
}

qa_write_plugin_matrix <- function(matrix) {
  report_dir <- file.path(QA_ROOT, "QualityAssurance", "reports")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(matrix, file.path(report_dir, "plugin_matrix.csv"), row.names = FALSE, na = "")

  text <- c(
    "Triple_A - Automatic Plugin Validation",
    paste0("Generated: ", Sys.time()),
    paste0("Plugins tested: ", nrow(matrix)),
    paste0("Passed: ", sum(matrix$status == "PASS")),
    paste0("Failed: ", sum(matrix$status == "FAIL")),
    "",
    apply(matrix, 1, function(x) paste(sprintf("%s=%s", names(x), x), collapse = " | "))
  )
  writeLines(text, file.path(report_dir, "plugin_matrix.txt"), useBytes = TRUE)

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(matrix, file.path(report_dir, "plugin_matrix.json"), pretty = TRUE, na = "null")
  }
  invisible(matrix)
}

qa_run_plugin_matrix <- function() {
  source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Triple_A.R"))
  triple_a_load(QA_ROOT, verbose = FALSE)

  dataset <- qa_plugin_fixture_dataset()
  ids <- qa_plugin_ids()
  qa_expect_true(length(ids) > 0L, "No plugins were discovered")

  rows <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    id <- ids[[i]]
    plugin <- qa_get_plugin_definition(id)
    parameters <- qa_plugin_defaults(plugin)
    plugin_test_cfg <- if (exists("aaa_plugin_test_configuration", mode = "function")) {
      aaa_plugin_test_configuration(id, dataset, context = list(testing = TRUE))
    } else {
      list(parameters = parameters, workflow_arguments = list())
    }
    parameters <- utils::modifyList(parameters, plugin_test_cfg$parameters %||% list())
    testing_root <- getOption("triplea.testing.root")
    if (is.null(testing_root) || length(testing_root) != 1L ||
        is.na(testing_root) || !nzchar(as.character(testing_root))) {
      testing_root <- file.path(QA_ROOT, "QualityAssurance", "tmp")
    }
    testing_root <- normalizePath(
      as.character(testing_root),
      winslash = "/",
      mustWork = FALSE
    )
    plugin_root <- file.path(testing_root, "plugin_matrix", id)
    unlink(plugin_root, recursive = TRUE, force = TRUE)
    dir.create(plugin_root, recursive = TRUE, showWarnings = FALSE)
    qa_expect_true(
      dir.exists(plugin_root),
      paste0("Could not create plugin test directory: ", plugin_root)
    )

    selected_outputs <- if (!is.null(plugin) && is.list(plugin$outputs)) {
      names(plugin$outputs)
    } else {
      character()
    }

    metadata_path <- file.path(
      QA_ROOT, "QualityAssurance", "fixtures", "valid", "environmental_metadata.csv"
    )

    workflow_arguments <- utils::modifyList(
      list(
        abundance_type = "counts",
        output_dir = plugin_root,
        outputs = selected_outputs,
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
      plugin_test_cfg$workflow_arguments %||% list()
    )

    context <- list(
      dataset = dataset,
      parameters = parameters,
      params = parameters,
      config = parameters,
      workflow_arguments = workflow_arguments,
      output_dir = plugin_root,
      project_dir = plugin_root,
      analysis_name = paste0("QualityAssurance_", id),
      environmental_file = metadata_path,
      testing = TRUE,
      use_cache = TRUE,
      verbose = FALSE,
      progress_verbosity = "standard"
    )

    validation_status <- "PASS"
    run_status <- "PASS"
    result_status <- "PASS"
    cache_status <- "UNVERIFIED"
    details <- ""
    elapsed_1 <- NA_real_
    elapsed_2 <- NA_real_

    plugin_warnings <- character()
    first <- tryCatch({
      validation_result <- qa_validate_one_plugin(id, context)
      if (is.data.frame(validation_result) &&
          "Blocking" %in% names(validation_result) &&
          any(isTRUE(validation_result$Blocking) | validation_result$Blocking %in% TRUE)) {
        validation_status <<- "FAIL"
        stop(
          paste(validation_result$Message[validation_result$Blocking], collapse = "; "),
          call. = FALSE
        )
      }
      started <- proc.time()[["elapsed"]]
      value <- withCallingHandlers(
        qa_run_one_plugin(id, context),
        warning = function(w) {
          plugin_warnings <<- unique(c(plugin_warnings, conditionMessage(w)))
          invokeRestart("muffleWarning")
        }
      )
      elapsed_1 <- proc.time()[["elapsed"]] - started
      value
    }, error = function(e) {
      run_status <<- "FAIL"
      details <<- conditionMessage(e)
      NULL
    })

    if (identical(run_status, "PASS") && !qa_result_has_content(first, plugin_root)) {
      result_status <- "FAIL"
      details <- "Plugin returned NULL and generated no result files"
    }

    if (identical(run_status, "PASS") && identical(result_status, "PASS")) {
      second <- tryCatch({
        started <- proc.time()[["elapsed"]]
        value <- withCallingHandlers(
          qa_run_one_plugin(id, context),
          warning = function(w) {
            plugin_warnings <<- unique(c(plugin_warnings, conditionMessage(w)))
            invokeRestart("muffleWarning")
          }
        )
        elapsed_2 <- proc.time()[["elapsed"]] - started
        value
      }, error = function(e) {
        details <<- paste(c(details, paste0("Second execution: ", conditionMessage(e))), collapse = " | ")
        NULL
      })
      if (qa_cache_evidence(second)) cache_status <- "PASS"
    }

    if (length(plugin_warnings) > 0L) {
      warning_text <- paste(plugin_warnings, collapse = " | ")
      details <- paste(c(details[nzchar(details)], paste0("Non-blocking warning: ", warning_text)), collapse = " | ")
    }

    overall <- if ("FAIL" %in% c(validation_status, run_status, result_status)) "FAIL" else "PASS"
    rows[[i]] <- data.frame(
      plugin = id,
      validation = validation_status,
      execution = run_status,
      results = result_status,
      cache = cache_status,
      first_run_seconds = round(elapsed_1, 3),
      second_run_seconds = round(elapsed_2, 3),
      status = overall,
      details = details,
      stringsAsFactors = FALSE
    )
  }

  matrix <- do.call(rbind, rows)
  qa_write_plugin_matrix(matrix)
  matrix
}
