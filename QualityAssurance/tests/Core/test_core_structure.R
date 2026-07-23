qa_register_test(
  "CORE_001", "smoke", "critical",
  "Required application structure exists and production code is independent of QualityAssurance",
  function() {
    qa_expect_files(file.path(QA_ROOT, c("AAApp", "AAApp/Biological", "AAApp/FASTQ", "AAApp/Launcher", "AAApp/FunctionBuilder", "AAApp/Common", "AAApp/Common/Engine", "Resources", "Resources/Documentation", "Plugins", "Cache", "Results")))
    # Only runtime production code must be independent of the test suite.
    # Tools/Diagnostics/Verify_Installation.R is intentionally excluded because its job is
    # to validate the QualityAssurance framework and therefore it legitimately
    # references files below QualityAssurance/.
    production_roots <- c(
      file.path(QA_ROOT, "AAApp", "Launcher"),
      file.path(QA_ROOT, "AAApp", "Biological"),
      file.path(QA_ROOT, "AAApp", "FASTQ"),
      file.path(QA_ROOT, "AAApp", "Common"),
      file.path(QA_ROOT, "Plugins")
    )
    files <- unlist(lapply(production_roots, function(path) {
      if (!dir.exists(path)) return(character())
      list.files(path, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
    }), use.names = FALSE)
    text <- paste(vapply(files, function(x) paste(readLines(x, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
    qa_expect_true(!grepl("QualityAssurance", text, fixed = TRUE), "Runtime production source references QualityAssurance")
    TRUE
  }
)

qa_register_test(
  "CORE_002", "smoke", "critical", "All R files parse",
  function() {
    files <- list.files(QA_ROOT, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
    parse_one <- function(file, attempts = 3L) {
      last_error <- NULL
      for (attempt in seq_len(attempts)) {
        result <- tryCatch({ parse(file = file, keep.source = FALSE); TRUE }, error = function(e) e)
        if (isTRUE(result)) return(TRUE)
        last_error <- result
        Sys.sleep(0.15 * attempt)
      }
      stop("Cannot parse ", file, ": ", conditionMessage(last_error), call. = FALSE)
    }
    invisible(vapply(files, parse_one, logical(1)))
    TRUE
  }
)

qa_register_test(
  "CORE_003", "smoke", "high", "Declared method images exist",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Triple_A.R"))
    triple_a_load(QA_ROOT, verbose = FALSE)
    methods <- triple_a_list_methods()
    images <- unique(stats::na.omit(methods$Example_image))
    qa_expect_files(file.path(QA_ROOT, "AAApp", "Biological", "www", "method_examples", images))
    TRUE
  }
)
