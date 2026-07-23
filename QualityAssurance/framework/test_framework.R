# Minimal self-contained test framework for Triple_A.
# This file is always sourced before test discovery by Run_All_Tests.R.

qa_registry <- new.env(parent = emptyenv())
qa_registry$tests <- list()

qa_register_test <- function(id,
                             level = c("smoke", "functional", "regression", "release"),
                             severity = c("critical", "high", "medium", "low"),
                             description,
                             run) {
  level <- match.arg(level)
  severity <- match.arg(severity)
  if (!is.character(id) || length(id) != 1L || !nzchar(id)) {
    stop("Test id must be a non-empty string.", call. = FALSE)
  }
  if (!is.function(run)) stop("Test run argument must be a function.", call. = FALSE)
  qa_registry$tests[[id]] <- list(
    id = id,
    level = level,
    severity = severity,
    description = as.character(description),
    run = run
  )
  invisible(id)
}

# Concise registration helper retained for all test files.
# `type` is accepted as a compatibility alias for `level` because some
# historical and third-party tests used qa_test(type = "regression", ...).
qa_test <- function(id,
                    level = NULL,
                    severity,
                    description,
                    code,
                    type = NULL,
                    ...) {
  dots <- list(...)
  if (length(dots)) {
    stop("Unused qa_test arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
  }

  if (is.null(level) && !is.null(type)) level <- type
  if (!is.null(level) && !is.null(type) && !identical(as.character(level), as.character(type))) {
    stop("qa_test received conflicting values for 'level' and compatibility alias 'type'.", call. = FALSE)
  }
  if (is.null(level)) level <- "regression"

  expr <- substitute(code)
  env <- parent.frame()

  # Support both historical expression-style tests:
  #   qa_test(..., { ... })
  # and explicit function-style tests:
  #   qa_test(..., function() { ... })
  # Wrapping a function expression in another function causes the test to
  # return the function object itself (reported as "function ()") instead of
  # executing its body.
  run <- if (is.call(expr) && identical(expr[[1L]], as.name("function"))) {
    eval(expr, envir = env)
  } else {
    as.function(c(alist(), expr), envir = env)
  }

  qa_register_test(
    id = id,
    level = level,
    severity = severity,
    description = description,
    run = run
  )
}

qa_discover_tests <- function(dir) {
  if (!dir.exists(dir)) stop("Test directory not found: ", dir, call. = FALSE)
  files <- sort(list.files(dir, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE))
  for (file in files) {
    sys.source(file, envir = globalenv(), keep.source = FALSE)
  }
  invisible(files)
}

qa_run_registered <- function(levels = c("smoke", "functional", "regression", "release")) {
  rows <- list()
  tests <- qa_registry$tests
  for (id in names(tests)) {
    test <- tests[[id]]
    if (!test$level %in% levels) next
    started <- Sys.time()
    status <- "PASS"
    details <- ""
    tryCatch(
      {
        value <- test$run()
        if (!isTRUE(value)) {
          status <- "FAIL"
          details <- paste(capture.output(str(value)), collapse = " ")
        }
      },
      warning = function(w) {
        status <<- "WARN"
        details <<- conditionMessage(w)
      },
      error = function(e) {
        status <<- "FAIL"
        details <<- conditionMessage(e)
      }
    )
    rows[[length(rows) + 1L]] <- data.frame(
      id = id,
      level = test$level,
      severity = test$severity,
      description = test$description,
      status = status,
      duration_seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 3),
      details = details,
      stringsAsFactors = FALSE
    )
    cat(sprintf("[%s] %s | %s | %s | %s\n", status, id, test$level, test$severity, test$description))
    if (identical(status, "FAIL") && identical(test$severity, "critical")) break
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}
