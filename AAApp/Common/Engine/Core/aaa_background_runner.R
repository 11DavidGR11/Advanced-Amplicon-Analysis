# =============================================================================
# Background workflow controller
#
# The workflow is executed in a separate R process so the Shiny event loop
# remains responsive and can stop a long-running analysis.
# =============================================================================

aaa_start_background_workflow <- function(
  config,
  project_root,
  verbose = FALSE
) {
  if (!requireNamespace(
    "callr",
    quietly = TRUE
  )) {
    stop(
      "Stopping a running workflow from Shiny requires the R package ",
      "'callr'. Install it with install.packages(\"callr\")."
    )
  }

  session_directory <- tempfile(
    pattern = "triple_a_background_"
  )

  dir.create(
    session_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  paths <- list(
    directory = session_directory,
    config = file.path(
      session_directory,
      "config.rds"
    ),
    progress = file.path(
      session_directory,
      "progress.rds"
    ),
    result = file.path(
      session_directory,
      "result.rds"
    ),
    error = file.path(
      session_directory,
      "error.rds"
    )
  )

  saveRDS(
    config,
    paths$config
  )

  process <- callr::r_bg(
    func = function(project_root,
                    config_path,
                    progress_path,
                    result_path,
                    error_path,
                    verbose) {
      atomic_save <- function(object,
                              path) {
        temporary <- tempfile(
          pattern = paste0(
            basename(path),
            "_"
          ),
          tmpdir = dirname(path)
        )

        saveRDS(
          object,
          temporary
        )

        if (file.exists(path)) {
          unlink(path, force = TRUE)
        }

        if (!file.rename(
          temporary,
          path
        )) {
          file.copy(
            temporary,
            path,
            overwrite = TRUE
          )
          unlink(
            temporary,
            force = TRUE
          )
        }

        invisible(path)
      }

      tryCatch(
        {
          source(
            file.path(
              project_root,
              "AAApp", "Common", "Engine", "Triple_A.R"
            ),
            local = globalenv()
          )

          triple_a_load(
            project_root = project_root,
            install_missing = FALSE,
            verbose = FALSE
          )

          config <- readRDS(
            config_path
          )

          progress_callback <- function(stage,
                                        detail,
                                        completed,
                                        total) {
            atomic_save(
              list(
                stage = stage,
                detail = detail,
                completed = completed,
                total = total,
                timestamp = Sys.time()
              ),
              progress_path
            )
          }

          result <- triple_a_run_config(
            config,
            verbose = verbose,
            progress_callback =
              progress_callback
          )

          atomic_save(
            result,
            result_path
          )

          progress_callback(
            stage = "completed",
            detail = "Workflow completed.",
            completed = 1,
            total = 1
          )
        },
        error = function(error_condition) {
          error_call <- conditionCall(
            error_condition
          )

          atomic_save(
            list(
              message = conditionMessage(
                error_condition
              ),
              call = if (is.null(error_call)) {
                NULL
              } else {
                paste(
                  deparse(error_call),
                  collapse = " "
                )
              },
              timestamp = Sys.time()
            ),
            error_path
          )

          quit(
            save = "no",
            status = 1L,
            runLast = FALSE
          )
        }
      )

      invisible(TRUE)
    },
    args = list(
      project_root = normalizePath(
        project_root,
        winslash = "/",
        mustWork = TRUE
      ),
      config_path = paths$config,
      progress_path = paths$progress,
      result_path = paths$result,
      error_path = paths$error,
      verbose = isTRUE(verbose)
    ),
    supervise = TRUE,
    stdout = "|",
    stderr = "|"
  )

  c(
    list(process = process),
    paths
  )
}

aaa_stop_background_workflow <- function(
  background,
  timeout = 2
) {
  if (is.null(background) ||
    is.null(background$process)) {
    return(invisible(FALSE))
  }

  process <- background$process

  if (!process$is_alive()) {
    return(invisible(TRUE))
  }

  try(
    process$interrupt(),
    silent = TRUE
  )

  deadline <- Sys.time() +
    as.difftime(
      timeout,
      units = "secs"
    )

  while (
    process$is_alive() &&
      Sys.time() < deadline
  ) {
    Sys.sleep(0.05)
  }

  if (process$is_alive()) {
    try(
      process$kill(),
      silent = TRUE
    )
  }

  invisible(!process$is_alive())
}

aaa_cleanup_background_workflow <- function(
  background
) {
  if (is.null(background)) {
    return(invisible(FALSE))
  }

  directory <- background$directory %||%
    NULL

  if (!is.null(directory) &&
    dir.exists(directory)) {
    unlink(
      directory,
      recursive = TRUE,
      force = TRUE
    )
  }

  invisible(TRUE)
}
