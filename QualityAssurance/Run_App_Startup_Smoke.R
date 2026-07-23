# Starts each Triple_A Shiny component in a child R process, verifies that its
# local HTTP port opens, then terminates the process. Run from the project root.
source(file.path("QualityAssurance", "framework", "test_helpers.R"))
root <- qa_find_root()

required <- c("callr", "httpuv", "shiny")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing smoke-test packages: ", paste(missing, collapse = ", "))

port_open <- function(port, timeout = 30) {
  deadline <- Sys.time() + timeout
  repeat {
    con <- suppressWarnings(tryCatch(socketConnection("127.0.0.1", port, open = "r+", timeout = 1), error = function(e) NULL))
    if (!is.null(con)) { close(con); return(TRUE) }
    if (Sys.time() >= deadline) return(FALSE)
    Sys.sleep(0.25)
  }
}

start_app <- function(name, app_dir, option_name, required_packages = character()) {
  absent <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    cat(sprintf("[SKIP] %s: missing %s\n", name, paste(absent, collapse = ", ")))
    return(invisible(NA))
  }
  port <- as.integer(httpuv::randomPort())
  log <- tempfile(pattern = paste0("triple_a_", gsub("[^a-z]+", "_", tolower(name))), fileext = ".log")
  process <- callr::r_bg(function(app_dir, option_name, port) {
    x <- list(); x[[option_name]] <- app_dir; do.call(options, x)
    shiny::runApp(app_dir, host = "127.0.0.1", port = port, launch.browser = FALSE)
  }, args = list(normalizePath(app_dir, winslash = "/"), option_name, port), stdout = log, stderr = log, supervise = TRUE)
  on.exit(if (process$is_alive()) process$kill(), add = TRUE)
  ok <- port_open(port, 30)
  if (!ok) {
    details <- if (file.exists(log)) paste(readLines(log, warn = FALSE), collapse = "\n") else "No log"
    stop(name, " did not open its port.\n", details, call. = FALSE)
  }
  process$kill()
  cat(sprintf("[PASS] %s opened port %d\n", name, port))
  invisible(TRUE)
}

start_app("Biological Analysis", file.path(root, "AAApp", "Biological"), "triple_a_app_dir", c("bslib", "DT", "shinyjs"))
start_app("FASTQ", file.path(root, "AAApp", "FASTQ"), "triple_a_fastq_app_dir", c("bslib", "dada2", "jsonlite"))
start_app("MultiAmplicon", file.path(root, "AAApp", "MultiAmplicon"), "triple_a_multiamplicon_dir", c("bslib", "readxl"))
start_app("Cache Manager", file.path(root, "AAApp", "CacheManager"), "triple_a_cache_manager_dir", c("bslib", "DBI", "RSQLite"))
start_app("Function Builder", file.path(root, "AAApp", "FunctionBuilder"), "triple_a_function_builder_dir", c("bslib", "jsonlite"))
cat("Startup smoke test completed.\n")
