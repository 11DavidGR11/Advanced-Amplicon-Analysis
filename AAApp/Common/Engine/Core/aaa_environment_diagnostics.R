# Triple_A environment diagnostics ------------------------------------
aaa_environment_diagnostics <- function(project_dir = getwd(), analyses = character()) {
  core_packages <- c("shiny", "bslib", "DT", "jsonlite")
  analysis_packages <- tryCatch(aaa_required_packages(analyses = analyses),
    error = function(e) character()
  )
  packages <- unique(c(core_packages, analysis_packages))
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  writable <- dir.exists(project_dir) && file.access(project_dir, 2L) == 0L
  excel_ok <- requireNamespace("readxl", quietly = TRUE)
  background_ok <- requireNamespace("callr", quietly = TRUE)
  r_ok <- getRversion() >= "4.2.0"
  rows <- data.frame(
    Check = c(
      "R version", "Required packages", "Writable project directory",
      "Excel import", "Background execution"
    ),
    Status = c(
      if (r_ok) "Valid" else "Warning",
      if (!length(missing)) "Valid" else "Error",
      if (writable) "Valid" else "Error",
      if (excel_ok) "Valid" else "Warning",
      if (background_ok) "Valid" else "Warning"
    ),
    Message = c(
      as.character(getRversion()),
      if (!length(missing)) "Available" else paste("Missing:", paste(missing, collapse = ", ")),
      normalizePath(project_dir, winslash = "/", mustWork = FALSE),
      if (excel_ok) "XLS/XLSX readers available" else "readxl is not installed",
      if (background_ok) "callr is available" else "Execution will use the foreground session"
    ),
    stringsAsFactors = FALSE
  )
  rows
}
