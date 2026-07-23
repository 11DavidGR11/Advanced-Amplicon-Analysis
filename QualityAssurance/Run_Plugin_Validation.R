qa_resolve_script_dir <- function() {
  script_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
    error = function(e) NA_character_
  )
  if (!is.na(script_file) && nzchar(script_file) && file.exists(script_file)) {
    return(dirname(script_file))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  candidates <- unique(c(
    wd,
    file.path(wd, "QualityAssurance"),
    dirname(wd),
    file.path(dirname(wd), "QualityAssurance")
  ))
  valid <- candidates[file.exists(file.path(candidates, "framework", "test_helpers.R"))]
  if (!length(valid)) {
    stop(
      paste0(
        "Cannot locate the QualityAssurance folder from working directory: ", wd,
        ". Run the script from the project root or from inside QualityAssurance."
      ),
      call. = FALSE
    )
  }
  valid[[1]]
}

script_dir <- qa_resolve_script_dir()
source(file.path(script_dir, "framework", "test_helpers.R"))
source(file.path(script_dir, "framework", "test_framework.R"))
QA_ROOT <- qa_find_root(dirname(script_dir))
assign("QA_ROOT", QA_ROOT, envir = globalenv())
options(
  triplea.testing = TRUE,
  triplea.testing.root = file.path(tempdir(), "Triple_A_QualityAssurance")
)
dir.create(getOption("triplea.testing.root"), recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "framework", "plugin_validation.R"))
matrix <- qa_run_plugin_matrix()
print(matrix, row.names = FALSE)
cat("\nReports:\n")
cat(file.path(QA_ROOT, "QualityAssurance", "reports", "plugin_matrix.csv"), "\n")
cat(file.path(QA_ROOT, "QualityAssurance", "reports", "plugin_matrix.txt"), "\n")
if (any(matrix$status == "FAIL")) stop("One or more plugins failed validation. Review plugin_matrix.txt", call. = FALSE)
invisible(matrix)
