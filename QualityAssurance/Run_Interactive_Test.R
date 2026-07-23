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
checklist <- file.path(script_dir, "manual_checklists", "Complete_Functional_Test.md")
cat(paste(readLines(checklist, warn=FALSE, encoding="UTF-8"), collapse="\n"), "\n")
cat("\nOpen the application in another R session with:\nsource(\"Tools/Diagnostics/Run_Biological.R\")\n")
invisible(checklist)
