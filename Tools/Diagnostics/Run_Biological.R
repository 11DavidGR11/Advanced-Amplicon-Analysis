# Triple_A Biological developer launcher.
# Launches AAApp/Biological directly, bypassing the Launcher. For
# development and manual QA only - end users should use Run_Triple_A.R.

required_packages <- c("shiny", "bslib", "DT", "shinyjs")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Triple_A Biological cannot start. Install these packages first: ",
    paste(missing_packages, collapse = ", "),
    ".",
    call. = FALSE
  )
}

get_this_file <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(sub("^--file=", "", file_arg[[1L]]))
  }

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    candidate <- frames[[i]]$ofile
    if (!is.null(candidate) && nzchar(candidate)) return(candidate)
  }

  stop(
    "Could not determine the location of Run_Biological.R. ",
    "Run it with source('Tools/Diagnostics/Run_Biological.R') or Rscript Tools/Diagnostics/Run_Biological.R.",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
read_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

script_dir <- dirname(normalizePath(get_this_file(), winslash = "/", mustWork = TRUE))
app_dir <- normalizePath(file.path(script_dir, "..", "..", "AAApp", "Biological"), winslash = "/", mustWork = TRUE)
port_value <- suppressWarnings(as.integer(read_arg("port", NA_character_)))
if (is.na(port_value)) port_value <- NULL
launch_browser <- !identical(tolower(read_arg("browser", "true")), "false")

options(triple_a_app_dir = app_dir)
shiny::runApp(
  appDir = app_dir,
  host = "127.0.0.1",
  port = port_value,
  launch.browser = launch_browser,
  display.mode = "normal"
)
