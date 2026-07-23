# Triple_A single distribution entry point

required <- c("shiny", "bslib", "callr", "httpuv")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
 stop(
    "Install these launcher packages first: ",
    paste(missing, collapse = ", "),
    ".",
    call. = FALSE
  )
}

get_this_file <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) return(sub("^--file=", "", file_arg[[1L]]))

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    candidate <- frames[[i]]$ofile
    if (!is.null(candidate) && nzchar(candidate)) return(candidate)
  }

  candidate <- file.path(getwd(), "Run_Triple_A.R")
  if (file.exists(candidate)) return(candidate)
  stop("Could not determine the Triple_A project folder.", call. = FALSE)
}

root <- dirname(normalizePath(get_this_file(), winslash = "/", mustWork = TRUE))
launcher_dir <- file.path(root, "AAApp", "Launcher")
if (!file.exists(file.path(launcher_dir, "app.R"))) {
  stop("AAApp/Launcher/app.R was not found inside the Triple_A project.", call. = FALSE)
}

options(triple_a_root = root, triple_a_launcher_dir = launcher_dir)
shiny::runApp(
  appDir = launcher_dir,
  host = "127.0.0.1",
  launch.browser = TRUE,
  display.mode = "normal"
)
