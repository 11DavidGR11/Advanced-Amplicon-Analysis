# Triple_A MultiAmplicon developer launcher.
# Launches AAApp/MultiAmplicon directly, bypassing the Launcher. For
# development and manual QA only - end users should use Run_Triple_A.R.
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
app_dir <- normalizePath(file.path(script_dir, "..", "..", "AAApp", "MultiAmplicon"), winslash = "/", mustWork = TRUE)
options(triple_a_multiamplicon_dir = app_dir)
shiny::runApp(app_dir, launch.browser = TRUE)
