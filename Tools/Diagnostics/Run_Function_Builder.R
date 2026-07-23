# Triple_A Function Builder developer launcher.
# Launches AAApp/FunctionBuilder directly, bypassing the Launcher. For
# development and manual QA only - end users should use Run_Triple_A.R.
script_dir <- normalizePath(dirname(sys.frame(1)$ofile), winslash = "/", mustWork = TRUE)
app_dir <- normalizePath(file.path(script_dir, "..", "..", "AAApp", "FunctionBuilder"), winslash = "/", mustWork = TRUE)
options(triple_a_function_builder_dir = app_dir)
shiny::runApp(app_dir)
