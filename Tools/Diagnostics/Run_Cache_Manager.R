# Triple_A Cache Manager developer launcher.
# Launches AAApp/CacheManager directly, bypassing the Launcher. For
# development and manual QA only - end users should use Run_Triple_A.R.
script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
app_dir <- normalizePath(file.path(script_dir, "..", "..", "AAApp", "CacheManager"), winslash = "/", mustWork = TRUE)
options(triple_a_cache_manager_dir = app_dir)
shiny::runApp(app_dir, launch.browser = TRUE)
