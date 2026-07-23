# =============================================================================
# Triple_A AAApp Shiny application
# =============================================================================

library(shiny)
library(bslib)
library(DT)
library(shinyjs)

app_dir <- normalizePath(
  getOption("triple_a_app_dir", getwd()),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(dirname(app_dir), "Common", "paths.R"), local = TRUE)
project_root <- normalizePath(
  getOption("triple_a_root", aaa_distribution_root(app_dir)),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(
  project_root,
  "AAApp", "Common", "Engine", "Triple_A.R"
))

source(file.path(project_root, "AAApp", "Common", "help.R"), local = TRUE)

triple_a_load(
  project_root = project_root,
  install_missing = FALSE,
  verbose = FALSE
)

# Load the stabilized application modules in a fixed order.
module_files <- c(
  "00_shared.R",
  "10_ui.R",
  "20_server.R"
)

for (module_file in module_files) {
  sys.source(
    file.path(app_dir, "modules", module_file),
    envir = environment()
  )
}

shinyApp(ui, server)
