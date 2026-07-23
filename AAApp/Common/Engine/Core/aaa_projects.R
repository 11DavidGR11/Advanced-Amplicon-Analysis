# Triple_A project compatibility layer ---------------------------------------------
aaa_project_file <- function(project_dir) file.path(project_dir, "project.yaml")
aaa_project_runtime_dirs <- function() c("Runs", "Cache")

aaa_create_project <- function(project_dir, name = basename(project_dir), create_runtime = FALSE) {
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(create_runtime)) invisible(lapply(file.path(project_dir, aaa_project_runtime_dirs()), dir.create, recursive = TRUE, showWarnings = FALSE))
  project <- list(
    schema_version = "2.0", name = name,
    created = format(Sys.time(), tz = "UTC", usetz = TRUE),
    modified = format(Sys.time(), tz = "UTC", usetz = TRUE),
    application = "Triple_A",
    storage_policy = "reference-inputs_content-addressed-cache"
  )
  aaa_write_yaml(project, aaa_project_file(project_dir))
  project$directory <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  project
}

aaa_open_project <- function(project_dir) {
  file <- aaa_project_file(project_dir)
  if (!file.exists(file)) stop("Not a compatible Triple_A project: ", project_dir)
  project <- aaa_read_yaml(file)
  if (!identical(as.character(project$schema_version), "2.0")) stop("Unsupported project schema: ", project$schema_version)
  project$directory <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  project
}

aaa_write_yaml <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("yaml", quietly = TRUE)) yaml::write_yaml(x, file) else jsonlite::write_json(x, file, pretty = TRUE, auto_unbox = TRUE)
  invisible(file)
}
aaa_read_yaml <- function(file) if (requireNamespace("yaml", quietly = TRUE)) yaml::read_yaml(file) else jsonlite::read_json(file, simplifyVector = TRUE)
