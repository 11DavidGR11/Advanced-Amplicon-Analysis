# Triple_A installation, architecture and source validation.

get_this_file <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) return(sub("^--file=", "", file_arg[[1L]]))
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    candidate <- frames[[i]]$ofile
    if (!is.null(candidate) && nzchar(candidate)) return(candidate)
  }
  candidates <- c(
    file.path(getwd(), "Tools", "Diagnostics", "Verify_Installation.R"),
    file.path(getwd(), "Verify_Installation.R")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) return(existing[[1L]])
  stop("Could not determine the Triple_A project directory.", call. = FALSE)
}

find_project_root <- function(start) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "Run_Triple_A.R")) &&
        dir.exists(file.path(current, "AAApp"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate the Triple_A project root.", call. = FALSE)
}

script_file <- normalizePath(get_this_file(), winslash = "/", mustWork = TRUE)
root <- find_project_root(dirname(script_file))
options(triple_a_root = root)

required_paths <- c(
  "Run_Triple_A.R",
  "AAApp/Launcher/app.R",
  "AAApp/Biological/app.R",
  "AAApp/FASTQ/app.R",
  "AAApp/CacheManager/app.R",
  "AAApp/FunctionBuilder/app.R",
  "AAApp/MultiAmplicon/app.R",
  "AAApp/Common/help.R",
  "AAApp/Common/paths.R",
  "AAApp/Common/Engine/Triple_A.R",
  "Tools/Diagnostics/Verify_Installation.R",
  "Tools/Diagnostics/Run_Biological.R",
  "Tools/Diagnostics/Run_FASTQ.R",
  "Tools/Diagnostics/Run_Cache_Manager.R",
  "Tools/Diagnostics/Run_Function_Builder.R",
  "Tools/Diagnostics/Run_MultiAmplicon.R",
  "Resources/Documentation",
  "Resources/FunctionalDB",
  "Cache/GenomeCache.sqlite",
  "Cache/GFF",
  "Results",
  "Plugins",
  "QualityAssurance/framework/test_framework.R",
  "QualityAssurance/framework/test_helpers.R"
)
missing_paths <- required_paths[!file.exists(file.path(root, required_paths)) & !dir.exists(file.path(root, required_paths))]

cat("Triple_A validation\n")
cat("Project:", root, "\n\n")

if (length(missing_paths)) {
  cat("[ERROR] Missing required paths:\n")
  cat(paste0(" - ", missing_paths, collapse = "\n"), "\n")
} else {
  cat("[OK] Distribution structure is complete.\n")
}

forbidden_paths <- c(
  "AAApp/Run_Triple_A.R", "AAApp/Verify_Installation.R", "Run_MultiAmplicon.R", "Build_Clean_Release.R",
  "AAApp/Biological/Run_Biological.R", "AAApp/FASTQ/Run_FASTQ.R", "AAApp/CacheManager/Run_Cache_Manager.R",
  "AAApp/FunctionBuilder/Run_Function_Builder.R", "AAApp/MultiAmplicon/Run_MultiAmplicon.R",
  "Developer", "Database", "Documentation", "Cache/Proteins"
)
present_forbidden <- forbidden_paths[file.exists(file.path(root, forbidden_paths)) | dir.exists(file.path(root, forbidden_paths))]
if (length(present_forbidden)) {
  cat("[ERROR] Legacy or unsupported paths remain:\n")
  cat(paste0(" - ", present_forbidden, collapse = "\n"), "\n")
} else {
  cat("[OK] No legacy distribution folders were found.\n")
}

r_files <- list.files(root, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
parse_errors <- character()
for (file in r_files) {
  error <- tryCatch({ parse(file = file, keep.source = FALSE); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(error)) parse_errors[file] <- error
}
if (!length(parse_errors)) {
  cat("[OK] All", length(r_files), "R files passed syntax parsing.\n")
} else {
  cat("[ERROR] Syntax errors detected:\n")
  for (file in names(parse_errors)) cat(" -", file, ":", parse_errors[[file]], "\n")
}

source_text <- paste(unlist(lapply(r_files, readLines, warn = FALSE, encoding = "UTF-8")), collapse = "\n")
stale_tokens <- c("Biological_Analysis", "AAApp/Run_Triple_A.R")
stale_found <- stale_tokens[vapply(stale_tokens, grepl, logical(1), x = source_text, fixed = TRUE)]
# The architecture test intentionally mentions the forbidden launcher, so only inspect AAApp runtime files.
runtime_files <- list.files(file.path(root, "AAApp"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
runtime_files <- runtime_files[basename(runtime_files) != "Verify_Installation.R"]
runtime_text <- paste(unlist(lapply(runtime_files, readLines, warn = FALSE, encoding = "UTF-8")), collapse = "\n")
stale_runtime <- c("Biological_Analysis", "Developer/", "Database/")
stale_runtime <- stale_runtime[vapply(stale_runtime, grepl, logical(1), x = runtime_text, fixed = TRUE)]
if (length(stale_runtime)) {
  cat("[ERROR] Obsolete runtime path references:", paste(stale_runtime, collapse = ", "), "\n")
} else {
  cat("[OK] Runtime code contains no obsolete folder references.\n")
}

launcher_packages <- c("shiny", "bslib", "DT", "shinyjs", "callr", "httpuv")
core_packages <- c("dplyr", "tidyr", "tibble", "readr", "readxl", "stringr", "ggplot2", "ggrepel", "openxlsx", "jsonlite", "DBI", "RSQLite")
optional_packages <- c("rentrez", "pheatmap", "vegan", "pls", "dada2")
all_packages <- unique(c(launcher_packages, core_packages, optional_packages))
installed <- vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)

cat("\nPackage status:\n")
for (package in all_packages) {
  group <- if (package %in% launcher_packages) "launcher" else if (package %in% core_packages) "core" else "optional"
  cat(sprintf(" %-12s %-9s %s\n", package, group, if (installed[[package]]) "OK" else "MISSING"))
}
required <- unique(c(launcher_packages, core_packages))
required_missing <- required[!installed[required]]

ok <- !length(missing_paths) && !length(present_forbidden) && !length(parse_errors) && !length(stale_runtime) && !length(required_missing)
if (!ok) {
  cat("\nValidation completed with problems.\n")
  if (length(required_missing)) {
    cat("Install required packages with:\n")
    cat("install.packages(c(", paste(sprintf("'%s'", required_missing), collapse = ", "), "))\n", sep = "")
  }
} else {
  cat("\n[OK] Triple_A core installation is ready.\n")
}
invisible(ok)
