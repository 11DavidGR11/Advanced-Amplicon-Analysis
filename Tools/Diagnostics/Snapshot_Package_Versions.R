# Records the exact version of every R package Triple_A depends on, plus R
# itself and the active Bioconductor release, into
# Tools/Diagnostics/package_versions.json.
#
# Run this after confirming the application and QualityAssurance suite work
# correctly, to freeze a known-good combination. Re-run it whenever
# dependencies are deliberately upgraded. See package_versions.json for how
# to reinstall a pinned version if a future upstream release breaks
# compatibility (this already happened once: newer CVXR releases renamed the
# exported function ANCOMBC::ancombc2() depends on).

get_this_file <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) return(sub("^--file=", "", file_arg[[1L]]))
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    candidate <- frames[[i]]$ofile
    if (!is.null(candidate) && nzchar(candidate)) return(candidate)
  }
  candidate <- file.path(getwd(), "Tools", "Diagnostics", "Snapshot_Package_Versions.R")
  if (file.exists(candidate)) return(candidate)
  stop("Could not determine the Triple_A project folder.", call. = FALSE)
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

root <- find_project_root(dirname(normalizePath(get_this_file(), winslash = "/", mustWork = TRUE)))
source(file.path(root, "AAApp", "Common", "Engine", "Core", "aaa_dependencies.R"))

all_packages <- sort(unique(c(
  TRIPLE_A_DEPENDENCIES$launcher,
  TRIPLE_A_DEPENDENCIES$core,
  unlist(TRIPLE_A_DEPENDENCIES$analyses, use.names = FALSE),
  unlist(TRIPLE_A_DEPENDENCIES$input_formats, use.names = FALSE)
)))
all_packages <- all_packages[nzchar(all_packages)]

installed <- installed.packages()[, c("Package", "Version")]
versions <- setNames(
  vapply(all_packages, function(pkg) {
    if (pkg %in% installed[, "Package"]) installed[installed[, "Package"] == pkg, "Version"][1]
    else NA_character_
  }, character(1)),
  all_packages
)

missing <- names(versions)[is.na(versions)]
if (length(missing) > 0) {
  cat("Warning: not installed, skipped from the snapshot:", paste(missing, collapse = ", "), "\n")
}

snapshot <- list(
  generated = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = as.character(getRversion()),
  bioc_version = if (requireNamespace("BiocManager", quietly = TRUE)) {
    tryCatch(as.character(BiocManager::version()), error = function(e) NA_character_)
  } else NA_character_,
  packages = as.list(versions[!is.na(versions)]),
  known_incompatibilities = list(
    list(
      package = "CVXR",
      issue = paste(
        "CVXR versions that renamed the exported solve() generic to psolve()",
        "(seen starting around CVXR 1.8.x) break ANCOMBC::ancombc2(), which",
        "still calls CVXR::solve() internally."
      ),
      fix = paste(
        "install.packages('CVXR',",
        "repos='https://packagemanager.posit.co/cran/2024-07-01')",
        "to get a pre-rename binary (CVXR 1.0-14 as of this snapshot)."
      )
    )
  )
)

target <- file.path(root, "Tools", "Diagnostics", "package_versions.json")
jsonlite::write_json(snapshot, target, pretty = TRUE, auto_unbox = TRUE, null = "null")
cat("Wrote", length(versions[!is.na(versions)]), "package versions to", target, "\n")
