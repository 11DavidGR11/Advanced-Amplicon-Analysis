# Central path service for Triple_A.

.AAA_REQUIRED_ROOT_MARKERS <- c(
  file.path("AAApp", "Launcher", "app.R"),
  file.path("AAApp", "Biological", "app.R"),
  file.path("AAApp", "FASTQ", "app.R"),
  file.path("AAApp", "CacheManager", "app.R"),
  file.path("AAApp", "FunctionBuilder", "app.R"),
  file.path("AAApp", "Common", "Engine", "Triple_A.R"),
  "Run_Triple_A.R"
)

aaa_distribution_root <- function(start = getOption("triple_a_root", getwd())) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    valid <- all(file.exists(file.path(current, .AAA_REQUIRED_ROOT_MARKERS)))
    if (valid) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Triple_A distribution root not found from: ", start, call. = FALSE)
    }
    current <- parent
  }
}

aaa_path <- function(..., root = getOption("triple_a_root", NULL), mustWork = FALSE) {
  if (is.null(root)) root <- aaa_distribution_root()
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  path <- file.path(root, ...)
  if (mustWork && !file.exists(path) && !dir.exists(path)) {
    stop("Triple_A path does not exist: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

aaa_cache_path <- function(...) aaa_path("Cache", ...)
aaa_genome_cache_db <- function() aaa_cache_path("GenomeCache.sqlite")
aaa_cached_gff_path <- function(...) aaa_cache_path("GFF", ...)
aaa_cache_backups_path <- function(...) aaa_cache_path("Backups", ...)
