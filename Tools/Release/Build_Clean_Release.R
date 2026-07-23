# Build a minimal Triple_A distribution archive from any project subdirectory.
#
# The archive is an allow-list of what an end user actually needs to run the
# application. Developer-only trees (QualityAssurance, Tools) are never copied
# in, rather than copied and then deleted, so a new development folder cannot
# silently leak into a release.
find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "Run_Triple_A.R")) &&
        dir.exists(file.path(current, "AAApp"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Run this script from inside a Triple_A project.", call. = FALSE)
}

# A beta-testing build keeps the accumulated reference cache (assembly
# accessions and downloaded GFF annotations) so testers don't have to re-run
# every NCBI/GFF lookup a maintainer has already performed. A plain clean
# release resets it to an empty, portable starting point instead. Toggle
# with: TRIPLE_A_KEEP_CACHE=TRUE Rscript Tools/Release/Build_Clean_Release.R
keep_cache <- isTRUE(as.logical(Sys.getenv("TRIPLE_A_KEEP_CACHE", "FALSE")))
release_label <- if (keep_cache) "BetaTesting" else "clean"

project_root <- find_project_root()

# The distribution folder name is fixed rather than taken from basename(project_root):
# a maintainer's working folder may be named anything, and that name must not
# become part of what the user unzips.
dist_name <- "Triple_A"
out_file <- file.path(dirname(project_root), paste0(dist_name, "_", release_label, ".zip"))
stage_parent <- tempfile("triple_a_release_")
stage_root <- file.path(stage_parent, dist_name)
dir.create(stage_root, recursive = TRUE)
on.exit(unlink(stage_parent, recursive = TRUE, force = TRUE), add = TRUE)

# --- Allow-list: everything an end user needs, and nothing else ---------------
# Cache/ is handled separately below because its contents depend on the flavour.
runtime_dirs <- c(
  "AAApp",       # launcher, the five applications and the analysis engine
  "Plugins",     # analysis plugins discovered at runtime
  "Resources"    # contextual documentation shown in-app, plus FunctionalDB
)
runtime_files <- c(
  "Run_Triple_A.R",
  "LICENSE",
  "Triple_A_User_Manual.docx"
)

for (d in runtime_dirs) {
  src <- file.path(project_root, d)
  if (!dir.exists(src)) stop("Required distribution directory is missing: ", d, call. = FALSE)
  if (!file.copy(src, stage_root, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)) {
    stop("Could not stage directory: ", d, call. = FALSE)
  }
}
for (f in runtime_files) {
  src <- file.path(project_root, f)
  if (file.exists(src)) file.copy(src, file.path(stage_root, f), copy.date = TRUE)
}

# Results/ must exist and be empty: the application writes runs into it and the
# directory contract is part of what the engine expects on first launch.
dir.create(file.path(stage_root, "Results", "Logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_root, "Results", "Runs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_root, "Results", "FASTQ"), recursive = TRUE, showWarnings = FALSE)
for (p in c("Results", "Results/Logs", "Results/Runs", "Results/FASTQ")) {
  file.create(file.path(stage_root, p, ".gitkeep"), showWarnings = FALSE)
}

# --- Reference cache ---------------------------------------------------------
cache_src <- file.path(project_root, "Cache")
cache_dst <- file.path(stage_root, "Cache")
dir.create(file.path(cache_dst, "GFF"), recursive = TRUE, showWarnings = FALSE)

if (keep_cache) {
  # Ship the accumulated GenomeCache.sqlite and Cache/GFF as-is so testers
  # start with the reference lookups already resolved. Backups/ is maintainer
  # housekeeping and is deliberately not copied.
  db <- file.path(cache_src, "GenomeCache.sqlite")
  if (file.exists(db)) file.copy(db, file.path(cache_dst, "GenomeCache.sqlite"), copy.date = TRUE)
  gff <- list.files(file.path(cache_src, "GFF"), full.names = TRUE, recursive = TRUE)
  if (length(gff)) {
    file.copy(gff, file.path(cache_dst, "GFF"), copy.date = TRUE)
  }
} else {
  # Reset the portable cache while keeping its required structure.
  file.create(file.path(cache_dst, "GFF", ".gitkeep"), showWarnings = FALSE)
  cache_db <- file.path(cache_dst, "GenomeCache.sqlite")
  if (requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {
    connection <- DBI::dbConnect(RSQLite::SQLite(), cache_db)
    DBI::dbDisconnect(connection)
  } else {
    file.create(cache_db)
  }
}

# --- Transient and development leftovers anywhere in the staged tree ----------
all_files <- list.files(stage_root, recursive = TRUE, all.files = TRUE,
                        full.names = TRUE, no.. = TRUE)
transient <- grepl(
  "\\.(RData|Rhistory|Rproj\\.user|log|tmp|bak)$|(^|/)\\.(claude|vscode|idea|Renviron|DS_Store)(/|$)",
  gsub("\\\\", "/", all_files), ignore.case = TRUE, perl = TRUE
)
if (any(transient)) unlink(all_files[transient], recursive = TRUE, force = TRUE)

# --- Package ------------------------------------------------------------------
if (file.exists(out_file)) unlink(out_file)
old <- setwd(stage_parent)
on.exit(setwd(old), add = TRUE)

if (nzchar(Sys.which("zip"))) {
  status <- system2("zip", c("-qr9", shQuote(out_file), shQuote(dist_name)))
  if (!identical(status, 0L)) stop("The ZIP archive could not be created.", call. = FALSE)
} else if (.Platform$OS.type == "windows" && nzchar(Sys.which("powershell"))) {
  # A bare Windows install has no 'zip' on PATH. Fall back to PowerShell's
  # built-in Compress-Archive via a temporary script, instead of failing
  # outright and instead of building a fragile nested-quoted -Command string.
  escape_ps <- function(x) gsub("'", "''", x, fixed = TRUE)
  ps_script <- tempfile(fileext = ".ps1")
  # -ErrorAction Stop makes Compress-Archive's errors terminating, so a
  # failure inside the cmdlet actually produces a non-zero exit code instead
  # of the script silently reporting success with no archive on disk.
  ps_content <- sprintf(
    "Compress-Archive -Path '%s' -DestinationPath '%s' -Force -ErrorAction Stop",
    escape_ps(dist_name), escape_ps(out_file)
  )
  # Windows PowerShell 5.1 falls back to the system codepage for .ps1 files
  # without a byte-order mark, which corrupts non-ASCII characters. Write an
  # explicit UTF-8 BOM so the script is decoded correctly regardless of locale.
  con <- file(ps_script, open = "wb")
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(charToRaw(enc2utf8(ps_content)), con)
  close(con)
  status <- system2("powershell", c("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", shQuote(ps_script)))
  unlink(ps_script)
  if (!identical(status, 0L)) stop("The ZIP archive could not be created (Compress-Archive fallback failed).", call. = FALSE)
} else {
  stop(
    "No 'zip' utility was found on PATH, and no Windows/PowerShell fallback applies. ",
    "Install zip (e.g. via Git for Windows) and re-run.",
    call. = FALSE
  )
}

# Both branches above can report a zero exit status without actually having
# written the archive (this happened once with Compress-Archive before
# -ErrorAction Stop was added); verify the file is really there before
# declaring success.
if (!file.exists(out_file) || file.info(out_file)$size <= 0) {
  stop("The ZIP archive was not found on disk after packaging: ", out_file, call. = FALSE)
}
cat(
  if (keep_cache) "Beta-testing release (cache retained) created: " else "Clean release created: ",
  out_file, " (", round(file.info(out_file)$size / 1024^2, 1), " MB)\n", sep = ""
)
