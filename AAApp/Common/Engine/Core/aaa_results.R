# Result manifest -------------------------------------------------------------
aaa_build_result_manifest <- function(run_dir) {
  files <- list.files(run_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  ext <- tolower(tools::file_ext(files))
  type <- ifelse(ext %in% c("png", "svg", "tif", "tiff", "jpg", "jpeg"), "figure",
    ifelse(ext %in% c("csv", "tsv", "xlsx", "xls"), "table",
      ifelse(ext %in% c("pdf", "docx", "html"), "document", "file")
    )
  )
  rel <- substring(files, nchar(run_dir) + 2)
  data.frame(
    id = vapply(rel, aaa_hash_object, character(1)), title = basename(files), path = rel, type = type, format = ext,
    size_bytes = file.info(files)$size, modified = format(file.info(files)$mtime, tz = "UTC", usetz = TRUE), stringsAsFactors = FALSE
  )
}
aaa_write_result_manifest <- function(run_dir) {
  manifest <- aaa_build_result_manifest(run_dir)
  file <- file.path(run_dir, "manifest.json")
  jsonlite::write_json(list(schema_version = "2.7", run_directory = normalizePath(run_dir, winslash = "/"), results = manifest), file, pretty = TRUE, auto_unbox = TRUE, na = "null")
  file
}
