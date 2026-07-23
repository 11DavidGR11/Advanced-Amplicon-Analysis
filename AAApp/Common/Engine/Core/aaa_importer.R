# =============================================================================
# Triple_A canonical table importer
#
# Public functions:
#   aaa_supported_input_formats()
#   aaa_detect_table_format()
#   aaa_import_table()
#   aaa_list_workbook_sheets()
#
# The importer detects the real file format from its contents. Extensions are
# treated only as hints. This preserves compatibility with legacy pipelines
# that write tab-delimited text using an .xls suffix.
# =============================================================================

aaa_supported_input_formats <- function() {
  c("csv", "tsv", "txt", "xls", "xlsx")
}

aaa_normalize_extension <- function(path = NULL, original_name = NULL) {
  candidate <- original_name
  if (is.null(candidate) || length(candidate) != 1L || is.na(candidate) || !nzchar(candidate)) {
    candidate <- path
  }
  if (is.null(candidate) || length(candidate) != 1L || is.na(candidate)) candidate <- ""
  tolower(tools::file_ext(candidate))
}

aaa_read_signature <- function(path, n = 16L) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "raw", n = n)
}

aaa_raw_starts_with <- function(value, signature) {
  length(value) >= length(signature) && identical(value[seq_along(signature)], signature)
}

aaa_file_looks_textual <- function(path, sample_bytes = 65536L) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = sample_bytes)
  if (length(bytes) == 0L) {
    return(FALSE)
  }

  numeric_bytes <- as.integer(bytes)
  # NUL bytes are a strong binary-file signal. Ordinary control characters
  # are allowed only for tab, carriage return, line feed and form feed.
  if (any(numeric_bytes == 0L)) {
    return(FALSE)
  }
  allowed_controls <- c(9L, 10L, 12L, 13L)
  invalid_controls <- numeric_bytes < 32L & !numeric_bytes %in% allowed_controls
  mean(invalid_controls) < 0.01
}

aaa_detect_zip_table_format <- function(path) {
  members <- tryCatch(utils::unzip(path, list = TRUE)$Name, error = function(e) character())
  if (length(members) == 0L) {
    return("zip_unknown")
  }
  lower <- tolower(members)
  if (any(lower == "xl/workbook.xml") || any(grepl("^xl/worksheets/", lower))) {
    return("xlsx")
  }
  if (any(lower == "mimetype") && any(grepl("content\\.xml$", lower))) {
    return("ods")
  }
  "zip_unknown"
}

aaa_detect_table_format <- function(path, original_name = NULL) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop("A readable input-file path is required for format detection.", call. = FALSE)
  }

  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0L) {
    return("empty")
  }

  signature <- aaa_read_signature(path, 16L)
  ole_signature <- as.raw(c(0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1))
  zip_signature <- as.raw(c(0x50, 0x4B, 0x03, 0x04))
  gzip_signature <- as.raw(c(0x1F, 0x8B))

  if (aaa_raw_starts_with(signature, ole_signature)) {
    return("xls")
  }
  if (aaa_raw_starts_with(signature, zip_signature)) {
    return(aaa_detect_zip_table_format(path))
  }
  if (aaa_raw_starts_with(signature, gzip_signature)) {
    return("gzip")
  }
  if (aaa_file_looks_textual(path)) {
    return("text")
  }

  declared <- aaa_normalize_extension(path, original_name)
  if (declared %in% c("xls", "xlsx")) {
    return(paste0("unrecognized_", declared))
  }
  "unknown"
}

aaa_stage_excel_file <- function(path, extension = NULL) {
  if (!file.exists(path)) stop("Input file does not exist: ", path, call. = FALSE)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) stop("The uploaded Excel file is empty or cannot be read.", call. = FALSE)

  if (is.null(extension) || length(extension) != 1L || is.na(extension) || !nzchar(extension)) extension <- tools::file_ext(path)
  extension <- tolower(extension)
  if (!extension %in% c("xls", "xlsx")) stop("Excel staging requires an .xls or .xlsx format.", call. = FALSE)

  staged_dir <- file.path(tempdir(), "Triple_A_excel_uploads")
  if (!dir.exists(staged_dir)) dir.create(staged_dir, recursive = TRUE, showWarnings = FALSE)

  token <- paste0(format(Sys.time(), "%Y%m%d%H%M%OS6"), "_", Sys.getpid(), "_", sample.int(.Machine$integer.max, 1L))
  token <- gsub("[^0-9A-Za-z_]", "", token)
  staged_path <- file.path(staged_dir, paste0("triple_a_", token, ".", extension))

  copied <- file.copy(path, staged_path, overwrite = TRUE, copy.mode = FALSE)
  if (!isTRUE(copied) || !file.exists(staged_path)) stop("Triple_A could not stage the uploaded Excel file for reading.", call. = FALSE)
  staged_info <- file.info(staged_path)
  if (is.na(staged_info$size) || staged_info$size != info$size) {
    unlink(staged_path, force = TRUE)
    stop("The uploaded Excel file could not be copied completely.", call. = FALSE)
  }
  normalizePath(staged_path, winslash = "/", mustWork = TRUE)
}

aaa_read_excel_table <- function(path, extension, sheet = NULL) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read genuine .xls and .xlsx workbooks. Install it with install.packages('readxl').", call. = FALSE)
  }
  staged_path <- aaa_stage_excel_file(path, extension)
  on.exit(unlink(staged_path, force = TRUE), add = TRUE)

  args <- list(path = staged_path, .name_repair = "minimal")
  if (!is.null(sheet) && length(sheet) == 1L && !is.na(sheet) && nzchar(as.character(sheet))) args$sheet <- sheet
  reader <- if (identical(extension, "xls")) readxl::read_xls else readxl::read_xlsx

  tryCatch(do.call(reader, args), error = function(e) {
    stop("Triple_A could not open the genuine Excel workbook. It may be damaged or password-protected. Original error: ", conditionMessage(e), call. = FALSE)
  })
}

aaa_decode_text_lines <- function(path, n = 20L) {
  encodings <- c("UTF-8-BOM", "UTF-8", "latin1")
  for (encoding in encodings) {
    value <- tryCatch(readLines(path, n = n, warn = FALSE, encoding = encoding), error = function(e) NULL)
    if (!is.null(value) && length(value) > 0L) {
      return(value)
    }
  }
  character()
}

aaa_detect_text_delimiter <- function(path, fallback = "\t") {
  lines <- aaa_decode_text_lines(path, n = 20L)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0L) {
    return(fallback)
  }

  candidates <- c(tab = "\t", comma = ",", semicolon = ";", pipe = "|")
  scores <- vapply(candidates, function(delim) {
    fields <- lengths(strsplit(lines, delim, fixed = TRUE))
    useful <- fields > 1L
    if (!any(useful)) {
      return(-Inf)
    }
    # Prefer delimiters yielding many columns consistently across lines.
    stats::median(fields[useful]) - stats::mad(fields[useful], constant = 1) - mean(!useful)
  }, numeric(1))

  if (!any(is.finite(scores))) {
    return(fallback)
  }
  unname(candidates[[which.max(scores)]])
}

aaa_import_text_table <- function(path, delimiter = "auto") {
  if (!requireNamespace("readr", quietly = TRUE)) stop("Package 'readr' is required to read delimited text tables.", call. = FALSE)
  resolved <- if (identical(delimiter, "auto")) aaa_detect_text_delimiter(path) else delimiter
  imported <- suppressWarnings(readr::read_delim(
    path,
    delim = resolved,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    name_repair = "minimal",
    locale = readr::locale(encoding = "UTF-8"),
    trim_ws = TRUE,
    progress = FALSE
  ))

  imported <- as.data.frame(imported, check.names = FALSE, stringsAsFactors = FALSE)
  imported[] <- lapply(imported, function(column) {
    utils::type.convert(
      column,
      as.is = TRUE,
      na.strings = c("", "NA", "NaN", "null", "NULL")
    )
  })
  imported
}

aaa_format_mismatch_warning <- function(declared, detected) {
  if (!nzchar(declared)) {
    return(invisible(FALSE))
  }
  expected <- switch(declared,
    csv = "text",
    tsv = "text",
    txt = "text",
    xls = "xls",
    xlsx = "xlsx",
    declared
  )
  if (identical(expected, detected)) {
    return(invisible(FALSE))
  }
  if (declared %in% c("xls", "xlsx") && identical(detected, "text")) {
    warning(
      "The file is named .", declared, " but contains delimited text. ",
      "Triple_A imported it as a text table for compatibility with legacy files.",
      call. = FALSE, immediate. = TRUE
    )
    return(invisible(TRUE))
  }
  warning("The filename extension suggests .", declared, " but the content was detected as ", detected, ". Triple_A used the detected content format.", call. = FALSE, immediate. = TRUE)
  invisible(TRUE)
}

aaa_list_workbook_sheets <- function(path, original_name = NULL) {
  if (!file.exists(path)) stop("Input file does not exist: ", path, call. = FALSE)
  detected <- aaa_detect_table_format(path, original_name)
  if (!detected %in% c("xls", "xlsx")) {
    return(character())
  }
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required to inspect Excel workbooks.", call. = FALSE)

  staged_path <- aaa_stage_excel_file(path, detected)
  on.exit(unlink(staged_path, force = TRUE), add = TRUE)
  tryCatch(readxl::excel_sheets(staged_path), error = function(e) {
    stop("Triple_A could not inspect the Excel workbook. Original error: ", conditionMessage(e), call. = FALSE)
  })
}

aaa_import_table <- function(
  path,
  sheet = NULL,
  delimiter = "auto",
  check_minimum_columns = TRUE,
  original_name = NULL
) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("A single input-file path is required.", call. = FALSE)
  if (!file.exists(path)) stop("Input file does not exist: ", path, call. = FALSE)

  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0L) stop("The input table is empty or unreadable.", call. = FALSE)

  declared <- aaa_normalize_extension(path, original_name)
  detected <- aaa_detect_table_format(path, original_name)
  aaa_format_mismatch_warning(declared, detected)

  imported <- switch(detected,
    text = aaa_import_text_table(path, delimiter),
    xls = aaa_read_excel_table(path, "xls", sheet),
    xlsx = aaa_read_excel_table(path, "xlsx", sheet),
    ods = stop("ODS content was detected. Save the table as XLSX, XLS, CSV, TSV or TXT before importing.", call. = FALSE),
    gzip = stop("A gzip-compressed file was detected. Decompress it before importing.", call. = FALSE),
    empty = stop("The input table is empty.", call. = FALSE),
    zip_unknown = stop("A ZIP-based file was detected, but it is not a recognized XLSX workbook.", call. = FALSE),
    unrecognized_xls = stop("The file is neither a genuine XLS workbook nor readable delimited text.", call. = FALSE),
    unrecognized_xlsx = stop("The file is neither a genuine XLSX workbook nor readable delimited text.", call. = FALSE),
    stop("Triple_A could not recognize the table format from its contents. Supported table contents are delimited text, XLS and XLSX.", call. = FALSE)
  )

  imported <- as.data.frame(imported, check.names = FALSE, stringsAsFactors = FALSE)
  if (isTRUE(check_minimum_columns) && ncol(imported) < 2L) {
    stop("The table was imported as a single column. Select the correct delimiter or verify the file contents.", call. = FALSE)
  }
  imported
}
