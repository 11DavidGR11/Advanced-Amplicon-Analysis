# =============================================================================
# Triple_A MultiAmplicon table validation and integration
# =============================================================================

multiamplicon_clean_names <- function(x) {
  x <- as.character(x)
  sub("^\\ufeff", "", x)
}

multiamplicon_validate_headers <- function(tables, file_names = names(tables)) {
  if (!is.list(tables) || length(tables) < 2L) {
    stop("Select at least two amplicon count tables.", call. = FALSE)
  }
  if (is.null(file_names) || length(file_names) != length(tables)) {
    file_names <- paste0("file_", seq_along(tables))
  }
  reference <- multiamplicon_clean_names(names(tables[[1L]]))
  if (!length(reference) || any(!nzchar(reference)) || anyDuplicated(reference)) {
    stop("The first table has empty or duplicated column headers.", call. = FALSE)
  }
  for (i in seq_along(tables)) {
    current <- multiamplicon_clean_names(names(tables[[i]]))
    if (any(!nzchar(current)) || anyDuplicated(current)) {
      stop(sprintf("%s has empty or duplicated column headers.", file_names[[i]]), call. = FALSE)
    }
    if (!setequal(current, reference)) {
      missing <- setdiff(reference, current)
      extra <- setdiff(current, reference)
      detail <- paste0(
        if (length(missing)) paste0(" Missing: ", paste(missing, collapse = ", "), ".") else "",
        if (length(extra)) paste0(" Extra: ", paste(extra, collapse = ", "), ".") else ""
      )
      stop(sprintf("Header mismatch in %s.%s", file_names[[i]], detail), call. = FALSE)
    }
  }
  invisible(reference)
}

multiamplicon_as_count <- function(x, column_name, file_name) {
  original <- x
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x[trimws(x) == ""] <- NA_character_
    suppressWarnings(value <- as.numeric(x))
  } else {
    value <- suppressWarnings(as.numeric(x))
  }
  invalid_text <- !is.na(original) & is.na(value)
  invalid_value <- !is.na(value) & (!is.finite(value) | value < 0 | abs(value - round(value)) > 1e-8)
  if (any(invalid_text | invalid_value)) {
    rows <- which(invalid_text | invalid_value)
    shown <- paste(head(rows, 8L), collapse = ", ")
    stop(sprintf(
      "%s: column '%s' is not a non-negative integer count column (problematic rows: %s%s).",
      file_name, column_name, shown, if (length(rows) > 8L) ", ..." else ""
    ), call. = FALSE)
  }
  value[is.na(value)] <- 0
  as.numeric(round(value))
}



# Infer likely count columns from a bounded sample instead of scanning every
# cell during upload. Full validation is still performed at integration time.
multiamplicon_guess_count_columns <- function(tables, sample_rows = 250L) {
  if (!is.list(tables) || !length(tables)) return(character())
  headers <- names(tables[[1L]])
  sample_rows <- max(1L, as.integer(sample_rows[[1L]]))
  headers[vapply(headers, function(column) {
    all(vapply(tables, function(tab) {
      value <- tab[[column]]
      if (length(value) > sample_rows) value <- value[seq_len(sample_rows)]
      tryCatch({
        multiamplicon_as_count(value, column, "column preview")
        TRUE
      }, error = function(e) FALSE)
    }, logical(1)))
  }, logical(1))]
}

multiamplicon_validate_and_combine <- function(
    tables,
    count_columns,
    file_names = names(tables),
    aggregate_duplicates = TRUE) {
  headers <- multiamplicon_validate_headers(tables, file_names)
  if (is.null(file_names) || length(file_names) != length(tables)) {
    file_names <- paste0("file_", seq_along(tables))
  }
  count_columns <- as.character(count_columns)
  if (!length(count_columns)) stop("Select at least one count column.", call. = FALSE)
  if (anyDuplicated(count_columns) || !all(count_columns %in% headers)) {
    stop("The selected count columns must be unique and present in every table.", call. = FALSE)
  }
  descriptor_columns <- setdiff(headers, count_columns)
  if (!length(descriptor_columns)) {
    stop("At least one non-count identifier or taxonomy column is required.", call. = FALSE)
  }

  normalized <- Map(function(tab, nm) {
    tab <- as.data.frame(tab, check.names = FALSE, stringsAsFactors = FALSE)
    names(tab) <- multiamplicon_clean_names(names(tab))

    # Tables may contain the same columns in a different order. Align every
    # table to the first file before validating counts or combining rows.
    tab <- tab[, headers, drop = FALSE]

    for (column in count_columns) {
      tab[[column]] <- multiamplicon_as_count(tab[[column]], column, nm)
    }
    tab
  }, tables, file_names)

  combined <- do.call(rbind, normalized)
  rownames(combined) <- NULL
  if (!isTRUE(aggregate_duplicates)) return(combined)

  # Aggregate only rows whose complete descriptor tuple is identical. This
  # prevents taxa with similar labels but different lineages from being merged.
  key_data <- combined[descriptor_columns]
  key_data[] <- lapply(key_data, function(z) {
    z <- as.character(z)
    z[is.na(z)] <- ""
    z
  })
  key <- do.call(paste, c(key_data, sep = "\r"))
  groups <- match(key, unique(key))
  first_rows <- match(seq_len(max(groups)), groups)
  result <- combined[first_rows, descriptor_columns, drop = FALSE]
  for (column in count_columns) {
    result[[column]] <- as.numeric(rowsum(combined[[column]], groups, reorder = FALSE)[, 1L])
  }
  result <- result[headers]
  rownames(result) <- NULL
  result
}
