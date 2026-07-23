# =============================================================================
# Triple_A canonical data model
# =============================================================================

# `%||%` is defined once in aaa_globals.R.

# Normalize sample identifiers only for matching. Original abundance-column names
# remain canonical in the dataset and in generated outputs. Common separators are
# treated as equivalent so identifiers such as S-1, S.1 and S_1 can be aligned.
aaa_normalize_sample_id <- function(x) {
  x <- trimws(enc2utf8(as.character(x)))
  x <- tolower(x)
  gsub("[[:space:]_.-]+", "", x, perl = TRUE)
}

aaa_match_sample_ids <- function(reference_ids, metadata_ids) {
  reference_ids <- as.character(reference_ids)
  metadata_ids <- as.character(metadata_ids)
  reference_keys <- aaa_normalize_sample_id(reference_ids)
  metadata_keys <- aaa_normalize_sample_id(metadata_ids)

  invalid_reference <- is.na(reference_keys) | !nzchar(reference_keys)
  invalid_metadata <- is.na(metadata_keys) | !nzchar(metadata_keys)
  if (any(invalid_reference)) {
    stop("Abundance sample identifiers cannot be missing or empty after normalisation.", call. = FALSE)
  }
  if (any(invalid_metadata)) {
    stop("Metadata sample identifiers cannot be missing or empty after normalisation.", call. = FALSE)
  }

  duplicated_reference <- unique(reference_keys[duplicated(reference_keys) | duplicated(reference_keys, fromLast = TRUE)])
  if (length(duplicated_reference)) {
    labels <- reference_ids[reference_keys %in% duplicated_reference]
    stop(
      "Abundance sample identifiers are ambiguous after normalisation: ",
      paste(unique(labels), collapse = ", "),
      call. = FALSE
    )
  }

  duplicated_metadata <- unique(metadata_keys[duplicated(metadata_keys) | duplicated(metadata_keys, fromLast = TRUE)])
  if (length(duplicated_metadata)) {
    labels <- metadata_ids[metadata_keys %in% duplicated_metadata]
    stop(
      "Metadata sample identifiers are ambiguous after normalisation: ",
      paste(unique(labels), collapse = ", "),
      call. = FALSE
    )
  }

  positions <- match(reference_keys, metadata_keys)
  data.frame(
    reference_id = reference_ids,
    metadata_id = metadata_ids[positions],
    metadata_position = positions,
    matched = !is.na(positions),
    exact = !is.na(positions) & reference_ids == metadata_ids[positions],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Build a sample design from a metadata column
#'
#' Each sample is assigned to the group named in `group_column`, matched by
#' sample identifier rather than by column position. Groups may therefore be
#' interleaved in the abundance table and may contain different numbers of
#' replicates.
#'
#' @param sample_columns Abundance sample columns, in table order.
#' @param metadata Metadata table, one row per sample.
#' @param sample_id_column Metadata column holding the sample identifiers.
#' @param group_column Metadata column holding the experimental group.
#' @param group_levels Optional explicit group order; defaults to order of
#'   first appearance in the abundance table.
#' @return A data frame with Sample_column, Treatment and Replicate.
aaa_sample_design_from_metadata <- function(sample_columns, metadata,
                                            sample_id_column, group_column,
                                            group_levels = NULL) {
  sample_columns <- as.character(sample_columns)
  if (!length(sample_columns)) stop("No sample columns were supplied.", call. = FALSE)
  if (is.null(metadata) || !is.data.frame(metadata) || !nrow(metadata)) {
    stop("A metadata table is required to build the sample design.", call. = FALSE)
  }
  for (column in c(sample_id_column, group_column)) {
    if (is.null(column) || !nzchar(column) || !column %in% names(metadata)) {
      stop("Metadata column was not found: ", column %||% "<none>", call. = FALSE)
    }
  }
  if (identical(sample_id_column, group_column)) {
    stop(
      "The grouping column must be different from the sample-identifier column.",
      call. = FALSE
    )
  }

  matches <- aaa_match_sample_ids(
    sample_columns, as.character(metadata[[sample_id_column]])
  )
  if (any(!matches$matched)) {
    stop(
      "Metadata are missing for these abundance samples: ",
      paste(matches$reference_id[!matches$matched], collapse = ", "),
      call. = FALSE
    )
  }

  groups <- trimws(as.character(metadata[[group_column]][matches$metadata_position]))
  missing_group <- is.na(groups) | !nzchar(groups)
  if (any(missing_group)) {
    stop(
      "These samples have no value in the grouping column '", group_column, "': ",
      paste(sample_columns[missing_group], collapse = ", "),
      call. = FALSE
    )
  }

  # Group ORDER is taken from the metadata file, not from the order in which the
  # samples happen to appear in the abundance table. Otherwise reordering the
  # abundance columns would silently change which group acts as the reference
  # level, flipping the direction of every pairwise comparison and the order of
  # legends and heatmap columns. The metadata is a deliberate declaration; the
  # column order of the abundance table is incidental.
  if (is.null(group_levels)) {
    group_levels <- if (is.factor(metadata[[group_column]])) {
      levels(droplevels(metadata[[group_column]]))
    } else {
      declared <- trimws(as.character(metadata[[group_column]]))
      unique(declared[!is.na(declared) & nzchar(declared)])
    }
  }
  group_levels <- as.character(group_levels)
  # Any group present in the samples but absent from the declared order is
  # appended rather than dropped, so no sample can lose its label.
  group_levels <- c(group_levels, setdiff(unique(groups), group_levels))

  data.frame(
    Sample_column = sample_columns,
    # Kept as a factor so the declared order survives into samples_name and
    # from there into comparison directions, legends and column ordering.
    Treatment = factor(groups, levels = group_levels),
    # Replicate is an ordinal within its own group, not a position in the table.
    Replicate = stats::ave(seq_along(groups), groups, FUN = seq_along),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

aaa_new_dataset <- function(abundance, sample_design, metadata = NULL,
                            metadata_roles = NULL, source = list()) {
  abundance <- as.data.frame(abundance,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_design <- as.data.frame(sample_design,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(abundance) == 0L) stop("The abundance table contains no rows.")
  taxonomy_column <- intersect(c("Taxonomy", "taxonomy"), names(abundance))[1]
  if (is.na(taxonomy_column)) stop("The abundance table requires a Taxonomy column.")
  if (taxonomy_column != "Taxonomy") names(abundance)[names(abundance) == taxonomy_column] <- "Taxonomy"

  required_design <- c("Sample_column", "Treatment", "Replicate")
  missing_design <- setdiff(required_design, names(sample_design))
  if (length(missing_design)) {
    stop("Sample design is missing: ", paste(missing_design, collapse = ", "))
  }
  if (nrow(sample_design) == 0L) stop("The sample design contains no samples.")
  if (anyDuplicated(sample_design$Sample_column)) stop("Sample identifiers must be unique.")
  if (anyNA(sample_design$Sample_column) || any(!nzchar(trimws(sample_design$Sample_column)))) {
    stop("Sample identifiers cannot be missing or empty.")
  }

  missing_samples <- setdiff(sample_design$Sample_column, names(abundance))
  if (length(missing_samples)) {
    stop(
      "Sample columns are absent from the abundance table: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  non_numeric <- sample_design$Sample_column[
    !vapply(abundance[sample_design$Sample_column], is.numeric, logical(1))
  ]
  if (length(non_numeric)) stop("Non-numeric sample columns: ", paste(non_numeric, collapse = ", "))
  if (any(vapply(abundance[sample_design$Sample_column], function(x) any(x < 0, na.rm = TRUE), logical(1)))) {
    stop("Abundance values cannot be negative.")
  }

  if (!"FeatureID" %in% names(abundance)) {
    abundance$FeatureID <- sprintf("Feature_%06d", seq_len(nrow(abundance)))
  }
  abundance <- abundance[, c("FeatureID", "Taxonomy", setdiff(names(abundance), c("FeatureID", "Taxonomy"))), drop = FALSE]

  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(metadata_roles) && nrow(metadata_roles)) {
      metadata_roles <- as.data.frame(metadata_roles, stringsAsFactors = FALSE)
      id_cols <- metadata_roles$Column[metadata_roles$Role == "identifier"]
      if (length(id_cols) != 1L) stop("Metadata must have exactly one sample identifier role.")
      id_col <- id_cols[[1]]
      if (!id_col %in% names(metadata)) stop("Metadata identifier column was not found: ", id_col)
      metadata_ids <- as.character(metadata[[id_col]])
      sample_matches <- aaa_match_sample_ids(sample_design$Sample_column, metadata_ids)
      if (any(!sample_matches$matched)) {
        missing_metadata <- sample_matches$reference_id[!sample_matches$matched]
        stop("Metadata are missing samples: ", paste(missing_metadata, collapse = ", "))
      }

      # Canonicalise only matched metadata identifiers to the exact abundance
      # column names. This keeps all downstream joins strict and deterministic
      # while preserving unrelated metadata rows unchanged.
      metadata[[id_col]][sample_matches$metadata_position] <- sample_matches$reference_id
    }
  }

  structure(list(
    abundance = abundance,
    sample_design = sample_design,
    metadata = metadata,
    metadata_roles = metadata_roles,
    source = source,
    schema_version = "2.0"
  ), class = c("Triple_A_dataset", "list"))
}

aaa_validate_dataset <- function(dataset) {
  if (!inherits(dataset, "Triple_A_dataset")) {
    stop("'dataset' must be a Triple_A_dataset created by aaa_new_dataset().")
  }
  invisible(aaa_new_dataset(
    dataset$abundance, dataset$sample_design,
    dataset$metadata, dataset$metadata_roles,
    dataset$source
  ))
}

aaa_dataset_environment <- function(dataset, variables = NULL) {
  aaa_validate_dataset(dataset)
  if (is.null(dataset$metadata)) stop("No metadata table is attached to this dataset.")
  roles <- dataset$metadata_roles
  if (is.null(roles) || !nrow(roles)) stop("Metadata roles have not been defined.")
  id_col <- roles$Column[roles$Role == "identifier"]
  if (length(id_col) != 1L) stop("Exactly one metadata identifier is required.")
  if (is.null(variables)) variables <- roles$Column[roles$Role == "environmental_variable"]
  missing <- setdiff(variables, names(dataset$metadata))
  if (length(missing)) stop("Environmental variables not found: ", paste(missing, collapse = ", "))
  out <- dataset$metadata[, c(id_col, variables), drop = FALSE]
  names(out)[1] <- "Sample_column"
  out
}
