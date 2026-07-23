# Triple_A pre-flight validation -----------------------------------------
aaa_validation_row <- function(check, status, message, blocking = status %in% c("Error", "Pending")) {
  data.frame(Check = check, Status = status, Message = message, Blocking = blocking, stringsAsFactors = FALSE)
}

#' @param group_sizes Optional named vector of samples per group. When supplied,
#'   the design is validated from it rather than from the replicate keyword,
#'   which only describes balanced designs built from consecutive blocks.
aaa_validate_preflight <- function(input_file = NULL, data = NULL, sample_identifier = "", sample_columns = NULL, treatments = character(),
                                   replicates = "none", analyses = character(), outputs = character(),
                                   environmental_file = NULL, group_sizes = NULL) {
  rows <- list()
  rows[[1]] <- if (is.null(input_file) || !file.exists(input_file)) {
    aaa_validation_row("Input file", "Pending", "Load a supported abundance table.")
  } else {
    aaa_validation_row("Input file", "Valid", paste("Readable:", basename(input_file)))
  }
  rows[[2]] <- if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
    aaa_validation_row("Excel/data table", "Pending", "No parsed data are available.")
  } else {
    aaa_validation_row("Excel/data table", "Valid", sprintf("%d rows and %d columns loaded.", nrow(data), ncol(data)))
  }
  sample_cols <- if (!is.null(sample_columns)) {
    intersect(as.character(sample_columns), names(data))
  } else if (is.null(data) || !nzchar(sample_identifier)) {
    character()
  } else {
    grep(sample_identifier, names(data), value = TRUE, fixed = TRUE)
  }
  rows[[3]] <- if (!length(sample_cols)) {
    aaa_validation_row("Sample columns", "Error", "No sample columns match the configured identifier.")
  } else {
    aaa_validation_row("Sample columns", "Valid", sprintf("%d sample columns detected.", length(sample_cols)))
  }
  rows[[4]] <- if (length(treatments) < 2) {
    aaa_validation_row("Treatments", "Error", "At least two treatments are required.")
  } else {
    aaa_validation_row("Treatments", "Valid", paste(length(treatments), "treatments detected."))
  }
  rows[[5]] <- if (!is.null(group_sizes) && length(group_sizes)) {
    # An explicit design assigns each sample to a group by name, so there is no
    # divisibility requirement and groups may differ in size.
    summary_text <- paste(
      sprintf("%s=%d", names(group_sizes), as.integer(group_sizes)),
      collapse = ", "
    )
    if (any(as.integer(group_sizes) < 1L)) {
      aaa_validation_row("Replicates", "Error", paste("Empty treatment group:", summary_text))
    } else if (length(unique(as.integer(group_sizes))) == 1L) {
      aaa_validation_row("Replicates", "Valid", paste("Balanced design:", summary_text))
    } else {
      aaa_validation_row("Replicates", "Valid", paste("Unbalanced design (supported):", summary_text))
    }
  } else {
    expected <- c(none = 1, duplicate = 2, triplicate = 3, quadruplicate = 4, quintuplicate = 5)[replicates]
    if (is.na(expected)) {
      aaa_validation_row("Replicates", "Error", "Unknown replicate design.")
    } else if (length(sample_cols) %% expected != 0) {
      aaa_validation_row("Replicates", "Warning", "The number of sample columns is not divisible by the selected replicate count.", FALSE)
    } else {
      aaa_validation_row("Replicates", "Valid", paste("Design:", replicates))
    }
  }
  missing_count <- if (!length(sample_cols) || is.null(data)) NA_integer_ else sum(is.na(data[sample_cols]))
  rows[[6]] <- if (is.na(missing_count)) {
    aaa_validation_row("Missing values", "Pending", "Waiting for sample data.")
  } else if (missing_count > 0) {
    aaa_validation_row("Missing values", "Warning", paste(missing_count, "missing abundance values will require handling."), FALSE)
  } else {
    aaa_validation_row("Missing values", "Valid", "No missing abundance values detected.")
  }
  tax_candidates <- if (is.null(data)) character() else grep("tax|lineage|classification", names(data), ignore.case = TRUE, value = TRUE)
  rows[[7]] <- if (!length(tax_candidates)) {
    aaa_validation_row("Taxonomy format", "Warning", "No taxonomy-like column name was detected; verify the first identifier column.", FALSE)
  } else {
    aaa_validation_row("Taxonomy format", "Valid", paste("Detected:", tax_candidates[1]))
  }
  rows[[8]] <- if ("rda" %in% analyses && (is.null(environmental_file) || !file.exists(environmental_file))) {
    aaa_validation_row("Environmental variables", "Error", "RDA requires an environmental metadata file.")
  } else if (!is.null(environmental_file) && file.exists(environmental_file)) {
    aaa_validation_row("Environmental variables", "Valid", paste("Readable:", basename(environmental_file)))
  } else {
    aaa_validation_row("Environmental variables", "Valid", "Not required by the selected analyses.")
  }
  rows[[9]] <- if (!length(analyses)) {
    aaa_validation_row("Analyses", "Error", "Select at least one analysis.")
  } else {
    aaa_validation_row("Analyses", "Valid", paste(length(analyses), "analysis modules selected."))
  }
  rows[[10]] <- if (!length(analyses)) {
    aaa_validation_row("Automatic outputs", "Pending", "Outputs will be generated after selecting an analysis.")
  } else {
    aaa_validation_row(
      "Automatic outputs",
      "Valid",
      paste(length(outputs), "standard outputs will be generated automatically."),
      FALSE
    )
  }
  missing_packages <- tryCatch(aaa_missing_packages(analyses = analyses, input_files = c(input_file, environmental_file), include_core = TRUE), error = function(e) conditionMessage(e))
  rows[[11]] <- if (length(missing_packages)) {
    aaa_validation_row("Packages", "Error", paste("Missing:", paste(missing_packages, collapse = ", ")))
  } else {
    aaa_validation_row("Packages", "Valid", "All required packages are available.")
  }
  diagnostics <- tryCatch(aaa_environment_diagnostics(getwd(), analyses), error = function(e) NULL)
  rows[[12]] <- if (is.null(diagnostics)) {
    aaa_validation_row("Environment", "Warning", "Environment diagnostics could not be completed.", FALSE)
  } else if (any(diagnostics$Status == "Error")) {
    aaa_validation_row("Environment", "Error", paste(diagnostics$Message[diagnostics$Status == "Error"], collapse = "; "))
  } else {
    aaa_validation_row("Environment", "Valid", "R, packages, file access and execution services are ready.", FALSE)
  }
  do.call(rbind, rows)
}

# Scientific input and analysis eligibility checks ----------------------
#
# Every check returns the same shape. `reason` states what is wrong in one
# sentence; `guidance` lists the concrete steps that unlock the analysis, so the
# interface can tell the user what to do instead of only what failed.
aaa_eligibility_result <- function(available, reason = NULL, details = list(),
                                   guidance = character()) {
  list(
    available = isTRUE(available),
    reason = if (isTRUE(available)) NULL else as.character(reason)[1L],
    guidance = if (isTRUE(available)) character() else as.character(guidance),
    details = details
  )
}

# PLS-DA requires declared biological replication. Keeping this rule in the
# engine makes the availability decision testable independently of Shiny.
# `sample_columns` and `treatments` are optional so existing callers that only
# know the replicate design keep working; when supplied, the check mirrors every
# stop() inside aaa_plsda_analysis() so the interface blocks what would fail.
#' @param group_sizes Optional named vector of samples per group. When supplied
#'   it takes precedence over `replicates`, which can only describe a balanced
#'   design built from consecutive blocks of columns.
aaa_check_plsda_eligibility <- function(replicates = "none",
                                        sample_columns = NULL,
                                        treatments = NULL,
                                        label = "PLS-DA",
                                        group_sizes = NULL) {
  if (!is.null(group_sizes) && length(group_sizes)) {
    group_sizes <- stats::setNames(as.integer(group_sizes), names(group_sizes))
    if (length(group_sizes) < 2L) {
      return(aaa_eligibility_result(
        FALSE,
        paste0(label, " is disabled: at least two treatment groups are required."),
        guidance = c(
          "Declare at least two distinct groups in the experimental design.",
          paste0(label, " is a supervised classifier: with a single group there is nothing to discriminate.")
        )
      ))
    }
    small <- names(group_sizes)[group_sizes < 2L]
    if (length(small)) {
      return(aaa_eligibility_result(
        FALSE,
        paste0(
          label, " is disabled: these groups have a single sample: ",
          paste(small, collapse = ", "), "."
        ),
        guidance = c(
          "Give every group at least two samples, or exclude the listed groups.",
          "Cross-validation must be able to hold out a sample and still leave the group represented.",
          "Group sizes may differ from each other.",
          "Descriptive ordinations (PCA, PCoA, NMDS) stay available with unreplicated designs."
        )
      ))
    }
    if (sum(group_sizes) < 3L) {
      return(aaa_eligibility_result(
        FALSE,
        paste0(label, " is disabled: at least three samples are required."),
        guidance = c("Select at least three sample columns.")
      ))
    }
    return(aaa_eligibility_result(
      TRUE,
      details = list(groups = length(group_sizes), replicates_per_group = group_sizes)
    ))
  }

  replicates <- tolower(trimws(as.character(replicates %||% "none")[1L]))
  replicated <- replicates %in% c("duplicate", "triplicate", "quadruplicate", "quintuplicate")
  if (!replicated) {
    return(aaa_eligibility_result(
      FALSE,
      paste0(label, " is disabled: the design has no biological replicates."),
      details = list(replicates = replicates),
      guidance = c(
        "Open 'Data and metadata' and set the replicate design to duplicate or higher.",
        paste0(label, " estimates a cross-validated classifier, so each group needs at least two samples to hold one out."),
        "Descriptive ordinations (PCA, PCoA, NMDS) stay available with unreplicated designs."
      )
    ))
  }

  if (!is.null(treatments)) {
    groups <- unique(trimws(as.character(treatments)))
    groups <- groups[nzchar(groups)]
    if (length(groups) < 2L) {
      return(aaa_eligibility_result(
        FALSE,
        paste0(label, " is disabled: at least two treatment groups are required."),
        details = list(groups = length(groups)),
        guidance = c(
          "Declare at least two distinct treatment names in 'Data and metadata'.",
          paste0(label, " is a supervised classifier: with a single group there is nothing to discriminate.")
        )
      ))
    }
  }

  if (!is.null(sample_columns)) {
    n_samples <- length(as.character(sample_columns))
    if (n_samples < 3L) {
      return(aaa_eligibility_result(
        FALSE,
        paste0(label, " is disabled: at least three samples are required."),
        details = list(samples = n_samples),
        guidance = c(
          "Select at least three sample columns in 'Data and metadata'.",
          "Cross-validation cannot form training and validation folds with fewer samples."
        )
      ))
    }
  }

  aaa_eligibility_result(TRUE, details = list(replicates = replicates))
}

# Ordination and diversity need enough samples to place points in two
# dimensions; the group-based tests inside the module need enough groups. The
# module is only blocked outright when no ordination at all could be produced.
aaa_check_community_structure_eligibility <- function(sample_columns,
                                                      treatments = NULL) {
  sample_columns <- as.character(sample_columns %||% character())

  if (!length(sample_columns)) {
    return(aaa_eligibility_result(
      FALSE,
      "Community structure is disabled: no valid sample columns are selected.",
      guidance = c(
        "Import an abundance table in 'Data and metadata'.",
        "Confirm the sample-column selection so at least three samples are included."
      )
    ))
  }

  if (length(sample_columns) < 3L) {
    return(aaa_eligibility_result(
      FALSE,
      paste0(
        "Community structure is disabled: two-dimensional ordination needs at least three samples, but ",
        length(sample_columns), " are selected."
      ),
      details = list(samples = length(sample_columns)),
      guidance = c(
        "Select at least three sample columns in 'Data and metadata'.",
        "PCA, PCoA and NMDS place samples on two axes, which is undefined with fewer than three points.",
        "With one or two samples, use Top abundance for a descriptive summary instead."
      )
    ))
  }

  groups <- unique(trimws(as.character(treatments %||% character())))
  groups <- groups[nzchar(groups)]

  aaa_eligibility_result(
    TRUE,
    details = list(samples = length(sample_columns), groups = length(groups))
  )
}

# Validates a metadata-driven experimental design before it is built, so the
# interface can explain what is wrong instead of surfacing a raw engine error.
aaa_check_sample_design_eligibility <- function(sample_columns, metadata,
                                                sample_id_column, group_column) {
  sample_columns <- as.character(sample_columns %||% character())

  if (!length(sample_columns)) {
    return(aaa_eligibility_result(
      FALSE, "The experimental design is incomplete: no sample columns are selected.",
      guidance = c("Import an abundance table and confirm the sample-column selection.")
    ))
  }
  if (is.null(metadata) || !is.data.frame(metadata) || !nrow(metadata)) {
    return(aaa_eligibility_result(
      FALSE, "The experimental design is incomplete: no metadata file is loaded.",
      guidance = c(
        "Upload a metadata file with one row per sample.",
        "It needs a sample-identifier column and a column naming the group of each sample.",
        "Alternatively, switch the design source back to 'consecutive blocks'."
      )
    ))
  }
  if (is.null(group_column) || !nzchar(group_column) || !group_column %in% names(metadata)) {
    return(aaa_eligibility_result(
      FALSE, "The experimental design is incomplete: no grouping column is selected.",
      guidance = c(
        "Choose the metadata column that names the experimental group of each sample.",
        "Typical names are Treatment, Group or Condition."
      )
    ))
  }
  if (identical(group_column, sample_id_column)) {
    return(aaa_eligibility_result(
      FALSE, "The grouping column and the sample-identifier column are the same.",
      guidance = c(
        "Select a different column for the group.",
        "The identifier names each sample; the grouping column says which samples belong together."
      )
    ))
  }

  design <- tryCatch(
    aaa_sample_design_from_metadata(
      sample_columns, metadata, sample_id_column, group_column
    ),
    error = function(e) e
  )
  if (inherits(design, "error")) {
    return(aaa_eligibility_result(
      FALSE,
      paste0("The experimental design could not be built: ", conditionMessage(design)),
      guidance = c(
        "Add one metadata row per selected sample, with a non-empty value in the grouping column.",
        "Identifier differences in case, spaces, '-', '.' and '_' are tolerated; anything else must match."
      )
    ))
  }

  # Groups declared in the metadata but absent from the selected samples are not
  # part of this design; counting them would report zero-size groups and turn
  # the smallest group size into 0.
  counts <- aaa_treatment_counts(design$Treatment)
  if (length(counts) < 1L) {
    return(aaa_eligibility_result(
      FALSE, "The experimental design has no groups.",
      guidance = c("Check that the grouping column contains values.")
    ))
  }

  # A column with one distinct value per sample is almost always the wrong
  # choice: it is usually a second identifier or a continuous measurement, and
  # it leaves every sample alone in its own group, which disables every analysis
  # that needs replication. Caught here rather than as a puzzling downstream
  # "no replicates" message.
  if (length(counts) == length(sample_columns) && length(sample_columns) > 1L) {
    return(aaa_eligibility_result(
      FALSE,
      paste0(
        "The grouping column '", group_column, "' has a different value for every sample, ",
        "so each sample would form its own group of one."
      ),
      guidance = c(
        "Choose the column that names the experimental condition shared by several samples, not one that identifies each sample individually.",
        "A continuous measurement is not a grouping column: those belong in the environmental-variable role.",
        "If the metadata has no such column, add one, for example Group or Condition, with one value per experimental condition."
      )
    ))
  }

  aaa_eligibility_result(
    TRUE,
    details = list(
      groups = names(counts),
      replicates = as.integer(counts),
      balanced = length(unique(as.integer(counts))) == 1L,
      minimum_replicates = as.integer(min(counts))
    )
  )
}

aaa_check_top_abundance_eligibility <- function(sample_columns) {
  sample_columns <- as.character(sample_columns %||% character())
  if (!length(sample_columns)) {
    return(aaa_eligibility_result(
      FALSE,
      "Top abundance is disabled: no valid sample columns are selected.",
      guidance = c(
        "Import an abundance table in 'Data and metadata'.",
        "Check the sample-identifier pattern so the sample columns are detected."
      )
    ))
  }
  aaa_eligibility_result(TRUE, details = list(samples = length(sample_columns)))
}

aaa_check_functional_potential_eligibility <- function(sample_columns,
                                                       selected_functions) {
  sample_columns <- as.character(sample_columns %||% character())
  selected_functions <- as.character(selected_functions %||% character())

  if (!length(sample_columns)) {
    return(aaa_eligibility_result(
      FALSE,
      "Functional potential is disabled: no valid sample columns are selected.",
      guidance = c("Import an abundance table in 'Data and metadata' before selecting functions.")
    ))
  }
  if (!length(selected_functions)) {
    return(aaa_eligibility_result(
      FALSE,
      "Functional potential is disabled: no biological function is selected.",
      guidance = c(
        "Pick at least one function in the 'Functional analysis' panel, or apply one of the presets.",
        "Use the search box to find a pathway by name if the thematic list is long."
      )
    ))
  }
  aaa_eligibility_result(
    TRUE,
    details = list(functions = length(selected_functions))
  )
}

# ANCOM-BC2 and MaAsLin2 model taxa against metadata columns, so they need a
# metadata table with at least one usable predictor matched to the samples.
aaa_check_taxon_association_eligibility <- function(sample_columns,
                                                    metadata,
                                                    sample_id_column,
                                                    variables) {
  sample_columns <- as.character(sample_columns %||% character())
  variables <- as.character(variables %||% character())

  if (!length(sample_columns)) {
    return(aaa_eligibility_result(
      FALSE,
      "Taxon-association analyses are disabled: no valid abundance sample columns are selected.",
      guidance = c("Import an abundance table in 'Data and metadata'.")
    ))
  }
  if (is.null(metadata) || !is.data.frame(metadata) || !nrow(metadata)) {
    return(aaa_eligibility_result(
      FALSE,
      "Taxon-association analyses are disabled: no sample metadata are loaded.",
      guidance = c(
        "Upload a metadata file in 'Data and metadata', with one row per sample.",
        "ANCOM-BC2 and MaAsLin2 model each taxon against metadata columns, so a predictor is mandatory."
      )
    ))
  }
  if (is.null(sample_id_column) || !nzchar(sample_id_column) ||
    !sample_id_column %in% names(metadata)) {
    return(aaa_eligibility_result(
      FALSE,
      "Taxon-association analyses are disabled: no valid metadata sample-ID column is selected.",
      guidance = c(
        "Assign the 'identifier' role to the metadata column holding the sample names.",
        "Its values must match the abundance sample columns (differences in case, spaces, '-', '.' and '_' are tolerated)."
      )
    ))
  }
  if (!length(intersect(variables, names(metadata)))) {
    return(aaa_eligibility_result(
      FALSE,
      "Taxon-association analyses are disabled: no metadata column is classified as an experimental factor or environmental variable.",
      guidance = c(
        "Give at least one metadata column the 'experimental_factor' role (categorical predictor) or 'environmental_variable' role (continuous predictor).",
        "The identifier column alone is not a predictor."
      )
    ))
  }

  metadata_ids <- trimws(as.character(metadata[[sample_id_column]]))
  sample_matches <- tryCatch(
    aaa_match_sample_ids(sample_columns, metadata_ids),
    error = function(e) e
  )
  if (inherits(sample_matches, "error")) {
    return(aaa_eligibility_result(
      FALSE,
      paste0("Taxon-association analyses are disabled: ", conditionMessage(sample_matches)),
      guidance = c("Make the sample identifiers unique and non-empty in both the abundance table and the metadata.")
    ))
  }
  if (any(!sample_matches$matched)) {
    missing_ids <- sample_matches$reference_id[!sample_matches$matched]
    return(aaa_eligibility_result(
      FALSE,
      paste0(
        "Taxon-association analyses are disabled: metadata are missing for these abundance samples: ",
        paste(missing_ids, collapse = ", "), "."
      ),
      guidance = c(
        "Add one metadata row per listed sample, or deselect those samples.",
        "Check for typographical differences between the two files."
      )
    ))
  }

  aaa_eligibility_result(
    TRUE,
    details = list(
      variables = intersect(variables, names(metadata)),
      matched_samples = length(sample_columns)
    )
  )
}

# envfit, partial RDA, dbRDA and variance partitioning accept factors as well as
# continuous variables, so they are less restrictive than RDA proper.
aaa_check_constrained_ordination_eligibility <- function(sample_columns,
                                                          metadata,
                                                          sample_id_column,
                                                          variables,
                                                          minimum_samples = 4L,
                                                          label = "Constrained ordination") {
  base <- aaa_check_taxon_association_eligibility(
    sample_columns = sample_columns,
    metadata = metadata,
    sample_id_column = sample_id_column,
    variables = variables
  )
  if (!isTRUE(base$available)) {
    return(aaa_eligibility_result(
      FALSE,
      sub("^Taxon-association analyses are disabled", paste0(label, " is disabled"), base$reason),
      guidance = base$guidance
    ))
  }

  if (length(as.character(sample_columns)) < minimum_samples) {
    return(aaa_eligibility_result(
      FALSE,
      paste0(
        label, " is disabled: at least ", minimum_samples,
        " matched samples are required, but ", length(as.character(sample_columns)), " are selected."
      ),
      guidance = c(
        paste0("Select at least ", minimum_samples, " sample columns with complete metadata."),
        "Constrained ordination fits environmental predictors and needs residual degrees of freedom."
      )
    ))
  }

  aaa_eligibility_result(TRUE, details = base$details)
}

aaa_check_variance_partitioning_eligibility <- function(sample_columns,
                                                         metadata,
                                                         sample_id_column,
                                                         environmental_variables,
                                                         experimental_factors) {
  base <- aaa_check_constrained_ordination_eligibility(
    sample_columns = sample_columns,
    metadata = metadata,
    sample_id_column = sample_id_column,
    variables = unique(c(environmental_variables, experimental_factors)),
    label = "Variance partitioning"
  )
  if (!isTRUE(base$available)) {
    return(base)
  }

  environmental_variables <- intersect(as.character(environmental_variables %||% character()), names(metadata))
  experimental_factors <- intersect(as.character(experimental_factors %||% character()), names(metadata))

  if (!length(environmental_variables) || !length(experimental_factors)) {
    return(aaa_eligibility_result(
      FALSE,
      "Variance partitioning is disabled: it needs one environmental variable set AND one experimental factor set.",
      guidance = c(
        "Give at least one metadata column the 'environmental_variable' role (continuous predictor).",
        "Give at least one different metadata column the 'experimental_factor' role (experimental design predictor).",
        "The analysis splits explained variation between those two sets, so both must be non-empty."
      )
    ))
  }

  aaa_eligibility_result(
    TRUE,
    details = list(
      environmental_variables = environmental_variables,
      experimental_factors = experimental_factors
    )
  )
}

aaa_abundance_numeric_matrix <- function(data, sample_columns) {
  if (is.null(data) || !is.data.frame(data)) {
    return(NULL)
  }
  sample_columns <- intersect(as.character(sample_columns), names(data))
  if (!length(sample_columns)) {
    return(NULL)
  }

  converted <- lapply(data[sample_columns], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  matrix <- as.matrix(as.data.frame(converted, check.names = FALSE))
  colnames(matrix) <- sample_columns
  matrix
}

aaa_validate_abundance_nature <- function(
  data,
  sample_columns,
  abundance_type = c("proportion", "percentage", "counts"),
  proportion_sum_tolerance = 0.02,
  percentage_sum_tolerance = 2
) {
  abundance_type <- match.arg(abundance_type)
  values <- aaa_abundance_numeric_matrix(data, sample_columns)

  if (is.null(values) || !length(values)) {
    return(aaa_eligibility_result(
      FALSE,
      "Abundance format cannot be verified until valid sample columns are selected."
    ))
  }

  if (anyNA(values)) {
    return(aaa_eligibility_result(
      FALSE,
      "Abundance format cannot be verified because one or more selected sample columns contain missing or non-numeric values."
    ))
  }

  if (any(!is.finite(values))) {
    return(aaa_eligibility_result(
      FALSE,
      "Abundance values must be finite numeric values."
    ))
  }

  if (any(values < 0)) {
    return(aaa_eligibility_result(
      FALSE,
      "Abundance values cannot be negative."
    ))
  }

  totals <- colSums(values)
  maximum <- max(values, na.rm = TRUE)

  if (identical(abundance_type, "proportion")) {
    invalid_range <- maximum > 1 + sqrt(.Machine$double.eps)
    invalid_totals <- abs(totals - 1) > proportion_sum_tolerance

    if (invalid_range || any(invalid_totals)) {
      bad <- names(totals)[invalid_totals]
      reason <- paste0(
        "The table was declared as proportions, but values must be between 0 and 1 and each complete sample column must sum to approximately 1."
      )
      if (length(bad)) {
        reason <- paste0(reason, " Incompatible sample columns: ", paste(bad, collapse = ", "), ".")
      }
      return(aaa_eligibility_result(FALSE, reason, list(column_totals = totals)))
    }
  }

  if (identical(abundance_type, "percentage")) {
    invalid_range <- maximum > 100 + sqrt(.Machine$double.eps)
    invalid_totals <- abs(totals - 100) > percentage_sum_tolerance

    if (invalid_range || any(invalid_totals)) {
      bad <- names(totals)[invalid_totals]
      reason <- paste0(
        "The table was declared as percentages, but values must be between 0 and 100 and each complete sample column must sum to approximately 100."
      )
      if (length(bad)) {
        reason <- paste0(reason, " Incompatible sample columns: ", paste(bad, collapse = ", "), ".")
      }
      return(aaa_eligibility_result(FALSE, reason, list(column_totals = totals)))
    }
  }

  if (identical(abundance_type, "counts")) {
    integer_like <- abs(values - round(values)) <= sqrt(.Machine$double.eps)
    if (!all(integer_like)) {
      return(aaa_eligibility_result(
        FALSE,
        "The table was declared as raw read counts, but it contains non-integer values."
      ))
    }
    if (any(totals <= 0)) {
      return(aaa_eligibility_result(
        FALSE,
        "The table was declared as raw read counts, but at least one selected sample has a total count of zero."
      ))
    }
  }

  aaa_eligibility_result(
    TRUE,
    details = list(
      abundance_type = abundance_type,
      minimum = min(values),
      maximum = maximum,
      column_totals = totals
    )
  )
}

#' @param group_sizes Optional named vector of samples per group. When supplied
#'   it takes precedence over `replicates`, which can only describe a balanced
#'   design built from consecutive blocks of columns.
aaa_check_differential_abundance_eligibility <- function(
  sample_columns,
  treatments,
  replicates,
  minimum_replicates = 2L,
  group_sizes = NULL
) {
  sample_columns <- as.character(sample_columns %||% character())
  treatments <- as.character(treatments %||% character())

  if (!is.null(group_sizes) && length(group_sizes)) {
    group_sizes <- stats::setNames(as.integer(group_sizes), names(group_sizes))
    small <- names(group_sizes)[group_sizes < minimum_replicates]
    if (length(small)) {
      return(aaa_eligibility_result(
        FALSE,
        paste0(
          "Differential abundance is disabled: these groups have fewer than ",
          minimum_replicates, " samples: ",
          paste(sprintf("%s (%d)", small, group_sizes[small]), collapse = ", "), "."
        ),
        guidance = c(
          paste0("Give every group at least ", minimum_replicates, " samples, or exclude the listed groups."),
          "Wilcoxon and t tests estimate within-group variability, which is undefined with one sample per group.",
          "Group sizes may differ from each other; only the minimum matters."
        )
      ))
    }
    if (length(group_sizes) < 2L) {
      return(aaa_eligibility_result(
        FALSE,
        "Differential abundance is disabled: at least two treatment groups are required.",
        guidance = c("Declare at least two distinct groups in the experimental design.")
      ))
    }
    return(aaa_eligibility_result(
      TRUE,
      details = list(groups = length(group_sizes), replicates_per_group = group_sizes)
    ))
  }

  replicate_n <- tryCatch(aaa_n_replicates(replicates), error = function(e) NA_integer_)

  if (!length(sample_columns)) {
    return(aaa_eligibility_result(
      FALSE,
      "Differential abundance is disabled: no valid sample columns are selected.",
      guidance = c(
        "Import an abundance table in 'Data and metadata'.",
        "Check the sample-identifier pattern so the sample columns are detected."
      )
    ))
  }
  if (length(treatments) < 2L) {
    return(aaa_eligibility_result(
      FALSE,
      "Differential abundance is disabled: at least two treatment groups are required.",
      guidance = c(
        "Declare at least two treatment names in 'Data and metadata'.",
        "The analysis compares taxa between pairs of groups, so a single group has nothing to compare."
      )
    ))
  }
  if (is.na(replicate_n) || replicate_n < minimum_replicates) {
    return(aaa_eligibility_result(
      FALSE,
      paste0(
        "Differential abundance is disabled: the design has no biological replicates. At least ",
        minimum_replicates, " samples per treatment are required."
      ),
      guidance = c(
        "Set the replicate design to duplicate or higher in 'Data and metadata'.",
        "Wilcoxon and t tests estimate within-group variability, which is undefined with one sample per group.",
        "Descriptive abundance summaries (Top abundance) stay available."
      )
    ))
  }
  if (length(sample_columns) != length(treatments) * replicate_n) {
    return(aaa_eligibility_result(
      FALSE,
      paste0(
        "Differential abundance is disabled: ", length(sample_columns),
        " sample columns are selected but the declared design needs ",
        length(treatments) * replicate_n, " (", length(treatments),
        " treatments x ", replicate_n, " replicates)."
      ),
      guidance = c(
        "Make the number of selected sample columns equal to treatments x replicates.",
        "Either adjust the replicate design, the treatment names, or the sample selection.",
        "Remember that samples are assigned to treatments by position, in consecutive blocks."
      )
    ))
  }

  aaa_eligibility_result(
    TRUE,
    details = list(groups = length(treatments), replicates_per_group = replicate_n)
  )
}

aaa_check_rda_eligibility <- function(
  sample_columns,
  metadata,
  sample_id_column,
  environmental_variables
) {
  sample_columns <- as.character(sample_columns %||% character())
  environmental_variables <- as.character(environmental_variables %||% character())

  if (!length(sample_columns)) {
    return(aaa_eligibility_result(
      FALSE, "RDA is disabled: no valid abundance sample columns are selected.",
      guidance = c("Import an abundance table in 'Data and metadata'.")
    ))
  }
  if (is.null(metadata) || !is.data.frame(metadata) || !nrow(metadata)) {
    return(aaa_eligibility_result(
      FALSE, "RDA is disabled: no environmental metadata file is loaded.",
      guidance = c(
        "Upload an environmental metadata file in 'Data and metadata'.",
        "It needs one row per sample, a sample-identifier column, and one column per measured variable.",
        "RDA explains community composition using those measured gradients, so the file is mandatory."
      )
    ))
  }
  if (is.null(sample_id_column) || !nzchar(sample_id_column) || !sample_id_column %in% names(metadata)) {
    return(aaa_eligibility_result(
      FALSE, "RDA is disabled: no valid environmental sample-ID column is selected.",
      guidance = c(
        "Assign the 'identifier' role to the metadata column holding the sample names.",
        "Its values must match the abundance sample columns (differences in case, spaces, '-', '.' and '_' are tolerated)."
      )
    ))
  }

  environmental_variables <- intersect(environmental_variables, names(metadata))
  if (!length(environmental_variables)) {
    return(aaa_eligibility_result(
      FALSE, "RDA is disabled: no metadata column is classified as an environmental variable.",
      guidance = c(
        "Give the 'environmental_variable' role to at least one numeric metadata column.",
        "RDA only accepts numeric predictors: for categorical variables use envfit, partial RDA or variance partitioning."
      )
    ))
  }

  metadata_ids <- trimws(as.character(metadata[[sample_id_column]]))
  sample_matches <- tryCatch(
    aaa_match_sample_ids(sample_columns, metadata_ids),
    error = function(e) e
  )
  if (inherits(sample_matches, "error")) {
    return(aaa_eligibility_result(
      FALSE, paste0("RDA is disabled: ", conditionMessage(sample_matches)),
      guidance = c("Make the sample identifiers unique and non-empty in both the abundance table and the metadata.")
    ))
  }
  matched <- sample_matches$metadata_position
  if (any(!sample_matches$matched)) {
    missing_ids <- sample_matches$reference_id[!sample_matches$matched]
    return(aaa_eligibility_result(
      FALSE,
      paste0("RDA is disabled: metadata are missing for these abundance samples: ", paste(missing_ids, collapse = ", "), "."),
      guidance = c(
        "Add one metadata row per listed sample, or deselect those samples.",
        "Check for typographical differences between the two files."
      )
    ))
  }

  design <- metadata[matched, environmental_variables, drop = FALSE]
  design[] <- lapply(design, function(x) suppressWarnings(as.numeric(as.character(x))))
  complete <- stats::complete.cases(design)
  if (sum(complete) < 3L) {
    return(aaa_eligibility_result(
      FALSE,
      paste0(
        "RDA is disabled: only ", sum(complete),
        " matched samples have complete environmental values; at least three are required."
      ),
      guidance = c(
        "Fill in the missing environmental values, or deselect the variables that have gaps.",
        "Non-numeric entries (such as 'n.d.' or 'ND') count as missing: replace them with numbers or leave the cell empty.",
        "Samples with any missing value among the selected variables are dropped from the model."
      )
    ))
  }

  variable_ok <- vapply(design[complete, , drop = FALSE], function(x) {
    length(unique(x[is.finite(x)])) >= 2L && stats::var(x, na.rm = TRUE) > 0
  }, logical(1))
  if (!all(variable_ok)) {
    return(aaa_eligibility_result(
      FALSE,
      paste0("RDA is disabled: these environmental variables have no usable variation: ", paste(names(variable_ok)[!variable_ok], collapse = ", "), "."),
      guidance = c(
        "Remove the 'environmental_variable' role from the listed columns, or replace them with variables that actually vary.",
        "A predictor with the same value in every sample cannot explain any difference between samples."
      )
    ))
  }

  model_matrix <- tryCatch(
    stats::model.matrix(~., data = design[complete, , drop = FALSE]),
    error = function(e) NULL
  )
  if (is.null(model_matrix)) {
    return(aaa_eligibility_result(
      FALSE, "RDA is disabled: the environmental design matrix could not be constructed.",
      guidance = c("Check that every selected environmental column contains plain numeric values.")
    ))
  }
  design_rank <- qr(model_matrix)$rank
  if (design_rank >= sum(complete)) {
    return(aaa_eligibility_result(
      FALSE,
      paste0(
        "RDA is disabled: ", length(environmental_variables), " environmental variables against ",
        sum(complete), " complete samples leaves no residual degrees of freedom."
      ),
      guidance = c(
        paste0("Select at most ", max(1L, sum(complete) - 2L), " environmental variables for this number of samples."),
        "Alternatively, add more samples with complete metadata.",
        "A model with as many predictors as samples fits perfectly and cannot be tested."
      )
    ))
  }

  aaa_eligibility_result(
    TRUE,
    details = list(
      matched_samples = length(sample_columns),
      normalised_matches = sum(sample_matches$matched & !sample_matches$exact),
      complete_samples = sum(complete),
      variables = environmental_variables,
      design_rank = design_rank
    )
  )
}
