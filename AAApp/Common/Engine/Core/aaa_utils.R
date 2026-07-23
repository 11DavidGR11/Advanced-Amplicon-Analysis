#' Shared utilities for Advanced_Amplicon_Analysis
#'
#' Internal helpers used by all analysis modules.
#' @keywords internal
NULL

# `%||%` is defined once in aaa_globals.R.

aaa_n_replicates <- function(replicates) {
  replicates <- match.arg(
    replicates,
    c("none", "duplicate", "triplicate", "quadruplicate", "quintuplicate")
  )
  switch(replicates,
    none = 1L,
    duplicate = 2L,
    triplicate = 3L,
    quadruplicate = 4L,
    quintuplicate = 5L
  )
}

#' Groups present in a design, in their declared order
#'
#' The two ways of declaring a design produce different types: consecutive-block
#' mode yields a character vector with no inherent order beyond first
#' appearance, while a metadata-driven design yields a factor whose levels
#' encode the order declared in the metadata file.
#'
#' Levels with no sample are dropped. A metadata file commonly describes more
#' groups than the samples currently selected for analysis, and a group with no
#' sample is not part of the design: counting it makes `table()` report zero-size
#' groups, which turns the smallest group size into 0 and the group count into
#' the wrong number.
#'
#' @param treatment Treatment column of a sample design.
#' @return Character vector of group names, in declared order.
aaa_treatment_levels <- function(treatment) {
  if (is.factor(treatment)) {
    levels(droplevels(treatment))
  } else {
    unique(as.character(treatment))
  }
}

#' Number of samples per group, in declared order and without empty groups
#' @param treatment Treatment column of a sample design.
#' @return Named integer-valued table.
aaa_treatment_counts <- function(treatment) {
  table(factor(
    as.character(treatment),
    levels = aaa_treatment_levels(treatment)
  ))
}

#' Create the Advanced_Amplicon_Analysis project structure
#'
#' @param project_dir Root project directory.
#' @param analysis_name Analysis subdirectory.
#' @return Named list of paths.
aaa_results_root <- function(path) {
  if (!is.character(path) ||
    length(path) != 1L ||
    is.na(path) ||
    !nzchar(trimws(path))) {
    stop("'path' must be one non-empty character value.")
  }

  normalized <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )

  components <- strsplit(
    normalized,
    "/",
    fixed = TRUE
  )[[1]]

  runs_position <- which(
    tolower(components) == "runs"
  )

  if (length(runs_position) > 0) {
    position <- runs_position[length(runs_position)]

    if (position > 1L) {
      root <- paste(
        components[seq_len(position - 1L)],
        collapse = "/"
      )

      if (grepl("^[A-Za-z]:$", components[1])) {
        root <- paste0(
          components[1],
          "/",
          paste(
            components[2:(position - 1L)],
            collapse = "/"
          )
        )
      }

      return(root)
    }
  }

  if (tolower(basename(normalized)) == "results") {
    return(normalized)
  }

  file.path(
    normalized,
    "Results"
  )
}


aaa_create_project_structure <- function(
  project_dir,
  analysis_name = NULL,
  results_root = NULL
) {
  if (!is.character(project_dir) || length(project_dir) != 1L ||
    is.na(project_dir) || !nzchar(trimws(project_dir))) {
    stop("'project_dir' must be a single non-empty character string.")
  }

  if (!is.null(analysis_name)) {
    if (!is.character(analysis_name) || length(analysis_name) != 1L ||
      is.na(analysis_name) || !nzchar(trimws(analysis_name))) {
      stop("'analysis_name' must be NULL or one non-empty character string.")
    }
    if (grepl("[<>:\"/\\\\|?*]", analysis_name)) {
      stop("'analysis_name' contains invalid filename characters.")
    }
  }

  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

  if (is.null(results_root)) results_root <- aaa_results_root(project_dir)
  dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
  results_root <- normalizePath(results_root, winslash = "/", mustWork = TRUE)

  # storage contract: user outputs remain under the central Results tree,
  # while replaceable genome resources are shared by all projects in root/Cache.
  runs_dir <- file.path(results_root, "Runs")
  distribution_root <- getOption("triple_a_root", NULL)
  if (is.null(distribution_root) || !dir.exists(distribution_root)) {
    distribution_root <- dirname(results_root)
  }
  cache_dir <- file.path(distribution_root, "Cache")
  dir.create(runs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cache_dir, "GFF"), recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    project = project_dir,
    results_root = results_root,
    runs = runs_dir,
    cache = cache_dir,
    gff = file.path(cache_dir, "GFF"),
    reference_cache = cache_dir,
    # Metadata and logs are files within the run rather than empty top-level
    # directories.  Keeping this alias avoids touching stable writers.
    metadata = project_dir,
    analyses = project_dir
  )

  paths$analysis <- if (is.null(analysis_name)) {
    NULL
  } else {
    file.path(project_dir, analysis_name)
  }

  if (!is.null(paths$analysis)) {
    dir.create(paths$analysis, recursive = TRUE, showWarnings = FALSE)
  }

  # Cache subdirectories are created lazily only when data are written.
  paths
}


aaa_result_cache_status <- function(x, max_depth = 5L) {
  statuses <- character()
  visit <- function(object, depth = 0L) {
    if (is.null(object) || depth > max_depth) {
      return(invisible(NULL))
    }
    status <- attr(object, "triple_a_cache_status", exact = TRUE)
    if (!is.null(status)) statuses <<- c(statuses, toupper(as.character(status)))
    if (is.list(object) && !is.data.frame(object)) {
      nms <- names(object)
      if (!is.null(nms) && "metadata" %in% nms && is.list(object[["metadata"]])) {
        meta <- object[["metadata"]]
        if ("cache_status" %in% names(meta)) statuses <<- c(statuses, toupper(as.character(meta[["cache_status"]])))
      }
      for (item in object) visit(item, depth + 1L)
    }
    invisible(NULL)
  }
  visit(x)
  statuses <- statuses[statuses %in% c("HIT", "MISS")]
  if ("HIT" %in% statuses) "HIT" else "MISS"
}

# Collapsing cache activity to a single HIT/MISS badge hides how effective
# the cache actually was (a run with 1 lookup reused out of 500 shows the
# same badge as a fully cached run); sum the real per-analysis counts too.
aaa_result_cache_counts <- function(x, max_depth = 5L) {
  hits <- 0L
  misses <- 0L
  visit <- function(object, depth = 0L) {
    if (is.null(object) || depth > max_depth) {
      return(invisible(NULL))
    }
    if (is.list(object) && !is.data.frame(object)) {
      nms <- names(object)
      if (!is.null(nms) && "metadata" %in% nms && is.list(object[["metadata"]])) {
        meta <- object[["metadata"]]
        if (is.numeric(meta[["cache_hits"]])) hits <<- hits + meta[["cache_hits"]]
        if (is.numeric(meta[["cache_misses"]])) misses <<- misses + meta[["cache_misses"]]
      }
      for (item in object) visit(item, depth + 1L)
    }
    invisible(NULL)
  }
  visit(x)
  list(hits = hits, misses = misses)
}

# Run-state helpers ---------------------------------------------------------
aaa_run_status_file <- function(run_dir) file.path(run_dir, "run_status.json")

aaa_write_run_status <- function(run_dir, status, detail = NULL, extra = list()) {
  allowed <- c("running", "completed", "failed", "cancelled", "stopped")
  status <- tolower(as.character(status)[1L])
  if (!status %in% allowed) stop("Unsupported run status: ", status, call. = FALSE)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  record <- c(list(
    status = status,
    updated = format(Sys.time(), tz = "UTC", usetz = TRUE),
    detail = detail %||% NA_character_
  ), extra)
  jsonlite::write_json(record, aaa_run_status_file(run_dir),
    pretty = TRUE,
    auto_unbox = TRUE, null = "null"
  )
  invisible(record)
}

aaa_read_run_status <- function(run_dir) {
  file <- aaa_run_status_file(run_dir)
  if (!file.exists(file)) {
    return(list(status = "unknown", updated = NA_character_, detail = NA_character_))
  }
  tryCatch(jsonlite::read_json(file, simplifyVector = TRUE),
    error = function(e) {
      list(
        status = "unknown", updated = NA_character_,
        detail = conditionMessage(e)
      )
    }
  )
}

aaa_cache_status <- function(project_dir) {
  if (!is.character(project_dir) || length(project_dir) != 1L ||
    is.na(project_dir) || !nzchar(trimws(project_dir))) {
    stop("'project_dir' must be a single non-empty character string.")
  }

  project_dir <- normalizePath(
    project_dir,
    winslash = "/",
    mustWork = FALSE
  )
  candidate_root <- project_dir
  if (tolower(basename(candidate_root)) == "results") candidate_root <- dirname(candidate_root)
  cache_dir <- file.path(candidate_root, "Cache")

  cache_files <- if (dir.exists(cache_dir)) {
    list.files(
      cache_dir,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE
    )
  } else {
    character()
  }

  cache_files <- cache_files[file.exists(cache_files)]
  file_info <- if (length(cache_files) > 0L) file.info(cache_files) else NULL
  size_bytes <- if (is.null(file_info)) 0 else sum(file_info$size, na.rm = TRUE)
  last_update <- if (is.null(file_info) || nrow(file_info) == 0L) {
    "Not yet used"
  } else {
    format(max(file_info$mtime, na.rm = TRUE), "%Y-%m-%d %H:%M:%S")
  }

  cached_extensions <- c(
    "rds", "rda", "sqlite", "db", "json", "yaml", "yml"
  )
  cached_objects <- sum(
    tolower(tools::file_ext(cache_files)) %in% cached_extensions
  )

  data.frame(
    Metric = c(
      "Cache status",
      "Cache backend",
      "Reusable cached objects",
      "Cache size (MB)",
      "Last update"
    ),
    Value = c(
      if (dir.exists(cache_dir)) "Active" else "Ready (created on first use)",
      "SQLite + files",
      cached_objects,
      round(size_bytes / 1024^2, 2),
      last_update
    ),
    stringsAsFactors = FALSE
  )
}


aaa_acquire_cache_lock <- function(
  lock_file,
  timeout = 120,
  stale_after = 900
) {
  dir.create(
    dirname(lock_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  started <- Sys.time()

  repeat {
    if (file.exists(lock_file)) {
      age <- as.numeric(
        difftime(
          Sys.time(),
          file.info(lock_file)$mtime,
          units = "secs"
        )
      )

      if (is.finite(age) &&
        age > stale_after) {
        unlink(
          lock_file,
          force = TRUE
        )
      }
    }

    created <- file.create(
      lock_file,
      showWarnings = FALSE
    )

    if (isTRUE(created)) {
      writeLines(
        c(
          paste0("pid=", Sys.getpid()),
          paste0(
            "created=",
            format(
              Sys.time(),
              "%Y-%m-%d %H:%M:%S"
            )
          )
        ),
        lock_file
      )

      return(invisible(TRUE))
    }

    elapsed <- as.numeric(
      difftime(
        Sys.time(),
        started,
        units = "secs"
      )
    )

    if (elapsed >= timeout) {
      stop(
        "Timed out while waiting for cache lock: ",
        lock_file
      )
    }

    Sys.sleep(0.2)
  }
}


aaa_release_cache_lock <- function(lock_file) {
  if (file.exists(lock_file)) {
    unlink(
      lock_file,
      force = TRUE
    )
  }

  invisible(TRUE)
}


aaa_write_session_information <- function(path) {
  metadata <- aaa_project_metadata()

  information <- capture.output({
    cat("Project:", metadata$title, "\n")
    cat("Short name:", metadata$subtitle, "\n")
    cat("Triple_A release identifier:", metadata$version, "\n")
    cat("Author(s):", aaa_author_display(" ; "), "\n")
    cat("License:", metadata$license, "\n")

    if (nzchar(metadata$repository)) {
      cat("Repository:", metadata$repository, "\n")
    }

    if (nzchar(metadata$doi)) {
      cat("DOI:", metadata$doi, "\n")
    }

    cat("Suggested citation:", metadata$citation, "\n")
    cat("Generated:", format(Sys.time()), "\n\n")
    print(sessionInfo())
  })

  writeLines(
    information,
    con = path,
    useBytes = TRUE
  )

  invisible(path)
}

aaa_append_log <- function(path, stage, detail) {
  line <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    stage,
    " | ",
    detail
  )

  directory <- dirname(path)
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  written <- tryCatch(
    {
      cat(
        line,
        "\n",
        file = path,
        append = TRUE,
        sep = ""
      )
      TRUE
    },
    error = function(e) {
      warning(
        "Triple_A could not write to the execution log: ",
        conditionMessage(e),
        call. = FALSE
      )
      FALSE
    }
  )

  invisible(structure(
    line,
    written = written,
    path = path
  ))
}


#' Parse taxonomy strings from common amplicon/metagenomic formats
#'
#' Supports prefixed lineages (d__/k__/p__/.../s__), unprefixed lineages,
#' and separators commonly produced by QIIME 2, SILVA, GTDB, Kraken/Bracken,
#' MetaPhlAn and spreadsheet exports.
aaa_get_tax_info <- function(x) {
  # This function runs once per taxon inside a lapply()+bind_rows() call in
  # aaa_prepare_amplicon_data(), i.e. potentially thousands of times per
  # analysis. Returning plain list()s instead of tibble::tibble() avoids
  # tibble's per-call quasiquotation/glue overhead, which profiling showed
  # dominates aaa_prepare_amplicon_data()'s runtime. dplyr::bind_rows() binds
  # a list of lists identically to a list of one-row tibbles.
  empty <- function() {
    list(
      Taxonomy = NA_character_, Genus = NA_character_, Tax_level = "Unknown"
    )
  }
  if (length(x) == 0L || is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(empty())
  }

  value <- trimws(as.character(x))
  value <- gsub("\\|", ";", value)
  value <- gsub(",(?=\\s*(?:[dkpcofgs](?:__|_|:)|[A-Za-z]+[:=]))", ";", value, perl = TRUE)
  taxa <- trimws(strsplit(value, ";", fixed = TRUE)[[1]])
  taxa <- taxa[nzchar(taxa)]
  taxa <- taxa[!tolower(taxa) %in% c("na", "none", "unknown", "unclassified", "uncultured")]
  if (!length(taxa)) {
    return(empty())
  }

  rank_letter <- function(z) {
    z <- trimws(z)
    if (grepl("^[dkpcofgs](?:__|_|:)", z, ignore.case = TRUE, perl = TRUE)) {
      return(tolower(substr(z, 1, 1)))
    }
    key <- tolower(sub("[:=].*$", "", z))
    map <- c(
      domain = "d", superkingdom = "d", kingdom = "k", phylum = "p",
      class = "c", order = "o", family = "f", genus = "g", species = "s"
    )
    if (key %in% names(map)) {
      return(unname(map[[key]]))
    }
    NA_character_
  }

  clean_rank <- function(z) {
    z <- trimws(z)
    z <- sub("^[A-Za-z]+[:=]\\s*", "", z)
    z <- sub("^[dkpcofgs](?:__|_|:)\\s*", "", z, ignore.case = TRUE, perl = TRUE)
    trimws(z)
  }

  cleaned <- vapply(taxa, clean_rank, character(1))
  invalid_values <- c("na", "none", "unknown", "unclassified", "uncultured", "")
  valid <- !tolower(cleaned) %in% invalid_values
  taxa <- taxa[valid]
  cleaned <- cleaned[valid]
  if (!length(cleaned)) {
    return(empty())
  }

  letters <- vapply(taxa, rank_letter, character(1))
  level_map <- c(
    d = "Domain", k = "Kingdom", p = "Phylum", c = "Class",
    o = "Order", f = "Family", g = "Genus", s = "Species"
  )

  # Use the deepest explicitly labelled genus/species entry when available.
  usable_idx <- which(letters %in% c("g", "s"))
  if (length(usable_idx)) {
    lowest_index <- max(usable_idx)
    lowest_letter <- letters[lowest_index]
    taxonomy <- cleaned[lowest_index]
    genus_idx <- which(letters == "g" & seq_along(letters) <= lowest_index)
    genus <- if (length(genus_idx)) cleaned[max(genus_idx)] else taxonomy
    level <- unname(level_map[[lowest_letter]])
    return(list(Taxonomy = taxonomy, Genus = genus, Tax_level = level))
  }

  # Infer common unprefixed lineage layouts. Six ranks generally end at genus;
  # seven or more generally end at species with genus immediately before it.
  if (all(is.na(letters))) {
    if (length(cleaned) >= 7L) {
      return(list(
        Taxonomy = cleaned[length(cleaned)],
        Genus = cleaned[length(cleaned) - 1L],
        Tax_level = "Species"
      ))
    }
    if (length(cleaned) == 6L) {
      return(list(
        Taxonomy = cleaned[length(cleaned)],
        Genus = cleaned[length(cleaned)],
        Tax_level = "Genus"
      ))
    }
    # A single unprefixed value is commonly a genus/species name selected by
    # the user as the taxonomic column. Preserve it rather than discarding the
    # whole dataset. Ambiguous multi-rank lineages shorter than six are not
    # promoted above genus.
    if (length(cleaned) == 1L) {
      return(list(
        Taxonomy = cleaned[1], Genus = cleaned[1], Tax_level = "Genus"
      ))
    }
  }

  empty()
}

# In-memory memoization for the expensive taxonomy-parsing/pivot step below.
# A single workflow run typically calls aaa_prepare_amplicon_data() once per
# selected analysis (functional_potential, top_abundance, differential_abundance,
# functional_abundance, community_structure, plsda, rda) with the SAME dataset
# and abundance_type, differing only in analysis_name/project_dir (which only
# affect the cheap aaa_create_project_structure() side effect) and filter_genus
# (which has exactly two possible values in practice). Caching wide/long by a
# content hash of the inputs that actually affect them avoids redundant taxonomy
# parsing across analyses in the same run. Bounded to a handful of entries so a
# long-running session moving between projects/datasets does not grow unbounded.
.triple_a_prepared_cache <- new.env(parent = emptyenv())
.triple_a_prepared_cache$order <- character()
.triple_a_prepared_cache$max_entries <- 4L

aaa_prepared_cache_get <- function(key) {
  if (is.na(key) || !exists(key, envir = .triple_a_prepared_cache, inherits = FALSE)) {
    return(NULL)
  }
  get(key, envir = .triple_a_prepared_cache, inherits = FALSE)
}

aaa_prepared_cache_set <- function(key, value) {
  if (is.na(key)) {
    return(invisible(NULL))
  }
  assign(key, value, envir = .triple_a_prepared_cache)
  .triple_a_prepared_cache$order <- c(setdiff(.triple_a_prepared_cache$order, key), key)
  while (length(.triple_a_prepared_cache$order) > .triple_a_prepared_cache$max_entries) {
    stale_key <- .triple_a_prepared_cache$order[1]
    if (exists(stale_key, envir = .triple_a_prepared_cache, inherits = FALSE)) {
      rm(list = stale_key, envir = .triple_a_prepared_cache)
    }
    .triple_a_prepared_cache$order <- .triple_a_prepared_cache$order[-1]
  }
  invisible(NULL)
}

aaa_prepare_amplicon_data <- function(
  dataset,
  abundance_type = c("proportion", "percentage", "counts"),
  project_dir, analysis_name, filter_genus = TRUE
) {
  abundance_type <- match.arg(abundance_type)
  aaa_validate_dataset(dataset)
  project <- aaa_create_project_structure(project_dir, analysis_name)

  raw <- dataset$abundance
  sample_map <- dataset$sample_design
  sample_columns <- as.character(sample_map$Sample_column)

  abundance_validation <- aaa_validate_abundance_nature(
    data = raw,
    sample_columns = sample_columns,
    abundance_type = abundance_type
  )
  if (!isTRUE(abundance_validation$available)) {
    stop(abundance_validation$reason, call. = FALSE)
  }

  samples_name <- aaa_treatment_levels(sample_map$Treatment)
  # Treatments are identified by label, not by column position, so groups may
  # have different numbers of replicates. `n_replicates` is the SMALLEST group
  # size, because that is what the analyses that require replication need to
  # check; `has_replicates` says whether averaging into group means is
  # meaningful at all.
  replicate_counts <- aaa_treatment_counts(sample_map$Treatment)
  n_replicates <- as.integer(min(replicate_counts))
  has_replicates <- any(replicate_counts > 1L)
  balanced <- length(unique(as.integer(replicate_counts))) == 1L
  replicates <- if (!has_replicates) {
    "none"
  } else if (balanced) {
    switch(as.character(n_replicates),
      `2` = "duplicate",
      `3` = "triplicate",
      `4` = "quadruplicate",
      `5` = "quintuplicate",
      paste0(n_replicates, " replicates per treatment")
    )
  } else {
    paste0(
      "unbalanced (",
      paste(sprintf("%s=%d", names(replicate_counts), as.integer(replicate_counts)),
        collapse = ", "
      ),
      ")"
    )
  }

  cache_key <- tryCatch(
    aaa_hash_object(list(
      raw = raw[, c("Taxonomy", sample_columns), drop = FALSE],
      sample_map = sample_map,
      abundance_type = abundance_type,
      filter_genus = filter_genus
    )),
    error = function(e) NA_character_
  )

  cached <- aaa_prepared_cache_get(cache_key)

  # Library sizes are only knowable while the table still holds raw counts.
  # They are captured here and carried forward because `wide` is rescaled to
  # percentages below, and analyses that genuinely need counts (ANCOM-BC2)
  # cannot reconstruct sequencing depth from a percentage table.
  library_sizes <- NULL

  if (!is.null(cached)) {
    wide <- cached$wide
    long <- cached$long
    library_sizes <- cached$library_sizes
  } else {
    tax_info <- dplyr::bind_rows(lapply(raw$Taxonomy, aaa_get_tax_info))
    wide <- dplyr::bind_cols(tax_info, raw[, sample_columns, drop = FALSE]) |>
      dplyr::group_by(Taxonomy, Genus, Tax_level) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(sample_columns), \(x) sum(x, na.rm = TRUE)),
        .groups = "drop"
      )

    if (filter_genus) {
      wide <- wide |> dplyr::filter(!is.na(Genus), nzchar(Genus), Genus != "None")
    }
    if (nrow(wide) == 0L) {
      examples <- unique(utils::head(as.character(raw$Taxonomy), 5L))
      stop(paste0(
        "No usable genus- or species-level taxa remained after filtering. ",
        "Check the selected taxonomy column/rank mapping. Example imported values: ",
        paste(examples, collapse = " | ")
      ))
    }

    if (abundance_type == "proportion") {
      wide <- wide |> dplyr::mutate(dplyr::across(dplyr::all_of(sample_columns), \(x) x * 100))
    } else if (abundance_type == "counts") {
      totals <- colSums(wide[sample_columns], na.rm = TRUE)
      if (any(totals <= 0)) stop("Count data contain sample columns with a total of zero.")
      library_sizes <- totals
      wide[sample_columns] <- sweep(wide[sample_columns], 2, totals, "/") * 100
    }

    long <- wide |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(sample_columns),
        names_to = "Sample_column", values_to = "Abundance"
      ) |>
      dplyr::left_join(sample_map, by = "Sample_column")

    aaa_prepared_cache_set(
      cache_key,
      list(wide = wide, long = long, library_sizes = library_sizes)
    )
  }

  list(
    wide = wide, long = long, raw = raw,
    sample_columns = sample_columns, sample_map = sample_map,
    samples_name = samples_name, replicates = replicates,
    n_replicates = n_replicates, has_replicates = has_replicates,
    balanced = balanced, replicate_counts = replicate_counts,
    abundance_type = abundance_type,
    library_sizes = library_sizes,
    metadata = dataset$metadata, metadata_roles = dataset$metadata_roles,
    project = project
  )
}

#' Calculate means, SDs and display labels for replicate groups
#' Calculate means, SDs and display labels for replicate groups
#'
#' Columns are assigned to treatments by looking each column name up in
#' `sample_map`, never by their position in the table. The previous positional
#' version sliced the columns into consecutive blocks of equal size, which
#' silently mislabelled interleaved designs (C1, T1, C2, T2...) and made equal
#' group sizes a hard requirement of the whole application.
#'
#' @param data_frame Taxa x samples table. Column names must be sample IDs.
#' @param sample_map Data frame with `Sample_column` and `Treatment`.
#' @param group_levels Optional treatment order for the output columns.
aaa_replicate_summary <- function(data_frame, sample_map, group_levels = NULL) {
  sample_map <- as.data.frame(sample_map, stringsAsFactors = FALSE)
  if (!all(c("Sample_column", "Treatment") %in% names(sample_map))) {
    stop("'sample_map' must contain Sample_column and Treatment.")
  }

  columns <- names(data_frame)
  positions <- match(columns, as.character(sample_map$Sample_column))
  if (anyNA(positions)) {
    stop(
      "These sample columns are absent from the sample design: ",
      paste(columns[is.na(positions)], collapse = ", ")
    )
  }
  treatments <- as.character(sample_map$Treatment)[positions]

  if (is.null(group_levels)) group_levels <- unique(treatments)
  group_levels <- as.character(group_levels)
  group_levels <- group_levels[group_levels %in% treatments]
  if (!length(group_levels)) {
    stop("No treatment group matches the supplied sample columns.")
  }

  mean_data <- data.frame(row.names = rownames(data_frame))
  sd_data <- data.frame(row.names = rownames(data_frame))

  for (name in group_levels) {
    cols <- columns[treatments == name]
    mean_data[[name]] <- rowMeans(data_frame[, cols, drop = FALSE],
      na.rm = TRUE
    )
    # A single-sample group has no dispersion to report; apply() would return
    # NA anyway, but stating it keeps the label formatting predictable.
    sd_data[[name]] <- if (length(cols) < 2L) {
      rep(NA_real_, nrow(data_frame))
    } else {
      apply(data_frame[, cols, drop = FALSE], 1, stats::sd, na.rm = TRUE)
    }
  }

  # Groups with a single sample print the mean alone instead of "x ± NA%".
  mean_values <- as.matrix(mean_data)
  sd_values <- as.matrix(sd_data)
  labels <- matrix(
    ifelse(
      is.finite(sd_values),
      sprintf("%.1f%% ± %.1f%%", mean_values, sd_values),
      sprintf("%.1f%%", mean_values)
    ),
    nrow = nrow(mean_data), dimnames = dimnames(mean_data)
  ) |>
    as.data.frame(stringsAsFactors = FALSE)

  list(mean = mean_data, sd = sd_data, labels = labels)
}

aaa_safe_name <- function(x) {
  gsub("[^A-Za-z0-9_-]", "_", trimws(x))
}

aaa_result <- function(tables = list(), plots = list(), files = character(),
                       output_dir, metadata = list()) {
  if (is.null(tables$summary)) {
    tables <- c(list(summary = data.frame(
      Metric = c("Output directory", "Tables generated", "Figures generated", "Files generated"),
      Value = c(
        normalizePath(output_dir, winslash = "/", mustWork = FALSE),
        length(tables), length(plots), length(files)
      ),
      stringsAsFactors = FALSE
    )), tables)
  }
  structure(
    list(
      tables = tables,
      plots = plots,
      files = files,
      output_dir = output_dir,
      metadata = metadata
    ),
    class = c("Triple_A_result", "list")
  )
}

print.Triple_A_result <- function(x, ...) {
  cat("<Advanced_Amplicon_Analysis result>\n")
  cat("Output directory:", x$output_dir, "\n")
  cat("Tables:", paste(names(x$tables), collapse = ", "), "\n")
  cat("Plots:", paste(names(x$plots), collapse = ", "), "\n")
  invisible(x)
}


aaa_collect_files <- function(x, prefix = NULL) {
  collected <- list()

  walk <- function(value, path) {
    if (is.null(value)) {
      return(invisible(NULL))
    }

    if (is.character(value) &&
      length(value) >= 1 &&
      all(file.exists(value) | is.na(value))) {
      valid <- value[!is.na(value)]
      for (i in seq_along(valid)) {
        name <- paste(
          c(path, names(valid)[i] %||% i),
          collapse = " / "
        )
        collected[[name]] <<- valid[i]
      }
      return(invisible(NULL))
    }

    if (is.list(value)) {
      child_names <- names(value)
      if (is.null(child_names)) {
        child_names <- as.character(seq_along(value))
      }

      for (i in seq_along(value)) {
        walk(
          value[[i]],
          c(path, child_names[i])
        )
      }
    }

    invisible(NULL)
  }

  if (is.list(x)) {
    walk(x, prefix %||% character())
  }

  collected
}

aaa_collect_tables <- function(x, prefix = NULL) {
  collected <- list()

  walk <- function(value, path) {
    if (is.null(value)) {
      return(invisible(NULL))
    }

    if (is.data.frame(value) ||
      is.matrix(value)) {
      name <- paste(path, collapse = " / ")
      collected[[name]] <<- value
      return(invisible(NULL))
    }

    if (is.list(value)) {
      child_names <- names(value)
      if (is.null(child_names)) {
        child_names <- as.character(seq_along(value))
      }

      for (i in seq_along(value)) {
        walk(
          value[[i]],
          c(path, child_names[i])
        )
      }
    }

    invisible(NULL)
  }

  walk(x, prefix %||% character())
  collected
}
