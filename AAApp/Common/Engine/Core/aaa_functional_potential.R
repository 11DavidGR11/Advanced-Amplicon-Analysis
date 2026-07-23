# Return TRUE only when a cached taxonomy record contains enough information
# to evaluate reference-genome markers. Failed or incomplete lookups are not
# permanent cache hits and must be retried on a later run.
aaa_reference_cache_record_is_reusable <- function(cache, index) {
  if (is.null(cache) || !is.data.frame(cache) || length(index) != 1L ||
    is.na(index) || index < 1L || index > nrow(cache)) {
    return(FALSE)
  }

  required <- c("TaxID", "Reference_genome")
  if (!all(required %in% names(cache))) {
    return(FALSE)
  }

  taxid <- as.character(cache[["TaxID"]][index])
  genome <- as.character(cache[["Reference_genome"]][index])

  length(taxid) == 1L && !is.na(taxid) && nzchar(trimws(taxid)) &&
    length(genome) == 1L && !is.na(genome) && nzchar(trimws(genome))
}


# Insert or update one taxonomy/reference row without creating duplicate
# taxonomy records. Existing marker columns are preserved when possible.
aaa_upsert_reference_cache_record <- function(
  cache, taxonomy, genus, tax_level, taxid, reference_genome
) {
  new_values <- list(
    TaxID = taxid,
    Taxonomy = taxonomy,
    Genus = genus,
    Tax_level = tax_level,
    Reference_genome = reference_genome
  )

  if (is.null(cache) || !is.data.frame(cache) || nrow(cache) == 0L) {
    return(tibble::as_tibble(new_values))
  }

  missing_columns <- setdiff(names(new_values), names(cache))
  for (column in missing_columns) cache[[column]] <- NA_character_

  index <- match(taxonomy, cache$Taxonomy)
  if (is.na(index)) {
    row <- as.list(rep(NA, ncol(cache)))
    names(row) <- names(cache)
    for (column in names(new_values)) row[[column]] <- new_values[[column]]
    return(dplyr::bind_rows(cache, tibble::as_tibble(row)))
  }

  for (column in names(new_values)) cache[[column]][index] <- new_values[[column]]
  cache
}

#' Infer functional potential from amplicon-derived taxonomy
#'
#' Links taxonomic assignments to NCBI Taxonomy and RefSeq assemblies, searches
#' GFF annotations for selected genes and applies a user-supplied classifier.
#'
#' @return A Triple_A_result object.
aaa_functional_potential <- function(
  dataset, genes, graph_main,
  classification_function, graph_note = NULL,
  abundance_type = c("proportion", "percentage", "counts"),
  project_dir, analysis_name,
  biological_function_id = NULL,
  verbose = TRUE,
  progress_callback = NULL,
  progress_verbosity = c(
    "standard", "detailed", "developer"
  )
) {
  abundance_type <- match.arg(abundance_type)
  progress_verbosity <- match.arg(progress_verbosity)
  aaa_check_packages(analyses = "functional_potential")

  verbosity_rank <- c(
    standard = 1L,
    detailed = 2L,
    developer = 3L
  )

  report_progress <- function(stage,
                              detail,
                              completed,
                              total,
                              level = "standard") {
    if (
      verbosity_rank[[level]] <=
        verbosity_rank[[progress_verbosity]] &&
        is.function(progress_callback)
    ) {
      progress_callback(
        stage = stage,
        detail = detail,
        completed = completed,
        total = total
      )
    }

    invisible(NULL)
  }

  report_progress(
    "input",
    "Validating the functional-analysis parameters.",
    0,
    1
  )

  if (!is.character(genes) || length(genes) == 0 ||
    any(is.na(genes)) || any(!nzchar(trimws(genes)))) {
    stop("'genes' must be a non-empty character vector.")
  }
  genes <- trimws(genes)
  if (anyDuplicated(genes)) stop("'genes' contains duplicated entries.")
  if (!is.function(classification_function)) {
    stop("'classification_function' must be a function.")
  }
  test <- tryCatch(
    classification_function(stats::setNames(rep(FALSE, length(genes)), genes)),
    error = function(e) NULL
  )
  if (!is.character(test) || length(test) != 1 || is.na(test) ||
    !nzchar(trimws(test))) {
    stop("'classification_function' must return one non-empty character value.")
  }

  report_progress(
    "input",
    "Reading and preparing the amplicon abundance table.",
    0.15,
    1
  )

  prepared <- aaa_prepare_amplicon_data(
    dataset = dataset,
    abundance_type = abundance_type,
    project_dir = project_dir,
    analysis_name = analysis_name
  )
  project <- prepared$project
  taxa <- prepared$wide |>
    dplyr::select(Taxonomy, Genus, Tax_level) |>
    dplyr::distinct()

  report_progress(
    "cache",
    "Loading the taxonomic and reference-annotation cache.",
    0.25,
    1
  )

  cache_info <- aaa_load_reference_cache(project)
  cache <- cache_info$cache
  cache_file <- cache_info$file
  cache_updated <- FALSE
  cache_hits <- 0L
  cache_misses <- 0L

  # Hash-based Taxonomy -> row index lookup. A plain match(x, cache$Taxonomy)
  # inside the per-taxon loop below rebuilds a hash of the whole cache on
  # every call (O(n_taxa * n_cache)). Building the hash once up front and
  # keeping it in sync as new rows are appended keeps lookups O(1) while
  # preserving the original behaviour of reusing rows inserted earlier in
  # the same run.
  taxonomy_lookup <- new.env(parent = emptyenv())
  if (nrow(cache) > 0L) {
    for (row_index in seq_len(nrow(cache))) {
      taxonomy_lookup[[cache$Taxonomy[row_index]]] <- row_index
    }
  }

  report_progress(
    "cache",
    paste0(
      "Shared SQLite cache loaded: ",
      nrow(cache),
      " reference record(s) and ",
      length(list.files(
        project$gff,
        pattern = "\\.gff\\.gz$",
        full.names = TRUE
      )),
      " cached GFF annotation(s)."
    ),
    0.30,
    1
  )

  taxa$TaxID <- NA_character_
  taxa$Reference_genome <- NA_character_

  n_taxa <- nrow(taxa)

  for (i in seq_len(n_taxa)) {
    report_progress(
      "functional_monitor",
      paste(
        "reference",
        i,
        n_taxa,
        taxa$Taxonomy[i],
        analysis_name,
        0,
        length(genes),
        sep = "|||"
      ),
      i - 1,
      max(1, n_taxa),
      level = "standard"
    )

    report_progress(
      "ncbi_reference",
      paste0(
        "Resolving taxon ",
        i,
        " of ",
        n_taxa,
        ": ",
        taxa$Taxonomy[i]
      ),
      i - 1,
      max(1, n_taxa),
      level = "detailed"
    )

    idx <- taxonomy_lookup[[taxa$Taxonomy[i]]]
    if (is.null(idx)) idx <- NA_integer_
    reusable_reference <- aaa_reference_cache_record_is_reusable(cache, idx)

    if (isTRUE(reusable_reference)) {
      taxa$TaxID[i] <- as.character(cache$TaxID[idx])
      taxa$Reference_genome[i] <- as.character(cache$Reference_genome[idx])
      cache_hits <- cache_hits + 1L

      report_progress(
        "ncbi_reference",
        paste0(
          "Reused cached NCBI reference for ",
          taxa$Taxonomy[i],
          "."
        ),
        i,
        max(1, n_taxa),
        level = "detailed"
      )
    } else {
      # Incomplete rows commonly result from temporary NCBI/network failures.
      # They are deliberately retried instead of being treated as permanent
      # cache hits, which previously made every marker remain NA/Unknown.
      cache_misses <- cache_misses + 1L
      query_candidates <- aaa_functional_taxon_candidates(
        taxonomy = taxa$Taxonomy[i],
        genus = taxa$Genus[i],
        tax_level = taxa$Tax_level[i]
      )

      taxid <- NA_character_
      genome <- NA_character_
      resolved_query <- NA_character_
      resolved_rank <- NA_character_

      if (nrow(query_candidates) > 0L) {
        for (candidate_index in seq_len(nrow(query_candidates))) {
          candidate_query <- query_candidates$Query[candidate_index]
          candidate_rank <- query_candidates$Rank[candidate_index]
          candidate_taxid <- aaa_get_taxid(candidate_query, verbose = verbose)

          if (is.na(candidate_taxid) || !nzchar(trimws(candidate_taxid))) {
            next
          }

          candidate_genome <- aaa_get_reference_genome(
            candidate_taxid,
            verbose = verbose
          )

          if (is.na(candidate_genome) || !nzchar(trimws(candidate_genome))) {
            next
          }

          taxid <- candidate_taxid
          genome <- candidate_genome
          resolved_query <- candidate_query
          resolved_rank <- candidate_rank
          break
        }
      }

      taxa$TaxID[i] <- taxid
      taxa$Reference_genome[i] <- genome
      if (!"Reference_query" %in% names(taxa)) taxa$Reference_query <- NA_character_
      if (!"Reference_rank" %in% names(taxa)) taxa$Reference_rank <- NA_character_
      taxa$Reference_query[i] <- resolved_query
      taxa$Reference_rank[i] <- resolved_rank
      cache <- aaa_upsert_reference_cache_record(
        cache = cache,
        taxonomy = taxa$Taxonomy[i],
        genus = taxa$Genus[i],
        tax_level = taxa$Tax_level[i],
        taxid = taxid,
        reference_genome = genome
      )
      cache_updated <- TRUE
      # Keep the hash lookup in sync so a taxonomy string that repeats
      # later in this same run (e.g. same Taxonomy under a different
      # Genus/Tax_level) is treated as a cache hit instead of re-querying NCBI.
      taxonomy_lookup[[taxa$Taxonomy[i]]] <- nrow(cache)

      report_progress(
        "ncbi_reference",
        paste0(
          "Resolved NCBI TaxID and reference genome for ",
          taxa$Taxonomy[i],
          "."
        ),
        i,
        max(1, n_taxa),
        level = "detailed"
      )
    }
  }

  report_progress(
    "gene_search",
    paste0(
      "Preparing ",
      length(genes),
      " curated gene markers."
    ),
    0,
    1
  )

  for (g in genes) {
    taxa[[g]] <- NA
    if (!g %in% names(cache)) cache[[g]] <- NA
  }

  valid_taxid <- !is.na(taxa$TaxID) & nzchar(taxa$TaxID)
  cache_idx <- rep(NA_integer_, nrow(taxa))
  cache_idx[valid_taxid] <- match(taxa$TaxID[valid_taxid], cache$TaxID)
  gff_memory <- list()

  for (i in seq_len(n_taxa)) {
    report_progress(
      "functional_monitor",
      paste(
        "genome",
        i,
        n_taxa,
        taxa$Taxonomy[i],
        analysis_name,
        0,
        length(genes),
        sep = "|||"
      ),
      i - 1,
      max(1, n_taxa),
      level = "standard"
    )

    report_progress(
      "gff",
      paste0(
        "Checking GFF annotations for taxon ",
        i,
        " of ",
        n_taxa,
        ": ",
        taxa$Taxonomy[i]
      ),
      i - 1,
      max(1, n_taxa),
      level = "detailed"
    )

    idx <- cache_idx[i]
    if (is.na(idx)) next

    missing_genes <- genes[
      is.na(unlist(cache[idx, genes, drop = FALSE], use.names = FALSE))
    ]

    cached_gene_count <- length(genes) -
      length(missing_genes)
    cache_hits <- cache_hits + cached_gene_count
    cache_misses <- cache_misses + length(missing_genes)

    report_progress(
      "functional_monitor",
      paste(
        "genes",
        i,
        n_taxa,
        taxa$Taxonomy[i],
        analysis_name,
        cached_gene_count,
        length(genes),
        sep = "|||"
      ),
      i - 1,
      max(1, n_taxa),
      level = "standard"
    )

    if (length(missing_genes) > 0 &&
      !is.na(taxa$Reference_genome[i]) &&
      nzchar(taxa$Reference_genome[i])) {
      accession <- taxa$Reference_genome[i]
      if (!accession %in% names(gff_memory)) {
        report_progress(
          "gff",
          paste0(
            "Checking GFF annotation ",
            accession,
            " for ",
            taxa$Taxonomy[i],
            "."
          ),
          i - 0.75,
          max(1, n_taxa),
          level = "detailed"
        )

        gff_memory[[accession]] <- aaa_get_gff(
          accession,
          project,
          verbose = verbose
        )

        gff_status <- attr(
          gff_memory[[accession]],
          "triple_a_cache_status"
        )

        report_progress(
          "gff",
          if (identical(gff_status, "cache_hit")) {
            paste0("Using cached GFF annotation ", accession, ".")
          } else if (identical(gff_status, "downloaded")) {
            paste0("Downloaded and cached GFF annotation ", accession, ".")
          } else {
            paste0("GFF annotation could not be loaded for ", accession, ".")
          },
          i - 0.50,
          max(1, n_taxa),
          level = "standard"
        )
      }
      gff <- gff_memory[[accession]]
      if (!is.null(gff)) {
        for (gene_index in seq_along(missing_genes)) {
          g <- missing_genes[gene_index]

          report_progress(
            "gene_search",
            paste0(
              "Searching ",
              g,
              " in ",
              accession,
              " (",
              gene_index,
              "/",
              length(missing_genes),
              ")."
            ),
            gene_index - 1,
            max(1, length(missing_genes)),
            level = "developer"
          )

          cache[idx, g] <- aaa_search_gene(
            gff,
            g,
            verbose = verbose
          )
          cache_updated <- TRUE

          report_progress(
            "functional_monitor",
            paste(
              "genes",
              i,
              n_taxa,
              taxa$Taxonomy[i],
              analysis_name,
              cached_gene_count + gene_index,
              length(genes),
              sep = "|||"
            ),
            i - 1,
            max(1, n_taxa),
            level = "standard"
          )
        }
      }
    }

    for (g in genes) {
      taxa[[g]][i] <- as.logical(cache[idx, g])
    }

    report_progress(
      "functional_monitor",
      paste(
        "complete",
        i,
        n_taxa,
        taxa$Taxonomy[i],
        analysis_name,
        length(genes),
        length(genes),
        sep = "|||"
      ),
      i,
      max(1, n_taxa),
      level = "standard"
    )

    report_progress(
      "gff",
      paste0(
        "Completed GFF and gene evaluation for ",
        taxa$Taxonomy[i],
        "."
      ),
      i,
      max(1, n_taxa),
      level = "detailed"
    )
  }

  if (cache_updated) {
    report_progress(
      "cache",
      "Updating the reference-genome and gene cache.",
      0.70,
      1
    )
    cache <- aaa_update_reference_cache(cache, cache_file)
  }

  report_progress(
    "classification",
    "Applying the curated biological classification rules.",
    0.80,
    1
  )

  taxa$Potential <- vapply(
    seq_len(nrow(taxa)),
    function(i) {
      x <- vapply(genes, function(g) {
        value <- taxa[[g]][i]
        if (length(value) == 0L || is.na(value)) NA else isTRUE(as.logical(value))
      }, logical(1))
      names(x) <- genes
      classification_function(x)
    },
    character(1)
  )

  evaluated_markers <- rowSums(!is.na(taxa[, genes, drop = FALSE]))
  detected_markers <- rowSums(as.data.frame(lapply(taxa[, genes, drop = FALSE], function(z) z %in% TRUE)))
  taxa$Markers_evaluated <- evaluated_markers
  taxa$Markers_detected <- detected_markers
  if (all(evaluated_markers == 0L)) {
    warning(
      "No functional markers could be evaluated. Triple_A found no reusable ",
      "reference/GFF evidence for the supplied taxonomy. Check internet access, ",
      "taxonomy resolution and legacy cache migration; classifications are not ",
      "interpretable until at least one marker is evaluated.",
      call. = FALSE
    )
  }

  evidence_table <- NULL

  if (!is.null(biological_function_id)) {
    evidence_table <- aaa_summarize_function_evidence(
      taxa,
      biological_function_id
    )
    taxa <- evidence_table
  }

  report_progress(
    "abundance_summary",
    "Combining predicted functional potential with taxonomic abundances.",
    0.87,
    1
  )

  joined <- prepared$wide |>
    dplyr::left_join(taxa, by = c("Taxonomy", "Genus", "Tax_level"))

  functional_abundance <- joined |>
    dplyr::group_by(Potential) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(prepared$sample_columns),
        \(x) sum(x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  report_progress(
    "figures",
    "Preparing the functional-potential heatmap and summary tables.",
    0.93,
    1
  )

  heatmap_matrix <- functional_abundance |>
    tibble::column_to_rownames("Potential") |>
    as.data.frame()

  if (isTRUE(prepared$has_replicates)) {
    summary <- aaa_replicate_summary(
      heatmap_matrix, prepared$sample_map, prepared$samples_name
    )
    plot_matrix <- summary$mean
    plot_labels <- summary$labels
  } else {
    plot_matrix <- heatmap_matrix
    names(plot_matrix) <- prepared$samples_name
    plot_labels <- matrix(
      sprintf("%.1f%%", as.matrix(plot_matrix)),
      nrow = nrow(plot_matrix), dimnames = dimnames(plot_matrix)
    )
  }

  heatmap_file <- file.path(
    project$analysis,
    paste0(
      aaa_safe_name(graph_main),
      if (!isTRUE(prepared$has_replicates)) "_by_samples" else "", ".png"
    )
  )

  heatmap_values <- as.matrix(plot_matrix)
  finite_heatmap_values <- heatmap_values[is.finite(heatmap_values)]

  if (length(finite_heatmap_values) == 0L) {
    heatmap_breaks <- seq(0, 1, length.out = 51L)
  } else {
    heatmap_range <- range(finite_heatmap_values, na.rm = TRUE)

    if (!all(is.finite(heatmap_range))) {
      heatmap_breaks <- seq(0, 1, length.out = 51L)
    } else if (isTRUE(all.equal(heatmap_range[1], heatmap_range[2]))) {
      centre <- heatmap_range[1]
      tolerance <- max(abs(centre) * 1e-8, 1e-8)
      heatmap_breaks <- seq(
        centre - tolerance,
        centre + tolerance,
        length.out = 51L
      )
    } else {
      heatmap_breaks <- seq(
        heatmap_range[1],
        heatmap_range[2],
        length.out = 51L
      )
    }
  }

  # A fixed width clips long row labels (classification categories can be
  # arbitrarily long descriptive strings) instead of the plot growing to fit
  # them; measure the actual rendered label width and size around it.
  heatmap_width <- max(3, 0.55 * ncol(heatmap_values)) +
    aaa_measure_text_width_inches(rownames(heatmap_values), 10) + 1.6

  heatmap_height <- max(4, 0.5 * nrow(plot_matrix) + 1.5)
  heatmap_args <- list(
    mat = heatmap_values,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = graph_main,
    angle_col = 315,
    display_numbers = plot_labels,
    fontsize_number = 10,
    number_color = "black",
    color = grDevices::colorRampPalette(
      c("white", "yellow", "orange", "purple")
    )(50),
    breaks = heatmap_breaks,
    border_color = "grey70"
  )

  if (is.null(graph_note) || !nzchar(trimws(graph_note))) {
    do.call(
      pheatmap::pheatmap,
      c(
        heatmap_args,
        list(
          filename = heatmap_file,
          width = heatmap_width,
          height = heatmap_height
        )
      )
    )
  } else {
    # pheatmap has no native footer. Add an optional registry-provided note to
    # the returned gtable; the default path remains unchanged for all analyses.
    heatmap_object <- do.call(
      pheatmap::pheatmap,
      c(heatmap_args, list(silent = TRUE))
    )
    footer_height <- grid::unit(0.45, "in")
    heatmap_gtable <- gtable::gtable_add_rows(
      heatmap_object$gtable,
      footer_height,
      pos = -1
    )
    heatmap_gtable <- gtable::gtable_add_grob(
      heatmap_gtable,
      grid::textGrob(
        graph_note,
        x = grid::unit(0.5, "npc"),
        just = "centre",
        gp = grid::gpar(fontsize = 8, col = "grey30")
      ),
      t = nrow(heatmap_gtable),
      l = 1,
      r = ncol(heatmap_gtable),
      name = "functional_heatmap_note"
    )
    grDevices::png(
      heatmap_file,
      width = heatmap_width,
      height = heatmap_height + 0.45,
      units = "in",
      res = 150
    )
    grid::grid.newpage()
    grid::grid.draw(heatmap_gtable)
    grDevices::dev.off()
  }

  summary_file <- file.path(project$analysis, "aaa_functional_potential_summary.xlsx")
  openxlsx::write.xlsx(
    list(
      Taxon_results = taxa,
      aaa_functional_abundance = functional_abundance,
      Heatmap_values = tibble::rownames_to_column(
        plot_matrix,
        "Potential"
      ),
      Registry_gene_roles = if (
        !is.null(biological_function_id)
      ) {
        aaa_registry_gene_catalogue(
          biological_function_id
        )
      } else {
        data.frame()
      },
      Cache_summary = data.frame(
        Metric = c(
          "Taxon/gene lookups reused from cache",
          "New NCBI/GFF lookups performed"
        ),
        Count = c(cache_hits, cache_misses),
        stringsAsFactors = FALSE
      )
    ),
    summary_file,
    overwrite = TRUE
  )

  metadata <- list(
    genes = genes,
    classification_function = deparse(substitute(classification_function)),
    abundance_type = abundance_type,
    replicate_design = prepared$replicates,
    n_replicates = prepared$n_replicates,
    cache_status = if (cache_hits > 0L) "HIT" else "MISS",
    cache_hits = cache_hits,
    cache_misses = cache_misses,
    evaluated_markers = sum(evaluated_markers),
    detected_markers = sum(detected_markers)
  )

  report_progress(
    "completed",
    "Functional-potential analysis completed.",
    1,
    1
  )

  result <- aaa_result(
    tables = list(
      taxa = taxa,
      functional_abundance = functional_abundance,
      heatmap_values = plot_matrix
    ),
    plots = list(heatmap = heatmap_file),
    files = c(summary = summary_file, heatmap = heatmap_file),
    output_dir = project$analysis,
    metadata = metadata
  )
  attr(result, "triple_a_cache_status") <- tolower(metadata$cache_status)
  result
}


aaa_functional_taxon_candidates <- function(
  taxonomy,
  genus = NA_character_,
  tax_level = NA_character_
) {
  clean_taxon <- function(x) {
    if (length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
      return(NA_character_)
    }

    value <- trimws(as.character(x))
    value <- sub("^[A-Za-z]__", "", value)
    value <- gsub("_", " ", value, fixed = TRUE)
    value <- trimws(gsub("\\s+", " ", value, perl = TRUE))

    invalid_exact <- c(
      "unknown", "unclassified", "uncultured", "uncultured bacterium",
      "uncultured archaeon", "bacterium", "archaeon", "metagenome",
      "environmental sample", "na", "none"
    )
    lower_value <- tolower(value)

    if (!nzchar(value) || lower_value %in% invalid_exact) {
      return(NA_character_)
    }

    if (grepl("(^|\\s)(sp\\.?|bacterium|archaeon|uncultured|unclassified)(\\s|$)",
      lower_value,
      perl = TRUE
    )) {
      return(NA_character_)
    }

    value
  }

  genus_value <- clean_taxon(genus)
  taxonomy_value <- clean_taxon(taxonomy)
  level_value <- tolower(trimws(as.character(tax_level)[1L]))

  queries <- character()
  ranks <- character()

  # The taxonomy stored in the prepared dataset is the lowest assigned rank.
  # Use it as a species query only when the recorded rank is Species and the
  # cleaned value is a genuine binomial. Otherwise do not infer a species.
  if (!is.na(taxonomy_value) && identical(level_value, "species")) {
    words <- strsplit(taxonomy_value, " ", fixed = TRUE)[[1L]]
    words <- words[nzchar(words)]
    if (length(words) >= 2L) {
      queries <- c(queries, taxonomy_value)
      ranks <- c(ranks, "Species")
    }
  }

  if (!is.na(genus_value)) {
    queries <- c(queries, genus_value)
    ranks <- c(ranks, "Genus")
  }

  if (!length(queries)) {
    return(data.frame(
      Query = character(),
      Rank = character(),
      stringsAsFactors = FALSE
    ))
  }

  keep <- !duplicated(tolower(queries))
  data.frame(
    Query = queries[keep],
    Rank = ranks[keep],
    stringsAsFactors = FALSE
  )
}


# Backward-compatible helper: return the preferred candidate only.
aaa_functional_taxon_query <- function(
  taxonomy,
  genus = NA_character_,
  tax_level = NA_character_
) {
  candidates <- aaa_functional_taxon_candidates(
    taxonomy = taxonomy,
    genus = genus,
    tax_level = tax_level
  )
  if (!nrow(candidates)) NA_character_ else candidates$Query[[1L]]
}


# NCBI E-utilities allow 3 requests/second without an API key, or 10/second
# with one (https://www.ncbi.nlm.nih.gov/books/NBK25497/). On a cold cache,
# taxon resolution issues 1-2 sequential requests per new taxon, so a large
# taxon table can take a while; registering a free NCBI API key (from the
# user's NCBI account "Settings" page) via
# options(triple_a_ncbi_api_key = "...") or the NCBI_API_KEY environment
# variable roughly triples that throughput without violating NCBI's rate
# limit, unlike running requests in parallel would.
aaa_ncbi_request_delay <- function() {
  key <- getOption("triple_a_ncbi_api_key", Sys.getenv("NCBI_API_KEY", ""))
  if (is.character(key) && nzchar(trimws(key))) {
    if (requireNamespace("rentrez", quietly = TRUE)) {
      rentrez::set_entrez_key(trimws(key))
    }
    return(0.11)
  }
  0.35
}


aaa_get_taxid <- function(name, verbose = FALSE) {
  if (is.na(name) || !nzchar(trimws(name))) {
    return(NA_character_)
  }
  original <- trimws(name)
  query_name <- gsub("_", " ", original)
  Sys.sleep(aaa_ncbi_request_delay())

  search <- tryCatch(
    rentrez::entrez_search(
      db = "taxonomy",
      term = paste0(query_name, "[Scientific Name]")
    ),
    error = function(e) {
      if (verbose) message("Taxonomy search failed for ", original, ": ", e$message)
      NULL
    }
  )
  if (is.null(search) || length(search$ids) == 0) {
    return(NA_character_)
  }
  as.character(search$ids[1])
}


aaa_get_reference_genome <- function(taxid, verbose = FALSE) {
  if (is.na(taxid) || !nzchar(as.character(taxid))) {
    return(NA_character_)
  }
  Sys.sleep(aaa_ncbi_request_delay())
  search <- tryCatch(
    rentrez::entrez_search(
      db = "assembly",
      term = paste0("txid", taxid, "[Organism] AND latest_refseq[filter]"),
      retmax = 100
    ),
    error = function(e) {
      if (verbose) {
        message(
          "Assembly search failed for TaxID ", taxid, ": ",
          e$message
        )
      }
      NULL
    }
  )
  if (is.null(search) || length(search$ids) == 0) {
    return(NA_character_)
  }

  summaries <- tryCatch(
    rentrez::entrez_summary(db = "assembly", id = search$ids),
    error = function(e) NULL
  )
  if (is.null(summaries)) {
    return(NA_character_)
  }

  accession <- rentrez::extract_from_esummary(summaries, "assemblyaccession")
  category <- rentrez::extract_from_esummary(summaries, "refseq_category")
  level <- rentrez::extract_from_esummary(summaries, "assemblystatus")

  df <- data.frame(
    accession = accession,
    category = ifelse(is.na(category), "", category),
    level = ifelse(is.na(level), "", level),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$accession) & nzchar(df$accession), , drop = FALSE]
  if (nrow(df) == 0) {
    return(NA_character_)
  }

  df$priority <- dplyr::case_when(
    df$category == "reference genome" ~ 1L,
    df$category == "representative genome" ~ 2L,
    df$level == "Complete Genome" ~ 3L,
    TRUE ~ 4L
  )
  df$level_priority <- match(
    df$level,
    c("Complete Genome", "Chromosome", "Scaffold", "Contig")
  )
  df$level_priority[is.na(df$level_priority)] <- 5L
  df <- df[order(df$priority, df$level_priority), , drop = FALSE]
  as.character(df$accession[1])
}


aaa_find_cached_gff <- function(
  accession,
  gff_directory
) {
  if (is.na(accession) ||
    !nzchar(accession) ||
    !dir.exists(gff_directory)) {
    return(character())
  }

  # A previous version of this pattern double-escaped the literal dots in
  # "_genomic.gff.gz$" (\\\\. instead of \\.), which made the regex require
  # a literal backslash before "gff"/"gz" in the filename. No real cached
  # file could ever match, so aaa_get_gff() always treated every accession
  # as uncached and re-downloaded it from NCBI on every run, even when the
  # GFF was already sitting in Cache/GFF from a previous analysis.
  files <- list.files(
    gff_directory,
    pattern = "_genomic\\.gff\\.gz$",
    full.names = TRUE
  )

  # Require the accession to be followed by "_" (as in
  # "<accession>_<assembly name>_genomic.gff.gz") so an accession that is a
  # string prefix of another (e.g. GCF_000005825.2 vs GCF_000005825.20)
  # cannot match the wrong cached file.
  matches <- files[startsWith(basename(files), paste0(accession, "_"))]
  matches[file.exists(matches)]
}


aaa_read_cached_gff <- function(
  path,
  accession
) {
  gff <- tryCatch(
    utils::read.delim(
      gzfile(path),
      sep = "\t",
      header = FALSE,
      comment.char = "#",
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      warning(
        "Failed to read cached GFF for ",
        accession,
        ": ",
        e$message
      )

      NULL
    }
  )

  if (is.null(gff) ||
    ncol(gff) != 9) {
    return(NULL)
  }

  names(gff) <- c(
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
  )

  attr(
    gff,
    "triple_a_cache_status"
  ) <- "cache_hit"

  gff
}


aaa_get_gff <- function(
  accession,
  project,
  verbose = FALSE
) {
  if (is.na(accession) ||
    !nzchar(accession)) {
    return(NULL)
  }

  # Cache directories are intentionally lazy, but they must exist before
  # locks, temporary downloads or GFF files are created.
  dir.create(project$cache, recursive = TRUE, showWarnings = FALSE)
  dir.create(project$gff, recursive = TRUE, showWarnings = FALSE)

  cached <- aaa_find_cached_gff(
    accession,
    project$gff
  )

  if (length(cached) > 0) {
    return(aaa_read_cached_gff(
      cached[1],
      accession
    ))
  }

  lock_file <- file.path(
    project$cache,
    paste0(
      gsub(
        "[^A-Za-z0-9_.-]",
        "_",
        accession
      ),
      ".gff.lock"
    )
  )

  aaa_acquire_cache_lock(
    lock_file
  )

  on.exit(
    aaa_release_cache_lock(
      lock_file
    ),
    add = TRUE
  )

  cached <- aaa_find_cached_gff(
    accession,
    project$gff
  )

  if (length(cached) > 0) {
    return(aaa_read_cached_gff(
      cached[1],
      accession
    ))
  }

  search <- tryCatch(
    rentrez::entrez_search(
      db = "assembly",
      term = accession,
      retmax = 1
    ),
    error = function(e) NULL
  )

  if (is.null(search) ||
    length(search$ids) == 0) {
    return(NULL)
  }

  summary <- tryCatch(
    rentrez::entrez_summary(
      db = "assembly",
      id = search$ids[1]
    ),
    error = function(e) NULL
  )

  if (is.null(summary)) {
    return(NULL)
  }

  ftp <- summary$ftppath_refseq

  if (is.null(ftp) ||
    !nzchar(ftp)) {
    return(NULL)
  }

  ftp <- sub(
    "^ftp://",
    "https://",
    ftp
  )

  assembly <- basename(ftp)
  filename <- paste0(
    assembly,
    "_genomic.gff.gz"
  )

  destination <- file.path(
    project$gff,
    filename
  )

  temporary <- tempfile(
    pattern = paste0(
      filename,
      "_"
    ),
    tmpdir = project$gff,
    fileext = ".download"
  )

  if (verbose) {
    message(
      "Downloading ",
      filename
    )
  }

  ok <- tryCatch(
    {
      utils::download.file(
        paste0(
          ftp,
          "/",
          filename
        ),
        temporary,
        mode = "wb",
        quiet = !verbose
      )

      isTRUE(
        file.info(temporary)$size > 0
      )
    },
    error = function(e) {
      warning(
        "Failed to download GFF for ",
        accession,
        ": ",
        e$message
      )

      FALSE
    }
  )

  if (!isTRUE(ok)) {
    unlink(
      temporary,
      force = TRUE
    )

    return(NULL)
  }

  if (!file.rename(
    temporary,
    destination
  )) {
    file.copy(
      temporary,
      destination,
      overwrite = TRUE
    )

    unlink(
      temporary,
      force = TRUE
    )
  }

  gff <- aaa_read_cached_gff(
    destination,
    accession
  )

  if (!is.null(gff)) {
    attr(
      gff,
      "triple_a_cache_status"
    ) <- "downloaded"
  }

  gff
}


aaa_search_gene <- function(gff, gene, verbose = FALSE) {
  gene <- trimws(gene)
  if (!nzchar(gene) || is.null(gff) || nrow(gff) == 0) {
    return(FALSE)
  }
  if (!"attributes" %in% names(gff)) {
    stop("Invalid GFF object: missing 'attributes' column.")
  }
  if (!exists("gene_aliases", inherits = TRUE)) {
    stop("'gene_aliases' was not found. Source dictionary.R first.")
  }

  # gene_aliases[[gene]] comes from a dictionary that can be extended at
  # runtime with user-supplied custom aliases (Resources/FunctionalDB/
  # CustomGeneAliases.json); coerce explicitly so a malformed or unexpected
  # entry there cannot reach nzchar() as anything other than character.
  aliases <- as.character(unique(c(gene, gene_aliases[[gene]])))
  aliases <- aliases[!is.na(aliases) & nzchar(aliases)]

  found <- any(vapply(
    aliases,
    function(term) {
      escaped <- stringr::str_replace_all(
        term, "([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1"
      )
      pattern <- paste0(
        "(^|;)(gene|name|gene_name|gene_synonym|locus_tag|old_locus_tag)=",
        escaped, "([;,]|$)|product=[^;]*\\b", escaped, "\\b"
      )
      any(grepl(pattern, gff$attributes,
        perl = TRUE,
        ignore.case = TRUE
      ), na.rm = TRUE)
    },
    logical(1)
  ))

  if (verbose) {
    message(
      "Gene ", gene, ": ",
      ifelse(found, "found", "not found")
    )
  }
  found
}


infer_functional_potential <- aaa_functional_potential
