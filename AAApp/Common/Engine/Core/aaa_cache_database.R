# =============================================================================
# Persistent SQLite cache
#
# Stores taxonomic references and marker-gene search results with stable types.
# The database is shared by every Run_* execution.
# =============================================================================

aaa_reference_cache_database_path <- function(project) {
  file.path(
    project$reference_cache,
    "GenomeCache.sqlite"
  )
}


aaa_require_cache_database_packages <- function() {
  missing <- c("DBI", "RSQLite")[
    !vapply(
      c("DBI", "RSQLite"),
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing) > 0) {
    stop(
      "The cache requires R packages: ",
      paste(missing, collapse = ", "),
      ". Install them before running the workflow."
    )
  }

  invisible(TRUE)
}


aaa_reference_cache_connect <- function(
  project = NULL,
  database_file = NULL
) {
  aaa_require_cache_database_packages()

  if (is.null(database_file)) {
    if (is.null(project)) {
      stop(
        "Provide either 'project' or 'database_file'."
      )
    }

    database_file <-
      aaa_reference_cache_database_path(project)
  }

  dir.create(
    dirname(database_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  connection <- DBI::dbConnect(
    RSQLite::SQLite(),
    database_file
  )

  DBI::dbExecute(
    connection,
    "PRAGMA busy_timeout = 30000"
  )

  try(
    DBI::dbExecute(
      connection,
      "PRAGMA journal_mode = WAL"
    ),
    silent = TRUE
  )

  DBI::dbExecute(
    connection,
    "PRAGMA foreign_keys = ON"
  )

  aaa_reference_cache_initialize(
    connection
  )

  connection
}


aaa_reference_cache_initialize <- function(connection) {
  DBI::dbExecute(
    connection,
    paste(
      "CREATE TABLE IF NOT EXISTS reference_genomes (",
      "taxonomy TEXT PRIMARY KEY NOT NULL,",
      "taxid TEXT,",
      "genus TEXT,",
      "tax_level TEXT,",
      "reference_genome TEXT,",
      "updated_at TEXT NOT NULL",
      ")"
    )
  )

  DBI::dbExecute(
    connection,
    paste(
      "CREATE TABLE IF NOT EXISTS gene_results (",
      "taxonomy TEXT NOT NULL,",
      "gene TEXT NOT NULL,",
      "found INTEGER NOT NULL CHECK(found IN (0, 1)),",
      "updated_at TEXT NOT NULL,",
      "PRIMARY KEY (taxonomy, gene),",
      "FOREIGN KEY (taxonomy)",
      "REFERENCES reference_genomes(taxonomy)",
      "ON DELETE CASCADE",
      ")"
    )
  )

  DBI::dbExecute(
    connection,
    paste(
      "CREATE TABLE IF NOT EXISTS cache_metadata (",
      "key TEXT PRIMARY KEY NOT NULL,",
      "value TEXT,",
      "updated_at TEXT NOT NULL",
      ")"
    )
  )

  DBI::dbExecute(
    connection,
    paste(
      "CREATE INDEX IF NOT EXISTS",
      "idx_reference_genomes_taxid",
      "ON reference_genomes(taxid)"
    )
  )

  DBI::dbExecute(
    connection,
    paste(
      "CREATE INDEX IF NOT EXISTS",
      "idx_gene_results_gene",
      "ON gene_results(gene)"
    )
  )

  invisible(connection)
}


aaa_cache_as_text <- function(x) {
  result <- as.character(x)
  result[is.na(x)] <- NA_character_
  result
}


aaa_cache_as_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  normalized <- toupper(
    trimws(
      as.character(x)
    )
  )

  result <- rep(
    NA,
    length(normalized)
  )

  result[
    normalized %in%
      c("TRUE", "T", "1", "YES")
  ] <- TRUE

  result[
    normalized %in%
      c("FALSE", "F", "0", "NO")
  ] <- FALSE

  result
}


aaa_normalize_reference_cache <- function(cache) {
  required <- c(
    "TaxID",
    "Taxonomy",
    "Genus",
    "Tax_level",
    "Reference_genome"
  )

  if (is.null(cache) ||
    !is.data.frame(cache)) {
    cache <- data.frame(
      TaxID = character(),
      Taxonomy = character(),
      Genus = character(),
      Tax_level = character(),
      Reference_genome = character(),
      stringsAsFactors = FALSE
    )
  }

  for (column in setdiff(
    required,
    names(cache)
  )) {
    cache[[column]] <- NA_character_
  }

  for (column in required) {
    cache[[column]] <-
      aaa_cache_as_text(
        cache[[column]]
      )
  }

  gene_columns <- setdiff(
    names(cache),
    required
  )

  for (column in gene_columns) {
    cache[[column]] <-
      aaa_cache_as_logical(
        cache[[column]]
      )
  }

  valid_taxonomy <- !is.na(cache$Taxonomy) &
    nzchar(trimws(cache$Taxonomy))

  cache <- cache[
    valid_taxonomy,
    c(required, gene_columns),
    drop = FALSE
  ]

  if (nrow(cache) > 0) {
    cache <- cache[
      !duplicated(
        cache$Taxonomy,
        fromLast = TRUE
      ), ,
      drop = FALSE
    ]
  }

  rownames(cache) <- NULL
  cache
}


aaa_reference_cache_upsert <- function(
  connection,
  cache
) {
  cache <- aaa_normalize_reference_cache(
    cache
  )

  if (nrow(cache) == 0) {
    return(invisible(cache))
  }

  timestamp <- format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )

  base_columns <- c(
    "TaxID",
    "Taxonomy",
    "Genus",
    "Tax_level",
    "Reference_genome"
  )

  gene_columns <- setdiff(
    names(cache),
    base_columns
  )

  DBI::dbWithTransaction(
    connection,
    {
      reference_sql <- paste(
        "INSERT INTO reference_genomes",
        "(taxonomy, taxid, genus, tax_level,",
        "reference_genome, updated_at)",
        "VALUES (?, ?, ?, ?, ?, ?)",
        "ON CONFLICT(taxonomy) DO UPDATE SET",
        "taxid = excluded.taxid,",
        "genus = excluded.genus,",
        "tax_level = excluded.tax_level,",
        "reference_genome = excluded.reference_genome,",
        "updated_at = excluded.updated_at"
      )

      gene_sql <- paste(
        "INSERT INTO gene_results",
        "(taxonomy, gene, found, updated_at)",
        "VALUES (?, ?, ?, ?)",
        "ON CONFLICT(taxonomy, gene) DO UPDATE SET",
        "found = excluded.found,",
        "updated_at = excluded.updated_at"
      )

      if (nrow(cache) > 0) {
        reference_statement <- DBI::dbSendStatement(
          connection,
          reference_sql
        )

        DBI::dbBind(
          reference_statement,
          params = list(
            cache$Taxonomy,
            cache$TaxID,
            cache$Genus,
            cache$Tax_level,
            cache$Reference_genome,
            rep(timestamp, nrow(cache))
          )
        )

        DBI::dbClearResult(reference_statement)
      }

      if (length(gene_columns) > 0) {
        # Long-format table of (taxonomy, gene, found) rows across every
        # gene column, skipping cells that were never evaluated (NA).
        # A single dbBind() call then executes all rows as one batch
        # instead of issuing one dbExecute() per (row, gene) combination.
        gene_rows <- lapply(gene_columns, function(gene) {
          values <- cache[[gene]]
          keep <- !is.na(values)

          if (!any(keep)) {
            return(NULL)
          }

          data.frame(
            taxonomy = cache$Taxonomy[keep],
            gene = gene,
            found = as.integer(values[keep]),
            stringsAsFactors = FALSE
          )
        })

        gene_rows <- do.call(
          rbind,
          gene_rows[!vapply(gene_rows, is.null, logical(1))]
        )

        if (!is.null(gene_rows) && nrow(gene_rows) > 0) {
          gene_statement <- DBI::dbSendStatement(
            connection,
            gene_sql
          )

          DBI::dbBind(
            gene_statement,
            params = list(
              gene_rows$taxonomy,
              gene_rows$gene,
              gene_rows$found,
              rep(timestamp, nrow(gene_rows))
            )
          )

          DBI::dbClearResult(gene_statement)
        }
      }

      DBI::dbExecute(
        connection,
        paste(
          "INSERT INTO cache_metadata",
          "(key, value, updated_at)",
          "VALUES ('schema_version', '1', ?)",
          "ON CONFLICT(key) DO UPDATE SET",
          "value = excluded.value,",
          "updated_at = excluded.updated_at"
        ),
        params = list(timestamp)
      )
    }
  )

  invisible(cache)
}


aaa_reference_cache_read <- function(
  connection
) {
  references <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT",
      "taxid AS TaxID,",
      "taxonomy AS Taxonomy,",
      "genus AS Genus,",
      "tax_level AS Tax_level,",
      "reference_genome AS Reference_genome",
      "FROM reference_genomes",
      "ORDER BY taxonomy"
    )
  )

  references <- aaa_normalize_reference_cache(
    references
  )

  genes <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT taxonomy, gene, found",
      "FROM gene_results",
      "ORDER BY gene, taxonomy"
    )
  )

  if (nrow(references) == 0 ||
    nrow(genes) == 0) {
    return(references)
  }

  gene_names <- unique(
    genes$gene
  )

  for (gene in gene_names) {
    references[[gene]] <- NA
  }

  row_index <- match(
    genes$taxonomy,
    references$Taxonomy
  )

  matched <- !is.na(row_index)
  genes <- genes[matched, , drop = FALSE]
  row_index <- row_index[matched]

  # One vectorized column assignment per gene instead of one data.frame
  # cell assignment per (taxon, gene) row; equivalent result, far fewer
  # (and far cheaper) assignments when the cache holds many taxa/genes.
  for (gene in gene_names) {
    rows_for_gene <- genes$gene == gene
    references[[gene]][row_index[rows_for_gene]] <- as.logical(
      genes$found[rows_for_gene]
    )
  }

  aaa_normalize_reference_cache(
    references
  )
}


aaa_load_reference_cache <- function(project) {
  database_file <-
    aaa_reference_cache_database_path(
      project
    )

  connection <-
    aaa_reference_cache_connect(
      database_file = database_file
    )

  on.exit(
    DBI::dbDisconnect(connection),
    add = TRUE
  )

  list(
    cache =
      aaa_reference_cache_read(
        connection
      ),
    file = database_file,
    backend = "SQLite"
  )
}


aaa_update_reference_cache <- function(
  cache,
  cache_file
) {
  connection <-
    aaa_reference_cache_connect(
      database_file = cache_file
    )

  on.exit(
    DBI::dbDisconnect(connection),
    add = TRUE
  )

  normalized <-
    aaa_normalize_reference_cache(
      cache
    )

  aaa_reference_cache_upsert(
    connection,
    normalized
  )

  aaa_reference_cache_read(
    connection
  )
}
