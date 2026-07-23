# =============================================================================
# Triple_A portable reference-annotation cache manager
# =============================================================================

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

aaa_cache_manager_require_packages <- function() {
  missing <- c("DBI", "RSQLite")[!vapply(c("DBI", "RSQLite"), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Cache management requires R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

aaa_cache_root_from_database <- function(database_file) {
  normalizePath(dirname(database_file), winslash = "/", mustWork = FALSE)
}

aaa_cache_normalize_relative_path <- function(path, cache_root) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(trimws(path))) {
    return(NA_character_)
  }
  value <- gsub("\\\\", "/", trimws(as.character(path)))
  root <- sub("/+$", "", normalizePath(cache_root, winslash = "/", mustWork = FALSE))
  absolute <- normalizePath(value, winslash = "/", mustWork = FALSE)
  if (startsWith(tolower(absolute), paste0(tolower(root), "/"))) {
    value <- substring(absolute, nchar(root) + 2L)
  }
  value <- sub("^\\./", "", value)
  value <- sub("^Cache/", "", value, ignore.case = TRUE)
  if (grepl("(^|/)\\.\\.(/|$)", value)) stop("Unsafe cache path: ", path, call. = FALSE)
  value
}

aaa_cache_manager_initialize <- function(connection) {
  # Single source of truth for the core cache schema. The reference-genome /
  # gene-results / cache-metadata tables are created only by
  # aaa_reference_cache_initialize() in aaa_cache_database.R, which both the
  # analysis engine and the CacheManager app load. This replaces a duplicated
  # inline CREATE TABLE fallback that could silently drift from the canonical.
  if (!exists("aaa_reference_cache_initialize", mode = "function")) {
    stop(
      "aaa_reference_cache_initialize() is not loaded; source ",
      "aaa_cache_database.R before aaa_genome_cache_manager.R.",
      call. = FALSE
    )
  }
  aaa_reference_cache_initialize(connection)
  # cache_operations is an audit log specific to the CacheManager layer.
  DBI::dbExecute(connection, "CREATE TABLE IF NOT EXISTS cache_operations (id INTEGER PRIMARY KEY AUTOINCREMENT, operation TEXT NOT NULL, source TEXT, summary TEXT, created_at TEXT NOT NULL)")
  DBI::dbExecute(connection, "CREATE INDEX IF NOT EXISTS idx_reference_genomes_taxonomy ON reference_genomes(taxonomy)")
  DBI::dbExecute(connection, "CREATE INDEX IF NOT EXISTS idx_reference_genomes_genus ON reference_genomes(genus)")
  DBI::dbExecute(connection, "CREATE INDEX IF NOT EXISTS idx_reference_genomes_taxid_v12 ON reference_genomes(taxid)")
  now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  DBI::dbExecute(connection, "INSERT INTO cache_metadata(key,value,updated_at) VALUES('schema_version','2',?) ON CONFLICT(key) DO UPDATE SET value=excluded.value, updated_at=excluded.updated_at", params = list(now))
  invisible(connection)
}

aaa_cache_manager_connect <- function(database_file) {
  aaa_cache_manager_require_packages()
  dir.create(dirname(database_file), recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), database_file)
  DBI::dbExecute(con, "PRAGMA busy_timeout=30000")
  DBI::dbExecute(con, "PRAGMA foreign_keys=ON")
  try(DBI::dbExecute(con, "PRAGMA journal_mode=WAL"), silent = TRUE)
  aaa_cache_manager_initialize(con)
  con
}

aaa_cache_backup <- function(database_file, backup_dir = file.path(dirname(database_file), "Backups")) {
  if (!file.exists(database_file)) {
    return(NA_character_)
  }
  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  destination <- file.path(backup_dir, paste0("GenomeCache_", stamp, ".sqlite"))
  if (!file.copy(database_file, destination, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)) {
    stop("Could not create cache backup: ", destination, call. = FALSE)
  }
  normalizePath(destination, winslash = "/", mustWork = TRUE)
}

aaa_cache_list_files <- function(directory) {
  if (!dir.exists(directory)) {
    return(character())
  }
  files <- list.files(directory, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  normalizePath(files[file.exists(files) & !dir.exists(files)], winslash = "/", mustWork = FALSE)
}

aaa_cache_file_inventory <- function(cache_root) {
  files <- aaa_cache_list_files(file.path(cache_root, "GFF"))
  data.frame(
    kind = rep("gff", length(files)),
    path = files,
    stringsAsFactors = FALSE
  )
}

aaa_cache_statistics <- function(database_file) {
  cache_root <- aaa_cache_root_from_database(database_file)
  con <- aaa_cache_manager_connect(database_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  refs <- DBI::dbGetQuery(con, "SELECT COUNT(*) n, COUNT(DISTINCT genus) genera, COUNT(DISTINCT taxid) taxids FROM reference_genomes")
  genes <- DBI::dbGetQuery(con, "SELECT COUNT(*) n, COUNT(DISTINCT gene) genes FROM gene_results")
  inventory <- aaa_cache_file_inventory(cache_root)
  gff_sizes <- if (nrow(inventory)) file.info(inventory$path)$size else numeric()
  gff_bytes <- sum(gff_sizes, na.rm = TRUE)

  database_candidates <- unique(c(
    database_file,
    paste0(database_file, "-wal"),
    paste0(database_file, "-shm")
  ))
  database_candidates <- database_candidates[file.exists(database_candidates)]
  database_bytes <- if (length(database_candidates)) {
    sum(file.info(database_candidates)$size, na.rm = TRUE)
  } else {
    0
  }

  list(
    references = as.integer(refs$n[1]),
    genera = as.integer(refs$genera[1]),
    taxids = as.integer(refs$taxids[1]),
    gene_results = as.integer(genes$n[1]),
    genes = as.integer(genes$genes[1]),
    gff_files = sum(inventory$kind == "gff"),
    gff_bytes = as.numeric(gff_bytes),
    database_bytes = as.numeric(database_bytes),
    bytes = as.numeric(gff_bytes + database_bytes),
    database_file = database_file
  )
}

aaa_cache_reference_table <- function(database_file) {
  con <- aaa_cache_manager_connect(database_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbGetQuery(con, "SELECT taxonomy, taxid, genus, tax_level, reference_genome, updated_at FROM reference_genomes ORDER BY taxonomy")
}

aaa_cache_operation_log <- function(database_file, limit = 100L) {
  con <- aaa_cache_manager_connect(database_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbGetQuery(con, "SELECT operation, source, summary, created_at FROM cache_operations ORDER BY id DESC LIMIT ?", params = list(as.integer(limit)))
}

aaa_cache_record_operation <- function(connection, operation, source = "", summary = "") {
  DBI::dbExecute(connection, "INSERT INTO cache_operations(operation,source,summary,created_at) VALUES(?,?,?,?)", params = list(operation, source, summary, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

aaa_cache_file_checksum <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(path)[[1L]])
}

aaa_cache_copy_tree <- function(source_dir, destination_dir, conflict = c("keep_existing", "replace_newer")) {
  conflict <- match.arg(conflict)

  added <- 0L
  skipped <- 0L
  replaced <- 0L
  deduplicated <- 0L
  conflicts <- character()

  if (!dir.exists(source_dir)) {
    return(list(
      added = 0L,
      skipped = 0L,
      replaced = 0L,
      deduplicated = 0L,
      conflicts = character()
    ))
  }

  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  source_files <- aaa_cache_list_files(source_dir)

  for (src in source_files) {
    dst <- file.path(destination_dir, basename(src))

    if (!file.exists(dst)) {
      copied <- file.copy(
        src,
        dst,
        overwrite = FALSE,
        copy.mode = TRUE,
        copy.date = TRUE
      )

      if (!isTRUE(copied)) {
        stop("Could not copy cache file: ", src, call. = FALSE)
      }

      added <- added + 1L
      next
    }

    src_checksum <- aaa_cache_file_checksum(src)
    dst_checksum <- aaa_cache_file_checksum(dst)

    if (
      !is.na(src_checksum) &&
        !is.na(dst_checksum) &&
        identical(src_checksum, dst_checksum)
    ) {
      skipped <- skipped + 1L
      deduplicated <- deduplicated + 1L
      next
    }

    conflicts <- c(conflicts, basename(src))

    replace <- identical(conflict, "replace_newer") &&
      isTRUE(file.info(src)$mtime > file.info(dst)$mtime)

    if (replace) {
      copied <- file.copy(
        src,
        dst,
        overwrite = TRUE,
        copy.mode = TRUE,
        copy.date = TRUE
      )

      if (!isTRUE(copied)) {
        stop("Could not replace cache file: ", dst, call. = FALSE)
      }

      replaced <- replaced + 1L
    } else {
      skipped <- skipped + 1L
    }
  }

  list(
    added = added,
    skipped = skipped,
    replaced = replaced,
    deduplicated = deduplicated,
    conflicts = unique(conflicts)
  )
}

aaa_cache_merge <- function(target_database, source_database, source_cache_root = dirname(source_database), conflict = c("keep_existing", "replace_newer"), create_backup = TRUE) {
  conflict <- match.arg(conflict)
  if (!file.exists(source_database)) stop("Source GenomeCache.sqlite was not found.", call. = FALSE)
  target_root <- aaa_cache_root_from_database(target_database)
  source_root <- normalizePath(source_cache_root, winslash = "/", mustWork = TRUE)
  dir.create(file.path(target_root, "GFF"), recursive = TRUE, showWarnings = FALSE)
  backup <- if (isTRUE(create_backup)) aaa_cache_backup(target_database) else NA_character_
  target <- aaa_cache_manager_connect(target_database)
  on.exit(DBI::dbDisconnect(target), add = TRUE)
  source <- aaa_cache_manager_connect(source_database)
  on.exit(DBI::dbDisconnect(source), add = TRUE)
  source_refs <- DBI::dbGetQuery(source, "SELECT taxonomy,taxid,genus,tax_level,reference_genome,updated_at FROM reference_genomes")
  source_genes <- DBI::dbGetQuery(source, "SELECT taxonomy,gene,found,updated_at FROM gene_results")
  existing <- DBI::dbGetQuery(target, "SELECT taxonomy,updated_at FROM reference_genomes")
  added <- updated <- skipped <- 0L
  DBI::dbWithTransaction(target, {
    # Decide which reference rows to write with vectorized lookups, then
    # issue one batched dbBind() for all of them instead of one dbExecute()
    # per row (matches the batching pattern used by aaa_reference_cache_upsert()).
    idx <- match(source_refs$taxonomy, existing$taxonomy)
    is_new <- is.na(idx)
    is_newer <- !is_new &
      identical(conflict, "replace_newer") &
      as.character(source_refs$updated_at) > as.character(existing$updated_at[idx])
    should_write <- is_new | is_newer
    skipped <- sum(!should_write)
    added <- sum(is_new)
    updated <- sum(is_newer)

    write_rows <- source_refs[should_write, , drop = FALSE]
    if (nrow(write_rows) > 0) {
      relative_paths <- vapply(
        write_rows$reference_genome,
        aaa_cache_normalize_relative_path,
        character(1),
        source_root
      )
      reference_statement <- DBI::dbSendStatement(
        target,
        "INSERT INTO reference_genomes(taxonomy,taxid,genus,tax_level,reference_genome,updated_at) VALUES(?,?,?,?,?,?) ON CONFLICT(taxonomy) DO UPDATE SET taxid=excluded.taxid,genus=excluded.genus,tax_level=excluded.tax_level,reference_genome=excluded.reference_genome,updated_at=excluded.updated_at"
      )
      DBI::dbBind(reference_statement, params = list(
        write_rows$taxonomy, write_rows$taxid, write_rows$genus,
        write_rows$tax_level, relative_paths, write_rows$updated_at
      ))
      DBI::dbClearResult(reference_statement)
    }

    if (nrow(source_genes) > 0) {
      gene_statement <- DBI::dbSendStatement(
        target,
        "INSERT INTO gene_results(taxonomy,gene,found,updated_at) VALUES(?,?,?,?) ON CONFLICT(taxonomy,gene) DO UPDATE SET found=CASE WHEN excluded.updated_at > gene_results.updated_at THEN excluded.found ELSE gene_results.found END, updated_at=MAX(gene_results.updated_at,excluded.updated_at)"
      )
      DBI::dbBind(gene_statement, params = list(
        source_genes$taxonomy, source_genes$gene,
        source_genes$found, source_genes$updated_at
      ))
      DBI::dbClearResult(gene_statement)
    }
  })
  gff <- aaa_cache_copy_tree(file.path(source_root, "GFF"), file.path(target_root, "GFF"), conflict)
  summary <- list(
    references_added = added, references_updated = updated, references_skipped = skipped,
    gff_added = gff$added,
    duplicate_files_ignored = gff$deduplicated,
    files_replaced = gff$replaced,
    file_conflicts = gff$conflicts, backup = backup
  )
  aaa_cache_record_operation(target, "merge", source_database, paste(names(summary)[1:8], unlist(summary[1:8]), collapse = "; "))
  summary
}

aaa_cache_verify <- function(database_file) {
  root <- aaa_cache_root_from_database(database_file)
  con <- aaa_cache_manager_connect(database_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  integrity <- DBI::dbGetQuery(con, "PRAGMA integrity_check")[[1L]]
  refs <- DBI::dbGetQuery(con, "SELECT taxonomy,reference_genome FROM reference_genomes ORDER BY taxonomy")
  inventory <- aaa_cache_file_inventory(root)
  accessions <- trimws(as.character(refs$reference_genome))
  accessions[!nzchar(accessions)] <- NA_character_

  # reference_genome stores an NCBI assembly accession, not a local path.
  # A reference is backed locally when a GFF filename contains that
  # accession. This preserves compatibility with the existing GFF downloader.
  # Build the accession x filename substring-match matrix once and derive
  # both directions (local coverage, orphan files) from it, instead of two
  # independent O(accessions x files) scans.
  filenames <- if (nrow(inventory)) basename(inventory$path) else character(0)
  has_accession <- which(!is.na(accessions))
  match_matrix <- matrix(FALSE, nrow = length(accessions), ncol = length(filenames))
  if (length(has_accession) > 0 && length(filenames) > 0) {
    match_matrix[has_accession, ] <- t(vapply(
      accessions[has_accession],
      function(accession) grepl(accession, filenames, fixed = TRUE),
      logical(length(filenames))
    ))
  }

  has_local_file <- if (length(filenames)) rowSums(match_matrix) > 0 else rep(FALSE, length(accessions))
  missing <- refs[!is.na(accessions) & !has_local_file, , drop = FALSE]

  referenced_file <- if (length(filenames)) colSums(match_matrix) > 0 else logical(0)
  orphans <- inventory[!referenced_file, , drop = FALSE]
  valid <- all(tolower(integrity) == "ok")
  list(valid = valid, sqlite_integrity = integrity, missing_files = missing, orphan_files = orphans, references = nrow(refs), local_coverage = sum(has_local_file), coverage_warnings = nrow(missing))
}

aaa_cache_clean_orphans <- function(database_file, dry_run = TRUE) {
  check <- aaa_cache_verify(database_file)
  files <- check$orphan_files$path
  removed <- 0L
  if (!dry_run && length(files)) {
    aaa_cache_backup(database_file)
    removed <- sum(vapply(files, unlink, integer(1), force = TRUE) == 0L)
    con <- aaa_cache_manager_connect(database_file)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    aaa_cache_record_operation(con, "clean_orphans", "", paste("removed", removed, "files"))
  }
  list(dry_run = dry_run, candidates = length(files), removed = removed, files = files)
}

# Reference-genome assignments are cached indefinitely (NCBI's "reference
# genome" for a species can change over time), so give users a way to see
# and manually clear old entries. This is deliberately manual, not an
# automatic re-check on every run: re-verifying every cached taxon against
# NCBI on each analysis is exactly the always-re-download regression that
# was fixed earlier, just moved from the GFF layer to the taxonomy layer.
aaa_cache_stale_references <- function(database_file, max_age_days) {
  con <- aaa_cache_manager_connect(database_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  cutoff <- format(Sys.time() - max_age_days * 86400, "%Y-%m-%d %H:%M:%S")
  DBI::dbGetQuery(
    con,
    "SELECT taxonomy, taxid, genus, tax_level, reference_genome, updated_at FROM reference_genomes WHERE updated_at < ? ORDER BY updated_at",
    params = list(cutoff)
  )
}

aaa_cache_prune_stale_references <- function(database_file, max_age_days, dry_run = TRUE) {
  candidates <- aaa_cache_stale_references(database_file, max_age_days)
  removed <- 0L
  if (!dry_run && nrow(candidates)) {
    aaa_cache_backup(database_file)
    con <- aaa_cache_manager_connect(database_file)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    statement <- DBI::dbSendStatement(con, "DELETE FROM reference_genomes WHERE taxonomy = ?")
    DBI::dbBind(statement, params = list(candidates$taxonomy))
    DBI::dbClearResult(statement)
    removed <- nrow(candidates)
    aaa_cache_record_operation(
      con, "prune_stale_references", "",
      paste0("removed ", removed, " reference(s) older than ", max_age_days, " day(s)")
    )
  }
  list(dry_run = dry_run, max_age_days = max_age_days, candidates = nrow(candidates), removed = removed, references = candidates)
}

aaa_cache_export <- function(database_file, destination_zip) {
  root <- aaa_cache_root_from_database(database_file)
  check <- aaa_cache_verify(database_file)
  stage <- tempfile("TripleA_GenomeCache_")
  dir.create(stage, recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  file.copy(database_file, file.path(stage, "GenomeCache.sqlite"), copy.mode = TRUE, copy.date = TRUE)
  if (dir.exists(file.path(root, "GFF"))) file.copy(file.path(root, "GFF"), stage, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
  manifest <- c("Triple_A reference annotation cache package", "cache_format_version=2", paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0("integrity=", if (check$valid) "valid" else "warnings"), paste0("references=", check$references))
  writeLines(manifest, file.path(stage, "cache_manifest.txt"))
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(stage)
  destination_zip <- normalizePath(destination_zip, winslash = "/", mustWork = FALSE)
  if (file.exists(destination_zip)) unlink(destination_zip)
  utils::zip(destination_zip, files = list.files(stage, all.files = FALSE), flags = "-r9X")
  normalizePath(destination_zip, winslash = "/", mustWork = TRUE)
}


aaa_cache_package_preview <- function(package_file) {
  if (!file.exists(package_file)) stop("Cache package not found: ", package_file, call. = FALSE)

  extension <- tolower(tools::file_ext(package_file))
  package_size <- unname(file.info(package_file)$size)

  if (identical(extension, "sqlite")) {
    return(list(
      package_type = "sqlite",
      package_bytes = package_size,
      database_files = 1L,
      gff_files = 0L,
      total_entries = 1L,
      portable = FALSE,
      warning = paste(
        "The selected file contains only the SQLite index.",
        "GFF annotation files cannot be imported from this upload."
      )
    ))
  }

  if (!identical(extension, "zip")) {
    stop("Select a .zip portable cache or a .sqlite recovery file.", call. = FALSE)
  }

  listing <- utils::unzip(package_file, list = TRUE)
  names <- chartr("\\", "/", listing$Name)
  is_directory <- grepl("/$", names)
  file_names <- names[!is_directory]

  database_match <- grepl("(^|/)(GenomeCache|ReferenceGenomeCache)\\.sqlite$", file_names, ignore.case = TRUE)
  gff_match <- grepl("(^|/)GFF/", file_names, ignore.case = TRUE)

  database_files <- sum(database_match)
  gff_files <- sum(gff_match)
  portable <- database_files == 1L && gff_files > 0L

  warning <- NULL
  if (database_files == 0L) warning <- "No GenomeCache.sqlite database was found in the ZIP."
  if (database_files > 1L) warning <- "More than one cache database was found in the ZIP."
  if (database_files == 1L && gff_files == 0L) {
    warning <- "The ZIP contains a database but no GFF/ annotation files."
  }

  list(
    package_type = "zip",
    package_bytes = package_size,
    database_files = as.integer(database_files),
    gff_files = as.integer(gff_files),
    total_entries = as.integer(length(file_names)),
    portable = portable,
    warning = warning
  )
}

aaa_cache_import_package <- function(target_database, package_file, conflict = c("keep_existing", "replace_newer")) {
  conflict <- match.arg(conflict)
  if (!file.exists(package_file)) stop("Cache package not found.", call. = FALSE)
  stage <- tempfile("TripleA_cache_import_")
  dir.create(stage, recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  if (grepl("\\.zip$", package_file, ignore.case = TRUE)) {
    utils::unzip(package_file, exdir = stage)
    candidates <- list.files(stage, pattern = "GenomeCache\\.sqlite$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (!length(candidates)) stop("The ZIP does not contain GenomeCache.sqlite.", call. = FALSE)
    source_db <- candidates[1L]
    source_root <- dirname(source_db)
  } else {
    source_db <- package_file
    source_root <- dirname(package_file)
  }
  result <- aaa_cache_merge(target_database, source_db, source_root, conflict = conflict, create_backup = TRUE)
  result$source_format <- if (grepl("\\.zip$", package_file, ignore.case = TRUE)) "portable_zip" else "sqlite_only"
  result$warning <- if (identical(result$source_format, "sqlite_only")) "Only the SQLite index was supplied. GFF annotation files were copied only when a companion GFF folder was present beside it." else NULL
  result
}
