library(shiny)
library(bslib)
library(DT)

# The cache manager is a local application and may need to receive large portable
# packages. The limit can be changed before launching with:
# options(triple_a_max_upload_bytes = <bytes>)
max_upload_bytes <- getOption("triple_a_max_upload_bytes", 10 * 1024^3)
options(shiny.maxRequestSize = max_upload_bytes)

app_dir <- normalizePath(getOption("triple_a_cache_manager_dir", getwd()), winslash = "/", mustWork = TRUE)
source(file.path(dirname(app_dir), "Common", "paths.R"), local = TRUE)
root <- normalizePath(getOption("triple_a_root", aaa_distribution_root(app_dir)), winslash = "/", mustWork = TRUE)
source(file.path(root, "AAApp", "Common", "help.R"), local = TRUE)
source(file.path(root, "AAApp", "Common", "Engine", "Core", "aaa_cache_database.R"), local = TRUE)
source(file.path(root, "AAApp", "Common", "Engine", "Core", "aaa_genome_cache_manager.R"), local = TRUE)
database_file <- aaa_genome_cache_db()
dir.create(aaa_cached_gff_path(), recursive = TRUE, showWarnings = FALSE)
dir.create(aaa_cache_backups_path(), recursive = TRUE, showWarnings = FALSE)

format_bytes <- function(x) {
  if (length(x) == 0L || !is.finite(x) || is.na(x)) return("0 B")
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- 1L
  x <- as.numeric(x)
  while (x >= 1024 && i < length(units)) {
    x <- x / 1024
    i <- i + 1L
  }
  sprintf(if (i == 1L) "%.0f %s" else "%.2f %s", x, units[[i]])
}

stat_card <- function(label, value, detail, icon) {
  tags$div(
    class = "cache-stat-card",
    tags$div(class = "cache-stat-icon", icon),
    tags$div(
      class = "cache-stat-content",
      tags$div(class = "cache-stat-label", label),
      tags$div(class = "cache-stat-number", value),
      tags$div(class = "cache-stat-detail", detail)
    )
  )
}

ui <- page_navbar(
  title = "Cache (Triple A)",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#5B3F7A"),
  header = tags$head(
    tags$style(HTML("\
      .cache-summary-grid {margin-bottom: .85rem;}\
      .cache-stat-card {\
        min-height: 108px; padding: .9rem 1rem; border: 1px solid #d9dee5;\
        border-radius: .75rem; background: #fff; box-shadow: 0 2px 8px rgba(0,0,0,.07);\
        display: flex; gap: .8rem; align-items: center; overflow: hidden;\
      }\
      .cache-stat-icon {font-size: 1.65rem; line-height: 1; flex: 0 0 auto;}\
      .cache-stat-content {min-width: 0; flex: 1 1 auto;}\
      .cache-stat-label {font-size: .78rem; font-weight: 700; color: #52606d; text-transform: uppercase; letter-spacing: .03em;}\
      .cache-stat-number {font-size: clamp(1.55rem, 2.6vw, 2.15rem); font-weight: 800; line-height: 1.08; color: #243447; overflow-wrap: anywhere;}\
      .cache-stat-detail {font-size: .82rem; color: #687786; line-height: 1.2; white-space: normal;}\
      .cache-path-block {padding: .8rem 1rem; border-radius: .5rem; background: #f7f8fa; border: 1px solid #e1e5ea;}\
      .cache-path {display: block; margin-top: .25rem; overflow-wrap: anywhere; word-break: break-word; font-family: var(--bs-font-monospace); font-size: .86rem;}\
      .cache-table-wrap {overflow-x: auto; overflow-y: auto; max-height: calc(100vh - 350px); min-height: 280px;}\
      .cache-compact-card .card-body {padding: 1rem 1.15rem;}\
      .cache-actions {display: flex; flex-wrap: wrap; gap: .75rem; align-items: center;}\
      .cache-actions .btn {width: auto; min-width: 180px;}\
      .cache-section-divider {border: 0; border-top: 1px solid #d8dde3; margin: 1.1rem 0;}\
      .package-preview {border-left: 5px solid #5B3F7A;}\
      .package-metric {font-size: .98rem; margin-bottom: .3rem;}\
      .upload-limit {font-size: .88rem; color: #52606d; margin-bottom: .5rem;}\
      .cache-note {font-size: .84rem; color: #687786;}\
      @media (max-width: 768px) {\
        .cache-stat-card {min-height: 96px;}\
        .cache-stat-number {font-size: 1.6rem;}\
        .cache-table-wrap {max-height: none; min-height: 240px;}\
        .cache-actions .btn {width: 100%;}\
      }\
    "))
  ),
  nav_panel(
    "Overview",
    aaa_help_button(),
    uiOutput("cache_summary"),
    card(
      class = "cache-compact-card",
      card_header("Active cache"),
      tags$div(
        class = "cache-path-block",
        tags$strong("Active database"),
        tags$span(class = "cache-path", database_file)
      ),
      tags$p(class = "cache-note mt-2 mb-2", "Assembly accessions are indexed in GenomeCache.sqlite. GFF annotations are downloaded lazily and stored once in Cache/GFF/."),
      tags$div(class = "cache-actions", actionButton("refresh_overview", "Refresh overview", class = "btn-outline-primary"))
    ),
    card(
      class = "cache-compact-card",
      card_header("Stored references"),
      tags$div(class = "cache-table-wrap", DTOutput("references"))
    )
  ),
  nav_panel(
    "Import / Merge",
    card(
      class = "cache-compact-card",
      card_header("Import a portable cache"),
      tags$p("Select a Triple_A portable ZIP. It must include the SQLite index and, when present, the GFF/ folder annotations."),
      tags$div(
        class = "alert alert-warning",
        "A standalone GenomeCache.sqlite file lets you recover the index, but does not carry the GFF annotations. Using the portable ZIP is recommended."
      ),
      fileInput(
        "cache_file",
        "Cache package",
        accept = c(".zip", ".sqlite"),
        width = "100%",
        buttonLabel = "Browse...",
        placeholder = "No file selected"
      ),
      tags$p(class = "upload-limit", paste("Maximum upload size:", format_bytes(max_upload_bytes))),
      uiOutput("package_preview"),
      tags$hr(class = "cache-section-divider"),
      radioButtons(
        "conflict",
        "Conflict policy",
        choices = c("Keep existing records" = "keep_existing", "Use the newest records" = "replace_newer"),
        inline = TRUE
      ),
      tags$div(class = "cache-actions", actionButton("merge", "Import and merge", class = "btn-primary")),
      tags$hr(class = "cache-section-divider"),
      verbatimTextOutput("merge_report")
    )
  ),
  nav_panel(
    "Export",
    card(
      class = "cache-compact-card",
      card_header("Create a portable cache"),
      tags$p("The ZIP contains GenomeCache.sqlite, GFF/, and a manifest. It can be imported on another machine without manually copying files."),
      tags$div(class = "cache-actions", downloadButton("export_cache", "Export portable cache", class = "btn-primary"))
    )
  ),
  nav_panel(
    "Integrity",
    layout_columns(
      col_widths = c(4, 4, 4),
      card(class = "cache-compact-card", card_header("Verification"), actionButton("verify", "Verify integrity", class = "btn-primary"), verbatimTextOutput("verify_report")),
      card(
        class = "cache-compact-card",
        card_header("Orphan file cleanup"),
        tags$p("Only files not referenced by the SQLite index are considered candidates."),
        tags$div(class = "cache-actions", actionButton("preview_clean", "Preview cleanup"), actionButton("clean", "Delete orphan files", class = "btn-danger")),
        verbatimTextOutput("clean_report")
      ),
      card(
        class = "cache-compact-card",
        card_header("Stale reference cleanup"),
        tags$p("A taxon's NCBI reference genome is cached indefinitely. Remove entries not refreshed in a while to force them to be re-resolved on the next run."),
        numericInput("stale_days", "Older than (days)", value = 180, min = 1, step = 1),
        tags$div(class = "cache-actions", actionButton("preview_stale", "Preview cleanup"), actionButton("clean_stale", "Delete stale references", class = "btn-danger")),
        verbatimTextOutput("stale_report")
      )
    )
  ),
  nav_panel("History", card(class = "cache-compact-card", card_header("Cache operations"), tags$div(class = "cache-table-wrap", DTOutput("history"))))
)

server <- function(input, output, session) {
  aaa_register_context_help(input, output, session, root, "CacheManager")
  refresh <- reactiveVal(0L)
  merge_result <- reactiveVal(NULL)
  verify_result <- reactiveVal(NULL)
  clean_result <- reactiveVal(NULL)
  stale_result <- reactiveVal(NULL)
  package_preview <- reactiveVal(NULL)

  stats <- reactive({
    refresh()
    aaa_cache_statistics(database_file)
  })

  output$cache_summary <- renderUI({
    x <- stats()
    layout_columns(
      class = "cache-summary-grid",
      col_widths = c(3, 3, 3, 3),
      stat_card("References", format(x$references, big.mark = ","), paste(format(x$genera, big.mark = ","), "genera"), "🧬"),
      stat_card("GFF annotations", format(x$gff_files, big.mark = ","), paste("Local files ·", format_bytes(x$gff_bytes)), "📄"),
      stat_card("Gene records", format(x$gene_results, big.mark = ","), paste(format(x$genes, big.mark = ","), "genes evaluated"), "🧬"),
      stat_card("Total size", format_bytes(x$bytes), paste("SQLite:", format_bytes(x$database_bytes)), "💾")
    )
  })

  output$references <- renderDT({
    refresh()
    DT::datatable(
      head(aaa_cache_reference_table(database_file), 5000L),
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE, filter = "top"
    )
  })

  output$history <- renderDT({
    refresh()
    DT::datatable(
      aaa_cache_operation_log(database_file, 500L),
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE
    )
  })

  observeEvent(input$refresh_overview, {
    refresh(refresh() + 1L)
  })

  observeEvent(input$cache_file, {
    req(input$cache_file$datapath)
    preview <- tryCatch(
      aaa_cache_package_preview(input$cache_file$datapath),
      error = function(e) e
    )
    package_preview(preview)
  }, ignoreInit = TRUE)

  output$package_preview <- renderUI({
    x <- package_preview()
    if (is.null(x)) return(NULL)
    if (inherits(x, "error")) {
      return(tags$div(class = "alert alert-danger", conditionMessage(x)))
    }

    warning_ui <- if (!is.null(x$warning) && nzchar(x$warning)) {
      tags$div(class = "alert alert-warning mt-3", x$warning)
    }

    card(
      class = "package-preview mt-3 mb-3",
      card_header("Detected content"),
      tags$div(class = "package-metric", tags$strong("File: "), input$cache_file$name),
      tags$div(class = "package-metric", tags$strong("Size: "), format_bytes(x$package_bytes)),
      tags$div(class = "package-metric", tags$strong("SQLite databases: "), x$database_files),
      tags$div(class = "package-metric", tags$strong("GFF: "), x$gff_files),
      warning_ui
    )
  })

  observeEvent(input$merge, {
    req(input$cache_file$datapath)
    preview <- package_preview()
    if (inherits(preview, "error")) {
      showNotification("The package could not be processed.", type = "error")
      return()
    }
    if (!is.null(preview$warning)) {
      showNotification(preview$warning, type = "warning", duration = 8)
    }

    result <- withProgress(message = "Importing the cache", value = 0, {
      incProgress(0.10, detail = "Preparing the package...")
      uploaded <- tempfile(fileext = paste0("_", basename(input$cache_file$name)))
      copied <- file.copy(input$cache_file$datapath, uploaded, overwrite = TRUE)
      if (!isTRUE(copied)) stop("The selected file could not be prepared.", call. = FALSE)
      on.exit(unlink(uploaded, force = TRUE), add = TRUE)

      incProgress(0.25, detail = "Validating the content...")
      aaa_cache_package_preview(uploaded)
      incProgress(0.20, detail = "Merging records and files...")
      answer <- tryCatch(
        aaa_cache_import_package(database_file, uploaded, input$conflict),
        error = function(e) e
      )
      incProgress(0.45, detail = "Updating the overview...")
      answer
    })

    merge_result(result)
    if (!inherits(result, "error")) {
      refresh(refresh() + 1L)
      showNotification("The cache was imported successfully.", type = "message")
    }
  })

  output$merge_report <- renderPrint({
    x <- merge_result()
    if (is.null(x)) {
      cat("No import has been performed yet in this session.\n")
    } else if (inherits(x, "error")) {
      cat("ERROR:", conditionMessage(x), "\n")
    } else {
      cat("Import completed\n")
      cat("References added:", x$references_added, "\n")
      cat("References updated:", x$references_updated, "\n")
      cat("GFF files added:", x$gff_added, "\n")
      cat("Duplicates ignored:", x$duplicate_files_ignored, "\n")
      cat("Files replaced:", x$files_replaced, "\n")
      cat("Conflicts:", length(x$file_conflicts), "\n")
      cat("Backup:", x$backup, "\n")
      if (!is.null(x$warning)) cat("WARNING:", x$warning, "\n")
    }
  })

  output$export_cache <- downloadHandler(
    filename = function() paste0("TripleA_ReferenceAnnotationCache_", format(Sys.Date(), "%Y-%m-%d"), ".zip"),
    content = function(file) aaa_cache_export(database_file, file),
    contentType = "application/zip"
  )

  observeEvent(input$verify, {
    verify_result(tryCatch(aaa_cache_verify(database_file), error = function(e) e))
    refresh(refresh() + 1L)
  })

  output$verify_report <- renderPrint({
    x <- verify_result()
    if (is.null(x)) cat("Click Verify integrity.\n")
    else if (inherits(x, "error")) cat("ERROR:", conditionMessage(x), "\n")
    else {
      cat("SQLite:", paste(x$sqlite_integrity, collapse = ", "), "\n")
      cat("References:", x$references, "\n")
      cat("Local coverage:", x$local_coverage, "\n")
      cat("Missing referenced files:", nrow(x$missing_files), "\n")
      cat("Orphan files:", nrow(x$orphan_files), "\n")
      cat("Status:", if (x$valid) "VALID" else "NEEDS REVIEW", "\n")
    }
  })

  observeEvent(input$preview_clean, {
    clean_result(tryCatch(aaa_cache_clean_orphans(database_file, TRUE), error = function(e) e))
  })

  observeEvent(input$clean, {
    showModal(modalDialog(
      title = "Delete orphan files?",
      tags$p("A backup of the database will be created. Only files not referenced by the index will be deleted."),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_clean", "Delete", class = "btn-danger"))
    ))
  })

  observeEvent(input$confirm_clean, {
    removeModal()
    clean_result(tryCatch(aaa_cache_clean_orphans(database_file, FALSE), error = function(e) e))
    refresh(refresh() + 1L)
  })

  output$clean_report <- renderPrint({
    x <- clean_result()
    if (is.null(x)) cat("Preview the cleanup before deleting files.\n")
    else if (inherits(x, "error")) cat("ERROR:", conditionMessage(x), "\n")
    else print(x)
  })

  observeEvent(input$preview_stale, {
    req(input$stale_days)
    stale_result(tryCatch(aaa_cache_prune_stale_references(database_file, input$stale_days, TRUE), error = function(e) e))
  })

  observeEvent(input$clean_stale, {
    req(input$stale_days)
    showModal(modalDialog(
      title = "Delete stale references?",
      tags$p("A backup of the database will be created. Removed taxa keep their downloaded GFF files and will be re-resolved from NCBI the next time they are analysed."),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_clean_stale", "Delete", class = "btn-danger"))
    ))
  })

  observeEvent(input$confirm_clean_stale, {
    removeModal()
    stale_result(tryCatch(aaa_cache_prune_stale_references(database_file, input$stale_days, FALSE), error = function(e) e))
    refresh(refresh() + 1L)
  })

  output$stale_report <- renderPrint({
    x <- stale_result()
    if (is.null(x)) cat("Preview the cleanup before deleting references.\n")
    else if (inherits(x, "error")) cat("ERROR:", conditionMessage(x), "\n")
    else {
      cat(if (x$dry_run) "Preview only, nothing deleted.\n" else "Cleanup completed.\n")
      cat("Older than:", x$max_age_days, "day(s)\n")
      cat("Candidates:", x$candidates, "\n")
      cat("Removed:", x$removed, "\n")
    }
  })
}

shinyApp(ui, server)
