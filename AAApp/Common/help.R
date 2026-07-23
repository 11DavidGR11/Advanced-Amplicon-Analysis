# Shared contextual documentation utilities.

aaa_help_title <- function(path) {
  lines <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  heading <- grep("^\\s*#\\s+.+", lines, value = TRUE)
  if (length(heading)) {
    title <- trimws(sub("^\\s*#\\s+", "", heading[[1]]))
    if (nzchar(title)) return(title)
  }
  trimws(gsub("[_-]+", " ", tools::file_path_sans_ext(basename(path))))
}

aaa_discover_context_help <- function(project_root, component) {
  folder <- file.path(project_root, "Resources", "Documentation", component)
  if (!dir.exists(folder)) return(setNames(character(), character()))
  files <- list.files(folder, pattern = "\\.[Mm][Dd]$", full.names = TRUE, recursive = TRUE)
  files <- files[!file.info(files)$isdir]
  if (!length(files)) return(setNames(character(), character()))
  titles <- vapply(files, aaa_help_title, character(1))
  idx <- order(tolower(titles), tolower(files))
  stats::setNames(files[idx], make.unique(titles[idx], sep = " — "))
}

aaa_help_button <- function(id = "context_help", label = "Help") {
  shiny::actionButton(id, label, icon = shiny::icon("circle-question"), class = "btn-outline-secondary w-100")
}

aaa_register_context_help <- function(input, output, session, project_root, component, id = "context_help") {
  shiny::observeEvent(input[[id]], ignoreInit = TRUE, {
    topics <- aaa_discover_context_help(project_root, component)
    if (!length(topics)) {
      shiny::showModal(shiny::modalDialog(
        title = "Documentation unavailable",
        shiny::tags$p("No Markdown documentation was found for this application."),
        shiny::tags$code(file.path("Resources", "Documentation", component)),
        easyClose = TRUE, footer = shiny::modalButton("Close")
      ))
      return()
    }
    ns <- session$ns
    selector_id <- ns(paste0(id, "_topic"))
    content_id <- ns(paste0(id, "_content"))
    shiny::showModal(shiny::modalDialog(
      title = "Application help",
      shiny::selectInput(selector_id, "Topic", choices = stats::setNames(unname(topics), names(topics))),
      shiny::uiOutput(content_id),
      size = "l", easyClose = TRUE, footer = shiny::modalButton("Close")
    ))
    output[[paste0(id, "_content")]] <- shiny::renderUI({
      selected <- input[[paste0(id, "_topic")]]
      shiny::req(selected, selected %in% unname(topics), file.exists(selected))
      shiny::includeMarkdown(selected)
    })
  })
  invisible(TRUE)
}
