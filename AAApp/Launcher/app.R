library(shiny)
library(bslib)

launcher_dir <- normalizePath(getOption("triple_a_launcher_dir", getwd()), winslash = "/", mustWork = TRUE)
paths_file <- file.path(dirname(launcher_dir), "Common", "paths.R")
if (!file.exists(paths_file)) stop("AAApp/Common/paths.R was not found.", call. = FALSE)
source(paths_file, local = TRUE)
root <- normalizePath(getOption("triple_a_root", aaa_distribution_root(launcher_dir)), winslash = "/", mustWork = TRUE)
source(file.path(root, "AAApp", "Common", "help.R"), local = TRUE)

required_launcher_packages <- c("shiny", "bslib", "callr", "httpuv")
missing_launcher_packages <- required_launcher_packages[
  !vapply(required_launcher_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_launcher_packages) > 0L) {
  stop(
    "The Triple_A launcher requires: ",
    paste(missing_launcher_packages, collapse = ", "),
    ".",
    call. = FALSE
  )
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

# Presence check that does NOT load the package namespace. requireNamespace()
# would load dada2 (a heavy Bioconductor package) into the single-threaded
# launcher process, freezing its UI for several seconds and making the FASTQ
# button feel unresponsive. system.file() only resolves the install path.
package_installed <- function(pkg) nzchar(system.file(package = pkg))

new_local_port <- function() {
  as.integer(httpuv::randomPort())
}

wait_for_local_app <- function(process, port, timeout = 20) {
  deadline <- Sys.time() + timeout
  repeat {
    if (!process$is_alive()) return(FALSE)
    connection <- suppressWarnings(tryCatch(
      socketConnection(
        host = "127.0.0.1",
        port = port,
        open = "r+",
        blocking = TRUE,
        timeout = 1
      ),
      error = function(e) NULL
    ))
    if (!is.null(connection)) {
      try(close(connection), silent = TRUE)
      return(TRUE)
    }
    if (Sys.time() >= deadline) return(FALSE)
    Sys.sleep(0.25)
  }
}

child_app_function <- function(app_dir, option_name, port, log_file, max_upload_bytes) {
  options_list <- list(
    shiny.maxRequestSize = max_upload_bytes,
    triple_a_max_upload_bytes = max_upload_bytes
  )
  options_list[[option_name]] <- app_dir
  do.call(options, options_list)

  cat(sprintf("[%s] Starting app in %s on port %s\n", Sys.time(), app_dir, port))
  shiny::runApp(
    appDir = app_dir,
    host = "127.0.0.1",
    port = port,
    launch.browser = FALSE,
    display.mode = "normal"
  )
}

launcher_card <- function(title, description, button_id, button_label, icon_name = "circle") {
  tags$article(
    class = "launcher-card",
    tags$div(
      class = "launcher-card-icon",
      icon(icon_name)
    ),
    tags$div(
      class = "launcher-card-body",
      tags$h3(class = "launcher-card-title", title),
      tags$p(class = "launcher-card-description", description),
      actionButton(button_id, button_label, class = "btn-primary launcher-card-button")
    )
  )
}

ui <- page_fluid(
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#5B3F7A"),
  tags$head(
    tags$title("Advanced Amplicon Analysis (Triple A)"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML(
      ".launcher-shell{max-width:1280px;margin:0 auto;padding:2rem 1rem 3rem;}\n       .launcher-header{margin-bottom:2rem;text-align:center;}\n       .launcher-brand{display:flex;flex-direction:column;align-items:center;}\n       .launcher-logo{display:block;width:min(100%,520px);height:auto;margin:0 auto 1rem;}\n       .launcher-header h1{font-weight:700;margin-bottom:.35rem;}\n       .launcher-introduction{max-width:960px;margin:0 auto;text-align:left;}
       .launcher-summary{max-width:900px;margin:1rem auto 0;line-height:1.6;}
       .launcher-summary ul{margin-bottom:0;}\n       .launcher-section{margin-top:2rem;}\n       .launcher-section-title{font-size:1.45rem;font-weight:650;margin-bottom:1rem;}\n       .launcher-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1rem;}\n       .launcher-grid-single{grid-template-columns:minmax(0,1fr);}\n       .launcher-card{display:flex;align-items:flex-start;gap:1rem;min-width:0;height:100%;padding:1.25rem;border:1px solid var(--bs-border-color);border-radius:.75rem;background:var(--bs-body-bg);box-shadow:0 .15rem .5rem rgba(0,0,0,.06);}\n       .launcher-card-icon{display:flex;align-items:center;justify-content:center;flex:0 0 3rem;height:3rem;border-radius:.65rem;background:rgba(91,63,122,.12);color:#5B3F7A;font-size:1.35rem;}\n       .launcher-card-body{display:flex;flex-direction:column;align-items:flex-start;min-width:0;height:100%;}\n       .launcher-card-title{font-size:1.15rem;font-weight:650;margin:0 0 .4rem;}\n       .launcher-card-description{margin:0 0 1rem;overflow-wrap:anywhere;}\n       .launcher-card-button{margin-top:auto;white-space:normal;}\n       .launcher-actions{display:flex;gap:.6rem;flex-wrap:wrap;}\n       .workflow-figure{width:100%;height:auto;border:1px solid var(--bs-border-color);border-radius:.5rem;}\n       .launcher-status{margin-top:2rem;}\n       @media (max-width:760px){.launcher-grid{grid-template-columns:minmax(0,1fr)}.launcher-shell{padding-top:1.25rem}.launcher-card{padding:1rem}}"
    )),
    tags$script(HTML(
      "Shiny.addCustomMessageHandler('tripleAOpenUrl', function(message) {\n         var win = window.open(message.url, '_blank');\n         if (win) { win.focus(); }\n       });"
    ))
  ),
  tags$main(
    class = "launcher-shell",
    tags$header(
      class = "launcher-header",
      tags$div(
        class = "launcher-brand",
        tags$img(
          src = "triple_a_logo.png",
          class = "launcher-logo",
          alt = "Triple A — Advanced Amplicon Analysis Platform"
        ),
        tags$h1("Triple A", class = "visually-hidden")
      ),
      tags$div(
        class = "launcher-introduction",
        tags$p(class = "lead", HTML("<strong>Triple A</strong> is a comprehensive platform for the analysis of amplicon sequencing data, integrating taxonomic profiling, functional prediction, biological-function inference and ecological statistics within a single reproducible workflow.")),
        tags$p("Choose the workflow corresponding to your starting data."),
        tags$div(
          class = "launcher-actions",
          tags$div(style = "max-width:220px;", aaa_help_button()),
          tags$div(
            style = "max-width:220px;",
            actionButton(
              "show_workflow", "Workflow",
              icon = icon("diagram-project"),
              class = "btn-outline-secondary w-100"
            )
          )
        )
      ),
    ),
    tags$section(
      class = "launcher-section",
      tags$h2(class = "launcher-section-title", "Analysis"),
      tags$div(
        class = "launcher-grid launcher-grid-single",
        launcher_card(
          title = "Biological Analysis (Triple A)",
          description = "Import abundance data and run taxonomic, statistical and functional analyses.",
          button_id = "open_analysis",
          button_label = "Open Biological Analysis",
          icon_name = "chart-bar"
        )
      )
    ),
    tags$section(
      class = "launcher-section",
      tags$h2(class = "launcher-section-title", "Data preparation"),
      tags$div(
        class = "launcher-grid",
        launcher_card(
          title = "Amplicon Integrator (Triple A)",
          description = "Validate and combine compatible amplicon count tables.",
          button_id = "open_multiamplicon",
          button_label = "Open MultiAmplicon",
          icon_name = "layer-group"
        ),
        launcher_card(
          title = "FASTQ Pipeline (Triple A)",
          description = "Run DADA2 preprocessing and export a compatible abundance table.",
          button_id = "open_fastq",
          button_label = "Open FASTQ Pipeline",
          icon_name = "dna"
        )
      )
    ),
    tags$section(
      class = "launcher-section",
      tags$h2(class = "launcher-section-title", "Resources and extension"),
      tags$div(
        class = "launcher-grid",
        launcher_card(
          title = "Cache (Triple A)",
          description = "Inspect, merge, export and verify cached reference annotations.",
          button_id = "open_cache_manager",
          button_label = "Open Cache Manager",
          icon_name = "database"
        ),
        launcher_card(
          title = "Biological Function Builder (Triple A)",
          description = "Create and validate declarative biological functions without editing R scripts.",
          button_id = "open_function_builder",
          button_label = "Open Function Builder",
          icon_name = "puzzle-piece"
        )
      )
    ),
    tags$div(class = "launcher-status", uiOutput("launcher_status")),
    tags$p(
      class = "launcher-note",
      HTML(
        "Each workflow runs in an independent R process. Startup logs are saved in <code>Results/Logs</code>.<br><br>
    <strong>Triple A</strong><br>
    Developed by David Garrido Rodríguez"
      )
    )
  )
)

server <- function(input, output, session) {
  aaa_register_context_help(input, output, session, root, "General")

  # The workflow diagram is a static asset in www/, so it is addressed by name
  # and opens full size in a new tab for reading the small print.
  workflow_figure <- "Triple_A_complete_workflow.png"

  observeEvent(input$show_workflow, ignoreInit = TRUE, {
    # Resolved from launcher_dir rather than the working directory, which is not
    # guaranteed to be the app folder.
    if (!file.exists(file.path(launcher_dir, "www", workflow_figure))) {
      showModal(modalDialog(
        title = "Workflow diagram unavailable",
        tags$p("The workflow figure was not found in the launcher's www folder."),
        tags$code(file.path("AAApp", "Launcher", "www", workflow_figure)),
        easyClose = TRUE, footer = modalButton("Close")
      ))
      return(invisible(NULL))
    }
    showModal(modalDialog(
      title = "Triple A — complete project workflow",
      tags$p(
        class = "text-muted mb-3",
        "From raw sequencing reads to reproducible, documented biological results. ",
        "Open it full size to read the detail."
      ),
      tags$img(src = workflow_figure, class = "workflow-figure",
               alt = "Complete Triple A project workflow"),
      easyClose = TRUE,
      size = "xl",
      footer = tagList(
        tags$a(
          href = workflow_figure, target = "_blank", class = "btn btn-outline-secondary",
          icon("up-right-from-square"), " Open full size"
        ),
        modalButton("Close")
      )
    ))
  })

  child_processes <- reactiveVal(list())
  status_text <- reactiveVal("Launcher ready.")
  # Persists the most recently started module so the launcher can always show a
  # visible, clickable link. The automatic window.open() in show_launch_result()
  # is fired from an asynchronous server message rather than the click gesture,
  # so browsers routinely block it as a pop-up; without a persistent link the
  # user is left staring at the launcher with the module running but unreachable.
  active_link <- reactiveVal(NULL)
  pending_fastq_launch <- reactiveVal(FALSE)

  output$launcher_status <- renderUI({
    link <- active_link()
    tagList(
      tags$div(class = "alert alert-light py-2", status_text()),
      if (!is.null(link)) tags$div(
        class = "alert alert-success py-2 d-flex align-items-center flex-wrap gap-2",
        tags$span(paste0(link$app_name, " is running.")),
        tags$a(
          href = link$url, target = "_blank", rel = "noopener noreferrer",
          class = "btn btn-sm btn-primary",
          icon("up-right-from-square"), paste("Open", link$app_name)
        ),
        tags$span(class = "text-muted small", link$url)
      )
    )
  })

  register_process <- function(process, app_name) {
    processes <- child_processes()
    processes[[paste0(safe_name(app_name), "_", as.integer(Sys.time()))]] <- process
    child_processes(processes)
  }

  show_launch_result <- function(app_name, url, log_file) {
    status_text(paste(app_name, "started at", url))
    active_link(list(app_name = app_name, url = url))
    session$sendCustomMessage("tripleAOpenUrl", list(url = url))
    showModal(modalDialog(
      title = paste(app_name, "started"),
      tags$p("The application is running in a separate local R process."),
      tags$p(
        "If the new browser tab did not open automatically, use this link ",
        "(it also stays available in the launcher status panel below): ",
        tags$a(
          href = url, target = "_blank", rel = "noopener noreferrer",
          class = "btn btn-sm btn-primary",
          paste("Open", app_name)
        )
      ),
      tags$p(class = "text-muted small", paste("Startup log:", log_file)),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  }

  launch_child_app <- function(app_dir, option_name, app_name, log_name) {
    app_dir <- normalizePath(app_dir, winslash = "/", mustWork = FALSE)
    app_file <- file.path(app_dir, "app.R")

    if (!file.exists(app_file)) {
      showModal(modalDialog(
        title = paste("Unable to open", app_name),
        tags$p("The application file was not found."),
        tags$pre(app_file),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return(invisible(FALSE))
    }

    port <- new_local_port()
    url <- sprintf("http://127.0.0.1:%d", port)
    log_dir <- file.path(root, "Results", "Logs")
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- normalizePath(
      file.path(log_dir, paste0(log_name, ".log")),
      winslash = "/",
      mustWork = FALSE
    )

    cat(sprintf("\n[%s] Launcher requested %s\n", Sys.time(), app_name),
        file = log_file, append = TRUE)

    status_text(paste("Starting", app_name, "..."))

    process <- tryCatch(
      callr::r_bg(
        func = child_app_function,
        args = list(
          app_dir = app_dir,
          option_name = option_name,
          port = port,
          log_file = log_file,
          max_upload_bytes = getOption("triple_a_max_upload_bytes", 10 * 1024^3)
        ),
        supervise = FALSE,
        stdout = log_file,
        stderr = log_file,
        system_profile = TRUE,
        user_profile = TRUE
      ),
      error = function(e) e
    )

    if (inherits(process, "error")) {
      status_text(paste("Could not start", app_name, "."))
      showModal(modalDialog(
        title = paste("Unable to open", app_name),
        tags$p(conditionMessage(process)),
        tags$p("Diagnostic log:"),
        tags$pre(log_file),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return(invisible(FALSE))
    }

    register_process(process, app_name)

    if (!wait_for_local_app(process, port, timeout = 20)) {
      status_text(paste(app_name, "did not finish starting."))
      log_tail <- if (file.exists(log_file)) {
        paste(tail(readLines(log_file, warn = FALSE), 30), collapse = "\n")
      } else {
        "No startup log was created."
      }
      showModal(modalDialog(
        title = paste("Unable to open", app_name),
        tags$p("The child R process stopped or did not open its local port."),
        tags$p("Startup diagnostics:"),
        tags$pre(style = "max-height: 320px; overflow-y: auto;", log_tail),
        tags$p(class = "text-muted small", paste("Full log:", log_file)),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return(invisible(FALSE))
    }

    show_launch_result(app_name, url, log_file)
    invisible(TRUE)
  }

  launch_analysis <- function() {
    launch_child_app(
      app_dir = file.path(root, "AAApp", "Biological"),
      option_name = "triple_a_app_dir",
      app_name = "Biological Analysis (Triple A)",
      log_name = "analysis_startup"
    )
  }

  launch_multiamplicon <- function() {
    launch_child_app(
      app_dir = file.path(root, "AAApp", "MultiAmplicon"),
      option_name = "triple_a_multiamplicon_dir",
      app_name = "Amplicon Integrator (Triple A)",
      log_name = "multiamplicon_startup"
    )
  }

  launch_fastq <- function() {
    launch_child_app(
      app_dir = file.path(root, "AAApp", "FASTQ"),
      option_name = "triple_a_fastq_app_dir",
      app_name = "FASTQ Pipeline (Triple A)",
      log_name = "fastq_startup"
    )
  }

  launch_cache_manager <- function() {
    launch_child_app(
      app_dir = file.path(root, "AAApp", "CacheManager"),
      option_name = "triple_a_cache_manager_dir",
      app_name = "Cache (Triple A)",
      log_name = "cache_manager_startup"
    )
  }

  launch_function_builder <- function() {
    launch_child_app(
      app_dir = file.path(root, "AAApp", "FunctionBuilder"),
      option_name = "triple_a_function_builder_dir",
      app_name = "Biological Function Builder (Triple A)",
      log_name = "function_builder_startup"
    )
  }


  ask_to_install_dada2 <- function() {
    showModal(modalDialog(
      title = "Install the FASTQ dependency?",
      tags$p("The FASTQ workflow requires the Bioconductor package DADA2, which is not installed."),
      tags$p("Triple_A can install DADA2 and its required Bioconductor dependencies now."),
      tags$div(
        class = "alert alert-warning",
        "The download is relatively large and the installation may require an internet connection and permission to write to your R library."
      ),
      easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("install_dada2", "Install DADA2", class = "btn-primary")
      )
    ))
  }

  observeEvent(input$open_analysis, ignoreInit = TRUE, {
    launch_analysis()
  })

  observeEvent(input$open_multiamplicon, ignoreInit = TRUE, {
    launch_multiamplicon()
  })

  observeEvent(input$open_cache_manager, ignoreInit = TRUE, {
    launch_cache_manager()
  })

  observeEvent(input$open_function_builder, ignoreInit = TRUE, {
    launch_function_builder()
  })


  observeEvent(input$open_fastq, ignoreInit = TRUE, {
    missing_fastq_support <- c("jsonlite")[
      !vapply(c("jsonlite"), requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_fastq_support) > 0L) {
      showModal(modalDialog(
        title = "FASTQ support package missing",
        tags$p("Install the following CRAN package before opening the FASTQ module:"),
        tags$pre(paste(missing_fastq_support, collapse = ", ")),
        tags$pre(sprintf('install.packages(c(%s))', paste(sprintf('"%s"', missing_fastq_support), collapse = ", "))),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return()
    }

    if (package_installed("dada2")) {
      launch_fastq()
    } else {
      pending_fastq_launch(TRUE)
      ask_to_install_dada2()
    }
  })

  observeEvent(input$install_dada2, ignoreInit = TRUE, {
    removeModal()
    status_text("Installing DADA2 and its Bioconductor dependencies...")

    installation_error <- NULL
    withProgress(message = "Installing DADA2", value = 0, {
      incProgress(0.1, detail = "Preparing Bioconductor")
      tryCatch({
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        incProgress(0.25, detail = "Downloading and installing packages")
        BiocManager::install("dada2", ask = FALSE, update = FALSE)
        incProgress(0.6, detail = "Verifying installation")
      }, error = function(e) {
        installation_error <<- e
      })
    })

    if (!is.null(installation_error) || !requireNamespace("dada2", quietly = TRUE)) {
      status_text("DADA2 installation did not complete.")
      message_text <- if (is.null(installation_error)) {
        "R did not report a specific error, but the dada2 package is still unavailable."
      } else {
        conditionMessage(installation_error)
      }
      showModal(modalDialog(
        title = "DADA2 installation failed",
        tags$p(message_text),
        tags$p("You can retry later by pressing Open FASTQ Pipeline."),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      pending_fastq_launch(FALSE)
      return()
    }

    status_text("DADA2 installed successfully.")
    showNotification("DADA2 was installed successfully.", type = "message", duration = 6)
    if (isTRUE(pending_fastq_launch())) {
      pending_fastq_launch(FALSE)
      launch_fastq()
    }
  })
}

shinyApp(ui, server)
