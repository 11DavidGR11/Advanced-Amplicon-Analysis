server <- function(input, output, session) {
  aaa_register_context_help(input, output, session, project_root, "Biological")

  function_group_input_ids <- stats::setNames(
    paste0("function_group_", gsub("[^A-Za-z0-9]+", "_", names(function_groups))),
    names(function_groups)
  )

  observe({
    selected <- unique(unlist(lapply(function_group_input_ids, function(id) input[[id]] %||% character()), use.names=FALSE))
    updateCheckboxGroupInput(session, "functional_functions", selected=selected)
  })

  update_function_groups <- function(selected) {
    selected <- intersect(unique(selected), names(biological_function_registry))
    for (category in names(function_groups)) {
      updateCheckboxGroupInput(session, function_group_input_ids[[category]],
        selected = intersect(selected, unname(function_groups[[category]])))
    }
  }

  observeEvent(input$functions_all, {
    update_function_groups(names(biological_function_registry))
  })
  observeEvent(input$functions_none, {
    update_function_groups(character())
  })
  observeEvent(input$function_preset, {
    req(nzchar(input$function_preset))
    update_function_groups(function_presets[[input$function_preset]] %||% character())
  }, ignoreInit=TRUE)
  observeEvent(input$function_search, {
    req(nzchar(input$function_search))
    selected <- unique(c(input$functional_functions %||% character(), input$function_search))
    update_function_groups(selected)
    updateSelectizeInput(session, "function_search", selected="")
  }, ignoreInit=TRUE)

  output$function_selection_summary <- renderUI({
    ids <- input$functional_functions %||% character()
    if (!length(ids)) return(tags$div(class="alert alert-warning py-2", "No biological functions selected."))
    cats <- unique(function_catalogue$Category[match(ids, function_catalogue$ID)])
    tags$div(class="alert alert-light border py-2 small",
      tags$strong(length(ids), " function(s) selected"),
      tags$span(" — ", paste(stats::na.omit(cats), collapse=", ")))
  })

  # Whether the design is replicated is a property of the design itself, not of
  # the replicate keyword: that control is hidden in metadata mode and keeps
  # whatever value it last held, which would silently drop the supervised
  # analyses from a perfectly well replicated metadata design.
  design_is_replicated <- reactive({
    sizes <- design_group_sizes()
    if (!is.null(sizes) && length(sizes)) {
      return(length(sizes) >= 2L && all(sizes >= 2L))
    }
    !identical(input$replicates %||% "none", "none")
  })

  selected_analyses <- reactive({
    c(
      if (isTRUE(input$use_functional_potential)) "functional_potential",
      if (isTRUE(input$use_top_abundance)) "top_abundance",
      if (isTRUE(input$use_differential_abundance)) "differential_abundance",
      if (isTRUE(input$use_functional_abundance)) "functional_abundance",
      if (isTRUE(input$analysis_differential_functions)) "differential_functions",
      if (isTRUE(input$analysis_functional_enrichment)) "functional_enrichment",
      if (isTRUE(input$use_community_structure)) "community_structure",
      if (isTRUE(input$use_plsda) && design_is_replicated()) "plsda",
      if (isTRUE(input$use_splsda) && design_is_replicated()) "splsda",
      if (isTRUE(input$use_rda)) "rda",
      if (isTRUE(input$analysis_envfit)) "envfit",
      if (isTRUE(input$analysis_partial_rda)) "partial_rda",
      if (isTRUE(input$analysis_dbrda)) "dbrda",
      if (isTRUE(input$analysis_partial_dbrda)) "partial_dbrda",
      if (isTRUE(input$analysis_varpart)) "variance_partitioning",
      if (isTRUE(input$analysis_ancombc2)) "ancombc2",
      if (isTRUE(input$analysis_maaslin)) "maaslin"
    )
  })

  # Keep the hidden legacy input synchronized for compatibility with any
  # third-party code, but all UI and workflow logic uses the reactive
  # selection above as the single source of truth.
  observeEvent(selected_analyses(), {
    updateCheckboxGroupInput(session, "analyses", selected = selected_analyses())
  }, ignoreInit = FALSE)


  v3_community_selected <- reactive({
    any(vapply(c(
      "analysis_alpha", "analysis_beta", "analysis_pca", "analysis_pcoa",
      "analysis_nmds", "analysis_clustering", "analysis_permanova",
      "analysis_anosim", "analysis_pairwise_permanova", "analysis_beta_dispersion"
    ), function(id) isTRUE(input[[id]]), logical(1)))
  })

  observe({
    updateCheckboxInput(session, "use_community_structure", value = v3_community_selected())
  })


  output$taxon_association_availability_ui <- renderUI({
    if (is.null(state$data)) {
      return(tags$div(class="alert alert-warning py-2 small", "Import an abundance table before configuring taxon associations."))
    }
    if (is.null(state$environmental_data)) {
      return(tags$div(class="alert alert-warning py-2 small", "ANCOM-BC2 and MaAsLin require sample metadata defining a factor or environmental predictor."))
    }
    tags$div(class="alert alert-success py-2 small", "Taxon-association inputs are available. Package dependencies are checked again before execution.")
  })

  selected_outputs <- reactive({
    selected_modules <- selected_analyses()
    automatic <- output_catalogue$ID[
      output_catalogue$Module %in% selected_modules &
        output_catalogue$Type != "figure"
    ]

    selected_figures <- unlist(list(
      if ("functional_potential" %in% selected_modules) input$figures_functional_potential,
      if ("functional_abundance" %in% selected_modules) input$figures_functional_abundance,
      if ("top_abundance" %in% selected_modules) input$figures_top_abundance,
      if ("differential_abundance" %in% selected_modules) input$figures_differential_abundance,
      if ("community_structure" %in% selected_modules) input$figures_community_structure,
      if ("community_structure" %in% selected_modules) input$figures_community_comparison,
      if ("plsda" %in% selected_modules) input$figures_plsda,
      if ("splsda" %in% selected_modules) input$figures_splsda,
      if ("rda" %in% selected_modules) input$figures_rda
    ), use.names = FALSE)

    ids <- unique(c(automatic, selected_figures %||% character()))

    if (length(input$functional_functions %||% character()) < 2L) {
      ids <- ids[ids != "functional_contributors"]
    }

    ids
  })


  state <- reactiveValues(
    data = NULL,
    results = NULL,
    log = character(),
    tables = list(),
    output_dir = NULL,
    workflow_started = NULL,
    environmental_data = NULL,
    metadata_roles = NULL,
    background = NULL,
    running = FALSE,
    cancel_requested = FALSE,
    progress_signature = NULL,
    selected_result_path = NULL,
    selected_result_data = NULL
  )

  dependency_revision <- reactiveVal(0L)

  # Apply a real HTML lock to unavailable analysis panels. A disabled
  # fieldset also disables controls created later inside dynamic uiOutput()
  # elements, so the lock cannot be bypassed after a reactive re-render.
  set_analysis_controls_enabled <- function(container_id, enabled) {
    disabled_js <- if (isTRUE(enabled)) "false" else "true"
    shinyjs::runjs(sprintf(
      paste0(
        "(function(){",
        "var root=document.getElementById('%s');",
        "if(!root){return;}",
        "root.disabled=%s;",
        "root.querySelectorAll('input,select,textarea,button').forEach(function(el){el.disabled=%s;});",
        "root.classList.toggle('triple-a-analysis-disabled', %s);",
        "})();"
      ),
      container_id, disabled_js, disabled_js, disabled_js
    ))
  }

  # Lock a whole accordion panel and, within it, individual analyses that have
  # stricter requirements than the panel as a whole. Container and per-item
  # states are applied in one place so the panel-wide pass cannot re-enable an
  # item that its own gate has just disabled.
  set_group_availability <- function(container_id, container_available, items = list()) {
    container_available <- isTRUE(container_available)

    resolved <- vapply(
      names(items),
      function(input_id) container_available && isTRUE(items[[input_id]]),
      logical(1)
    )

    # An analysis that just became unavailable must also leave the run.
    for (input_id in names(resolved)) {
      if (!resolved[[input_id]] && isTRUE(input[[input_id]])) {
        updateCheckboxInput(session, input_id, value = FALSE)
      }
    }

    apply_state <- function() {
      set_analysis_controls_enabled(container_id, container_available)
      for (input_id in names(resolved)) {
        shinyjs::toggleState(input_id, condition = resolved[[input_id]])
      }
    }

    # Reapply after the current reactive flush because Shiny/selectize may
    # replace input nodes while processing the same invalidating event.
    apply_state()
    session$onFlushed(apply_state, once = TRUE)
  }

  # A blocked analysis must say what is wrong AND what to do about it, so the
  # notice renders the eligibility reason followed by its guidance steps.
  analysis_availability_notice <- function(status, available_message) {
    if (isTRUE(status$available)) {
      return(tags$div(class = "alert alert-success py-2 small", available_message))
    }
    guidance <- as.character(status$guidance %||% character())
    guidance <- guidance[nzchar(trimws(guidance))]
    tags$div(
      class = "alert alert-warning py-2 small",
      tags$div(tags$strong("Locked. "), status$reason),
      if (length(guidance)) {
        tags$div(
          class = "mt-2",
          tags$strong("How to unlock:"),
          tags$ul(
            class = "mb-0 ps-3",
            lapply(guidance, function(step) tags$li(step))
          )
        )
      }
    )
  }

  # ---------------------------------------------------------------------
  # Analysis gating. Every module has an eligibility reactive backed by an
  # engine-level check, so the interface refuses exactly what the engine would
  # stop on, and every locked module explains how to unlock itself.
  # ---------------------------------------------------------------------

  top_abundance_eligibility <- reactive({
    aaa_check_top_abundance_eligibility(sample_columns())
  })

  output$top_abundance_availability_ui <- renderUI({
    analysis_availability_notice(
      top_abundance_eligibility(),
      "Top abundance is available: the selected sample columns can be ranked and summarised."
    )
  })

  observe({
    set_group_availability(
      "top_abundance_controls",
      top_abundance_eligibility()$available,
      list(use_top_abundance = TRUE)
    )
  })

  community_structure_eligibility <- reactive({
    aaa_check_community_structure_eligibility(
      sample_columns = sample_columns(),
      treatments = treatment_names()
    )
  })

  # Group-based comparisons need groups; the ordinations in the same module do
  # not, so they are gated separately instead of blocking the whole module.
  community_comparison_eligibility <- reactive({
    structure_status <- community_structure_eligibility()
    if (!isTRUE(structure_status$available)) {
      return(structure_status)
    }
    groups <- unique(trimws(treatment_names()))
    groups <- groups[nzchar(groups)]
    if (length(groups) < 2L) {
      return(aaa_eligibility_result(
        FALSE,
        "Community comparison is disabled: at least two treatment groups are required.",
        guidance = c(
          "Declare at least two distinct treatment names in 'Data and metadata'.",
          "PERMANOVA, ANOSIM and beta dispersion all compare groups against each other.",
          "Pairwise PERMANOVA additionally needs three or more groups.",
          "The ordinations in 'Community structure' remain available with a single group."
        )
      ))
    }
    aaa_eligibility_result(TRUE, details = list(groups = length(groups)))
  })

  output$community_structure_availability_ui <- renderUI({
    analysis_availability_notice(
      community_structure_eligibility(),
      "Diversity and ordination are available: enough samples are selected to build a two-dimensional ordination."
    )
  })

  output$community_comparison_availability_ui <- renderUI({
    analysis_availability_notice(
      community_comparison_eligibility(),
      "Community comparison is available: the declared groups can be tested against each other."
    )
  })

  # The individual checkboxes have to be listed so a locked panel also clears
  # them: use_community_structure is derived from them, so leaving them ticked
  # would keep the module in the run even though its panel is disabled.
  observe({
    available <- community_structure_eligibility()$available
    set_group_availability("diversity_controls", available, list(
      analysis_alpha = TRUE, analysis_beta = TRUE
    ))
    set_group_availability("community_structure_controls", available, list(
      analysis_pca = TRUE, analysis_pcoa = TRUE,
      analysis_nmds = TRUE, analysis_clustering = TRUE
    ))
    set_group_availability(
      "community_comparison_controls",
      community_comparison_eligibility()$available,
      list(
        analysis_permanova = TRUE, analysis_anosim = TRUE,
        analysis_pairwise_permanova = TRUE, analysis_beta_dispersion = TRUE
      )
    )
  })

  plsda_eligibility <- reactive({
    design_status <- sample_design_eligibility()
    if (!isTRUE(design_status$available)) return(design_status)
    aaa_check_plsda_eligibility(
      replicates = input$replicates %||% "none",
      sample_columns = sample_columns(),
      treatments = treatment_names(),
      label = "PLS-DA",
      group_sizes = design_group_sizes()
    )
  })

  output$plsda_availability_ui <- renderUI({
    analysis_availability_notice(
      plsda_eligibility(),
      "PLS-DA and sPLS-DA are available: the declared design has replicated groups to classify."
    )
  })

  observe({
    available <- plsda_eligibility()$available
    set_group_availability("plsda_controls", available, list(use_plsda = TRUE))
    set_group_availability("splsda_controls", available, list(use_splsda = TRUE))
  })

  functional_potential_eligibility <- reactive({
    aaa_check_functional_potential_eligibility(
      sample_columns = sample_columns(),
      selected_functions = input$functional_functions %||% character()
    )
  })

  output$functional_potential_availability_ui <- renderUI({
    analysis_availability_notice(
      functional_potential_eligibility(),
      "Functional potential is available: reference-genome evidence will be resolved for the selected functions."
    )
  })

  observe({
    set_group_availability(
      "functional_potential_controls",
      functional_potential_eligibility()$available,
      list(use_functional_potential = TRUE)
    )
  })

  functional_abundance_eligibility <- reactive({
    selected <- input$functional_functions %||% character()
    if (!length(selected)) {
      return(aaa_eligibility_result(
        FALSE,
        "Pathway-abundance analysis is disabled: no biological function is selected.",
        guidance = c(
          "Pick at least one function in the 'Functional analysis' panel, or apply one of the presets.",
          "This analysis aggregates the abundance of the taxa assigned to each selected function, so it needs at least one."
        )
      ))
    }
    aaa_eligibility_result(
      TRUE,
      details = list(selected_functions = length(selected))
    )
  })

  output$functional_abundance_availability_ui <- renderUI({
    analysis_availability_notice(
      functional_abundance_eligibility(),
      "Pathway-abundance analysis is available because at least one biological function is selected."
    )
  })

  observe({
    set_group_availability(
      "functional_abundance_controls",
      functional_abundance_eligibility()$available,
      list(use_functional_abundance = TRUE)
    )
  })

  differential_abundance_eligibility <- reactive({
    design_status <- sample_design_eligibility()
    if (!isTRUE(design_status$available)) return(design_status)
    aaa_check_differential_abundance_eligibility(
      sample_columns = sample_columns(),
      treatments = treatment_names(),
      replicates = input$replicates %||% "none",
      minimum_replicates = 2L,
      group_sizes = design_group_sizes()
    )
  })

  taxon_association_eligibility <- reactive({
    aaa_check_taxon_association_eligibility(
      sample_columns = sample_columns(),
      metadata = state$environmental_data,
      sample_id_column = input$environmental_sample_id %||% "",
      variables = c(selected_environmental_variables(), selected_experimental_factors())
    )
  })

  output$taxon_association_availability_ui <- renderUI({
    analysis_availability_notice(
      taxon_association_eligibility(),
      "Taxon-association analyses are available: metadata predictors are matched to the selected samples. Package dependencies are checked again before execution."
    )
  })

  observe({
    set_group_availability(
      "differential_abundance_controls",
      taxon_association_eligibility()$available,
      list(
        analysis_ancombc2 = TRUE,
        analysis_maaslin = TRUE,
        use_differential_abundance = differential_abundance_eligibility()$available
      )
    )
  })

  rda_eligibility <- reactive({
    aaa_check_rda_eligibility(
      sample_columns = sample_columns(),
      metadata = state$environmental_data,
      sample_id_column = input$environmental_sample_id %||% "",
      environmental_variables = selected_environmental_variables()
    )
  })

  # envfit, partial RDA and dbRDA also accept categorical predictors, so the
  # panel as a whole is gated by the looser constrained-ordination check and
  # only the RDA checkbox itself carries the numeric-predictor requirement.
  constrained_ordination_eligibility <- reactive({
    aaa_check_constrained_ordination_eligibility(
      sample_columns = sample_columns(),
      metadata = state$environmental_data,
      sample_id_column = input$environmental_sample_id %||% "",
      variables = c(selected_environmental_variables(), selected_experimental_factors()),
      label = "Environmental analysis"
    )
  })

  variance_partitioning_eligibility <- reactive({
    aaa_check_variance_partitioning_eligibility(
      sample_columns = sample_columns(),
      metadata = state$environmental_data,
      sample_id_column = input$environmental_sample_id %||% "",
      environmental_variables = selected_environmental_variables(),
      experimental_factors = selected_experimental_factors()
    )
  })

  output$constrained_ordination_availability_ui <- renderUI({
    analysis_availability_notice(
      constrained_ordination_eligibility(),
      "Environmental analyses are available: metadata predictors are matched to the selected samples."
    )
  })

  # Only shown once the panel itself is usable, so a single missing metadata
  # file does not produce two stacked warnings saying the same thing.
  output$rda_availability_ui <- renderUI({
    if (!isTRUE(constrained_ordination_eligibility()$available)) {
      return(NULL)
    }
    status <- rda_eligibility()
    if (isTRUE(status$available)) {
      return(NULL)
    }
    analysis_availability_notice(status, "")
  })

  output$variance_partitioning_availability_ui <- renderUI({
    if (!isTRUE(constrained_ordination_eligibility()$available)) {
      return(NULL)
    }
    status <- variance_partitioning_eligibility()
    if (isTRUE(status$available)) {
      return(NULL)
    }
    analysis_availability_notice(status, "")
  })

  observe({
    set_group_availability(
      "rda_controls",
      constrained_ordination_eligibility()$available,
      list(
        analysis_envfit = TRUE,
        analysis_partial_rda = TRUE,
        analysis_dbrda = TRUE,
        analysis_partial_dbrda = TRUE,
        use_rda = rda_eligibility()$available,
        analysis_varpart = variance_partitioning_eligibility()$available
      )
    )
  })

  pairwise_choices <- reactive({
    groups <- unique(trimws(treatment_names()))
    groups <- groups[nzchar(groups)]
    if (length(groups) < 2L) return(character())
    pairs <- utils::combn(groups, 2, simplify=FALSE)
    stats::setNames(vapply(pairs, function(z) paste(z, collapse="||"), character(1)),
                    vapply(pairs, function(z) paste(z[2], "vs", z[1]), character(1)))
  })
  output$pairwise_comparisons_ui <- renderUI({
    choices <- pairwise_choices()
    if (!length(choices)) return(tags$div(class="small text-muted", "Define at least two treatments to select pairwise comparisons."))
    checkboxGroupInput("differential_comparisons", "Pairwise comparisons to run", choices=choices, selected=unname(choices))
  })
  cache_revision <- reactiveVal(0L)
  rda_was_selected <- reactiveVal(FALSE)

  # Build the result inventory only when the active run changes. Selecting a
  # file must not rescan the filesystem or rebuild the complete tree.
  result_inventory <- reactive({
    req(state$output_dir)

    run_root <- normalizePath(
      state$output_dir,
      winslash = "/",
      mustWork = TRUE
    )

    files <- list.files(
      run_root,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE
    )

    if (!length(files)) {
      return(data.frame(
        path = character(),
        relative = character(),
        stringsAsFactors = FALSE
      ))
    }

    file_information <- file.info(files)
    files <- files[!is.na(file_information$isdir) & !file_information$isdir]

    if (!length(files)) {
      return(data.frame(
        path = character(),
        relative = character(),
        stringsAsFactors = FALSE
      ))
    }

    normalized <- normalizePath(
      files,
      winslash = "/",
      mustWork = TRUE
    )

    data.frame(
      path = normalized,
      relative = substring(normalized, nchar(run_root) + 2L),
      stringsAsFactors = FALSE
    )
  })

  # Keep non-reactive references for session cleanup. Shiny session-end
  # callbacks do not run inside a reactive consumer and therefore must not
  # read from or write to reactiveValues directly.
  cleanup_state <- new.env(parent = emptyenv())
  cleanup_state$background <- NULL
  cleanup_state$resource_alias <- NULL

  session$onSessionEnded(function() {
    background_process <- cleanup_state$background
    resource_alias <- cleanup_state$resource_alias

    if (!is.null(background_process)) {
      try(
        aaa_stop_background_workflow(background_process),
        silent = TRUE
      )
      try(
        aaa_cleanup_background_workflow(background_process),
        silent = TRUE
      )
    }

    if (
      is.character(resource_alias) &&
      length(resource_alias) == 1L &&
      !is.na(resource_alias) &&
      nzchar(resource_alias)
    ) {
      try(
        shiny::removeResourcePath(resource_alias),
        silent = TRUE
      )
    }

    cleanup_state$background <- NULL
    cleanup_state$resource_alias <- NULL
  })

  session$onFlushed(
    function() {
      shinyjs::disable("stop_run")
    },
    once = TRUE
  )

  update_progress_dom <- function(
      stage = NULL,
      detail = NULL,
      percent = NULL,
      elapsed = NULL,
      log_text = NULL,
      genome_current = NULL,
      genome_total = NULL,
      gene_current = NULL,
      gene_total = NULL,
      current_genome = NULL,
      current_pathway = NULL) {

    updates <- list(
      stage = stage,
      detail = detail,
      percent = percent,
      elapsed = elapsed,
      log_text = log_text,
      genome_current = genome_current,
      genome_total = genome_total,
      gene_current = gene_current,
      gene_total = gene_total,
      current_genome = current_genome,
      current_pathway = current_pathway
    )

    json <- jsonlite::toJSON(
      updates,
      auto_unbox = TRUE,
      null = "null",
      na = "null"
    )

    shinyjs::runjs(
      paste0(
        "(function(x){",
        "if(x.stage!==null){",
        "document.getElementById('triple_a_progress_stage').textContent=x.stage;",
        "}",
        "if(x.detail!==null){",
        "document.getElementById('triple_a_progress_detail').textContent=x.detail;",
        "}",
        "if(x.percent!==null){",
        "document.getElementById('triple_a_progress_percent').textContent=x.percent;",
        "}",
        "if(x.elapsed!==null){",
        "document.getElementById('triple_a_progress_elapsed').textContent=x.elapsed;",
        "}",
        "if(x.log_text!==null){",
        "var el=document.getElementById('triple_a_progress_log');",
        "if(el){el.textContent=x.log_text;el.scrollTop=el.scrollHeight;}",
        "}",
        "if(x.genome_current!==null&&x.genome_total!==null){",
        "var gc=document.getElementById('triple_a_genome_count');",
        "var gb=document.getElementById('triple_a_genome_bar');",
        "var gp=x.genome_total>0?100*x.genome_current/x.genome_total:0;",
        "if(gc){gc.textContent=x.genome_current+' / '+x.genome_total;}",
        "if(gb){gb.style.width=Math.max(0,Math.min(100,gp))+'%';}",
        "}",
        "if(x.gene_current!==null&&x.gene_total!==null){",
        "var gec=document.getElementById('triple_a_gene_count');",
        "var geb=document.getElementById('triple_a_gene_bar');",
        "var gep=x.gene_total>0?100*x.gene_current/x.gene_total:0;",
        "if(gec){gec.textContent=x.gene_current+' / '+x.gene_total;}",
        "if(geb){geb.style.width=Math.max(0,Math.min(100,gep))+'%';}",
        "}",
        "if(x.current_genome!==null){",
        "var cg=document.getElementById('triple_a_current_genome');",
        "if(cg){cg.textContent=x.current_genome;}",
        "}",
        "if(x.current_pathway!==null){",
        "var cp=document.getElementById('triple_a_current_pathway');",
        "if(cp){cp.textContent=x.current_pathway;}",
        "}",
        "})(",
        json,
        ");"
      )
    )
  }

  format_elapsed <- function(start_time) {
    if (is.null(start_time)) {
      return("00:00:00")
    }

    seconds <- max(
      0,
      as.numeric(
        difftime(
          Sys.time(),
          start_time,
          units = "secs"
        )
      )
    )

    sprintf(
      "%02d:%02d:%02d",
      floor(seconds / 3600),
      floor((seconds %% 3600) / 60),
      floor(seconds %% 60)
    )
  }

  append_log <- function(text) {
    state$log <- c(
      state$log,
      paste0(
        format(Sys.time(), "%H:%M:%S"),
        " | ",
        text
      )
    )

    update_progress_dom(
      elapsed = format_elapsed(
        state$workflow_started
      ),
      log_text = paste(
        state$log,
        collapse = "\n"
      )
    )

    invisible(text)
  }

  observeEvent(input$input_file, {
    req(input$input_file)

    state$data <- tryCatch({
      path <- input$input_file$datapath
      extension <- tolower(tools::file_ext(input$input_file$name))

      imported <- aaa_import_table(
        path = path,
        original_name = input$input_file$name
      )

      imported <- as.data.frame(
        imported,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      if (ncol(imported) < 2L) {
        stop(
          paste0(
            "The abundance table was imported as a single column. ",
            "Check that the delimiter matches the file extension (. ",
            extension,
            ")."
          ),
          call. = FALSE
        )
      }

      imported
    }, error = function(e) {
      showNotification(
        conditionMessage(e),
        type = "error",
        duration = NULL
      )
      NULL
    })
  })


  output$environmental_tab_ui <- renderUI({
    tagList(
      layout_columns(
        class = "metadata-layout",
        card(
          card_header("Optional sample metadata"),
          tags$div(
            class = "p-3",
            tags$div(class = "alert alert-light border",
              tags$strong("Format: "),
              "one row per sample. The identifier column must match abundance-table sample-column names. Metadata can be added or replaced without changing the abundance table."
            ),
            fileInput("environmental_file", "Sample variables and experimental factors",
              accept=c(".csv", ".tsv", ".txt", ".xls", ".xlsx"), width="100%"),
            tags$p(class="small text-muted mb-0",
              "Triple_A does not restrict column names. Assign the analytical role of every column after import."
            )
          )
        ),
        card(
          card_header("Column roles and values"),
          tags$div(class="p-3", uiOutput("environmental_controls_ui"))
        ),
        col_widths=c(5,7)
      ),
      card(
        card_header("Environmental-model defaults"),
        tags$div(class="p-3",
          layout_columns(
            numericInput("rda_permutations", "Permutation count", value=999,min=99,step=100),
            numericInput("rda_alpha", "Significance threshold", value=.05,min=.001,max=.99,step=.01),
            col_widths=c(6,6)
          ),
          tags$p(class="small text-muted mb-0",
            "These defaults are shared by RDA, dbRDA, envfit and related permutation tests."
          )
        )
      )
    )
  })

  persist_metadata_role_inputs <- function() {
    if (is.null(state$metadata_roles) || nrow(state$metadata_roles) == 0L) return(invisible(NULL))

    mapping <- state$metadata_roles
    if (!"Selected_role" %in% names(mapping)) {
      mapping$Selected_role <- mapping$Suggested_role
    }

    for (i in seq_len(nrow(mapping))) {
      value <- input[[paste0("metadata_role_", i)]]
      if (!is.null(value) && length(value) == 1L && nzchar(value)) {
        mapping$Selected_role[i] <- value
      }
    }

    selected_id <- input$environmental_sample_id
    if (!is.null(selected_id) && length(selected_id) == 1L && nzchar(selected_id)) {
      state$environmental_sample_id_value <- selected_id
    }

    state$metadata_roles <- mapping
    invisible(NULL)
  }

  observeEvent(selected_analyses(), {
    # Preserve manually assigned metadata roles whenever analysis controls
    # invalidate. Metadata now live in the first workflow tab, so no automatic
    # tab navigation is required or permitted.
    persist_metadata_role_inputs()

    selected <- "rda" %in% selected_analyses()
    if (selected && !isTRUE(rda_was_selected())) {
      showNotification(
        "RDA selected. Verify the metadata roles in '1. Data and metadata'.",
        type = "message",
        duration = 6
      )
    }
    rda_was_selected(selected)
  }, ignoreInit = TRUE)

  infer_metadata_type <- function(values) {
    non_missing <- values[!is.na(values) & trimws(as.character(values)) != ""]
    if (length(non_missing) == 0L) return("empty")

    text <- as.character(non_missing)
    numeric_values <- suppressWarnings(as.numeric(text))
    if (all(is.finite(numeric_values))) return("numeric")

    logical_values <- tolower(text)
    if (all(logical_values %in% c("true", "false", "yes", "no", "0", "1"))) {
      return("logical")
    }

    is_unambiguous_date <- function(x) {
      date_patterns <- c(
        "^\\d{4}-\\d{2}-\\d{2}$",       # YYYY-MM-DD
        "^\\d{4}/\\d{2}/\\d{2}$",       # YYYY/MM/DD
        "^\\d{2}/\\d{2}/\\d{4}$",       # DD/MM/YYYY or MM/DD/YYYY
        "^\\d{2}-\\d{2}-\\d{4}$",       # DD-MM-YYYY or MM-DD-YYYY
        "^\\d{4}-\\d{2}-\\d{2}[ T]\\d{2}:\\d{2}(:\\d{2})?$"
      )
      all(vapply(x, function(value) {
        any(vapply(date_patterns, grepl, logical(1), x = value, perl = TRUE))
      }, logical(1)))
    }

    if (is_unambiguous_date(text)) return("date")

    unique_count <- length(unique(text))
    if (unique_count <= max(10L, ceiling(0.25 * length(text)))) return("factor")
    "text"
  }

  suggest_metadata_role <- function(column_name, detected_type, position = 1L) {
    key <- tolower(gsub("[^a-z0-9]", "", column_name))
    if (position == 1L || key %in% c("sample", "sampleid", "sampleidentifier", "id")) {
      return("identifier")
    }
    if (grepl("treatment|group|condition|block|batch|farm|site|season|diet|timepoint", key)) {
      return("experimental_factor")
    }
    if (identical(detected_type, "numeric")) return("environmental_variable")
    if (identical(detected_type, "factor") || identical(detected_type, "logical")) {
      return("experimental_factor")
    }
    "ignore"
  }

  observeEvent(input$environmental_file, {
    req(input$environmental_file)
    path <- input$environmental_file$datapath
    extension <- tolower(tools::file_ext(input$environmental_file$name))
    state$environmental_data <- tryCatch(
      aaa_import_table(
        path = path,
        original_name = input$environmental_file$name
      ),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
        NULL
      }
    )

    if (!is.null(state$environmental_data)) {
      columns <- names(state$environmental_data)
      types <- vapply(state$environmental_data, infer_metadata_type, character(1))
      roles <- vapply(seq_along(columns), function(i) {
        suggest_metadata_role(columns[i], types[i], i)
      }, character(1))
      state$metadata_roles <- data.frame(
        Column = columns,
        Detected_type = unname(types),
        Suggested_role = unname(roles),
        Selected_role = unname(roles),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      identifier_index <- which(roles == "identifier")
      if (length(identifier_index) == 0L) identifier_index <- 1L
      state$environmental_sample_id_value <- columns[identifier_index[1]]
    }
  })

  metadata_role_map <- reactive({
    if (is.null(state$metadata_roles)) return(data.frame())
    mapping <- state$metadata_roles
    if (!"Selected_role" %in% names(mapping)) {
      mapping$Selected_role <- mapping$Suggested_role
    }
    mapping$Role <- vapply(seq_len(nrow(mapping)), function(i) {
      input[[paste0("metadata_role_", i)]] %||% mapping$Selected_role[i]
    }, character(1))

    selected_id <- input$environmental_sample_id %||%
      state$environmental_sample_id_value %||% mapping$Column[1]
    mapping$Role[mapping$Role == "identifier" & mapping$Column != selected_id] <- "ignore"
    mapping$Role[mapping$Column == selected_id] <- "identifier"
    mapping
  })

  selected_environmental_variables <- reactive({
    mapping <- metadata_role_map()
    if (nrow(mapping) == 0L) return(character())
    mapping$Column[
      mapping$Role == "environmental_variable" &
        mapping$Detected_type == "numeric"
    ]
  })

  selected_experimental_factors <- reactive({
    mapping <- metadata_role_map()
    if (nrow(mapping) == 0L) return(character())
    mapping$Column[mapping$Role == "experimental_factor"]
  })

  output$environmental_controls_ui <- renderUI({
    if (is.null(state$environmental_data) || is.null(state$metadata_roles)) {
      return(tags$div(
        class = "alert alert-secondary mb-0",
        "Upload a metadata file to activate these controls."
      ))
    }

    mapping <- state$metadata_roles

    # "Sample identifier" used to be offered here as well as in the Sample-ID
    # selector below. It never did anything: metadata_role_map() forces the
    # identifier role onto whichever column that selector names and demotes any
    # other column claiming it to "ignore". Two controls for one concept, one of
    # them inert, is what made the identifier and the grouping column look like
    # the same setting. The identifier is now chosen in exactly one place.
    role_choices <- c(
      "Experimental factor" = "experimental_factor",
      "Environmental variable" = "environmental_variable",
      "Ignore" = "ignore"
    )

    identifier_index <- which(mapping$Suggested_role == "identifier")
    if (length(identifier_index) == 0L) identifier_index <- 1L
    current_id <- input$environmental_sample_id %||%
      state$environmental_sample_id_value %||% mapping$Column[identifier_index[1]]
    current_group <- input$design_group_column %||% ""

    badge <- function(text, colour) {
      tags$div(
        class = "small fw-semibold",
        style = paste0("color:", colour, ";"),
        text
      )
    }

    role_rows <- lapply(seq_len(nrow(mapping)), function(i) {
      column <- mapping$Column[i]
      missing_fraction <- mean(is.na(state$environmental_data[[column]]))
      is_identifier <- identical(column, current_id)
      is_group <- identical(column, current_group)

      control <- if (is_identifier) {
        # The identifier carries no analytical role, so there is nothing to
        # choose: state what it is doing instead of showing a dead dropdown.
        badge("Sample ID — matched to the abundance table", "#1F4E79")
      } else {
        tagList(
          if (is_group) {
            badge("Grouping column — defines the experimental design", "#2E6B4F")
          },
          selectInput(
            paste0("metadata_role_", i),
            label = NULL,
            choices = role_choices,
            selected = {
              stored <- if ("Selected_role" %in% names(mapping)) {
                mapping$Selected_role[i]
              } else {
                mapping$Suggested_role[i]
              }
              if (stored %in% role_choices) stored else "ignore"
            },
            width = "100%"
          )
        )
      }

      tags$div(
        class = "border rounded p-2 mb-2",
        layout_columns(
          tags$div(
            tags$strong(column),
            tags$div(
              class = "small text-muted",
              paste0(
                "Detected type: ", mapping$Detected_type[i],
                " · Missing: ", round(100 * missing_fraction, 1), "%"
              )
            )
          ),
          control,
          col_widths = c(7, 5)
        )
      )
    })

    # The two columns that are not analytical variables but structural
    # declarations are chosen together, and next to the per-column roles, which
    # is where users look for "what does each column mean".
    tagList(
      selectInput(
        "environmental_sample_id",
        "Sample ID column",
        choices = mapping$Column,
        selected = current_id
      ),
      tags$p(
        class = "small text-muted mt-n2",
        "Its values must match the sample column names of the abundance table. ",
        "Differences in case, spaces, '-', '.' and '_' are tolerated."
      ),
      uiOutput("design_group_column_ui"),
      tags$hr(class = "my-2"),
      tags$div(class = "small text-muted mb-2", "Analytical role of the remaining columns:"),
      tags$div(role_rows),
      tags$div(
        class = "alert alert-light border mt-2 mb-0",
        tags$strong(nrow(state$environmental_data)), " metadata rows loaded. ",
        tags$span(textOutput("metadata_role_summary", inline = TRUE)),
        tags$br(),
        tags$span(textOutput("environmental_sample_match_summary", inline = TRUE))
      )
    )
  })

  output$metadata_role_summary <- renderText({
    mapping <- metadata_role_map()
    if (nrow(mapping) == 0L) return("")
    paste0(
      sum(mapping$Role == "experimental_factor"), " experimental factor(s); ",
      sum(mapping$Role == "environmental_variable"), " environmental variable(s); ",
      sum(mapping$Role == "ignore"), " ignored column(s)."
    )
  })

  output$environmental_sample_match_summary <- renderText({
    if (is.null(state$environmental_data)) return("")
    id_column <- input$environmental_sample_id %||% ""
    if (!nzchar(id_column) || !id_column %in% names(state$environmental_data)) {
      return("Select a valid sample-ID column to check sample matching.")
    }

    matches <- tryCatch(
      aaa_match_sample_ids(
        sample_columns(),
        state$environmental_data[[id_column]]
      ),
      error = function(e) e
    )
    if (inherits(matches, "error")) return(conditionMessage(matches))

    unmatched <- matches$reference_id[!matches$matched]
    normalised <- sum(matches$matched & !matches$exact)
    if (length(unmatched)) {
      return(paste0(
        "Matched abundance samples: ", sum(matches$matched), "/", nrow(matches),
        ". Missing metadata: ", paste(unmatched, collapse = ", "), "."
      ))
    }
    paste0(
      "Matched abundance samples: ", nrow(matches), "/", nrow(matches),
      if (normalised > 0L) paste0(" (", normalised, " matched after identifier normalisation).") else "."
    )
  })


  history_data <- reactiveVal(data.frame())
  history_aliases <- reactiveVal(character())
  refresh_history_data <- function() {
    runs_dir <- file.path(project_root, "Results", "Runs")
    runs <- if (dir.exists(runs_dir)) list.dirs(runs_dir, recursive = FALSE, full.names = TRUE) else character()
    old_aliases <- shiny::isolate(history_aliases())
    if (length(old_aliases)) for (alias in old_aliases) try(shiny::removeResourcePath(alias), silent = TRUE)
    aliases <- character()
    history <- if (length(runs) == 0) data.frame() else do.call(rbind, lapply(seq_along(runs), function(i) {
      run <- runs[[i]]
      metadata_file <- file.path(run, "Run_metadata.json")
      metadata <- if (file.exists(metadata_file)) tryCatch(jsonlite::read_json(metadata_file, simplifyVector = TRUE), error = function(e) NULL) else NULL
      status_record <- aaa_read_run_status(run)
      status <- status_record$status %||% metadata$status %||% "unknown"
      alias <- paste0("triple_a_history_", i, "_", abs(sum(utf8ToInt(basename(run)))))
      try(shiny::addResourcePath(alias, normalizePath(run, winslash = "/", mustWork = TRUE)), silent = TRUE)
      aliases <<- c(aliases, alias)
      report <- file.path(run, "Triple_A_report.html")
      open_link <- if (file.exists(report)) sprintf("<a href='/%s/%s' target='_blank' rel='noopener'>Open run</a>", alias, basename(report)) else ""
      cache_hits <- metadata$cache_hits %||% status_record$cache_hits
      cache_misses <- metadata$cache_misses %||% status_record$cache_misses
      cache_total <- (cache_hits %||% 0) + (cache_misses %||% 0)
      cache_label <- if (is.numeric(cache_hits) && is.numeric(cache_misses) && cache_total > 0) {
        sprintf("%d/%d cached", cache_hits, cache_total)
      } else {
        metadata$cache_status %||% status_record$cache_status %||% "UNVERIFIED"
      }
      data.frame(
        Run = basename(run),
        Status = tools::toTitleCase(status),
        Cache = cache_label,
        Date = metadata$started %||% as.character(file.info(run)$mtime),
        Duration_seconds = metadata$elapsed_seconds %||% NA,
        Analyses = if (!is.null(metadata$analyses)) paste(metadata$analyses, collapse = ", ") else "",
        Open = open_link,
        Directory = normalizePath(run, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    history_aliases(aliases)
    history_data(history)
  }
  observeEvent(input$refresh_history, refresh_history_data())
  session$onFlushed(
    function() {
      refresh_history_data()
    },
    once = TRUE
  )
  output$history_table <- renderDT({
    DT::datatable(history_data(), escape = FALSE, selection = "single", options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
  })

  selected_history_run <- reactive({
    idx <- input$history_table_rows_selected
    data <- history_data()
    if (length(idx) != 1L || !nrow(data)) return(NULL)
    data[idx, , drop = FALSE]
  })

  observeEvent(input$open_history_run, {
    run <- selected_history_run()
    req(run)
    report <- file.path(run$Directory[[1]], "Triple_A_report.html")
    if (!file.exists(report)) {
      showNotification("The selected run has no saved HTML report.", type = "warning")
      return(invisible(NULL))
    }
    browseURL(normalizePath(report, winslash = "/", mustWork = TRUE))
  })

  observeEvent(input$reload_history_run, {
    run <- selected_history_run()
    req(run)
    snapshot <- file.path(run$Directory[[1]], "Run_snapshot.rds")
    if (!file.exists(snapshot)) {
      showNotification(
        "This older run has no reloadable snapshot. New runs save one automatically.",
        type = "warning", duration = 9
      )
      return(invisible(NULL))
    }
    saved <- tryCatch(readRDS(snapshot), error = function(e) e)
    if (inherits(saved, "error")) {
      showNotification(paste("Unable to reload run:", conditionMessage(saved)), type = "error")
      return(invisible(NULL))
    }
    state$loaded_project_snapshot <- saved
    state$results <- saved$results %||% NULL
    if (!is.null(saved$dataset$abundance)) {
      state$data <- saved$dataset$abundance
      state$loaded_dataset <- saved$dataset
    }
    if (!is.null(saved$config$analyses)) {
      updateCheckboxInput(session, "use_community_structure", value = "community_structure" %in% saved$config$analyses)
      updateCheckboxInput(session, "use_rda", value = "rda" %in% saved$config$analyses)
      updateCheckboxInput(session, "use_plsda", value = isTRUE(saved$has_biological_replicates) && "plsda" %in% saved$config$analyses)
    }
    showNotification(
      "The saved dataset and analysis configuration were restored. Review sample selection before rerunning.",
      type = "message", duration = 9
    )
  })

  taxonomy_rank_labels <- c(
    Domain = "domain", Kingdom = "kingdom", Phylum = "phylum",
    Class = "class", Order = "order", Family = "family",
    Genus = "genus", Species = "species"
  )

  normalized_column_name <- function(x) {
    tolower(gsub("[^a-z0-9]", "", x))
  }

  guess_taxonomy_columns <- function(data) {
    columns <- names(data)
    keys <- normalized_column_name(columns)
    rank_keys <- c(
      "taxonomy", "taxon", "lineage", "classification",
      "domain", "superkingdom", "kingdom", "phylum", "class",
      "order", "family", "genus", "species"
    )
    columns[vapply(keys, function(key) {
      any(vapply(rank_keys, function(rank) grepl(rank, key, fixed = TRUE), logical(1)))
    }, logical(1))]
  }

  guess_rank_column <- function(data, rank) {
    columns <- names(data)
    keys <- normalized_column_name(columns)
    aliases <- switch(
      rank,
      domain = c("domain", "superkingdom", "d"),
      kingdom = c("kingdom", "k"),
      phylum = c("phylum", "p"),
      class = c("class", "taxclass", "c"),
      order = c("order", "taxorder", "o"),
      family = c("family", "f"),
      genus = c("genus", "g"),
      species = c("species", "s"),
      rank
    )
    exact <- which(keys %in% aliases)
    if (length(exact)) return(columns[exact[1]])
    contains <- which(vapply(keys, function(key) any(vapply(aliases[nchar(aliases) > 1], function(a) grepl(a, key, fixed = TRUE), logical(1))), logical(1)))
    if (length(contains)) columns[contains[1]] else ""
  }

  selected_rank_mapping <- reactive({
    ranks <- unname(taxonomy_rank_labels)
    ids <- paste0("taxonomy_rank_", ranks)
    values <- vapply(ids, function(id) input[[id]] %||% "", character(1))
    stats::setNames(values[nzchar(values)], ranks[nzchar(values)])
  })

  selected_rank_columns <- reactive({
    unname(selected_rank_mapping())
  })

  output$abundance_import_controls_ui <- renderUI({
    if (is.null(state$data)) {
      return(tags$div(
        class = "alert alert-light border small",
        tags$strong("Import workflow"),
        tags$br(),
        "Load an abundance table. Triple_A will then show a preview and ask how its columns should be interpreted."
      ))
    }

    columns <- names(state$data)
    numeric_columns <- columns[vapply(state$data, function(x) {
      values <- suppressWarnings(as.numeric(as.character(x)))
      length(values) > 0 && mean(!is.na(values)) >= 0.8
    }, logical(1))]
    tax_candidates <- guess_taxonomy_columns(state$data)
    single_candidates <- tax_candidates[
      grepl("tax|lineage|classification", tax_candidates, ignore.case = TRUE)
    ]
    default_tax <- if (length(single_candidates)) single_candidates[1] else if (length(tax_candidates)) tax_candidates[1] else columns[1]
    detected_rank_count <- sum(vapply(unname(taxonomy_rank_labels), function(rank) nzchar(guess_rank_column(state$data, rank)), logical(1)))
    default_layout <- if (detected_rank_count >= 2L) "separate" else "single"
    none_choice <- stats::setNames("", "Not present")
    column_choices_optional <- c(none_choice, stats::setNames(columns, columns))

    tagList(
      tags$div(
        class = "alert alert-info py-2 small",
        icon("wand-magic-sparkles"),
        paste0(
          "Detected ", ncol(state$data), " columns, including ",
          length(numeric_columns), " numeric candidate sample columns. Confirm the mapping below."
        )
      ),
      tags$div(
        class = "border rounded p-3 mb-3 bg-light",
        tags$div(class = "fw-semibold mb-2", "1. Sample columns"),
        radioButtons(
          "sample_selection_mode", NULL,
          choices = c(
            "Analyze all numeric columns not assigned to taxonomy" = "all_numeric",
            "Select sample columns manually" = "manual",
            "Select columns containing a shared text pattern" = "pattern"
          ),
          selected = "all_numeric"
        ),
        conditionalPanel(
          condition = "input.sample_selection_mode == 'manual'",
          selectizeInput(
            "sample_columns_manual", "Columns containing abundance values",
            choices = columns,
            selected = setdiff(numeric_columns, tax_candidates),
            multiple = TRUE,
            options = list(plugins = list("remove_button"))
          )
        ),
        conditionalPanel(
          condition = "input.sample_selection_mode == 'pattern'",
          textInput(
            "sample_identifier", "Text shared by sample-column names",
            value = "Sample", placeholder = "e.g. Sample, Treatment or _R"
          ),
          tags$p(
            class = "small text-muted",
            "Matching is case-insensitive. Several alternatives may be separated with commas or semicolons."
          )
        ),
        uiOutput("preview_status"),
        uiOutput("selected_samples_ui"),
        radioButtons(
          "identifier_column_mode", "Additional identifier column",
          choices = c(
            "No additional identifier column" = "none",
            "Exclude an identifier column from abundance values" = "column"
          ),
          selected = "none"
        ),
        conditionalPanel(
          condition = "input.identifier_column_mode == 'column'",
          selectInput(
            "sample_id_column", "Identifier column",
            choices = columns,
            selected = columns[1]
          )
        ),
        tags$p(
          class = "small text-muted mb-0",
          "The identifier is optional. When omitted, Triple_A analyzes all selected sample columns."
        )
      ),
      tags$div(
        class = "border rounded p-3 mb-3 bg-light",
        tags$div(class = "fw-semibold mb-2", "2. Taxonomic information"),
        radioButtons(
          "taxonomy_layout", NULL,
          choices = c(
            "The complete lineage is stored in one column" = "single",
            "Each taxonomic rank is stored in a separate column" = "separate"
          ),
          selected = default_layout
        ),
        conditionalPanel(
          condition = "input.taxonomy_layout == 'single'",
          selectInput(
            "taxonomy_column", "Column containing the taxonomic lineage",
            choices = columns, selected = default_tax
          ),
          tags$p(
            class = "small text-muted",
            "Accepted examples include d__/p__/g__ lineages, rank labels such as genus=..., and unprefixed semicolon- or pipe-separated lineages."
          )
        ),
        conditionalPanel(
          condition = "input.taxonomy_layout == 'separate'",
          tags$p(
            class = "small text-muted",
            "Assign only the ranks present in the file. Missing ranks may remain as ‘Not present’."
          ),
          lapply(names(taxonomy_rank_labels), function(label) {
            rank <- taxonomy_rank_labels[[label]]
            selectInput(
              paste0("taxonomy_rank_", rank), label,
              choices = column_choices_optional,
              selected = guess_rank_column(state$data, rank)
            )
          })
        )
      ),
      tags$div(
        class = "border rounded p-3 mb-3 bg-light",
        tags$div(class = "fw-semibold mb-2", "3. Experimental design"),
        radioButtons(
          "design_source", "How samples are assigned to groups",
          choices = c(
            "By column order, in consecutive blocks" = "blocks",
            "From a metadata column" = "metadata"
          ),
          selected = "blocks"
        ),
        conditionalPanel(
          condition = "input.design_source == 'blocks'",
          selectInput(
            "replicates", "Replicates per treatment",
            choices = c(
              "One sample per treatment" = "none",
              "Two replicates" = "duplicate",
              "Three replicates" = "triplicate",
              "Four replicates" = "quadruplicate",
              "Five replicates" = "quintuplicate"
            ),
            selected = "triplicate"
          ),
          tags$p(
            class = "small text-muted",
            paste(
              "Samples are split into consecutive blocks of the selected size,",
              "in the order the sample columns appear. All groups must therefore",
              "have the same number of replicates, and the column order must",
              "already be grouped by treatment."
            )
          ),
          conditionalPanel(
            condition = "input.replicates == 'none'",
            checkboxInput(
              "keep_original_sample_names",
              "Keep original sample-column names in final results",
              value = TRUE
            )
          )
        ),
        conditionalPanel(
          condition = "input.design_source == 'metadata'",
          tags$p(
            class = "small text-muted",
            paste(
              "Each sample takes the group named in the metadata column selected",
              "under \"Grouping column\", in the Column roles and values panel below.",
              "Column order is irrelevant and groups may have different numbers",
              "of replicates."
            )
          )
        ),
        uiOutput("sample_design_summary_ui"),
        selectInput(
          "abundance_type", "Abundance values",
          choices = c(
            "Proportions (0–1)" = "proportion",
            "Percentages (0–100)" = "percentage",
            "Raw read counts" = "counts"
          ),
          selected = "proportion"
        )
      )
    )
  })

  canonical_abundance_data <- reactive({
    req(state$data)
    data <- as.data.frame(state$data, check.names = FALSE, stringsAsFactors = FALSE)

    if (identical(input$taxonomy_layout %||% "single", "separate")) {
      rank_mapping <- selected_rank_mapping()
      rank_mapping <- rank_mapping[rank_mapping %in% names(data)]
      rank_columns <- unname(rank_mapping)
      rank_names <- names(rank_mapping)
      validate(need(length(rank_columns) > 0, "Select at least one taxonomic-rank column."))

      prefix_map <- c(
        domain = "d__", kingdom = "k__", phylum = "p__", class = "c__",
        order = "o__", family = "f__", genus = "g__", species = "s__"
      )
      taxonomy <- vapply(seq_len(nrow(data)), function(i) {
        pieces <- vapply(seq_along(rank_columns), function(j) {
          value <- trimws(as.character(data[[rank_columns[j]]][i]))
          if (is.na(value) || !nzchar(value) || tolower(value) %in% c("na", "none", "unclassified", "unknown")) return(NA_character_)
          if (grepl("^[a-zA-Z]__", value)) value else paste0(prefix_map[[rank_names[j]]], value)
        }, character(1))
        paste(pieces[!is.na(pieces)], collapse = ";")
      }, character(1))
      data$Taxonomy <- taxonomy
    } else {
      taxonomy_column <- input$taxonomy_column %||% ""
      validate(need(taxonomy_column %in% names(data), "Select the taxonomy column."))
      data$Taxonomy <- as.character(data[[taxonomy_column]])
    }

    data
  })

  sample_columns <- reactive({
    req(state$data)
    data <- canonical_abundance_data()
    columns <- names(data)

    excluded <- unique(c(
      "Taxonomy",
      input$taxonomy_column %||% character(),
      selected_rank_columns(),
      if (identical(input$identifier_column_mode %||% "none", "column")) input$sample_id_column %||% character()
    ))

    mode <- input$sample_selection_mode %||% "all_numeric"
    selected <- switch(
      mode,
      manual = input$sample_columns_manual %||% character(),
      pattern = {
        identifier <- trimws(input$sample_identifier %||% "")
        patterns <- trimws(unlist(strsplit(identifier, "[,;]", perl = TRUE)))
        patterns <- patterns[nzchar(patterns)]
        if (!length(patterns)) {
          character()
        } else {
          matched <- vapply(columns, function(column_name) {
            any(vapply(patterns, function(pattern) {
              grepl(pattern, column_name, ignore.case = TRUE, fixed = TRUE)
            }, logical(1)))
          }, logical(1))
          columns[matched]
        }
      },
      all_numeric = columns[vapply(data, function(x) {
        values <- suppressWarnings(as.numeric(as.character(x)))
        mean(!is.na(values)) >= 0.8
      }, logical(1))],
      character()
    )

    selected <- setdiff(selected, excluded)
    selected[selected %in% columns]
  })

  preview_data <- reactive({
    req(state$data)
    data <- canonical_abundance_data()
    visible_columns <- unique(c("Taxonomy", sample_columns()))
    data[, visible_columns[visible_columns %in% names(data)], drop = FALSE]
  })

  output$preview_status <- renderUI({
    req(state$data)

    matched <- sample_columns()
    hidden <- setdiff(names(state$data), names(preview_data()))

    tags$div(
      class = if (length(matched) > 0) "alert alert-light border py-2" else "alert alert-warning py-2",
      tags$strong(length(matched)),
      " sample column(s) selected. ",
      length(hidden),
      " non-analysis column(s) are hidden from the preview."
    )
  })

  output$selected_samples_ui <- renderUI({
    req(state$data)
    selected <- sample_columns()
    if (!length(selected)) {
      return(tags$div(
        class = "alert alert-warning py-2 mb-0",
        icon("triangle-exclamation"),
        " No sample columns match the current selection."
      ))
    }

    tags$div(
      class = "border rounded p-2 bg-white",
      tags$div(class = "fw-semibold mb-2", paste0("Selected samples (", length(selected), ")")),
      tags$div(
        style = "max-height: 150px; overflow-y: auto;",
        lapply(selected, function(sample_name) {
          tags$span(class = "badge text-bg-light border me-1 mb-1", sample_name)
        })
      )
    )
  })


  treatment_count <- reactive({
    columns <- sample_columns()
    n_reps <- aaa_n_replicates(
      input$replicates %||% "triplicate"
    )

    if (length(columns) == 0 ||
        length(columns) %% n_reps != 0) {
      return(NA_integer_)
    }

    length(columns) / n_reps
  })

  output$treatments_ui <- renderUI({
    # In metadata mode the group names come from the metadata column, so there
    # is nothing to type here.
    if (identical(input$design_source %||% "blocks", "metadata")) {
      return(NULL)
    }

    count <- treatment_count()

    if (is.na(count)) {
      return(tags$div(
        class = "alert alert-warning",
        "Sample columns and replicate design are not compatible."
      ))
    }

    no_replicates <- identical(input$replicates %||% "triplicate", "none")
    keep_original <- no_replicates && isTRUE(input$keep_original_sample_names %||% TRUE)

    if (keep_original) {
      return(tags$div(
        class = "alert alert-success py-2",
        icon("circle-check"),
        " Final sample names will be preserved: ",
        paste(sample_columns(), collapse = ", ")
      ))
    }

    textInput(
      "treatments",
      paste0("Treatment names — ", count, " required"),
      value = paste0("Treatment_", seq_len(count), collapse = ", ")
    )
  })

  # input$treatments is a free-text field: without debouncing, every
  # keystroke recomputes treatment_names() and its entire downstream chain
  # (pairwise_choices() -> pairwise UI, validation_report() -> validation
  # UI, selection_summary). Debouncing collapses rapid keystrokes into a
  # single recomputation once typing pauses, without changing behaviour.
  treatments_input_debounced <- shiny::debounce(
    reactive(input$treatments),
    millis = 500
  )

  # Group names for consecutive-block mode. Kept separate from treatment_names()
  # so that sample_design_table() never depends on a reactive that can itself
  # read the design back: treatment_names() -> sample_design_table() ->
  # treatment_names() is a cycle the moment the two branch differently.
  block_treatment_names <- reactive({
    no_replicates <- identical(input$replicates %||% "triplicate", "none")
    keep_original <- no_replicates && isTRUE(input$keep_original_sample_names %||% TRUE)
    if (keep_original) return(sample_columns())

    req(treatments_input_debounced())
    names <- trimws(strsplit(treatments_input_debounced(), ",", fixed = TRUE)[[1]])
    names[nzchar(names)]
  })

  treatment_names <- reactive({
    if (identical(input$design_source %||% "blocks", "metadata")) {
      design <- sample_design_table()
      if (is.null(design)) return(character())
      # Declared order, matching what the engine derives for samples_name, so
      # the interface and the results agree on which group comes first.
      return(if (is.factor(design$Treatment)) {
        levels(droplevels(design$Treatment))
      } else {
        unique(as.character(design$Treatment))
      })
    }
    block_treatment_names()
  })

  # ---------------------------------------------------------------------
  # Experimental design. This reactive is the single source of truth for
  # which sample belongs to which group; everything downstream (eligibility,
  # validation, the canonical dataset) reads it instead of re-deriving the
  # design from the column order.
  # ---------------------------------------------------------------------

  metadata_group_column <- reactive({
    if (!identical(input$design_source %||% "blocks", "metadata")) return("")
    selected <- input$design_group_column %||% ""
    if (nzchar(selected)) return(selected)
    # Fall back to the first column carrying the experimental-factor role, so
    # a typical metadata file works without extra configuration.
    factors <- selected_experimental_factors()
    if (length(factors)) factors[[1L]] else ""
  })

  # Rendered inside the metadata card, next to the sample-ID selector, because
  # both declare what a column IS rather than what it measures. Only relevant
  # when the design is read from the metadata.
  output$design_group_column_ui <- renderUI({
    if (!identical(input$design_source %||% "blocks", "metadata")) {
      return(tags$p(
        class = "small text-muted",
        tags$strong("Grouping column: "),
        "not used. The experimental design is currently taken from the column order; ",
        "switch it to \"From a metadata column\" under Experimental design to choose one here."
      ))
    }
    if (is.null(state$environmental_data)) {
      return(tags$div(
        class = "alert alert-warning py-2 small mb-2",
        "Upload a metadata file to use it as the design source."
      ))
    }
    columns <- setdiff(
      names(state$environmental_data),
      input$environmental_sample_id %||% ""
    )
    if (!length(columns)) {
      return(tags$div(
        class = "alert alert-warning py-2 small mb-2",
        "The metadata file has no column other than the sample identifier."
      ))
    }
    tagList(
      selectInput(
        "design_group_column", "Grouping column",
        choices = columns,
        selected = if ((input$design_group_column %||% "") %in% columns) {
          input$design_group_column
        } else {
          metadata_group_column()
        }
      ),
      tags$p(
        class = "small text-muted mt-n2",
        "It names the experimental group of each sample and defines the design. ",
        "Groups may be interleaved and may have different numbers of replicates."
      )
    )
  })

  sample_design_eligibility <- reactive({
    if (!identical(input$design_source %||% "blocks", "metadata")) {
      return(aaa_eligibility_result(TRUE))
    }
    aaa_check_sample_design_eligibility(
      sample_columns = sample_columns(),
      metadata = state$environmental_data,
      sample_id_column = input$environmental_sample_id %||% "",
      group_column = metadata_group_column()
    )
  })

  sample_design_table <- reactive({
    columns <- sample_columns()
    if (!length(columns)) return(NULL)

    if (identical(input$design_source %||% "blocks", "metadata")) {
      if (!isTRUE(sample_design_eligibility()$available)) return(NULL)
      return(aaa_sample_design_from_metadata(
        sample_columns = columns,
        metadata = state$environmental_data,
        sample_id_column = input$environmental_sample_id %||% "",
        group_column = metadata_group_column()
      ))
    }

    # Consecutive-block mode, kept for existing workflows.
    replicate_n <- aaa_n_replicates(input$replicates %||% "triplicate")
    names <- block_treatment_names()
    if (!length(names) || length(columns) != length(names) * replicate_n) return(NULL)
    data.frame(
      Sample_column = columns,
      Treatment = rep(names, each = replicate_n),
      Replicate = rep(seq_len(replicate_n), times = length(names)),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  })

  design_group_sizes <- reactive({
    design <- sample_design_table()
    if (is.null(design) || !nrow(design)) return(NULL)
    counts <- table(factor(design$Treatment, levels = unique(design$Treatment)))
    stats::setNames(as.integer(counts), names(counts))
  })

  output$sample_design_summary_ui <- renderUI({
    status <- sample_design_eligibility()
    if (!isTRUE(status$available)) {
      return(analysis_availability_notice(status, ""))
    }
    sizes <- design_group_sizes()
    if (is.null(sizes)) {
      if (identical(input$design_source %||% "blocks", "metadata")) return(NULL)
      return(tags$div(
        class = "alert alert-warning py-2 small mb-0",
        "Sample columns and the declared replicate design are not compatible yet."
      ))
    }
    balanced <- length(unique(sizes)) == 1L
    tags$div(
      class = "alert alert-success py-2 small mb-0",
      tags$strong(length(sizes), " group(s): "),
      paste(sprintf("%s (n=%d)", names(sizes), sizes), collapse = ", "),
      if (!balanced) {
        tags$div(
          class = "mt-1",
          "Unbalanced design. Group means and dispersion are computed per group, so this is supported; ",
          "note that tests with very different group sizes have less power for the smaller groups."
        )
      }
    )
  })


  methodology_parameters <- reactive({
    list(
      differential_abundance = list(
        method = input$method %||% "wilcox",
        paired = FALSE,
        min_prevalence = input$min_prevalence %||% 0.20,
        min_mean_abundance = input$min_mean_abundance %||% 0.01,
        alpha = input$alpha %||% 0.05,
        log2fc_threshold = input$log2fc %||% 1
      ),
      community_structure = list(
        transformation = input$community_transformation %||% "hellinger",
        distance_method = input$community_distance %||% "bray",
        permutations = input$permutations %||% 999,
        significance_alpha = input$community_alpha %||% 0.05,
        nmds_trymax = 50
      ),
      supervised_multivariate = list(
        transformation = input$community_transformation %||% "hellinger",
        plsda_components = input$plsda_components %||% 2,
        plsda_cv_folds = input$plsda_folds %||% 5,
        splsda_components = input$splsda_components %||% 2,
        splsda_cv_folds = input$splsda_folds %||% 5,
        splsda_repeats = input$splsda_repeats %||% 10,
        splsda_keepx = suppressWarnings(as.integer(trimws(strsplit(input$splsda_keepx %||% "5,10,20", ",")[[1]]))),
        splsda_tune = isTRUE(input$splsda_tune),
        rda_permutations = input$rda_permutations %||% 999,
        rda_alpha = input$rda_alpha %||% 0.05
      )
    )
  })

  dependency_input_files <- reactive({
    files <- character()

    if (!is.null(input$input_file)) {
      files <- c(
        files,
        input$input_file$name %||%
          input$input_file$datapath
      )
    }

    if (!is.null(input$environmental_file)) {
      files <- c(
        files,
        input$environmental_file$name %||%
          input$environmental_file$datapath
      )
    }

    files
  })

  selected_dependency_status <- reactive({
    dependency_revision()

    analyses <- selected_analyses() %||%
      character()

    aaa_dependency_status(
      analyses = analyses,
      input_files = dependency_input_files()
    )
  })

  missing_selected_packages <- reactive({
    status <- selected_dependency_status()

    status$Package[
      !status$Installed
    ]
  })

  output$dependency_status_ui <- renderUI({
    status <- selected_dependency_status()

    if (nrow(status) == 0) {
      return(NULL)
    }

    missing <- status[
      !status$Installed,
      ,
      drop = FALSE
    ]

    if (nrow(missing) == 0) {
      return(tags$div(
        class = paste(
          "alert alert-success",
          "small py-2"
        ),
        icon("circle-check"),
        " All dependencies required by the selected analyses are installed."
      ))
    }

    tags$div(
      class = paste(
        "alert alert-warning",
        "small py-2"
      ),
      icon("triangle-exclamation"),
      tags$strong(
        " Missing dependencies: "
      ),
      paste(
        missing$Package,
        collapse = ", "
      ),
      tags$br(),
      "Run will remain blocked until these packages are installed."
    )
  })

  show_dependency_modal <- function(
      missing_status) {

    package_items <- lapply(
      seq_len(nrow(missing_status)),
      function(i) {
        tags$li(
          tags$code(
            missing_status$Package[i]
          ),
          if (nzchar(
            missing_status$Required_by[i]
          )) {
            paste0(
              " — required by ",
              missing_status$Required_by[i]
            )
          }
        )
      }
    )

    showModal(
      modalDialog(
        title = "Required R packages are missing",
        tags$p(
          paste(
            "Triple_A cannot run the selected analyses",
            "until the following packages are installed:"
          )
        ),
        tags$ul(package_items),
        tags$p(
          class = "text-muted small",
          paste(
            "Installation uses your configured CRAN repository",
            "and may require write permission for the R library."
          )
        ),
        tags$div(
          id = "dependency_install_message"
        ),
        easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "install_selected_dependencies",
            "Install packages",
            icon = icon("download"),
            class = "btn-primary"
          )
        )
      )
    )
  }


  selected_methodology <- reactive({
    analyses <- selected_analyses()
    outputs <- selected_outputs()

    aaa_methodology_table(
      analyses = analyses,
      outputs = outputs,
      parameters = methodology_parameters()
    )
  })

  output$selected_methods_ui <- renderUI({
    methodology_cards(selected_methodology())
  })

  output$outputs_ui <- renderUI({
    # Retained only as a hidden compatibility bridge. Outputs are generated
    # automatically from the selected analyses and are not user-selectable.
    checkboxGroupInput(
      "outputs",
      NULL,
      choices = selected_outputs(),
      selected = selected_outputs()
    )
  })

  output$function_details_ui <- renderUI({
    ids <- input$functional_functions

    if (length(ids) == 0) return(NULL)

    catalogue <- triple_a_list_functions()
    selected <- catalogue[
      catalogue$ID %in% ids,
      ,
      drop = FALSE
    ]

    tags$div(
      lapply(
        seq_len(nrow(selected)),
        function(i) {
          tags$details(
            tags$summary(
              selected$Name[i]
            ),
            tags$p(
              selected$Description[i]
            ),
            tags$p(
              tags$strong("Diagnostic genes: "),
              selected$Diagnostic_genes[i]
            ),
            tags$p(
              tags$strong("Supporting genes: "),
              selected$Supporting_genes[i]
            )
          )
        }
      )
    )
  })

  abundance_nature_validation <- reactive({
    if (is.null(state$data)) {
      return(aaa_eligibility_result(
        FALSE,
        "Abundance format verification is waiting for an uploaded table."
      ))
    }
    aaa_validate_abundance_nature(
      data = canonical_abundance_data(),
      sample_columns = sample_columns(),
      abundance_type = input$abundance_type %||% "proportion"
    )
  })

  validation_report <- reactive({
    report <- aaa_validate_preflight(
      input_file = if (is.null(input$input_file)) NULL else input$input_file$datapath,
      data = if (is.null(state$data)) NULL else canonical_abundance_data(),
      sample_identifier = input$sample_identifier %||% "",
      sample_columns = sample_columns(),
      treatments = treatment_names(),
      replicates = input$replicates %||% "none",
      analyses = selected_analyses(),
      outputs = selected_outputs(),
      environmental_file = if (is.null(input$environmental_file)) NULL else input$environmental_file$datapath,
      group_sizes = design_group_sizes()
    )

    nature <- abundance_nature_validation()
    abundance_row <- aaa_validation_row(
      "Abundance format",
      if (isTRUE(nature$available)) "Valid" else if (is.null(state$data)) "Pending" else "Error",
      if (isTRUE(nature$available)) {
        paste0(
          "Values are compatible with the declared format: ",
          input$abundance_type %||% "proportion",
          "."
        )
      } else {
        nature$reason
      }
    )

    rbind(report, abundance_row)
  })

  output$validation_status_ui <- renderUI({
    report <- validation_report()
    tags$div(lapply(seq_len(nrow(report)), function(i) {
      status <- report$Status[i]
      cls <- switch(status, Valid="success", Warning="warning", Error="danger", Pending="secondary", "secondary")
      tags$div(class=paste("alert py-2 mb-2 alert", cls, sep="-"),
        tags$strong(paste0(report$Check[i], ": ")), report$Message[i])
    }))
  })

  output$preview <- renderDT({
    req(state$data)

    DT::datatable(
      prepare_dt_table(
        preview_data(),
        max_rows = 20
      ),
      options = list(
        scrollX = TRUE,
        pageLength = 10
      ),
      rownames = FALSE
    )
  }, server = FALSE)

  shiny::outputOptions(
    output,
    "preview",
    suspendWhenHidden = FALSE
  )

  output$selection_summary <- renderText({
    paste0(
      "Detected sample columns: ",
      length(sample_columns()),
      "\nTreatments: ",
      paste(treatment_names(), collapse = ", "),
      "\nProgress detail: ",
      input$progress_verbosity %||% "standard",
      "\nAnalyses: ",
      paste(selected_analyses(), collapse = ", "),
      "\nOutputs: ",
      paste(selected_outputs(), collapse = ", "),
      "\nFunctions: ",
      paste(
        input$functional_functions,
        collapse = ", "
      )
    )
  })

  observeEvent(
    input$install_selected_dependencies,
    {
      analyses <- selected_analyses() %||%
        character()

      missing_before <- missing_selected_packages()

      if (length(missing_before) == 0) {
        removeModal()
        return(invisible(NULL))
      }

      shinyjs::disable(
        "install_selected_dependencies"
      )

      shinyjs::html(
        id = "dependency_install_message",
        html = paste0(
          "<div class='alert alert-info'>",
          "Installing ",
          paste(
            missing_before,
            collapse = ", "
          ),
          "…</div>"
        )
      )

      installation_error <- NULL

      withProgress(
        message = "Installing R packages",
        value = 0,
        {
          for (i in seq_along(
            missing_before
          )) {
            package <- missing_before[i]

            setProgress(
              value = (i - 1) /
                length(missing_before),
              detail = paste(
                "Installing",
                package
              )
            )

            tryCatch(
              utils::install.packages(
                package
              ),
              error = function(e) {
                installation_error <<-
                  conditionMessage(e)
              }
            )
          }

          setProgress(
            value = 1,
            detail = "Checking installation"
          )
        }
      )

      missing_after <- aaa_missing_packages(
        analyses = analyses,
        input_files =
          dependency_input_files(),
        include_core = TRUE
      )

      if (length(missing_after) > 0) {
        shinyjs::enable(
          "install_selected_dependencies"
        )

        shinyjs::html(
          id = "dependency_install_message",
          html = paste0(
            "<div class='alert alert-danger'>",
            "<strong>Installation was not completed.</strong><br>",
            if (!is.null(
              installation_error
            )) {
              paste0(
                installation_error,
                "<br>"
              )
            } else {
              ""
            },
            "Still missing: ",
            paste(
              missing_after,
              collapse = ", "
            ),
            "<br><br>Install manually with:<br>",
            "<code>install.packages(c(",
            paste(
              sprintf(
                "&quot;%s&quot;",
                missing_after
              ),
              collapse = ", "
            ),
            "))</code>",
            "</div>"
          )
        )

        return(invisible(NULL))
      }

      dependency_revision(
        dependency_revision() + 1L
      )

      removeModal()

      showNotification(
        paste(
          "Dependencies installed successfully.",
          "Triple_A will now start the analysis."
        ),
        type = "message",
        duration = 5
      )

      session$onFlushed(
        function() {
          shinyjs::click("run")
        },
        once = TRUE
      )
    }
  )

  select_result_file <- function(path, notify_errors = TRUE) {
    if (is.null(path) || !length(path) || !nzchar(path[[1L]]) ||
        !file.exists(path[[1L]])) {
      state$selected_result_path <- NULL
      state$selected_result_data <- NULL
      return(invisible(FALSE))
    }

    candidate <- normalizePath(path[[1L]], winslash = "/", mustWork = TRUE)
    root <- normalizePath(state$output_dir, winslash = "/", mustWork = TRUE)

    if (!startsWith(tolower(candidate), paste0(tolower(root), "/"))) {
      if (isTRUE(notify_errors)) {
        showNotification(
          "The selected result is outside the current run directory.",
          type = "error"
        )
      }
      return(invisible(FALSE))
    }

    state$selected_result_path <- candidate
    state$selected_result_data <- NULL

    extension <- tolower(tools::file_ext(candidate))
    if (extension %in% c("csv", "tsv", "txt", "xls", "xlsx")) {
      file_size <- file.info(candidate)$size %||% 0
      preview_limit <- if (extension %in% c("xls", "xlsx")) {
        25 * 1024^2
      } else {
        50 * 1024^2
      }

      if (is.finite(file_size) && file_size > preview_limit) {
        if (isTRUE(notify_errors)) {
          showNotification(
            paste(
              "The selected table is too large for a safe in-app preview.",
              "Use Open to inspect it without blocking the Shiny session."
            ),
            type = "message",
            duration = 8
          )
        }
        return(invisible(TRUE))
      }

      state$selected_result_data <- tryCatch(
        aaa_import_table(
          path = candidate,
          original_name = basename(candidate)
        ),
        error = function(e) {
          if (isTRUE(notify_errors)) {
            showNotification(
              paste(
                "Unable to preview the selected table:",
                conditionMessage(e)
              ),
              type = "warning",
              duration = 8
            )
          }
          NULL
        }
      )
    }

    invisible(TRUE)
  }

  finalize_workflow_result <- function(result) {
    state$results <- result

    if (!is.null(result)) {
      append_log("Analysis completed.")

      update_progress_dom(
        stage = "Completed",
        detail = paste0(
          "All selected analyses finished. Results are available in ",
          state$output_dir,
          "."
        ),
        percent = "100%",
        elapsed = format_elapsed(
          state$workflow_started
        )
      )

      shinyjs::runjs(
        paste0(
          "['triple_a_genome_bar','triple_a_gene_bar'].forEach(",
          "function(id){",
          "var el=document.getElementById(id);",
          "if(el){",
          "el.classList.remove('progress-bar-animated');",
          "}",
          "});"
        )
      )

      if (!is.null(result$metadata$run_directory)) state$output_dir <- result$metadata$run_directory
      refresh_history_data()
      cache_revision(
        cache_revision() + 1L
      )
      state$tables <- aaa_collect_tables(result)

      collected_files <- aaa_collect_files(result)

      all_files <- unlist(
        collected_files,
        recursive = TRUE,
        use.names = TRUE
      )

      if (length(all_files) == 0) {
        all_files <- character()
      } else {
        all_files <- as.character(all_files)

        valid_paths <- !is.na(all_files) &
          nzchar(all_files)

        all_files <- all_files[valid_paths]

        if (length(all_files) > 0) {
          existing_paths <- vapply(
            all_files,
            function(path) {
              is.character(path) &&
                length(path) == 1L &&
                !is.na(path) &&
                nzchar(path) &&
                file.exists(path)
            },
            logical(1)
          )

          all_files <- all_files[existing_paths]
        }
      }

      image_extensions <- c(
        "png", "jpg", "jpeg", "webp", "gif"
      )

      if (length(all_files) == 0) {
        state$image_files <- character()
      } else {
        image_paths <- tolower(
          tools::file_ext(all_files)
        ) %in% image_extensions

        image_files <- all_files[image_paths]

        if (length(image_files) > 0) {
          normalized_images <- normalizePath(
            image_files,
            winslash = "/",
            mustWork = TRUE
          )

          keep_unique <- !duplicated(
            tolower(normalized_images)
          )

          image_files <- image_files[keep_unique]
          names(image_files) <- make.unique(
            basename(image_files)
          )
        }

        state$image_files <- image_files
      }

      if (!is.null(state$resource_alias)) {
        try(
          shiny::removeResourcePath(state$resource_alias),
          silent = TRUE
        )
        cleanup_state$resource_alias <- NULL
      }

      state$resource_alias <- paste0(
        "triple_a_results_",
        as.integer(Sys.time())
      )
      cleanup_state$resource_alias <- state$resource_alias

      shiny::addResourcePath(
        state$resource_alias,
        normalizePath(
          state$output_dir,
          winslash = "/",
          mustWork = TRUE
        )
      )

      table_choices <- names(state$tables)
      summary_first <- grepl("(^| / )summary$", table_choices, ignore.case = TRUE)
      table_choices <- c(table_choices[summary_first], table_choices[!summary_first])
      updateSelectInput(
        session,
        "table_name",
        choices = table_choices,
        selected = table_choices[1] %||% character()
      )

      # Select and load the first Summary file when available. Only one
      # result is previewed at a time.
      summary_candidates <- all_files[
        grepl("summary", basename(all_files), ignore.case = TRUE)
      ]
      default_result <- if (length(summary_candidates)) {
        summary_candidates[[1L]]
      } else if (length(all_files)) {
        all_files[[1L]]
      } else {
        NULL
      }
      select_result_file(default_result, notify_errors = FALSE)

      updateTabsetPanel(
        session,
        "main_tabs",
        selected = "Results"
      )
    }
  }

  finish_workflow_controls <- function() {
    state$running <- FALSE
    shinyjs::enable("run")
    shinyjs::disable("stop_run")
  }

  apply_background_progress <- function(event) {
    if (is.null(event) ||
        is.null(event$stage)) {
      return(invisible(NULL))
    }

    stage <- event$stage
    detail <- event$detail %||% ""
    completed <- event$completed %||% 0
    total <- event$total %||% 1

    monitor_event <- identical(
      stage,
      "functional_monitor"
    ) ||
      endsWith(
        stage,
        ".functional_monitor"
      )

    if (monitor_event) {
      monitor_detail <- sub(
        "^[^:]+:\\s*",
        "",
        detail
      )

      monitor <- strsplit(
        monitor_detail,
        "|||",
        fixed = TRUE
      )[[1]]

      if (length(monitor) >= 7) {
        update_progress_dom(
          stage =
            "Potential metabolic pathways analysis",
          detail = if (identical(
            monitor[1],
            "reference"
          )) {
            paste0(
              "Resolving cached taxonomic reference ",
              monitor[2],
              " of ",
              monitor[3],
              "."
            )
          } else {
            paste0(
              "Evaluating ",
              monitor[5],
              " in reference genomes."
            )
          },
          elapsed = format_elapsed(
            state$workflow_started
          ),
          genome_current =
            as.integer(monitor[2]),
          genome_total =
            as.integer(monitor[3]),
          current_genome = monitor[4],
          current_pathway = monitor[5],
          gene_current =
            as.integer(monitor[6]),
          gene_total =
            as.integer(monitor[7])
        )
      }

      return(invisible(NULL))
    }

    value <- max(
      0,
      min(
        1,
        completed / max(1, total)
      )
    )

    stage_labels <- c(
      functional_potential =
        "Potential metabolic pathways analysis",
      functional_potential.figures =
        "Potential metabolic pathways: figures",
      functional_potential.completed =
        "Potential metabolic pathways: completed",
      functional_abundance =
        "Potential metabolomic pathways abundance",
      functional_abundance.definitions =
        "Potential metabolomic pathways: definitions",
      functional_abundance.summary =
        "Potential metabolomic pathways: abundance summary",
      functional_abundance.outputs =
        "Potential metabolomic pathways: outputs",
      initialization = "Preparing workflow",
      functional_potential.cache =
        "Loading shared cache",
      functional_potential.ncbi_reference =
        "Resolving reference genomes",
      functional_potential.gene_search =
        "Evaluating marker genes",
      functional_potential.gff =
        "Reading cached GFF annotations",
      completed = "Completed"
    )

    stage_label <- if (
      stage %in% names(stage_labels)
    ) {
      unname(stage_labels[[stage]])
    } else {
      tools::toTitleCase(
        gsub("[._]", " ", stage)
      )
    }

    percentage <- paste0(
      round(100 * value),
      "%"
    )

    update_progress_dom(
      stage = stage_label,
      detail = detail,
      percent = percentage,
      elapsed = format_elapsed(
        state$workflow_started
      )
    )

    if (!identical(
      input$progress_verbosity,
      "standard"
    )) {
      append_log(
        paste0(
          percentage,
          " | ",
          stage_label,
          " | ",
          detail
        )
      )
    }

    invisible(NULL)
  }

  read_background_rds <- function(path, attempts = 3L, delay = 0.03) {
    for (attempt in seq_len(attempts)) {
      if (!file.exists(path)) {
        if (attempt < attempts) Sys.sleep(delay)
        next
      }
      value <- tryCatch(
        suppressWarnings(readRDS(path)),
        error = function(e) NULL
      )
      if (!is.null(value)) return(value)
      if (attempt < attempts) Sys.sleep(delay)
    }
    NULL
  }

  background_poll <- reactiveTimer(
    350,
    session
  )

  observe({
    background_poll()

    if (!isTRUE(state$running) ||
        is.null(state$background)) {
      return(invisible(NULL))
    }

    background <- state$background

    if (file.exists(
      background$progress
    )) {
      info <- file.info(
        background$progress
      )

      signature <- paste(
        info$size,
        as.numeric(info$mtime),
        sep = "|"
      )

      if (!identical(
        signature,
        state$progress_signature
      )) {
        event <- read_background_rds(
          background$progress
        )

        if (!is.null(event)) {
          state$progress_signature <-
            signature

          apply_background_progress(
            event
          )
        }
      }
    }

    if (background$process$is_alive()) {
      return(invisible(NULL))
    }

    if (isTRUE(
      state$cancel_requested
    )) {
      return(invisible(NULL))
    }

    result <- if (file.exists(
      background$result
    )) {
      read_background_rds(
        background$result
      )
    } else {
      NULL
    }

    if (!is.null(result)) {
      finalize_workflow_result(
        result
      )

      finish_workflow_controls()

      aaa_cleanup_background_workflow(
        background
      )

      state$background <- NULL
      cleanup_state$background <- NULL
      return(invisible(NULL))
    }

    error_information <- if (
      file.exists(background$error)
    ) {
      read_background_rds(
        background$error
      )
    } else {
      NULL
    }

    error_detail <- if (
      is.null(error_information)
    ) {
      paste(
        "The background workflow stopped before",
        "returning a result."
      )
    } else if (
      is.null(error_information$call)
    ) {
      error_information$message
    } else {
      paste0(
        error_information$message,
        " | Call: ",
        error_information$call
      )
    }

    append_log(
      paste(
        "ERROR:",
        error_detail
      )
    )

    if (!is.null(state$output_dir) && dir.exists(state$output_dir)) {
      try(aaa_write_run_status(state$output_dir, "failed", error_detail), silent = TRUE)
    }

    update_progress_dom(
      stage = "Error",
      detail = error_detail,
      elapsed = format_elapsed(
        state$workflow_started
      )
    )

    showNotification(
      error_detail,
      type = "error",
      duration = NULL
    )

    finish_workflow_controls()

    aaa_cleanup_background_workflow(
      background
    )

    state$background <- NULL
    cleanup_state$background <- NULL

    invisible(NULL)
  })

  observeEvent(input$stop_run, {
    if (!isTRUE(state$running) ||
        is.null(state$background)) {
      return(invisible(NULL))
    }

    state$cancel_requested <- TRUE

    update_progress_dom(
      stage = "Stopping",
      detail = paste(
        "Stopping the active workflow.",
        "Files already written will be preserved."
      ),
      elapsed = format_elapsed(
        state$workflow_started
      )
    )

    stopped <- aaa_stop_background_workflow(
      state$background
    )

    if (!is.null(state$output_dir) && dir.exists(state$output_dir)) {
      try(aaa_write_run_status(state$output_dir, "cancelled",
                               "Execution cancelled by the user."), silent = TRUE)
    }

    append_log(
      if (isTRUE(stopped)) {
        "Execution cancelled by the user."
      } else {
        paste(
          "A stop request was sent, but the",
          "background process did not exit immediately."
        )
      }
    )

    update_progress_dom(
      stage = "Cancelled",
      detail = paste(
        "The analysis was stopped by the user.",
        "Partial files already written have been preserved."
      ),
      elapsed = format_elapsed(
        state$workflow_started
      )
    )

    showNotification(
      "Triple_A execution cancelled.",
      type = "warning",
      duration = 6
    )

    aaa_cleanup_background_workflow(
      state$background
    )

    state$background <- NULL
    cleanup_state$background <- NULL
    finish_workflow_controls()

    invisible(NULL)
  })

  observeEvent(input$run, {
    validation <- validation_report()
    blocking <- validation$Status %in% c("Error", "Pending")
    if (any(blocking)) {
      showModal(modalDialog(
        title = "Run blocked by validation",
        tags$p("Correct the following items before running:"),
        tags$ul(lapply(which(blocking), function(i) tags$li(paste(validation$Check[i], validation$Message[i], sep=": ")))),
        easyClose = TRUE, footer = modalButton("Close")
      ))
      return(invisible(NULL))
    }

    req(input$input_file)
    req(length(selected_analyses()) > 0)
    req(length(selected_outputs()) > 0)

    # Final gate. The panels are locked reactively, but a stale browser state or
    # a programmatic input update could still submit an ineligible module, so
    # every selected analysis is re-checked here against the same engine rules
    # and the user is told how to unlock whatever failed.
    design_status <- sample_design_eligibility()
    if (!isTRUE(design_status$available)) {
      guidance <- as.character(design_status$guidance %||% character())
      showModal(modalDialog(
        title = "Run blocked: the experimental design is not usable",
        tags$div(tags$strong(design_status$reason)),
        if (length(guidance)) {
          tags$ul(class = "mt-2 ps-3", lapply(guidance, function(step) tags$li(step)))
        },
        easyClose = TRUE, footer = modalButton("Close")
      ))
      return(invisible(NULL))
    }

    active_analyses <- selected_analyses()
    gates <- list(
      top_abundance = top_abundance_eligibility,
      community_structure = community_structure_eligibility,
      differential_abundance = differential_abundance_eligibility,
      plsda = plsda_eligibility,
      splsda = plsda_eligibility,
      rda = rda_eligibility,
      envfit = constrained_ordination_eligibility,
      partial_rda = constrained_ordination_eligibility,
      dbrda = constrained_ordination_eligibility,
      partial_dbrda = constrained_ordination_eligibility,
      variance_partitioning = variance_partitioning_eligibility,
      ancombc2 = taxon_association_eligibility,
      maaslin = taxon_association_eligibility,
      functional_potential = functional_potential_eligibility,
      functional_abundance = functional_abundance_eligibility
    )

    blocked <- Filter(Negate(is.null), lapply(
      intersect(active_analyses, names(gates)),
      function(id) {
        status <- gates[[id]]()
        if (isTRUE(status$available)) NULL else list(id = id, status = status)
      }
    ))

    if (length(blocked)) {
      showModal(modalDialog(
        title = "Run blocked: some selected analyses do not meet their requirements",
        tags$p(
          "The analyses below cannot run with the current data and design. ",
          "Follow the steps to unlock them, or deselect them and run the rest."
        ),
        lapply(blocked, function(entry) {
          guidance <- as.character(entry$status$guidance %||% character())
          guidance <- guidance[nzchar(trimws(guidance))]
          tags$div(
            class = "mb-3",
            tags$div(tags$strong(entry$status$reason)),
            if (length(guidance)) {
              tags$ul(class = "mb-0 ps-3", lapply(guidance, function(step) tags$li(step)))
            }
          )
        }),
        easyClose = TRUE, footer = modalButton("Close")
      ))
      return(invisible(NULL))
    }

    abundance_status <- abundance_nature_validation()
    validate(need(isTRUE(abundance_status$available), abundance_status$reason))

    missing_status <- selected_dependency_status()
    missing_status <- missing_status[
      !missing_status$Installed,
      ,
      drop = FALSE
    ]

    if (nrow(missing_status) > 0) {
      show_dependency_modal(
        missing_status
      )

      return(invisible(NULL))
    }

    # RDA metadata requirements are covered by the eligibility gate above, which
    # reports them in the modal together with the steps to fix them instead of
    # silently aborting the observer.

    # Block mode is the only one that constrains the number of treatment names,
    # because it derives the groups from the column order.
    if (!identical(input$design_source %||% "blocks", "metadata")) {
      count <- treatment_count()
      validate(
        need(!is.na(count), "Invalid sample/replicate design."),
        need(
          length(treatment_names()) == count,
          paste0("Exactly ", count, " treatment names are required.")
        )
      )
    }

    updateTabsetPanel(
      session,
      "main_tabs",
      selected = "Run"
    )

    state$log <- character()
    state$workflow_started <- Sys.time()

    update_progress_dom(
      stage = "Preparing analysis",
      detail = "Validating the uploaded table and selected parameters.",
      percent = "0%",
      elapsed = "00:00:00",
      log_text = "",
      genome_current = 0,
      genome_total = 0,
      gene_current = 0,
      gene_total = 0,
      current_genome = "—",
      current_pathway = "—"
    )

    append_log("Preparing Triple_A configuration.")

    output_dir <- file.path(
      project_root,
      "Results"
    )
    state$output_dir <- output_dir

    selected_samples <- sample_columns()
    # Built by sample_design_table(), which handles both the consecutive-block
    # mode and the metadata-column mode. Deriving it again here is what used to
    # tie the whole application to the order of the sample columns.
    sample_design <- sample_design_table()
    if (is.null(sample_design) || !nrow(sample_design)) {
      showModal(modalDialog(
        title = "Run blocked: the experimental design is incomplete",
        tags$p("Triple_A could not determine which sample belongs to which group."),
        tags$ul(
          tags$li("In block mode, the number of sample columns must equal treatments x replicates."),
          tags$li("In metadata mode, every selected sample needs a row with a non-empty group value.")
        ),
        easyClose = TRUE, footer = modalButton("Close")
      ))
      return(invisible(NULL))
    }
    append_log(paste0(
      "Experimental design (",
      if (identical(input$design_source %||% "blocks", "metadata")) {
        paste0("from metadata column '", metadata_group_column(), "'")
      } else {
        "consecutive blocks"
      },
      "): ",
      paste(sprintf("%s=%d", names(design_group_sizes()), design_group_sizes()), collapse = ", ")
    ))
    canonical_dataset <- aaa_new_dataset(
      abundance = canonical_abundance_data(),
      sample_design = sample_design,
      metadata = state$environmental_data,
      metadata_roles = metadata_role_map(),
      source = list(
        file = input$input_file$name,
        environmental_file = if (is.null(input$environmental_file)) NULL else input$environmental_file$name
      )
    )

    config <- triple_a_config(
      dataset = canonical_dataset,
      abundance_type =
        input$abundance_type,
      output_dir = output_dir,
      progress_verbosity =
        input$progress_verbosity %||% "standard",
      analyses = selected_analyses(),
      outputs = selected_outputs(),
      functional_functions =
        input$functional_functions,
      top_abundance = list(
        top_n = input$top_n %||% 20
      ),
      differential_abundance = list(
        method = input$method %||% "wilcox",
        min_prevalence =
          input$min_prevalence %||% 0.20,
        min_mean_abundance =
          input$min_mean_abundance %||% 0.01,
        alpha = input$alpha %||% 0.05,
        log2fc_threshold =
          input$log2fc %||% 1,
        max_labels =
          input$max_labels %||% 10,
        label_only_significant = FALSE,
        colour_by = "log2FC",
        point_size = "mean_abundance",
        x_limit = 6,
        comparisons = input$differential_comparisons %||% character()
      ),
      functional_abundance = list(
        top_taxa_per_pathway =
          input$top_contributors %||% 5
      ),
      community_structure = list(
        transformation = input$community_transformation %||% "hellinger",
        distance_method = input$community_distance %||% "bray",
        permutations = input$permutations %||% 999,
        significance_alpha = input$community_alpha %||% 0.05,
        nmds_trymax = 50,
        show_sample_labels = input$ordination_labels %||% FALSE
      ),
      supervised_multivariate = list(
        transformation = input$community_transformation %||% "hellinger",
        plsda_components = input$plsda_components %||% 2,
        plsda_cv_folds = input$plsda_folds %||% 5,
        plsda_seed = 123,
        splsda_components = input$splsda_components %||% 2,
        splsda_cv_folds = input$splsda_folds %||% 5,
        splsda_repeats = input$splsda_repeats %||% 10,
        splsda_keepx = suppressWarnings(as.integer(trimws(strsplit(input$splsda_keepx %||% "5,10,20", ",")[[1]]))),
        splsda_tune = isTRUE(input$splsda_tune),
        splsda_seed = 123,
        rda_permutations = input$rda_permutations %||% 999,
        rda_alpha = input$rda_alpha %||% 0.05,
        show_sample_labels = input$ordination_labels %||% FALSE
      ),
      environmental = list(
        distance = input$environmental_distance %||% "bray",
        sample_id_column = input$environmental_sample_id %||% NULL,
        variables = selected_environmental_variables(),
        experimental_factors = selected_experimental_factors(),
        column_roles = metadata_role_map()
      )
    )

    state$cancel_requested <- FALSE
    state$progress_signature <- NULL
    state$running <- TRUE

    shinyjs::disable("run")
    shinyjs::enable("stop_run")

    append_log("Starting Triple_A in a background R process.")

    background <- tryCatch(
      aaa_start_background_workflow(
        config = config,
        project_root = project_root,
        verbose = identical(
          input$progress_verbosity,
          "developer"
        )
      ),
      error = function(e) {
        error_detail <- conditionMessage(e)

        append_log(
          paste(
            "ERROR:",
            error_detail
          )
        )

        update_progress_dom(
          stage = "Error",
          detail = error_detail,
          elapsed = format_elapsed(
            state$workflow_started
          )
        )

        showNotification(
          error_detail,
          type = "error",
          duration = NULL
        )

        finish_workflow_controls()

        NULL
      }
    )

    if (is.null(background)) {
      return(invisible(NULL))
    }

    state$background <- background
    cleanup_state$background <- background

    update_progress_dom(
      stage = "Starting background workflow",
      detail = paste(
        "The analysis is running independently.",
        "You may stop it without closing Shiny."
      ),
      percent = "0%",
      elapsed = format_elapsed(
        state$workflow_started
      )
    )

    invisible(NULL)

  })

  output$result_tree_ui <- renderUI({
    req(state$results, state$output_dir)

    inventory <- result_inventory()
    files <- inventory$path
    rel <- inventory$relative

    if (!length(files)) {
      return(tags$div(
        class = "alert alert-info",
        "No results are available."
      ))
    }

    tree_file_link <- function(path, relative_path) {
      js_value <- jsonlite::toJSON(
        relative_path,
        auto_unbox = TRUE
      )

      tags$a(
        href = "#",
        class = "result-tree-link d-block py-1 px-2 rounded text-decoration-none",
        onclick = paste0(
          "document.querySelectorAll('.result-tree-link.active').forEach(",
          "function(el){el.classList.remove('active','bg-primary','text-white');});",
          "this.classList.add('active','bg-primary','text-white');",
          "Shiny.setInputValue('result_tree_selection', ",
          js_value,
          ", {priority: 'event'}); return false;"
        ),
        icon("file"),
        tags$span(class = "ms-1 text-break", basename(path))
      )
    }

    classify_result <- function(path) {
      name <- tolower(basename(path))
      ext <- tolower(tools::file_ext(path))

      if (grepl("summary|overview", name)) {
        return("Summary")
      }
      if (ext %in% c("png", "jpg", "jpeg", "webp", "gif", "svg")) {
        return("Figures")
      }
      if (ext %in% c("csv", "tsv", "xls", "xlsx")) {
        return("Tables")
      }
      if (ext %in% c("html", "htm", "pdf", "docx", "pptx")) {
        return("Reports")
      }
      "Technical files"
    }

    first_component <- vapply(
      strsplit(rel, "/", fixed = TRUE),
      function(parts) {
        if (length(parts) > 1L) parts[1L] else "Run summary"
      },
      character(1)
    )
    analyses <- split(seq_along(files), first_component)

    category_order <- c(
      "Summary", "Figures", "Tables", "Reports", "Technical files"
    )

    tags$div(
      class = "result-tree",
      tags$p(
        class = "small text-muted mb-2",
        "Expand an analysis and select one result to display it."
      ),
      lapply(seq_along(analyses), function(i) {
        idx <- analyses[[i]]
        analysis_name <- names(analyses)[i]
        categories <- split(idx, vapply(files[idx], classify_result, character(1)))
        categories <- categories[intersect(category_order, names(categories))]

        tags$details(
          open = i == 1L,
          tags$summary(
            class = "fw-semibold py-1",
            icon("folder"),
            tags$span(
              class = "ms-1",
              gsub("_", " ", analysis_name)
            )
          ),
          tags$div(
            class = "ps-2 pb-2",
            lapply(names(categories), function(category_name) {
              category_idx <- categories[[category_name]]
              tags$details(
                open = identical(category_name, "Summary"),
                tags$summary(
                  class = "small fw-semibold py-1",
                  category_name,
                  tags$span(
                    class = "badge text-bg-light ms-1",
                    length(category_idx)
                  )
                ),
                tags$div(
                  class = "small ps-2 pb-1",
                  lapply(category_idx, function(j) {
                    tree_file_link(files[j], rel[j])
                  })
                )
              )
            })
          )
        )
      })
    )
  })

  observeEvent(input$result_tree_selection, {
    req(state$output_dir, input$result_tree_selection)

    relative_path <- as.character(input$result_tree_selection)[1L]
    candidate <- normalizePath(
      file.path(state$output_dir, relative_path),
      winslash = "/",
      mustWork = FALSE
    )
    root <- normalizePath(
      state$output_dir,
      winslash = "/",
      mustWork = TRUE
    )

    # Prevent selecting paths outside the current run directory.
    if (!startsWith(tolower(candidate), paste0(tolower(root), "/")) ||
        !file.exists(candidate)) {
      showNotification(
        "The selected result is not available in the current run.",
        type = "error"
      )
      return(invisible(NULL))
    }

    select_result_file(candidate, notify_errors = TRUE)

  }, ignoreInit = TRUE)

  selected_result_url <- reactive({
    req(state$selected_result_path, state$output_dir, state$resource_alias)

    output_root <- normalizePath(
      state$output_dir,
      winslash = "/",
      mustWork = TRUE
    )
    selected <- normalizePath(
      state$selected_result_path,
      winslash = "/",
      mustWork = TRUE
    )
    relative <- substring(selected, nchar(output_root) + 2L)
    encoded <- paste(
      vapply(
        strsplit(relative, "/", fixed = TRUE)[[1]],
        utils::URLencode,
        character(1),
        reserved = TRUE
      ),
      collapse = "/"
    )

    paste0(
      state$resource_alias,
      "/",
      encoded,
      "?v=",
      as.integer(file.info(selected)$mtime)
    )
  })

  output$selected_result_ui <- renderUI({
    if (is.null(state$selected_result_path) ||
        !file.exists(state$selected_result_path)) {
      return(tags$div(
        class = "alert alert-info",
        "Select a file from the result tree to preview it."
      ))
    }

    path <- state$selected_result_path
    extension <- tolower(tools::file_ext(path))
    title <- basename(path)

    header <- tags$div(
      class = "d-flex justify-content-between align-items-center mb-3",
      tags$strong(class = "text-break", title),
      tags$a(
        class = "btn btn-outline-secondary btn-sm",
        href = selected_result_url(),
        target = "_blank",
        icon("up-right-from-square"),
        " Open"
      )
    )

    if (extension %in% c("csv", "tsv", "txt", "xls", "xlsx")) {
      return(header)
    }

    if (extension %in% c("png", "jpg", "jpeg", "webp", "gif", "svg")) {
      return(tags$div(
        header,
        tags$img(
          src = selected_result_url(),
          style = paste(
            "display:block; width:100%; height:auto;",
            "max-height:72vh; object-fit:contain; background:white;"
          )
        )
      ))
    }

    if (extension %in% c("html", "htm", "pdf")) {
      return(tags$div(
        header,
        tags$iframe(
          src = selected_result_url(),
          style = "width:100%; height:72vh; border:1px solid #dee2e6;"
        )
      ))
    }

    if (extension %in% c("json", "yaml", "yml", "log", "md")) {
      text <- tryCatch(
        paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
        error = function(e) conditionMessage(e)
      )
      return(tags$div(
        header,
        tags$pre(
          style = "max-height:72vh; overflow:auto; white-space:pre-wrap;",
          text
        )
      ))
    }

    tags$div(
      header,
      class = "alert alert-light border",
      "Preview is not available for this file type. Use Open to inspect the file."
    )
  })

  output$selected_result_table <- renderDT({
    req(state$selected_result_path)

    extension <- tolower(tools::file_ext(state$selected_result_path))
    validate(
      need(
        extension %in% c("csv", "tsv", "txt", "xls", "xlsx"),
        ""
      ),
      need(
        !is.null(state$selected_result_data),
        "The selected table could not be previewed."
      )
    )

    DT::datatable(
      prepare_dt_table(
        as.data.frame(state$selected_result_data, check.names = FALSE),
        max_rows = 10000
      ),
      extensions = c("Buttons", "Scroller"),
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel"),
        scrollX = TRUE,
        scrollY = "62vh",
        scroller = TRUE,
        deferRender = TRUE,
        pageLength = 25
      ),
      filter = "top",
      rownames = FALSE
    )
  }, server = FALSE)

  # Results are rendered exclusively on demand through the tree selector.
  # Global figure/table galleries were intentionally removed to avoid eager
  # rendering and excessive reactive invalidation.

  output$cache_status_ui <- renderUI({
    cache_revision()

    status <- aaa_cache_status(project_root)

    tags$div(
      class = "table-responsive mb-3",
      tags$table(
        class = "table table-sm align-middle",
        tags$tbody(
          lapply(
            seq_len(nrow(status)),
            function(i) {
              tags$tr(
                tags$th(
                  status$Metric[i],
                  scope = "row"
                ),
                tags$td(
                  status$Value[i]
                )
              )
            }
          )
        )
      )
    )
  })

  observeEvent(input$refresh_cache_status, {
    cache_revision(
      cache_revision() + 1L
    )
  })

  observeEvent(
    input$request_clear_cache,
    {
      if (isTRUE(state$running)) {
        showNotification(
          "Stop the active workflow before clearing the cache.",
          type = "warning",
          duration = 6
        )

        return(invisible(NULL))
      }

      showModal(
        modalDialog(
          title = "Clear cache?",
          tags$p(
            paste(
              "This removes downloaded GFF annotations,",
              "NCBI references and cached gene-search results."
            )
          ),
          tags$p(
            class = "text-danger fw-semibold",
            "Future analyses will need to download and evaluate them again."
          ),
          easyClose = TRUE,
          footer = tagList(
            modalButton("Cancel"),
            actionButton(
              "confirm_clear_cache",
              "Clear cache",
              class = "btn-danger"
            )
          )
        )
      )
    }
  )

  observeEvent(
    input$confirm_clear_cache,
    {
      removeModal()

      results_root <- file.path(
        project_root,
        "Results"
      )

      paths <- aaa_create_project_structure(
        project_dir = results_root,
        analysis_name = NULL,
        results_root = results_root
      )

      unlink(
        paths$cache,
        recursive = TRUE,
        force = TRUE
      )

      cache_revision(
        cache_revision() + 1L
      )

      showNotification(
        "The cache was cleared.",
        type = "warning",
        duration = 6
      )
    }
  )


  output$download_results <- downloadHandler(
    filename = function() {
      paste0(
        "Triple_A_results_",
        format(Sys.Date(), "%Y%m%d"),
        ".zip"
      )
    },
    content = function(file) {
      req(dir.exists(state$output_dir))

      old <- setwd(dirname(state$output_dir))
      on.exit(setwd(old), add = TRUE)

      utils::zip(
        zipfile = file,
        files = basename(state$output_dir),
        flags = "-r9X"
      )
    }
  )
  output$download_report <- downloadHandler(
    filename = function() paste0("Triple_A_report_", format(Sys.Date(), "%Y%m%d"), ".html"),
    content = function(file) {
      req(state$results)
      report <- state$results$metadata$report %||% file.path(state$output_dir, "Triple_A_report.html")
      validate(need(file.exists(report), "No report has been generated for the current run."))
      file.copy(report, file, overwrite = TRUE)
    }
  )

}
