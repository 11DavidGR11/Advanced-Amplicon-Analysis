analysis_info_box <- function(question, description) {
  tags$div(
    class = "alert alert-light border small mb-3",
    tags$strong(question),
    tags$br(),
    description
  )
}

analysis_option <- function(id, label, description) {
  tags$div(
    class = "analysis-option border rounded p-3 mb-2",
    checkboxInput(id, label, FALSE),
    tags$div(class = "small text-muted ms-4", description)
  )
}

ui <- page_sidebar(
  shinyjs::useShinyjs(),
  tags$style(HTML("
    html, body { min-height: 100%; overflow-x: hidden; }
    .bslib-page-fill, .bslib-page-fill > body { height: auto !important; min-height: 100vh; }
    .triple-a-analysis-disabled { opacity: 0.58; }
    .triple-a-analysis-disabled input,
    .triple-a-analysis-disabled select,
    .triple-a-analysis-disabled textarea,
    .triple-a-analysis-disabled button { cursor: not-allowed !important; }
    .analysis-question { border-left: 4px solid #5B3F7A; }
    .card { min-width: 0; height: auto !important; }
    .card-body { min-height: 0 !important; overflow: visible !important; }
    .bslib-grid { align-items: start; }
    .accordion { --bs-accordion-btn-padding-y: .8rem; }
    .accordion-body { padding: 1rem; overflow: visible; }
    .analysis-option { background: #fff; transition: background-color .15s ease; }
    .analysis-option:hover { background: #f8f6fa; }
    .analysis-option .form-check { margin-bottom: .15rem; }
    .metadata-layout .card { width: 100%; }
    @media (max-width: 991.98px) {
      .bslib-sidebar-layout > .main { min-width: 0; }
    }
  ")),

  title = "Biological Analysis (Triple A)",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#5B3F7A"),
  sidebar = sidebar(
    width = 320,
    aaa_help_button(),
    tags$hr(),
    tags$h5("Triple_A workflow"),
    tags$p(class="small text-muted",
      "1. Import and describe the data. 2. Select compatible analyses. 3. Validate and run. 4. Inspect or download the results."),
    uiOutput("dependency_status_ui"),
    selectInput(
      "progress_verbosity", "Progress detail",
      choices = c(
        "Standard — concise user messages" = "standard",
        "Detailed — one message per taxon/genome" = "detailed",
        "Developer — include every gene search" = "developer"
      ), selected = "standard"
    ),
    tags$hr(),
    tags$div(class="small text-muted",
      tags$strong(TRIPLE_A_PROJECT_SUBTITLE),
      tags$br(), "Developed by ", paste(TRIPLE_A_AUTHORS$name, collapse=", "),
      tags$br(), TRIPLE_A_LICENSE, " License"
    )
  ),

  navset_card_tab(
    id = "main_tabs",

    nav_panel(
      "1. Data and metadata",
      tags$div(class="alert alert-info",
        tags$strong("Define the analytical input before selecting methods. "),
        "The abundance table is required. Sample metadata are optional and unlock comparisons, environmental models and multivariable taxon associations."
      ),
      layout_columns(
        class = "metadata-layout",
        card(
          card_header("Abundance and taxonomy"),
          tags$div(class="p-3",
            fileInput("input_file", "Amplicon abundance table",
              accept=c(".csv", ".tsv", ".txt", ".xls", ".xlsx"), width="100%"),
            uiOutput("abundance_import_controls_ui"),
            uiOutput("treatments_ui")
          )
        ),
        card(
          card_header("Input preview"),
          tags$div(class="p-2", DTOutput("preview"))
        ),
        col_widths=c(5,7)
      ),
      uiOutput("environmental_tab_ui"),
      layout_columns(
        card(card_header("Current data definition"), verbatimTextOutput("selection_summary")),
        card(card_header("Input validation"), uiOutput("validation_status_ui")),
        col_widths=c(6,6)
      )
    ),

    nav_panel(
      "2. Analyses and parameters",
      tags$div(class="alert alert-light border",
        tags$strong("Analyses are unlocked from the declared data roles. "),
        "Unavailable methods remain disabled and explain which input is missing."
      ),
      tags$div(style="display:none;", checkboxGroupInput("analyses", NULL, choices=analysis_choices, selected=character(0))),
      accordion(
        id="analysis_accordion", open=FALSE,
        accordion_panel(
          "Community composition — Which taxa are present and most abundant?",
          analysis_info_box("Biological question", "Which microorganisms compose the community, in what proportions, and which taxa dominate the samples?"),
          uiOutput("top_abundance_availability_ui"),
          tags$fieldset(id="top_abundance_controls", style="border:0;padding:0;margin:0;min-width:0;",
            analysis_option("use_top_abundance", "Top abundance", "Ranks the most abundant taxa and generates composition, distribution, heatmap and lollipop visualisations."),
            numericInput("top_n", "Number of top taxa", value=20, min=1),
            checkboxGroupInput("figures_top_abundance", "Figures to generate", choices=figure_choices_for("top_abundance"), selected=unname(figure_choices_for("top_abundance")))
          )
        ),
        accordion_panel(
          "Diversity — How diverse is each sample and how different are samples?",
          analysis_info_box("Biological question", "How many taxa occur within samples, how evenly distributed are they, and how dissimilar are communities between samples?"),
          uiOutput("community_structure_availability_ui"),
          tags$fieldset(id="diversity_controls", style="border:0;padding:0;margin:0;min-width:0;",
            analysis_option("analysis_alpha", "Alpha diversity", "Observed richness, Shannon, Simpson and inverse Simpson indices for each sample."),
            analysis_option("analysis_beta", "Beta diversity", "Pairwise community dissimilarity using Bray–Curtis or presence/absence Jaccard distances."),
            selectInput("community_distance", "Beta-diversity distance", choices=c("Bray-Curtis"="bray", "Jaccard (presence/absence)"="jaccard"), selected="bray")
          )
        ),
        accordion_panel(
          "Community structure — What are the main patterns among samples?",
          analysis_info_box("Biological question", "Do samples form gradients, clusters or recurring compositional patterns without imposing an environmental model?"),
          uiOutput("community_structure_availability_ui"),
          tags$fieldset(id="community_structure_controls", style="border:0;padding:0;margin:0;min-width:0;",
            analysis_option("analysis_pca", "PCA", "Linear ordination of transformed abundances; emphasizes dominant variance gradients."),
            analysis_option("analysis_pcoa", "PCoA", "Ordination of the selected beta-diversity distance matrix."),
            analysis_option("analysis_nmds", "NMDS", "Rank-based ordination suitable for ecological dissimilarities; reports stress as a fit diagnostic."),
            analysis_option("analysis_clustering", "Hierarchical clustering", "Groups samples according to the selected community dissimilarity."),
            selectInput("community_transformation", "Ordination transformation", choices=c("Hellinger (recommended)"="hellinger", "Relative abundance"="relative", "Log1p relative abundance"="log1p"), selected="hellinger"),
            checkboxInput("ordination_labels", "Show sample labels", FALSE),
            tags$div(style="display:none;", checkboxInput("use_community_structure", "Run community structure", FALSE)),
            checkboxGroupInput("figures_community_structure", "Figures to generate", choices=figure_choices_for_panel("community_structure", exclude=community_comparison_figure_ids), selected=unname(figure_choices_for_panel("community_structure", exclude=community_comparison_figure_ids)))
          )
        ),
        accordion_panel(
          "Community comparison — Are whole communities different between groups?",
          analysis_info_box("Biological question", "Does community composition differ between experimental groups, and are significant results compatible with similar within-group dispersion?"),
          uiOutput("community_comparison_availability_ui"),
          tags$fieldset(id="community_comparison_controls", style="border:0;padding:0;margin:0;min-width:0;",
            analysis_option("analysis_permanova", "PERMANOVA", "Permutation test of group or multivariable effects on a distance matrix; reports effect size as R²."),
            analysis_option("analysis_anosim", "ANOSIM", "Rank-based test comparing between-group and within-group dissimilarities."),
            analysis_option("analysis_pairwise_permanova", "Pairwise PERMANOVA", "Post-hoc PERMANOVA for every pair of levels, with multiple-testing correction."),
            analysis_option("analysis_beta_dispersion", "Beta dispersion", "Tests homogeneity of multivariate dispersion, an essential companion diagnostic for PERMANOVA."),
            numericInput("permutations", "Permutations", value=999, min=99, step=100),
            numericInput("community_alpha", "Significance threshold", value=.05, min=.001, max=.99, step=.01),
            checkboxGroupInput("figures_community_comparison", "Figures to generate", choices=figure_choices_for_panel("community_structure", ids=community_comparison_figure_ids), selected=unname(figure_choices_for_panel("community_structure", ids=community_comparison_figure_ids)))
          )
        ),
        accordion_panel(
          "Community discrimination — Can composition classify predefined groups?",
          analysis_info_box("Biological question", "Can the amplicon abundance profile distinguish or predict the experimental group assigned to each sample?"),
          uiOutput("plsda_availability_ui"),
          analysis_option("use_plsda", "PLS-DA", "Supervised multivariate classification that constructs components maximizing separation among declared groups."),
          tags$fieldset(id="plsda_controls", style="border:0;padding:0;margin:0;min-width:0;",
            numericInput("plsda_components", "Components", 2, min=1, max=5),
            numericInput("plsda_folds", "Cross-validation folds", 5, min=2, max=10),
            checkboxGroupInput("figures_plsda", "Figures", choices=figure_choices_for("plsda"), selected=unname(figure_choices_for("plsda")))
          ),
          analysis_option("use_splsda", "sPLS-DA", "Sparse supervised classification that selects a compact microbial signature while separating declared groups."),
          tags$fieldset(id="splsda_controls", style="border:0;padding:0;margin:0;min-width:0;",
            numericInput("splsda_components", "Maximum components", 2, min=1, max=5),
            numericInput("splsda_folds", "Cross-validation folds", 5, min=2, max=10),
            numericInput("splsda_repeats", "Repeated validations", 10, min=1, max=100),
            checkboxInput("splsda_tune", "Tune the number of selected taxa automatically", TRUE),
            textInput("splsda_keepx", "Candidate taxa per component", "5, 10, 20"),
            checkboxGroupInput("figures_splsda", "Figures", choices=figure_choices_for("splsda"), selected=unname(figure_choices_for("splsda")))
          )
        ),
        accordion_panel(
          "Environmental analysis — Which measured variables explain community composition?",
          analysis_info_box("Biological question", "Which environmental gradients are associated with the community, how much variation do they explain, and what remains after controlling for confounders?"),
          uiOutput("constrained_ordination_availability_ui"),
          tags$fieldset(id="rda_controls", style="border:0;padding:0;margin:0;min-width:0;",
            analysis_option("analysis_envfit", "envfit", "Fits environmental vectors or factor centroids onto PCA, PCoA and NMDS using permutation tests."),
            uiOutput("rda_availability_ui"),
            analysis_option("use_rda", "RDA", "Constrained linear ordination of Hellinger-transformed abundances using continuous environmental predictors."),
            analysis_option("analysis_partial_rda", "Partial RDA", "RDA after conditioning on declared blocking or confounding variables."),
            analysis_option("analysis_dbrda", "dbRDA", "Distance-based constrained ordination for Bray–Curtis or Jaccard dissimilarities."),
            analysis_option("analysis_partial_dbrda", "Partial dbRDA", "dbRDA after conditioning on blocking or confounding variables."),
            uiOutput("variance_partitioning_availability_ui"),
            analysis_option("analysis_varpart", "Variance partitioning", "Partitions explained variation among environmental and experimental predictor sets."),
            selectInput("environmental_distance", "dbRDA distance", choices=c("Bray-Curtis"="bray","Jaccard"="jaccard")),
            checkboxGroupInput("figures_rda", "RDA figures", choices=figure_choices_for("rda"), selected=unname(figure_choices_for("rda")))
          )
        ),
        accordion_panel(
          "Taxon associations — Which taxa differ or track metadata variables?",
          analysis_info_box("Biological question", "Which taxa are more abundant in one group, or change along environmental gradients after adjustment for other variables?"),
          uiOutput("taxon_association_availability_ui"),
          tags$fieldset(id="differential_abundance_controls", style="border:0;padding:0;margin:0;min-width:0;",
            analysis_option("analysis_ancombc2", "ANCOM-BC2", "Compositional differential-abundance analysis with bias correction, covariates and multi-group contrasts."),
            analysis_option("analysis_maaslin", "MaAsLin", "Multivariable taxon–metadata association models for categorical and continuous predictors."),
            sliderInput("min_prevalence", "Minimum prevalence", min=0,max=1,value=.20,step=.05),
            numericInput("alpha", "Adjusted P-value threshold", value=.05,min=.001,max=.99),
            tags$div(style="display:none;", checkboxInput("use_differential_abundance", "Legacy differential abundance", FALSE)),
            tags$div(style="display:none;", selectInput("method", "Legacy test", choices=c("Wilcoxon"="wilcox", "Student t-test"="t_test"))),
            tags$div(style="display:none;", numericInput("min_mean_abundance", "Minimum mean abundance", .01)),
            tags$div(style="display:none;", numericInput("log2fc", "log2FC", 1)),
            tags$div(style="display:none;", numericInput("max_labels", "Labels", 10)),
            uiOutput("pairwise_comparisons_ui"),
            checkboxGroupInput("figures_differential_abundance", "Figures", choices=figure_choices_for("differential_abundance"), selected=unname(figure_choices_for("differential_abundance")))
          )
        ),
        accordion_panel(
          "Functional analysis — What biological capabilities are represented?",
          analysis_info_box("Biological question", "Which curated biological functions are potentially represented, how abundant are their contributing taxa, and which functions differ between conditions?"),
          uiOutput("functional_potential_availability_ui"),
          tags$fieldset(id="functional_potential_controls", style="border:0;padding:0;margin:0;min-width:0;",
            analysis_option("use_functional_potential", "Functional potential", "Infers curated functional potential from taxonomic assignments and reference-genome/GFF evidence.")
          ),
          tags$div(style="display:none;", checkboxGroupInput("functional_functions", NULL, choices=function_choices, selected=c("Methanogenesis","Methanogenesis_subtype","Methanotrophy"))),
          selectizeInput("function_search", "Find a biological function", choices=c("Select a function..."="", function_choices), multiple=FALSE, options=list(placeholder="Type a function, pathway or process...")),
          layout_columns(actionButton("functions_all", "Select all", class="btn-outline-primary btn-sm w-100"), actionButton("functions_none", "Clear", class="btn-outline-secondary btn-sm w-100"), selectInput("function_preset", "Preset", choices=c("Custom"="", names(function_presets))), col_widths=c(3,3,6), gap="0.4rem"),
          function_group_accordion,
          uiOutput("function_selection_summary"), uiOutput("function_details_ui"),
          checkboxGroupInput("figures_functional_potential", "Figures", choices=figure_choices_for("functional_potential"), selected=unname(figure_choices_for("functional_potential"))),
          tags$hr(),
          analysis_option("use_functional_abundance", "Functional abundance", "Summarises the relative abundance of taxa contributing to selected functions."),
          uiOutput("functional_abundance_availability_ui"),
          tags$fieldset(id="functional_abundance_controls", style="border:0;padding:0;margin:0;min-width:0;",
            numericInput("top_contributors", "Top contributors per function", value=5,min=1),
            checkboxGroupInput("figures_functional_abundance", "Figures", choices=figure_choices_for("functional_abundance"), selected=unname(figure_choices_for("functional_abundance")))
          ),
          analysis_option("analysis_differential_functions", "Differential functions", "Compares inferred functional profiles between declared experimental groups."),
          analysis_option("analysis_functional_enrichment", "Functional enrichment", "Tests whether selected functions are over-represented among significant or high-abundance features.")
        )
      ),
      tags$div(style="display:none;", uiOutput("outputs_ui")),
      layout_columns(
        actionButton("run", "Run selected analyses", icon=icon("play"), class="btn-primary w-100"),
        actionButton("stop_run", "Stop", icon=icon("stop"), class="btn-danger w-100"),
        col_widths=c(9,3), gap="0.5rem"
      )
    ),

    nav_panel(
      "Methods",
      card(
        card_header("Methods, parameters and interpretation"),
        tags$div(
          class = "alert alert-light border",
          "This panel describes the selected analyses using the current parameters. PCA, PCoA and NMDS are exploratory and do not produce significance P-values; PERMANOVA and differential abundance use explicit statistical criteria."
        ),
        uiOutput("selected_methods_ui")
      )
    ),

    nav_panel(
      "Run",
      layout_columns(
        card(
          card_header("Current stage"),
          tags$div(
            class = "p-3",
            tags$div(
              id = "triple_a_progress_stage",
              class = "fs-4 fw-semibold",
              "Waiting to start"
            ),
            tags$div(
              id = "triple_a_progress_detail",
              class = "text-muted mt-2",
              "Select the analyses and press Run."
            )
          )
        ),
        card(
          card_header("Execution status"),
          tags$div(
            class = "p-3",
            tags$div(
              tags$strong("Progress: "),
              tags$span(
                id = "triple_a_progress_percent",
                "0%"
              )
            ),
            tags$div(
              class = "mt-2",
              tags$strong("Elapsed time: "),
              tags$span(
                id = "triple_a_progress_elapsed",
                "00:00:00"
              )
            )
          )
        ),
        col_widths = c(8, 4)
      ),

      conditionalPanel(
        condition = "input.progress_verbosity == 'standard'",
        card(
          card_header(
            "Potential metabolic pathways analysis"
          ),
          tags$div(
            class = "p-3",
            tags$div(
              class = "mb-4",
              tags$div(
                class = "d-flex justify-content-between",
                tags$strong("Reference annotations"),
                tags$span(
                  id = "triple_a_genome_count",
                  "0 / 0"
                )
              ),
              tags$div(
                class = "progress mt-2",
                style = "height: 18px;",
                tags$div(
                  id = "triple_a_genome_bar",
                  class = paste(
                    "progress-bar",
                    "progress-bar-striped",
                    "progress-bar-animated"
                  ),
                  role = "progressbar",
                  style = paste0(
                    "width: 0%;",
                    "background-color: ",
                    "#5B3F7A;"
                  )
                )
              )
            ),

            tags$div(
              class = "mb-4",
              tags$div(
                class = "d-flex justify-content-between",
                tags$strong("Genes"),
                tags$span(
                  id = "triple_a_gene_count",
                  "0 / 0"
                )
              ),
              tags$div(
                class = "progress mt-2",
                style = "height: 18px;",
                tags$div(
                  id = "triple_a_gene_bar",
                  class = paste(
                    "progress-bar",
                    "progress-bar-striped",
                    "progress-bar-animated"
                  ),
                  role = "progressbar",
                  style = paste0(
                    "width: 0%;",
                    "background-color: ",
                    "#2A9D8F;"
                  )
                )
              )
            ),

            layout_columns(
              tags$div(
                class = "border rounded p-3 h-100",
                tags$div(
                  class = "small text-muted",
                  "Current genome"
                ),
                tags$div(
                  id = "triple_a_current_genome",
                  class = "fw-semibold mt-1",
                  "—"
                )
              ),
              tags$div(
                class = "border rounded p-3 h-100",
                tags$div(
                  class = "small text-muted",
                  "Current pathway"
                ),
                tags$div(
                  id = "triple_a_current_pathway",
                  class = "fw-semibold mt-1",
                  "—"
                )
              ),
              col_widths = c(6, 6)
            )
          )
        )
      ),

      conditionalPanel(
        condition = "input.progress_verbosity != 'standard'",
        card(
          card_header("Live execution log"),
          tags$pre(
            id = "triple_a_progress_log",
            style = paste(
              "min-height: 420px;",
              "max-height: 620px;",
              "overflow-y: auto;",
              "white-space: pre-wrap;",
              "background: #111827;",
              "color: #f9fafb;",
              "padding: 1rem;",
              "border-radius: 0.4rem;"
            ),
            "No workflow has been started."
          )
        )
      )
    ),

    nav_panel(
      "Results",
      layout_columns(
        card(
          card_header("Result tree"),
          uiOutput("result_tree_ui")
        ),
        card(
          card_header("Selected result"),
          tags$div(
            class = "p-2",
            uiOutput("selected_result_ui"),
            DTOutput("selected_result_table")
          )
        ),
        col_widths = c(3, 9)
      )
    ),

    nav_panel(
      "Downloads",
      card(
        card_header("Exports and reports"),
        tags$div(class = "p-3 d-flex gap-2 flex-wrap",
          downloadButton("download_results", "Download current run ZIP"),
          downloadButton("download_report", "Download report")
        )
      )
    ),

    nav_panel(
      "Cache",
      card(
        card_header("Shared cache"),
        tags$div(
          class = "p-3",
          tags$p(
            paste(
              "Triple_A reuses compatible results automatically.",
              "The cache is created only when needed and can be cleared safely from this panel."
            )
          ),
          uiOutput("cache_status_ui"),
          layout_columns(
            actionButton(
              "refresh_cache_status",
              "Refresh status",
              icon = icon("rotate"),
              class = "btn-outline-primary"
            ),
            actionButton(
              "request_clear_cache",
              "Clear cache",
              icon = icon("trash"),
              class = "btn-outline-danger"
            ),
            col_widths = c(6, 6)
          )
        )
      )
    ),

    nav_panel(
      "History",
      card(
        card_header("Execution history"),
        tags$div(class = "p-2 d-flex gap-2 flex-wrap",
          actionButton("refresh_history", "Refresh history", icon = icon("rotate")),
          actionButton("open_history_run", "Open selected results", icon = icon("folder-open")),
          actionButton("reload_history_run", "Reload selected project", icon = icon("rotate-left"))
        ),
        tags$p(class="small text-muted px-2", "Select one run. Open displays its saved results; Reload restores the dataset and analysis configuration when a run snapshot is available."),
        DTOutput("history_table")
      )
    ),

    nav_panel(
      "About",
      card(
        card_header(
          paste0(
            TRIPLE_A_PROJECT_TITLE,
            " — ",
            TRIPLE_A_PROJECT_SUBTITLE
          )
        ),
        tags$div(
          class = "p-3",
          tags$h5("Author and development"),
          tags$div(
            lapply(
              seq_len(nrow(TRIPLE_A_AUTHORS)),
              function(i) {
                tags$div(
                  class = "mb-3",
                  tags$strong(
                    TRIPLE_A_AUTHORS$name[i]
                  ),
                  tags$br(),
                  TRIPLE_A_AUTHORS$role[i],
                  tags$br(),
                  TRIPLE_A_AUTHORS$institution[i]
                )
              }
            )
          ),
          tags$p(
            tags$strong("License: "),
            TRIPLE_A_LICENSE
          ),
          tags$p(
            tags$strong("Suggested citation: "),
            aaa_project_metadata()$citation
          )
        )
      )
    )
  )
)

