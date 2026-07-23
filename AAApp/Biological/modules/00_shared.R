
required_background_functions <- c(
  "aaa_start_background_workflow",
  "aaa_stop_background_workflow",
  "aaa_cleanup_background_workflow"
)

missing_background_functions <- required_background_functions[
  !vapply(
    required_background_functions,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )
]

if (length(missing_background_functions) > 0) {
  stop(
    paste0(
      "Triple_A could not load the background workflow controller. ",
      "Missing function(s): ",
      paste(
        missing_background_functions,
        collapse = ", "
      ),
      ". Confirm that AAApp/Common/Engine/Core/aaa_background_runner.R exists ",
      "and that AAApp/Common/Engine/Triple_A.R includes it in triple_a_load()."
    ),
    call. = FALSE
  )
}

function_catalogue <- triple_a_list_functions()
function_choices <- stats::setNames(
  function_catalogue$ID,
  function_catalogue$Name
)
function_groups <- split(
  stats::setNames(function_catalogue$ID, function_catalogue$Name),
  function_catalogue$Category
)

function_presets <- list(
  `Methane and biogas` = intersect(names(biological_function_registry), c(
    "Methanogenesis","Methanogenesis_subtype","Hydrogenotrophic_methanogenesis",
    "Acetoclastic_methanogenesis","Methylotrophic_methanogenesis",
    "Methanotrophy","Homoacetogenesis","Carboxydotrophic_hydrogenogenesis")),
  `Wastewater treatment` = intersect(names(biological_function_registry), c(
    "Nitrification","Denitrification","DNRA","Nitrous_oxide_reduction",
    "Anammox","Polyphosphate_metabolism","Ureolysis","Sulfate_reduction",
    "Sulfide_oxidation")),
  `CO, CO2 and hydrogen` = intersect(names(biological_function_registry), c(
    "Homoacetogenesis","Acetogenesis","Aerobic_CO_oxidation",
    "Carboxydotrophic_hydrogenogenesis","Fermentative_H2_production",
    "Calvin_cycle","Reverse_TCA_cycle","Hydrogenotrophic_methanogenesis"))
)

function_group_accordion <- do.call(
  bslib::accordion,
  c(
    list(id = "biological_function_groups", open = FALSE),
    lapply(names(function_groups), function(category) {
      safe_id <- paste0(
        "function_group_",
        gsub("[^A-Za-z0-9]+", "_", category)
      )
      bslib::accordion_panel(
        category,
        shiny::checkboxGroupInput(
          safe_id,
          NULL,
          choices = function_groups[[category]],
          selected = intersect(
            function_groups[[category]],
            c("Methanogenesis", "Methanogenesis_subtype", "Methanotrophy")
          )
        )
      )
    })
  )
)

analysis_choices <- stats::setNames(
  triple_a_list_analyses()$ID,
  triple_a_list_analyses()$Name
)

output_catalogue <- triple_a_list_outputs()
output_choices <- stats::setNames(
  output_catalogue$ID,
  output_catalogue$Name
)

figure_catalogue <- output_catalogue[
  output_catalogue$Type == "figure",
  ,
  drop = FALSE
]

automatic_output_ids <- output_catalogue$ID[
  output_catalogue$Type != "figure"
]

figure_choices_for <- function(module) {
  rows <- figure_catalogue$Module == module
  stats::setNames(
    figure_catalogue$ID[rows],
    figure_catalogue$Name[rows]
  )
}

# The community-structure engine computes ordination/diversity figures and
# community-comparison statistics in a single pass, but the interface presents
# them as two separate questions ("Community structure" vs "Community
# comparison"). These IDs route the comparison figures to the matching panel so
# the figure checkboxes line up with the analysis they belong to.
community_comparison_figure_ids <- c(
  "pairwise_permanova_heatmap", "beta_dispersion_boxplot",
  "anosim_plot", "permanova_variance_plot"
)

figure_choices_for_panel <- function(module, ids = NULL, exclude = NULL) {
  choices <- figure_choices_for(module)
  if (!is.null(ids)) choices <- choices[choices %in% ids]
  if (!is.null(exclude)) choices <- choices[!choices %in% exclude]
  choices
}


prepare_dt_table <- function(x, max_rows = NULL) {
  x <- as.data.frame(
    x,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  x[] <- lapply(
    x,
    function(column) {
      if (is.list(column)) {
        return(vapply(
          column,
          function(value) {
            paste(
              as.character(value),
              collapse = " | "
            )
          },
          character(1)
        ))
      }

      if (is.factor(column) ||
          inherits(column, c("Date", "POSIXct", "POSIXlt"))) {
        return(as.character(column))
      }

      if (is.numeric(column)) {
        column[!is.finite(column)] <- NA_real_
      }

      column
    }
  )

  if (!is.null(max_rows) &&
      nrow(x) > max_rows) {
    x <- utils::head(x, max_rows)
  }

  x
}


methodology_cards <- function(methods) {
  if (is.null(methods) || nrow(methods) == 0) {
    return(tags$div(
      class = "alert alert-info",
      "Select at least one analysis and output to display its methodology."
    ))
  }

  tags$div(
    lapply(
      seq_len(nrow(methods)),
      function(i) {
        row <- methods[i, , drop = FALSE]

        tags$details(
          class = "card mb-3",
          open = i == 1,
          tags$summary(
            class = "card-header fw-semibold",
            paste0(row$Analysis, " — ", row$Output)
          ),
          tags$div(
            class = "card-body",
            layout_columns(
              tags$div(
                class = "h-100",
                tags$h5("Method and criteria"),
                tags$p(tags$strong("Objective: "), row$Objective),
                tags$p(tags$strong("Data preparation: "), row$Data_preparation),
                tags$p(tags$strong("Method: "), row$Method),
                tags$p(tags$strong("Statistical test: "), row$Statistical_test),
                tags$p(tags$strong("Significance criterion: "), row$Significance_criterion),
                tags$p(tags$strong("Interpretation: "), row$Interpretation),
                tags$p(tags$strong("Assumptions and limitations: "), row$Assumptions_and_limitations),
                tags$p(tags$strong("Implementation: "), row$Implementation),
                tags$p(tags$strong("Current parameters: "), tags$code(row$Run_parameters))
              ),
              tags$div(
                class = "h-100 border-start ps-3",
                tags$h5("Illustrative result"),
                if (
                  !is.na(row$Example_image) &&
                  nzchar(row$Example_image)
                ) {
                  tags$img(
                    src = file.path(
                      "method_examples",
                      row$Example_image
                    ),
                    alt = paste0(
                      "Illustrative ",
                      row$Output
                    ),
                    style = paste(
                      "width:100%;",
                      "height:auto;",
                      "max-height:520px;",
                      "object-fit:contain;",
                      "background:white;",
                      "border:1px solid #dee2e6;",
                      "border-radius:0.4rem;",
                      "padding:0.5rem;"
                    )
                  )
                },
                tags$div(
                  class = "alert alert-warning mt-3 mb-0 small",
                  tags$strong("Illustrative example. "),
                  row$Example_caption
                )
              ),
              col_widths = c(6, 6)
            )
          )
        )
      }
    )
  )
}


