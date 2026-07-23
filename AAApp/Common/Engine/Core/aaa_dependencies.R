# =============================================================================
# Triple_A dependency registry
#
# Dependencies are divided into:
# - launcher: required before the Shiny application can open;
# - core: required for the standard analytical engine;
# - optional: required only by selected analyses or input formats.
# =============================================================================

TRIPLE_A_DEPENDENCIES <- list(
  launcher = c(
    "shiny",
    "bslib",
    "DT",
    "shinyjs"
  ),
  core = c(
    "dplyr",
    "tidyr",
    "tibble",
    "readr",
    "readxl",
    "stringr",
    "ggplot2",
    "ggrepel",
    "openxlsx",
    "jsonlite",
    "callr",
    "DBI",
    "RSQLite",
    "base64enc"
  ),
  analyses = list(
    functional_potential = c(
      "rentrez",
      "pheatmap"
    ),
    top_abundance = c(
      "pheatmap"
    ),
    differential_abundance = character(),
    functional_abundance = c(
      "pheatmap"
    ),
    community_structure = c(
      "vegan"
    ),
    plsda = c(
      "pls"
    ),
    splsda = c(
      "mixOmics"
    ),
    rda = c("vegan"),
    envfit = c("vegan"),
    partial_rda = c("vegan"),
    dbrda = c("vegan"),
    partial_dbrda = c("vegan"),
    variance_partitioning = c("vegan"),
    ancombc2 = c("ANCOMBC", "TreeSummarizedExperiment", "S4Vectors"),
    maaslin = c("Maaslin2"),
    # Derived analyses that consume results already computed by the functional
    # and differential-abundance engines. They add no packages of their own
    # (base stats + core ggplot2/openxlsx), but they must still be listed here:
    # aaa_required_packages() rejects any selected analysis missing from this
    # map as "Unknown analyses", which blocked runs that selected them.
    differential_functions = character(),
    functional_enrichment = character()
  ),
  input_formats = list(
    xls = "readxl",
    xlsx = "readxl"
  )
)

aaa_required_packages <- function(
  analyses = character(),
  input_files = character(),
  include_launcher = FALSE,
  include_core = TRUE
) {
  analyses <- unique(
    analyses[
      !is.na(analyses) &
        nzchar(analyses)
    ]
  )

  unknown <- setdiff(
    analyses,
    names(TRIPLE_A_DEPENDENCIES$analyses)
  )

  if (length(unknown) > 0) {
    stop(
      "Unknown analyses when resolving dependencies: ",
      paste(unknown, collapse = ", ")
    )
  }

  packages <- character()

  if (isTRUE(include_launcher)) {
    packages <- c(
      packages,
      TRIPLE_A_DEPENDENCIES$launcher
    )
  }

  if (isTRUE(include_core)) {
    packages <- c(
      packages,
      TRIPLE_A_DEPENDENCIES$core
    )
  }

  if (length(analyses) > 0) {
    packages <- c(
      packages,
      unlist(
        TRIPLE_A_DEPENDENCIES$analyses[
          analyses
        ],
        use.names = FALSE
      )
    )
  }

  input_files <- input_files[
    !is.na(input_files) &
      nzchar(input_files)
  ]

  if (length(input_files) > 0) {
    extensions <- unique(
      tolower(
        tools::file_ext(input_files)
      )
    )

    format_packages <- unlist(
      TRIPLE_A_DEPENDENCIES$input_formats[
        intersect(
          extensions,
          names(
            TRIPLE_A_DEPENDENCIES$input_formats
          )
        )
      ],
      use.names = FALSE
    )

    packages <- c(
      packages,
      format_packages
    )
  }

  unique(packages[nzchar(packages)])
}

aaa_missing_packages <- function(
  analyses = character(),
  input_files = character(),
  include_launcher = FALSE,
  include_core = TRUE
) {
  required <- aaa_required_packages(
    analyses = analyses,
    input_files = input_files,
    include_launcher = include_launcher,
    include_core = include_core
  )

  required[
    !vapply(
      required,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
}

aaa_dependency_status <- function(
  analyses = character(),
  input_files = character()
) {
  required <- aaa_required_packages(
    analyses = analyses,
    input_files = input_files,
    include_core = TRUE
  )

  installed <- vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )

  data.frame(
    Package = required,
    Installed = unname(installed),
    Required_by = vapply(
      required,
      function(package) {
        relevant <- names(
          Filter(
            function(packages) {
              package %in% packages
            },
            TRIPLE_A_DEPENDENCIES$analyses[
              analyses
            ]
          )
        )

        labels <- aaa_analysis_catalogue()$Name[
          match(
            relevant,
            aaa_analysis_catalogue()$ID
          )
        ]

        labels <- labels[
          !is.na(labels)
        ]

        if (
          package %in%
            TRIPLE_A_DEPENDENCIES$core
        ) {
          labels <- c(
            "Triple_A core",
            labels
          )
        }

        paste(
          unique(labels),
          collapse = ", "
        )
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

aaa_check_packages <- function(
  analyses = character(),
  input_files = character(),
  install_missing = FALSE,
  include_launcher = FALSE,
  include_core = TRUE
) {
  required <- aaa_required_packages(
    analyses = analyses,
    input_files = input_files,
    include_launcher = include_launcher,
    include_core = include_core
  )

  missing <- aaa_missing_packages(
    analyses = analyses,
    input_files = input_files,
    include_launcher = include_launcher,
    include_core = include_core
  )

  if (length(missing) > 0 &&
    isTRUE(install_missing)) {
    utils::install.packages(missing)

    missing <- aaa_missing_packages(
      analyses = analyses,
      input_files = input_files,
      include_launcher = include_launcher,
      include_core = include_core
    )
  }

  if (length(missing) > 0) {
    stop(
      "Missing R packages required for the selected workflow: ",
      paste(missing, collapse = ", "),
      ". Install them with install.packages(c(",
      paste(
        sprintf('"%s"', missing),
        collapse = ", "
      ),
      "))."
    )
  }

  invisible(required)
}
