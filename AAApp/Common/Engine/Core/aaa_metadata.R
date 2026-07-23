# =============================================================================
# Triple_A project identity and authorship
#
# Edit this file to update authors, institutional affiliation, contact details,
# repository URL or DOI. Shiny and the generated metadata read these values.
# =============================================================================

TRIPLE_A_PROJECT_TITLE <- "Advanced Amplicon Analysis"
TRIPLE_A_PROJECT_SUBTITLE <- "Triple_A"
TRIPLE_A_LICENSE <- "MIT"

TRIPLE_A_AUTHORS <- data.frame(
  name = "David Garrido Rodríguez",
  institution = "Universidad de Valladolid, Institute of Sustainable Processes (ISP)",
  email = "david.garrido23@uva.es",
  orcid = "0000-0001-5180-0006",
  role = "Creator and lead developer",
  stringsAsFactors = FALSE
)

TRIPLE_A_REPOSITORY <- ""
TRIPLE_A_DOI <- ""

aaa_project_metadata <- function() {
  authors <- TRIPLE_A_AUTHORS

  authors$display <- vapply(
    seq_len(nrow(authors)),
    function(i) {
      values <- c(
        authors$name[i],
        authors$institution[i],
        if (nzchar(authors$orcid[i])) {
          paste0("ORCID: ", authors$orcid[i])
        } else {
          character()
        },
        if (nzchar(authors$email[i])) {
          authors$email[i]
        } else {
          character()
        }
      )

      paste(
        values[nzchar(values)],
        collapse = " | "
      )
    },
    character(1)
  )

  list(
    project_name = TRIPLE_A_NAME,
    short_name = TRIPLE_A_SHORT_NAME,
    title = TRIPLE_A_PROJECT_TITLE,
    subtitle = TRIPLE_A_PROJECT_SUBTITLE,
    license = TRIPLE_A_LICENSE,
    authors = authors,
    repository = TRIPLE_A_REPOSITORY,
    doi = TRIPLE_A_DOI,
    citation = paste0(
      paste(authors$name, collapse = ", "),
      ". ",
      TRIPLE_A_PROJECT_TITLE,
      " (",
      TRIPLE_A_PROJECT_SUBTITLE,
      ")."
    )
  )
}

aaa_author_display <- function(separator = " · ") {
  paste(
    aaa_project_metadata()$authors$display,
    collapse = separator
  )
}
