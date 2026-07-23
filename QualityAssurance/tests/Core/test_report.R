# The analysis report must be self-contained: figures embedded as base64
# data URIs (so the file can be moved/emailed), plus the exact run parameters
# and the methodological summary. PDF rendering is best-effort (headless
# Chrome) and is disabled here for a fast, deterministic check.
qa_register_test(
  "CORE_008", "regression", "high",
  "Automatic analysis report embeds figures, parameters and methodology in a self-contained file",
  function() {
    old <- options(triple_a_render_pdf = FALSE)
    on.exit(options(old), add = TRUE)

    fig <- tempfile(fileext = ".png")
    ggplot2::ggsave(
      fig,
      ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
        ggplot2::geom_point(),
      width = 4, height = 3
    )

    methodology <- data.frame(
      Analysis = "Community structure", Output = "PERMANOVA table",
      Objective = "Test whether centroids differ", Method = "vegan::adonis2",
      Statistical_test = "PERMANOVA", Significance_criterion = "P <= 0.05",
      Interpretation = "A significant result supports differences in centroids",
      Run_parameters = "permutations=999", stringsAsFactors = FALSE
    )
    settings <- data.frame(
      Parameter = c("distance", "permutations"),
      Value = c("bray", "999"), stringsAsFactors = FALSE
    )
    functions <- data.frame(
      Function = "Arsenate_respiration", Category = "Arsenic metabolism",
      stringsAsFactors = FALSE
    )
    packages <- data.frame(Package = "vegan", Version = "2.6-4", stringsAsFactors = FALSE)

    report_dir <- tempfile("triplea_report_")
    dir.create(report_dir, recursive = TRUE)
    res <- aaa_write_analysis_report(
      report_dir = report_dir, run_id = "TEST_RUN", elapsed_seconds = 1.0,
      analyses = "community_structure", methodology_table = methodology,
      settings_table = settings, function_table = functions,
      package_versions = packages, figure_files = fig
    )

    qa_expect_true(file.exists(res$html), "The analysis report HTML was not generated.")
    qa_expect_true(is.na(res$pdf), "PDF must be skipped when triple_a_render_pdf is FALSE.")

    html <- paste(readLines(res$html, warn = FALSE), collapse = "\n")
    qa_expect_true(
      grepl("data:image/png;base64,", html, fixed = TRUE),
      "Report figures were not embedded as self-contained data URIs."
    )
    qa_expect_true(
      !grepl("<img src='[A-Za-z]:", html),
      "Report still references figures by absolute file path instead of embedding them."
    )
    qa_expect_true(
      grepl("bray", html, fixed = TRUE) && grepl("PERMANOVA", html, fixed = TRUE),
      "Report is missing the run-parameter table or the methodological summary."
    )
    TRUE
  }
)
