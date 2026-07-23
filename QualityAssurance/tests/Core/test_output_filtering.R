# aaa_filter_selected_outputs() had no branch for results$splsda at all, so
# deselecting the sPLS-DA score/confusion/selected-feature plots in the
# workflow config had no effect: those plots and files were always kept.
# No test exercised this function directly, which is how the gap went
# unnoticed alongside the working PLS-DA/RDA branches beside it.

qa_fake_analysis_result <- function(file_keys, plot_keys) {
  files <- stats::setNames(
    vapply(file_keys, function(k) tempfile(fileext = ".png"), character(1)),
    file_keys
  )
  for (path in files) file.create(path)
  plots <- stats::setNames(as.list(seq_along(plot_keys)), plot_keys)
  list(files = files, plots = plots)
}

qa_register_test(
  "CORE_005", "regression", "critical",
  "Deselecting a PLS-DA output removes its plot and file (established working pattern)",
  function() {
    results <- list(
      plsda = qa_fake_analysis_result(
        c("plot", "confusion", "vip", "summary"),
        c("plsda", "confusion_matrix", "vip")
      )
    )
    filtered <- aaa_filter_selected_outputs(results, character())
    qa_expect_true(
      is.null(filtered$plsda$plots$plsda) && !"plot" %in% names(filtered$plsda$files),
      "Deselecting plsda_plot did not remove the PLS-DA score plot/file."
    )
    TRUE
  }
)

qa_register_test(
  "CORE_006", "regression", "critical",
  "Deselecting an sPLS-DA output removes its plot and file",
  function() {
    results <- list(
      splsda = qa_fake_analysis_result(
        c("plot", "confusion", "selected", "summary"),
        c("splsda", "confusion_matrix", "selected_features")
      )
    )
    all_selected <- c("splsda_plot", "splsda_confusion_plot", "splsda_selected_features")

    kept <- aaa_filter_selected_outputs(results, all_selected)
    qa_expect_true(
      identical(names(kept$splsda$plots), c("splsda", "confusion_matrix", "selected_features")),
      "Selected sPLS-DA outputs were removed even though they were requested."
    )

    filtered <- aaa_filter_selected_outputs(results, setdiff(all_selected, "splsda_plot"))
    qa_expect_true(
      is.null(filtered$splsda$plots$splsda) && !"plot" %in% names(filtered$splsda$files),
      "Deselecting splsda_plot did not remove the sPLS-DA score plot/file."
    )

    filtered2 <- aaa_filter_selected_outputs(results, setdiff(all_selected, "splsda_confusion_plot"))
    qa_expect_true(
      is.null(filtered2$splsda$plots$confusion_matrix) && !"confusion" %in% names(filtered2$splsda$files),
      "Deselecting splsda_confusion_plot did not remove the sPLS-DA confusion plot/file."
    )

    filtered3 <- aaa_filter_selected_outputs(results, setdiff(all_selected, "splsda_selected_features"))
    qa_expect_true(
      is.null(filtered3$splsda$plots$selected_features) && !"selected" %in% names(filtered3$splsda$files),
      "Deselecting splsda_selected_features did not remove the sPLS-DA feature plot/file."
    )
    TRUE
  }
)
