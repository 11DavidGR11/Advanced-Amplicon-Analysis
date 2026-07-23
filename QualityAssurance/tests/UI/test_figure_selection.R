qa_register_test(
  "UI_007", "regression", "high",
  "Each analysis exposes selectable figures while tables remain automatic",
  function() {
    app <- qa_read_app_source(QA_ROOT)
    required <- c("figures_functional_potential", "figures_functional_abundance", "figures_top_abundance",
      "figures_differential_abundance", "figures_community_structure", "figures_community_comparison",
      "figures_plsda", "figures_rda")
    qa_expect_true(all(vapply(required, grepl, logical(1), x = app, fixed = TRUE)),
      "One or more per-analysis figure selectors are missing")
  }
)
