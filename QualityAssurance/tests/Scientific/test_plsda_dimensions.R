qa_register_test(
  "PLSDA_002", "regression", "critical",
  "PLS-DA response arrays are normalized to sample-by-class matrices",
  function() {
    qa_expect_true(exists("aaa_plsda_component_matrix", mode = "function"),
                   "aaa_plsda_component_matrix() is unavailable")
    arr <- array(seq_len(12L * 3L * 2L), dim = c(12L, 3L, 2L))
    out <- aaa_plsda_component_matrix(arr, 2L, 12L, 3L, "test predictions")
    qa_expect_true(identical(dim(out), c(12L, 3L)), "Array normalization returned wrong dimensions")
    mat <- matrix(seq_len(36L), nrow = 12L, ncol = 3L)
    out2 <- aaa_plsda_component_matrix(mat, 1L, 12L, 3L, "test fitted values")
    qa_expect_true(identical(dim(out2), c(12L, 3L)), "Matrix normalization returned wrong dimensions")
    TRUE
  }
)

qa_register_test(
  "PLUGIN_004", "regression", "high",
  "All plugin runners and validators accept the canonical extensible signature",
  function() {
    ids <- qa_plugin_ids()
    problems <- character()
    for (id in ids) {
      p <- qa_get_plugin_definition(id)
      runner_formals <- names(formals(p$runner))
      validator_formals <- names(formals(p$validator))
      if (!all(c("context", "parameters", "...") %in% runner_formals)) {
        problems <- c(problems, paste0(id, " runner signature"))
      }
      if (!all(c("context", "...") %in% validator_formals)) {
        problems <- c(problems, paste0(id, " validator signature"))
      }
    }
    qa_expect_true(!length(problems), paste(problems, collapse = "; "))
    TRUE
  }
)
