qa_register_test(
  "PLSDA_001",
  "regression",
  "high",
  "PLS-DA preprocessing removes non-finite and near-constant predictors",
  function() {
    source(
      file.path(QA_ROOT, "AAApp", "Common", "Engine", "Core", "aaa_supervised_analysis.R"),
      local = FALSE
    )

    x <- cbind(
      informative_1 = c(1, 2, 3, 4, 5, 6),
      informative_2 = c(6, 4, 5, 2, 3, 1),
      contains_inf = c(1, 2, Inf, 4, 5, 6),
      constant = rep(1, 6),
      near_constant = 1 + seq_len(6) * 1e-18
    )

    prepared <- aaa_prepare_plsda_matrix(x)

    qa_expect_true(
      is.matrix(prepared) && all(is.finite(prepared)),
      "Prepared PLS-DA matrix contains non-finite values."
    )
    qa_expect_true(
      ncol(prepared) >= 2L,
      "Informative predictors were incorrectly removed."
    )
    qa_expect_true(
      !"constant" %in% colnames(prepared) &&
        !"near_constant" %in% colnames(prepared),
      "Constant or numerically unstable predictors were retained."
    )
  }
)
