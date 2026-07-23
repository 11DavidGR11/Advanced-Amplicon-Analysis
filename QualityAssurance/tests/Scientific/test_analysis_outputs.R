qa_register_test("ANALYSIS_001", "regression", "high", "PLS-DA, RDA, PCA and differential analysis expose standard summaries and diagnostics", function() {
  sup <- paste(readLines(file.path(QA_ROOT,"AAApp","Common","Engine","Core","aaa_supervised_analysis.R"),warn=FALSE),collapse="\n")
  com <- paste(readLines(file.path(QA_ROOT,"AAApp","Common","Engine","Core","aaa_community_analysis.R"),warn=FALSE),collapse="\n")
  dif <- paste(readLines(file.path(QA_ROOT,"AAApp","Common","Engine","Core","aaa_extended_analyses.R"),warn=FALSE),collapse="\n")
  qa_expect_true(all(vapply(c("VIP_scores","Confusion_matrix","Balanced accuracy","Permutation p-value"), function(x) grepl(x,sup,fixed=TRUE), logical(1))),"PLS-DA diagnostics missing")
  qa_expect_true(all(vapply(c("Axis_tests","Adjusted R2","VIF"), function(x) grepl(x,sup,fixed=TRUE), logical(1))),"RDA diagnostics missing")
  qa_expect_true(all(vapply(c("PCA_loadings","PCA taxon contributions"), function(x) grepl(x,com,fixed=TRUE), logical(1))),"PCA diagnostics missing")
  qa_expect_true(grepl("QQ plot",dif,fixed=TRUE),"Differential QQ plot missing")
  TRUE
})

qa_register_test(
  "ANALYSIS_002", "regression", "critical",
  "Reported multivariate statistics are numerically valid, not just present",
  function() {
    # ANALYSIS_001 only greps the source for column names, so it stayed green
    # while R2Y/Q2 were always negative and the RDA constrained-variance metric
    # was always 100. This asserts the values themselves.
    abundance <- utils::read.csv(
      file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv"),
      check.names = FALSE
    )
    samples <- grep("^Sample_", names(abundance), value = TRUE)
    metadata <- utils::read.csv(
      file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "environmental_metadata.csv"),
      stringsAsFactors = FALSE
    )
    roles <- data.frame(
      Column = c("Sample", "Treatment", "Temperature", "pH", "Moisture", "NH4"),
      Role = c("identifier", "experimental_factor", rep("environmental_variable", 4)),
      stringsAsFactors = FALSE
    )
    make <- function(groups, per_group) {
      cols <- samples[seq_len(groups * per_group)]
      design <- data.frame(
        Sample_column = cols,
        Treatment = rep(LETTERS[seq_len(groups)], each = per_group),
        Replicate = rep(seq_len(per_group), times = groups),
        stringsAsFactors = FALSE
      )
      aaa_new_dataset(abundance, design, metadata = metadata, metadata_roles = roles)
    }
    metric <- function(table, name) table$Value[table$Metric == name]

    # --- PLS-DA: R2Y and Q2 must be consistent with the classification the same
    # model achieves. Leaving the dummy response matrix uncentred removed the
    # intercept, so the residuals absorbed the class means and both statistics
    # came out at or below zero even for a model classifying at 92%. Bounds
    # alone are too weak to catch that (R2Y can still land above zero), so the
    # assertion ties Q2 to the out-of-fold accuracy actually obtained.
    check_plsda <- function(result, classes, label) {
      summary_table <- result$tables$summary
      r2y <- metric(summary_table, "R2Y")
      q2 <- metric(summary_table, "Q2")
      accuracy <- metric(summary_table, "Cross-validated accuracy")
      chance <- 1 / classes

      qa_expect_true(
        is.finite(r2y) && r2y > 0 && r2y <= 1,
        paste0(label, ": R2Y is not a valid goodness-of-fit value: ", r2y)
      )
      qa_expect_true(
        is.finite(q2) && q2 <= 1,
        paste0(label, ": Q2 exceeds 1: ", q2)
      )
      qa_expect_true(
        r2y >= q2,
        paste0(label, ": cross-validated Q2 cannot exceed the fitted R2Y.")
      )
      # The fixture separates its groups cleanly, so this branch is the one
      # that actually runs; it is guarded so the test degrades to the bounds
      # checks instead of failing spuriously if the fixture ever changes.
      if (is.finite(accuracy) && accuracy > chance + 0.25) {
        qa_expect_true(
          q2 > 0,
          paste0(
            label, ": Q2 is ", round(q2, 4), " while the model classifies ",
            round(100 * accuracy, 1), "% of held-out samples correctly (chance = ",
            round(100 * chance, 1), "%). A predictive model cannot have a negative Q2."
          )
        )
      }
    }

    check_plsda(
      aaa_plsda_analysis(
        make(3, 4), "counts", n_components = 2, cv_folds = 3,
        permutation_repetitions = 9, project_dir = tempfile("triplea_plsda_num_")
      ),
      classes = 3, label = "PLS-DA (3 groups)"
    )

    # --- PLS-DA is not restricted to (classes - 1) components: with two groups
    # the second component used to collapse to a constant zero column.
    two_group <- aaa_plsda_analysis(
      make(2, 4), "counts", n_components = 2, cv_folds = 2,
      permutation_repetitions = 0, project_dir = tempfile("triplea_plsda_2g_")
    )
    check_plsda(two_group, classes = 2, label = "PLS-DA (2 groups)")
    qa_expect_true(
      metric(two_group$tables$summary, "Components") >= 2,
      "A two-group PLS-DA was capped at a single component."
    )
    qa_expect_true(
      length(unique(round(two_group$tables$scores$Component2, 9))) > 1L,
      "The second PLS-DA component is constant, so the score plot is degenerate."
    )

    # --- RDA: constrained variance must be the share of TOTAL inertia, which
    # equals R2. Summing the per-axis percentages always gave exactly 100.
    rda <- aaa_rda_analysis(
      make(3, 4), "counts", c("Temperature", "pH"),
      permutations = 99, project_dir = tempfile("triplea_rda_num_")
    )
    constrained <- metric(rda$tables$summary, "Constrained variance (% of total inertia)")
    r2 <- metric(rda$tables$summary, "R2")
    qa_expect_true(
      is.finite(constrained) && constrained > 0 && constrained < 100,
      paste0("RDA constrained variance is not a share of total inertia: ", constrained)
    )
    qa_expect_true(
      abs(constrained - 100 * r2) < 1e-6,
      "RDA constrained variance does not agree with the model R2."
    )
    qa_expect_true(
      abs(sum(rda$tables$explained_variance$Percent) - 100) < 1e-6,
      "Per-axis RDA percentages should still be relative to the constrained part."
    )
    TRUE
  }
)

qa_register_test(
  "ANALYSIS_003", "regression", "critical",
  "Taxon-association engines receive real sequencing depth for count input",
  function() {
    # prepared$wide is rescaled to percentages, so rounding it handed ANCOM-BC2
    # a table whose sample totals were ~100 instead of the real library size.
    abundance <- utils::read.csv(
      file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv"),
      check.names = FALSE
    )
    samples <- grep("^Sample_", names(abundance), value = TRUE)
    # Scale to a realistic sequencing depth so a percentage table is clearly
    # distinguishable from a count table.
    abundance[samples] <- round(sweep(
      as.matrix(abundance[samples]), 2,
      colSums(abundance[samples]), "/"
    ) * 50000)
    metadata <- utils::read.csv(
      file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "environmental_metadata.csv"),
      stringsAsFactors = FALSE
    )
    design <- data.frame(
      Sample_column = samples,
      Treatment = rep(c("A", "B", "C"), each = 4),
      Replicate = rep(1:4, 3),
      stringsAsFactors = FALSE
    )
    roles <- data.frame(
      Column = c("Sample", "Treatment", "Temperature", "pH", "Moisture", "NH4"),
      Role = c("identifier", "experimental_factor", rep("environmental_variable", 4)),
      stringsAsFactors = FALSE
    )
    dataset <- aaa_new_dataset(abundance, design, metadata = metadata, metadata_roles = roles)

    expected <- colSums(abundance[samples])
    prepared <- aaa_taxon_association_input(
      dataset, "counts", c("Treatment", "Temperature"),
      tempfile("triplea_taxassoc_"), "counts_depth"
    )
    observed <- rowSums(prepared$counts)

    qa_expect_true(
      all(abs(observed[names(expected)] - expected) <= 0.01 * expected),
      paste0(
        "Count input did not reach the taxon-association engines at its real depth: ",
        "expected ~", round(mean(expected)), " reads per sample, got ~", round(mean(observed)), "."
      )
    )
    qa_expect_true(
      max(prepared$counts) > 100,
      "The reconstructed count table looks like rounded percentages."
    )
    TRUE
  }
)
