# Runtime coverage for the envfit / partial RDA / dbRDA / partial dbRDA /
# variance partitioning / ANCOM-BC2 / MaAsLin2 engines in
# aaa_advanced_environmental_taxon.R. These analyses previously had no
# automated execution coverage (only static UI/wiring checks elsewhere),
# which is how real runtime bugs in the ANCOM-BC2 and MaAsLin2 engines went
# undetected until they were exercised by hand against installed packages.

qa_advanced_environmental_taxon_dataset <- function() {
  abundance <- utils::read.csv(
    file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv"),
    check.names = FALSE
  )
  samples <- grep("^Sample_", names(abundance), value = TRUE)
  design <- data.frame(
    Sample_column = samples,
    Treatment = rep(c("A", "B", "C"), each = 4),
    Replicate = rep(1:4, 3),
    stringsAsFactors = FALSE
  )
  metadata <- utils::read.csv(
    file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "environmental_metadata.csv"),
    stringsAsFactors = FALSE
  )
  metadata_roles <- data.frame(
    Column = c("Sample", "Treatment", "Temperature", "pH", "Moisture", "NH4"),
    Role = c("identifier", "experimental_factor", "environmental_variable",
             "environmental_variable", "environmental_variable", "environmental_variable"),
    stringsAsFactors = FALSE
  )
  aaa_new_dataset(abundance, design, metadata = metadata, metadata_roles = metadata_roles)
}

qa_register_test(
  "SCI_016", "regression", "critical",
  "envfit fits environmental vectors onto an unconstrained ordination",
  function() {
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_envfit_analysis(
      dataset, "counts", c("Temperature", "pH"),
      transformation = "hellinger", permutations = 49,
      project_dir = tempfile("triplea_envfit_")
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "envfit did not return a Triple_A_result.")
    qa_expect_true(nrow(result$tables$Vectors) == 2L, "envfit did not fit both environmental vectors.")
    TRUE
  }
)

qa_register_test(
  "SCI_016B", "regression", "critical",
  "envfit reports a real P-value for a categorical (factor) environmental variable",
  function() {
    dataset <- qa_advanced_environmental_taxon_dataset()
    # vegan::envfit() names factors$centroids rows "<Variable><Level>" (no
    # separator) and pvals by variable only; a previous version of this
    # engine stripped everything up to "=" before matching against pvals,
    # which is a no-op against that format and always produced NA here.
    result <- aaa_envfit_analysis(
      dataset, "counts", "Treatment",
      transformation = "hellinger", permutations = 49,
      project_dir = tempfile("triplea_envfit_factor_")
    )
    centroids <- result$tables$Factor_centroids
    qa_expect_true(nrow(centroids) == 3L, "envfit did not report all three Treatment level centroids.")
    qa_expect_true(
      all(is.finite(centroids$P_value)),
      "envfit factor P-values were NA instead of the real permutation P-value."
    )
    qa_expect_true(
      length(unique(centroids$P_value)) == 1L,
      "All levels of the same factor should share one P-value, but they differ."
    )
    TRUE
  }
)

qa_register_test(
  "SCI_017", "regression", "critical",
  "Partial RDA conditions ordination on a declared experimental factor",
  function() {
    dataset <- qa_advanced_environmental_taxon_dataset()
    # Moisture varies within each Treatment group but not between groups, so
    # it stays informative once Treatment is conditioned out; pH is a
    # deterministic function of Treatment in this fixture and would be
    # perfectly aliased against the conditioning variable.
    result <- aaa_constrained_analysis(
      dataset, "counts", environmental_variables = "Moisture",
      conditioning_variables = "Treatment", transformation = "hellinger",
      distance = NULL, permutations = 49,
      project_dir = tempfile("triplea_prda_"), analysis_name = "Partial_RDA"
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "Partial RDA did not return a Triple_A_result.")
    qa_expect_true(
      grepl("Partial", result$tables$Model_variables$Conditioning[1]) ||
        nzchar(result$tables$Model_variables$Conditioning[1]),
      "Partial RDA did not record the conditioning variable."
    )
    TRUE
  }
)

qa_register_test(
  "SCI_018", "regression", "critical",
  "dbRDA exposes a Bray-Curtis constrained ordination",
  function() {
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_constrained_analysis(
      dataset, "counts", environmental_variables = "pH",
      conditioning_variables = character(), transformation = "hellinger",
      distance = "bray", permutations = 49,
      project_dir = tempfile("triplea_dbrda_"), analysis_name = "dbRDA"
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "dbRDA did not return a Triple_A_result.")
    qa_expect_true(nrow(result$tables$Variance) > 0L, "dbRDA did not report ordination axis variance.")
    TRUE
  }
)

qa_register_test(
  "SCI_019", "regression", "critical",
  "Partial dbRDA combines a distance measure with a conditioning factor",
  function() {
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_constrained_analysis(
      dataset, "counts", environmental_variables = "Moisture",
      conditioning_variables = "Treatment", transformation = "hellinger",
      distance = "bray", permutations = 49,
      project_dir = tempfile("triplea_pdbrda_"), analysis_name = "Partial_dbRDA"
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "Partial dbRDA did not return a Triple_A_result.")
    TRUE
  }
)

qa_register_test(
  "SCI_020", "regression", "critical",
  "Variance partitioning requires both an environmental and an experimental predictor set",
  function() {
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_variance_partitioning(
      dataset, "counts", environmental_variables = "Moisture",
      experimental_factors = "Treatment", transformation = "hellinger",
      project_dir = tempfile("triplea_varpart_")
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "Variance partitioning did not return a Triple_A_result.")
    qa_expect_true(nrow(result$tables$Fractions) > 0L, "Variance partitioning did not report variance fractions.")
    TRUE
  }
)

qa_register_test(
  "SCI_021", "regression", "critical",
  "ANCOM-BC2 runs end-to-end on a canonical count dataset",
  function() {
    if (!requireNamespace("ANCOMBC", quietly = TRUE) ||
        !requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
      return(TRUE)
    }
    dataset <- qa_advanced_environmental_taxon_dataset()
    # ANCOMBC2 always warns about unstable variance estimation below 5
    # replicates per group; this fixture intentionally stays small for
    # speed, so the warning is expected and not a signal of a defect here.
    result <- suppressWarnings(aaa_ancombc2_analysis(
      dataset, "counts", variables = "Treatment", group_variable = "Treatment",
      project_dir = tempfile("triplea_ancombc_"), analysis_name = "ANCOM_BC2"
    ))
    qa_expect_true(inherits(result, "Triple_A_result"), "ANCOM-BC2 did not return a Triple_A_result.")
    qa_expect_true(nrow(result$tables$Results) > 0L, "ANCOM-BC2 did not return per-taxon results.")
    qa_expect_true(
      !is.null(result$files[["volcano"]]) && file.exists(result$files[["volcano"]]),
      "ANCOM-BC2 did not save its volcano plot."
    )
    TRUE
  }
)

qa_register_test(
  "SCI_022", "regression", "critical",
  "MaAsLin2 runs end-to-end and picks a reference level for a 3+ level fixed effect",
  function() {
    if (!requireNamespace("Maaslin2", quietly = TRUE)) {
      return(TRUE)
    }
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_maaslin_analysis(
      dataset, "counts", variables = "Treatment",
      project_dir = tempfile("triplea_maaslin_"), analysis_name = "MaAsLin2"
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "MaAsLin2 did not return a Triple_A_result.")
    qa_expect_true(is.data.frame(result$tables$results), "MaAsLin2 did not return a results table.")
    qa_expect_true(
      !is.null(result$files[["forest"]]) && file.exists(result$files[["forest"]]),
      "MaAsLin2 did not save its coefficient forest plot."
    )
    TRUE
  }
)

qa_register_test(
  "SCI_026", "regression", "high",
  "Community structure saves the pairwise PERMANOVA, beta-dispersion, ANOSIM and PERMANOVA-variance figures",
  function() {
    # The graphical summaries draw only from the pairwise
    # PERMANOVA / ANOSIM / beta-dispersion tables the engine already computes.
    # No test previously ran aaa_community_analysis() end-to-end, so this
    # locks in that the four new figures are produced for a 3-group design.
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_community_analysis(
      dataset, abundance_type = "counts", transformation = "hellinger",
      distance_method = "bray", permutations = 99, nmds_trymax = 5,
      project_dir = tempfile("triplea_community_")
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "Community structure did not return a Triple_A_result.")
    figures <- c("pairwise_permanova", "beta_dispersion", "anosim", "permanova_variance")
    for (key in figures) {
      qa_expect_true(
        !is.null(result$files[[key]]) && file.exists(result$files[[key]]),
        paste0("Community structure did not save the '", key, "' figure.")
      )
      qa_expect_true(
        !is.null(result$plots[[key]]),
        paste0("Community structure did not expose the '", key, "' plot object.")
      )
    }
    TRUE
  }
)

qa_register_test(
  "SCI_034", "regression", "high",
  "Hierarchical clustering produces a dendrogram and a cluster table from the ordination distances",
  function() {
    # "Hierarchical clustering" was offered as a selectable analysis in the
    # interface while no implementation existed: the checkbox produced nothing
    # and no output was catalogued. This locks in that it now runs, that it
    # reuses the same distance matrix as the ordinations rather than computing
    # its own, and that both of its outputs reach disk.
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_community_analysis(
      dataset, abundance_type = "counts", transformation = "hellinger",
      distance_method = "bray", permutations = 99, nmds_trymax = 5,
      project_dir = tempfile("triplea_clustering_")
    )

    qa_expect_true(
      !is.null(result$files[["dendrogram"]]) && file.exists(result$files[["dendrogram"]]),
      "Hierarchical clustering did not save the dendrogram figure."
    )
    qa_expect_true(
      !is.null(result$plots[["dendrogram"]]),
      "Hierarchical clustering did not expose the dendrogram plot object."
    )

    clusters <- result$tables$hierarchical_clustering
    qa_expect_true(
      is.data.frame(clusters) && nrow(clusters) > 0L,
      "Hierarchical clustering did not return a cluster table."
    )
    required_columns <- c("Sample_column", "Treatment", "Cluster", "Dendrogram_order")
    qa_expect_true(
      all(required_columns %in% names(clusters)),
      paste(
        "Cluster table is missing columns:",
        paste(setdiff(required_columns, names(clusters)), collapse = ", ")
      )
    )

    # Every sample must appear exactly once, or the dendrogram is not a
    # partition of the dataset and the cluster column cannot be trusted.
    design_samples <- as.character(dataset$sample_design$Sample_column)
    qa_expect_true(
      setequal(clusters$Sample_column, design_samples) &&
        !anyDuplicated(clusters$Sample_column),
      "The cluster table does not contain each sample exactly once."
    )
    qa_expect_true(
      all(is.finite(clusters$Cluster)) && all(clusters$Cluster >= 1L),
      "Cluster membership contains non-finite or non-positive values."
    )

    # The cut is made at the number of declared groups, so the number of
    # clusters can never exceed it.
    declared <- length(unique(as.character(dataset$sample_design$Treatment)))
    qa_expect_true(
      length(unique(clusters$Cluster)) <= max(2L, declared),
      "Hierarchical clustering produced more clusters than treatment groups."
    )

    qa_expect_true(
      identical(result$metadata$clustering_linkage, "average"),
      "Hierarchical clustering did not use average linkage."
    )

    # Both outputs must be catalogued and documented, or they cannot be
    # selected in the interface nor described in the methods workbook.
    catalogue <- aaa_output_catalogue()
    for (id in c("dendrogram_plot", "hierarchical_clustering_table")) {
      qa_expect_true(
        id %in% catalogue$ID,
        paste0("Output '", id, "' is not in the output catalogue.")
      )
    }
    TRUE
  }
)

qa_register_test(
  "SCI_023", "regression", "critical",
  "sPLS-DA runs end-to-end and predicts fitted classes",
  function() {
    # A prior version called mixOmics::predict(...); predict() is a base R
    # generic that mixOmics only provides an S3 method for, so the
    # namespace-qualified call failed at runtime with "'predict' is not an
    # exported object from 'namespace:mixOmics'" the first time a user
    # actually ran sPLS-DA. No test executed the full function, only its
    # matrix-preprocessing helpers, so this went unnoticed until reported.
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_splsda_analysis(
      dataset, "counts", n_components = 2, tune = FALSE,
      project_dir = tempfile("triplea_splsda_")
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "sPLS-DA did not return a Triple_A_result.")
    qa_expect_true(
      is.data.frame(result$tables$performance) && "Predicted" %in% names(result$tables$performance),
      "sPLS-DA did not return predicted classes."
    )
    qa_expect_true(nrow(result$tables$confusion_matrix) > 0L, "sPLS-DA did not return a confusion matrix.")
    TRUE
  }
)

qa_register_test(
  "SCI_024", "regression", "critical",
  "PLS-DA runs end-to-end and predicts fitted classes",
  function() {
    # ANALYSIS_001 only greps the source file for diagnostic column names; it
    # never actually calls aaa_plsda_analysis(), which is exactly the kind of
    # gap that let the sibling sPLS-DA function ship with a runtime-only bug
    # (see SCI_023). This exercises the real function end-to-end instead.
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_plsda_analysis(
      dataset, "counts", n_components = 2,
      project_dir = tempfile("triplea_plsda_")
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "PLS-DA did not return a Triple_A_result.")
    qa_expect_true(
      is.data.frame(result$tables$performance) && "Predicted" %in% names(result$tables$performance),
      "PLS-DA did not return predicted classes."
    )
    qa_expect_true(is.data.frame(result$tables$vip_scores), "PLS-DA did not return VIP scores.")
    TRUE
  }
)

qa_register_test(
  "SCI_025", "regression", "critical",
  "Top abundance runs end-to-end and saves its plots",
  function() {
    # aaa_top_abundance() stores its Taxon column as a factor (for plot
    # ordering) and later measures those same labels via
    # aaa_flipped_axis_plot_width(). aaa_measure_text_width_inches() called
    # nzchar() directly on whatever it was given, which crashed with
    # "'nzchar()' requires a character vector" (dispatch through
    # [.factor -> NextMethod("[")) the first time a real user ran this
    # analysis. No test executed the full function end-to-end, only the
    # underlying width helper with plain character input (see CORE_013).
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_top_abundance(
      dataset, top_n = 20, abundance_type = "counts",
      project_dir = tempfile("triplea_top_abundance_")
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "Top abundance did not return a Triple_A_result.")
    qa_expect_true(
      !is.null(result$files[["lollipop"]]) && file.exists(result$files[["lollipop"]]),
      "Top abundance did not save its lollipop plot (this is the plot whose width computation crashed)."
    )
    TRUE
  }
)

qa_register_test(
  "SCI_027", "regression", "critical",
  "Differential abundance runs end-to-end and saves its volcano/MA plots",
  function() {
    # No test previously executed aaa_differential_abundance() end-to-end; that
    # gap is exactly how the ggplot2 >= 4.0 ggsave create.dir regression reached
    # the volcano/MA save path unnoticed until a full workflow was run by hand.
    dataset <- qa_advanced_environmental_taxon_dataset()
    result <- aaa_differential_abundance(
      dataset, abundance_type = "counts", method = "wilcox",
      min_prevalence = 0.1, min_mean_abundance = 0,
      project_dir = tempfile("triplea_diff_")
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "Differential abundance did not return a Triple_A_result.")
    volcano_files <- grep("^volcano_", names(result$files), value = TRUE)
    ma_files <- grep("^ma_", names(result$files), value = TRUE)
    qa_expect_true(length(volcano_files) >= 1L, "Differential abundance produced no volcano plot.")
    qa_expect_true(
      all(file.exists(unlist(result$files[volcano_files]))) &&
        all(file.exists(unlist(result$files[ma_files]))),
      "Differential abundance did not save its volcano/MA plot files."
    )
    qa_expect_true(is.data.frame(result$tables$combined), "Differential abundance did not return a combined results table.")
    TRUE
  }
)

qa_register_test(
  "SCI_028", "regression", "critical",
  "Potential metabolomic pathways abundance runs end-to-end and saves its figures",
  function() {
    # Exercises aaa_functional_abundance() end-to-end with a synthetic pathway
    # built from the prepared abundance table, covering the heatmap, bar and dot
    # plots that previously had no runtime test.
    dataset <- qa_advanced_environmental_taxon_dataset()
    prepared <- aaa_prepare_amplicon_data(dataset, "counts", tempfile("triplea_prep_"), "probe")
    taxa <- unique(prepared$wide[, c("Taxonomy", "Genus", "Tax_level")])
    taxa$Potential <- "Test pathway"
    pathways <- list(`Test pathway` = list(results = taxa, include = "Test pathway"))
    result <- aaa_functional_abundance(
      dataset, pathways, abundance_type = "counts",
      project_dir = tempfile("triplea_func_")
    )
    qa_expect_true(inherits(result, "Triple_A_result"), "Functional abundance did not return a Triple_A_result.")
    for (key in c("heatmap", "barplot", "dotplot")) {
      qa_expect_true(
        !is.null(result$files[[key]]) && file.exists(result$files[[key]]),
        paste0("Functional abundance did not save its '", key, "' figure.")
      )
    }
    TRUE
  }
)
