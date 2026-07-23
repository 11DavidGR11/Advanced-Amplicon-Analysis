# =============================================================================
# Triple_A centralized settings
# =============================================================================

TRIPLE_A_NAME <- "Advanced_Amplicon_Analysis"
TRIPLE_A_SHORT_NAME <- "Triple_A"

aaa_get_defaults <- function() {
  list(
    output_dir = "Results",
    progress_verbosity = "standard",
    analyses = c(
      "functional_potential",
      "top_abundance",
      "differential_abundance",
      "functional_abundance",
      "community_structure",
      "plsda",
      "splsda",
      "rda", "envfit", "partial_rda", "dbrda", "partial_dbrda", "variance_partitioning", "ancombc2", "maaslin", "differential_functions", "functional_enrichment"
    ),
    outputs = c(
      "functional_heatmap",
      "top_heatmap",
      "top_lollipop",
      "top_composition",
      "top_distribution",
      "volcano_plot",
      "ma_plot",
      "functional_abundance_heatmap",
      "functional_abundance_barplot",
      "functional_abundance_dotplot",
      "functional_contributors",
      "pca_plot",
      "pca_scree",
      "pcoa_plot",
      "nmds_plot",
      "permanova_table",
      "alpha_diversity_plot",
      "alpha_diversity_table",
      "beta_distance_table",
      "pairwise_permanova_heatmap",
      "beta_dispersion_boxplot",
      "anosim_plot",
      "permanova_variance_plot",
      "dendrogram_plot",
      "hierarchical_clustering_table",
      "plsda_plot",
      "plsda_performance",
      "plsda_loadings",
      "splsda_plot",
      "splsda_confusion_plot",
      "splsda_selected_features",
      "splsda_performance",
      "splsda_stability",
      "rda_plot",
      "rda_anova",
      "rda_variance",
      "differential_functions_plot",
      "differential_functions_table",
      "functional_enrichment_plot",
      "functional_enrichment_table"
    ),
    top_abundance = list(
      top_n = 20
    ),
    differential_abundance = list(
      method = "wilcox",
      paired = FALSE,
      pseudocount = 1e-6,
      min_prevalence = 0.20,
      min_mean_abundance = 0.01,
      alpha = 0.05,
      log2fc_threshold = 1,
      top_n_table = 25,
      max_labels = 10,
      label_only_significant = FALSE,
      colour_by = "log2FC",
      point_size = "mean_abundance",
      x_limit = 6,
      comparisons = NULL,
      dynamic_ylim = TRUE
    ),
    functional_abundance = list(
      top_taxa_per_pathway = 5
    ),
    community_structure = list(
      transformation = "hellinger",
      distance_method = "bray",
      permutations = 999,
      significance_alpha = 0.05,
      nmds_trymax = 50,
      show_sample_labels = FALSE
    ),
    supervised_multivariate = list(
      transformation = "hellinger",
      plsda_components = 2,
      plsda_cv_folds = 5,
      plsda_seed = 123,
      splsda_components = 2,
      splsda_cv_folds = 5,
      splsda_repeats = 10,
      splsda_keepx = c(5, 10, 20),
      splsda_tune = TRUE,
      splsda_seed = 123,
      rda_permutations = 999,
      rda_alpha = 0.05,
      show_sample_labels = FALSE
    ),
    graphics = list(
      width = 8,
      height = 6,
      dpi = 300,
      legend_position = "bottom"
    )
  )
}

aaa_analysis_catalogue <- function() {
  data.frame(
    ID = c(
      "functional_potential",
      "top_abundance",
      "differential_abundance",
      "functional_abundance",
      "community_structure",
      "plsda",
      "splsda",
      "rda", "envfit", "partial_rda", "dbrda", "partial_dbrda", "variance_partitioning", "ancombc2", "maaslin", "differential_functions", "functional_enrichment"
    ),
    Name = c(
      "Functional potential",
      "Top abundance",
      "Differential abundance",
      "Potential metabolomic pathways abundance",
      "Community structure and diversity",
      "PLS-DA",
      "sPLS-DA",
      "RDA with environmental variables", "envfit", "Partial RDA", "dbRDA", "Partial dbRDA", "Variance partitioning", "ANCOM-BC2", "MaAsLin2", "Differential functions", "Functional enrichment"
    ),
    Description = c(
      "Reference-genome gene search and pathway classification.",
      "Heatmap and abundance summaries for the most abundant taxa.",
      "Pairwise differential testing with optional volcano and MA plots.",
      "Abundance summaries for selected putative metabolic pathways.",
      "PCA, PCoA, NMDS, PERMANOVA, alpha diversity and beta distances.",
      "Supervised discrimination of treatment groups with cross-validated PLS-DA.",
      "Sparse supervised discrimination with cross-validated taxon selection.",
      "Constrained ordination relating community composition to selected environmental variables.",
      "Permutation-based environmental vector fitting.", "Constrained ordination controlling for covariates.", "Distance-based constrained ordination.", "Distance-based constrained ordination controlling for covariates.", "Partitioning of explained community variation between predictor sets.", "Bias-corrected compositional differential abundance.", "Multivariable taxon-metadata association modelling.", "Comparison of inferred functional profiles between experimental groups.", "Over-representation of curated functions among differentially abundant taxa."
    ),
    stringsAsFactors = FALSE
  )
}

aaa_output_catalogue <- function() {
  data.frame(
    ID = c(
      "functional_heatmap",
      "top_heatmap",
      "top_lollipop",
      "top_composition",
      "top_distribution",
      "volcano_plot",
      "ma_plot",
      "qq_plot",
      "functional_abundance_heatmap",
      "functional_abundance_barplot",
      "functional_abundance_dotplot",
      "functional_contributors",
      "pca_plot",
      "pca_scree",
      "pcoa_plot",
      "nmds_plot",
      "permanova_table",
      "alpha_diversity_plot",
      "alpha_diversity_table",
      "beta_distance_table",
      "pairwise_permanova_heatmap",
      "beta_dispersion_boxplot",
      "anosim_plot",
      "permanova_variance_plot",
      "dendrogram_plot",
      "hierarchical_clustering_table",
      "plsda_plot",
      "plsda_confusion_plot",
      "plsda_vip_plot",
      "plsda_performance",
      "plsda_loadings",
      "splsda_plot",
      "splsda_confusion_plot",
      "splsda_selected_features",
      "splsda_performance",
      "splsda_stability",
      "rda_plot",
      "rda_variance_plot",
      "rda_anova",
      "rda_variance",
      "differential_functions_plot",
      "differential_functions_table",
      "functional_enrichment_plot",
      "functional_enrichment_table"
    ),
    Module = c(
      "functional_potential",
      rep("top_abundance", 4),
      rep("differential_abundance", 3),
      rep("functional_abundance", 4),
      rep("community_structure", 14),
      rep("plsda", 5),
      rep("splsda", 5),
      rep("rda", 4),
      rep("differential_functions", 2),
      rep("functional_enrichment", 2)
    ),
    Name = c(
      "Functional-potential heatmap",
      "Top-taxa heatmap",
      "Top-taxa lollipop plot",
      "Top-taxa composition plot",
      "Top-taxa distribution plot",
      "Volcano plot",
      "MA plot",
      "QQ plot",
      "Potential metabolomic pathways abundance heatmap",
      "Potential metabolomic pathways abundance bar plot",
      "Potential metabolomic pathways abundance dot plot",
      "Potential pathway contributors",
      "PCA sample plot",
      "PCA scree plot",
      "PCoA sample plot",
      "NMDS sample plot",
      "PERMANOVA table",
      "Alpha-diversity plot",
      "Alpha-diversity table",
      "Beta-distance matrix",
      "Pairwise PERMANOVA heatmap",
      "Beta-dispersion boxplot",
      "ANOSIM R gauge",
      "PERMANOVA explained-variance plot",
      "Hierarchical clustering dendrogram",
      "Hierarchical clustering table",
      "PLS-DA score plot",
      "PLS-DA confusion-matrix plot",
      "PLS-DA VIP plot",
      "PLS-DA cross-validation performance",
      "PLS-DA taxon loadings",
      "sPLS-DA score plot",
      "sPLS-DA confusion-matrix plot",
      "sPLS-DA selected-feature plot",
      "sPLS-DA cross-validation performance",
      "sPLS-DA feature-selection stability",
      "RDA ordination plot",
      "RDA explained-variance plot",
      "RDA permutation tests",
      "RDA explained variance",
      "Differential functions plot",
      "Differential functions table",
      "Functional enrichment plot",
      "Functional enrichment table"
    ),
    Type = c(
      rep("figure", 12),
      rep("figure", 4),
      "table", "figure", "table", "table",
      rep("figure", 5), "table",
      rep("figure", 3), "table", "table",
      rep("figure", 3), "table", "table",
      rep("figure", 2), "table", "table",
      "figure", "table", "figure", "table"
    ),
    stringsAsFactors = FALSE
  )
}


aaa_validate_selections <- function(analyses, outputs) {
  valid_analyses <- aaa_analysis_catalogue()$ID
  valid_outputs <- aaa_output_catalogue()$ID

  unknown_analyses <- setdiff(analyses, valid_analyses)
  unknown_outputs <- setdiff(outputs, valid_outputs)

  if (length(unknown_analyses) > 0) {
    stop(
      "Unknown analyses: ",
      paste(unknown_analyses, collapse = ", ")
    )
  }

  if (length(unknown_outputs) > 0) {
    stop(
      "Unknown outputs: ",
      paste(unknown_outputs, collapse = ", ")
    )
  }

  output_catalogue <- aaa_output_catalogue()
  selected_modules <- output_catalogue$Module[
    match(outputs, output_catalogue$ID)
  ]
  incompatible <- outputs[
    !selected_modules %in% analyses
  ]

  if (length(incompatible) > 0) {
    stop(
      "Outputs selected without their parent analysis: ",
      paste(incompatible, collapse = ", ")
    )
  }

  invisible(TRUE)
}
