# =============================================================================
# Triple_A analysis-method registry
#
# This file is the single source of truth for the objective, data preparation,
# statistical method, significance criteria and interpretation of each output.
# Shiny, metadata workbooks and documentation can read the same definitions.
# =============================================================================

aaa_analysis_method_registry_base <- function() {
  list(
    functional_potential = list(
      analysis_name = "Functional potential",
      outputs = list(
        functional_heatmap = list(
          output_name = "Functional-potential heatmap",
          example_image = "functional_heatmap.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = paste(
            "Summarise the abundance of taxa assigned to curated",
            "putative biological functions."
          ),
          preparation = paste(
            "Taxonomy is linked to a representative NCBI RefSeq genome;",
            "curated genes are searched in its GFF annotation."
          ),
          method = paste(
            "Rule-based classification using diagnostic, supporting and",
            "accessory gene evidence from the local biological registry."
          ),
          statistical_test = "No inferential statistical test.",
          significance = paste(
            "Not applicable. Categories represent curated genomic evidence,",
            "not statistical significance."
          ),
          interpretation = paste(
            "Positive categories indicate putative potential in a",
            "representative reference genome, not confirmed activity in the",
            "sampled strain. The Taxon_results table additionally reports the",
            "diagnostic-module completeness and a confidence tier (High/Medium/",
            "Low/Insufficient evidence) so partial and fully supported calls can",
            "be distinguished."
          ),
          assumptions = paste(
            "The taxonomic assignment and representative genome are suitable",
            "proxies for the sampled organism."
          ),
          implementation = paste(
            "NCBI taxonomy/RefSeq retrieval, GFF attribute search and local",
            "Triple_A classification functions."
          )
        )
      )
    ),
    top_abundance = list(
      analysis_name = "Top abundance",
      outputs = list(
        top_heatmap = list(
          output_name = "Top-taxa heatmap",
          example_image = "top_heatmap.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Compare the relative abundance of the most abundant taxa among samples.",
          preparation = "Taxa are ranked by mean abundance and the selected top N are retained.",
          method = "Descriptive heatmap of the selected abundance matrix.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable; the plot is descriptive.",
          interpretation = "Colour intensity represents relative abundance, not statistical evidence.",
          assumptions = "Abundance values are comparable among samples after the selected input scaling.",
          implementation = "Triple_A abundance ranking and heatmap export."
        ),
        top_lollipop = list(
          output_name = "Top-taxa lollipop plot",
          example_image = "top_lollipop.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Rank the most abundant taxa using their mean abundance.",
          preparation = "Taxa are ranked by mean abundance and the selected top N are retained.",
          method = "Descriptive lollipop chart.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable; the plot is descriptive.",
          interpretation = "Longer stems indicate greater mean abundance.",
          assumptions = "Mean abundance is an appropriate summary for ranking taxa.",
          implementation = "Triple_A descriptive abundance summary with ggplot2."
        ),
        top_composition = list(
          output_name = "Top-taxa composition plot",
          example_image = "top_composition.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Display the proportional composition of selected taxa across treatments or samples.",
          preparation = "The top N taxa are retained and displayed as stacked abundance components.",
          method = "Descriptive stacked composition chart.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable; the plot is descriptive.",
          interpretation = "Segment sizes represent composition and should not be interpreted as significant differences.",
          assumptions = "The displayed abundance scale is appropriate for compositional comparison.",
          implementation = "Triple_A abundance aggregation with ggplot2."
        ),
        top_distribution = list(
          output_name = "Top-taxa distribution plot",
          example_image = "top_distribution.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Visualise abundance variability of selected taxa across treatments.",
          preparation = "The top N taxa are retained and sample-level abundances are grouped by treatment.",
          method = "Descriptive sample-distribution plot.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable unless a separate differential-abundance analysis is run.",
          interpretation = "Overlap or separation is descriptive and does not by itself demonstrate significance.",
          assumptions = "Replicate labels and treatment assignments are correct.",
          implementation = "Triple_A grouped abundance display with ggplot2."
        )
      )
    ),
    differential_abundance = list(
      analysis_name = "Differential abundance",
      outputs = list(
        volcano_plot = list(
          output_name = "Volcano plot",
          example_image = "volcano_plot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Identify taxa differing between each pair of treatments.",
          preparation = paste(
            "Taxa are filtered by minimum prevalence and mean abundance;",
            "a pseudocount is used for log2 fold-change calculation."
          ),
          method = "Pairwise taxon-level hypothesis testing and log2 fold-change estimation.",
          statistical_test = "Configured at run time: Wilcoxon rank-sum test or Student t-test.",
          significance = paste(
            "Benjamini-Hochberg adjusted P value <= alpha AND absolute",
            "log2 fold change >= the configured threshold."
          ),
          interpretation = paste(
            "A taxon is labelled significant only when it meets both the",
            "adjusted-P and effect-size criteria."
          ),
          assumptions = paste(
            "Samples are independent unless paired=TRUE; the selected test",
            "and replicate design are appropriate."
          ),
          implementation = "stats::wilcox.test or stats::t.test; stats::p.adjust(method='BH')."
        ),
        ma_plot = list(
          output_name = "MA plot",
          example_image = "ma_plot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Display effect size against mean abundance for a pairwise comparison.",
          preparation = "Uses the same filtering, testing and fold-change calculations as the volcano plot.",
          method = "Mean abundance versus log2 fold-change display.",
          statistical_test = "The underlying significance classification uses the configured pairwise test.",
          significance = paste(
            "Benjamini-Hochberg adjusted P value <= alpha AND absolute",
            "log2 fold change >= the configured threshold."
          ),
          interpretation = "The MA plot is a visualisation; significance comes from the underlying test table.",
          assumptions = "The same assumptions as the selected differential-abundance test apply.",
          implementation = "Triple_A differential table with ggplot2."
        ),
        qq_plot = list(
          output_name = "QQ plot",
          example_image = "volcano_plot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Check whether the pairwise test's P-value distribution departs from the uniform expectation under the null.",
          preparation = "Uses the same filtering and testing as the volcano plot; expected quantiles come from the theoretical uniform(0,1) distribution of P-values.",
          method = "Observed versus expected -log10(P) quantile-quantile comparison, one plot per pairwise comparison.",
          statistical_test = "Diagnostic only; no additional hypothesis test is performed here.",
          significance = "Not applicable; use as a diagnostic for P-value inflation or deflation, not as a significance test itself.",
          interpretation = "Points following the diagonal support a well-calibrated test; systematic deviation suggests inflation, deflation, or model misspecification.",
          assumptions = "P-values from the underlying pairwise test are what is being diagnosed; degenerate designs can distort the expected quantiles.",
          implementation = "Triple_A P-value quantile comparison with ggplot2."
        )
      )
    ),
    functional_abundance = list(
      analysis_name = "Potential metabolomic pathways abundance",
      outputs = list(
        functional_abundance_heatmap = list(
          output_name = "Potential metabolomic pathways abundance heatmap",
          example_image = "potential_metabolomic_pathways_abundance_heatmap.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Compare the abundance of taxa assigned to selected putative functions.",
          preparation = "Taxa are selected using the curated functional classifier and aggregated by function.",
          method = "Descriptive heatmap of potential metabolomic pathway abundance.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable; the heatmap is descriptive.",
          interpretation = "Values combine taxonomic abundance with predicted reference-genome potential.",
          assumptions = "Functional assignment is an appropriate proxy for the sampled taxon.",
          implementation = "Triple_A functional registry, abundance aggregation and heatmap export."
        ),
        functional_abundance_barplot = list(
          output_name = "Potential metabolomic pathways abundance bar plot",
          example_image = "potential_metabolomic_pathways_abundance_barplot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Summarise potential metabolomic pathway abundance across treatments.",
          preparation = "Taxa classified within each selected function are aggregated.",
          method = "Descriptive bar chart.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable; apparent differences are descriptive.",
          interpretation = "Bar height represents inferred abundance, not measured pathway activity.",
          assumptions = "Reference-genome predictions adequately represent the reported taxa.",
          implementation = "Triple_A aggregation with ggplot2."
        ),
        functional_abundance_dotplot = list(
          output_name = "Potential metabolomic pathways abundance dot plot",
          example_image = "potential_metabolomic_pathways_abundance_dotplot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Compare inferred functions across treatments using compact abundance markers.",
          preparation = "Taxa classified within selected functions are aggregated.",
          method = "Descriptive dot plot.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable; the plot is descriptive.",
          interpretation = "Marker position and size reflect inferred abundance according to the plot legend.",
          assumptions = "The curated classifier is sufficiently specific for the biological question.",
          implementation = "Triple_A aggregation with ggplot2."
        ),
        functional_contributors = list(
          output_name = "Potential pathway contributors",
          example_image = "potential_metabolomic_pathway_contributors.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Identify the principal taxa contributing to each inferred function.",
          preparation = "Contributors are ranked by abundance within each selected function.",
          method = "Descriptive ranked table.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable; ranking does not imply significance.",
          interpretation = "The table identifies abundant predicted contributors, not causal organisms.",
          assumptions = "Taxonomic abundance and function assignment are suitable for contributor ranking.",
          implementation = "Triple_A functional contributor ranking."
        )
      )
    ),
    community_structure = list(
      analysis_name = "Community structure and diversity",
      outputs = list(
        pca_plot = list(
          output_name = "PCA sample plot",
          example_image = "pca_plot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Summarise the largest linear sources of variation among community profiles.",
          preparation = "Configured transformation is applied; taxa with zero variance are removed.",
          method = "Principal component analysis on the transformed abundance matrix.",
          statistical_test = "No hypothesis test; PCA is exploratory.",
          significance = "No P-value or significance threshold is used for PCA separation.",
          interpretation = paste(
            "Nearby samples have more similar transformed profiles. Axis",
            "labels report explained variance, but visual separation is not",
            "evidence of statistical significance."
          ),
          assumptions = "PCA represents linear variation and is sensitive to transformation and dominant taxa.",
          implementation = "stats::prcomp(center=TRUE, scale.=FALSE)."
        ),
        pca_scree = list(
          output_name = "PCA scree plot",
          example_image = "pca_scree.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Report the proportion of variance explained by successive principal components.",
          preparation = "Uses the eigenvalues from the fitted PCA.",
          method = "Explained-variance and cumulative-variance summary.",
          statistical_test = "No hypothesis test.",
          significance = "No universal significance threshold is applied.",
          interpretation = "Components explaining more variance contribute more strongly to the ordination.",
          assumptions = "The selected transformation provides an informative variance structure.",
          implementation = "Variance calculated from stats::prcomp singular values."
        ),
        pcoa_plot = list(
          output_name = "PCoA sample plot",
          example_image = "pcoa_plot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Represent pairwise beta-diversity distances in two dimensions.",
          preparation = "Relative abundances are converted to the configured Bray-Curtis or Jaccard distance.",
          method = "Principal coordinates analysis of the beta-diversity matrix.",
          statistical_test = "No hypothesis test; PCoA is exploratory.",
          significance = "No P-value is assigned to visual group separation.",
          interpretation = "Closer samples have smaller beta-diversity distances under the selected metric.",
          assumptions = "The selected distance metric reflects the ecological differences of interest.",
          implementation = "vegan::vegdist followed by stats::cmdscale(add=TRUE)."
        ),
        nmds_plot = list(
          output_name = "NMDS sample plot",
          example_image = "nmds_plot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Represent the rank order of beta-diversity distances in a low-dimensional space.",
          preparation = "Relative abundances are converted to the configured Bray-Curtis or Jaccard distance.",
          method = "Two-dimensional non-metric multidimensional scaling.",
          statistical_test = "No hypothesis test; NMDS is exploratory.",
          significance = "No P-value is used. Model adequacy is described by stress.",
          interpretation = paste(
            "Lower stress indicates a better two-dimensional representation;",
            "as a guide, <0.10 is good, 0.10-0.20 is usable with caution,",
            "and >0.20 indicates a weak representation."
          ),
          assumptions = "Rank-order distances can be represented adequately in two dimensions.",
          implementation = "vegan::metaMDS(k=2, autotransform=FALSE)."
        ),
        permanova_table = list(
          output_name = "PERMANOVA table",
          example_image = "permanova_table.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Test whether multivariate community centroids differ among treatments.",
          preparation = "Uses the same beta-diversity distance matrix as PCoA and NMDS.",
          method = "Permutational multivariate analysis of variance.",
          statistical_test = "PERMANOVA with the configured number of permutations.",
          significance = "Permutation P value <= the configured community alpha threshold.",
          interpretation = paste(
            "A significant result supports differences in multivariate",
            "centroids. It should be interpreted alongside group-dispersion",
            "diagnostics because unequal dispersion can affect PERMANOVA."
          ),
          assumptions = "Exchangeability under permutation and an appropriate experimental design.",
          implementation = "vegan::adonis2(distance ~ Treatment)."
        ),
        alpha_diversity_plot = list(
          output_name = "Alpha-diversity plot",
          example_image = "alpha_diversity_plot.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Display within-sample richness and diversity across treatments.",
          preparation = "Calculated from relative-abundance profiles for each sample.",
          method = "Observed taxa, Shannon, Simpson and inverse Simpson indices.",
          statistical_test = "No inferential test is currently applied to alpha-diversity groups.",
          significance = "Not applicable; the boxplots and points are descriptive.",
          interpretation = "Apparent group differences require a separate validated statistical test before being called significant.",
          assumptions = "Sampling depth and preprocessing permit meaningful diversity comparison.",
          implementation = "vegan::diversity and observed non-zero taxa counts."
        ),
        alpha_diversity_table = list(
          output_name = "Alpha-diversity table",
          example_image = "alpha_diversity_table.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Provide sample-level richness and diversity values.",
          preparation = "Calculated from relative-abundance profiles.",
          method = "Observed taxa, Shannon, Simpson and inverse Simpson indices.",
          statistical_test = "No inferential test.",
          significance = "Not applicable; the table contains descriptive indices.",
          interpretation = "Values should be compared in the context of sequencing depth and preprocessing.",
          assumptions = "The retained taxa and abundance values are suitable for diversity estimation.",
          implementation = "vegan::diversity."
        ),
        beta_distance_table = list(
          output_name = "Beta-distance matrix",
          example_image = "beta_distance_table.png",
          example_caption = paste(
            "Illustrative example generated with simulated data.",
            "It is not based on user data and is not an analytical result."
          ),
          objective = "Provide all pairwise ecological distances among samples.",
          preparation = "Relative abundances are converted using the configured distance metric.",
          method = "Bray-Curtis dissimilarity or binary Jaccard distance.",
          statistical_test = "No hypothesis test in the distance table itself.",
          significance = "Not applicable; significance is assessed separately by PERMANOVA.",
          interpretation = "Smaller values indicate more similar samples under the selected metric.",
          assumptions = "The selected distance metric is appropriate for the abundance data and biological question.",
          implementation = "vegan::vegdist."
        ),
        pairwise_permanova_heatmap = list(
          output_name = "Pairwise PERMANOVA heatmap",
          example_image = "",
          example_caption = paste(
            "Graphical summary of the pairwise PERMANOVA table.",
            "No additional statistical test is performed for the figure."
          ),
          objective = "Show, at a glance, which pairs of treatments differ in community composition.",
          preparation = "Uses the pairwise PERMANOVA table already computed for every pair of groups.",
          method = "Symmetric heatmap of the BH-adjusted permutation P-values (with R2 annotations).",
          statistical_test = "None for the figure; values come from vegan::adonis2 run per pair of groups.",
          significance = "Interpret against the configured community alpha; darker cells are more strongly separated.",
          interpretation = "Requires at least three treatment groups; otherwise an informative placeholder is drawn.",
          assumptions = "Same exchangeability assumptions as the underlying pairwise PERMANOVA.",
          implementation = "ggplot2 tile plot of the Pairwise_PERMANOVA table."
        ),
        beta_dispersion_boxplot = list(
          output_name = "Beta-dispersion boxplot",
          example_image = "",
          example_caption = paste(
            "Graphical summary of the beta-dispersion analysis.",
            "No additional statistical test is performed for the figure."
          ),
          objective = "Assess within-group multivariate dispersion, an assumption of PERMANOVA.",
          preparation = "Uses per-sample distances to the group centroid from the same betadisper fit that produces the test.",
          method = "Boxplot of distance-to-centroid per group with individual points and group means.",
          statistical_test = "The permutation test is reported in the Beta_dispersion table; the figure is descriptive.",
          significance = "Unequal dispersion can confound PERMANOVA; read alongside the PERMANOVA result.",
          interpretation = "Similar spread across groups supports the PERMANOVA homogeneity assumption.",
          assumptions = "Distances to centroid are computed in the PCoA space of the selected distance metric.",
          implementation = "vegan::betadisper distances rendered with ggplot2."
        ),
        anosim_plot = list(
          output_name = "ANOSIM R gauge",
          example_image = "",
          example_caption = paste(
            "Graphical summary of the ANOSIM statistic.",
            "No additional statistical test is performed for the figure."
          ),
          objective = "Communicate the strength of community separation from ANOSIM.",
          preparation = "Uses the ANOSIM R statistic and permutation P-value already computed.",
          method = "Gauge placing the R statistic on the -1 to 1 scale, coloured by effect size.",
          statistical_test = "vegan::anosim with the configured number of permutations (reported in the ANOSIM table).",
          significance = "Permutation P value <= the configured community alpha.",
          interpretation = "R near 1 indicates strong separation; R near 0 indicates little separation.",
          assumptions = "ANOSIM is rank-based and sensitive to unequal group sizes.",
          implementation = "ggplot2 gauge of the ANOSIM table."
        ),
        permanova_variance_plot = list(
          output_name = "PERMANOVA explained-variance plot",
          example_image = "",
          example_caption = paste(
            "Graphical summary of the PERMANOVA variance decomposition.",
            "No additional statistical test is performed for the figure."
          ),
          objective = "Show the proportion of community variation explained by the model term(s).",
          preparation = "Uses the R2 column from the PERMANOVA table.",
          method = "Horizontal bar chart of R2 by model term (Total excluded).",
          statistical_test = "None for the figure; R2 comes from vegan::adonis2.",
          significance = "Use the PERMANOVA permutation P-value for inference; the bars are descriptive.",
          interpretation = "Longer bars indicate more community variation attributed to that term.",
          assumptions = "Same assumptions as the underlying PERMANOVA model.",
          implementation = "ggplot2 bar chart of the PERMANOVA R2 values."
        ),
        dendrogram_plot = list(
          output_name = "Hierarchical clustering dendrogram",
          example_image = "",
          example_caption = paste(
            "Graphical summary of the hierarchical clustering.",
            "No additional statistical test is performed for the figure."
          ),
          objective = "Show which samples join first when the community distances are grouped hierarchically.",
          preparation = "Uses the same transformed abundance matrix and beta-diversity distance as the ordinations, so the dendrogram and the PCoA cannot disagree.",
          method = "Agglomerative hierarchical clustering with average linkage (UPGMA); leaves are coloured by declared treatment.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable. Clustering is exploratory and always produces groups, whether or not they are real.",
          interpretation = "Samples joining at a low height are compositionally similar. Read the height at which branches merge, not the left-to-right order of the leaves, which is arbitrary among equivalent rotations. Agreement between the branches and the declared treatments is suggestive, never a test: use PERMANOVA for that.",
          assumptions = "Average linkage does not assume Euclidean geometry, which is why it suits ecological dissimilarities. A different linkage can produce a different tree from the same distances.",
          implementation = "stats::hclust on the vegan::vegdist distance matrix, rendered with ggplot2."
        ),
        hierarchical_clustering_table = list(
          output_name = "Hierarchical clustering table",
          example_image = "",
          example_caption = "Cluster membership derived from the dendrogram.",
          objective = "State, per sample, the cluster it falls into and the treatment it was declared as.",
          preparation = "Cuts the dendrogram at the number of declared treatment groups.",
          method = "stats::cutree at k equal to the number of declared groups, bounded by the number of samples.",
          statistical_test = "No inferential statistical test.",
          significance = "Not applicable; the table is descriptive.",
          interpretation = "Compare the Cluster column against the Treatment column: agreement indicates that an unsupervised grouping recovers the experimental design, disagreement that it does not. Neither outcome is evidence on its own.",
          assumptions = "Cutting at the number of declared groups assumes that is a sensible number of clusters, which the data may not support.",
          implementation = "stats::cutree over the stats::hclust result."
        )
      )
    )
  )
}

aaa_resolve_method_parameters <- function(
  analysis_id,
  output_id,
  parameters = list()
) {
  defaults <- aaa_get_defaults()

  diff_cfg <- utils::modifyList(
    defaults$differential_abundance,
    parameters$differential_abundance %||% list()
  )

  community_cfg <- utils::modifyList(
    defaults$community_structure,
    parameters$community_structure %||% list()
  )

  if (analysis_id == "differential_abundance") {
    test_name <- if (identical(diff_cfg$method, "wilcox")) {
      "Wilcoxon rank-sum test"
    } else {
      "Student t-test"
    }

    return(list(
      statistical_test = paste0(
        test_name,
        if (isTRUE(diff_cfg$paired)) " (paired)" else " (unpaired)"
      ),
      significance = paste0(
        "Benjamini-Hochberg adjusted P <= ",
        format(diff_cfg$alpha),
        " and |log2FC| >= ",
        format(diff_cfg$log2fc_threshold),
        "."
      ),
      run_parameters = paste0(
        "method=", diff_cfg$method,
        "; paired=", diff_cfg$paired,
        "; min_prevalence=", diff_cfg$min_prevalence,
        "; min_mean_abundance=", diff_cfg$min_mean_abundance,
        "; alpha=", diff_cfg$alpha,
        "; log2fc_threshold=", diff_cfg$log2fc_threshold
      )
    ))
  }

  if (analysis_id == "community_structure") {
    test <- if (output_id == "permanova_table") {
      paste0(
        "PERMANOVA (vegan::adonis2; ",
        as.integer(community_cfg$permutations),
        " permutations)"
      )
    } else if (output_id %in% c("pca_plot", "pca_scree")) {
      "No hypothesis test (exploratory PCA)"
    } else if (output_id == "nmds_plot") {
      "No hypothesis test (ordination quality assessed by stress)"
    } else if (output_id == "pcoa_plot") {
      "No hypothesis test (exploratory PCoA)"
    } else {
      "No inferential test"
    }

    criterion <- if (output_id == "permanova_table") {
      paste0(
        "Permutation P <= ",
        format(community_cfg$significance_alpha),
        "."
      )
    } else {
      "Not applicable; this output is descriptive or exploratory."
    }

    return(list(
      statistical_test = test,
      significance = criterion,
      run_parameters = paste0(
        "transformation=", community_cfg$transformation,
        "; distance=", community_cfg$distance_method,
        "; permutations=", community_cfg$permutations,
        "; significance_alpha=", community_cfg$significance_alpha,
        "; nmds_trymax=", community_cfg$nmds_trymax
      )
    ))
  }

  list(
    statistical_test = NULL,
    significance = NULL,
    run_parameters = "No additional inferential parameters."
  )
}

aaa_methodology_table <- function(
  analyses,
  outputs,
  parameters = list()
) {
  registry <- aaa_analysis_method_registry()
  output_catalogue <- aaa_output_catalogue()

  selected <- output_catalogue[
    output_catalogue$ID %in% outputs &
      output_catalogue$Module %in% analyses, ,
    drop = FALSE
  ]

  rows <- lapply(
    seq_len(nrow(selected)),
    function(i) {
      analysis_id <- selected$Module[i]
      output_id <- selected$ID[i]
      definition <- registry[[analysis_id]]$outputs[[output_id]]

      if (is.null(definition)) {
        return(NULL)
      }

      resolved <- aaa_resolve_method_parameters(
        analysis_id = analysis_id,
        output_id = output_id,
        parameters = parameters
      )

      data.frame(
        Analysis_ID = analysis_id,
        Output_ID = output_id,
        Analysis = registry[[analysis_id]]$analysis_name,
        Output = definition$output_name,
        Example_image = definition$example_image %||% "",
        Example_caption = definition$example_caption %||% "",
        Objective = definition$objective,
        Data_preparation = definition$preparation,
        Method = definition$method,
        Statistical_test = resolved$statistical_test %||%
          definition$statistical_test,
        Significance_criterion = resolved$significance %||%
          definition$significance,
        Interpretation = definition$interpretation,
        Assumptions_and_limitations = definition$assumptions,
        Implementation = definition$implementation,
        Run_parameters = resolved$run_parameters,
        stringsAsFactors = FALSE
      )
    }
  )

  dplyr::bind_rows(rows)
}

aaa_write_methodology_files <- function(
  methodology,
  metadata_dir
) {
  workbook <- file.path(
    metadata_dir,
    "Analysis_methods.xlsx"
  )

  text_file <- file.path(
    metadata_dir,
    "Analysis_methods.txt"
  )

  openxlsx::write.xlsx(
    list(Analysis_methods = methodology),
    workbook,
    overwrite = TRUE
  )

  lines <- unlist(lapply(
    seq_len(nrow(methodology)),
    function(i) {
      row <- methodology[i, , drop = FALSE]
      c(
        paste0("ANALYSIS: ", row$Analysis),
        paste0("OUTPUT: ", row$Output),
        paste0("OBJECTIVE: ", row$Objective),
        paste0("DATA PREPARATION: ", row$Data_preparation),
        paste0("METHOD: ", row$Method),
        paste0("STATISTICAL TEST: ", row$Statistical_test),
        paste0("SIGNIFICANCE: ", row$Significance_criterion),
        paste0("INTERPRETATION: ", row$Interpretation),
        paste0("ASSUMPTIONS/LIMITATIONS: ", row$Assumptions_and_limitations),
        paste0("IMPLEMENTATION: ", row$Implementation),
        paste0("RUN PARAMETERS: ", row$Run_parameters),
        strrep("-", 80),
        ""
      )
    }
  ))

  writeLines(lines, text_file, useBytes = TRUE)

  list(
    workbook = workbook,
    text = text_file
  )
}


aaa_analysis_method_registry <- function() {
  registry <- aaa_analysis_method_registry_base()
  registry$plsda <- list(
    analysis_name = "PLS-DA",
    outputs = list(
      plsda_plot = list(
        output_name = "PLS-DA score plot", example_image = "plsda_plot.png",
        example_caption = "Illustrative supervised ordination generated with simulated data.",
        objective = "Assess whether treatment groups can be discriminated from transformed community profiles.",
        preparation = "Abundances are transformed, centred and scaled; treatment labels define the response classes.",
        method = "Partial least-squares discriminant analysis with dummy-coded classes.",
        statistical_test = "Supervised model; predictive performance is assessed by cross-validation.",
        significance = "No P-value is assigned to visual separation. Report cross-validated accuracy and avoid interpretation when performance is weak.",
        interpretation = "Samples close in score space have similar supervised latent profiles; separation must be supported by validation.",
        assumptions = "Groups are predefined, sample size is adequate, and cross-validation is representative.",
        implementation = "pls::plsr with k-fold cross-validation."
      ),
      plsda_confusion_plot = list(
        output_name = "PLS-DA confusion matrix", example_image = "plsda_performance.png",
        example_caption = "Illustrative confusion matrix generated with simulated data.",
        objective = "Visualise how predicted classes compare with the declared groups across held-out folds.",
        preparation = "Predictions are pooled across cross-validation folds.",
        method = "Observed-versus-predicted class contingency table rendered as a confusion matrix.",
        statistical_test = "Cross-validated classification agreement; no P-value is computed.",
        significance = "Read alongside overall accuracy; a clean diagonal alone does not establish significance.",
        interpretation = "Off-diagonal cells identify which groups the model confuses most often.",
        assumptions = "Folds preserve representative class structure and are not leaked between training and validation.",
        implementation = "pls validation predictions rendered as a confusion matrix."
      ),
      plsda_vip_plot = list(
        output_name = "PLS-DA variable importance (VIP)", example_image = "plsda_loadings.png",
        example_caption = "Illustrative VIP ranking generated with simulated data.",
        objective = "Identify taxa contributing strongly to the supervised latent components.",
        preparation = "VIP scores are computed from the fitted model's X weights and explained Y variance.",
        method = "Standard Variable Importance in Projection (VIP) approximation, ranked in decreasing order.",
        statistical_test = "No hypothesis test; VIP > 1 is a common, non-formal screening threshold.",
        significance = "VIP indicates model contribution, not univariate significance.",
        interpretation = "Taxa with high VIP scores are the strongest drivers of the supervised separation.",
        assumptions = "The fitted model is stable and not overfit.", implementation = "Computed from pls::plsr X weights and Y-variance explained."
      ),
      plsda_performance = list(
        output_name = "PLS-DA cross-validation performance", example_image = "plsda_performance.png",
        example_caption = "Illustrative confusion matrix generated with simulated data.",
        objective = "Quantify out-of-fold classification performance.", preparation = "Predictions are obtained from cross-validation.",
        method = "Observed-versus-predicted class comparison.", statistical_test = "Cross-validated classification accuracy.",
        significance = "Accuracy is descriptive; compare with class balance and chance expectation.", interpretation = "High accuracy supports predictive discrimination, but external validation is preferable.",
        assumptions = "Folds preserve representative class structure.", implementation = "pls validation predictions."
      ),
      plsda_loadings = list(
        output_name = "PLS-DA taxon loadings", example_image = "plsda_loadings.png",
        example_caption = "Illustrative loading ranking generated with simulated data.",
        objective = "Identify taxa contributing strongly to supervised latent components.", preparation = "Loadings are extracted from the fitted PLS model.",
        method = "Ranking by Euclidean loading magnitude across the first two components.", statistical_test = "No hypothesis test.",
        significance = "Loadings indicate model contribution, not univariate significance.", interpretation = "Large absolute loadings identify influential taxa that require biological validation.",
        assumptions = "The fitted model is stable and not overfit.", implementation = "pls::loadings."
      )
    )
  )
  registry$splsda <- list(
    analysis_name = "sPLS-DA",
    outputs = list(
      splsda_plot = list(
        output_name = "sPLS-DA score plot", example_image = "plsda_plot.png",
        example_caption = "Illustrative sparse supervised ordination.",
        objective = "Assess whether a reduced microbial signature discriminates predefined groups.",
        preparation = "The canonical transformed abundance matrix is filtered, centred and scaled before sparse modelling.",
        method = "Sparse partial least-squares discriminant analysis with component-specific keepX selection.",
        statistical_test = "Repeated stratified M-fold cross-validation reports balanced error rate.",
        significance = "Visual separation is not a significance test; interpret only alongside cross-validated BER.",
        interpretation = "Separated samples and a low BER support discrimination by the selected microbial signature.",
        assumptions = "Groups are predefined, classes contain biological replicates and validation folds are representative.",
        implementation = "mixOmics::splsda, mixOmics::tune.splsda and mixOmics::perf."
      ),
      splsda_confusion_plot = list(
        output_name = "sPLS-DA confusion matrix", example_image = "plsda_performance.png",
        example_caption = "Illustrative confusion matrix generated with simulated data.",
        objective = "Visualise how predicted classes compare with the declared groups across held-out folds.",
        preparation = "Predictions are pooled across repeated stratified M-fold cross-validation.",
        method = "Observed-versus-predicted class contingency table rendered as a confusion matrix.",
        statistical_test = "Cross-validated classification agreement; no P-value is computed.",
        significance = "Read alongside the balanced error rate; a clean diagonal alone does not establish significance.",
        interpretation = "Off-diagonal cells identify which groups the sparse signature confuses most often.",
        assumptions = "Folds preserve representative class structure and are not leaked between training and validation.",
        implementation = "mixOmics::perf confusion-matrix summaries."
      ),
      splsda_performance = list(
        output_name = "sPLS-DA validation performance", example_image = "plsda_performance.png",
        example_caption = "Illustrative classification performance.",
        objective = "Quantify predictive error across repeated held-out folds.", preparation = "Class-stratified repeated M-fold validation.",
        method = "Balanced error rate and observed-versus-predicted summaries.", statistical_test = "Cross-validated prediction error.",
        significance = "BER is descriptive and should be compared with chance and class imbalance.", interpretation = "Lower BER indicates better balanced prediction.",
        assumptions = "No sample leakage across folds.", implementation = "mixOmics::perf."
      ),
      splsda_selected_features = list(
        output_name = "sPLS-DA selected taxa", example_image = "plsda_loadings.png",
        example_caption = "Illustrative sparse loading ranking.",
        objective = "Identify taxa retained in the sparse discriminant signature.", preparation = "Features are selected separately on each latent component.",
        method = "Absolute sparse loading ranking.", statistical_test = "No univariate hypothesis test.",
        significance = "Selection indicates multivariate model contribution, not taxon-level significance.", interpretation = "Repeatedly selected taxa are stronger signature candidates.",
        assumptions = "The tuned model is stable.", implementation = "mixOmics::selectVar."
      ),
      splsda_stability = list(
        output_name = "sPLS-DA feature-selection stability", example_image = "plsda_loadings.png",
        example_caption = "Illustrative selection-stability summary.",
        objective = "Summarise how consistently taxa contribute across model components.", preparation = "Selected taxa are collated from component-specific sparse loadings.",
        method = "Component selection frequency.", statistical_test = "No P-value.", significance = "Use as a robustness diagnostic.",
        interpretation = "Higher frequency indicates a more persistent multivariate contribution.", assumptions = "Components are interpretable and tuning is adequate.", implementation = "Derived from mixOmics::selectVar outputs."
      )
    )
  )
  registry$rda <- list(
    analysis_name = "RDA with environmental variables",
    outputs = list(
      rda_plot = list(
        output_name = "RDA ordination plot", example_image = "rda_plot.png",
        example_caption = "Illustrative constrained ordination generated with simulated data.",
        objective = "Relate community composition to selected environmental variables.",
        preparation = "Environmental rows are matched by sample ID; selected numeric variables are standardised.",
        method = "Redundancy analysis, a constrained linear ordination.", statistical_test = "Permutation tests are reported separately.",
        significance = "Visual alignment is descriptive; significance is determined by permutation tests.",
        interpretation = "Arrow direction indicates increasing environmental values; arrow length reflects association strength.",
        assumptions = "Relationships are approximately linear and environmental predictors are not excessively collinear.", implementation = "vegan::rda and vegan::scores."
      ),
      rda_anova = list(
        output_name = "RDA permutation tests", example_image = "rda_anova.png",
        example_caption = "Illustrative permutation-test table generated with simulated data.",
        objective = "Test the overall constrained model and individual environmental terms.", preparation = "Uses the fitted RDA model.",
        method = "Permutation ANOVA for constrained ordination.", statistical_test = "vegan::anova.cca overall and by term.",
        significance = "Permutation P <= configured alpha.", interpretation = "A significant term explains community variation beyond random permutation expectation.",
        assumptions = "Exchangeability under permutation and an appropriate environmental model.", implementation = "vegan::anova.cca."
      ),
      rda_variance = list(
        output_name = "RDA explained variance", example_image = "rda_variance.png",
        example_caption = "Illustrative constrained-variance summary generated with simulated data.",
        objective = "Summarise the variance represented by constrained axes.", preparation = "Eigenvalues are extracted from the constrained component.",
        method = "Percentage of constrained inertia by RDA axis.", statistical_test = "No additional test.",
        significance = "Use permutation tests for inference.", interpretation = "Higher percentages indicate more constrained variation represented by an axis.",
        assumptions = "The RDA model is appropriate.", implementation = "vegan constrained eigenvalues."
      ),
      rda_variance_plot = list(
        output_name = "RDA explained variance (plot)", example_image = "rda_variance.png",
        example_caption = "Illustrative constrained-variance bar plot generated with simulated data.",
        objective = "Visualise the variance represented by each constrained axis.", preparation = "Eigenvalues are extracted from the constrained component.",
        method = "Percentage of constrained inertia by RDA axis, rendered as a bar plot.", statistical_test = "No additional test.",
        significance = "Use permutation tests for inference.", interpretation = "Taller bars indicate more constrained variation represented by that axis.",
        assumptions = "The RDA model is appropriate.", implementation = "vegan constrained eigenvalues."
      )
    )
  )

  registry$differential_functions <- list(
      analysis_name = "Differential functions",
      outputs = list(
        differential_functions_plot = list(
          output_name = "Differential functions plot",
          example_image = "", example_caption = "Graphical summary of the differential-function test.",
          objective = "Show which inferred functional profiles differ between experimental groups.",
          preparation = "Uses the pathway-by-sample matrix produced by the potential metabolomic pathways abundance analysis; no abundance is recomputed.",
          method = "Bar chart of the log2 fold change per pathway and comparison, coloured by significance.",
          statistical_test = "Wilcoxon rank-sum test per pathway for every pair of treatments.",
          significance = "Benjamini-Hochberg adjusted P below the configured alpha.",
          interpretation = "A significant pathway means its inferred genomic potential is differentially represented between the groups. It does not mean the pathway is more active: no expression is measured.",
          assumptions = "Pathway abundances inherit the skew of the taxon abundances they are summed from, which is why a non-parametric test is used. Pathways that are constant across all samples are reported as untested rather than given a P value.",
          implementation = "stats::wilcox.test with stats::p.adjust(method = 'BH'), rendered with ggplot2."
        ),
        differential_functions_table = list(
          output_name = "Differential functions table",
          example_image = "", example_caption = "Per-pathway comparison statistics.",
          objective = "Report the full comparison statistics for every pathway and pair of groups.",
          preparation = "Same pathway-by-sample matrix as the figure.",
          method = "One row per pathway and comparison with group means, difference, log2 fold change and adjusted P.",
          statistical_test = "Wilcoxon rank-sum test per pathway.",
          significance = "Benjamini-Hochberg adjusted P below the configured alpha.",
          interpretation = "Read the adjusted P together with the difference in means: a small P on a negligible difference is not biologically meaningful.",
          assumptions = "The Tested column records whether a pathway could be tested at all.",
          implementation = "Triple_A comparison table exported to the analysis workbook."
        )
      )
    )

  registry$functional_enrichment <- list(
      analysis_name = "Functional enrichment",
      outputs = list(
        functional_enrichment_plot = list(
          output_name = "Functional enrichment plot",
          example_image = "", example_caption = "Graphical summary of the enrichment test.",
          objective = "Show which curated functions are over-represented among the differentially abundant taxa.",
          preparation = "Crosses the per-taxon functional calls with the taxa flagged as significant by the differential-abundance analysis.",
          method = "Dot plot of adjusted significance per function, sized by the number of overlapping taxa.",
          statistical_test = "Fisher's exact test per function.",
          significance = "Benjamini-Hochberg adjusted P below the configured alpha, with an odds ratio above one.",
          interpretation = "An enriched function means the taxa carrying it are over-represented among those that changed, which suggests a functional theme rather than isolated taxa. It remains inferred potential, not measured activity.",
          assumptions = "The background is every taxon that received a functional call, not every taxon in the dataset: a taxon with no reference genome could never have been called positive, and including it would inflate the apparent enrichment.",
          implementation = "stats::fisher.test with stats::p.adjust(method = 'BH'), rendered with ggplot2."
        ),
        functional_enrichment_table = list(
          output_name = "Functional enrichment table",
          example_image = "", example_caption = "Per-function enrichment statistics.",
          objective = "Report the contingency counts and test statistics behind each enrichment call.",
          preparation = "Same functional calls and differential-abundance results as the figure.",
          method = "One row per function with the number of taxa carrying it, the background size, the overlap with significant taxa, the expected overlap, the odds ratio and the adjusted P.",
          statistical_test = "Fisher's exact test per function.",
          significance = "Benjamini-Hochberg adjusted P below the configured alpha.",
          interpretation = "Compare Overlap against Expected: an odds ratio above one with few overlapping taxa is fragile evidence, however small the P value.",
          assumptions = "Enrichment depends entirely on how the background is defined; the background used here is stated in the preparation field.",
          implementation = "Triple_A enrichment table exported to the analysis workbook."
        )
      )
  )

  registry
}
