# Differential functions and functional enrichment were offered as selectable
# analyses in the interface while no implementation existed: both checkboxes
# produced nothing and neither had a catalogued output. These tests exercise the
# engines directly with synthetic inputs, so they do not depend on the reference
# cache or on network access.

qa_register_test(
  "SCI_035", "regression", "high",
  "Differential functions compares pathway profiles between groups and adjusts for multiple testing",
  function() {
    set.seed(11)
    samples <- c(paste0("Ctrl_", 1:4), paste0("Trt_", 1:4))
    design <- data.frame(
      Sample_column = samples,
      Treatment = rep(c("Control", "Treated"), each = 4L),
      Replicate = rep(1:4, times = 2L),
      stringsAsFactors = FALSE
    )
    # One pathway carries a clear group difference, one is pure noise and one is
    # constant: the third must be reported as untested rather than given a P value.
    pathways <- data.frame(
      Pathway = c("Methanogenesis", "Noise_pathway", "Constant_pathway"),
      stringsAsFactors = FALSE
    )
    values <- rbind(
      c(runif(4, 0.05, 0.10), runif(4, 0.60, 0.70)),
      runif(8, 0.20, 0.30),
      rep(0.10, 8)
    )
    colnames(values) <- samples
    pathway_by_sample <- cbind(pathways, as.data.frame(values, stringsAsFactors = FALSE))

    result <- aaa_differential_functions(
      pathway_by_sample = pathway_by_sample,
      sample_design = design,
      alpha = 0.05,
      project_dir = tempfile("triplea_difffun_")
    )

    qa_expect_true(inherits(result, "Triple_A_result"), "Differential functions did not return a Triple_A_result.")
    qa_expect_true(
      !is.null(result$files[["differential_functions"]]) &&
        file.exists(result$files[["differential_functions"]]),
      "Differential functions did not save its figure."
    )

    tb <- result$tables$comparisons
    qa_expect_true(is.data.frame(tb) && nrow(tb) == 3L, "Expected one row per pathway for a two-group design.")
    qa_expect_true(
      all(c("Pathway", "log2FC", "P_value", "Adjusted_P", "Significance", "Tested") %in% names(tb)),
      "The comparison table is missing required columns."
    )

    constant <- tb[tb$Pathway == "Constant_pathway", , drop = FALSE]
    qa_expect_true(
      isFALSE(constant$Tested[1]) && is.na(constant$P_value[1]),
      "A pathway that is constant across all samples must be reported as untested, not given a P value."
    )
    # BH can only ever raise a P value, never lower it.
    tested <- tb[isTRUE(tb$Tested) | tb$Tested, , drop = FALSE]
    tested <- tested[!is.na(tested$P_value), , drop = FALSE]
    qa_expect_true(
      all(tested$Adjusted_P >= tested$P_value - 1e-12),
      "Adjusted P values are smaller than the raw P values, so the correction was not applied."
    )
    qa_expect_true(
      tb$log2FC[tb$Pathway == "Methanogenesis"] > 0,
      "The pathway that is more abundant in the second group should have a positive log2 fold change."
    )
    TRUE
  }
)

qa_register_test(
  "SCI_036", "regression", "high",
  "Functional enrichment tests curated functions against a background of taxa with a functional call",
  function() {
    # Two taxa carry the function and both are differentially abundant; the
    # remaining background taxa are neither. The overlap should therefore be 2.
    taxa_a <- data.frame(
      Taxonomy = paste0("Taxon_", 1:6),
      Potential = c(
        "Methanogenesis potential", "Methanogenesis potential",
        "No methanogenesis potential", "No methanogenesis potential",
        "Insufficient evidence", "No methanogenesis potential"
      ), stringsAsFactors = FALSE
    )
    functional_potential <- list(
      Methanogenesis = list(tables = list(taxa = taxa_a))
    )
    differential <- list(tables = list(comparisons = data.frame(
      Taxonomy = paste0("Taxon_", 1:6),
      Significance = c("Higher in Treated", "Higher in Treated",
                       "Not significant", "Not significant",
                       "Not significant", "Not significant"),
      stringsAsFactors = FALSE
    )))

    result <- aaa_functional_enrichment(
      functional_potential = functional_potential,
      differential_abundance = differential,
      alpha = 0.05,
      project_dir = tempfile("triplea_enrich_")
    )

    qa_expect_true(inherits(result, "Triple_A_result"), "Functional enrichment did not return a Triple_A_result.")
    qa_expect_true(
      !is.null(result$files[["functional_enrichment"]]) &&
        file.exists(result$files[["functional_enrichment"]]),
      "Functional enrichment did not save its figure."
    )

    tb <- result$tables$enrichment
    qa_expect_true(is.data.frame(tb) && nrow(tb) == 1L, "Expected one row per curated function.")
    qa_expect_true(
      all(c("Function", "Taxa_with_function", "Taxa_in_background", "Overlap",
            "Odds_ratio", "P_value", "Adjusted_P", "Enriched") %in% names(tb)),
      "The enrichment table is missing required columns."
    )
    qa_expect_true(tb$Overlap[1] == 2L, "The overlap between positive and significant taxa was not counted correctly.")
    qa_expect_true(tb$Taxa_with_function[1] == 2L, "Only calls that are neither negative nor insufficient count as positive.")
    # The background is every taxon with a functional call, including the
    # insufficient-evidence one, which could not have been called positive.
    qa_expect_true(
      tb$Taxa_in_background[1] == 6L,
      "The background must be every taxon that received a functional call."
    )
    qa_expect_true(is.finite(tb$P_value[1]), "Fisher's exact test did not produce a P value.")
    TRUE
  }
)
