# Every selectable output ID in aaa_output_catalogue() must have a matching
# help/documentation entry in aaa_analysis_method_registry(); otherwise the
# methodology table silently drops that row instead of erroring, which is
# how 4 gaps (splsda_confusion_plot, plsda_confusion_plot, plsda_vip_plot,
# rda_variance_plot, differential_abundance's qq_plot) went unnoticed.
qa_register_test(
  "CORE_007", "regression", "high",
  "Every catalogued output has a matching methodology-registry documentation entry",
  function() {
    catalogue <- aaa_output_catalogue()
    registry <- aaa_analysis_method_registry()

    missing <- character()
    for (i in seq_len(nrow(catalogue))) {
      module_id <- catalogue$Module[i]
      output_id <- catalogue$ID[i]
      definition <- registry[[module_id]]$outputs[[output_id]]
      if (is.null(definition)) {
        missing <- c(missing, paste0(module_id, "::", output_id))
      }
    }

    qa_expect_true(
      length(missing) == 0L,
      paste0("Catalogued outputs with no registry documentation entry: ", paste(missing, collapse = ", "))
    )
    TRUE
  }
)

# Every selectable analysis in aaa_analysis_catalogue() must have an entry in the
# dependency map, and vice versa. aaa_required_packages() rejects any selected
# analysis missing from TRIPLE_A_DEPENDENCIES$analyses as "Unknown analyses",
# which blocked runs that selected differential_functions or functional_enrichment
# (they were catalogued and runnable but absent from the dependency map).
qa_register_test(
  "CORE_014", "regression", "critical",
  "Every catalogued analysis is resolvable by the dependency registry",
  function() {
    catalogue_ids <- aaa_analysis_catalogue()$ID
    dependency_ids <- names(TRIPLE_A_DEPENDENCIES$analyses)

    missing_from_deps <- setdiff(catalogue_ids, dependency_ids)
    orphan_in_deps <- setdiff(dependency_ids, catalogue_ids)

    qa_expect_true(
      length(missing_from_deps) == 0L,
      paste0("Catalogued analyses absent from the dependency map (would fail with 'Unknown analyses'): ",
             paste(missing_from_deps, collapse = ", "))
    )
    qa_expect_true(
      length(orphan_in_deps) == 0L,
      paste0("Dependency-map analyses not present in the catalogue: ", paste(orphan_in_deps, collapse = ", "))
    )

    # The resolver must actually run for the full catalogue without erroring.
    resolved <- tryCatch(
      { aaa_required_packages(analyses = catalogue_ids); TRUE },
      error = function(e) conditionMessage(e)
    )
    qa_expect_true(isTRUE(resolved),
      paste0("aaa_required_packages() failed for the full catalogue: ", resolved))
    TRUE
  }
)
