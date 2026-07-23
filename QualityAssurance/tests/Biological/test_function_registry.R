qa_register_test(
  "BIO_001", "regression", "critical",
  "Biological functions have complete and valid registry definitions",
  function() {
    registry <- aaa_registry()
    required_fields <- c("display_name", "category", "description", "gene_roles",
                         "genes", "classifier", "graph_main", "analysis_name",
                         "pathway_selector")
    qa_expect_true(length(registry) >= 20L, "Too few registered functions")
    qa_expect_true(length(unique(names(registry))) == length(registry), "Duplicated registry IDs")
    qa_expect_true(all(vapply(registry, function(x) all(required_fields %in% names(x)), logical(1))), "Missing registry fields")
    qa_expect_true(all(vapply(registry, function(x) is.function(x$classifier), logical(1))), "Invalid classifier")
    qa_expect_true(all(vapply(registry, function(x) length(x$genes) > 0L, logical(1))), "Function without genes")
    validation <- aaa_validate_biological_registry(registry, strict = FALSE)
    qa_expect_true(isTRUE(validation$valid), paste(validation$issues, collapse="; "))
    TRUE
  }
)

qa_register_test(
  "BIO_002", "regression", "critical",
  "Functional evidence summary adds module completeness and a confidence tier",
  function() {
    # A generic, additive confidence layer sits over the existing Potential
    # classification. Diagnostic genes are the essential module; completeness is
    # the fraction detected, and confidence also accounts for how much of the
    # module could be evaluated. This exercises the five decision boundaries.
    taxa <- data.frame(
      Taxonomy = paste0("Tax", 1:5), Genus = paste0("G", 1:5), Tax_level = "Genus",
      Potential = NA_character_,
      mcrA = c(TRUE,  TRUE,  FALSE, TRUE, NA),
      mcrB = c(TRUE,  FALSE, FALSE, NA,   NA),
      mcrG = c(TRUE,  FALSE, FALSE, NA,   NA),
      stringsAsFactors = FALSE
    )
    ev <- aaa_summarize_function_evidence(taxa, "Methanogenesis")
    qa_expect_true(
      all(c("Diagnostic_completeness", "Confidence") %in% names(ev)),
      "Evidence summary is missing the new completeness/confidence columns."
    )
    qa_expect_true(
      isTRUE(all.equal(
        ev$Diagnostic_completeness, round(c(1, 1/3, 0, 1/3, 0), 3), tolerance = 1e-6
      )),
      "Diagnostic completeness was not the detected fraction of the diagnostic module."
    )
    qa_expect_equal(
      ev$Confidence,
      c("High confidence", "Low confidence", "No evidence detected",
        "Insufficient evidence", "Insufficient evidence"),
      "Confidence tiers did not match the expected module-completeness boundaries."
    )
    TRUE
  }
)

qa_register_test(
  "BIO_003", "regression", "critical",
  "Biogeochemical cycles are registered with conservative classifiers",
  function() {
    registry <- aaa_registry()
    new_functions <- c(
      "Phosphonate_degradation", "Perchlorate_reduction", "Organohalide_respiration",
      "Dissimilatory_iron_reduction", "Iron_oxidation", "Arsenate_respiration",
      "Arsenite_oxidation", "Arsenic_detoxification", "Periplasmic_nitrate_reduction"
    )
    qa_expect_true(all(new_functions %in% names(registry)),
      "One or more new biogeochemical cycles are missing from the registry.")

    marker <- function(genes, present) stats::setNames(genes %in% present, genes)

    # Iron reduction/oxidation must never over-call a complete pathway from
    # lineage-specific cytochromes; the strongest call is "Probable".
    fe_red <- registry[["Dissimilatory_iron_reduction"]]
    qa_expect_true(
      grepl("^Probable", fe_red$classifier(marker(fe_red$genes, c("omcA", "omcS")))),
      "Two metal-reducing cytochromes should read as Probable iron reduction.")
    qa_expect_true(
      grepl("^No detected", fe_red$classifier(marker(fe_red$genes, character()))),
      "Absent cytochromes should read as No detected iron reduction.")
    qa_expect_false(
      any(grepl("^Complete", vapply(list(c("omcA","omcB","omcS","omcZ")),
        function(g) fe_red$classifier(marker(fe_red$genes, g)), character(1)))),
      "Iron reduction must not report a complete pathway from these markers.")

    fe_ox <- registry[["Iron_oxidation"]]
    qa_expect_true(
      grepl("^Probable", fe_ox$classifier(marker(fe_ox$genes, "cyc2"))),
      "cyc2 should read as Probable iron oxidation.")

    # Chlorite dismutase is the diagnostic hallmark of (per)chlorate reduction.
    perc <- registry[["Perchlorate_reduction"]]
    qa_expect_true(
      grepl("^Complete", perc$classifier(marker(perc$genes, c("pcrA", "cld")))),
      "pcrA plus cld should read as Complete (per)chlorate reduction.")

    # arrA (respiratory) and aioA (oxidative) must resolve to distinct labels.
    ar_resp <- registry[["Arsenate_respiration"]]$classifier(
      marker(registry[["Arsenate_respiration"]]$genes, "arrA"))
    ai_ox <- registry[["Arsenite_oxidation"]]$classifier(
      marker(registry[["Arsenite_oxidation"]]$genes, "aioA"))
    qa_expect_true(grepl("arsenate respiration", ar_resp) && grepl("arsenite oxidation", ai_ox),
      "Arsenate respiration and arsenite oxidation classifiers were not distinct.")
    TRUE
  }
)

qa_register_test(
  "BIO_004", "regression", "critical",
  "No functional classifier leaks a negative/insufficient call through its pathway_selector",
  function() {
    # Regression guard for the selector-leak bug: the methanogenesis/WLP
    # family previously counted "Insufficient evidence…" and "No … assigned"
    # outputs as positives in the functional-abundance selection. The unified
    # aaa_positive_pathway_call() must exclude every negative vocabulary.
    qa_expect_true(
      all(!vapply(c(
        "No detected acetogenic potential", "Non-methanogen",
        "Insufficient evidence to classify metabolic potential",
        "No homoacetogenesis assigned: methanogenic pathway detected",
        "Unknown methanotrophic potential", "Excluded: methanogenic pathway detected"
      ), aaa_positive_pathway_call, logical(1))),
      "aaa_positive_pathway_call() let a negative/insufficient/unknown call count as positive."
    )
    qa_expect_true(
      all(vapply(c(
        "Complete perchlorate reduction potential", "Probable iron-oxidation potential",
        "Ectoine producer", "Hydrogenotrophic methanogenesis potential",
        "Partial ectoine biosynthetic potential"
      ), aaa_positive_pathway_call, logical(1))),
      "aaa_positive_pathway_call() dropped a genuine positive call."
    )

    # Fuzz every function: no negative-sounding classifier output may pass its
    # own pathway_selector.
    set.seed(4)
    registry <- aaa_registry()
    negative_rx <- "insufficient|unknown|^no |^non-|assigned|excluded|not detected"
    apply_sel <- function(sel, s) {
      if (is.function(sel)) isTRUE(as.logical(sel(s)))
      else if (is.character(sel)) s %in% sel else FALSE
    }
    leaks <- character()
    for (id in names(registry)) {
      genes <- registry[[id]]$genes
      cl <- registry[[id]]$classifier
      sel <- registry[[id]]$pathway_selector
      for (frac in c(0, 0.3, 0.6, 1)) for (r in 1:25) {
        v <- runif(length(genes)) < frac
        v[runif(length(genes)) < 0.15] <- NA
        out <- tryCatch(cl(stats::setNames(as.logical(v), genes)), error = function(e) NA_character_)
        if (!is.na(out) && grepl(negative_rx, out, ignore.case = TRUE) && apply_sel(sel, out)) {
          leaks <- c(leaks, paste0(id, ": ", out))
        }
      }
    }
    qa_expect_true(
      length(unique(leaks)) == 0L,
      paste0("Pathway selectors leak negative calls as positive: ",
             paste(unique(leaks), collapse = " | "))
    )
    TRUE
  }
)
