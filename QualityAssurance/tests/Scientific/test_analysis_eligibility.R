qa_register_test(
  "SCI_014", "regression", "critical",
  "Abundance declarations are checked against the supplied numeric data",
  function() {
    proportion_data <- data.frame(
      Taxonomy = c("g__A", "g__B"),
      S1 = c(0.25, 0.75),
      S2 = c(0.40, 0.60),
      check.names = FALSE
    )
    percentage_data <- transform(
      proportion_data,
      S1 = S1 * 100,
      S2 = S2 * 100
    )
    count_data <- transform(
      proportion_data,
      S1 = c(10, 30),
      S2 = c(20, 30)
    )

    qa_expect_true(
      aaa_validate_abundance_nature(proportion_data, c("S1", "S2"), "proportion")$available,
      "Valid proportions were rejected."
    )
    qa_expect_true(
      aaa_validate_abundance_nature(percentage_data, c("S1", "S2"), "percentage")$available,
      "Valid percentages were rejected."
    )
    qa_expect_true(
      aaa_validate_abundance_nature(count_data, c("S1", "S2"), "counts")$available,
      "Valid counts were rejected."
    )
    qa_expect_true(
      !aaa_validate_abundance_nature(percentage_data, c("S1", "S2"), "proportion")$available,
      "Percentages incorrectly passed as proportions."
    )
    qa_expect_true(
      !aaa_validate_abundance_nature(proportion_data, c("S1", "S2"), "counts")$available,
      "Non-integer proportions incorrectly passed as counts."
    )
  }
)

qa_register_test(
  "SCI_015", "regression", "critical",
  "RDA and differential abundance eligibility checks reject invalid designs",
  function() {
    da_invalid <- aaa_check_differential_abundance_eligibility(
      sample_columns = c("A", "B"),
      treatments = c("Control", "Treatment"),
      replicates = "none"
    )
    da_valid <- aaa_check_differential_abundance_eligibility(
      sample_columns = c("C1", "C2", "T1", "T2"),
      treatments = c("Control", "Treatment"),
      replicates = "duplicate"
    )
    qa_expect_true(!da_invalid$available, "Unreplicated differential abundance was not blocked.")
    qa_expect_true(da_valid$available, "A replicated two-group differential design was rejected.")

    metadata <- data.frame(
      Sample = c("S1", "S2", "S3", "S4"),
      pH = c(6.5, 7.0, 7.5, 8.0),
      Constant = 1,
      stringsAsFactors = FALSE
    )
    rda_valid <- aaa_check_rda_eligibility(
      sample_columns = metadata$Sample,
      metadata = metadata,
      sample_id_column = "Sample",
      environmental_variables = "pH"
    )
    rda_invalid <- aaa_check_rda_eligibility(
      sample_columns = metadata$Sample,
      metadata = metadata,
      sample_id_column = "Sample",
      environmental_variables = "Constant"
    )
    qa_expect_true(rda_valid$available, "A valid RDA design was rejected.")
    qa_expect_true(!rda_invalid$available, "A constant RDA predictor was not blocked.")

    separator_metadata <- data.frame(
      Sample = c("S-1", "S-2", "S-3", "S-4"),
      pH = c(6.5, 7.0, 7.5, 8.0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    separator_rda <- aaa_check_rda_eligibility(
      sample_columns = c("S.1", "S.2", "S.3", "S.4"),
      metadata = separator_metadata,
      sample_id_column = "Sample",
      environmental_variables = "pH"
    )
    qa_expect_true(
      separator_rda$available,
      "Equivalent sample identifiers using hyphens and periods were not aligned."
    )
    qa_expect_true(
      identical(separator_rda$details$normalised_matches, 4L),
      "Normalised RDA sample matches were not reported correctly."
    )

    abundance <- data.frame(
      Taxonomy = c("g__A", "g__B"),
      `S.1` = c(1, 2), `S.2` = c(2, 3),
      `S.3` = c(3, 4), `S.4` = c(4, 5),
      check.names = FALSE
    )
    design <- data.frame(
      Sample_column = c("S.1", "S.2", "S.3", "S.4"),
      Treatment = rep("Group", 4),
      Replicate = seq_len(4),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    roles <- data.frame(
      Column = c("Sample", "pH"),
      Role = c("identifier", "environmental_variable"),
      stringsAsFactors = FALSE
    )
    dataset <- aaa_new_dataset(abundance, design, separator_metadata, roles)
    environment <- aaa_dataset_environment(dataset, "pH")
    qa_expect_true(
      identical(environment$Sample_column, design$Sample_column),
      "Canonical dataset metadata were not aligned to abundance sample identifiers."
    )
  }
)

qa_register_test(
  "SCI_033", "regression", "critical",
  "Both ways of declaring the design agree, including groups with no samples",
  function() {
    # The two design sources produce different types for Treatment: character
    # from consecutive blocks, factor from a metadata column. Nothing downstream
    # may behave differently because of that, so the same experiment declared
    # both ways must give identical numbers.
    set.seed(23)
    groups <- c("CTL", "TRT", "REC")
    per_group <- 4L
    samples <- unlist(lapply(groups, function(g) paste0(g, "_", seq_len(per_group))))
    group_of <- rep(groups, each = per_group)

    taxa <- c(paste0("k__B;p__P;g__Common", 1:6), paste0("k__B;p__P;g__MARKER_", groups))
    counts <- matrix(stats::rpois(length(taxa) * length(samples), 300),
                     nrow = length(taxa), dimnames = list(taxa, samples))
    for (g in groups) {
      row <- paste0("k__B;p__P;g__MARKER_", g)
      counts[row, ] <- 0
      counts[row, group_of == g] <- stats::rpois(per_group, 4000)
    }
    abundance <- data.frame(Taxonomy = taxa, counts, check.names = FALSE)
    metadata <- data.frame(Sample = samples, Grupo = group_of, stringsAsFactors = FALSE)

    block_design <- data.frame(
      Sample_column = samples,
      Treatment = rep(groups, each = per_group),
      Replicate = rep(seq_len(per_group), times = length(groups)),
      stringsAsFactors = FALSE
    )
    by_blocks <- aaa_new_dataset(abundance, block_design)
    by_metadata <- aaa_new_dataset(
      abundance,
      aaa_sample_design_from_metadata(samples, metadata, "Sample", "Grupo")
    )
    qa_expect_true(
      is.character(by_blocks$sample_design$Treatment) &&
        is.factor(by_metadata$sample_design$Treatment),
      "The two design sources no longer produce the types this test is guarding."
    )

    out <- tempfile("triplea_modes_")
    modes <- list(blocks = by_blocks, metadata = by_metadata)

    top <- lapply(modes, function(ds)
      aaa_top_abundance(ds, top_n = 9, abundance_type = "counts", project_dir = out))
    for (name in names(top)) {
      values <- top[[name]]$tables$heatmap_values
      correct <- vapply(groups, function(g) {
        row <- grep(paste0("MARKER_", g), rownames(values))
        if (!length(row)) return(FALSE)
        v <- unlist(values[row, , drop = TRUE])
        v[[g]] > 5 && all(v[setdiff(groups, g)] < 0.01)
      }, logical(1))
      qa_expect_true(all(correct),
        paste0("Group means are wrong when the design comes from ", name, "."))
    }
    normalize <- function(result) {
      m <- result$tables$heatmap_values
      as.matrix(m[order(rownames(m)), order(names(m)), drop = FALSE])
    }
    qa_expect_true(
      isTRUE(all.equal(normalize(top$blocks), normalize(top$metadata))),
      "Group means differ between the two design sources."
    )

    differential <- lapply(modes, function(ds)
      aaa_differential_abundance(ds, "counts", project_dir = out))
    key <- function(result) {
      x <- result$tables$combined
      x <- x[order(x$Comparison, x$Taxonomy), c("Comparison", "Taxonomy", "log2FC", "Adjusted_P")]
      rownames(x) <- NULL
      x
    }
    qa_expect_true(
      isTRUE(all.equal(key(differential$blocks), key(differential$metadata))),
      "Differential abundance differs between the two design sources."
    )

    community <- suppressWarnings(lapply(modes, function(ds)
      aaa_community_analysis(ds, "counts", "hellinger", "bray",
                             permutations = 99, project_dir = out)))
    qa_expect_true(
      isTRUE(all.equal(community$blocks$tables$permanova$R2[1],
                       community$metadata$tables$permanova$R2[1])),
      "PERMANOVA R2 differs between the two design sources."
    )

    plsda <- lapply(modes, function(ds)
      aaa_plsda_analysis(ds, "counts", cv_folds = 3,
                         permutation_repetitions = 0, project_dir = out))
    metric <- function(r, m) r$tables$summary$Value[r$tables$summary$Metric == m]
    qa_expect_true(
      isTRUE(all.equal(metric(plsda$blocks, "R2Y"), metric(plsda$metadata, "R2Y"))) &&
        isTRUE(all.equal(metric(plsda$blocks, "Cross-validated accuracy"),
                         metric(plsda$metadata, "Cross-validated accuracy"))),
      "PLS-DA results differ between the two design sources."
    )

    # A metadata file routinely describes more groups than the samples selected
    # for a given run. The unused group must not be counted as an empty one:
    # that made table() report a zero-size group, so the smallest group size
    # became 0 and the group count was one too many.
    subset_samples <- samples[group_of %in% c("CTL", "TRT")]
    status <- aaa_check_sample_design_eligibility(
      subset_samples, metadata, "Sample", "Grupo"
    )
    qa_expect_true(
      isTRUE(status$available), "A design excluding one metadata group was rejected."
    )
    qa_expect_true(
      identical(length(status$details$groups), 2L),
      paste0("A group with no samples was counted: ", length(status$details$groups), " groups reported.")
    )
    qa_expect_true(
      identical(status$details$minimum_replicates, per_group),
      paste0("Smallest group size is wrong with an unused group: ",
             status$details$minimum_replicates)
    )

    subset_dataset <- aaa_new_dataset(
      abundance[, c("Taxonomy", subset_samples)],
      aaa_sample_design_from_metadata(subset_samples, metadata, "Sample", "Grupo")
    )
    prepared <- aaa_prepare_amplicon_data(
      subset_dataset, "counts", project_dir = out, analysis_name = "subset"
    )
    qa_expect_true(
      identical(prepared$samples_name, c("CTL", "TRT")),
      paste0("An unused group survived into samples_name: ",
             paste(prepared$samples_name, collapse = ", "))
    )
    qa_expect_true(
      identical(prepared$n_replicates, per_group),
      paste0("An unused group corrupted the replicate count: ", prepared$n_replicates)
    )
    TRUE
  }
)

qa_register_test(
  "SCI_032", "regression", "critical",
  "Results depend on the declared design, not on the order of the abundance columns",
  function() {
    # The strongest available check that nothing is still positional: run the
    # same experiment twice, changing only the column order, and require the
    # results to be identical. Marker taxa unique to each group verify
    # independently that group means are computed over the right samples,
    # rather than merely being self-consistent.
    set.seed(11)
    sizes <- c(CTL = 5L, TRT = 3L, REC = 4L)          # deliberately unbalanced
    samples <- unlist(lapply(names(sizes), function(g) paste0(g, "_", seq_len(sizes[[g]]))))
    group_of <- rep(names(sizes), sizes)

    taxa <- c(paste0("k__B;p__P;g__Common", 1:8),
              paste0("k__B;p__P;g__MARKER_", names(sizes)))
    counts <- matrix(stats::rpois(length(taxa) * length(samples), 300),
                     nrow = length(taxa), dimnames = list(taxa, samples))
    for (g in names(sizes)) {
      row <- paste0("k__B;p__P;g__MARKER_", g)
      counts[row, ] <- 0
      counts[row, group_of == g] <- stats::rpois(sizes[[g]], 4000)
    }
    abundance <- data.frame(Taxonomy = taxa, counts, check.names = FALSE)
    metadata <- data.frame(Sample = samples, Grupo = group_of, stringsAsFactors = FALSE)

    build <- function(table) {
      columns <- setdiff(names(table), "Taxonomy")
      design <- aaa_sample_design_from_metadata(columns, metadata, "Sample", "Grupo")
      aaa_new_dataset(table, design)
    }
    shuffled <- abundance[, c("Taxonomy", sample(samples))]
    reference <- build(abundance)
    permuted <- build(shuffled)

    qa_expect_true(
      identical(levels(permuted$sample_design$Treatment), names(sizes)),
      paste0(
        "Group order follows the abundance column order instead of the metadata: ",
        paste(levels(permuted$sample_design$Treatment), collapse = " -> ")
      )
    )

    out <- tempfile("triplea_invariance_")
    top <- lapply(list(reference = reference, permuted = permuted), function(ds)
      aaa_top_abundance(ds, top_n = 11, abundance_type = "counts", project_dir = out))

    for (name in names(top)) {
      values <- top[[name]]$tables$heatmap_values
      correct <- vapply(names(sizes), function(g) {
        row <- grep(paste0("MARKER_", g), rownames(values))
        if (!length(row)) return(FALSE)
        v <- unlist(values[row, , drop = TRUE])
        v[[g]] > 5 && all(v[setdiff(names(sizes), g)] < 0.01)
      }, logical(1))
      qa_expect_true(
        all(correct),
        paste0("Group means are not computed over the samples of each group (", name, ").")
      )
    }

    normalize <- function(result) {
      m <- result$tables$heatmap_values
      as.matrix(m[order(rownames(m)), order(names(m)), drop = FALSE])
    }
    qa_expect_true(
      isTRUE(all.equal(normalize(top$permuted), normalize(top$reference))),
      "Group means change when the abundance columns are reordered."
    )

    differential <- lapply(list(reference = reference, permuted = permuted), function(ds)
      aaa_differential_abundance(ds, "counts", project_dir = out))
    key <- function(result) {
      x <- result$tables$combined
      x <- x[order(x$Comparison, x$Taxonomy), c("Comparison", "Taxonomy", "log2FC", "Adjusted_P")]
      rownames(x) <- NULL
      x
    }
    qa_expect_true(
      isTRUE(all.equal(key(differential$permuted), key(differential$reference))),
      "Differential abundance changes when the abundance columns are reordered."
    )

    # The marker taxa make this synthetic matrix nearly degenerate, so cmdscale
    # reports negative squared distances. That is a property of the fixture, not
    # of the behaviour under test, and it would otherwise mask the result.
    community <- suppressWarnings(
      lapply(list(reference = reference, permuted = permuted), function(ds)
        aaa_community_analysis(ds, "counts", "hellinger", "bray",
                               permutations = 99, project_dir = out))
    )
    qa_expect_true(
      isTRUE(all.equal(community$permuted$tables$permanova$R2[1],
                       community$reference$tables$permanova$R2[1])),
      "PERMANOVA R2 changes when the abundance columns are reordered."
    )
    alpha <- lapply(community, function(result) {
      a <- result$tables$alpha_diversity
      a <- a[order(a$Sample_column), c("Sample_column", "Shannon")]
      rownames(a) <- NULL
      a
    })
    qa_expect_true(
      isTRUE(all.equal(alpha$permuted, alpha$reference)),
      "Per-sample alpha diversity changes when the abundance columns are reordered."
    )
    qa_expect_true(
      identical(
        as.character(community$reference$tables$alpha_diversity$Treatment),
        unname(group_of[match(community$reference$tables$alpha_diversity$Sample_column, samples)])
      ),
      "Alpha-diversity rows carry the wrong treatment label."
    )
    TRUE
  }
)

qa_register_test(
  "SCI_030", "regression", "critical",
  "Groups are assigned by label, so interleaved and unbalanced designs are correct",
  function() {
    # Treatments used to be sliced out of the sample columns in consecutive
    # blocks of equal size. An interleaved layout was therefore mislabelled
    # silently, and unequal group sizes were rejected outright.
    abundance <- utils::read.csv(
      file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv"),
      check.names = FALSE
    )
    samples <- grep("^Sample_", names(abundance), value = TRUE)
    groups <- c("Ctl", "Trt", "Ctl", "Trt", "Ctl", "Trt",
                "Ctl", "Trt", "Trt", "Trt", "Ctl", "Trt")

    design <- data.frame(
      Sample_column = samples,
      Treatment = groups,
      Replicate = stats::ave(seq_along(groups), groups, FUN = seq_along),
      stringsAsFactors = FALSE
    )
    dataset <- aaa_new_dataset(abundance, design)
    prepared <- aaa_prepare_amplicon_data(
      dataset, "counts",
      project_dir = tempfile("triplea_design_"), analysis_name = "design"
    )

    qa_expect_true(
      isTRUE(prepared$has_replicates) && !isTRUE(prepared$balanced),
      "An unbalanced design was not reported as replicated and unbalanced."
    )
    qa_expect_true(
      identical(prepared$n_replicates, 5L),
      paste0("n_replicates should be the smallest group size, got ", prepared$n_replicates)
    )

    values <- as.data.frame(prepared$wide[, prepared$sample_columns, drop = FALSE])
    summary <- aaa_replicate_summary(values, prepared$sample_map, prepared$samples_name)

    for (group in c("Ctl", "Trt")) {
      expected <- rowMeans(values[, samples[groups == group], drop = FALSE])
      qa_expect_true(
        isTRUE(all.equal(unname(summary$mean[[group]]), unname(expected))),
        paste0("Group mean for '", group, "' does not use the samples of that group.")
      )
    }

    # A positional split would have averaged the first six columns together,
    # which for this layout is three Ctl and three Trt samples.
    positional <- rowMeans(values[, samples[1:6], drop = FALSE])
    qa_expect_true(
      !isTRUE(all.equal(unname(summary$mean[["Ctl"]]), unname(positional))),
      "Group means still look like a positional block split."
    )
    TRUE
  }
)

qa_register_test(
  "SCI_031", "regression", "critical",
  "A sample design can be built from a metadata column regardless of row order",
  function() {
    samples <- c("S_1", "S_2", "S_3", "S_4", "S_5")
    metadata <- data.frame(
      # Different separator on purpose: identifier matching normalises it.
      Sample = c("s.4", "s.1", "s.5", "s.2", "s.3"),
      Grupo = c("B", "A", "B", "A", "B"),
      stringsAsFactors = FALSE
    )

    design <- aaa_sample_design_from_metadata(samples, metadata, "Sample", "Grupo")
    qa_expect_true(
      identical(design$Sample_column, samples),
      "The design does not follow the abundance-table sample order."
    )
    qa_expect_true(
      identical(as.character(design$Treatment), c("A", "A", "B", "B", "B")),
      paste0("Groups were not read from the metadata column: ",
             paste(design$Treatment, collapse = ", "))
    )
    # Group order comes from the metadata, which is where the user declares it.
    qa_expect_true(
      identical(levels(design$Treatment), c("B", "A")),
      paste0(
        "Group order does not follow the metadata row order; expected B, A and got: ",
        paste(levels(design$Treatment), collapse = ", ")
      )
    )
    qa_expect_true(
      identical(as.integer(design$Replicate), c(1L, 2L, 1L, 2L, 3L)),
      "Replicate numbers are not ordinals within their own group."
    )

    ok <- aaa_check_sample_design_eligibility(samples, metadata, "Sample", "Grupo")
    qa_expect_true(ok$available, "A valid metadata design was rejected.")
    qa_expect_true(!isTRUE(ok$details$balanced), "A 2/3 design was reported as balanced.")

    incomplete <- metadata
    incomplete$Grupo[2] <- NA
    blocked <- list(
      no_metadata = aaa_check_sample_design_eligibility(samples, NULL, "Sample", "Grupo"),
      same_column = aaa_check_sample_design_eligibility(samples, metadata, "Sample", "Sample"),
      missing_rows = aaa_check_sample_design_eligibility(samples, metadata[1:2, ], "Sample", "Grupo"),
      empty_group = aaa_check_sample_design_eligibility(samples, incomplete, "Sample", "Grupo")
    )
    for (name in names(blocked)) {
      qa_expect_true(
        !isTRUE(blocked[[name]]$available),
        paste0("An invalid metadata design was accepted: ", name)
      )
      qa_expect_true(
        length(blocked[[name]]$guidance) > 0L,
        paste0("A blocked metadata design gives no unlocking guidance: ", name)
      )
    }
    TRUE
  }
)

qa_register_test(
  "SCI_029", "regression", "critical",
  "Every analysis gate blocks invalid designs and explains how to unlock itself",
  function() {
    metadata <- data.frame(
      Sample = c("S1", "S2", "S3", "S4"),
      pH = c(6.5, 7.0, 7.5, 8.0),
      Site = c("a", "a", "b", "b"),
      stringsAsFactors = FALSE
    )
    samples <- metadata$Sample

    # A blocked analysis without instructions is only half a message, so every
    # gate is required to return actionable guidance alongside the reason.
    blocked <- list(
      top_abundance = aaa_check_top_abundance_eligibility(character()),
      community_structure = aaa_check_community_structure_eligibility(c("S1", "S2")),
      community_groups = aaa_check_plsda_eligibility(
        replicates = "duplicate", sample_columns = samples, treatments = "OnlyOne"
      ),
      plsda_unreplicated = aaa_check_plsda_eligibility("none"),
      functional_potential = aaa_check_functional_potential_eligibility(samples, character()),
      taxon_association = aaa_check_taxon_association_eligibility(
        samples, NULL, "Sample", "pH"
      ),
      constrained = aaa_check_constrained_ordination_eligibility(
        c("S1", "S2"), metadata, "Sample", "pH"
      ),
      varpart = aaa_check_variance_partitioning_eligibility(
        samples, metadata, "Sample",
        environmental_variables = "pH", experimental_factors = character()
      )
    )

    for (name in names(blocked)) {
      status <- blocked[[name]]
      qa_expect_true(!isTRUE(status$available), paste0("Invalid design was not blocked: ", name))
      qa_expect_true(
        is.character(status$reason) && nzchar(status$reason),
        paste0("Blocked analysis reports no reason: ", name)
      )
      qa_expect_true(
        length(status$guidance) > 0L && all(nzchar(status$guidance)),
        paste0("Blocked analysis reports no unlocking guidance: ", name)
      )
    }

    # Valid designs must stay unblocked, or the gates would be unusable.
    allowed <- list(
      top_abundance = aaa_check_top_abundance_eligibility(samples),
      community_structure = aaa_check_community_structure_eligibility(samples, c("A", "B")),
      plsda = aaa_check_plsda_eligibility(
        replicates = "duplicate", sample_columns = samples, treatments = c("A", "B")
      ),
      functional_potential = aaa_check_functional_potential_eligibility(samples, "Methanogenesis"),
      taxon_association = aaa_check_taxon_association_eligibility(
        samples, metadata, "Sample", c("pH", "Site")
      ),
      constrained = aaa_check_constrained_ordination_eligibility(
        samples, metadata, "Sample", c("pH", "Site")
      ),
      varpart = aaa_check_variance_partitioning_eligibility(
        samples, metadata, "Sample",
        environmental_variables = "pH", experimental_factors = "Site"
      )
    )

    for (name in names(allowed)) {
      qa_expect_true(
        isTRUE(allowed[[name]]$available),
        paste0("Valid design was incorrectly blocked: ", name)
      )
      qa_expect_true(
        length(allowed[[name]]$guidance) == 0L,
        paste0("An available analysis should carry no unlocking guidance: ", name)
      )
    }
    TRUE
  }
)
