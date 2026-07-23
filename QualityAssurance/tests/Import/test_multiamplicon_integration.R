qa_register_test(
  "IMP_004", "regression", "critical",
  "MultiAmplicon aligns equivalent headers and enforces integer counts",
  function() {
    source(file.path(QA_ROOT, "AAApp", "MultiAmplicon", "multiamplicon_core.R"), local = FALSE)

    a <- data.frame(
      Taxon = c("A", "B"),
      Kingdom = c("Bacteria", "Bacteria"),
      S1 = c(1, 2),
      S2 = c(0, 3),
      check.names = FALSE
    )
    b <- data.frame(
      Taxon = c("A", "C"),
      Kingdom = c("Bacteria", "Archaea"),
      S1 = c(4, 5),
      S2 = c(6, 7),
      check.names = FALSE
    )

    combined <- multiamplicon_validate_and_combine(
      list(a, b), c("S1", "S2"), c("a.tsv", "b.tsv"), TRUE
    )
    qa_expect_true(identical(names(combined), names(a)), "Output headers or order changed")
    qa_expect_true(nrow(combined) == 3L, "Identical descriptor rows were not aggregated")
    qa_expect_true(combined$S1[combined$Taxon == "A"] == 5, "Duplicate counts were not summed")

    # Equivalent headers in a different order must be accepted and aligned to
    # the first table before rows are combined.
    reordered <- b[c("S2", "Kingdom", "Taxon", "S1")]
    aligned <- multiamplicon_validate_and_combine(
      list(a, reordered), c("S1", "S2"), c("a.tsv", "reordered.tsv"), FALSE
    )
    qa_expect_true(
      identical(names(aligned), names(a)),
      "Differently ordered headers were not aligned to the reference table"
    )
    qa_expect_true(
      identical(as.numeric(aligned$S1), c(1, 2, 4, 5)) &&
        identical(as.numeric(aligned$S2), c(0, 3, 6, 7)),
      "Automatic header alignment changed the count values"
    )

    # Structural differences remain errors.
    missing_column <- b[, setdiff(names(b), "S2"), drop = FALSE]
    missing_error <- tryCatch({
      multiamplicon_validate_headers(list(a, missing_column), c("a.tsv", "missing.tsv"))
      FALSE
    }, error = function(e) TRUE)
    qa_expect_true(missing_error, "A table with a missing column was accepted")

    extra_column <- b
    extra_column$Unexpected <- 1L
    extra_error <- tryCatch({
      multiamplicon_validate_headers(list(a, extra_column), c("a.tsv", "extra.tsv"))
      FALSE
    }, error = function(e) TRUE)
    qa_expect_true(extra_error, "A table with an additional column was accepted")

    duplicated_header <- b
    names(duplicated_header)[2L] <- names(duplicated_header)[1L]
    duplicate_error <- tryCatch({
      multiamplicon_validate_headers(list(a, duplicated_header), c("a.tsv", "duplicate.tsv"))
      FALSE
    }, error = function(e) TRUE)
    qa_expect_true(duplicate_error, "A table with duplicated headers was accepted")

    bad_count <- b
    bad_count$S1[1] <- 1.5
    count_error <- tryCatch({
      multiamplicon_validate_and_combine(list(a, bad_count), c("S1", "S2"))
      FALSE
    }, error = function(e) TRUE)
    qa_expect_true(count_error, "Fractional abundance values were accepted as counts")

    guessed <- multiamplicon_guess_count_columns(list(a, b), sample_rows = 1L)
    qa_expect_true(all(c("S1", "S2") %in% guessed), "Integer count columns were not inferred from a bounded preview")
    qa_expect_true(!"Taxon" %in% guessed, "A descriptor column was inferred as a count column")

    one_file_error <- tryCatch({
      multiamplicon_validate_headers(list(a), "a.tsv")
      FALSE
    }, error = function(e) TRUE)
    qa_expect_true(one_file_error, "A single input file was accepted for multi-amplicon integration")
  }
)
