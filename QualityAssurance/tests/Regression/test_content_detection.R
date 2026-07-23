qa_register_test(
  "IMPORT_004",
  "regression",
  "critical",
  "Legacy tab-delimited text with an XLS suffix is imported by content",
  function() {
    fixture <- file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance_legacy_text.xls")
    qa_expect_true(file.exists(fixture), "Legacy text-with-XLS-suffix fixture is missing")

    detected <- aaa_detect_table_format(fixture, original_name = "featureTable.sample.total.relative.xls")
    qa_expect_true(identical(detected, "text"), paste("Expected text content, detected", detected))

    imported <- suppressWarnings(aaa_import_table(
      fixture,
      original_name = "featureTable.sample.total.relative.xls"
    ))
    qa_expect_true(nrow(imported) == 40L, "Legacy text-XLS fixture has an unexpected row count")
    qa_expect_true(ncol(imported) == 13L, "Legacy text-XLS fixture has an unexpected column count")
    qa_expect_true(identical(names(imported)[1L], "Taxonomy"), "Legacy text-XLS fixture lost its first column")
    TRUE
  }
)

qa_register_test(
  "IMPORT_005",
  "regression",
  "high",
  "Content detection overrides a misleading filename extension",
  function() {
    source_file <- file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv")
    temporary_file <- tempfile(fileext = ".xlsx")
    on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
    file.copy(source_file, temporary_file, overwrite = TRUE)

    detected <- aaa_detect_table_format(temporary_file, original_name = "misleading.xlsx")
    qa_expect_true(identical(detected, "text"), paste("Expected text content, detected", detected))
    imported <- suppressWarnings(aaa_import_table(temporary_file, original_name = "misleading.xlsx"))
    qa_expect_true(ncol(imported) == 13L, "Misleading extension prevented text import")
    TRUE
  }
)
