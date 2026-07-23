qa_register_test("DATA_001", "regression", "critical",
  "Canonical dataset is created without legacy OTU columns", function() {
    d <- utils::read.csv(file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv"), check.names = FALSE)
    samples <- grep("^Sample_", names(d), value = TRUE)
    design <- data.frame(
      Sample_column = samples,
      Treatment = rep(c("A", "B", "C"), each = 4),
      Replicate = rep(1:4, 3), stringsAsFactors = FALSE
    )
    x <- aaa_new_dataset(d, design)
    qa_expect_true(inherits(x, "Triple_A_dataset"))
    qa_expect_true("FeatureID" %in% names(x$abundance))
    qa_expect_true(!"#OTU_num" %in% names(x$abundance))
  })

qa_register_test("DATA_002", "regression", "critical",
  "Scientific preparation consumes only the canonical dataset", function() {
    d <- utils::read.csv(file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid", "minimal_abundance.csv"), check.names = FALSE)
    samples <- grep("^Sample_", names(d), value = TRUE)
    design <- data.frame(Sample_column = samples,
      Treatment = rep(c("A", "B", "C"), each = 4),
      Replicate = rep(1:4, 3), stringsAsFactors = FALSE)
    x <- aaa_new_dataset(d, design)
    out <- aaa_prepare_amplicon_data(x, abundance_type = "counts",
      project_dir = tempfile("triplea_data_test_"), analysis_name = "prepare", filter_genus = FALSE)
    qa_expect_true(length(out$sample_columns) == 12L)
    qa_expect_true(nrow(out$long) > 0L)
  })
