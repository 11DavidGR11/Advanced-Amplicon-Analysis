qa_register_test(
  "IMPORT_001",
  "regression",
  "critical",
  "CSV, TSV, TXT, XLS and XLSX create equivalent canonical abundance tables",
  function() {
    base <- file.path(QA_ROOT, "QualityAssurance", "fixtures", "valid")
    files <- file.path(base, paste0("minimal_abundance.", c("csv", "tsv", "txt", "xls", "xlsx")))
    qa_expect_true(all(file.exists(files)), paste("Missing import fixture(s):", paste(files[!file.exists(files)], collapse = ", ")))

    tables <- lapply(files, aaa_import_table)
    reference <- tables[[1L]]

    for (i in seq_along(tables)) {
      current <- tables[[i]]
      qa_expect_true(nrow(current) == 40L, paste(basename(files[[i]]), "did not import 40 rows"))
      qa_expect_true(ncol(current) == 13L, paste(basename(files[[i]]), "did not import 13 columns"))
      qa_expect_true(identical(names(current), names(reference)), paste(basename(files[[i]]), "has different column names"))
      qa_expect_true(identical(as.character(current$Taxonomy), as.character(reference$Taxonomy)), paste(basename(files[[i]]), "has different taxonomy values"))
      qa_expect_true(isTRUE(all.equal(as.matrix(current[-1]), as.matrix(reference[-1]), check.attributes = FALSE)), paste(basename(files[[i]]), "has different abundance values"))
    }

    TRUE
  }
)

qa_register_test(
  "IMPORT_002",
  "regression",
  "high",
  "The application selectors and messages declare XLS support",
  function() {
    app_file <- file.path(QA_ROOT, "AAApp", "Biological", "app.R")
    importer_file <- file.path(QA_ROOT, "AAApp", "Common", "Engine", "Core", "aaa_importer.R")
    app_text <- qa_read_app_source(QA_ROOT)
    importer_text <- paste(readLines(importer_file, warn = FALSE), collapse = "\n")

    qa_expect_true(grepl('"\\.xls"', app_text), "AAApp file selectors do not include .xls")
    qa_expect_true(grepl('c\\("csv", "tsv", "txt", "xls", "xlsx"\\)', importer_text), "Canonical importer does not declare all supported formats")
    qa_expect_true(grepl('aaa_import_table', app_text, fixed = TRUE), "AAApp does not use the canonical importer")
    TRUE
  }
)
