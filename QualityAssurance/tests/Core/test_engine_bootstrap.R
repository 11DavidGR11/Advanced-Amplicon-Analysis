qa_register_test(
  "CORE_004", "smoke", "critical",
  "Public API loads and the biological registry validates",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Triple_A.R"), local = globalenv())
    triple_a_load(QA_ROOT, install_missing = FALSE, envir = globalenv(), verbose = FALSE)
    qa_expect_true(exists("aaa_registry", mode="function"), "aaa_registry was not loaded")
    x <- aaa_validate_biological_registry(strict=FALSE)
    qa_expect_true(isTRUE(x$valid), paste(x$issues, collapse="; "))
    TRUE
  }
)
