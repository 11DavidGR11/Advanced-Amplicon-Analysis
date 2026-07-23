qa_register_test("FUNC_006", "regression", "critical",
  "Functional lookup uses a specific taxon and creates the GFF cache lazily",
  function() {
    core <- file.path(QA_ROOT, "AAApp", "Common", "Engine", "Core", "aaa_functional_potential.R")
    src <- paste(readLines(core, warn = FALSE), collapse = "\n")
    qa_expect_true(grepl("aaa_functional_taxon_query", src, fixed = TRUE))
    qa_expect_true(grepl("dir.create(project$gff", src, fixed = TRUE))
    qa_expect_equal(aaa_functional_taxon_query("k__Bacteria;p__Firmicutes;g__Methanosarcina", "Methanosarcina"), "Methanosarcina")
    TRUE
  })

qa_register_test("FUNC_009", "regression", "high",
  "An optional NCBI API key shortens the request delay without changing default (no-key) throttling",
  function() {
    old_option <- getOption("triple_a_ncbi_api_key")
    old_env <- Sys.getenv("NCBI_API_KEY", unset = NA)
    on.exit({
      options(triple_a_ncbi_api_key = old_option)
      if (is.na(old_env)) Sys.unsetenv("NCBI_API_KEY") else Sys.setenv(NCBI_API_KEY = old_env)
    }, add = TRUE)

    options(triple_a_ncbi_api_key = NULL)
    Sys.unsetenv("NCBI_API_KEY")
    qa_expect_equal(aaa_ncbi_request_delay(), 0.35)

    options(triple_a_ncbi_api_key = "dummy_test_key")
    qa_expect_true(
      aaa_ncbi_request_delay() < 0.35,
      "Configuring an NCBI API key did not shorten the inter-request delay"
    )
    TRUE
  })
