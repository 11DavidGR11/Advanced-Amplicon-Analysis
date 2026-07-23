qa_register_test(
  "FUNC_005",
  "regression",
  "critical",
  "Incomplete functional-reference cache rows are retried and updated",
  function() {
    cache <- data.frame(
      TaxID = NA_character_,
      Taxonomy = "k__Bacteria;p__Example;g__Example",
      Genus = "Example",
      Tax_level = "Genus",
      Reference_genome = NA_character_,
      mcrA = NA,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    qa_expect_true(!aaa_reference_cache_record_is_reusable(cache, 1L))

    updated <- aaa_upsert_reference_cache_record(
      cache = cache,
      taxonomy = cache$Taxonomy[1],
      genus = cache$Genus[1],
      tax_level = cache$Tax_level[1],
      taxid = "12345",
      reference_genome = "GCF_000000001.1"
    )

    qa_expect_true(nrow(updated) == 1L)
    qa_expect_true(aaa_reference_cache_record_is_reusable(updated, 1L))
    qa_expect_true("mcrA" %in% names(updated))
    TRUE
  }
)
