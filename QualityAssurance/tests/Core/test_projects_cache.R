qa_register_test(
  "PLUGIN_001", "functional", "critical", "Plugins are discovered with unique identifiers",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Triple_A.R"))
    triple_a_load(QA_ROOT, verbose = FALSE)
    plugins <- triple_a_list_plugins()
    qa_expect_true(nrow(plugins) >= 7L, "Expected at least seven plugins")
    qa_expect_true(!anyDuplicated(plugins$ID), "Duplicate plugin IDs")
    TRUE
  }
)

qa_register_test(
  "PROJECT_001", "functional", "high", "Project create-open round trip works",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Triple_A.R"))
    triple_a_load(QA_ROOT, verbose = FALSE)
    directory <- file.path(getOption("triplea.testing.root"), "project_roundtrip")
    unlink(directory, recursive = TRUE, force = TRUE)
    on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
    triple_a_create_project(directory, "Beta project")
    reopened <- triple_a_open_project(directory)
    qa_expect_equal(reopened$name, "Beta project")
    TRUE
  }
)

qa_register_test(
  "CACHE_001", "functional", "high", "Content hashes are stable and input-sensitive",
  function() {
    key_1 <- aaa_hash_object(list(analysis = "community_structure", distance = "bray"))
    key_2 <- aaa_hash_object(list(analysis = "community_structure", distance = "bray"))
    key_3 <- aaa_hash_object(list(analysis = "community_structure", distance = "jaccard"))
    qa_expect_true(identical(key_1, key_2) && !identical(key_1, key_3))
    TRUE
  }
)

qa_register_test(
  "CORE_012", "regression", "high",
  "aaa_result_cache_counts sums real hits/misses across nested results instead of collapsing to a binary flag",
  function() {
    nested <- list(
      functional_potential = list(
        methanogenesis = list(metadata = list(cache_status = "HIT", cache_hits = 38L, cache_misses = 4L)),
        homoacetogenesis = list(metadata = list(cache_status = "MISS", cache_hits = 0L, cache_misses = 12L))
      ),
      community_structure = list(metadata = list(note = "no cache fields here"))
    )
    counts <- aaa_result_cache_counts(nested)
    qa_expect_true(
      identical(counts$hits, 38L) && identical(counts$misses, 16L),
      paste0("Expected hits=38/misses=16, got hits=", counts$hits, "/misses=", counts$misses)
    )
    qa_expect_true(
      identical(aaa_result_cache_status(nested), "HIT"),
      "A run with at least one cache hit anywhere should still report the binary status as HIT"
    )
    qa_expect_true(
      identical(aaa_result_cache_counts(list())$hits, 0L),
      "An empty results tree should report zero hits, not an error"
    )
    TRUE
  }
)
