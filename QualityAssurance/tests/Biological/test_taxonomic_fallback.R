qa_register_test(
  "FUNC_007",
  "regression",
  "critical",
  "Functional reference resolution uses species first, genus second and never higher ranks",
  function() {
    species_candidates <- aaa_functional_taxon_candidates(
      taxonomy = "Methanosarcina barkeri",
      genus = "Methanosarcina",
      tax_level = "Species"
    )
    qa_expect_equal(species_candidates$Query, c("Methanosarcina barkeri", "Methanosarcina"))
    qa_expect_equal(species_candidates$Rank, c("Species", "Genus"))

    genus_candidates <- aaa_functional_taxon_candidates(
      taxonomy = "Methanosarcina",
      genus = "Methanosarcina",
      tax_level = "Genus"
    )
    qa_expect_equal(genus_candidates$Query, "Methanosarcina")
    qa_expect_equal(genus_candidates$Rank, "Genus")

    ambiguous_species <- aaa_functional_taxon_candidates(
      taxonomy = "Methanosarcina sp.",
      genus = "Methanosarcina",
      tax_level = "Species"
    )
    qa_expect_equal(ambiguous_species$Query, "Methanosarcina")
    qa_expect_equal(ambiguous_species$Rank, "Genus")

    family_only <- aaa_functional_taxon_candidates(
      taxonomy = "Methanosarcinaceae",
      genus = NA_character_,
      tax_level = "Family"
    )
    qa_expect_true(nrow(family_only) == 0L)
    TRUE
  }
)
