qa_test(
  "BIO_002", "regression", "critical",
  "Custom biological functions support dictionary aliases, complex rules and a dedicated interface category",
  {
    builder <- file.path(QA_ROOT, "AAApp", "FunctionBuilder", "app.R")
    loader <- file.path(QA_ROOT, "AAApp", "Common", "Engine", "Extensions", "custom_biological_functions.R")
    dictionary_file <- file.path(QA_ROOT, "AAApp", "Common", "Engine", "Extensions", "gene_dictionary.R")
    qa_expect_files(c(builder, loader, dictionary_file))

    temp_root <- file.path(tempdir(), paste0("triple_a_custom_", Sys.getpid()))
    custom_dir <- file.path(temp_root, "Resources", "FunctionalDB", "CustomFunctions")
    unlink(temp_root, recursive = TRUE, force = TRUE)
    dir.create(custom_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)

    alias_updates <- list(geneA = c("QA enzyme alpha", "historical gene A"), geneB = "QA enzyme beta", geneC = "QA regulator")
    aaa_save_custom_gene_aliases(alias_updates, temp_root)
    loaded_aliases <- aaa_load_custom_gene_aliases(temp_root)
    qa_expect_true(all(c("geneA", "geneB", "geneC") %in% names(loaded_aliases)), "Custom aliases were not persisted.")
    qa_expect_true("QA enzyme alpha" %in% loaded_aliases$geneA, "A custom synonym was lost.")

    conflict <- tryCatch({ aaa_merge_gene_aliases(list(geneA = "shared alias"), list(geneB = "shared alias")); FALSE }, error = function(e) TRUE)
    qa_expect_true(conflict, "A synonym assigned to two genes was not rejected.")

    definition <- list(
      id = "QA_custom_function", display_name = "QA custom function", category = "QA biology",
      description = "QA definition", diagnostic = c("geneA", "geneB"), supporting = "geneC", accessory = character(),
      rule = list(type = "expression", minimum = 1L, expression = "(geneA OR geneB) AND geneC"),
      evidence_note = "QA", references = "QA"
    )
    definition_file <- file.path(custom_dir, "QA_custom_function.json")
    jsonlite::write_json(definition, definition_file, auto_unbox = TRUE)

    dictionary <- aaa_merge_gene_aliases(gene_aliases, loaded_aliases)
    normalized <- aaa_validate_custom_function_definition(jsonlite::read_json(definition_file, simplifyVector = FALSE), dictionary, require_dictionary = TRUE)
    qa_expect_true(is.character(normalized$accessory) && length(normalized$accessory) == 0L, "Empty JSON arrays were not normalized.")

    loaded <- suppressWarnings(aaa_load_custom_biological_functions(temp_root, biological_function_registry, dictionary))
    qa_expect_true("QA_custom_function" %in% names(loaded), "Valid custom function was not loaded.")
    qa_expect_true(identical(loaded$QA_custom_function$category, "Custom functions"), "Custom function was not grouped under Custom functions.")
    qa_expect_true(identical(loaded$QA_custom_function$biological_category, "QA biology"), "Biological category was not preserved.")

    positive <- loaded$QA_custom_function$classifier(c(geneA = TRUE, geneB = FALSE, geneC = TRUE))
    partial <- loaded$QA_custom_function$classifier(c(geneA = TRUE, geneB = FALSE, geneC = FALSE))
    negative <- loaded$QA_custom_function$classifier(c(geneA = FALSE, geneB = FALSE, geneC = FALSE))
    qa_expect_true(grepl("potential", positive, fixed = TRUE), "Complex positive rule failed.")
    qa_expect_true(grepl("Partial", partial, fixed = TRUE), "Complex partial evidence was not reported.")
    qa_expect_true(grepl("No detected", negative, fixed = TRUE), "Complex negative rule failed.")

    invalid_gene <- definition
    invalid_gene$diagnostic <- c("missingGene")
    invalid_gene$rule$expression <- "missingGene AND geneC"
    rejected <- tryCatch({ aaa_validate_custom_function_definition(invalid_gene, dictionary, require_dictionary = TRUE); FALSE }, error = function(e) TRUE)
    qa_expect_true(rejected, "A gene absent from the dictionary was accepted.")

    invalid_expression <- definition
    invalid_expression$rule$expression <- "geneA AND (geneB OR)"
    rejected_expression <- tryCatch({ aaa_validate_custom_function_definition(invalid_expression, dictionary, require_dictionary = TRUE); FALSE }, error = function(e) TRUE)
    qa_expect_true(rejected_expression, "Malformed Boolean logic was accepted.")
  }
)
