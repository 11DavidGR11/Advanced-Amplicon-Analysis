# Canonical null-coalescing operator for the engine. Defined once here (the
# first core file loaded) instead of being repeated in several engine files.
# Semantics match the previous effective engine behaviour (NULL only; a
# zero-length value is returned as-is). Standalone apps (FASTQ, FunctionBuilder)
# and the CacheManager keep their own copy because they run in separate
# processes that do not load the engine core.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Declares data-frame column names used inside dplyr/tidyr non-standard
# evaluation so that static checks (R CMD check, codetools::checkUsage) do not
# flag them as undefined global variables. These are column identifiers, not
# real objects; the declaration has no runtime effect.
utils::globalVariables(c(
  "Abundance", "Adjusted_P", "base", "Component", "Cos2", "Count",
  "Expected", "Genus", "Inverse_Simpson", "log2FC", "MA_X", "Mean",
  "Mean_abundance", "Observed", "Observed_taxa", "original", "P_value",
  "Pathway", "Percent", "Potential", "Predicted", "RDA1", "RDA2",
  "Sample_column", "SD", "Shannon", "Significance", "Simpson", "Tax_level",
  "Taxon", "Taxonomy", "Tested", "Treatment", "Value", "Variable",
  "Variance_explained", "VIP", "Volcano_Y"
))
