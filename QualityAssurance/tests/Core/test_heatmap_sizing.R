# aaa_plot_heatmap() previously hardcoded width=8 regardless of row-label
# length or column count, so realistic full taxonomic lineages (routinely
# 100+ characters) got clipped: the heatmap body was squeezed to a sliver
# and most sample columns were pushed off-canvas. No test rendered a real
# PNG and checked its dimensions, which is how this went unnoticed.

qa_png_pixel_width <- function(path) {
  raw <- readBin(path, "raw", n = 24)
  bytes <- as.integer(raw[17:20])
  sum(bytes * c(16777216L, 65536L, 256L, 1L))
}

qa_register_test(
  "CORE_008", "regression", "high",
  "aaa_measure_text_width_inches scales with label length",
  function() {
    short <- aaa_measure_text_width_inches("ab", 11)
    long <- aaa_measure_text_width_inches(strrep("Lachnospiraceae_wexlerae_", 4), 11)
    qa_expect_true(long > short * 5, "Text-width measurement did not scale with string length.")
    qa_expect_true(identical(aaa_measure_text_width_inches(character(0), 11), 0), "Empty labels should measure as zero width.")
    TRUE
  }
)

qa_register_test(
  "CORE_013", "regression", "critical",
  "aaa_measure_text_width_inches accepts factor labels without erroring",
  function() {
    # aaa_top_abundance() intentionally stores its Taxon column as a factor
    # (to control plot ordering) and then measures those same labels via
    # aaa_flipped_axis_plot_width() -> aaa_measure_text_width_inches(). A
    # bare nzchar() call on a factor errors with "requires a character
    # vector" instead of coercing, which crashed the real top-abundance
    # analysis (nzchar dispatch through [.factor -> NextMethod("[")) the
    # first time a user ran it end-to-end after this helper was introduced.
    factor_labels <- factor(c("Taxon_A", "Taxon_B", NA))
    result <- tryCatch(aaa_measure_text_width_inches(factor_labels, 10), error = function(e) e)
    qa_expect_true(!inherits(result, "error"), "Factor labels crashed aaa_measure_text_width_inches().")
    qa_expect_true(is.numeric(result) && result > 0, "Factor labels did not produce a positive width.")

    result2 <- tryCatch(aaa_flipped_axis_plot_width(factor_labels), error = function(e) e)
    qa_expect_true(!inherits(result2, "error"), "Factor labels crashed aaa_flipped_axis_plot_width().")
    TRUE
  }
)

qa_register_test(
  "CORE_009", "regression", "high",
  "aaa_plot_heatmap auto-sizes wider for long row labels instead of clipping them",
  function() {
    short_labels <- c("Taxon_A", "Taxon_B", "Taxon_C")
    long_labels <- c(
      paste(rep("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales", 3), collapse = ";"),
      paste(rep("d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales", 3), collapse = ";"),
      paste(rep("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria", 3), collapse = ";")
    )
    make_matrix <- function(labels) {
      matrix(runif(length(labels) * 4, 0, 10), nrow = length(labels),
             dimnames = list(labels, paste0("Sample_", 1:4)))
    }
    make_labels <- function(m) {
      matrix(sprintf("%.1f", m), nrow = nrow(m), dimnames = dimnames(m))
    }

    short_matrix <- make_matrix(short_labels)
    long_matrix <- make_matrix(long_labels)
    short_file <- tempfile(fileext = ".png")
    long_file <- tempfile(fileext = ".png")
    aaa_plot_heatmap(short_matrix, short_file, "Short", labels = make_labels(short_matrix))
    aaa_plot_heatmap(long_matrix, long_file, "Long", labels = make_labels(long_matrix))

    qa_expect_true(file.exists(short_file) && file.exists(long_file), "Heatmap PNG(s) were not created.")

    short_width <- qa_png_pixel_width(short_file)
    long_width <- qa_png_pixel_width(long_file)
    qa_expect_true(
      long_width > short_width * 1.5,
      paste0(
        "A heatmap with much longer row labels did not render meaningfully wider ",
        "(short=", short_width, "px, long=", long_width, "px); long labels are likely being clipped again."
      )
    )
    TRUE
  }
)

qa_register_test(
  "CORE_010", "regression", "high",
  "aaa_plot_heatmap's own labels=NULL default does not crash pheatmap",
  function() {
    m <- matrix(runif(9, 0, 10), nrow = 3, dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3")))
    out_file <- tempfile(fileext = ".png")
    result <- tryCatch(aaa_plot_heatmap(m, out_file, "No labels"), error = function(e) e)
    qa_expect_true(!inherits(result, "error"), "Calling aaa_plot_heatmap() without labels (its own documented default) errored.")
    qa_expect_true(file.exists(out_file), "No heatmap file was produced when labels was left at its default (NULL).")
    TRUE
  }
)

# coord_flip() taxon-name bar/lollipop/loading plots (top abundance, PCA
# loadings, PLS-DA VIP, sPLS-DA signature, pathway contributors) previously
# saved at a fixed width regardless of label length, the same clipping
# failure mode as CORE_009 but for ggplot2 rather than pheatmap output.
qa_register_test(
  "CORE_011", "regression", "high",
  "aaa_flipped_axis_plot_width grows with label length so coord_flip() taxon labels are not clipped",
  function() {
    short_width <- aaa_flipped_axis_plot_width(c("Taxon_A", "Taxon_B"))
    long_width <- aaa_flipped_axis_plot_width(c(
      paste(rep("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales", 3), collapse = ";"),
      paste(rep("d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales", 3), collapse = ";")
    ))
    qa_expect_true(
      long_width - short_width > 5,
      "Flipped-axis plot width did not scale with taxon label length; long lineage labels are likely being clipped again."
    )
    qa_expect_true(
      identical(aaa_flipped_axis_plot_width(character(0)), 5),
      "With no labels the flipped-axis width should fall back to the base width."
    )
    TRUE
  }
)
