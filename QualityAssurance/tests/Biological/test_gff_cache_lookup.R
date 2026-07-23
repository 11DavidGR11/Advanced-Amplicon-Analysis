# Regression coverage for aaa_find_cached_gff(): a double-escaped regex
# (\\\\. instead of \\.) previously made this always return zero matches,
# so aaa_get_gff() treated every already-cached GFF as missing and
# re-downloaded it from NCBI on every run. No test exercised the function
# with a real cached file on disk, which is how the bug went unnoticed.
qa_register_test("FUNC_008", "regression", "critical",
  "A previously cached GFF file is found and reused without re-downloading",
  function() {
    gff_dir <- tempfile("qa_gff_cache_")
    dir.create(gff_dir, recursive = TRUE)

    accession <- "GCF_000005825.2"
    cached_name <- paste0(accession, "_ASM582v2_genomic.gff.gz")
    file.create(file.path(gff_dir, cached_name))

    found <- aaa_find_cached_gff(accession, gff_dir)
    qa_expect_true(
      length(found) == 1L && identical(basename(found), cached_name),
      "An already-cached GFF file was not found on disk."
    )

    # A different accession that is a string prefix of another must not
    # collide (e.g. GCF_000005825.2 vs GCF_000005825.20).
    file.create(file.path(gff_dir, paste0(accession, "0_OTHER_genomic.gff.gz")))
    found_again <- aaa_find_cached_gff(accession, gff_dir)
    qa_expect_true(
      length(found_again) == 1L && identical(basename(found_again), cached_name),
      "A cached-GFF lookup matched an unrelated accession sharing a prefix."
    )

    qa_expect_true(
      length(aaa_find_cached_gff("GCF_999999999.9", gff_dir)) == 0L,
      "An accession with no cached file incorrectly reported a match."
    )
    TRUE
  })
