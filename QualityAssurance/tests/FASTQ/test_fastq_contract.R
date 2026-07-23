qa_register_test(
  "FASTQ_001", "regression", "critical",
  "FASTQ module supports paired reads and remains preprocessing-only",
  function() {
    x <- paste(readLines(file.path(QA_ROOT,"AAApp","FASTQ","app.R"),warn=FALSE),collapse="\n")
    required <- c("Paired-end","mergePairs","fastp","cutadapt","removeBimeraDenovo","assignTaxonomy","purpose=\"preprocessing_only\"")
    qa_expect_true(all(vapply(required,grepl,logical(1),x=x,fixed=TRUE)),"FASTQ contract incomplete")
    forbidden <- c("PERMANOVA","alpha diversity","DESeq2","functional potential")
    qa_expect_true(!any(vapply(forbidden,grepl,logical(1),x=x,fixed=TRUE)),"FASTQ module contains analytical functionality")
    TRUE
  }
)

qa_register_test(
  "FASTQ_003", "regression", "critical",
  "FASTQ sample names keep trailing digits and pair correctly across naming conventions",
  function() {
    # Regression guard: sample_key() used to strip a SECOND trailing ".1"/".2"
    # after removing the read-direction suffix, eating part of the sample name.
    # Samples like D90.4.1 and D90.4.2 both collapsed to "D90.4", so pair_files()
    # aborted with "sample names are ambiguous" and multi-sample runs could not
    # start. Load the helpers straight from the app source and exercise them.
    src <- readLines(file.path(QA_ROOT, "AAApp", "FASTQ", "app.R"), warn = FALSE)
    s <- grep("^sample_key <- function", src)[1]
    e <- grep("^pair_files <- function", src)[1]
    qa_expect_true(!is.na(s) && !is.na(e), "sample_key/pair_files not found in the FASTQ app.")
    env <- new.env(parent = globalenv())
    eval(parse(text = paste(src[s:(e - 1)], collapse = "\n")), envir = env)
    eval(parse(text = paste(src[e:(e + 5)], collapse = "\n")), envir = env)

    expected <- c(
      "D90.4.1_1.fastq.gz"             = "D90.4.1",
      "D90.4.2_2.fastq.gz"             = "D90.4.2",
      "D90.4.1.raw_1.fastq.gz"         = "D90.4.1.raw",
      "D90.4.1_R1.fastq.gz"            = "D90.4.1",
      "Sample_S1_L001_R1_001.fastq.gz" = "Sample_S1_L001",
      "muestra_1.fastq.gz"             = "muestra",
      "muestra.1.fastq.gz"             = "muestra",
      "muestra-2.fq.gz"                = "muestra"
    )
    for (f in names(expected)) {
      qa_expect_equal(
        env$sample_key(f), expected[[f]],
        paste0("sample_key('", f, "') should be '", expected[[f]], "'.")
      )
    }

    # Sample names ending in 1/2 must stay distinct and pair without error.
    paired <- env$pair_files(
      c("D90.4.1_1.fastq.gz", "D90.4.2_1.fastq.gz", "D90.4.3_1.fastq.gz"),
      c("D90.4.1_2.fastq.gz", "D90.4.2_2.fastq.gz", "D90.4.3_2.fastq.gz")
    )
    qa_expect_equal(
      paired$samples, c("D90.4.1", "D90.4.2", "D90.4.3"),
      "Samples whose names end in 1/2 collapsed together instead of staying distinct."
    )

    # Reverse files supplied in a different order must still be matched.
    shuffled <- env$pair_files(
      c("A.1_R1.fastq.gz", "B.2_R1.fastq.gz"),
      c("B.2_R2.fastq.gz", "A.1_R2.fastq.gz")
    )
    qa_expect_equal(shuffled$reverse_order, c(2L, 1L),
      "Reverse reads were not re-ordered to match their forward partners.")
    TRUE
  }
)
