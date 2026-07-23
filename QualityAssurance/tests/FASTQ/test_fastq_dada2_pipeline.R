# FASTQ_001 only greps AAApp/FASTQ/app.R's source text for expected function
# calls; it never actually runs dada2 against real reads. This test replays
# the same sequence of dada2 calls the app's paired-end pipeline uses
# (filterAndTrim -> learnErrors -> dada -> mergePairs -> makeSequenceTable ->
# removeBimeraDenovo -> assignTaxonomy) against synthetic FASTQ files with a
# realistic per-position quality profile, and checks the three known input
# taxa are correctly recovered end to end.

qa_fastq_write_gz <- function(path, ids, seqs, quals) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  for (i in seq_along(seqs)) {
    writeLines(c(paste0("@", ids[i]), seqs[i], "+", quals[i]), con)
  }
}

qa_register_test(
  "FASTQ_002", "regression", "critical",
  "The paired-end dada2 pipeline recovers known taxa from synthetic reads",
  function() {
    if (!requireNamespace("dada2", quietly = TRUE)) return(TRUE)

    set.seed(42)
    work_dir <- tempfile("qa_fastq_")
    raw_dir <- file.path(work_dir, "raw"); dir.create(raw_dir, recursive = TRUE)
    filt_dir <- file.path(work_dir, "filtered"); dir.create(filt_dir)

    bases <- c("A", "C", "G", "T")
    random_seq <- function(n) paste(sample(bases, n, replace = TRUE), collapse = "")
    revcomp <- function(seq) {
      chars <- rev(strsplit(seq, "")[[1]])
      comp <- c(A = "T", C = "G", G = "C", T = "A")
      paste(comp[chars], collapse = "")
    }
    mutate <- function(seq, error_rate = 0.005) {
      chars <- strsplit(seq, "")[[1]]
      hits <- runif(length(chars)) < error_rate
      chars[hits] <- sample(bases, sum(hits), replace = TRUE)
      paste(chars, collapse = "")
    }
    phred_to_char <- function(q) rawToChar(as.raw(q + 33L))
    quality_profile <- function(read_len) {
      position <- seq_len(read_len)
      mean_q <- pmax(15, 38 - (position / read_len) * 20)
      q <- pmin(40, pmax(3, round(mean_q + rnorm(read_len, sd = 2))))
      paste(vapply(q, phred_to_char, character(1)), collapse = "")
    }

    true_seqs <- vapply(1:3, function(i) random_seq(250), character(1))
    species_names <- c("Species_alpha", "Species_beta", "Species_gamma")
    samples <- paste0("Sample_", 1:4)
    mix <- list(c(0.7, 0.2, 0.1), c(0.5, 0.3, 0.2), c(0.2, 0.6, 0.2), c(0.1, 0.1, 0.8))
    reads_per_sample <- 2000

    fwd_files <- character(); rev_files <- character()
    for (s in seq_along(samples)) {
      which_species <- sample(1:3, reads_per_sample, replace = TRUE, prob = mix[[s]])
      full_reads <- vapply(true_seqs[which_species], mutate, character(1))
      fwd_reads <- unname(substr(full_reads, 1, 150))
      rev_reads <- unname(vapply(full_reads, function(x) revcomp(substr(x, 101, 250)), character(1)))
      ids <- paste0("read", seq_len(reads_per_sample))
      ff <- file.path(raw_dir, paste0(samples[s], "_R1.fastq.gz"))
      rf <- file.path(raw_dir, paste0(samples[s], "_R2.fastq.gz"))
      qa_fastq_write_gz(ff, ids, fwd_reads, vapply(nchar(fwd_reads), quality_profile, character(1)))
      qa_fastq_write_gz(rf, ids, rev_reads, vapply(nchar(rev_reads), quality_profile, character(1)))
      fwd_files <- c(fwd_files, ff); rev_files <- c(rev_files, rf)
    }

    training_fasta <- file.path(work_dir, "training.fa.gz")
    con <- gzfile(training_fasta, "wt")
    for (i in 1:3) {
      header <- paste0(
        ">Kingdom_Bacteria;Phylum_Testota;Class_Testia;Order_Testales;",
        "Family_Testaceae;Genus_Testus;", species_names[i]
      )
      writeLines(c(header, true_seqs[i]), con)
    }
    close(con)

    filt_f <- file.path(filt_dir, paste0(samples, "_F.fastq.gz"))
    filt_r <- file.path(filt_dir, paste0(samples, "_R.fastq.gz"))

    filtering <- dada2::filterAndTrim(
      fwd_files, filt_f, rev_files, filt_r,
      truncLen = c(150, 150), maxEE = c(2, 2), minLen = 50,
      truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE
    )
    qa_expect_true(all(filtering[, 2] > 0), "A sample lost all reads during filtering.")

    errF <- dada2::learnErrors(filt_f, multithread = FALSE, verbose = 0)
    errR <- dada2::learnErrors(filt_r, multithread = FALSE, verbose = 0)

    derepF <- dada2::derepFastq(filt_f); names(derepF) <- samples
    derepR <- dada2::derepFastq(filt_r); names(derepR) <- samples
    ddF <- dada2::dada(derepF, err = errF, multithread = FALSE, verbose = 0)
    ddR <- dada2::dada(derepR, err = errR, multithread = FALSE, verbose = 0)

    mergers <- dada2::mergePairs(ddF, derepF, ddR, derepR, minOverlap = 12)
    merge_counts <- vapply(mergers, nrow, integer(1))
    qa_expect_true(all(merge_counts > 0), "A sample produced zero merged read pairs.")

    seqtab <- dada2::makeSequenceTable(mergers)
    seqtab_nochim <- dada2::removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = FALSE)
    qa_expect_true(ncol(seqtab_nochim) > 0L, "No ASVs survived chimera removal.")

    taxa <- dada2::assignTaxonomy(seqtab_nochim, training_fasta, multithread = FALSE, verbose = FALSE)
    qa_expect_true(all(!is.na(taxa[, "Kingdom"])), "Not every recovered ASV received a Kingdom-level assignment.")

    species_hits <- sort(unname(taxa[, "Species"]))
    qa_expect_true(
      identical(species_hits, sort(species_names)),
      "The three known synthetic species were not recovered exactly."
    )
    TRUE
  }
)
