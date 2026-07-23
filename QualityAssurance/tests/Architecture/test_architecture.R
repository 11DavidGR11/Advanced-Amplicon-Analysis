qa_test("ARCH_100", "smoke", "critical", "Single-entry architecture and central storage paths are coherent", {
  required <- c(
    "Run_Triple_A.R",
    "AAApp/Launcher/app.R",
    "AAApp/Biological/app.R",
    "AAApp/FASTQ/app.R",
    "AAApp/Common/paths.R",
    "AAApp/Common/Engine/Triple_A.R",
    "Resources", "Resources/Documentation", "Resources/FunctionalDB",
    "Cache", "Cache/GenomeCache.sqlite", "Cache/GFF",
    "Results", "Plugins", "Tools/Diagnostics/Verify_Installation.R", "Tools/Release/Build_Clean_Release.R"
  )
  qa_expect_files(file.path(QA_ROOT, required))
  qa_expect_false(file.exists(file.path(QA_ROOT, "AAApp", "Run_Triple_A.R")), "AAApp/Run_Triple_A.R must not exist")
  qa_expect_false(file.exists(file.path(QA_ROOT, "AAApp", "Verify_Installation.R")), "Diagnostics must live under Tools/Diagnostics")
  qa_expect_false(file.exists(file.path(QA_ROOT, "Run_MultiAmplicon.R")), "Only Run_Triple_A.R may remain at project root")
  qa_expect_false(file.exists(file.path(QA_ROOT, "Build_Clean_Release.R")), "Release tools must live under Tools/Release")
  qa_expect_false(dir.exists(file.path(QA_ROOT, "Developer")), "Legacy Developer directory must not exist")
  qa_expect_false(dir.exists(file.path(QA_ROOT, "Database")), "Legacy Database directory must not exist")
  qa_expect_false(dir.exists(file.path(QA_ROOT, "Documentation")), "Documentation must live under Resources")
  qa_expect_false(dir.exists(file.path(QA_ROOT, "Cache", "Proteins")), "Protein cache must not exist without a protein-download stage")
  root_entry <- paste(readLines(file.path(QA_ROOT, "Run_Triple_A.R"), warn=FALSE), collapse="\n")
  qa_expect_true(grepl('AAApp", "Launcher', root_entry, fixed=TRUE), "Root entry must launch AAApp/Launcher")
})
