qa_test("UI_104", "regression", "critical", "Each application discovers contextual documentation from its own folder", {
  help_module <- file.path(QA_ROOT, "AAApp", "Common", "help.R")
  qa_expect_true(file.exists(help_module), "Shared contextual-help module is missing.")
  text <- paste(readLines(help_module, warn=FALSE), collapse="\n")
  qa_expect_true(grepl("recursive = TRUE", text, fixed=TRUE), "Contextual help is not recursively discovered.")
  components <- c("Biological","MultiAmplicon","FASTQ","CacheManager","FunctionBuilder","General")
  for (component in components) {
    folder <- file.path(QA_ROOT,"Resources","Documentation",component)
    qa_expect_true(dir.exists(folder), paste("Documentation folder missing:",component))
    guides <- list.files(folder,pattern="\\.[Mm][Dd]$",recursive=TRUE,full.names=TRUE)
    qa_expect_true(length(guides)>0L, paste("No documentation installed for:",component))
    qa_expect_true(all(file.info(guides)$size>0), paste("Empty documentation in:",component))
  }
  launcher <- paste(readLines(file.path(QA_ROOT,"AAApp","Launcher","app.R"),warn=FALSE),collapse="\n")
  qa_expect_false(grepl("open_help_center",launcher,fixed=TRUE),"Launcher still exposes a standalone Help Center.")
  # The Launcher itself should surface the onboarding docs (What is Triple_A,
  # installation, workflow) via the same shared contextual-help mechanism the
  # 5 sub-apps use, not a bespoke Help Center.
  qa_expect_true(grepl("aaa_help_button()", launcher, fixed=TRUE), "Launcher does not expose a help button.")
  qa_expect_true(grepl('aaa_register_context_help(input, output, session, root, "General")', launcher, fixed=TRUE),
                 "Launcher does not register contextual help for the General documentation folder.")
})
