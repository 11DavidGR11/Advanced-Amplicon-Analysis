qa_register_test(
  "RELEASE_001", "release", "critical",
  "Distribution contains no bundled runtime results or cache inside production sources",
  function() {
    production_roots <- file.path(QA_ROOT, c("AAApp", "Plugins", "Resources"))
    production_roots <- production_roots[dir.exists(production_roots)]

    forbidden_dir_names <- c("Results", "Runs", "Cache", "cache", "tmp", "temp")
    forbidden_file_patterns <- c(
      "\\.sqlite$", "\\.sqlite3$", "\\.db$", "\\.rds$", "\\.RData$",
      "\\.Rhistory$", "\\.log$", "validation_report\\.(csv|json|rds|txt)$",
      "plugin_matrix\\.(csv|json|txt)$", "session_info\\.txt$"
    )

    bad_dirs <- character()
    bad_files <- character()
    for (root in production_roots) {
      dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
      bad_dirs <- c(bad_dirs, dirs[basename(dirs) %in% forbidden_dir_names])
      files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
      for (pattern in forbidden_file_patterns) {
        bad_files <- c(bad_files, files[grepl(pattern, basename(files), ignore.case = TRUE)])
      }
    }

    # QA output is generated locally by Run_All_Tests.R and must not be bundled
    # in a release archive. The placeholder is allowed.
    report_dir <- file.path(QA_ROOT, "QualityAssurance", "reports")
    if (dir.exists(report_dir)) {
      report_files <- list.files(report_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
      report_files <- report_files[basename(report_files) != ".gitkeep"]
      bad_files <- c(bad_files, report_files)
    }
    qa_history <- file.path(QA_ROOT, "QualityAssurance", ".Rhistory")
    if (file.exists(qa_history)) bad_files <- c(bad_files, qa_history)

    bad_dirs <- unique(bad_dirs)
    bad_files <- unique(bad_files)
    qa_expect_true(
      length(bad_dirs) == 0L && length(bad_files) == 0L,
      paste(
        c(
          if (length(bad_dirs)) paste("Bundled runtime folders:", paste(bad_dirs, collapse = ", ")),
          if (length(bad_files)) paste("Bundled runtime files:", paste(bad_files, collapse = ", "))
        ),
        collapse = " | "
      )
    )
    TRUE
  }
)

qa_register_test(
  "RELEASE_002", "release", "critical",
  "Project root contains only the main launcher script and essential non-script files",
  function() {
    root_r <- list.files(QA_ROOT, pattern = "\\.[Rr]$", full.names = FALSE)
    qa_expect_true(
      identical(sort(root_r), "Run_Triple_A.R"),
      paste("Unexpected root-level R scripts:", paste(root_r, collapse = ", "))
    )
    # A single README.md at the project root is the repository landing page.
    # READMEs anywhere else are legacy duplicates of Resources/Documentation/
    # and are still rejected, because two copies of the same guidance drift.
    readmes <- list.files(
      QA_ROOT, pattern = "^README\\.md$", recursive = TRUE,
      full.names = TRUE, ignore.case = TRUE
    )
    root_path <- normalizePath(QA_ROOT, winslash = "/", mustWork = FALSE)
    nested_readmes <- readmes[
      normalizePath(dirname(readmes), winslash = "/", mustWork = FALSE) != root_path
    ]
    qa_expect_true(
      length(nested_readmes) == 0L,
      paste(
        "Legacy README files remain outside the project root:",
        paste(nested_readmes, collapse = ", ")
      )
    )

    # Root-level markdown files are almost always maintainer working notes
    # (audit logs, changelog drafts); catch new ones here instead of relying on
    # a manual audit. README.md is the one exception: it is the repository
    # landing page and cannot reach the user, because Build_Clean_Release.R
    # copies an explicit allow-list rather than pruning a full copy.
    allowed_root_md <- "README.md"
    root_md <- setdiff(
      list.files(QA_ROOT, pattern = "\\.md$", full.names = FALSE),
      allowed_root_md
    )
    qa_expect_true(
      length(root_md) == 0L,
      paste(
        "Unexpected root-level markdown files (add to Build_Clean_Release.R",
        "exclusions if they are maintainer notes, or move end-user docs",
        "under Resources/Documentation/):", paste(root_md, collapse = ", ")
      )
    )
    TRUE
  }
)

