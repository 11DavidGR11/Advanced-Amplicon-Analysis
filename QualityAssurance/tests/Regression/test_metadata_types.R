qa_register_test(
  "META_002",
  "regression",
  "high",
  "Ordinary text metadata do not trigger date-conversion errors",
  function() {
    infer_type_for_test <- function(values) {
      non_missing <- values[!is.na(values) & trimws(as.character(values)) != ""]
      if (length(non_missing) == 0L) return("empty")
      text <- as.character(non_missing)
      numeric_values <- suppressWarnings(as.numeric(text))
      if (all(is.finite(numeric_values))) return("numeric")
      logical_values <- tolower(text)
      if (all(logical_values %in% c("true", "false", "yes", "no", "0", "1"))) return("logical")
      date_patterns <- c(
        "^\\d{4}-\\d{2}-\\d{2}$",
        "^\\d{4}/\\d{2}/\\d{2}$",
        "^\\d{2}/\\d{2}/\\d{4}$",
        "^\\d{2}-\\d{2}-\\d{4}$",
        "^\\d{4}-\\d{2}-\\d{2}[ T]\\d{2}:\\d{2}(:\\d{2})?$"
      )
      is_date <- all(vapply(text, function(value) {
        any(vapply(date_patterns, grepl, logical(1), x = value, perl = TRUE))
      }, logical(1)))
      if (is_date) return("date")
      unique_count <- length(unique(text))
      if (unique_count <= max(10L, ceiling(0.25 * length(text)))) return("factor")
      "text"
    }
    qa_expect_true(identical(infer_type_for_test(c("Control", "Treatment A", "Treatment B")), "factor"))
    qa_expect_true(identical(infer_type_for_test(c("2026-01-01", "2026-01-02")), "date"))
  }
)
