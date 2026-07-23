# Content hashing -------------------------------------------------------------
aaa_hash_object <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(x, algo = "sha256", serialize = TRUE))
  }
  raw <- serialize(x, NULL, version = 3)
  tf <- tempfile()
  writeBin(raw, tf)
  on.exit(unlink(tf))
  unname(tools::md5sum(tf))
}
