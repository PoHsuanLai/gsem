# Locate the repo-root `tests/fixtures/` directory that holds the R-reference
# JSON fixtures and synthetic inputs the Rust integration tests also use.
#
# These parity tests assert the gsemr (Rust-backed) wrappers reproduce the SAME
# values that R GenomicSEM produced when the fixtures were generated — catching
# binding-layer bugs (argument marshalling, list/matrix round-tripping, default
# drift) that the Rust-core tests can't see.
#
# When run under `R CMD check`, the working directory is the unpacked test dir.
# We walk upward to find the repo's `tests/fixtures/`. If it can't be found
# (e.g. installed from a tarball without the repo tree), the fixtures are
# unavailable and the parity tests skip rather than fail.

fixtures_dir <- function() {
  # Candidate roots: env override, then walk up from cwd.
  env <- Sys.getenv("GSEMR_FIXTURES_DIR", unset = NA)
  if (!is.na(env) && dir.exists(env)) {
    return(normalizePath(env))
  }
  start <- getwd()
  dir <- start
  for (i in 1:8) {
    cand <- file.path(dir, "tests", "fixtures")
    if (dir.exists(cand) && file.exists(file.path(cand, "ldsc_synth.json"))) {
      return(normalizePath(cand))
    }
    parent <- dirname(dir)
    if (parent == dir) break
    dir <- parent
  }
  NA_character_
}

have_fixtures <- function() !is.na(fixtures_dir())

skip_if_no_fixtures <- function() {
  testthat::skip_if_not(
    have_fixtures(),
    "repo-root tests/fixtures/ not found (set GSEMR_FIXTURES_DIR)"
  )
}

load_fixture <- function(name) {
  jsonlite::fromJSON(file.path(fixtures_dir(), paste0(name, ".json")))
}

fixture_path <- function(...) file.path(fixtures_dir(), ...)
