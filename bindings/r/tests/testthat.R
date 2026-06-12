# Test runner for gsemr (Rust-backed GenomicSEM).
# Run via `R CMD check` or `testthat::test_local()`.
library(testthat)
library(gsemr)

test_check("gsemr")
