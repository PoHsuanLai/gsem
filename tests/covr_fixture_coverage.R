#!/usr/bin/env Rscript
#
# R-side coverage of GenomicSEM, DRIVEN BY THE FIXTURE GENERATORS.
#
# Instruments the GenomicSEM source package and runs each generate_*.R as the
# "test code". R lines/branches the generators never execute are options that
# no committed fixture exercises -> the R-side gap, directly comparable to the
# Rust-side cargo-llvm-cov gap.
#
# Usage:  Rscript tests/covr_fixture_coverage.R   (run from repo root OR tests/)
# Output: bench/covr_by_function.txt + printed overall %

suppressMessages({ library(covr) })

here   <- normalizePath(dirname(sub("--file=", "",
            grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
repo   <- normalizePath(file.path(here, ".."))
pkg    <- file.path(repo, "GenomicSEM")
stopifnot(file.exists(file.path(pkg, "DESCRIPTION")))

gens <- file.path(here, c(
  "generate_synthetic_reference.R",
  "generate_gwas_fixture.R",
  "generate_sldsc_reference.R",
  "generate_enrich_reference.R",
  "generate_covstruc_reference.R"
))
gens <- gens[file.exists(gens)]

# Write a standalone driver that covr will source in its subprocess.
driver <- tempfile("covr_driver_", fileext = ".R")
writeLines(c(
  sprintf('setwd(%s)', deparse(here)),
  'for (g in c(',
  paste0('    ', vapply(gens, deparse, ""), collapse = ",\n"),
  ')) {',
  '  message("### covr running ", basename(g))',
  '  tryCatch(source(g, local = new.env()),',
  '           error = function(e) message("  (failed: ", conditionMessage(e), ")"))',
  '}'
), driver)

cov <- covr::package_coverage(path = pkg, type = "none", code = sprintf('source(%s)', deparse(driver)))

cat("\n==================== R-side coverage (covr) ====================\n")
cat("OVERALL:", sprintf("%.1f%%", covr::percent_coverage(cov)), "\n\n")
print(covr::coverage_to_list(cov))

outdir <- file.path(repo, "bench")
dir.create(outdir, showWarnings = FALSE)
saveRDS(cov, file.path(outdir, "covr_fixture_coverage.rds"))
writeLines(capture.output(print(covr::coverage_to_list(cov))),
           file.path(outdir, "covr_by_function.txt"))
cat("\nwrote", file.path(outdir, "covr_by_function.txt"), "\n")
