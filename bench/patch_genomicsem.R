#' Hot-patch GenomicSEM internals in the current R session.
#'
#' GenomicSEM v0.0.5 (commit 123ec96) has an unguarded `solve()` call in
#' `.commonfactorGWAS_main` that crashes on SNPs where lavaan converges to
#' a degenerate local minimum (see upstream issue TBD). Until the fix is
#' merged upstream, this helper swaps in a patched version of the internal
#' at runtime so the bench can time and compare R's commonfactorGWAS on
#' real PGC data without hitting the crash.
#'
#' The patched source lives at
#'   /Users/pohsuanlai/Documents/GenomicSEM-rs/GenomicSEM/R/commonfactorGWAS_main.R
#' on the `fix-commonfactorGWAS-unguarded-solve` branch of the upstream
#' clone. Applying the patch in-memory (rather than rebuilding the R
#' package) keeps the bench lightweight and doesn't pollute the user's
#' installed GenomicSEM.
#'
#' @return Invisibly TRUE if the patch was applied, FALSE if the patched
#'   file could not be found (in which case the caller should fall back
#'   to treating R commonfactorGWAS as broken on this input).
patch_genomicsem <- function(patched_file = NULL) {
  if (is.null(patched_file)) {
    # Default: the upstream clone on the patch branch
    here <- tryCatch({
      args <- commandArgs(trailingOnly = FALSE)
      dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
    }, error = function(e) ".")
    if (is.na(here) || !nzchar(here)) here <- "."
    candidates <- c(
      file.path(here, "..", "GenomicSEM", "R", "commonfactorGWAS_main.R"),
      file.path(here, "GenomicSEM", "R", "commonfactorGWAS_main.R"),
      "GenomicSEM/R/commonfactorGWAS_main.R",
      "../GenomicSEM/R/commonfactorGWAS_main.R"
    )
    patched_file <- candidates[file.exists(candidates)][1]
    if (is.na(patched_file)) {
      warning("patch_genomicsem: cannot locate patched commonfactorGWAS_main.R; ",
              "tried: ", paste(candidates, collapse = ", "))
      return(invisible(FALSE))
    }
  }

  # Verify the file contains the fix (look for the `bread_catch` guard that
  # only exists in the patched version)
  patched_src <- readLines(patched_file, warn = FALSE)
  if (!any(grepl("bread_catch", patched_src))) {
    warning(sprintf(
      "patch_genomicsem: %s does not contain the bread_catch guard; ",
      patched_file))
    warning("  is this really the patched branch?")
    return(invisible(FALSE))
  }

  # Load the patched closure into a private env and install it into the
  # GenomicSEM namespace, unlocking the binding first so assign() works.
  env <- new.env()
  source(patched_file, local = env)
  if (!exists(".commonfactorGWAS_main", envir = env)) {
    warning(sprintf(
      "patch_genomicsem: %s did not define .commonfactorGWAS_main",
      patched_file))
    return(invisible(FALSE))
  }
  patched_fn <- env$.commonfactorGWAS_main
  environment(patched_fn) <- asNamespace("GenomicSEM")

  ns <- asNamespace("GenomicSEM")
  tryCatch({
    unlockBinding(".commonfactorGWAS_main", ns)
    assign(".commonfactorGWAS_main", patched_fn, envir = ns)
    lockBinding(".commonfactorGWAS_main", ns)
  }, error = function(e) {
    stop("patch_genomicsem: failed to install patched binding: ",
         conditionMessage(e))
  })

  message(sprintf(
    "patch_genomicsem: installed patched .commonfactorGWAS_main from %s",
    normalizePath(patched_file, mustWork = FALSE)))
  invisible(TRUE)
}
