#' Munge GWAS Summary Statistics
#'
#' Filters and harmonizes GWAS summary statistics for use with
#' \code{\link{ldsc}}.
#'
#' @param files Character vector of GWAS file paths
#' @param hm3 Path to HapMap3 SNP reference file
#' @param trait.names Character vector of trait names
#' @param N Optional sample-size override. Accepts \code{NULL} (infer
#'   N from every file), a scalar (applied uniformly to all traits),
#'   or a length-\code{length(files)} vector with one entry per trait
#'   (use \code{NA} in any slot to infer N from that file). Matches
#'   \code{GenomicSEM::munge}'s per-trait override semantics.
#' @param info.filter INFO score filter threshold (default 0.9)
#' @param maf.filter MAF filter threshold (default 0.01)
#' @param log.name Log file name (ignored in gsemr)
#' @param column.names Named list of column name overrides
#' @param parallel Accepted for API compatibility; gsemr's munge is currently
#'   single-threaded so this has no effect.
#' @param cores Accepted for API compatibility; see \code{parallel}.
#' @param overwrite Overwrite existing files (ignored in gsemr)
#' @return Character vector of output file paths (invisibly)
#' @export
munge <- function(files, hm3, trait.names=NULL, N=NULL, info.filter=0.9, maf.filter=0.01,
                  log.name=NULL, column.names=list(), parallel=FALSE, cores=NULL, overwrite=TRUE) {

  # log.name: handled after computation below
  if (identical(parallel, TRUE) || (!is.null(cores) && is.numeric(cores) && cores > 1)) {
    message("Note: gsemr munge is single-threaded; 'parallel'/'cores' are no-ops")
  }

  if (is.null(trait.names)) {
    trait.names <- paste0("trait", seq_along(files))
  }

  # Build a per-trait N-override vector the Rust side reads slot by
  # slot. `NA_real_` in any slot means "infer N from the file for this
  # trait". Scalar N is broadcast across all traits. `NULL` drops to
  # an all-NA vector (no overrides). This matches the per-trait
  # semantics of R GenomicSEM::munge, which also accepts a vector.
  n_per_trait <- if (is.null(N)) {
    rep(NA_real_, length(files))
  } else if (length(N) == 1L) {
    rep(as.double(N), length(files))
  } else if (length(N) == length(files)) {
    as.double(N)
  } else {
    stop(sprintf("length(N)=%d must be 1 or length(files)=%d",
                 length(N), length(files)))
  }

  # Convert column.names list to a simple `{"key":"value",...}` string the
  # Rust side parses in-place (no jsonlite dependency).
  column_names_json <- if (length(column.names) > 0) {
    entries <- mapply(function(k, v) {
      paste0("\"", k, "\":\"", v, "\"")
    }, names(column.names), column.names, USE.NAMES = FALSE)
    paste0("{", paste(entries, collapse = ","), "}")
  } else {
    "{}"
  }

  # Output goes to current directory (matching R GenomicSEM behavior)
  out <- "."
  dir.create(out, showWarnings = FALSE, recursive = TRUE)

  # Skip existing files if overwrite=FALSE
  if (!identical(overwrite, TRUE)) {
    existing <- sapply(trait.names, function(tn) {
      file.exists(file.path(out, paste0(tn, ".sumstats.gz")))
    })
    if (any(existing)) {
      skip_names <- trait.names[existing]
      message("Skipping existing files (overwrite=FALSE): ", paste(skip_names, collapse=", "))
      files <- files[!existing]
      trait.names <- trait.names[!existing]
      n_per_trait <- n_per_trait[!existing]
      if (length(files) == 0) {
        message("All files already exist, nothing to munge")
        return(invisible(character(0)))
      }
    }
  }

  result <- .Call("wrap__munge_rust",
    as.character(files),
    as.character(hm3),
    as.character(trait.names),
    as.double(info.filter),
    as.double(maf.filter),
    as.character(out),
    as.double(n_per_trait),
    as.character(column_names_json)
  )

  if (length(result) == 0) {
    warning("munge produced no output files")
  }

  # Write log file if requested
  if (!is.null(log.name)) {
    sink(log.name)
    cat("gsemr Munge Results\n")
    cat("===================\n\n")
    cat("Input files:", paste(files, collapse=", "), "\n")
    cat("Trait names:", paste(trait.names, collapse=", "), "\n")
    cat("INFO filter:", info.filter, "\n")
    cat("MAF filter:", maf.filter, "\n")
    cat("Output files:", paste(result, collapse=", "), "\n")
    sink()
    message("Munge log written to: ", log.name)
  }

  invisible(result)
}
