#' Munge GWAS Summary Statistics
#'
#' Filters and harmonizes GWAS summary statistics for use with
#' \code{\link{ldsc}}.
#'
#' @param files Character vector of GWAS file paths
#' @param hm3 Path to HapMap3 SNP reference file
#' @param trait.names Character vector of trait names
#' @param N Optional sample size override (single value or vector)
#' @param info.filter INFO score filter threshold (default 0.9)
#' @param maf.filter MAF filter threshold (default 0.01)
#' @param log.name Log file name (ignored in gsemr)
#' @param column.names Named list of column name overrides
#' @param parallel Use parallel processing (ignored in gsemr)
#' @param cores Number of cores (ignored in gsemr)
#' @param overwrite Overwrite existing files (ignored in gsemr)
#' @return Character vector of output file paths (invisibly)
#' @export
munge <- function(files, hm3, trait.names=NULL, N=NULL, info.filter=0.9, maf.filter=0.01,
                  log.name=NULL, column.names=list(), parallel=FALSE, cores=NULL, overwrite=TRUE) {

  # Ignored params
  if (!is.null(log.name)) {
    message("Note: 'log.name' is ignored in gsemr -- logging is handled by the Rust runtime")
  }
  if (!identical(parallel, FALSE)) {
    message("Note: 'parallel' is ignored in gsemr -- Rust uses native parallelism automatically")
  }
  if (!is.null(cores)) {
    message("Note: 'cores' is ignored in gsemr -- Rust uses native parallelism automatically")
  }
  if (!identical(overwrite, TRUE)) {
    message("Note: 'overwrite' is ignored in gsemr -- files are always overwritten")
  }

  if (is.null(trait.names)) {
    trait.names <- paste0("trait", seq_along(files))
  }

  # Convert N to a single double (NA if NULL)
  n_override <- if (is.null(N)) NaN else as.double(N[1])

  # Convert column.names list to JSON string
  column_names_json <- if (length(column.names) > 0) {
    jsonlite::toJSON(column.names, auto_unbox = TRUE)
  } else {
    "{}"
  }

  # Output goes to current directory (matching R GenomicSEM behavior)
  out <- "."
  dir.create(out, showWarnings = FALSE, recursive = TRUE)

  result <- .Call("wrap__munge_rust",
    as.character(files),
    as.character(hm3),
    as.character(trait.names),
    as.double(info.filter),
    as.double(maf.filter),
    as.character(out),
    n_override,
    as.character(column_names_json)
  )

  if (length(result) == 0) {
    warning("munge produced no output files")
  }

  invisible(result)
}
