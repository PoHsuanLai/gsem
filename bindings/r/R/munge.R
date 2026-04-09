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
#' @param parallel Use parallel processing (default FALSE)
#' @param cores Number of cores (default NULL)
#' @param overwrite Overwrite existing files (ignored in gsemr)
#' @return Character vector of output file paths (invisibly)
#' @export
munge <- function(files, hm3, trait.names=NULL, N=NULL, info.filter=0.9, maf.filter=0.01,
                  log.name=NULL, column.names=list(), parallel=FALSE, cores=NULL, overwrite=TRUE) {

  # log.name: handled after computation below
  if (identical(parallel, FALSE)) {
    Sys.setenv(RAYON_NUM_THREADS = "1")
  } else if (!is.null(cores)) {
    Sys.setenv(RAYON_NUM_THREADS = as.character(cores))
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
    n_override,
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
