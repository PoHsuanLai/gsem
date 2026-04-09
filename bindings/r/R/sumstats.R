#' Merge GWAS Summary Statistics
#'
#' Reads multiple GWAS files, applies QC filters, aligns alleles,
#' and outputs a merged tab-delimited file for multivariate GWAS.
#'
#' @param files Character vector of GWAS file paths
#' @param ref Path to LD score reference directory
#' @param trait.names Character vector of trait names
#' @param se.logit Logical vector indicating which traits have logistic SEs
#' @param OLS Logical vector indicating which traits are from OLS regression
#' @param linprob Logical vector indicating which traits are linear probability
#' @param N Numeric vector or list of sample size overrides
#' @param betas Named list of beta column overrides (ignored in gsemr)
#' @param info.filter INFO score filter (default 0.6)
#' @param maf.filter MAF filter (default 0.01)
#' @param keep.indel Keep indels (default FALSE)
#' @param parallel Use parallel processing (ignored in gsemr)
#' @param cores Number of cores (ignored in gsemr)
#' @param ambig Keep ambiguous SNPs (ignored in gsemr)
#' @param direct.filter Apply direct filter (ignored in gsemr)
#' @return Path to the merged output file
#' @export
sumstats <- function(files, ref, trait.names=NULL, se.logit, OLS=NULL, linprob=NULL,
                     N=NULL, betas=NULL, info.filter=0.6, maf.filter=0.01,
                     keep.indel=FALSE, parallel=FALSE, cores=NULL, ambig=FALSE, direct.filter=FALSE) {

  # Ignored params
  if (!is.null(betas)) {
    message("Note: 'betas' is ignored in gsemr -- column detection is automatic")
  }
  if (!identical(parallel, FALSE)) {
    message("Note: 'parallel' is ignored in gsemr -- Rust uses native parallelism automatically")
  }
  if (!is.null(cores)) {
    message("Note: 'cores' is ignored in gsemr -- Rust uses native parallelism automatically")
  }
  if (!identical(ambig, FALSE)) {
    message("Note: 'ambig' is ignored in gsemr -- ambiguous SNPs are always removed")
  }
  if (!identical(direct.filter, FALSE)) {
    message("Note: 'direct.filter' is ignored in gsemr -- not implemented")
  }

  if (is.null(trait.names)) {
    trait.names <- tools::file_path_sans_ext(basename(files))
  }

  n_traits <- length(files)

  # Handle se.logit: default to FALSE for all traits if missing
  if (missing(se.logit) || is.null(se.logit)) {
    se.logit <- rep(FALSE, n_traits)
  }

  # Handle OLS: default to NULL -> FALSE for all traits
  if (is.null(OLS)) {
    OLS <- rep(FALSE, n_traits)
  }

  # Handle linprob: default to NULL -> FALSE for all traits
  if (is.null(linprob)) {
    linprob <- rep(FALSE, n_traits)
  }

  # Handle N overrides: convert to JSON
  n_overrides_json <- if (is.null(N)) {
    "[]"
  } else {
    jsonlite::toJSON(as.list(N), auto_unbox = TRUE)
  }

  out <- "merged_sumstats.tsv"

  json <- .Call("wrap__sumstats_rust",
    as.character(files), as.character(ref), as.character(trait.names),
    as.double(info.filter), as.double(maf.filter), as.logical(keep.indel),
    as.character(out),
    as.integer(se.logit), as.integer(OLS), as.integer(linprob),
    as.character(n_overrides_json)
  )
  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::sumstats error: ", result$error)
  result$path
}
