#' Merge GWAS Summary Statistics
#'
#' Reads multiple GWAS files, applies QC filters, aligns alleles,
#' and outputs a merged tab-delimited file for multivariate GWAS.
#'
#' @param files Character vector of GWAS file paths
#' @param ref Path to reference panel file (e.g., w_hm3.snplist)
#' @param trait.names Character vector of trait names
#' @param se.logit Logical vector indicating which traits have logistic SEs
#' @param OLS Logical vector indicating which traits are from OLS regression
#' @param linprob Logical vector indicating which traits are linear probability
#' @param N Numeric vector or list of sample size overrides
#' @param betas Named list of beta column overrides per trait (default NULL = auto-detect)
#' @param info.filter INFO score filter (default 0.6)
#' @param maf.filter MAF filter (default 0.01)
#' @param keep.indel Keep indels (default FALSE)
#' @param parallel Accepted for API compatibility; gsemr's sumstats is currently
#'   single-threaded so this has no effect.
#' @param cores Accepted for API compatibility; see \code{parallel}.
#' @param ambig Keep ambiguous SNPs (default FALSE)
#' @param direct.filter Apply MAF filter directly to GWAS file frequencies (default FALSE)
#' @return Path to the merged output file
#' @export
sumstats <- function(files, ref, trait.names=NULL, se.logit, OLS=NULL, linprob=NULL,
                     N=NULL, betas=NULL, info.filter=0.6, maf.filter=0.01,
                     keep.indel=FALSE, parallel=FALSE, cores=NULL, ambig=FALSE,
                     direct.filter=FALSE, out="merged_sumstats.tsv") {

  if (identical(parallel, TRUE) || (!is.null(cores) && is.numeric(cores) && cores > 1)) {
    message("Note: gsemr sumstats is single-threaded; 'parallel'/'cores' are no-ops")
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

  # Handle betas: convert named list to JSON object
  betas_json <- if (is.null(betas)) {
    "{}"
  } else {
    jsonlite::toJSON(as.list(betas), auto_unbox = TRUE)
  }

  json <- .Call("wrap__sumstats_rust",
    as.character(files), as.character(ref), as.character(trait.names),
    as.double(info.filter), as.double(maf.filter), as.logical(keep.indel),
    as.character(out),
    as.integer(se.logit), as.integer(OLS), as.integer(linprob),
    as.character(n_overrides_json),
    as.logical(ambig),
    as.character(betas_json),
    as.logical(direct.filter)
  )
  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::sumstats error: ", result$error)
  result$path
}
