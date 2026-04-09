#' Run LD Score Regression
#'
#' Estimates genetic covariance matrix (S), sampling covariance matrix (V),
#' and intercept matrix (I) from GWAS summary statistics using LD Score
#' Regression with block jackknife.
#'
#' @param traits Character vector of paths to .sumstats.gz files
#' @param sample.prev Numeric vector of sample prevalences (NA for continuous traits)
#' @param population.prev Numeric vector of population prevalences (NA for continuous)
#' @param ld Path to LD score directory (containing {chr}.l2.ldscore.gz files)
#' @param wld Path to weight LD score directory (defaults to \code{ld})
#' @param trait.names Character vector of trait names (defaults to V1, V2, ...)
#' @param sep_weights Use separate weight LD scores (ignored in gsemr)
#' @param chr Number of chromosomes (default 22)
#' @param n.blocks Number of jackknife blocks (default 200)
#' @param ldsc.log Log file path (ignored in gsemr)
#' @param stand Standardize output (default FALSE)
#' @param select Variable selection method (default FALSE)
#' @param chisq.max Maximum chi-square filter (default NA = auto)
#' @return A list with components:
#'   \item{S}{Genetic covariance matrix (k x k)}
#'   \item{V}{Sampling covariance matrix of S (k* x k*, where k* = k(k+1)/2)}
#'   \item{I}{LDSC intercept matrix (k x k)}
#'   \item{N}{Sample size vector}
#'   \item{m}{Number of SNPs used}
#' @export
ldsc <- function(traits, sample.prev, population.prev, ld, wld,
                 trait.names=NULL, sep_weights=FALSE, chr=22,
                 n.blocks=200, ldsc.log=NULL, stand=FALSE, select=FALSE, chisq.max=NA) {

  # Ignored params
  if (!identical(sep_weights, FALSE)) {
    message("Note: sep_weights is always enabled in gsemr -- weight LD scores are read from the wld directory")
  }
  # ldsc.log: handled after computation below

  if (is.null(trait.names)) {
    trait.names <- paste0("V", seq_along(traits))
  }

  # Convert select to string for Rust
  select_str <- if (is.logical(select) && !select) "FALSE" else as.character(select)

  # Convert chisq.max: NA means auto (pass NaN to Rust)
  chisq_max_val <- if (is.na(chisq.max)) NaN else as.double(chisq.max)

  json <- .Call("wrap__ldsc_rust",
    as.character(traits),
    as.double(sample.prev),
    as.double(population.prev),
    as.character(ld),
    as.character(wld),
    as.integer(n.blocks),
    as.integer(chr),
    chisq_max_val,
    as.logical(stand),
    as.character(select_str)
  )

  result <- jsonlite::fromJSON(json)

  if (!is.null(result$error)) {
    stop("gsemr::ldsc error: ", result$error)
  }

  S <- as.matrix(result$s)
  rownames(S) <- colnames(S) <- trait.names

  V <- as.matrix(result$v)
  I <- as.matrix(result$i_mat)
  rownames(I) <- colnames(I) <- trait.names

  out <- list(
    S = S,
    V = V,
    I = I,
    N = result$n_vec,
    m = result$m
  )

  # Write log file if requested
  if (!is.null(ldsc.log)) {
    sink(ldsc.log)
    cat("gsemr LDSC Results\n")
    cat("==================\n\n")
    cat("Traits:", paste(traits, collapse=", "), "\n")
    cat("N blocks:", n.blocks, "\n")
    cat("M (SNPs used):", result$m, "\n\n")
    cat("Genetic Covariance Matrix (S):\n")
    print(round(S, 4))
    cat("\nIntercept Matrix (I):\n")
    print(round(I, 4))
    cat("\nSample Sizes (N):\n")
    print(result$n_vec)
    sink()
    message("LDSC log written to: ", ldsc.log)
  }

  out
}
