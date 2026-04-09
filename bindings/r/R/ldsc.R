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
#' @param n.blocks Number of jackknife blocks (default 200)
#' @return A list with components:
#'   \item{S}{Genetic covariance matrix (k x k)}
#'   \item{V}{Sampling covariance matrix of S (k* x k*, where k* = k(k+1)/2)}
#'   \item{I}{LDSC intercept matrix (k x k)}
#'   \item{N}{Sample size vector}
#'   \item{m}{Number of SNPs used}
#' @export
ldsc <- function(traits,
                 sample.prev = rep(NA, length(traits)),
                 population.prev = rep(NA, length(traits)),
                 ld,
                 wld = ld,
                 trait.names = NULL,
                 n.blocks = 200L) {

  if (is.null(trait.names)) {
    trait.names <- paste0("V", seq_along(traits))
  }

  json <- .Call("wrap__ldsc_rust",
    as.character(traits),
    as.double(sample.prev),
    as.double(population.prev),
    as.character(ld),
    as.character(wld),
    as.integer(n.blocks)
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

  list(
    S = S,
    V = V,
    I = I,
    N = result$n_vec,
    m = result$m
  )
}
