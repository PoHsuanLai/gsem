#' High-Definition Likelihood Estimation
#'
#' Estimates genetic covariance using the HDL method, which uses
#' LD eigenvalue decomposition for more precise estimation than LDSC.
#'
#' @param traits Character vector of paths to .sumstats.gz files
#' @param sample.prev Numeric vector of sample prevalences (NA for continuous traits)
#' @param population.prev Numeric vector of population prevalences (NA for continuous)
#' @param trait.names Character vector of trait names (defaults to V1, V2, ...)
#' @param LD.path Path to HDL LD reference panel directory (text format)
#' @param Nref Reference panel sample size (default 335265 for UKB)
#' @param method HDL method: "piecewise" (default) or "jackknife"
#' @return A list with components:
#'   \item{S}{Genetic covariance matrix (k x k)}
#'   \item{V}{Sampling covariance matrix}
#'   \item{I}{Intercept matrix}
#'   \item{m}{Number of SNPs used}
#' @export
hdl <- function(traits, sample.prev=NA, population.prev=NA, trait.names=NULL,
                LD.path, Nref=335265, method="piecewise") {

  if (is.null(trait.names)) {
    trait.names <- paste0("V", seq_along(traits))
  }

  json <- .Call("wrap__hdl_rust",
    as.character(traits),
    as.double(sample.prev),
    as.double(population.prev),
    as.character(LD.path),
    as.double(Nref),
    as.character(method)
  )

  result <- jsonlite::fromJSON(json)

  if (!is.null(result$error)) {
    stop("gsemr::hdl error: ", result$error)
  }

  S <- as.matrix(result$s)
  rownames(S) <- colnames(S) <- trait.names

  V <- as.matrix(result$v)
  # V is the sampling covariance matrix of the vectorized S elements;
  # label with trait-pair names for interpretability
  k <- length(trait.names)
  vpairs <- c()
  for (i in seq_len(k)) {
    for (j in i:k) {
      vpairs <- c(vpairs, paste0(trait.names[i], "_", trait.names[j]))
    }
  }
  if (nrow(V) == length(vpairs)) {
    rownames(V) <- colnames(V) <- vpairs
  }

  I <- as.matrix(result$i_mat)
  rownames(I) <- colnames(I) <- trait.names

  list(
    S = S,
    V = V,
    I = I,
    m = result$m
  )
}
