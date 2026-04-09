#' Parallel Analysis on LDSC Results
#'
#' Determines the number of non-spurious latent factors via Monte Carlo
#' eigenvalue comparison against simulated null distributions.
#'
#' @param S Genetic covariance matrix (from LDSC output)
#' @param V Sampling covariance matrix (from LDSC output)
#' @param r Number of Monte Carlo simulations (default NULL = 500)
#' @param p P-value threshold (ignored in gsemr)
#' @param save.pdf Save plot to PDF (ignored in gsemr)
#' @param diag Use diagonal of V only (ignored in gsemr)
#' @param fa Factor analysis method (ignored in gsemr)
#' @param fm Factor method (ignored in gsemr)
#' @param nfactors Number of factors to extract (ignored in gsemr)
#' @return A list with components:
#'   \item{observed}{Observed eigenvalues (descending)}
#'   \item{simulated_95}{95th percentile of simulated eigenvalues}
#'   \item{n_factors}{Suggested number of factors}
#' @export
paLDSC <- function(S=S, V=V, r=NULL, p=NULL, save.pdf=FALSE, diag=FALSE, fa=FALSE,
                   fm=NULL, nfactors=NULL) {

  # Ignored params
  if (!is.null(p)) {
    message("Note: 'p' is ignored in gsemr -- uses default 95th percentile threshold")
  }
  if (!identical(save.pdf, FALSE)) {
    message("Note: 'save.pdf' is ignored in gsemr -- plotting is not supported in Rust backend")
  }
  if (!identical(diag, FALSE)) {
    message("Note: 'diag' is ignored in gsemr -- not implemented")
  }
  if (!identical(fa, FALSE)) {
    message("Note: 'fa' is ignored in gsemr -- not implemented")
  }
  if (!is.null(fm)) {
    message("Note: 'fm' is ignored in gsemr -- not implemented")
  }
  if (!is.null(nfactors)) {
    message("Note: 'nfactors' is ignored in gsemr -- not implemented")
  }

  # Default r to 500 if NULL
  if (is.null(r)) r <- 500L

  # Convert S and V matrices to JSON
  s_json <- jsonlite::toJSON(as.matrix(S), digits = 15)
  v_json <- jsonlite::toJSON(as.matrix(V), digits = 15)

  json <- .Call("wrap__pa_ldsc_rust", as.character(s_json), as.character(v_json), as.integer(r))
  jsonlite::fromJSON(json)
}
