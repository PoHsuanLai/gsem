#' Parallel Analysis on LDSC Results
#'
#' Determines the number of non-spurious latent factors via Monte Carlo
#' eigenvalue comparison against simulated null distributions.
#'
#' @param covstruc LDSC result (list or JSON string)
#' @param r Number of Monte Carlo simulations (default 500)
#' @return A list with components:
#'   \item{observed}{Observed eigenvalues (descending)}
#'   \item{simulated_95}{95th percentile of simulated eigenvalues}
#'   \item{n_factors}{Suggested number of factors}
#' @export
paLDSC <- function(covstruc, r = 500) {
  if (is.list(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S, v = covstruc$V, i_mat = covstruc$I,
      n_vec = covstruc$N, m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

  json <- .Call("wrap__pa_ldsc_rust", as.character(covstruc_json), as.integer(r))
  jsonlite::fromJSON(json)
}
