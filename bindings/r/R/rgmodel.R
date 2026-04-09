#' Model-Implied Genetic Correlation Matrix
#'
#' Fits a common factor model and computes the model-implied genetic
#' correlation matrix R and its sampling covariance V_R.
#'
#' @param covstruc LDSC result (list or JSON string)
#' @param estimation Estimation method: "DWLS" (default) or "ML"
#' @return A list with components:
#'   \item{R}{Model-implied genetic correlation matrix}
#'   \item{V_R}{Sampling covariance of vech(R)}
#' @export
rgmodel <- function(covstruc, estimation = "DWLS") {
  if (is.list(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S, v = covstruc$V, i_mat = covstruc$I,
      n_vec = covstruc$N, m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

  json <- .Call("wrap__rgmodel_rust", as.character(covstruc_json), as.character(estimation))
  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::rgmodel error: ", result$error)

  list(
    R = as.matrix(result$R),
    V_R = as.matrix(result$V_R)
  )
}
