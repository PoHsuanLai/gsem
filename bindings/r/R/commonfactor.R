#' Fit Common Factor Model
#'
#' Fits a 1-factor CFA model to the genetic covariance structure,
#' with sandwich-corrected standard errors and fit indices.
#'
#' @param covstruc LDSC result (list with S, V, I components, or JSON string)
#' @param estimation Estimation method: "DWLS" (default) or "ML"
#' @return A list with components:
#'   \item{results}{Data frame of parameter estimates (lhs, op, rhs, est, se, z, p)}
#'   \item{modelfit}{Named vector of fit indices (chisq, df, p_chisq, aic, cfi, srmr)}
#' @export
commonfactor <- function(covstruc, estimation = "DWLS") {
  if (is.list(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S,
      v = covstruc$V,
      i_mat = covstruc$I,
      n_vec = covstruc$N,
      m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

  json <- .Call("wrap__commonfactor_rust", as.character(covstruc_json), as.character(estimation))
  result <- jsonlite::fromJSON(json)

  if (!is.null(result$error)) stop("gsemr::commonfactor error: ", result$error)

  list(
    results = as.data.frame(result$parameters),
    modelfit = c(
      chisq = result$chisq, df = result$df, p_chisq = result$p_chisq,
      AIC = result$aic, CFI = result$cfi, SRMR = result$srmr
    )
  )
}
