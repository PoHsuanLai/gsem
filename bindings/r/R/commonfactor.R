#' Fit Common Factor Model
#'
#' Fits a 1-factor CFA model to the genetic covariance structure,
#' with sandwich-corrected standard errors and fit indices.
#'
#' @param covstruc LDSC result (named list with S, V, I, N, m components)
#' @param estimation Estimation method: "DWLS" (default) or "ML"
#' @return A list with components:
#'   \item{results}{Data frame of parameter estimates (lhs, op, rhs, est, se, z, p)}
#'   \item{modelfit}{Named vector of fit indices (chisq, df, p_chisq, aic, cfi, srmr)}
#' @examples
#' # Synthetic 3-trait covariance structure (normally from `ldsc()`).
#' covstruc <- list(
#'   S = matrix(c(0.60, 0.42, 0.35,
#'                0.42, 0.50, 0.30,
#'                0.35, 0.30, 0.40), 3, 3,
#'              dimnames = list(c("V1", "V2", "V3"), c("V1", "V2", "V3"))),
#'   V = diag(6) * 0.001,
#'   I = diag(3),
#'   N = c(1e5, 1e5, 1e5),
#'   m = 1e6
#' )
#' cf <- commonfactor(covstruc, estimation = "DWLS")
#' cf$results   # parameter estimates
#' cf$modelfit  # chisq, df, p_chisq, AIC, CFI, SRMR
#' @export
commonfactor <- function(covstruc, estimation = "DWLS") {
  result <- .Call("wrap__commonfactor_rust",
    .covstruc_as_list(covstruc), as.character(estimation))

  if (!is.null(result$error)) stop("gsemr::commonfactor error: ", result$error)

  list(
    results = as.data.frame(result$parameters, stringsAsFactors = FALSE),
    modelfit = c(
      chisq = result$chisq, df = result$df, p_chisq = result$p_chisq,
      AIC = result$aic, CFI = result$cfi, SRMR = result$srmr
    )
  )
}
