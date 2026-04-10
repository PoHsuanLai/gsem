#' Model-Implied Genetic Correlation Matrix
#'
#' Fits a common factor model and computes the model-implied genetic
#' correlation matrix R and its sampling covariance V_R.
#'
#' @param LDSCoutput LDSC result (named list with S, V, I, N, m components)
#' @param model lavaan-style model syntax (optional, uses common factor if missing)
#' @param std.lv Standardize latent variables (default TRUE)
#' @param estimation Estimation method: TRUE means "DWLS", "ML" for ML
#' @param sub Subset of output (default NULL = all)
#' @param ... Additional arguments (ignored)
#' @return A list with components:
#'   \item{R}{Model-implied genetic correlation matrix}
#'   \item{V_R}{Sampling covariance of vech(R)}
#' @export
rgmodel <- function(LDSCoutput, model, std.lv=TRUE, estimation=TRUE, sub=NULL, ...) {

  # Convert estimation: TRUE means DWLS, character is used directly
  est_str <- if (is.logical(estimation) && estimation) {
    "DWLS"
  } else if (is.character(estimation)) {
    estimation
  } else {
    "DWLS"
  }

  # Convert sub: NULL means no subset
  sub_str <- if (is.null(sub)) "" else as.character(sub)

  # Handle model: if missing, default to empty string (Rust will use common factor)
  if (missing(model)) model <- ""

  result <- .Call("wrap__rgmodel_rust",
    .covstruc_as_list(LDSCoutput), as.character(est_str),
    as.character(model), as.logical(std.lv), as.character(sub_str))

  if (!is.null(result$error)) stop("gsemr::rgmodel error: ", result$error)

  list(
    R = result$R,
    V_R = result$V_R
  )
}
