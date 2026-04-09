#' Model-Implied Genetic Correlation Matrix
#'
#' Fits a common factor model and computes the model-implied genetic
#' correlation matrix R and its sampling covariance V_R.
#'
#' @param LDSCoutput LDSC result (list or JSON string)
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

  if (is.list(LDSCoutput)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = LDSCoutput$S, v = LDSCoutput$V, i_mat = LDSCoutput$I,
      n_vec = LDSCoutput$N, m = LDSCoutput$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- LDSCoutput
  }

  json <- .Call("wrap__rgmodel_rust",
    as.character(covstruc_json), as.character(est_str),
    as.character(model), as.logical(std.lv), as.character(sub_str))
  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::rgmodel error: ", result$error)

  list(
    R = as.matrix(result$R),
    V_R = as.matrix(result$V_R)
  )
}
