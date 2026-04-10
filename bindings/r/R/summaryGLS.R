#' Generalized Least Squares
#'
#' Runs GLS regression on genetic parameters.
#' Port of R GenomicSEM's \code{summaryGLS()}.
#'
#' @param OBJECT LDSC result object (used to extract Y and V_Y if provided)
#' @param Y Response vector (vech of genetic covariance matrix)
#' @param V_Y Covariance matrix of Y
#' @param PREDICTORS Predictor matrix (n x p)
#' @param INTERCEPT Include intercept (default TRUE)
#' @return A data frame with beta, se, z, and p columns
#' @export
summaryGLS <- function(OBJECT = NULL, Y = NULL, V_Y = NULL, PREDICTORS, INTERCEPT = TRUE) {
  # If OBJECT is provided, extract Y and V_Y from it
  if (!is.null(OBJECT)) {
    if (is.null(Y) && !is.null(OBJECT$subS)) {
      Y <- OBJECT$subS
    }
    if (is.null(V_Y) && !is.null(OBJECT$subV)) {
      V_Y <- OBJECT$subV
    }
  }

  if (is.null(Y)) stop("Y (response vector) must be provided")
  if (is.null(V_Y)) stop("V_Y (covariance matrix) must be provided")

  x_mat <- matrix(as.numeric(as.matrix(PREDICTORS)), nrow = nrow(as.matrix(PREDICTORS)))
  y_vec <- as.numeric(Y)
  v_mat <- matrix(as.numeric(as.matrix(V_Y)), nrow = nrow(as.matrix(V_Y)))

  result <- .Call("wrap__summary_gls_rust",
    x_mat,
    as.numeric(y_vec),
    v_mat,
    as.logical(INTERCEPT)
  )

  if (!is.null(result$error)) stop("gsemr::summaryGLS error: ", result$error)
  as.data.frame(result, stringsAsFactors = FALSE)
}
