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

  x_mat <- as.matrix(PREDICTORS)
  y_vec <- as.numeric(Y)
  v_mat <- as.matrix(V_Y)

  x_json <- jsonlite::toJSON(x_mat, digits = 15)
  v_json <- jsonlite::toJSON(v_mat, digits = 15)

  json <- .Call("wrap__summary_gls_rust",
    as.character(x_json),
    as.numeric(y_vec),
    as.character(v_json),
    as.logical(INTERCEPT)
  )

  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::summaryGLS error: ", result$error)
  as.data.frame(result)
}
