#' Generalized Least Squares
#'
#' Runs GLS regression on genetic parameters. Drop-in replacement for
#' \code{GenomicSEM::summaryGLS()}: returns a numeric matrix with
#' columns \code{betas, pvals, SE, Z} and row names \code{b0, b1, ...}
#' (intercept first when \code{INTERCEPT = TRUE}), and prints it before
#' returning.
#'
#' @param OBJECT LDSC result object (used to extract Y and V_Y if provided)
#' @param Y Response vector (vech of genetic covariance matrix)
#' @param V_Y Covariance matrix of Y
#' @param PREDICTORS Predictor matrix (n x p)
#' @param INTERCEPT Include intercept (default TRUE)
#' @return A numeric matrix with columns \code{betas, pvals, SE, Z},
#'   matching \code{GenomicSEM::summaryGLS}.
#' @examples
#' # GLS regression of a synthetic response on two predictors.
#' Y <- c(0.50, 0.30, 0.20, 0.10)
#' V_Y <- diag(4) * 0.01
#' X <- matrix(c(1, 2, 3, 4,
#'               0, 1, 0, 1), nrow = 4, ncol = 2)
#' fit <- summaryGLS(Y = Y, V_Y = V_Y, PREDICTORS = X, INTERCEPT = TRUE)
#' fit   # matrix with rows b0, b1, b2 and columns betas, pvals, SE, Z
#' @export
summaryGLS <- function(OBJECT = NULL, Y = NULL, V_Y = NULL, PREDICTORS, INTERCEPT = TRUE) {
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

  # Assemble R-compatible shape: matrix(betas, pvals, SE, Z) with
  # row names bN inherited from X's columns (b0 = intercept when
  # requested, otherwise b1..bk from the predictors).
  beta <- as.numeric(result$beta)
  se   <- as.numeric(result$se)
  z    <- as.numeric(result$z)
  p    <- as.numeric(result$p)
  out <- cbind(beta, p, se, z)
  colnames(out) <- c("betas", "pvals", "SE", "Z")
  k <- length(beta)
  rownames(out) <- if (isTRUE(INTERCEPT)) {
    c("b0", if (k >= 2L) paste0("b", seq_len(k - 1L)) else character(0))
  } else {
    paste0("b", seq_len(k))
  }
  print(out)
  invisible(out)
}
