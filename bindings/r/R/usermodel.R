#' Fit a User-Specified SEM Model
#'
#' Fits a structural equation model to the genetic covariance structure
#' estimated by \code{\link{ldsc}}.
#'
#' @param covstruc Output from \code{\link{ldsc}}
#' @param estimation Estimation method: \code{"DWLS"} (default) or \code{"ML"}
#' @param model lavaan-style model syntax string
#' @param CFIcalc Compute CFI (ignored in gsemr, always computed)
#' @param std.lv Standardize latent variables (default FALSE)
#' @param imp_cov Use model-implied covariance (ignored in gsemr)
#' @param fix_resid Fix residual variances to be positive (default TRUE)
#' @param toler Tolerance (accepted; convergence controlled by L-BFGS internally)
#' @param Q_Factor Compute Q factor (ignored in gsemr)
#' @return A list with components:
#'   \item{results}{Data frame of parameter estimates}
#'   \item{modelfit}{Data frame of fit indices}
#'   \item{converged}{Logical indicating convergence}
#' @export
usermodel <- function(covstruc, estimation="DWLS", model="", CFIcalc=TRUE,
                      std.lv=FALSE, imp_cov=FALSE, fix_resid=TRUE, toler=NULL, Q_Factor=FALSE) {

  # Ignored params
  if (!identical(CFIcalc, TRUE)) {
    message("Note: 'CFIcalc' is ignored in gsemr -- CFI is always computed")
  }
  if (!identical(imp_cov, FALSE)) {
    message("Note: 'imp_cov' is ignored in gsemr -- not implemented")
  }
  # Note: 'toler' is accepted but convergence tolerance is controlled by the
  # L-BFGS optimizer internally in the Rust backend.
  if (!identical(Q_Factor, FALSE)) {
    message("Note: 'Q_Factor' is ignored in gsemr -- not implemented")
  }

  # Convert covstruc to JSON for Rust
  covstruc_json <- jsonlite::toJSON(list(
    s = covstruc$S,
    v = covstruc$V,
    i_mat = covstruc$I,
    n_vec = covstruc$N,
    m = covstruc$m
  ), auto_unbox = TRUE, digits = 15)

  json <- .Call("wrap__usermodel_rust",
    as.character(covstruc_json),
    as.character(model),
    as.character(estimation),
    as.logical(std.lv),
    as.logical(fix_resid)
  )

  result <- jsonlite::fromJSON(json)

  if (!is.null(result$error)) {
    stop("gsemr::usermodel error: ", result$error)
  }

  params <- as.data.frame(result$parameters, stringsAsFactors = FALSE)
  if (nrow(params) > 0) {
    params$est <- as.numeric(params$est)
  }

  list(
    results = params,
    modelfit = data.frame(
      objective = result$objective,
      converged = result$converged
    ),
    converged = result$converged
  )
}

#' Fit a Common Factor Model
#'
#' Convenience wrapper that fits a single common factor model to the
#' genetic covariance structure.
#'
#' @param covstruc Output from \code{\link{ldsc}}
#' @param estimation Estimation method: \code{"DWLS"} (default) or \code{"ML"}
#' @return Same as \code{\link{usermodel}}
#' @export
commonfactor <- function(covstruc, estimation = "DWLS") {
  k <- nrow(covstruc$S)
  trait_names <- rownames(covstruc$S)
  if (is.null(trait_names)) {
    trait_names <- paste0("V", seq_len(k))
  }

  # Build common factor model syntax
  loadings <- paste0("NA*", trait_names, collapse = " + ")
  model_lines <- c(
    paste0("F1 =~ ", loadings),
    "F1 ~~ 1*F1",
    paste0(trait_names, " ~~ ", trait_names)
  )
  model <- paste(model_lines, collapse = "\n")

  usermodel(covstruc, estimation = estimation, model = model)
}
