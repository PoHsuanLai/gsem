#' Fit a User-Specified SEM Model
#'
#' Fits a structural equation model to the genetic covariance structure
#' estimated by \code{\link{ldsc}}.
#'
#' @param covstruc Output from \code{\link{ldsc}}
#' @param estimation Estimation method: \code{"DWLS"} (default) or \code{"ML"}
#' @param model lavaan-style model syntax string
#' @param CFIcalc Compute CFI (default TRUE; FALSE skips null model fitting)
#' @param std.lv Standardize latent variables (default FALSE)
#' @param imp_cov Return model-implied covariance matrix (default FALSE)
#' @param fix_resid Fix residual variances to be positive (default TRUE)
#' @param toler Gradient norm tolerance for L-BFGS optimizer (default NULL = 1e-6)
#' @param Q_Factor Compute Q factor heterogeneity statistic (default FALSE)
#' @return A list with components:
#'   \item{results}{Data frame of parameter estimates}
#'   \item{modelfit}{Data frame of fit indices (chisq, df, p_chisq, AIC, CFI, SRMR)}
#'   \item{converged}{Logical indicating convergence}
#'   \item{implied_cov}{Model-implied covariance matrix (if imp_cov=TRUE)}
#'   \item{Q_Factor}{Q factor results (if Q_Factor=TRUE)}
#' @export
usermodel <- function(covstruc, estimation="DWLS", model="", CFIcalc=TRUE,
                      std.lv=FALSE, imp_cov=FALSE, fix_resid=TRUE, toler=NULL, Q_Factor=FALSE) {

  # Convert toler: NULL/FALSE means auto (pass NaN), numeric means override
  toler_val <- if (is.null(toler) || identical(toler, FALSE)) NaN else as.double(toler)

  result <- .Call("wrap__usermodel_rust",
    .covstruc_as_list(covstruc),
    as.character(model),
    as.character(estimation),
    as.logical(std.lv),
    as.logical(fix_resid),
    as.logical(imp_cov),
    as.logical(Q_Factor),
    toler_val,
    as.logical(CFIcalc)
  )

  if (!is.null(result$error)) {
    stop("gsemr::usermodel error: ", result$error)
  }

  params <- as.data.frame(result$parameters, stringsAsFactors = FALSE)

  # Build modelfit data frame from fit indices
  fit_cols <- list(
    chisq = result$chisq,
    chisq_df = result$df,
    p_chisq = result$p_chisq,
    AIC = result$aic,
    SRMR = result$srmr
  )
  if (!is.null(result$cfi)) {
    fit_cols$CFI <- result$cfi
  }

  out <- list(
    results = params,
    modelfit = as.data.frame(fit_cols),
    converged = result$converged
  )

  if (imp_cov && !is.null(result$implied_cov)) {
    out$implied_cov <- result$implied_cov
  }

  if (Q_Factor && !is.null(result$Q_Factor)) {
    out$Q_Factor <- as.data.frame(result$Q_Factor, stringsAsFactors = FALSE)
  }

  out
}

# `commonfactor()` is defined in commonfactor.R as a direct Rust call.
