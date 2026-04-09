#' Parallel Analysis on LDSC Results
#'
#' Determines the number of non-spurious latent factors via Monte Carlo
#' eigenvalue comparison against simulated null distributions.
#'
#' When \code{fa=TRUE} or \code{fm} is specified, uses \pkg{psych} package's
#' factor analysis methods. Otherwise uses the fast Rust eigenvalue-based
#' parallel analysis.
#'
#' @param S Genetic covariance matrix (from LDSC output)
#' @param V Sampling covariance matrix (from LDSC output)
#' @param r Number of Monte Carlo simulations (default NULL = 500)
#' @param p Percentile threshold for null distribution (default NULL = 0.95)
#' @param save.pdf Save scree plot to PDF file (default FALSE; set to file path string to save)
#' @param diag Use diagonal of V only (default FALSE)
#' @param fa Use factor analysis instead of eigenvalue comparison (default FALSE; requires psych)
#' @param fm Factor method for psych::fa (e.g. "minres", "ml", "pa"; default NULL = "minres")
#' @param nfactors Number of factors to extract when using fa (default NULL = auto from eigenvalue analysis)
#' @return A list with components:
#'   \item{observed}{Observed eigenvalues (descending)}
#'   \item{simulated_95}{Simulated eigenvalues at the given percentile}
#'   \item{n_factors}{Suggested number of factors}
#' @export
paLDSC <- function(S=S, V=V, r=NULL, p=NULL, save.pdf=FALSE, diag=FALSE, fa=FALSE,
                   fm=NULL, nfactors=NULL) {

  # Default r to 500 if NULL
  if (is.null(r)) r <- 500L

  # Convert S and V matrices to JSON
  s_json <- jsonlite::toJSON(as.matrix(S), digits = 15)
  v_json <- jsonlite::toJSON(as.matrix(V), digits = 15)

  # Convert p: NULL means default 0.95
  p_val <- if (is.null(p)) NaN else as.double(p)

  json <- .Call("wrap__pa_ldsc_rust",
    as.character(s_json), as.character(v_json), as.integer(r),
    p_val, as.logical(diag))
  result <- jsonlite::fromJSON(json)

  # Factor analysis mode: delegate to psych::fa() for the actual extraction
  if (!identical(fa, FALSE) || !is.null(fm)) {
    if (!requireNamespace("psych", quietly = TRUE)) {
      stop("The 'psych' package is required for fa/fm options. Install with: install.packages('psych')")
    }

    # Convert S to correlation matrix
    s_mat <- as.matrix(S)
    sds <- sqrt(pmax(base::diag(s_mat), 1e-10))
    cor_mat <- s_mat / outer(sds, sds)
    base::diag(cor_mat) <- 1.0

    # Determine number of factors from eigenvalue analysis if not specified
    n_fa <- if (!is.null(nfactors)) nfactors else result$n_factors

    # Factor method
    fa_method <- if (!is.null(fm)) fm else "minres"

    fa_result <- tryCatch(
      psych::fa(cor_mat, nfactors = n_fa, fm = fa_method, rotate = "none"),
      error = function(e) {
        warning("psych::fa() failed: ", e$message, ". Falling back to eigenvalue analysis.")
        NULL
      }
    )

    if (!is.null(fa_result)) {
      result$fa <- list(
        loadings = unclass(fa_result$loadings),
        uniquenesses = fa_result$uniquenesses,
        fm = fa_method,
        nfactors = n_fa
      )
    }
  }

  # Generate scree plot PDF if requested
  if (!identical(save.pdf, FALSE)) {
    pdf_path <- if (is.character(save.pdf)) save.pdf else "paLDSC_scree.pdf"
    k <- length(result$observed)
    grDevices::pdf(pdf_path, width = 7, height = 5)
    graphics::plot(seq_len(k), result$observed, type = "b", pch = 16, col = "blue",
         xlab = "Factor", ylab = "Eigenvalue", main = "Parallel Analysis Scree Plot",
         ylim = range(c(result$observed, result$simulated_95)))
    graphics::lines(seq_len(k), result$simulated_95, type = "b", pch = 17, col = "red", lty = 2)
    graphics::legend("topright", legend = c("Observed", "Simulated"), col = c("blue", "red"),
           pch = c(16, 17), lty = c(1, 2))
    grDevices::dev.off()
    message("Scree plot saved to: ", pdf_path)
  }

  result
}
