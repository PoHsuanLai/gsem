#' Parallel Analysis on LDSC Results
#'
#' Determines the number of non-spurious latent factors via Monte Carlo
#' eigenvalue comparison against simulated null distributions.
#'
#' @param S Genetic covariance matrix (from LDSC output)
#' @param V Sampling covariance matrix (from LDSC output)
#' @param r Number of Monte Carlo simulations (default NULL = 500)
#' @param p Percentile threshold for null distribution (default NULL = 0.95)
#' @param save.pdf Save scree plot to PDF file (default FALSE; set to file path string to save)
#' @param diag Use diagonal of V only (default FALSE)
#' @param fa Factor analysis method (ignored in gsemr)
#' @param fm Factor method (ignored in gsemr)
#' @param nfactors Number of factors to extract (ignored in gsemr)
#' @return A list with components:
#'   \item{observed}{Observed eigenvalues (descending)}
#'   \item{simulated_95}{Simulated eigenvalues at the given percentile}
#'   \item{n_factors}{Suggested number of factors}
#' @export
paLDSC <- function(S=S, V=V, r=NULL, p=NULL, save.pdf=FALSE, diag=FALSE, fa=FALSE,
                   fm=NULL, nfactors=NULL) {

  if (!identical(fa, FALSE)) {
    message("Note: 'fa' is ignored in gsemr -- not implemented")
  }
  if (!is.null(fm)) {
    message("Note: 'fm' is ignored in gsemr -- not implemented")
  }
  if (!is.null(nfactors)) {
    message("Note: 'nfactors' is ignored in gsemr -- not implemented")
  }

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
