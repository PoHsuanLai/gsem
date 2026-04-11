#' Auto-Generate Model Syntax from Factor Loadings
#'
#' Generates lavaan-style model syntax from a factor loading matrix.
#'
#' @param Loadings Matrix of factor loadings (rows=phenotypes, cols=factors)
#' @param S_LD Not used (kept for R GenomicSEM compatibility)
#' @param cutoff Minimum absolute loading to include (default 0.3)
#' @param fix_resid Add positivity constraints on residual variances (default TRUE)
#' @param bifactor Generate bifactor model (default FALSE)
#' @param mustload Require all variables to load on at least one factor (default FALSE)
#' @param common Include a common factor (default FALSE)
#' @return Character string of lavaan model syntax
#' @examples
#' # Build a 1-factor loadings matrix and generate lavaan syntax.
#' loadings <- matrix(c(0.80, 0.70, 0.60), nrow = 3, ncol = 1,
#'                    dimnames = list(c("V1", "V2", "V3"), "F1"))
#' cat(write.model(loadings, S_LD = NULL, cutoff = 0.3))
#' @export
write.model <- function(Loadings, S_LD, cutoff, fix_resid=TRUE, bifactor=FALSE,
                        mustload=FALSE, common=FALSE) {

  if (!is.matrix(Loadings)) Loadings <- as.matrix(Loadings)
  names <- rownames(Loadings)
  if (is.null(names)) names <- paste0("V", seq_len(nrow(Loadings)))

  .Call("wrap__write_model_rust",
    as.double(as.vector(t(Loadings))),
    as.integer(nrow(Loadings)),
    as.integer(ncol(Loadings)),
    as.character(names),
    as.double(cutoff),
    as.logical(fix_resid),
    as.logical(bifactor),
    as.logical(mustload),
    as.logical(common)
  )
}
