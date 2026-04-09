#' Stratified LD Score Regression
#'
#' Extends standard LDSC to multiple annotation categories, estimating
#' per-annotation contributions to genetic covariance and heritability.
#'
#' @param traits Character vector of paths to .sumstats.gz files
#' @param sample.prev Numeric vector of sample prevalences (NULL for continuous traits)
#' @param population.prev Numeric vector of population prevalences (NULL for continuous)
#' @param ld Path to annotation LD score directory
#' @param wld Path to weight LD score directory
#' @param frq Path to frequency file directory
#' @param trait.names Character vector of trait names (defaults to V1, V2, ...)
#' @param n.blocks Number of jackknife blocks (default 200)
#' @param ldsc.log Log file path (ignored in gsemr)
#' @param exclude_cont Exclude continuous annotations (default TRUE)
#' @return A list with per-annotation S and V matrices
#' @export
s_ldsc <- function(traits, sample.prev=NULL, population.prev=NULL, ld, wld, frq,
                   trait.names=NULL, n.blocks=200, ldsc.log=NULL, exclude_cont=TRUE) {

  # Ignored params
  if (!is.null(ldsc.log)) {
    message("Note: 'ldsc.log' is ignored in gsemr -- logging is handled by the Rust runtime")
  }

  if (is.null(trait.names)) {
    trait.names <- paste0("V", seq_along(traits))
  }

  # Convert sample.prev/population.prev: NULL -> vector of NA
  if (is.null(sample.prev)) {
    sample.prev <- rep(NA, length(traits))
  }
  if (is.null(population.prev)) {
    population.prev <- rep(NA, length(traits))
  }

  json <- .Call("wrap__s_ldsc_rust",
    as.character(traits),
    as.double(sample.prev),
    as.double(population.prev),
    as.character(ld),
    as.character(wld),
    as.character(frq),
    as.integer(n.blocks),
    as.logical(exclude_cont)
  )

  result <- jsonlite::fromJSON(json)

  if (!is.null(result$error)) {
    stop("gsemr::s_ldsc error: ", result$error)
  }

  result
}
