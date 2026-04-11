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
#' @examples
#' \dontrun{
#' # Stratified LDSC requires annotation-specific LD, weight, and freq files.
#' result <- s_ldsc(
#'   traits = c("T1.sumstats.gz", "T2.sumstats.gz"),
#'   ld = "baseline_LD/",
#'   wld = "weights/",
#'   frq = "frq/",
#'   trait.names = c("V1", "V2")
#' )
#' # Feed into `enrich()` for annotation enrichment tests.
#' }
#' @export
s_ldsc <- function(traits, sample.prev=NULL, population.prev=NULL, ld, wld, frq,
                   trait.names=NULL, n.blocks=200, ldsc.log=NULL, exclude_cont=TRUE) {

  # Ignored params
  # ldsc.log: handled after computation below

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

  result <- .Call("wrap__s_ldsc_rust",
    as.character(traits),
    as.double(sample.prev),
    as.double(population.prev),
    as.character(ld),
    as.character(wld),
    as.character(frq),
    as.integer(n.blocks),
    as.logical(exclude_cont)
  )

  if (!is.null(result$error)) {
    stop("gsemr::s_ldsc error: ", result$error)
  }

  # Apply trait names to the intercept matrix and to each per-annotation
  # S matrix. V matrices are kstar x kstar so they get trait-pair names.
  if (!is.null(result$I)) {
    I <- result$I
    if (nrow(I) == length(trait.names)) {
      rownames(I) <- colnames(I) <- trait.names
    }
    result$I <- I
  }
  if (!is.null(result$S_annot)) {
    result$S_annot <- lapply(result$S_annot, function(S) {
      if (nrow(S) == length(trait.names)) {
        rownames(S) <- colnames(S) <- trait.names
      }
      S
    })
  }
  if (!is.null(result$V_annot)) {
    k <- length(trait.names)
    vpairs <- c()
    for (i in seq_len(k)) {
      for (j in i:k) {
        vpairs <- c(vpairs, paste0(trait.names[i], "_", trait.names[j]))
      }
    }
    result$V_annot <- lapply(result$V_annot, function(V) {
      if (nrow(V) == length(vpairs)) {
        rownames(V) <- colnames(V) <- vpairs
      }
      V
    })
  }

  # Write log file if requested
  if (!is.null(ldsc.log)) {
    sink(ldsc.log)
    cat("gsemr Stratified LDSC Results\n")
    cat("=============================\n\n")
    cat("Traits:", paste(traits, collapse=", "), "\n")
    if (!is.null(result$I)) {
      cat("\nIntercept Matrix (I):\n")
      print(round(result$I, 4))
    }
    if (!is.null(result$annotations)) {
      cat("\nAnnotations:", paste(result$annotations, collapse=", "), "\n")
    }
    sink()
    message("s_ldsc log written to: ", ldsc.log)
  }

  result
}
