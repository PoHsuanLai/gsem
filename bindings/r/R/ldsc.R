#' Run LD Score Regression
#'
#' Estimates genetic covariance matrix (S), sampling covariance matrix (V),
#' and intercept matrix (I) from GWAS summary statistics using LD Score
#' Regression with block jackknife.
#'
#' @param traits Character vector of paths to .sumstats.gz files
#' @param sample.prev Numeric vector of sample prevalences (NA for continuous traits)
#' @param population.prev Numeric vector of population prevalences (NA for continuous)
#' @param ld Path to LD score directory (containing {chr}.l2.ldscore.gz files)
#' @param wld Path to weight LD score directory (defaults to \code{ld})
#' @param trait.names Character vector of trait names (defaults to V1, V2, ...)
#' @param sep_weights Use separate weight LD scores (ignored in gsemr)
#' @param chr Number of chromosomes (default 22)
#' @param n.blocks Number of jackknife blocks (default 200)
#' @param ldsc.log Log file path (ignored in gsemr)
#' @param stand Standardize output (default FALSE)
#' @param select Variable selection method (default FALSE)
#' @param chisq.max Maximum chi-square filter (default NA = auto)
#' @param parallel Use a parallel rayon worker pool for the per-pair
#'   regression loop (default \code{TRUE}). Set to \code{FALSE} to force
#'   single-threaded execution.
#' @param cores Integer cap on the rayon pool size. When \code{NULL}
#'   (the default) rayon honours \code{RAYON_NUM_THREADS} if set, else
#'   it uses the number of logical cores reported by the OS. On
#'   many-core machines (32+) or when the underlying BLAS is
#'   multithreaded, set this explicitly to avoid oversubscribing CPUs
#'   with nested BLAS threads.
#' @return A list with components:
#'   \item{S}{Genetic covariance matrix (k x k)}
#'   \item{V}{Sampling covariance matrix of S (k* x k*, where k* = k(k+1)/2)}
#'   \item{I}{LDSC intercept matrix (k x k)}
#'   \item{N}{Sample size vector}
#'   \item{m}{Number of SNPs used}
#' @export
ldsc <- function(traits, sample.prev, population.prev, ld, wld,
                 trait.names=NULL, sep_weights=FALSE, chr=22,
                 n.blocks=200, ldsc.log=NULL, stand=FALSE, select=FALSE,
                 chisq.max=NA, parallel=TRUE, cores=NULL) {

  # sep_weights: when FALSE, use ld directory for weights too (ignore wld)
  if (identical(sep_weights, FALSE)) {
    wld <- ld
  }
  # ldsc.log: handled after computation below

  if (is.null(trait.names)) {
    trait.names <- paste0("V", seq_along(traits))
  }

  # Convert select to string for Rust
  select_str <- if (is.logical(select) && !select) "FALSE" else as.character(select)

  # Convert chisq.max: NA means auto (pass NaN to Rust)
  chisq_max_val <- if (is.na(chisq.max)) NaN else as.double(chisq.max)

  num_threads <- .resolve_num_threads(parallel, cores)

  result <- .Call("wrap__ldsc_rust",
    as.character(traits),
    as.double(sample.prev),
    as.double(population.prev),
    as.character(ld),
    as.character(wld),
    as.integer(n.blocks),
    as.integer(chr),
    chisq_max_val,
    as.logical(stand),
    as.character(select_str),
    num_threads
  )

  if (!is.null(result$error)) {
    stop("gsemr::ldsc error: ", result$error)
  }

  S <- result$s
  rownames(S) <- colnames(S) <- trait.names

  V <- result$v
  I <- result$i_mat
  rownames(I) <- colnames(I) <- trait.names

  out <- list(
    S = S,
    V = V,
    I = I,
    N = result$n_vec,
    m = result$m
  )

  if (!is.null(result$s_stand)) {
    s_stand <- result$s_stand
    rownames(s_stand) <- colnames(s_stand) <- trait.names
    out$S_Stand <- s_stand
  }
  if (!is.null(result$v_stand)) {
    out$V_Stand <- result$v_stand
  }

  # Write log file if requested
  if (!is.null(ldsc.log)) {
    sink(ldsc.log)
    cat("gsemr LDSC Results\n")
    cat("==================\n\n")
    cat("Traits:", paste(traits, collapse=", "), "\n")
    cat("N blocks:", n.blocks, "\n")
    cat("M (SNPs used):", result$m, "\n\n")
    cat("Genetic Covariance Matrix (S):\n")
    print(round(S, 4))
    cat("\nIntercept Matrix (I):\n")
    print(round(I, 4))
    cat("\nSample Sizes (N):\n")
    print(result$n_vec)
    sink()
    message("LDSC log written to: ", ldsc.log)
  }

  out
}
