#' Parallel Analysis on LDSC Results
#'
#' Determines the number of non-spurious latent factors via Monte Carlo
#' eigenvalue comparison against simulated null distributions. Drop-in
#' replacement for \code{GenomicSEM::paLDSC()}: prints the same
#' "Parallel Analysis suggests extracting N components" banner and
#' writes the same canonical PDF filenames (\code{PA_LDSC.pdf},
#' \code{Diagonalized_PA_LDSC.pdf}, \code{FA_PA_LDSC.pdf},
#' \code{FA_Diagonalized_PA_LDSC.pdf}) when \code{save.pdf = TRUE}. The
#' return value is returned invisibly so interactive callers don't see
#' the rich list dumped to the console, but if you assign it you get
#' the structured fields described below.
#'
#' When \code{fa=TRUE} or \code{fm} is specified, uses \pkg{psych}
#' package's factor analysis methods for the additional FA panels.
#'
#' @param S Genetic covariance matrix (from LDSC output)
#' @param V Sampling covariance matrix (from LDSC output)
#' @param r Number of Monte Carlo simulations (default NULL = 500)
#' @param p Percentile threshold for null distribution (default NULL = 0.95)
#' @param save.pdf When \code{TRUE}, write the canonical PDF filenames
#'   matching R \code{GenomicSEM::paLDSC} to the working directory.
#' @param diag Additionally run the analysis against \code{diag(V)} only
#'   and print/save the diagonalized variant (default FALSE).
#' @param fa Additionally run \pkg{psych} factor analysis on the
#'   correlation matrix (default FALSE; requires psych).
#' @param fm Factor method for psych::fa (e.g. "minres", "ml", "pa"; default NULL = "minres")
#' @param nfactors Number of factors to extract when using fa (default NULL = auto from eigenvalue analysis)
#' @param parallel Use a parallel rayon worker pool for the Monte
#'   Carlo simulation loop (default \code{TRUE}). Set to \code{FALSE}
#'   to force single-threaded execution.
#' @param cores Integer cap on the rayon pool size. When \code{NULL}
#'   (the default) rayon honours \code{RAYON_NUM_THREADS} if set, else
#'   it uses the number of logical cores reported by the OS. On
#'   many-core machines (32+) or when the underlying BLAS is
#'   multithreaded, set this explicitly to avoid oversubscribing CPUs
#'   with nested BLAS threads.
#' @return Invisibly, a list with components:
#'   \item{observed}{Observed eigenvalues (descending)}
#'   \item{simulated_95}{Simulated eigenvalues at the given percentile}
#'   \item{n_factors}{Suggested number of factors from the
#'     eigenvalue-based parallel analysis}
#'   \item{diag}{(when \code{diag=TRUE}) the same three fields for the
#'     diagonalized variant}
#'   \item{fa}{(when \code{fa=TRUE}) loadings + uniquenesses from
#'     \pkg{psych}}
#' @export
paLDSC <- function(S=S, V=V, r=NULL, p=NULL, save.pdf=FALSE, diag=FALSE, fa=FALSE,
                   fm=NULL, nfactors=NULL, parallel=TRUE, cores=NULL) {

  if (is.null(r)) r <- 500L

  s_mat <- matrix(as.numeric(as.matrix(S)), nrow = nrow(as.matrix(S)))
  v_mat <- matrix(as.numeric(as.matrix(V)), nrow = nrow(as.matrix(V)))

  p_val <- if (is.null(p)) NaN else as.double(p)

  num_threads <- .resolve_num_threads(parallel, cores)

  # --- Baseline (non-diagonal) run — always executed, its n_factors
  # is the top-level suggestion users see first.
  result <- .Call("wrap__pa_ldsc_rust",
    s_mat, v_mat, as.integer(r),
    p_val, FALSE, num_threads)

  if (!is.null(result$error)) {
    stop("gsemr::paLDSC error: ", result$error)
  }

  # --- Diagonalized variant: run a second time against diag(V) only,
  # matching R's behaviour of reporting both suggestions when diag=TRUE.
  result_diag <- NULL
  if (isTRUE(diag)) {
    result_diag <- .Call("wrap__pa_ldsc_rust",
      s_mat, v_mat, as.integer(r),
      p_val, TRUE, num_threads)
    if (!is.null(result_diag$error)) {
      stop("gsemr::paLDSC (diagonalized) error: ", result_diag$error)
    }
    result$diag <- list(
      observed     = result_diag$observed,
      simulated_95 = result_diag$simulated_95,
      n_factors    = result_diag$n_factors
    )
  }

  # --- Factor analysis: delegates to psych::fa() for the actual
  # extraction. Enabled when the user asks for fa=TRUE or specifies fm.
  fa_enabled <- !identical(fa, FALSE) || !is.null(fm)
  fa_result  <- NULL
  if (fa_enabled) {
    if (!requireNamespace("psych", quietly = TRUE)) {
      stop("The 'psych' package is required for fa/fm options. Install with: install.packages('psych')")
    }

    s_mat_full <- as.matrix(S)
    sds <- sqrt(pmax(base::diag(s_mat_full), 1e-10))
    cor_mat <- s_mat_full / outer(sds, sds)
    base::diag(cor_mat) <- 1.0

    n_fa <- if (!is.null(nfactors)) nfactors else result$n_factors
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
        loadings     = unclass(fa_result$loadings),
        uniquenesses = fa_result$uniquenesses,
        fm           = fa_method,
        nfactors     = n_fa
      )
    }
  }

  # --- Console output matching R GenomicSEM::paLDSC's cat() lines.
  .banner <- function() {
    cat(strrep("-", 72), "\n", sep = "")
  }
  .banner()
  cat("Parallel Analysis suggests extracting", result$n_factors, "components", "\n")
  .banner()
  if (!is.null(result_diag)) {
    cat("Diagonalized Parallel Analysis suggests extracting",
        result_diag$n_factors, "components", "\n")
    .banner()
  }
  if (fa_enabled) {
    cat("Parallel Analysis suggests extracting", result$n_factors, "factors", "\n")
    .banner()
    if (isTRUE(diag) && !is.null(result_diag)) {
      cat("Diagonalized Parallel Analysis suggests extracting",
          result_diag$n_factors, "factors", "\n")
      .banner()
    }
  }

  # --- PDF output using the same canonical filenames R writes. We use
  # a simple base-graphics scree plot instead of R's ggpubr/ggarrange
  # two-panel layout to avoid pulling in a ggplot2 + egg + ggpubr
  # dependency chain. The filenames and the existence of the files is
  # what drop-in callers depend on.
  .write_scree_pdf <- function(path, observed, simulated, title) {
    k <- length(observed)
    grDevices::pdf(path, width = 7, height = 5)
    on.exit(invisible(grDevices::dev.off()), add = TRUE)
    graphics::plot(seq_len(k), observed, type = "b", pch = 16, col = "blue",
                   xlab = "Factor", ylab = "Eigenvalue", main = title,
                   ylim = range(c(observed, simulated)))
    graphics::lines(seq_len(k), simulated, type = "b", pch = 17, col = "red", lty = 2)
    graphics::legend("topright", legend = c("Observed", "Simulated"),
                     col = c("blue", "red"), pch = c(16, 17), lty = c(1, 2))
  }

  if (isTRUE(save.pdf)) {
    .write_scree_pdf("PA_LDSC.pdf",
                     result$observed, result$simulated_95,
                     "Parallel Analysis (PA_LDSC)")
    if (!is.null(result_diag)) {
      .write_scree_pdf("Diagonalized_PA_LDSC.pdf",
                       result_diag$observed, result_diag$simulated_95,
                       "Diagonalized Parallel Analysis")
    }
    if (fa_enabled) {
      .write_scree_pdf("FA_PA_LDSC.pdf",
                       result$observed, result$simulated_95,
                       "Factor Analysis Parallel Analysis")
      if (isTRUE(diag) && !is.null(result_diag)) {
        .write_scree_pdf("FA_Diagonalized_PA_LDSC.pdf",
                         result_diag$observed, result_diag$simulated_95,
                         "Factor Analysis Diagonalized Parallel Analysis")
      }
    }
  }

  invisible(result)
}
