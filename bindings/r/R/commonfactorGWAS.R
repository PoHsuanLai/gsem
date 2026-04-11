# Session-local flag for the one-shot compatibility warning. Lives in the
# package namespace so it resets on reload and is not visible to users.
.gsemr_cfgwas_warned <- FALSE

#' @keywords internal
#' @noRd
.warn_commonfactorGWAS_semantics <- function() {
  if (isTRUE(getOption("gsemr.commonfactorGWAS.quiet", FALSE))) return(invisible())
  if (isTRUE(.gsemr_cfgwas_warned)) return(invisible())

  msg <- paste0(
    "gsemr::commonfactorGWAS fits a single-factor model using fixed-variance\n",
    "identification (F1 ~~ 1*F1, loadings free) with the fix_measurement\n",
    "baseline optimization. This is numerically stable and matches\n",
    "GenomicSEM::userGWAS on the equivalent model, but it does NOT\n",
    "numerically match GenomicSEM::commonfactorGWAS, which uses marker-\n",
    "indicator identification and refits the full model per SNP. On real\n",
    "GWAS data the two can disagree in sign and magnitude.\n\n",
    "If you need bit-for-bit parity with R GenomicSEM::commonfactorGWAS,\n",
    "there is no exact replacement currently. If you want stable single-\n",
    "factor GWAS with the same model R userGWAS would fit, this function\n",
    "is the right call. See ARCHITECTURE.md section 3.3 for the full\n",
    "rationale.\n\n",
    "Suppress this warning with: options(gsemr.commonfactorGWAS.quiet = TRUE)"
  )
  warning(msg, call. = FALSE)

  # Flip the session flag via namespace-local assign so subsequent calls
  # within the same session do not re-warn.
  ns <- asNamespace("gsemr")
  tryCatch({
    unlockBinding(".gsemr_cfgwas_warned", ns)
    assign(".gsemr_cfgwas_warned", TRUE, envir = ns)
    lockBinding(".gsemr_cfgwas_warned", ns)
  }, error = function(e) invisible())

  invisible()
}

#' Run Common Factor GWAS
#'
#' Runs multivariate GWAS with an auto-generated 1-factor model per SNP.
#'
#' @param covstruc LDSC result (named list with S, V, I, N, m components)
#' @param SNPs Path to merged summary statistics file
#' @param estimation Estimation method: "DWLS" (default) or "ML"
#' @param cores Integer cap on the rayon worker pool size used for the
#'   per-SNP fit loop. When \code{NULL} (the default) rayon honours
#'   \code{RAYON_NUM_THREADS} if set, else it uses the number of
#'   logical cores reported by the OS. On many-core machines (32+) or
#'   when the underlying BLAS is multithreaded, set this explicitly to
#'   avoid oversubscribing CPUs with nested BLAS threads.
#' @param toler Tolerance (accepted; convergence controlled by L-BFGS internally)
#' @param SNPSE SNP SE override (default FALSE = auto)
#' @param parallel Use a parallel rayon worker pool for the per-SNP
#'   fit loop (default \code{TRUE}). Set to \code{FALSE} to force
#'   single-threaded execution.
#' @param GC Genomic control: "standard" (default), "conservative", or "none"
#' @param MPI Use MPI (ignored in gsemr -- not applicable to Rust backend)
#' @param TWAS TWAS gene-level analysis mode (default FALSE)
#' @param smooth_check Check for non-positive-definite matrices (default FALSE)
#' @param identification Factor identification strategy: "fixed_variance"
#'   (default, fixes F1 variance to 1) or "marker_indicator" (fixes first
#'   loading to 1, matching R GenomicSEM's convention). Use "marker_indicator"
#'   for exact numerical parity with R GenomicSEM.
#' @return Data frame of per-SNP (or per-Gene if TWAS=TRUE) results
#' @note gsemr's \code{commonfactorGWAS} uses fixed-variance identification
#'   and is numerically equivalent to \code{\link{userGWAS}} on the same
#'   1-factor model. It does NOT match R \code{GenomicSEM::commonfactorGWAS},
#'   which uses marker-indicator identification. A first-use warning is
#'   emitted per session; suppress it with
#'   \code{options(gsemr.commonfactorGWAS.quiet = TRUE)}.
#'   See \code{ARCHITECTURE.md} section 3.3 for the full rationale.
#' @examples
#' \dontrun{
#' # Suppress the one-shot compatibility warning if desired:
#' options(gsemr.commonfactorGWAS.quiet = TRUE)
#'
#' # `covstruc` from `ldsc()`, merged SNP file from `sumstats()`:
#' result <- commonfactorGWAS(
#'   covstruc = covstruc,
#'   SNPs = "merged_sumstats.tsv",
#'   estimation = "DWLS",
#'   GC = "standard"
#' )
#' head(result)   # per-SNP est, se, z, p
#' }
#' @export
commonfactorGWAS <- function(covstruc=NULL, SNPs=NULL, estimation="DWLS", cores=NULL,
                             toler=FALSE, SNPSE=FALSE, parallel=TRUE, GC="standard",
                             MPI=FALSE, TWAS=FALSE, smooth_check=FALSE,
                             identification="fixed_variance") {

  .warn_commonfactorGWAS_semantics()

  if (!identical(MPI, FALSE)) {
    message("Note: 'MPI' is ignored in gsemr -- not applicable to Rust backend")
  }
  num_threads <- .resolve_num_threads(parallel, cores)

  # Convert SNPSE: FALSE means auto (pass NaN), numeric means override
  snp_se_val <- if (is.logical(SNPSE) && !SNPSE) NaN else as.double(SNPSE)

  # Accept both data frame (R-compatible) and file path
  if (is.data.frame(SNPs) || is.matrix(SNPs)) {
    tmp <- tempfile(fileext = ".tsv")
    write.table(SNPs, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
    snp_path <- tmp
    on.exit(unlink(tmp), add = TRUE)
  } else {
    snp_path <- as.character(SNPs)
  }

  result <- .Call("wrap__commonfactor_gwas_rust",
    .covstruc_as_list(covstruc), snp_path, as.character(GC),
    as.character(estimation), snp_se_val, as.logical(smooth_check),
    as.logical(TWAS), as.character(identification),
    num_threads)

  if (!is.null(result$error)) {
    stop("gsemr::commonfactorGWAS error: ", result$error)
  }

  raw <- .snp_columnar_to_df(result)

  # Pull the F1 ~ SNP (or F1 ~ Gene) effect row for every SNP via a
  # vectorized filter + hashed match(), then join back to the snps
  # table. `match()` returns NA for SNPs whose fit produced no effect
  # row, which propagates to NA est/se via vec[NA] semantics.
  id_col     <- if (isTRUE(TWAS)) "Gene" else "SNP"
  top_tbl    <- if (isTRUE(TWAS)) raw$genes else raw$snps
  target_rhs <- id_col

  pe       <- raw$params
  eff_mask <- pe$op == "~" & pe$rhs == target_rhs
  eff      <- pe[eff_mask, , drop = FALSE]
  m <- match(top_tbl[[id_col]], eff[[id_col]])

  out <- data.frame(
    setNames(list(top_tbl[[id_col]]), id_col),
    lhs = eff$lhs[m],
    op  = eff$op[m],
    rhs = eff$rhs[m],
    est = eff$est[m],
    se  = eff$se[m],
    chisq = top_tbl$chisq,
    df = top_tbl$df,
    converged = top_tbl$converged,
    stringsAsFactors = FALSE
  )

  if (!is.null(top_tbl$Q_chisq)) {
    out$Q_chisq <- top_tbl$Q_chisq
    out$Q_df    <- top_tbl$Q_df
    out$Q_pval  <- top_tbl$Q_pval
  }

  out$z <- ifelse(is.finite(out$se) & out$se > 0, out$est / out$se, NA_real_)
  out$p <- ifelse(is.finite(out$z), 2 * pnorm(abs(out$z), lower.tail = FALSE), NA_real_)
  out
}
