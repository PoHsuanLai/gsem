#' Run Common Factor GWAS
#'
#' Runs multivariate GWAS with an auto-generated 1-factor model per SNP.
#'
#' @param covstruc LDSC result (named list with S, V, I, N, m components)
#' @param SNPs Path to merged summary statistics file
#' @param estimation Estimation method: "DWLS" (default) or "ML"
#' @param cores Number of cores for Rayon thread pool (NULL = auto-detect)
#' @param toler Tolerance (accepted; convergence controlled by L-BFGS internally)
#' @param SNPSE SNP SE override (default FALSE = auto)
#' @param parallel Use parallel processing (default TRUE; FALSE sets single-threaded)
#' @param GC Genomic control: "standard" (default), "conservative", or "none"
#' @param MPI Use MPI (ignored in gsemr -- not applicable to Rust backend)
#' @param TWAS TWAS gene-level analysis mode (default FALSE)
#' @param smooth_check Check for non-positive-definite matrices (default FALSE)
#' @param identification Factor identification strategy: "fixed_variance"
#'   (default, fixes F1 variance to 1) or "marker_indicator" (fixes first
#'   loading to 1, matching R GenomicSEM's convention). Use "marker_indicator"
#'   for exact numerical parity with R GenomicSEM.
#' @return Data frame of per-SNP (or per-Gene if TWAS=TRUE) results
#' @export
# Session-local flag for the one-shot compatibility warning. Lives in the
# package namespace so it resets on reload and is not visible to users.
.gsemr_cfgwas_warned <- FALSE

#' @keywords internal
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

  raw <- if (isTRUE(TWAS)) {
    .snp_columnar_to_df(result, id_col = "Gene", extra_cols = c("Panel", "HSQ"))
  } else {
    .snp_columnar_to_df(result, id_col = "SNP")
  }

  # Extract the SNP/Gene effect row (F1 ~ SNP or F1 ~ Gene) from each
  # row's params, matching R GenomicSEM's commonfactorGWAS output format.
  id_col <- if (isTRUE(TWAS)) "Gene" else "SNP"
  target_rhs <- id_col
  n <- nrow(raw)
  out <- data.frame(
    SNP = raw[[id_col]],
    lhs = character(n), op = character(n), rhs = character(n),
    est = numeric(n), se = numeric(n),
    chisq = raw$chisq, df = raw$df,
    converged = raw$converged,
    stringsAsFactors = FALSE
  )
  if (id_col == "SNP") {
    names(out)[1] <- "SNP"
  } else {
    names(out)[1] <- "Gene"
  }

  for (i in seq_len(n)) {
    p <- raw$params[[i]]
    # Find the SNP/Gene regression row: op == "~" and rhs == target
    snp_row <- which(p$op == "~" & p$rhs == target_rhs)
    if (length(snp_row) > 0) {
      snp_row <- snp_row[1]
      out$lhs[i] <- p$lhs[snp_row]
      out$op[i]  <- p$op[snp_row]
      out$rhs[i] <- p$rhs[snp_row]
      out$est[i] <- p$est[snp_row]
      out$se[i]  <- p$se[snp_row]
    } else {
      out$est[i] <- NA
      out$se[i]  <- NA
    }
  }

  out$z <- ifelse(out$se > 0, out$est / out$se, NA)
  out$p <- ifelse(is.finite(out$z), 2 * pnorm(abs(out$z), lower.tail = FALSE), NA)
  out
}
