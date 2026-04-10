#' Run Common Factor GWAS
#'
#' Runs multivariate GWAS with an auto-generated 1-factor model per SNP.
#'
#' @param covstruc LDSC result (list or JSON string)
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
commonfactorGWAS <- function(covstruc=NULL, SNPs=NULL, estimation="DWLS", cores=NULL,
                             toler=FALSE, SNPSE=FALSE, parallel=TRUE, GC="standard",
                             MPI=FALSE, TWAS=FALSE, smooth_check=FALSE,
                             identification="fixed_variance") {

  if (!identical(parallel, TRUE)) {
    Sys.setenv(RAYON_NUM_THREADS = "1")
  }
  if (!identical(MPI, FALSE)) {
    message("Note: 'MPI' is ignored in gsemr -- not applicable to Rust backend")
  }

  # Set rayon thread count if cores is specified
  if (!is.null(cores) && cores > 0) {
    Sys.setenv(RAYON_NUM_THREADS = as.character(cores))
  }

  if (is.list(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S, v = covstruc$V, i_mat = covstruc$I,
      n_vec = covstruc$N, m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

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

  json <- .Call("wrap__commonfactor_gwas_rust",
    as.character(covstruc_json), snp_path, as.character(GC),
    as.character(estimation), snp_se_val, as.logical(smooth_check),
    as.logical(TWAS), as.character(identification))
  raw <- jsonlite::fromJSON(json)

  # Extract the SNP effect row (F1 ~ SNP) from each SNP's params,
  # matching R GenomicSEM's commonfactorGWAS output format.
  n <- nrow(raw)
  out <- data.frame(
    SNP = raw$SNP,
    lhs = character(n), op = character(n), rhs = character(n),
    est = numeric(n), se = numeric(n),
    chisq = raw$chisq, df = raw$df,
    converged = raw$converged,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n)) {
    p <- raw$params[[i]]
    # Find the SNP regression row: op == "~" and rhs == "SNP"
    snp_row <- which(p$op == "~" & p$rhs == "SNP")
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
