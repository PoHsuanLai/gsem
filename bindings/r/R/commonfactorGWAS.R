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
#' @param parallel Use parallel processing (ignored in gsemr)
#' @param GC Genomic control: "standard" (default), "conservative", or "none"
#' @param MPI Use MPI (ignored in gsemr)
#' @param TWAS TWAS mode (ignored in gsemr)
#' @param smooth_check Check for non-positive-definite matrices (default FALSE)
#' @return Data frame of per-SNP results
#' @export
commonfactorGWAS <- function(covstruc=NULL, SNPs=NULL, estimation="DWLS", cores=NULL,
                             toler=FALSE, SNPSE=FALSE, parallel=TRUE, GC="standard",
                             MPI=FALSE, TWAS=FALSE, smooth_check=FALSE) {

  if (!identical(parallel, TRUE)) {
    Sys.setenv(RAYON_NUM_THREADS = "1")
  }
  if (!identical(MPI, FALSE)) {
    message("Note: 'MPI' is ignored in gsemr -- not applicable to Rust backend")
  }
  if (!identical(TWAS, FALSE)) {
    message("Note: 'TWAS' is ignored in gsemr -- not implemented")
  }

  # Set rayon thread count if cores is specified
  if (!is.null(cores) && cores > 0) {
    Sys.setenv(RAYON_NUM_THREADS = as.character(cores))
  }

  # Note: 'toler' is accepted but convergence tolerance is controlled by the
  # L-BFGS optimizer internally in the Rust backend.

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

  json <- .Call("wrap__commonfactor_gwas_rust",
    as.character(covstruc_json), as.character(SNPs), as.character(GC),
    as.character(estimation), snp_se_val, as.logical(smooth_check))
  jsonlite::fromJSON(json)
}
