#' Run User-Specified GWAS
#'
#' Runs multivariate GWAS with a user-specified SEM model per SNP.
#'
#' @param covstruc LDSC result (list or JSON string)
#' @param SNPs Path to merged summary statistics file
#' @param estimation Estimation method: "DWLS" (default) or "ML"
#' @param model lavaan-style model syntax
#' @param printwarn Print warnings (ignored in gsemr)
#' @param sub Subset of results to return (default FALSE = all)
#' @param cores Number of cores (ignored in gsemr)
#' @param toler Tolerance (ignored in gsemr)
#' @param SNPSE SNP SE override (default FALSE = auto)
#' @param parallel Use parallel processing (ignored in gsemr)
#' @param GC Genomic control: "standard" (default), "conservative", or "none"
#' @param MPI Use MPI (ignored in gsemr)
#' @param smooth_check Check for non-positive-definite matrices (default FALSE)
#' @param TWAS TWAS mode (ignored in gsemr)
#' @param std.lv Standardize latent variables (default FALSE)
#' @param fix_measurement Fix measurement model (ignored in gsemr)
#' @param Q_SNP Compute Q_SNP heterogeneity statistic (ignored in gsemr)
#' @return Data frame of per-SNP results
#' @export
userGWAS <- function(covstruc=NULL, SNPs=NULL, estimation="DWLS", model="",
                     printwarn=TRUE, sub=FALSE, cores=NULL, toler=FALSE, SNPSE=FALSE,
                     parallel=TRUE, GC="standard", MPI=FALSE, smooth_check=FALSE,
                     TWAS=FALSE, std.lv=FALSE, fix_measurement=TRUE, Q_SNP=FALSE) {

  # Ignored params
  if (!identical(printwarn, TRUE)) {
    message("Note: 'printwarn' is ignored in gsemr -- warnings are always printed")
  }
  if (!is.null(cores)) {
    message("Note: 'cores' is ignored in gsemr -- Rust uses native parallelism automatically")
  }
  if (!identical(toler, FALSE)) {
    message("Note: 'toler' is ignored in gsemr -- Rust uses its own convergence criteria")
  }
  if (!identical(parallel, TRUE)) {
    message("Note: 'parallel' is ignored in gsemr -- Rust uses native parallelism automatically")
  }
  if (!identical(MPI, FALSE)) {
    message("Note: 'MPI' is ignored in gsemr -- not applicable to Rust backend")
  }
  if (!identical(TWAS, FALSE)) {
    message("Note: 'TWAS' is ignored in gsemr -- not implemented")
  }
  if (!identical(fix_measurement, TRUE)) {
    message("Note: 'fix_measurement' is ignored in gsemr -- not implemented")
  }
  if (!identical(Q_SNP, FALSE)) {
    message("Note: 'Q_SNP' is ignored in gsemr -- not implemented")
  }

  if (is.list(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S, v = covstruc$V, i_mat = covstruc$I,
      n_vec = covstruc$N, m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

  # Convert sub: FALSE means no subset, otherwise pass as string
  sub_str <- if (is.logical(sub) && !sub) "" else as.character(sub)

  # Convert SNPSE: FALSE means auto (pass NaN), numeric means override
  snp_se_val <- if (is.logical(SNPSE) && !SNPSE) NaN else as.double(SNPSE)

  json <- .Call("wrap__user_gwas_rust",
    as.character(covstruc_json), as.character(SNPs),
    as.character(model), as.character(estimation), as.character(GC),
    as.character(sub_str), snp_se_val, as.logical(smooth_check), as.logical(std.lv))
  jsonlite::fromJSON(json)
}
