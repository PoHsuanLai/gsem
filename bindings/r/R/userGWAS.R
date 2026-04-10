#' Run User-Specified GWAS
#'
#' Runs multivariate GWAS with a user-specified SEM model per SNP.
#'
#' @param covstruc LDSC result (named list with S, V, I, N, m components)
#' @param SNPs Path to merged summary statistics file
#' @param estimation Estimation method: "DWLS" (default) or "ML"
#' @param model lavaan-style model syntax
#' @param printwarn Print warnings (default TRUE; FALSE suppresses Rust warnings)
#' @param sub Subset of results to return (default FALSE = all)
#' @param cores Number of cores for Rayon thread pool (NULL = auto-detect)
#' @param toler Tolerance (accepted; convergence controlled by L-BFGS internally)
#' @param SNPSE SNP SE override (default FALSE = auto)
#' @param parallel Use parallel processing (default TRUE; FALSE sets single-threaded)
#' @param GC Genomic control: "standard" (default), "conservative", or "none"
#' @param MPI Use MPI (ignored in gsemr -- not applicable to Rust backend)
#' @param smooth_check Check for non-positive-definite matrices (default FALSE)
#' @param TWAS TWAS gene-level analysis mode (default FALSE)
#' @param std.lv Standardize latent variables (default FALSE)
#' @param fix_measurement Fix measurement model (default TRUE)
#' @param Q_SNP Compute Q_SNP heterogeneity statistic (default FALSE)
#' @return Data frame of per-SNP (or per-Gene if TWAS=TRUE) results
#' @export
userGWAS <- function(covstruc=NULL, SNPs=NULL, estimation="DWLS", model="",
                     printwarn=TRUE, sub=FALSE, cores=NULL, toler=FALSE, SNPSE=FALSE,
                     parallel=TRUE, GC="standard", MPI=FALSE, smooth_check=FALSE,
                     TWAS=FALSE, std.lv=FALSE, fix_measurement=TRUE, Q_SNP=FALSE) {

  # Suppress warnings if requested
  old_rust_log <- Sys.getenv("RUST_LOG", unset = NA)
  if (!identical(printwarn, TRUE)) {
    Sys.setenv(RUST_LOG = "error")
    on.exit(
      if (is.na(old_rust_log)) {
        Sys.unsetenv("RUST_LOG")
      } else {
        Sys.setenv(RUST_LOG = old_rust_log)
      },
      add = TRUE
    )
  }
  if (!identical(MPI, FALSE)) {
    message("Note: 'MPI' is ignored in gsemr -- not applicable to Rust backend")
  }
  num_threads <- .resolve_num_threads(parallel, cores)

  # Convert sub: FALSE means no subset, otherwise pass as string
  sub_str <- if (is.logical(sub) && !sub) "" else as.character(sub)

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

  result <- .Call("wrap__user_gwas_rust",
    .covstruc_as_list(covstruc), snp_path,
    as.character(model), as.character(estimation), as.character(GC),
    as.character(sub_str), snp_se_val, as.logical(smooth_check), as.logical(std.lv),
    as.logical(fix_measurement), as.logical(Q_SNP), as.logical(TWAS),
    num_threads)

  if (!is.null(result$error)) {
    stop("gsemr::userGWAS error: ", result$error)
  }

  if (isTRUE(TWAS)) {
    .snp_columnar_to_df(result, id_col = "Gene", extra_cols = c("Panel", "HSQ"))
  } else {
    .snp_columnar_to_df(result, id_col = "SNP")
  }
}
