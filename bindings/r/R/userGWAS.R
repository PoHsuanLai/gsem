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
#' @param smooth_check Check for non-positive-definite matrices (default FALSE)
#' @param TWAS TWAS gene-level analysis mode (default FALSE)
#' @param std.lv Standardize latent variables (default FALSE)
#' @param fix_measurement Fix measurement model (default TRUE)
#' @param Q_SNP Compute Q_SNP heterogeneity statistic (default FALSE)
#' @return A two-element list in normalized relational form:
#'   \describe{
#'     \item{\code{snps}}{One row per SNP with columns \code{SNP, CHR, BP,
#'       MAF, A1, A2, chisq, df, converged} (plus \code{Q_chisq, Q_df,
#'       Q_pval} when \code{Q_SNP = TRUE}).}
#'     \item{\code{params}}{Flat parameter table with columns \code{SNP,
#'       lhs, op, rhs, est, se, z, p}; join to \code{snps} on \code{SNP}.}
#'   }
#'   When \code{TWAS = TRUE} the first slot is named \code{genes} with
#'   columns \code{Gene, Panel, HSQ, chisq, df, converged}; the
#'   \code{params} table uses \code{Gene} in place of \code{SNP}.
#' @note This replaces the earlier one-row-per-SNP-with-list-column
#'   shape; the old layout made R reconstruction O(N) `data.frame()`
#'   calls and hung at multi-million-SNP scale. Example filters:
#'   `result$params[result$params$SNP == "rsXXXX", ]` for one SNP, or
#'   `result$params[result$params$op == "~" & result$params$rhs == "SNP", ]`
#'   for every SNP's regression effect.
#' @examples
#' \dontrun{
#' # `covstruc` from `ldsc()`, `merged_sumstats.tsv` from `sumstats()`.
#' result <- userGWAS(
#'   covstruc = covstruc,
#'   SNPs = "merged_sumstats.tsv",
#'   model = "F1 =~ NA*V1 + V2 + V3\nF1 ~ SNP\nF1 ~~ 1*F1"
#' )
#' # Per-SNP metadata:
#' head(result$snps)
#' # SNP regression effects (F1 ~ SNP) across all SNPs:
#' head(result$params[result$params$op == "~" & result$params$rhs == "SNP", ])
#' }
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
  .snp_columnar_to_df(result)
}
