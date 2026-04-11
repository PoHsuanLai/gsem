#' Multi-SNP Joint Analysis
#'
#' Fits a model with multiple SNPs simultaneously, accounting for LD between them.
#' Port of R GenomicSEM's \code{multiSNP()}.
#'
#' @param covstruc LDSC result (named list with S, V, I, N, m components)
#' @param SNPs Data frame of SNP summary statistics (with beta, se, var columns per SNP)
#' @param LD LD correlation matrix between SNPs
#' @param SNPSE SNP SE override (default FALSE = auto, numeric = override)
#' @param SNPlist Optional SNP list to subset (character vector)
#' @return A list with converged, chisq, df, and params
#' @examples
#' \dontrun{
#' # `covstruc` from `ldsc()`, `SNPs` from `sumstats()` subset to the
#' # target SNPs, and `LD` from a reference panel (e.g. LDlinkR or PLINK).
#' result <- multiSNP(
#'   covstruc = covstruc,
#'   SNPs = snp_df,       # data.frame with beta, se, var per SNP
#'   LD = ld_matrix,      # n_snps x n_snps correlation matrix
#'   SNPSE = 0.0005
#' )
#' result$results
#' }
#' @export
multiSNP <- function(covstruc, SNPs, LD, SNPSE = FALSE, SNPlist = NA) {

  # Build model: each SNP predicts all latent factors
  k <- nrow(as.matrix(covstruc$S))
  n_snps <- nrow(as.matrix(LD))
  snp_names <- if (!is.null(rownames(as.matrix(LD)))) {
    rownames(as.matrix(LD))
  } else {
    paste0("SNP", seq_len(n_snps))
  }

  # Extract beta, se, var from SNPs data frame and coerce to numeric matrices
  beta_mat <- matrix(as.numeric(as.matrix(SNPs$beta)), nrow = n_snps)
  se_mat <- matrix(as.numeric(as.matrix(SNPs$se)), nrow = n_snps)
  var_snp <- as.numeric(SNPs$var)

  model <- paste0("F1 =~ ", paste0("NA*V", seq_len(k), collapse = " + "),
                   "\n", paste0("F1 ~ ", snp_names, collapse = "\n"),
                   "\nF1 ~~ 1*F1")

  # Convert SNPSE: FALSE means auto (pass NaN), numeric means override
  snp_se_val <- if (is.logical(SNPSE) && !SNPSE) NaN else as.double(SNPSE)

  result <- .Call("wrap__multi_snp_rust",
    .covstruc_as_list(covstruc),
    as.character(model),
    "DWLS",
    beta_mat,
    se_mat,
    as.numeric(var_snp),
    matrix(as.numeric(as.matrix(LD)), nrow = nrow(as.matrix(LD))),
    as.character(snp_names),
    snp_se_val
  )

  if (!is.null(result$error)) stop("gsemr::multiSNP error: ", result$error)
  list(
    converged = result$converged,
    chisq = result$chisq,
    df = result$df,
    results = as.data.frame(result$params, stringsAsFactors = FALSE)
  )
}
