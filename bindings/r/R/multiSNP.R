#' Multi-SNP Joint Analysis
#'
#' Fits a model with multiple SNPs simultaneously, accounting for LD between them.
#' Port of R GenomicSEM's \code{multiSNP()}.
#'
#' @param covstruc LDSC result (list with S, V, I components, or JSON string)
#' @param SNPs Data frame of SNP summary statistics (with beta, se, var columns per SNP)
#' @param LD LD correlation matrix between SNPs
#' @param SNPSE SNP SE flag (ignored in gsemr)
#' @param SNPlist Optional SNP list to subset (character vector)
#' @return A list with converged, chisq, df, and params
#' @export
multiSNP <- function(covstruc, SNPs, LD, SNPSE = FALSE, SNPlist = NA) {
  if (!identical(SNPSE, FALSE)) {
    message("Note: 'SNPSE' is ignored in gsemr")
  }

  if (is.list(covstruc) && !is.character(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S,
      v = covstruc$V,
      i_mat = covstruc$I,
      n_vec = covstruc$N,
      m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

  # Build model: each SNP predicts all latent factors
  k <- nrow(as.matrix(covstruc$S))
  n_snps <- nrow(as.matrix(LD))
  snp_names <- if (!is.null(rownames(as.matrix(LD)))) {
    rownames(as.matrix(LD))
  } else {
    paste0("SNP", seq_len(n_snps))
  }

  # Extract beta, se, var from SNPs data frame
  beta_mat <- as.matrix(SNPs$beta)
  se_mat <- as.matrix(SNPs$se)
  var_snp <- as.numeric(SNPs$var)

  model <- paste0("F1 =~ ", paste0("NA*V", seq_len(k), collapse = " + "),
                   "\n", paste0("F1 ~ ", snp_names, collapse = "\n"),
                   "\nF1 ~~ 1*F1")

  beta_json <- jsonlite::toJSON(beta_mat, digits = 15)
  se_json <- jsonlite::toJSON(se_mat, digits = 15)
  ld_json <- jsonlite::toJSON(as.matrix(LD), digits = 15)

  json <- .Call("wrap__multi_snp_rust",
    as.character(covstruc_json),
    as.character(model),
    "DWLS",
    as.character(beta_json),
    as.character(se_json),
    as.numeric(var_snp),
    as.character(ld_json),
    as.character(snp_names)
  )

  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::multiSNP error: ", result$error)
  list(
    converged = result$converged,
    chisq = result$chisq,
    df = result$df,
    results = as.data.frame(result$params)
  )
}
