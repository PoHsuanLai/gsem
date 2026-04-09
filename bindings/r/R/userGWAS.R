#' Run User-Specified GWAS
#'
#' Runs multivariate GWAS with a user-specified SEM model per SNP.
#'
#' @param covstruc LDSC result (list or JSON string)
#' @param SNPs Path to merged summary statistics file
#' @param model lavaan-style model syntax
#' @param estimation Estimation method: "DWLS" (default) or "ML"
#' @param GC Genomic control: "standard" (default), "conservative", or "none"
#' @return Data frame of per-SNP results
#' @export
userGWAS <- function(covstruc, SNPs, model, estimation = "DWLS", GC = "standard") {
  if (is.list(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S, v = covstruc$V, i_mat = covstruc$I,
      n_vec = covstruc$N, m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

  json <- .Call("wrap__user_gwas_rust",
    as.character(covstruc_json), as.character(SNPs),
    as.character(model), as.character(estimation), as.character(GC))
  jsonlite::fromJSON(json)
}
