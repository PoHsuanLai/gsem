#' Run Common Factor GWAS
#'
#' Runs multivariate GWAS with an auto-generated 1-factor model per SNP.
#'
#' @param covstruc LDSC result (list or JSON string)
#' @param SNPs Path to merged summary statistics file
#' @param GC Genomic control: "standard" (default), "conservative", or "none"
#' @return Data frame of per-SNP results
#' @export
commonfactorGWAS <- function(covstruc, SNPs, GC = "standard") {
  if (is.list(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S, v = covstruc$V, i_mat = covstruc$I,
      n_vec = covstruc$N, m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

  json <- .Call("wrap__commonfactor_gwas_rust",
    as.character(covstruc_json), as.character(SNPs), as.character(GC))
  jsonlite::fromJSON(json)
}
