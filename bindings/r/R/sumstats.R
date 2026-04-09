#' Merge GWAS Summary Statistics
#'
#' Reads multiple GWAS files, applies QC filters, aligns alleles,
#' and outputs a merged tab-delimited file for multivariate GWAS.
#'
#' @param files Character vector of GWAS file paths
#' @param ref Path to LD score reference directory
#' @param trait.names Character vector of trait names
#' @param info.filter INFO score filter (default 0.6)
#' @param maf.filter MAF filter (default 0.01)
#' @param keep.indel Keep indels (default FALSE)
#' @param out Output file path (default "merged_sumstats.tsv")
#' @return Path to the merged output file
#' @export
sumstats <- function(files, ref, trait.names = NULL, info.filter = 0.6,
                     maf.filter = 0.01, keep.indel = FALSE,
                     out = "merged_sumstats.tsv") {
  if (is.null(trait.names)) {
    trait.names <- tools::file_path_sans_ext(basename(files))
  }

  json <- .Call("wrap__sumstats_rust",
    as.character(files), as.character(ref), as.character(trait.names),
    as.double(info.filter), as.double(maf.filter), as.logical(keep.indel),
    as.character(out)
  )
  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::sumstats error: ", result$error)
  result$path
}
