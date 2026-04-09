#' Munge GWAS Summary Statistics
#'
#' Filters and harmonizes GWAS summary statistics for use with
#' \code{\link{ldsc}}.
#'
#' @param files Character vector of GWAS file paths
#' @param hm3 Path to HapMap3 SNP reference file
#' @param trait.names Character vector of trait names
#' @param info.filter INFO score filter threshold (default 0.9)
#' @param maf.filter MAF filter threshold (default 0.01)
#' @param out Output directory for .sumstats.gz files
#' @return Character vector of output file paths (invisibly)
#' @export
munge <- function(files,
                  hm3,
                  trait.names = paste0("trait", seq_along(files)),
                  info.filter = 0.9,
                  maf.filter = 0.01,
                  out = ".") {

  dir.create(out, showWarnings = FALSE, recursive = TRUE)

  result <- .Call("wrap__munge_rust",
    as.character(files),
    as.character(hm3),
    as.character(trait.names),
    as.double(info.filter),
    as.double(maf.filter),
    as.character(out)
  )

  if (length(result) == 0) {
    warning("munge produced no output files")
  }

  invisible(result)
}
