#' Merge GWAS Summary Statistics
#'
#' Reads multiple GWAS files, applies QC filters, aligns alleles,
#' and outputs a merged tab-delimited file for multivariate GWAS.
#'
#' @param files Character vector of GWAS file paths
#' @param ref Path to reference panel file (e.g., w_hm3.snplist)
#' @param trait.names Character vector of trait names
#' @param se.logit Logical vector indicating which traits have logistic SEs
#' @param OLS Logical vector indicating which traits are from OLS regression
#' @param linprob Logical vector indicating which traits are linear probability
#' @param N Numeric vector or list of sample size overrides
#' @param betas Named list of beta column overrides per trait (default NULL = auto-detect)
#' @param info.filter INFO score filter (default 0.6)
#' @param maf.filter MAF filter (default 0.01)
#' @param keep.indel Keep indels (default FALSE)
#' @param parallel Use a parallel rayon worker pool to read the
#'   reference and GWAS files (default \code{TRUE}). Each input file is
#'   decompressed and parsed on its own worker thread. Set to
#'   \code{FALSE} to force single-threaded execution.
#' @param cores Integer cap on the rayon pool size. When \code{NULL}
#'   (the default) rayon honours \code{RAYON_NUM_THREADS} if set, else
#'   it uses the number of logical cores reported by the OS. Since the
#'   reads are parallelized across files, values above
#'   \code{length(files) + 1} don't help. On many-core machines (32+)
#'   or when the underlying BLAS is multithreaded, set this explicitly
#'   to avoid oversubscribing CPUs with nested BLAS threads.
#' @param ambig Keep ambiguous SNPs (default FALSE)
#' @param direct.filter Apply MAF filter directly to GWAS file frequencies (default FALSE)
#' @param out Output file path for the merged sumstats TSV (default "merged_sumstats.tsv")
#' @return Path to the merged output file
#' @examples
#' \dontrun{
#' # Merge munged sumstats into a single SNP x trait TSV for GWAS.
#' sumstats(
#'   files = c("T1.sumstats.gz", "T2.sumstats.gz"),
#'   ref = "eur_w_ld_chr/",
#'   trait.names = c("V1", "V2"),
#'   out = "merged_sumstats.tsv"
#' )
#' }
#' @export
sumstats <- function(files, ref, trait.names=NULL, se.logit, OLS=NULL, linprob=NULL,
                     N=NULL, betas=NULL, info.filter=0.6, maf.filter=0.01,
                     keep.indel=FALSE, parallel=TRUE, cores=NULL, ambig=FALSE,
                     direct.filter=FALSE, out="merged_sumstats.tsv") {

  num_threads <- .resolve_num_threads(parallel, cores)

  if (is.null(trait.names)) {
    trait.names <- tools::file_path_sans_ext(basename(files))
  }

  n_traits <- length(files)

  # Handle se.logit: default to FALSE for all traits if missing
  if (missing(se.logit) || is.null(se.logit)) {
    se.logit <- rep(FALSE, n_traits)
  }

  # Handle OLS: default to NULL -> FALSE for all traits
  if (is.null(OLS)) {
    OLS <- rep(FALSE, n_traits)
  }

  # Handle linprob: default to NULL -> FALSE for all traits
  if (is.null(linprob)) {
    linprob <- rep(FALSE, n_traits)
  }

  # Handle N overrides: serialize into the simple JSON-array format the
  # Rust side parses (e.g. `[1.0,null,2.5]`). Kept as a string because this
  # is a low-frequency config argument.
  n_overrides_json <- if (is.null(N)) {
    "[]"
  } else {
    vals <- sapply(N, function(x) {
      if (is.null(x) || (length(x) == 1 && is.na(x))) "null" else format(as.numeric(x), scientific = FALSE)
    })
    paste0("[", paste(vals, collapse = ","), "]")
  }

  # Handle betas: serialize into `{"trait1":"BETA_COL",...}` format.
  betas_json <- if (is.null(betas)) {
    "{}"
  } else {
    entries <- mapply(function(k, v) {
      paste0("\"", k, "\":\"", v, "\"")
    }, names(betas), betas, USE.NAMES = FALSE)
    paste0("{", paste(entries, collapse = ","), "}")
  }

  result <- .Call("wrap__sumstats_rust",
    as.character(files), as.character(ref), as.character(trait.names),
    as.double(info.filter), as.double(maf.filter), as.logical(keep.indel),
    as.character(out),
    as.integer(se.logit), as.integer(OLS), as.integer(linprob),
    as.character(n_overrides_json),
    as.logical(ambig),
    as.character(betas_json),
    as.logical(direct.filter),
    num_threads
  )
  if (!is.null(result$error)) stop("gsemr::sumstats error: ", result$error)
  result$path
}
