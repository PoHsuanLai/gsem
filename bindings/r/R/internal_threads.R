#' Normalize a `covstruc` argument into the exact list shape the Rust
#' binding expects.
#'
#' The Rust side accepts fields `s`, `v`, `i_mat`, `n_vec`, `m` (or their
#' capitalized `S`, `V`, `I`, `N`, `m` aliases). This helper passes the
#' output of \code{\link{ldsc}} through unchanged and converts any data
#' frames in the S/V/I slots into matrices.
#'
#' @param covstruc A named list from `ldsc()` (or a list with the same
#'   field names).
#' @return A named list ready to pass as a `covstruc` argument in a
#'   `.Call()` invocation.
#' @keywords internal
.covstruc_as_list <- function(covstruc) {
  if (!is.list(covstruc)) {
    stop("covstruc must be a list (output of ldsc())")
  }
  s <- covstruc$S; if (is.null(s)) s <- covstruc$s
  v <- covstruc$V; if (is.null(v)) v <- covstruc$v
  i <- covstruc$I; if (is.null(i)) i <- covstruc$i_mat
  n <- covstruc$N; if (is.null(n)) n <- covstruc$n_vec
  m <- covstruc$m; if (is.null(m)) m <- covstruc$M

  if (is.null(s)) stop("covstruc is missing 'S'")
  if (is.null(v)) stop("covstruc is missing 'V'")
  if (is.null(i)) stop("covstruc is missing 'I'")

  list(
    s = as.matrix(s),
    v = as.matrix(v),
    i_mat = as.matrix(i),
    n_vec = if (is.null(n)) numeric(0) else as.numeric(n),
    m = if (is.null(m)) 0 else as.numeric(m)
  )
}

#' Convert the columnar list returned by `snp_results_to_list` /
#' `twas_results_to_list` into the per-SNP data frame shape that matches
#' what `jsonlite::fromJSON()` previously produced: one top-level row per
#' SNP/gene with a list-column `params` holding a per-SNP data frame of
#' parameter estimates.
#'
#' The Rust side returns:
#'   list(SNP = chr[n], chisq = dbl[n], df = int[n], converged = lgl[n],
#'        [Q_chisq = dbl[n], Q_df = int[n], Q_pval = dbl[n],]
#'        params = list(snp_idx = int[m], lhs = chr[m], op = chr[m],
#'                      rhs = chr[m], est = dbl[m], se = dbl[m],
#'                      z = dbl[m], p = dbl[m]))
#'
#' where `snp_idx` is a 1-based index into the SNP columns. This helper
#' splits the flat params columns by `snp_idx` so each row gets its own
#' data frame under `$params[[i]]`, matching R GenomicSEM's historical
#' output and what downstream wrapper code (e.g. commonfactorGWAS.R)
#' expects.
#'
#' @keywords internal
.snp_columnar_to_df <- function(result, id_col = "SNP", extra_cols = character()) {
  n <- length(result[[id_col]])

  # Build the per-SNP params data frames by splitting on snp_idx.
  params_flat <- result$params
  per_snp <- vector("list", n)
  if (!is.null(params_flat) && length(params_flat$snp_idx) > 0) {
    for (i in seq_len(n)) {
      mask <- params_flat$snp_idx == i
      per_snp[[i]] <- data.frame(
        lhs = params_flat$lhs[mask],
        op = params_flat$op[mask],
        rhs = params_flat$rhs[mask],
        est = params_flat$est[mask],
        se = params_flat$se[mask],
        z = params_flat$z[mask],
        p = params_flat$p[mask],
        stringsAsFactors = FALSE
      )
    }
  } else {
    for (i in seq_len(n)) {
      per_snp[[i]] <- data.frame(
        lhs = character(0), op = character(0), rhs = character(0),
        est = numeric(0), se = numeric(0), z = numeric(0), p = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  }

  df <- data.frame(
    setNames(list(result[[id_col]]), id_col),
    stringsAsFactors = FALSE
  )
  for (col in extra_cols) {
    df[[col]] <- result[[col]]
  }
  df$chisq <- result$chisq
  df$df <- result$df
  df$converged <- result$converged
  if (!is.null(result$Q_chisq)) {
    df$Q_chisq <- result$Q_chisq
    df$Q_df <- result$Q_df
    df$Q_pval <- result$Q_pval
  }
  df$params <- per_snp

  df
}

#' Resolve `parallel`/`cores` arguments to an explicit thread count.
#'
#' Returns an integer suitable for passing to a Rust `.Call` that takes a
#' `num_threads` parameter. The Rust side treats `NA_integer_` as
#' "use rayon default" (= all available cores).
#'
#' Resolution order:
#'   1. If `cores` is a positive integer, return that.
#'   2. Else if `parallel` is `FALSE`, return 1L.
#'   3. Else return `NA_integer_` (let rayon pick).
#'
#' This is an internal helper. Not exported.
#'
#' @param parallel logical
#' @param cores integer or NULL
#' @return integer (possibly `NA_integer_`)
#' @keywords internal
.resolve_num_threads <- function(parallel, cores) {
  if (!is.null(cores) && is.numeric(cores) && length(cores) == 1L && cores > 0) {
    return(as.integer(cores))
  }
  if (identical(parallel, FALSE)) {
    return(1L)
  }
  NA_integer_
}
