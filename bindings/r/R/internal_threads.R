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

#' Wrap a two-table GWAS result from the Rust binding into R data frames.
#'
#' The Rust side returns `list(snps = ..., params = ...)` (or `list(genes
#' = ..., params = ...)` for TWAS) where each inner slot is already a
#' columnar list of equal-length vectors. Use `structure(..., class =
#' "data.frame", row.names = .set_row_names(n))` to skip the
#' `data.frame()` constructor's per-column type coercion, which is
#' prohibitive at 20M-row params tables.
#'
#' @param result Two-element list returned by the Rust binding.
#' @return A two-element named list (`snps` + `params`, or `genes` +
#'   `params`), each slot being a `data.frame`.
#' @keywords internal
.snp_columnar_to_df <- function(result) {
  as_df <- function(cols) {
    if (is.null(cols) || length(cols) == 0L) return(data.frame())
    n <- length(cols[[1]])
    structure(
      cols,
      class     = "data.frame",
      names     = names(cols),
      row.names = .set_row_names(n)
    )
  }
  top_name <- if (!is.null(result$snps)) "snps" else "genes"
  setNames(
    list(as_df(result[[top_name]]), as_df(result$params)),
    c(top_name, "params")
  )
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
