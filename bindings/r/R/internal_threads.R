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
