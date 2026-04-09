#' Enrichment Analysis
#'
#' Tests for annotation enrichment using stratified LDSC results.
#' Port of R GenomicSEM's \code{enrich()}.
#'
#' @param s_covstruc Stratified LDSC result (list with S_baseline, S_annot, V_annot, annotation_names, m_annot, m_total)
#' @param model Model syntax (ignored in gsemr -- enrichment uses basic test)
#' @param params Parameters to test (ignored in gsemr)
#' @param fix Which parameters to fix (ignored in gsemr)
#' @param std.lv Standardize latent variables (ignored in gsemr)
#' @param rm_flank Remove flanking regions (ignored in gsemr)
#' @param tau Use tau parameterization (ignored in gsemr)
#' @param base Include baseline (ignored in gsemr)
#' @param toler Tolerance (ignored in gsemr)
#' @param fixparam Fixed parameters (ignored in gsemr)
#' @return A data frame with annotation, enrichment, se, and p columns
#' @export
enrich <- function(s_covstruc, model = "", params = NULL, fix = "regressions",
                   std.lv = FALSE, rm_flank = TRUE, tau = FALSE, base = TRUE,
                   toler = NULL, fixparam = NULL) {
  # Notify about ignored params
  if (nchar(model) > 0) {
    message("Note: 'model' is accepted for API compatibility but enrichment uses the basic test in gsemr")
  }
  if (!is.null(params)) {
    message("Note: 'params' is ignored in gsemr")
  }
  if (!identical(fix, "regressions")) {
    message("Note: 'fix' is ignored in gsemr")
  }
  if (!identical(std.lv, FALSE)) {
    message("Note: 'std.lv' is ignored in gsemr")
  }
  if (!identical(rm_flank, TRUE)) {
    message("Note: 'rm_flank' is ignored in gsemr")
  }
  if (!identical(tau, FALSE)) {
    message("Note: 'tau' is ignored in gsemr")
  }
  if (!is.null(toler)) {
    message("Note: 'toler' is ignored in gsemr")
  }
  if (!is.null(fixparam)) {
    message("Note: 'fixparam' is ignored in gsemr")
  }

  # Extract components from stratified LDSC result
  s_baseline_json <- jsonlite::toJSON(as.matrix(s_covstruc$S_baseline), digits = 15)
  s_annot_json <- jsonlite::toJSON(lapply(s_covstruc$S_annot, as.matrix), digits = 15)
  v_annot_json <- jsonlite::toJSON(lapply(s_covstruc$V_annot, as.matrix), digits = 15)
  annotation_names <- as.character(s_covstruc$annotation_names)
  m_annot <- as.numeric(s_covstruc$m_annot)
  m_total <- as.numeric(s_covstruc$m_total)

  json <- .Call("wrap__enrich_rust",
    as.character(s_baseline_json),
    as.character(s_annot_json),
    as.character(v_annot_json),
    annotation_names,
    m_annot,
    m_total
  )

  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::enrich error: ", result$error)
  as.data.frame(result)
}
