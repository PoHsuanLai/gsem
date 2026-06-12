#' Enrichment Analysis
#'
#' Tests for annotation enrichment using stratified LDSC results.
#'
#' When \code{model} is provided, fits a SEM per annotation using
#' \code{\link{usermodel}} and tests for parameter differences vs baseline.
#' When \code{model} is empty, uses the fast Rust proportional enrichment test.
#'
#' @param s_covstruc Stratified LDSC result (list with S_baseline, S_annot, V_annot, annotation_names, m_annot, m_total)
#' @param model lavaan-style model syntax (default "" = use basic proportional enrichment test)
#' @param params Character vector of parameter names to test (default NULL = all free params)
#' @param fix Which parameters to fix: "regressions" (default), "loadings", or "none"
#' @param std.lv Standardize latent variables (default FALSE)
#' @param rm_flank Remove flanking regions (default TRUE; implemented in Rust s_ldsc engine)
#' @param tau Use tau parameterization (default FALSE; not yet used)
#' @param base Include baseline (default TRUE)
#' @param toler Gradient tolerance for optimizer (default NULL = auto)
#' @param fixparam Named list of parameters to fix at specific values (default NULL)
#' @return A data frame with enrichment results
#' @examples
#' \dontrun{
#' # Fast proportional enrichment test (default, model = "").
#' s_covstruc <- s_ldsc(
#'   traits = c("T1.sumstats.gz", "T2.sumstats.gz"),
#'   ld = "baseline_LD/", wld = "weights/", frq = "frq/"
#' )
#' result <- enrich(s_covstruc)
#' head(result)
#' }
#' @export
enrich <- function(s_covstruc, model = "", params = NULL, fix = "regressions",
                   std.lv = FALSE, rm_flank = TRUE, tau = FALSE, base = TRUE,
                   toler = NULL, fixparam = NULL) {

  # If no model specified, use the fast Rust proportional enrichment test
  if (!nzchar(model)) {
    as_num_matrix <- function(M) {
      M <- as.matrix(M)
      matrix(as.numeric(M), nrow = nrow(M))
    }
    s_baseline <- as_num_matrix(s_covstruc$S_baseline)
    s_annot_mats <- lapply(s_covstruc$S_annot, as_num_matrix)
    v_annot_mats <- lapply(s_covstruc$V_annot, as_num_matrix)
    annotation_names <- as.character(s_covstruc$annotation_names)
    m_annot <- as.numeric(s_covstruc$m_annot)
    m_total <- as.numeric(s_covstruc$m_total)

    result <- .Call("wrap__enrich_rust",
      s_baseline,
      s_annot_mats,
      v_annot_mats,
      annotation_names,
      m_annot,
      m_total
    )

    if (!is.null(result$error)) stop("gsemr::enrich error: ", result$error)
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }


  # SEM-based enrichment via the validated Rust implementation
  # (gsem_sem::enrich_model::model_enrichment), matching R GenomicSEM's
  # `enrich`: fit the model to the baseline annotation, fix the
  # regressions/loadings (per `fix`), re-fit each annotation freeing the
  # remaining parameters, and report enrichment = (est_annot/est_base)/Prop.
  as_num_matrix <- function(M) { M <- as.matrix(M); matrix(as.numeric(M), nrow = nrow(M)) }
  s_list <- lapply(s_covstruc$S_annot, as_num_matrix)
  v_list <- lapply(s_covstruc$V_annot, as_num_matrix)
  m_annot <- as.numeric(s_covstruc$m_annot)
  # Proportion of SNPs per annotation relative to the baseline (base = 1).
  prop <- m_annot / m_annot[1]
  annot_names <- as.character(s_covstruc$annotation_names)
  obs_names <- colnames(as.matrix(s_covstruc$S_annot[[1]]))
  if (is.null(obs_names)) obs_names <- paste0("V", seq_len(nrow(s_list[[1]])))
  if (is.null(params)) {
    stop("gsemr::enrich: 'params' must specify the target parameter(s), e.g. \"F1~~F1\".")
  }

  result <- .Call("wrap__enrich_model_rust",
    s_list, v_list, as.numeric(prop),
    as.character(obs_names), annot_names,
    as.character(model), as.character(params), as.character(fix)
  )
  if (!is.null(result$error)) stop("gsemr::enrich error: ", result$error)
  as.data.frame(result, stringsAsFactors = FALSE)
}
