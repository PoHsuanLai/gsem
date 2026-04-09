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
#' @export
enrich <- function(s_covstruc, model = "", params = NULL, fix = "regressions",
                   std.lv = FALSE, rm_flank = TRUE, tau = FALSE, base = TRUE,
                   toler = NULL, fixparam = NULL) {

  # If no model specified, use the fast Rust proportional enrichment test
  if (!nzchar(model)) {
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
    return(as.data.frame(result))
  }

  # SEM-based enrichment: fit model per annotation and compare to baseline
  n_annot <- length(s_covstruc$S_annot)
  annotation_names <- as.character(s_covstruc$annotation_names)
  m_annot <- as.numeric(s_covstruc$m_annot)

  # Tau parameterization: convert S_annot to per-SNP contribution
  # tau = S_annot / m_annot (removes M scaling, gives per-SNP genetic covariance)
  s_annot_use <- s_covstruc$S_annot
  if (tau) {
    s_annot_use <- lapply(seq_len(n_annot), function(a) {
      if (m_annot[a] > 0) {
        as.matrix(s_covstruc$S_annot[[a]]) / m_annot[a]
      } else {
        as.matrix(s_covstruc$S_annot[[a]])
      }
    })
  }

  # Fit baseline model
  baseline_covstruc <- list(
    S = as.matrix(s_covstruc$S_baseline),
    V = s_covstruc$V_annot[[1]],  # Use first annotation's V as baseline V
    I = diag(nrow(as.matrix(s_covstruc$S_baseline))),
    N = NULL,
    m = s_covstruc$m_total
  )
  baseline_fit <- usermodel(baseline_covstruc, model = model, std.lv = std.lv,
                            toler = toler, CFIcalc = FALSE)
  baseline_params <- baseline_fit$results

  # Apply fixparam if specified
  model_with_fixes <- model
  if (!is.null(fixparam)) {
    for (pname in names(fixparam)) {
      # fixparam adds constraints like "param_name == value"
      model_with_fixes <- paste0(model_with_fixes, "\n", pname, " == ", fixparam[[pname]])
    }
  }

  # Determine which baseline params to fix based on 'fix' argument
  if (identical(fix, "regressions")) {
    fix_rows <- baseline_params[baseline_params$op == "~", ]
  } else if (identical(fix, "loadings")) {
    fix_rows <- baseline_params[baseline_params$op == "=~", ]
  } else {
    fix_rows <- data.frame()
  }

  # Build per-annotation model with fixed parameters from baseline
  annot_model <- model_with_fixes
  if (nrow(fix_rows) > 0) {
    for (i in seq_len(nrow(fix_rows))) {
      constraint <- paste0(fix_rows$lhs[i], " ", fix_rows$op[i], " ",
                           format(fix_rows$est[i], digits = 10), "*", fix_rows$rhs[i])
      annot_model <- paste0(annot_model, "\n", constraint)
    }
  }

  # Fit model per annotation
  results_list <- list()
  for (a in seq_len(n_annot)) {
    annot_covstruc <- list(
      S = as.matrix(s_annot_use[[a]]),
      V = as.matrix(s_covstruc$V_annot[[a]]),
      I = diag(nrow(as.matrix(s_annot_use[[a]]))),
      N = NULL,
      m = if (tau) 1.0 else s_covstruc$m_annot[a]
    )

    annot_fit <- tryCatch(
      usermodel(annot_covstruc, model = annot_model, std.lv = std.lv,
                toler = toler, CFIcalc = FALSE),
      error = function(e) NULL
    )

    if (is.null(annot_fit)) {
      results_list[[a]] <- data.frame(
        annotation = annotation_names[a],
        stringsAsFactors = FALSE
      )
      next
    }

    annot_params <- annot_fit$results

    # Filter to requested params
    if (!is.null(params)) {
      param_keys <- paste0(annot_params$lhs, annot_params$op, annot_params$rhs)
      annot_params <- annot_params[param_keys %in% params, ]
    }

    if (nrow(annot_params) > 0) {
      annot_params$annotation <- annotation_names[a]
      results_list[[a]] <- annot_params
    } else {
      results_list[[a]] <- data.frame(
        annotation = annotation_names[a],
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results_list)
}
