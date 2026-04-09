#' Multi-Gene Joint Analysis
#'
#' Fits a model with multiple genes simultaneously, accounting for LD between them.
#' Port of R GenomicSEM's \code{multiGene()}. Reuses the multi-SNP engine.
#'
#' @param covstruc LDSC result (list with S, V, I components, or JSON string)
#' @param Genes Data frame of gene summary statistics (with beta, se, var columns per gene)
#' @param LD LD correlation matrix between genes
#' @param GeneSE Gene SE flag (ignored in gsemr)
#' @param Genelist Optional gene list to subset (character vector)
#' @return A list with converged, chisq, df, and params
#' @export
multiGene <- function(covstruc, Genes, LD, GeneSE = "F", Genelist = "F") {
  if (!identical(GeneSE, "F")) {
    message("Note: 'GeneSE' is ignored in gsemr")
  }
  if (!identical(Genelist, "F")) {
    message("Note: 'Genelist' is ignored in gsemr -- subsetting not yet supported")
  }

  if (is.list(covstruc) && !is.character(covstruc)) {
    covstruc_json <- jsonlite::toJSON(list(
      s = covstruc$S,
      v = covstruc$V,
      i_mat = covstruc$I,
      n_vec = covstruc$N,
      m = covstruc$m
    ), auto_unbox = TRUE)
  } else {
    covstruc_json <- covstruc
  }

  k <- nrow(as.matrix(covstruc$S))
  n_genes <- nrow(as.matrix(LD))
  gene_names <- if (!is.null(rownames(as.matrix(LD)))) {
    rownames(as.matrix(LD))
  } else {
    paste0("Gene", seq_len(n_genes))
  }

  beta_mat <- as.matrix(Genes$beta)
  se_mat <- as.matrix(Genes$se)
  var_gene <- as.numeric(Genes$var)

  model <- paste0("F1 =~ ", paste0("NA*V", seq_len(k), collapse = " + "),
                   "\n", paste0("F1 ~ ", gene_names, collapse = "\n"),
                   "\nF1 ~~ 1*F1")

  beta_json <- jsonlite::toJSON(beta_mat, digits = 15)
  se_json <- jsonlite::toJSON(se_mat, digits = 15)
  ld_json <- jsonlite::toJSON(as.matrix(LD), digits = 15)

  json <- .Call("wrap__multi_gene_rust",
    as.character(covstruc_json),
    as.character(model),
    "DWLS",
    as.character(beta_json),
    as.character(se_json),
    as.numeric(var_gene),
    as.character(ld_json),
    as.character(gene_names)
  )

  result <- jsonlite::fromJSON(json)
  if (!is.null(result$error)) stop("gsemr::multiGene error: ", result$error)
  list(
    converged = result$converged,
    chisq = result$chisq,
    df = result$df,
    results = as.data.frame(result$params)
  )
}
