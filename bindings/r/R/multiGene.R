#' Multi-Gene Joint Analysis
#'
#' Fits a model with multiple genes simultaneously, accounting for LD between them.
#' Port of R GenomicSEM's \code{multiGene()}. Reuses the multi-SNP engine.
#'
#' @param covstruc LDSC result (named list with S, V, I, N, m components)
#' @param Genes Data frame of gene summary statistics (with beta, se, var columns per gene)
#' @param LD LD correlation matrix between genes
#' @param GeneSE Gene SE override (default "F" = auto, numeric = override)
#' @param Genelist Optional gene list to subset (character vector)
#' @return A list with converged, chisq, df, and params
#' @export
multiGene <- function(covstruc, Genes, LD, GeneSE = "F", Genelist = "F") {

  k <- nrow(as.matrix(covstruc$S))
  n_genes <- nrow(as.matrix(LD))
  gene_names <- if (!is.null(rownames(as.matrix(LD)))) {
    rownames(as.matrix(LD))
  } else {
    paste0("Gene", seq_len(n_genes))
  }

  beta_mat <- matrix(as.numeric(as.matrix(Genes$beta)), nrow = n_genes)
  se_mat <- matrix(as.numeric(as.matrix(Genes$se)), nrow = n_genes)
  var_gene <- as.numeric(Genes$var)
  ld_mat <- matrix(as.numeric(as.matrix(LD)), nrow = nrow(as.matrix(LD)))

  model <- paste0("F1 =~ ", paste0("NA*V", seq_len(k), collapse = " + "),
                   "\n", paste0("F1 ~ ", gene_names, collapse = "\n"),
                   "\nF1 ~~ 1*F1")

  # Filter by Genelist if provided
  if (!identical(Genelist, "F") && !is.null(Genelist)) {
    keep <- gene_names %in% Genelist
    if (any(keep)) {
      gene_names <- gene_names[keep]
      beta_mat <- beta_mat[keep, , drop = FALSE]
      se_mat <- se_mat[keep, , drop = FALSE]
      var_gene <- var_gene[keep]
      ld_mat <- ld_mat[keep, keep, drop = FALSE]
    }
  }

  # Convert GeneSE: "F" means auto (pass NaN), numeric means override
  gene_se_val <- if (identical(GeneSE, "F")) NaN else as.double(GeneSE)

  result <- .Call("wrap__multi_gene_rust",
    .covstruc_as_list(covstruc),
    as.character(model),
    "DWLS",
    beta_mat,
    se_mat,
    as.numeric(var_gene),
    ld_mat,
    as.character(gene_names),
    gene_se_val
  )

  if (!is.null(result$error)) stop("gsemr::multiGene error: ", result$error)
  list(
    converged = result$converged,
    chisq = result$chisq,
    df = result$df,
    results = as.data.frame(result$params, stringsAsFactors = FALSE)
  )
}
