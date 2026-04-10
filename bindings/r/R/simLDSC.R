#' Simulate GWAS Summary Statistics
#'
#' Simulates multivariate GWAS summary statistics given a genetic covariance
#' structure. Port of R GenomicSEM's \code{simLDSC()}.
#'
#' @param covmat Population genetic covariance matrix (matrix or data.frame)
#' @param N Sample size matrix or vector of per-trait sample sizes
#' @param seed Random seed (ignored in gsemr -- Rust uses its own RNG)
#' @param ld Path to LD score directory
#' @param rPheno Phenotypic correlation matrix (default NULL = no environmental correlation)
#' @param int LDSC intercept matrix (default NULL = identity)
#' @param N_overlap Sample overlap proportion (default 0.99)
#' @param r Number of replications (ignored in gsemr -- returns 1 replication)
#' @param gzip_output Compress output (ignored in gsemr)
#' @param parallel Accepted for API compatibility; gsemr's simLDSC is currently
#'   single-threaded so this has no effect.
#' @param cores Accepted for API compatibility; see \code{parallel}.
#' @return A matrix of simulated Z-statistics (k traits x n_snps)
#' @export
simLDSC <- function(covmat, N, seed = 1234, ld, rPheno = NULL, int = NULL,
                    N_overlap = 0.99, r = 1, gzip_output = TRUE,
                    parallel = FALSE, cores = NULL) {
  if (r > 1) {
    message("Note: gsemr simLDSC returns 1 replication; 'r' > 1 is ignored")
  }
  if (identical(parallel, TRUE) || (!is.null(cores) && is.numeric(cores) && cores > 1)) {
    message("Note: gsemr simLDSC is single-threaded; 'parallel'/'cores' are no-ops")
  }

  s_mat <- as.matrix(covmat)
  k <- nrow(s_mat)

  # Extract per-trait N from diagonal if matrix, or use directly if vector
  if (is.matrix(N) || is.data.frame(N)) {
    n_per_trait <- diag(as.matrix(N))
  } else {
    n_per_trait <- rep(as.numeric(N), length.out = k)
  }

  # Read LD scores from the directory (use chr 1-22)
  ld_path <- as.character(ld)
  ld_scores <- c()
  m_total <- 0
  for (chr in 1:22) {
    f <- file.path(ld_path, paste0(chr, ".l2.ldscore.gz"))
    if (file.exists(f)) {
      d <- read.table(gzfile(f), header = TRUE)
      ld_scores <- c(ld_scores, d$L2)
    }
    mf <- file.path(ld_path, paste0(chr, ".l2.M_5_50"))
    if (file.exists(mf)) {
      m_total <- m_total + sum(as.numeric(readLines(mf)))
    }
  }

  if (length(ld_scores) == 0) {
    stop("No LD scores found in ", ld_path)
  }
  if (m_total == 0) m_total <- length(ld_scores)

  s_json <- jsonlite::toJSON(s_mat, digits = 15)

  # Convert intercept matrix to JSON
  int_json <- if (is.null(int)) "null" else jsonlite::toJSON(as.matrix(int), digits = 15)

  # Convert phenotypic correlation matrix to JSON
  r_pheno_json <- if (is.null(rPheno)) "null" else jsonlite::toJSON(as.matrix(rPheno), digits = 15)

  # N_overlap: use 0 if rPheno is NULL (no effect without it)
  n_overlap_val <- if (is.null(rPheno)) 0.0 else as.double(N_overlap)

  json <- .Call("wrap__sim_ldsc_rust",
    as.character(s_json),
    as.numeric(n_per_trait),
    as.numeric(ld_scores),
    as.numeric(m_total),
    as.character(int_json),
    as.character(r_pheno_json),
    n_overlap_val
  )

  jsonlite::fromJSON(json)
}
