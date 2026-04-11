#' Simulate GWAS Summary Statistics
#'
#' Simulates multivariate GWAS summary statistics given a genetic
#' covariance structure. Drop-in replacement for
#' \code{GenomicSEM::simLDSC()}: writes the same per-trait sumstats
#' files (\code{iter<r>GWAS<i>.sumstats[.gz]}) and \code{PopCovMat.RData}
#' to the working directory that R's upstream produces, and prints the
#' same "Writing sumstats for Phenotype i" banner. The simulated
#' \code{k x n_snps} Z-score matrix is returned invisibly as a gsemr
#' extension — assign the call to capture it, or ignore it to match R's
#' \code{NULL}-returning behaviour exactly.
#'
#' @param covmat Population genetic covariance matrix (matrix or data.frame)
#' @param N Sample size: either a scalar (applied to all traits via
#'   \code{diag(N, k, k)}), a length-k vector (per-trait), or a k x k
#'   sample-size matrix with sample overlap on the off-diagonal.
#'   Matches \code{GenomicSEM::simLDSC()}, which also accepts scalar or
#'   matrix N.
#' @param seed Random seed (ignored in gsemr -- Rust uses its own RNG)
#' @param ld Path to LD score directory
#' @param rPheno Phenotypic correlation matrix (default NULL = no environmental correlation)
#' @param int LDSC intercept matrix (default NULL = identity)
#' @param N_overlap Sample overlap proportion (default 0.99)
#' @param r Replication index — feeds into the output filename pattern
#'   \code{iter<r>GWAS<i>.sumstats}. gsemr returns a single replication;
#'   values \code{r > 1} are accepted for API compatibility but only
#'   the single replication is written (under the given \code{r}).
#' @param gzip_output Compress output with gzip (default TRUE)
#' @param parallel Accepted for API compatibility; gsemr's simLDSC is currently
#'   single-threaded so this has no effect.
#' @param cores Accepted for API compatibility; see \code{parallel}.
#' @return Invisibly, a matrix of simulated Z-statistics with shape
#'   \code{k traits x n_snps}. Side-effects: writes
#'   \code{iter<r>GWAS<i>.sumstats[.gz]} for each trait and
#'   \code{PopCovMat.RData} to the current working directory, matching
#'   R upstream.
#' @examples
#' \dontrun{
#' # Simulate per-trait sumstats from a target genetic covariance matrix.
#' covmat <- matrix(c(0.60, 0.42, 0.42, 0.50), 2, 2,
#'                  dimnames = list(c("T1","T2"), c("T1","T2")))
#' sim <- simLDSC(
#'   covmat = covmat,
#'   N = 1e5,
#'   ld = "eur_w_ld_chr/"
#' )
#' # Writes iter1GWAS1.sumstats.gz, iter1GWAS2.sumstats.gz, PopCovMat.RData.
#' }
#' @export
simLDSC <- function(covmat, N, seed = 1234, ld, rPheno = NULL, int = NULL,
                    N_overlap = 0.99, r = 1, gzip_output = TRUE,
                    parallel = FALSE, cores = NULL) {
  if (r > 1) {
    message("Note: gsemr simLDSC returns 1 replication; 'r' > 1 is accepted but only one file set is written")
  }
  if (identical(parallel, TRUE) || (!is.null(cores) && is.numeric(cores) && cores > 1)) {
    message("Note: gsemr simLDSC is single-threaded; 'parallel'/'cores' are no-ops")
  }

  s_mat <- as.matrix(covmat)
  k <- nrow(s_mat)

  # R GenomicSEM::simLDSC builds a K x K N matrix via `diag(N, k, k)`
  # when N is scalar, or uses the user-supplied matrix as-is, or
  # accepts a vector. Normalise to a per-trait vector for the Rust
  # core, and keep the full matrix around for PopCovMat.RData output.
  n_matrix <- if (is.matrix(N) || is.data.frame(N)) {
    as.matrix(N)
  } else if (length(N) > 1) {
    m <- base::diag(as.numeric(N), k, k)
    # Mirror R's off-diagonal overlap handling
    m[lower.tri(m, diag = FALSE)] <- m[1, 1] * N_overlap
    m[upper.tri(m, diag = FALSE)] <- m[1, 1] * N_overlap
    m
  } else {
    m <- base::diag(as.numeric(N), k, k)
    m[lower.tri(m, diag = FALSE)] <- m[1, 1] * N_overlap
    m[upper.tri(m, diag = FALSE)] <- m[1, 1] * N_overlap
    m
  }
  n_per_trait <- base::diag(n_matrix)

  # Load the LD score files once, keeping the full SNP metadata so we
  # can reconstruct R's per-trait sumstats file layout (SNP, CHR, BP,
  # ... Z, N, A1).
  ld_path <- as.character(ld)
  ld_df_list <- list()
  m_total <- 0
  for (chr in 1:22) {
    f <- file.path(ld_path, paste0(chr, ".l2.ldscore.gz"))
    if (file.exists(f)) {
      d <- read.table(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
      ld_df_list[[length(ld_df_list) + 1L]] <- d
    }
    mf <- file.path(ld_path, paste0(chr, ".l2.M_5_50"))
    if (file.exists(mf)) {
      m_total <- m_total + sum(as.numeric(readLines(mf)))
    }
  }
  if (length(ld_df_list) == 0L) {
    stop("No LD scores found in ", ld_path)
  }
  ld_df <- do.call(rbind, ld_df_list)
  if (m_total == 0) m_total <- nrow(ld_df)
  ld_scores <- ld_df$L2

  s_mat_num <- matrix(as.numeric(s_mat), nrow = nrow(s_mat))

  int_mat <- if (is.null(int)) {
    NULL
  } else {
    matrix(as.numeric(as.matrix(int)), nrow = nrow(as.matrix(int)))
  }
  r_pheno_mat <- if (is.null(rPheno)) {
    NULL
  } else {
    matrix(as.numeric(as.matrix(rPheno)), nrow = nrow(as.matrix(rPheno)))
  }

  n_overlap_val <- if (is.null(rPheno)) 0.0 else as.double(N_overlap)

  # Rust core returns a k x n_snps matrix of Z-statistics.
  z_mat <- .Call("wrap__sim_ldsc_rust",
    s_mat_num,
    as.numeric(n_per_trait),
    as.numeric(ld_scores),
    as.numeric(m_total),
    int_mat,
    r_pheno_mat,
    n_overlap_val
  )

  if (!is.matrix(z_mat) || nrow(z_mat) != k) {
    stop("gsemr::simLDSC: Rust core returned unexpected shape")
  }

  # Write the R-compatible per-trait sumstats files. Column layout
  # mirrors R upstream: start from the LD score data frame (SNP, CHR,
  # BP, MAF, L2, ...), add Z / N / A1. Phenotype names come from the
  # dimnames of covmat if present, otherwise Pheno_1..Pheno_k.
  pheno_names <- rownames(s_mat)
  if (is.null(pheno_names) || length(pheno_names) != k) {
    pheno_names <- paste0("Pheno_", seq_len(k))
  }

  cat(strrep("-", 20), "\n", sep = "")
  cat("\n", "Writing simulated sumstats", "\n", "\n", sep = "")
  for (i in seq_len(k)) {
    gwas_df <- ld_df
    # Guard against LD frames with more SNPs than we simulated (rare
    # but possible if the core trimmed).
    n_sim <- ncol(z_mat)
    if (nrow(gwas_df) != n_sim) {
      gwas_df <- gwas_df[seq_len(min(nrow(gwas_df), n_sim)), , drop = FALSE]
    }
    gwas_df$Z  <- z_mat[i, seq_len(nrow(gwas_df))]
    gwas_df$N  <- n_per_trait[i]
    gwas_df$A1 <- "A"

    cat("Writing sumstats for Phenotype", i, "\n", "\n", sep = "")
    out_path <- paste0("iter", r, "GWAS", i, ".sumstats")
    utils::write.table(gwas_df, file = out_path, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    if (isTRUE(gzip_output)) {
      # Use a pure-R gzip to avoid adding R.utils as a dependency.
      in_con  <- file(out_path, "rb")
      raw     <- readBin(in_con, what = "raw", n = file.info(out_path)$size)
      close(in_con)
      gz_con  <- gzfile(paste0(out_path, ".gz"), "wb")
      writeBin(raw, gz_con)
      close(gz_con)
      file.remove(out_path)
    }
  }
  cat(strrep("-", 37), "\n", sep = "")

  # Save the population covariance matrix R writes alongside the
  # per-trait files. Named `covMatrix` to match R upstream.
  covMatrix <- as.data.frame(s_mat)
  rownames(covMatrix) <- pheno_names
  colnames(covMatrix) <- pheno_names
  save(covMatrix, file = "PopCovMat.RData")

  invisible(z_mat)
}
