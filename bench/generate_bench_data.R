#!/usr/bin/env Rscript
# Generate realistic 3-trait GWAS data using simLDSC.
# Requires: GenomicSEM, LD scores in data/eur_w_ld_chr/

library(GenomicSEM)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

ld_path <- "data/eur_w_ld_chr/"
if (!dir.exists(ld_path)) stop("LD scores not found. Run setup.sh first.")

# 3-trait model: F1 =~ 0.7*V1 + 0.6*V2 + 0.5*V3
# True genetic covariance:
#   S = Lambda * Psi * Lambda' + Theta
#   Lambda = (0.7, 0.6, 0.5)', Psi = 0.10, Theta = diag(0.11, 0.14, 0.15)
Spop <- matrix(c(
  0.60, 0.42, 0.35,
  0.42, 0.50, 0.30,
  0.35, 0.30, 0.40
), 3, 3, byrow = TRUE)

# Sample sizes: 50k each, 99% overlap
N_mat <- matrix(50000, 3, 3)
diag(N_mat) <- 50000

# Phenotypic correlations
rPheno <- matrix(c(
  1.0, 0.3, 0.2,
  0.3, 1.0, 0.25,
  0.2, 0.25, 1.0
), 3, 3, byrow = TRUE)

# Intercepts (slight inflation from sample overlap)
intercepts <- c(1.02, 1.03, 1.01)

cat("Generating simulated GWAS data (3 traits, N=50k)...\n")

# Save ground truth
saveRDS(Spop, "data/ground_truth_S.rds")
write.csv(Spop, "data/ground_truth_S.csv", row.names = FALSE)

# Generate simulated data
simLDSC(
  covmat = Spop,
  N = N_mat,
  rPheno = rPheno,
  int = intercepts,
  ld = ld_path,
  r = 1,
  seed = 42,
  gzip_output = TRUE
)

# Move output files to data/
for (f in list.files(".", pattern = "iter1\\.GWAS.*\\.sumstats\\.gz")) {
  file.rename(f, file.path("data", f))
}

cat("Done. Files:\n")
cat(paste(" ", list.files("data", pattern = "sumstats")), sep = "\n")
