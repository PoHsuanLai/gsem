#!/usr/bin/env Rscript
# Run the R GenomicSEM pipeline and save outputs + timings.
# Usage: Rscript r_pipeline.R

library(GenomicSEM)
library(jsonlite)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

ld_path <- "data/eur_w_ld_chr/"
traits <- paste0("data/iter1GWAS", 1:3, ".sumstats.gz")
trait_names <- c("V1", "V2", "V3")
outdir <- "out_r"
dir.create(outdir, showWarnings = FALSE)

timings <- list()

# ── Step 1: LDSC ──────────────────────────────────────────────────────────────
cat("=== R: Running LDSC ===\n")
t0 <- proc.time()

ldsc_result <- ldsc(
  traits = traits,
  sample.prev = c(NA, NA, NA),
  population.prev = c(NA, NA, NA),
  ld = ld_path,
  wld = ld_path,
  trait.names = trait_names,
  n.blocks = 200
)

t1 <- proc.time()
timings$ldsc <- (t1 - t0)["elapsed"]
cat(sprintf("  LDSC time: %.2f s\n", timings$ldsc))

# Save LDSC outputs
S <- as.matrix(ldsc_result$S)
V <- as.matrix(ldsc_result$V)
I_mat <- as.matrix(ldsc_result$I)
N_vec <- as.numeric(ldsc_result$N)
m <- ldsc_result$m

write.csv(S, file.path(outdir, "S.csv"), row.names = FALSE)
write.csv(V, file.path(outdir, "V.csv"), row.names = FALSE)
write.csv(I_mat, file.path(outdir, "I.csv"), row.names = FALSE)

writeLines(toJSON(list(
  S = S,
  V = V,
  I = I_mat,
  N = N_vec,
  m = m
), auto_unbox = TRUE, digits = 15), file.path(outdir, "ldsc_result.json"))

cat(sprintf("  S matrix:\n"))
print(round(S, 4))
cat(sprintf("  Intercepts:\n"))
print(round(I_mat, 4))

# ── Step 2: SEM (common factor) ──────────────────────────────────────────────
cat("\n=== R: Running Common Factor SEM ===\n")
t0 <- proc.time()

sem_result <- tryCatch({
  commonfactor(ldsc_result, estimation = "DWLS")
}, error = function(e) {
  cat("  commonfactor() failed:", e$message, "\n")
  NULL
})

t1 <- proc.time()
timings$sem <- (t1 - t0)["elapsed"]
cat(sprintf("  SEM time: %.2f s\n", timings$sem))

if (!is.null(sem_result)) {
  write.csv(sem_result$results, file.path(outdir, "sem_results.csv"), row.names = FALSE)
  write.csv(sem_result$modelfit, file.path(outdir, "sem_fit.csv"), row.names = FALSE)
  cat("  Model fit:\n")
  print(sem_result$modelfit)
  cat("  Estimates:\n")
  est_cols <- intersect(c("lhs", "op", "rhs", "est", "STD_Genotype", "SE"), colnames(sem_result$results))
  print(sem_result$results[, est_cols])
}

# ── Save timings ──────────────────────────────────────────────────────────────
writeLines(toJSON(timings, auto_unbox = TRUE), file.path(outdir, "timings.json"))
cat(sprintf("\n=== R Pipeline Complete ===\n"))
cat(sprintf("  LDSC: %.2f s\n", timings$ldsc))
cat(sprintf("  SEM:  %.2f s\n", timings$sem))
