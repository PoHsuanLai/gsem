#!/usr/bin/env Rscript
#
# Generate GWAS test fixtures by running R GenomicSEM's commonfactorGWAS
# and userGWAS on real PGC data, saving per-SNP results as JSON.
#
# Prerequisites:
#   - bench/data/ must have GWAS files, LD scores, and reference.1000G.maf.0.005.txt.gz
#   - GenomicSEM and jsonlite packages installed
#
# Usage: cd tests && Rscript generate_gwas_fixture.R

library(jsonlite)
library(GenomicSEM)

outdir <- "fixtures"
dir.create(outdir, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Setup: use the bench data
# ---------------------------------------------------------------------------
bench_dir <- file.path("..", "bench")
setwd(bench_dir)

ld          <- "data/eur_w_ld_chr/"
ref_1000g   <- "data/reference.1000G.maf.0.005.txt.gz"
trait_names <- c("ANX", "OCD", "PTSD")
sample_prev <- c(0.5, 0.5, 0.5)
pop_prev    <- c(0.16, 0.02, 0.07)
n_overrides <- c(NA, 9725, 5831)
munged      <- paste0(trait_names, ".sumstats.gz")

ptsd_file <- Sys.glob("data/*PTSD*")[1]
raw_traits <- c(
  "data/anxiety.meta.full.cc.tbl.gz",
  "data/ocd_aug2017.gz",
  ptsd_file
)

if (!file.exists(ref_1000g)) {
  stop("Missing reference file: ", ref_1000g)
}

# Ensure munged files exist
if (!all(file.exists(munged))) {
  cat("Munging...\n")
  GenomicSEM::munge(
    files = raw_traits, hm3 = "data/w_hm3.snplist",
    trait.names = trait_names, N = n_overrides,
    info.filter = 0.9, maf.filter = 0.01
  )
}

# Run LDSC
cat("Running LDSC...\n")
r_cov <- GenomicSEM::ldsc(
  munged, sample_prev, pop_prev, ld, ld,
  trait.names = trait_names, n.blocks = 200
)

# Run sumstats with 1000G reference
cat("Running sumstats...\n")
r_ss <- GenomicSEM::sumstats(
  files = raw_traits, ref = ref_1000g, trait.names = trait_names,
  se.logit = c(FALSE, FALSE, FALSE),
  OLS = c(FALSE, FALSE, FALSE),
  linprob = c(FALSE, FALSE, FALSE),
  info.filter = 0.6, maf.filter = 0.01
)

cat("Total SNPs from sumstats:", nrow(r_ss), "\n")

# Use a small subset (20 SNPs) for the fixture
N_SNPS <- 20
snps_subset <- head(r_ss, N_SNPS)
cat("Using", nrow(snps_subset), "SNPs for GWAS fixture\n")

# ---------------------------------------------------------------------------
# commonfactorGWAS
# ---------------------------------------------------------------------------
cat("Running R commonfactorGWAS...\n")
cf_result <- tryCatch({
  GenomicSEM::commonfactorGWAS(
    covstruc = r_cov,
    SNPs = snps_subset,
    parallel = FALSE
  )
}, error = function(e) {
  cat("commonfactorGWAS error:", e$message, "\n")
  NULL
})

if (!is.null(cf_result)) {
  cat("commonfactorGWAS succeeded:", nrow(cf_result), "SNPs\n")
  cat("Columns:", colnames(cf_result), "\n")
  cat("Sample result:\n")
  print(head(cf_result, 3))
} else {
  stop("commonfactorGWAS failed — cannot generate fixture")
}

# ---------------------------------------------------------------------------
# userGWAS
# ---------------------------------------------------------------------------
cat("Running R userGWAS...\n")
user_model <- paste0(
  "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
  "F1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP"
)

user_result <- tryCatch({
  GenomicSEM::userGWAS(
    covstruc = r_cov,
    SNPs = snps_subset,
    model = user_model,
    parallel = FALSE
  )
}, error = function(e) {
  cat("userGWAS error:", e$message, "\n")
  NULL
})

if (!is.null(user_result)) {
  # userGWAS returns a list of data frames — one per SNP, each with all params
  if (is.list(user_result) && !is.data.frame(user_result)) {
    user_df <- do.call(rbind, user_result)
  } else {
    user_df <- user_result
  }
  cat("userGWAS raw result class:", class(user_result), "\n")
  cat("userGWAS raw result length:", length(user_result), "\n")
  if (is.list(user_result) && !is.data.frame(user_result)) {
    for (k in seq_along(user_result)) {
      cat(sprintf("  sub-model %d: %d rows, %d unique SNPs\n", k, nrow(user_result[[k]]), length(unique(user_result[[k]]$SNP))))
    }
  }
  cat("user_df rows:", nrow(user_df), "\n")
  cat("Columns:", colnames(user_df), "\n")
  cat("Unique SNPs:", length(unique(user_df$SNP)), "\n")
  cat("Sample result:\n")
  print(head(user_df, 3))
} else {
  stop("userGWAS failed — cannot generate fixture")
}

# ---------------------------------------------------------------------------
# Save fixtures
# ---------------------------------------------------------------------------
setwd(file.path("..", "tests"))

# Save the input data (S, V, I matrices and SNP subset) so Rust can reproduce
fixture <- list(
  # LDSC output
  s = lapply(1:nrow(r_cov$S), function(i) as.numeric(r_cov$S[i,])),
  v = lapply(1:nrow(r_cov$V), function(i) as.numeric(r_cov$V[i,])),
  i_mat = lapply(1:nrow(r_cov$I), function(i) as.numeric(r_cov$I[i,])),
  trait_names = trait_names,

  # SNP input data
  snps = lapply(1:nrow(snps_subset), function(i) {
    row <- snps_subset[i,]
    beta_cols <- grep("^beta\\.", colnames(snps_subset))
    se_cols <- grep("^se\\.", colnames(snps_subset))
    list(
      SNP = as.character(row$SNP),
      A1 = as.character(row$A1),
      A2 = as.character(row$A2),
      MAF = as.numeric(row$MAF),
      beta = as.numeric(row[, beta_cols]),
      se = as.numeric(row[, se_cols])
    )
  }),

  # commonfactorGWAS results — extract scalar values carefully
  commonfactor_gwas = lapply(1:nrow(cf_result), function(i) {
    row <- cf_result[i,]
    # chisq column name varies: "Q" in newer versions
    chisq_val <- NA
    if ("Q" %in% colnames(cf_result)) chisq_val <- as.numeric(row$Q)
    if ("chisq" %in% colnames(cf_result)) chisq_val <- as.numeric(row$chisq)
    df_val <- NA
    if ("Q_df" %in% colnames(cf_result)) df_val <- as.numeric(row$Q_df)
    if ("chisq_df" %in% colnames(cf_result)) df_val <- as.numeric(row$chisq_df)
    if ("df" %in% colnames(cf_result)) df_val <- as.numeric(row$df)
    z_val <- NA
    if ("Z_Estimate" %in% colnames(cf_result)) z_val <- as.numeric(row$Z_Estimate)
    if ("z" %in% colnames(cf_result)) z_val <- as.numeric(row$z)
    p_val <- NA
    if ("Pval_Estimate" %in% colnames(cf_result)) p_val <- as.numeric(row$Pval_Estimate)
    if ("p" %in% colnames(cf_result)) p_val <- as.numeric(row$p)

    list(
      SNP = as.character(row$SNP),
      est = as.numeric(row$est),
      se = as.numeric(row$se_c),
      z = z_val,
      p = p_val,
      chisq = chisq_val,
      chisq_df = df_val
    )
  }),

  # userGWAS results — extract the SNP effect row (F1 ~ SNP) per unique SNP
  user_gwas = {
    unique_snps <- unique(user_df$SNP)
    lapply(unique_snps, function(snp_id) {
      snp_data <- user_df[user_df$SNP == snp_id, ]
      # Find the F1 ~ SNP row (the free parameter)
      snp_effect <- snp_data[snp_data$op == "~" & snp_data$rhs == "SNP", ]
      if (nrow(snp_effect) == 0) {
        snp_effect <- snp_data[snp_data$free > 0, ]
      }
      if (nrow(snp_effect) == 0) return(NULL)
      row <- snp_effect[1, ]

      chisq_val <- if ("chisq" %in% colnames(row)) as.numeric(row$chisq) else NA
      df_val <- if ("chisq_df" %in% colnames(row)) as.numeric(row$chisq_df) else NA
      z_val <- if ("Z_Estimate" %in% colnames(row)) as.numeric(row$Z_Estimate) else NA
      p_val <- if ("Pval_Estimate" %in% colnames(row)) as.numeric(row$Pval_Estimate) else NA
      se_val <- if ("SE" %in% colnames(row)) as.numeric(row$SE) else NA

      list(
        SNP = as.character(snp_id),
        lhs = as.character(row$lhs),
        op = as.character(row$op),
        rhs = as.character(row$rhs),
        est = as.numeric(row$est),
        se = se_val,
        z = z_val,
        p = p_val,
        chisq = chisq_val,
        chisq_df = df_val
      )
    })
  }
)

path <- file.path(outdir, "gwas_per_snp.json")
writeLines(toJSON(fixture, auto_unbox = TRUE, digits = 15), path)
cat("\nWrote", path, "\n")
cat("Done.\n")
