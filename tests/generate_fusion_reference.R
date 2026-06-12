#!/usr/bin/env Rscript
#
# Reference fixture for read_fusion (R GenomicSEM's FUSION/TWAS .dat reader).
#
# Writes a pair of small FUSION .dat association files (committed under
# tests/fixtures/fusion/) and the merged output of the REAL
# GenomicSEM::read_fusion run on them (binary + continuous traits, NA dropping,
# inner join across traits). The Rust port (gsem::io::fusion_reader) must
# reproduce these numbers exactly.
#
# Usage: cd tests && Rscript generate_fusion_reference.R

suppressMessages(library(jsonlite))
suppressMessages(library(GenomicSEM))

outdir  <- "fixtures"
datadir <- file.path(outdir, "fusion")
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)

# Two FUSION .dat files. The FUSION format is whitespace-delimited with (among
# others) FILE (weight path), ID (gene), TWAS.Z, HSQ. G3 has a missing HSQ (".")
# so it is dropped by R's na.omit; trait 2 omits G3 entirely, so the inner join
# keeps {G1, G2}.
t1 <- c(
  "PANEL FILE ID TWAS.Z HSQ BEST.GWAS.ID",
  "x /path/to/GTEx_Brain/ENSG001.wgt.RDat G1 2.50 0.30 rs1",
  "x /path/to/GTEx_Brain/ENSG002.wgt.RDat G2 -1.20 0.18 rs2",
  "x /path/to/GTEx_Brain/ENSG003.wgt.RDat G3 0.80 . rs3"
)
t2 <- c(
  "PANEL FILE ID TWAS.Z HSQ BEST.GWAS.ID",
  "x /path/to/GTEx_Brain/ENSG001.wgt.RDat G1 1.10 0.30 rs1",
  "x /path/to/GTEx_Brain/ENSG002.wgt.RDat G2 3.40 0.18 rs2"
)
writeLines(t1, file.path(datadir, "trait1.dat"))
writeLines(t2, file.path(datadir, "trait2.dat"))

trait_names <- c("A", "B")
binary <- c(TRUE, FALSE)   # A: binary (liability conversion), B: continuous
n_vec  <- c(10000, 8000)

cat("=== read_fusion ===\n")
out <- GenomicSEM::read_fusion(
  files       = c(file.path(datadir, "trait1.dat"), file.path(datadir, "trait2.dat")),
  trait.names = trait_names,
  binary      = binary,
  N           = n_vec,
  perm        = FALSE
)
print(out)

fixture <- list(
  files       = c("fusion/trait1.dat", "fusion/trait2.dat"),
  trait_names = trait_names,
  binary      = binary,
  n           = n_vec,
  gene        = as.character(out$Gene),
  panel       = as.character(out$Panel),
  hsq         = as.numeric(out$HSQ),
  beta        = lapply(seq_len(nrow(out)), function(i)
                  as.numeric(c(out[i, "beta.A"], out[i, "beta.B"]))),
  se          = lapply(seq_len(nrow(out)), function(i)
                  as.numeric(c(out[i, "se.A"], out[i, "se.B"])))
)
writeLines(toJSON(fixture, auto_unbox = TRUE, digits = 15),
           file.path(outdir, "fusion_read.json"))
cat("Wrote fixtures/fusion_read.json\n")
