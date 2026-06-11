#!/usr/bin/env Rscript
#
# Reference fixture for enrich (R GenomicSEM's model-based functional
# enrichment). The canonical use is enrichment of a FACTOR VARIANCE: fit the
# measurement model to a baseline annotation, fix the loadings, then re-fit
# per annotation freeing the factor variance, and report
#   enrichment = (est_annot / est_baseline) / Prop  (null = 1).
#
# A clean, well-conditioned 4-trait covstruc (1-factor model over-identified;
# S PSD / V PD so no smoothing) makes the comparison numeric and stable.
#
# Requires: GenomicSEM, jsonlite. Usage: cd tests && Rscript generate_enrich_reference.R

suppressMessages({ library(jsonlite); library(GenomicSEM) })
set.seed(99)
outdir <- "fixtures"
mat_to_list <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))

k <- 4
traits <- paste0("T", 1:k)
L <- c(1.0, 0.8, 0.9, 0.7)              # loadings (marker variable T1 = 1)
resid <- c(0.30, 0.40, 0.35, 0.45)
mkS <- function(phi) {                  # scale the COMMON (factor-variance) part
  S <- outer(L, L) * phi
  diag(S) <- diag(S) + resid
  colnames(S) <- rownames(S) <- traits
  S
}
annot_names <- c("base", "A", "B")
S_list <- list(mkS(1.0), mkS(1.5), mkS(0.7))   # phi: 1.0, 1.5, 0.7
names(S_list) <- annot_names
kstar <- k * (k + 1) / 2                # 10
V_list <- lapply(seq_along(S_list), function(i) diag(kstar) * 1e-4)
names(V_list) <- annot_names
Prop <- data.frame(Prop = c(1.0, 0.30, 0.20))

covstruc <- list(
  S = S_list, V = V_list, S_Tau = S_list, V_Tau = V_list,
  I = diag(k), N = matrix(50000, 1, length(S_list)),
  m = matrix(c(1400, 420, 280), ncol = 1),
  Prop = Prop,
  Select = data.frame(V1 = c("A", "B"), V2 = c(1, 1))  # non-base annotations
)
rownames(covstruc$m) <- annot_names

model <- "F1 =~ T1 + T2 + T3 + T4\nF1 ~~ F1"
params <- c("F1~~F1")

res <- GenomicSEM::enrich(covstruc, model = model, params = params,
                          fix = "regressions", base = TRUE, tau = FALSE,
                          toler = 1e-60)
df <- res[[1]]
print(df[, c("Annotation", "lhs", "op", "rhs", "Enrichment",
             "Enrichment_SE", "Enrichment_p_value")])

out <- list(
  traits = traits, model = model, params = as.list(params),
  annot_names = annot_names,
  s = unname(lapply(S_list, function(M) mat_to_list(as.matrix(M)))),
  v = unname(lapply(V_list, function(M) mat_to_list(as.matrix(M)))),
  prop = as.numeric(Prop$Prop),
  enrichment    = as.numeric(df$Enrichment),
  enrichment_se = as.numeric(df$Enrichment_SE),
  enrichment_p  = as.numeric(df$Enrichment_p_value)
)
writeLines(toJSON(out, auto_unbox = TRUE, digits = 15),
           file.path(outdir, "enrich_synth.json"))
cat("\nWrote fixtures/enrich_synth.json\n")
unlink(list.files(".", pattern = "\\.log$", full.names = TRUE))
