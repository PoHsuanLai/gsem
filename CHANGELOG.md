# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Changes on `master` past the most recent release tag live here until the
next version is cut.

## [0.2.0] — 2026-06-13

### Added

- **Ported `read_fusion`** (`gsem::io::fusion_reader` + CLI `gsem
  read-fusion`). Reads raw FUSION TWAS `.dat` association files and converts
  each trait's `TWAS.Z` into a standardized gene-expression effect/SE (binary
  traits via the liability conversion `effect/√(effect²·HSQ + π²/3)`,
  continuous via `Z/√(N·HSQ)`), supports the permutation-Z path, then
  inner-joins traits on Gene+Panel into the merged TWAS format that
  `twas_reader` / `multiGene` / `userGWAS(TWAS=TRUE)` consume. This is the
  on-ramp from raw FUSION output into the Rust TWAS path (previously only
  pre-merged files could be read). Validated against live R `read_fusion`.

- **Ported `subSV`** (`gsem_matrix::vech::subset_sv`). Subsets `vech(S)` and
  the `V` sampling-covariance block by a set of 1-based vech positions, for
  both the with-diagonal (`TYPE="S"`/`"S_Stand"`) and off-diagonal
  (`TYPE="R"`) numbering. Validated against R to 1e-12. (R's matrix-input
  path has an undefined-`RMATRIX` bug; the port matches the bug-free
  `LDSC_OBJECT` path.)

- **Ported `summaryGLSbands`** numeric core
  (`gsem::stats::gls::summary_gls_bands`). GLS fit plus the confidence-band
  data — predictor grid, fitted line, and ±`BAND_SIZE`·SE envelope at
  `INTERVALS` points, with `INTERCEPT`/`QUAD`/`CONTROLVARS` support.
  Validated against R to 1e-9. The ggplot rendering is not ported.

- **Per-option R-equivalence coverage for the drop-in surface.** New
  real-package fixtures and tests for previously-untested option branches:
  `userGWAS` `estimation="ML"` / `GC={conserv,none}` / `Q_SNP` / `std.lv`;
  `ldsc` `stand=TRUE` (S_Stand/V_Stand) / `select="ODD"` / `chisq.max` /
  liability-scale; `sumstats` `ambig=TRUE` (full numeric parity); a
  non-degenerate misspecified-model CFI/chisq fixture; the `paLDSC` `diag`
  branch; and bindings parity suites (R testthat + Python pytest).

- **Coverage instrumentation in CI.** A `cargo-llvm-cov` coverage job with a
  ratcheting floor, Codecov upload + README badge, and an advisory R-side
  covr job.

- **R documentation site (pkgdown).** A published docs site for the `gsemr`
  binding (reference index + Compatibility/Architecture articles generated
  from the repo-root docs), deployed to GitHub Pages, with a navbar link to
  the upstream GenomicSEM project.

### Fixed

- **`std.lv` left the factor scale unidentified.** `parse_model(std_lv=true)`
  freed the auto-added latent variance instead of fixing it to 1 (lavaan's
  `std.lv=TRUE` convention: `auto.fix.first=FALSE` + latent variances fixed
  to 1). Any std.lv model — `userGWAS` / `usermodel` / `rgmodel(std_lv=)` —
  estimated an unidentified factor scale and was ~12% off R. Now the
  auto-add fixes the latent variance under `std_lv`. (`crates/gsem-sem/src/syntax.rs`)

- **Q_SNP heterogeneity statistic + LDSC intercept floor.** `compute_q_snp`
  computed the wrong quadratic form, and the LDSC intercept diagonal was not
  floored at 1.0 before building the per-SNP `V` (R's `userGWAS.R:156`).
  Both fixed; Q_SNP now matches R to ~1e-9.
  (`crates/gsem/src/gwas/{q_snp,gc_correction,user_gwas}.rs`)

---

The following entries accumulated on `master` after 0.1.3 and ship in 0.2.0:

### Added

- **R-equivalence coverage for `multiSNP`, `multiGene`, `simLDSC`, and
  `hdl`.** New live-R fixtures and tests validate the four previously
  untested functions, closing the last function-coverage gaps. `multiSNP`
  and `multiGene` reproduce R's augmented `S_Full`/`V_Full` matrices to
  1e-12/1e-14; `simLDSC`'s deterministic per-SNP Z covariance matches R's
  exact `varZ`/`covZ` algebra to 1e-9; `hdl` reproduces R's genetic-
  covariance matrix `S` to optimiser-level precision on a synthetic LD
  panel built in R's exact `.rda`/`.bim` format.

### Changed

- **`hdl` evaluated in the eigenspace (match R).** gsem's HDL likelihood
  was fed raw LD scores as `lam` and the raw per-SNP `bhat`, whereas R
  GenomicSEM/HDL evaluates the likelihood in the LD-block eigenspace:
  `lam` = eigenvalues of the block correlation matrix and
  `bstar = Vᵀ·bhat` (projection onto eigenvectors). `LdPiece` now carries
  the per-piece eigenvalues/eigenvectors (read from the `.rda` panel, or
  the new `chr*.eigen.tsv` text file emitted by `convert_hdl_panels`), the
  per-trait reference N uses R's `median(N)`, the likelihood floor matches
  R's `exp(-18)`, and the per-piece optimiser is a projected-gradient
  method mirroring R's `optim(L-BFGS-B)`. (`crates/gsem-ldsc/src/hdl.rs`)

### Fixed

- **`simLDSC` phenotypic-overlap term.** The off-diagonal environmental
  contribution to the per-SNP Z covariance was
  `rPheno·n_overlap·√(NᵢNⱼ)/n_snps`, carrying a spurious `√(NᵢNⱼ)/n_snps`
  factor. R GenomicSEM uses `rPheno_ij · N_ij/√(NᵢNⱼ)` (= `rPheno·n_overlap`
  for scalar sample overlap). Fixed and exposed as the deterministic
  `per_snp_z_cov`, validated against R's construction.
  (`crates/gsem/src/stats/simulation.rs`)

- **`multiSNP`/`multiGene` `V_Full` construction.** `run_multi_snp` built
  `V_Full` as a diagonal approximation that ignored the LDSC intercept
  matrix. The new `build_multi_snp_sv` reproduces R's full sampling
  covariance — SNP-trait variances `(SE·I_diag·varSNP)²`, cross-trait
  within-SNP, cross-SNP within-trait, and cross-SNP cross-trait blocks,
  plus the trait-trait `V_LD` block — matching R bit-for-bit.
  (`crates/gsem/src/gwas/multi_snp.rs`)

#### Documented R bugs (gsem implements the correct behavior)

- **`multiSNP`/`multiGene` constant cross-SNP cross-trait LD.** R weights
  every cross-SNP cross-trait sampling-covariance cell by a single
  constant LD value (`LD2[(f²−f)/2]`, the last lower-triangle entry)
  rather than the actual SNP-pair LD. gsem uses the correct per-pair
  `LD[a,b]`; the two coincide when the off-diagonal LD is constant (as in
  the equivalence fixtures). A unit test locks gsem's corrected path.

- **`multiGene` aborts for k ≥ 2 traits.** R's `multiGene` assigns into
  `V_SNP[y,x]` — a variable that does not exist in `multiGene` (a typo for
  `V_Gene`) — so the function errors with `object 'V_SNP' not found`
  whenever there is more than one trait. gsem implements the corrected
  algorithm (identical to `multiSNP` with gene heritabilities as the
  "variances"); the reference is generated from a minimally-patched
  `multiGene`.

- **`hdl` intercept weak identification.** The per-piece HDL intercept
  enters the likelihood only as `int·lam/N` and is weakly identified on
  small panels; R's `optim` drifts on some LD blocks while gsem's gradient
  method recovers the construction-true intercept. The genetic-covariance
  matrix `S` (HDL's primary output) is robustly identified and matches R.

- **`sumstats` beta standardization.** gsem's `sumstats` wrote the raw
  input betas/SEs instead of standardizing them like R GenomicSEM. Since
  `userGWAS`/`commonfactorGWAS` build `cov(SNP, trait) = varSNP · beta`,
  this propagated a `sqrt(varSNP)` scale error into per-SNP GWAS output.
  Now implements all of R's modes — `OLS` (`beta = Z/√(N·varSNP)`),
  `linprob`, `se.logit`, the default logistic "none" transform, and
  `betas` pass-through — reading the P and N columns and deriving Z from
  P. Validated formula-for-formula against R for every mode.
  (`crates/gsem/src/sumstats.rs`)

- **`sumstats` ambiguous-SNP default.** R GenomicSEM removes
  strand-ambiguous SNPs (A/T, C/G) only when `ambig=TRUE`; its default
  keeps them. The R and Python bindings (and the CLI) mapped
  `keep_ambig = ambig`, so their default dropped ~⅓ of SNPs where R keeps
  them. Fixed to `keep_ambig = !ambig` across the R binding, Python
  binding, and CLI (which gains an `--ambig` flag), matching R's default.

- **`rgmodel` V_R dimension.** `rgmodel` returned `V_R` as the full
  `kstar × kstar` sampling covariance of `vech(R)` (including the fixed
  diagonal-1 correlations, which have ~zero variance). R GenomicSEM returns
  `V_R` for the off-diagonal correlations only (`k(k−1)/2` square). Fixed to
  extract the off-diagonal vech block, matching R's shape and ordering (so
  the bindings now return R's shape too). (`crates/gsem-sem/src/rgmodel.rs`)

- **`enrich` model-based functional enrichment.** gsem's `enrich`
  (`enrichment_test`) computed the classic LDSC heritability-enrichment ratio
  (`prop_h²/prop_SNPs`, the original Finucane statistic), not R GenomicSEM's
  model-based `enrich`. R fits a SEM to a baseline annotation, fixes the
  regressions/loadings (per `fix`) to their baseline estimates, re-fits each
  annotation freeing the remaining parameters, and reports per-parameter
  `enrichment = (est_annot/est_baseline)/Prop` with SE and a 1-sided p-value.
  Added `gsem_sem::enrich_model::model_enrichment` (+ a shared `gsem_sem::fit`
  helper) implementing this with the regressions/covariances/variances fix
  modes, wired through the R binding (`enrich(s_covstruc, model, params, fix)`)
  and the Python binding (`model_enrichment`). Validated against live R
  GenomicSEM on a factor-variance covstruc. gsem's classic ratio is retained
  as `enrichment_test`. (`crates/gsem-sem/src/enrich_model.rs`)

- **`s_ldsc` overlap-weighted partitioned heritability.** Stratified LDSC
  returned per-annotation `S` as the raw coefficient contribution `tau · M`
  (R's `S_Tau`), missing R's overlap weighting. R computes each category's
  partitioned (co)heritability as `S = overlap · (tau · M)`, where
  `overlap[f,a] = M(f,a)/M_a` and `M(f,a)` is the number of SNPs in both
  annotations `f` and `a` (from the `.annot.gz` membership + `.frq` MAF
  filter). These agree only for disjoint annotations and diverge for the
  realistic overlapping baselineLD model. Added a `.annot.gz`/`.frq`
  cross-product reader and applied the overlap transform to both `S` and the
  jackknife-based `V`; the result now exposes `s_annot`/`v_annot` (R's
  `S`/`V`) and `s_tau`/`v_tau` (R's `S_Tau`/`V_Tau`). Threaded through the R
  binding, Python binding, and CLI (new `--frq` flag). Validated against R
  on synthetic overlapping annotations.
  (`crates/gsem-ldsc/src/stratified.rs`, `crates/gsem-ldsc/src/annot_reader.rs`)

### Added

- **In-CI R-equivalence coverage for the LDSC half** (`ldsc`, `munge`,
  `sumstats`) via small synthetic GWAS + LD-score fixtures generated by
  running R GenomicSEM (`tests/generate_synthetic_reference.R`).
  Previously only the data-dependent bench script covered these.

- **R-equivalence coverage for `summaryGLS`, `paLDSC` (observed eigenvalue
  spectrum), `rgmodel`, and `write.model`** (factor→indicator structure)
  via `tests/generate_covstruc_reference.R`, plus **ML-estimator and
  3-factor SEM** cases — broadening the function, estimator, and model-size
  combinations checked against R.

- **Tightened R-validation tolerances** to observed precision (SEM /
  commonfactor estimates now checked at 1e-6 vs the prior 5e-2; fit
  indices at 1e-10), with a `MEASURE=1` env-gated max-diff tracker.

## [0.1.3] — 2026-04-13

### Fixed

- **LDSC cross-trait allele alignment.** Z-scores were not flipped when
  the effect allele (A1) differed between traits in a cross-trait pair.
  This corrupted the Z₁·Z₂ product, producing wrong-sign genetic
  covariance and inflated cross-trait intercepts (e.g. 0.64 instead of
  0.05). Diagonal (heritability) estimates were unaffected.
  (`crates/gsem-ldsc/src/lib.rs`)

- **LDSC per-trait chi² filtering.** R GenomicSEM filters each trait's
  SNPs independently (threshold = `max(0.001·max(N), 80)`) before
  merging cross-trait pairs. gsemr was computing the threshold per-pair
  using `max(N_j, N_k)`, which kept SNPs that R would have excluded.
  Now matches R's pre-filter-then-merge behavior; cross-trait SNP
  counts agree exactly.

## [0.1.2] — 2026-04-13

### Changed

- **Publish CI R version pinned to `oldrel`** (currently R 4.4.x) instead
  of `release` (R 4.5.0). Pre-built binaries are now compiled against the
  previous R release so users who haven't upgraded to the latest R can
  install them directly.

### Fixed

- **Equality constraints in SEM models** (`a*V1 + a*V2` syntax) now
  correctly share a single free parameter across aliased cells. Previously
  each labelled term got its own independent parameter, silently ignoring
  the equality constraint. The Jacobian (analytical delta matrix) is
  updated to accumulate gradient contributions from aliased cells.

- **Labelled loadings treated as free parameters.** A label on the first
  indicator (e.g. `F1 =~ a*V1 + a*V2`) was incorrectly caught by the
  "fix first loading to 1" default, making it fixed instead of free.
  Labelled terms now bypass the first-loading rule, matching lavaan
  semantics.

- **Observed variable detection for covariance-only models.** Models like
  `RL ~~ MDD` (cross-covariance, no `=~` terms) failed with "model has
  0 observed variables". The fallback heuristic now collects all unique
  variable names from covariance and regression terms, not just
  self-covariances.

- **L-BFGS convergence flag on line-search exit.** When the line search
  exhausts all step sizes without improvement, the optimizer has reached
  a stationary point (or constrained boundary). This exit path now sets
  `converged = true`, matching the practical semantics of R GenomicSEM's
  `optim.force.converged = TRUE`. Estimates were already correct; only
  the flag was wrong.

### Added

- **Heywood case warnings.** After every SEM fit, negative diagonal
  entries in Theta (residual variances) or Psi (factor variances) are
  detected and logged at WARN level via the `log` crate. The warning
  names the offending variable(s) and suggests adding lower bounds or
  re-specifying the model. Surfaces in R, Python, and CLI.

- New tests for equality constraints: shared-parameter count reduction,
  value propagation through aliases, implied covariance correctness,
  and analytical-vs-numerical Jacobian agreement with labels.

- Negative variance detection test (`test_negative_variances_detected`).

## [0.1.1] — 2026-04-11

### Added

- **Pre-built R binary packages for Linux, macOS, and Windows** attached
  to every GitHub release. Users no longer need a Rust toolchain to
  install `gsemr` — pick the platform-native file from the release page
  and run `install.packages(<url>, repos = NULL)`. Source tarball +
  `remotes::install_github` remain as the "from source" fallback for
  unmatched R versions and dev builds.
- R package now ships with full Windows build support: `configure.win`,
  `src/Makevars.win`, and `-Wl,--export-all-symbols` in the MinGW link
  line so `R .Call("wrap__*_rust", ...)` resolves at runtime on
  Windows. `R CMD check` is green on ubuntu / macos / windows.
- `workflow_dispatch` trigger on `publish.yml` with a `release_tag`
  input, so the R-only job path can be manually re-run against an
  existing release (to ship a late R binary or fix without burning a
  crates.io / PyPI version slot).
- `readme = "README.md"` and `[project.urls]` in
  `bindings/python/pyproject.toml` so the PyPI project page renders a
  description with links to the repo, issue tracker, and architecture
  docs. (The 0.1.0 wheels shipped before this landed — the PyPI 0.1.0
  page will continue to show "no project description" until 0.1.1 is
  uploaded.)

### Changed

- **`bindings/python/README.md`**: absolute GitHub URLs instead of
  relative `../../API_COMPAT.md` / `../../ARCHITECTURE.md` links
  (which don't render on PyPI); function list expanded from 10 to all
  17 exports and grouped as *Core pipeline* vs *Advanced*.
- **Root `README.md`**: R install section restructured — binary
  packages as the primary path, source install demoted to a "From
  source" subsection. New Rust MSRV (1.88+) badge.
- **R package internal cleanup**:
  - Explicit `importFrom(stats, pnorm, setNames)` /
    `importFrom(utils, read.table, write.table)` in `NAMESPACE`,
    clearing the "no visible global function definition" notes.
  - `@param out` added to `sumstats.R` + regenerated `sumstats.Rd`
    (fixes an "Undocumented arguments in Rd file" WARNING).
  - `{chr}` in `ldsc.R`'s `@param ld` description escaped to
    `\code{<chr>.l2.ldscore.gz}` (fixes "Lost braces" NOTE).
  - `configure` / `configure.win` print `rustc --version` / `cargo
    --version` before the build (fixes "No rustc version reported"
    WARNING).
  - `configure` / `configure.win` now `rm -rf src/rust/target` after
    copying out `libgsemr.a`, so R CMD check doesn't scan a duplicate
    `rust/target/release/libgsemr.a` and double-report the Rust
    stdlib's `_exit` / `abort` / `exit` symbols.

### Removed

- **`PyLdscResult.to_json()` and `gsem.LdscResult.from_json(s)`** — the
  `json: String` field and the `conversions::json_to_ldsc` helper that
  backed them. These were dead compat from the pre-NumPy era; every hot
  path has been reading `s` / `v` / `i_mat` through NumPy getters for a
  while. **This is the only user-visible break in 0.1.1.** Callers who
  were serializing LDSC results to disk should use `pickle` or pass the
  result object directly into downstream functions (all 17 exported
  functions accept it).

### Fixed

- **`read_sumstats` aborted on NA tokens**, so `.sumstats.gz` files that
  worked in `GenomicSEM::ldsc` failed in `gsemr::ldsc` with a cryptic
  `gsemr::ldsc error: invalid N` long before `sample.prev` /
  `population.prev` were even evaluated. GenomicSEM tolerates these
  files because its reader calls `na.omit(read_delim(...))`; gsemr now
  matches that behaviour — rows whose `N` or `Z` field is an NA token
  (empty, `.`, `NA`, `NaN`, `N/A`, `NULL`; case-insensitive) are
  silently dropped, with a one-line INFO log naming the file and the
  dropped-row count. Genuine parse failures now report file path, line
  number, and the offending value, and `load_trait_data` wraps the
  error with the failing file path so multi-trait runs tell you which
  input went bad. (`crates/gsem/src/io/gwas_reader.rs`)
- **`publish.yml` skip cascade**: `check-r-package` and
  `upload-r-to-release` no longer inherit the implicit `if: success()`
  gate, which was evaluating to false on `workflow_dispatch` runs
  (because `publish-crates` is skipped there). They now use
  `if: always() && needs.<upstream>.result == 'success'` and proceed
  as long as their direct needs are green.
- **Versioned workspace deps**: root `Cargo.toml` internal deps now
  carry `version = "0.1.0"` alongside `path`, required by
  `cargo publish` for inter-crate references.
- **R binding Cargo patch**: `bindings/r/src/rust/Cargo.toml` declares
  internal crates at crates.io versions and uses `[patch.crates-io]`
  to redirect to the in-repo sources when available. `configure`
  strips the patch block when the local crates aren't reachable
  (tarball install), so cargo falls back to the published crates
  cleanly.

## [0.1.0] — 2026-04-11

First tagged release of GenomicSEM-rs. A Rust rewrite of R GenomicSEM,
shipped across four distribution channels:

### Added

- **Rust crates** on crates.io:
  - [`gsem-matrix`](https://crates.io/crates/gsem-matrix) — matrix
    utilities (nearest PD, half-vec, PSD smoothing).
  - [`gsem-ldsc`](https://crates.io/crates/gsem-ldsc) — LD Score
    Regression with block jackknife.
  - [`gsem-sem`](https://crates.io/crates/gsem-sem) — SEM engine
    (DWLS / ML, lavaan syntax, sandwich SEs).
  - [`gsem`](https://crates.io/crates/gsem) — pipeline + CLI binary.
- **R package `gsemr`** — drop-in compatible with R GenomicSEM on all
  18 user-facing functions (`ldsc`, `usermodel`, `commonfactor`,
  `commonfactorGWAS`, `userGWAS`, `munge`, `sumstats`, `paLDSC`,
  `write.model`, `rgmodel`, `hdl`, `s_ldsc`, `convert_hdl_panels`,
  `enrich`, `simLDSC`, `multiSNP`, `multiGene`, `summaryGLS`). Every
  exported function ships with a worked `@examples` block (inline
  where a synthetic covstruc suffices, `\dontrun{}` for file-dependent
  examples).
- **Python package `genomicsem`** on PyPI — PyO3 bindings exposing 17
  functions mirroring the R API. Linux / macOS / Windows wheels. Every
  `#[pyfunction]` carries a NumPy-style `Examples` section in its
  `__doc__` accessible via `help(gsem.ldsc)` etc.
- **`gsem` CLI** — 16 subcommands mirroring the R function set, each
  with an `EXAMPLE:` block in its `--help` long output.
- Significant wall-clock speedups vs R GenomicSEM on the PGC benchmark
  (Anxiety / OCD / PTSD). See the README comparison table and
  `bench/benchmark_plots.pdf` for the detailed breakdown.

### Known issues in 0.1.0

- **PyPI 0.1.0 project page is blank.** The 0.1.0 wheels were published
  before `readme = "README.md"` was added to `pyproject.toml`, so
  <https://pypi.org/project/genomicsem/> shows "The author of this
  package has not provided a project description". Fixed on `master`
  for 0.1.1 onwards; the 0.1.0 page can't be updated because PyPI
  doesn't allow re-uploading the same version.
- **R package `0.1.0` tarball originally shipped without Windows
  support** (no `configure.win`, no `Makevars.win`). Fixed post-tag:
  `gsemr_0.1.0.zip` and the rebuilt `gsemr_0.1.0.tar.gz` with Windows
  support are now attached to the v0.1.0 release via a manual publish
  re-trigger.

[Unreleased]: https://github.com/PoHsuanLai/gsem/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/PoHsuanLai/gsem/compare/v0.1.3...v0.2.0
[0.1.3]: https://github.com/PoHsuanLai/gsem/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/PoHsuanLai/gsem/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/PoHsuanLai/gsem/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/PoHsuanLai/gsem/releases/tag/v0.1.0
