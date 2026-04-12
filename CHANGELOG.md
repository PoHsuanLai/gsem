# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Changes on `master` past the most recent release tag live here until the
next version is cut.

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

[Unreleased]: https://github.com/PoHsuanLai/gsem/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/PoHsuanLai/gsem/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/PoHsuanLai/gsem/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/PoHsuanLai/gsem/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/PoHsuanLai/gsem/releases/tag/v0.1.0
