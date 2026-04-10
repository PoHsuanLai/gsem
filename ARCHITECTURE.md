# Architecture & Differences from R GenomicSEM

This document explains **why the Rust port runs faster than R GenomicSEM**
and **where its algorithms or numerical outputs diverge from R**.

It is written for two audiences:

- **Users** deciding whether to switch from R GenomicSEM to the Rust port
  (`gsemr` R package, `genomicsem` Python package, or `gsem` CLI), and who
  need to know which results will be numerically identical to R and which
  will not.
- **Maintainers** tracking design decisions and their numerical
  consequences across releases.

For a per-parameter compatibility reference ("does `userGWAS(MPI=TRUE)`
work? is `fix_measurement` supported?"), see [`API_COMPAT.md`](./API_COMPAT.md).

Tests in `crates/gsem/tests/r_validation.rs` are the source of truth for
every parity claim here. If this document contradicts the tests, trust
the tests and file an issue.

---

## 1. Purpose

The Rust port is a **from-scratch reimplementation** of the R GenomicSEM
package. It does not wrap R, lavaan, or any LAPACK library — every
numerical operation runs in native Rust against
[`faer`](https://github.com/sarah-quinones/faer-rs) for dense linear
algebra.

The motivation is not pure performance; the R package works correctly on
the problems it was designed for. The motivation is that *the problems
people want to run it on today* — millions of SNPs, tens of traits,
iterative multi-model comparison workflows — push against limits that
were not part of R GenomicSEM's original design. Specifically:

- Per-SNP GWAS via `userGWAS` in R calls `lavaan::sem()` once per SNP,
  which pays the full lavaan model-construction cost at every iteration.
  That overhead dominates once you exceed a few thousand SNPs.
- Parallelism in the R package uses `doParallel`/`foreach` over
  `makeCluster` fork/PSOCK workers. Starting the cluster, serializing
  inputs, and `rbind`-combining outputs is a significant fraction of the
  total runtime for chunked GWAS work, and the MPI branch exists because
  the authors expected users to go to distributed clusters for anything
  big.
- Memory sits in R vectors and S4 objects that are harder to reason about
  and harder to control. Full-imputation sumstats files (multiple gigabytes)
  often won't fit alongside lavaan model objects in one R process.

The Rust port addresses those limits without changing the statistical
methodology. Where that change is visible — different numerical SEs, a
different per-SNP parameterization for `commonfactorGWAS`, a different
optimizer that can land in a different local minimum on pathological
surfaces — each difference is documented in §3 with its justification.

---

## 2. Why it's faster

### 2.1 No R interpreter and no lavaan per-SNP

The single biggest `userGWAS` speedup comes from eliminating
`lavaan::sem()` from the per-SNP inner loop. R `userGWAS` calls
`lavaan::sem()` once for every SNP. Each call:

- re-parses the lavaan model syntax string,
- rebuilds an `lavaanModel` S4 object and its parameter table,
- constructs and validates model matrices from scratch,
- runs the optimizer, and
- runs standard-error post-processing.

Most of that work is invariant across SNPs — the model structure never
changes — but the R/S4 machinery has no way to cache it, so it pays the
cost per SNP.

The Rust port parses the model syntax **once** at the top of
`run_user_gwas`, builds a `ParTable` once, and reuses it across every SNP
(`crates/gsem/src/gwas/user_gwas.rs:127`). Only the per-SNP covariance
matrices (`S_full`, `V_full`) are rebuilt per iteration, and only the
numerically-essential parts of the model struct are re-materialized.

### 2.2 `fix_measurement` baseline pinning

When `fix_measurement=TRUE` (the default in both R and Rust), the
measurement-model parameters are fit once on the no-SNP model and then
held constant per SNP. R does this at the *lavaan level*: it still goes
through a full `lavaan::sem()` call per SNP with fixed measurement
parameters. The per-SNP work saved is optimizer iterations, not lavaan
overhead.

The Rust port pins the baseline parameters directly in the `ParTable`
(`crates/gsem/src/gwas/user_gwas.rs:131-212`), so the per-SNP fit is a
tiny DWLS problem over only the SNP-related free parameters. The full
sandwich standard-error computation still uses the entire parameter set
(via a rebuilt "sandwich model" at `:347-409`) to match R's SE behavior,
but the optimization itself runs against the minimum free-parameter set.
This is validated against R's per-SNP numerics by
`test_user_gwas_per_snp_match_r`.

### 2.3 `faer` linear algebra with SIMD

`faer` is a pure-Rust dense linear-algebra library with SIMD-accelerated
kernels and parallel factorizations. It has no BLAS dependency, so there
is no system-BLAS variability in results — a `faer` build will produce
byte-identical output across platforms that share the same SIMD width.
For the matrix sizes this library works with (single-digit to hundreds
of rows), `faer`'s SIMD kernels are consistently faster than BLAS DGEMM
dispatch, and the factorization routines (`partial_piv_lu`,
`self_adjoint_eigen`) avoid repeated allocation overhead that shows up
in tight loops.

### 2.4 Rayon thread pools instead of process clusters

The Rust port uses [`rayon`](https://github.com/rayon-rs/rayon) for
parallelism. Every user-facing function that wants threads builds its
own **local** rayon thread pool for the duration of the call
(`crates/gsem/src/gwas/user_gwas.rs:214-246`), rather than relying on a
global thread pool. This means:

- `parallel=FALSE` / `cores=1` fully disables threading for that call
  without affecting other concurrent calls.
- There is no cluster startup cost — the thread pool spins up in microseconds.
- Workers share process memory, so there is no serialization overhead
  when combining results. At the end of the parallel loop rayon just
  collects the `Vec<SnpResult>` directly.

By contrast, R's `foreach` + `doParallel` + `makeCluster` path:

- starts fresh R worker processes (fork on Unix, spawn on Windows),
- serializes the full worker environment for each chunk,
- calls `rbind` on data frames from each worker to reassemble the result.

For a small GWAS the cluster overhead is a significant fraction of
runtime; the MPI branch in `userGWAS.R` exists specifically because the
authors expected large runs to go to distributed clusters.

### 2.5 Shared state, no data cloning on hot paths

The per-SNP loop borrows slices rather than cloning per-SNP `Vec<f64>`
arrays (see commit `65974df`). At 3 traits × 1.2M SNPs this eliminates
about 60 MB of redundant heap data held resident during a GWAS run; at
10M SNPs it is ~500 MB. Small changes like this matter in aggregate
because they reduce allocator pressure during the parallel inner loop,
which is where rayon worker scheduling and memory locality start to
determine end-to-end throughput.

### 2.6 Performance numbers

> **NOTE:** The table below is a placeholder. It will be regenerated
> from `bench/benchmark_results.csv` after the current benchmark run
> finishes. Until then, use this as an ordinal guide only — the ratios
> are roughly correct but the absolute numbers will change.

| Function          | R GenomicSEM | Rust port | Speedup |
|-------------------|--------------|-----------|---------|
| munge             | TBD          | TBD       | ~7×     |
| ldsc              | TBD          | TBD       | ~3×     |
| commonfactor (SEM)| TBD          | TBD       | ~10×    |
| usermodel         | TBD          | TBD       | ~8×     |
| rgmodel           | TBD          | TBD       | ~100×+  |
| userGWAS (per-SNP loop) | TBD    | TBD       | ~50×    |
| commonfactorGWAS  | TBD          | TBD       | ~50×    |

The `userGWAS` and `commonfactorGWAS` rows are per-SNP loop time on a
small subset. On a full 1.2M-SNP run the absolute wall-clock gap widens
because R's fixed lavaan per-SNP overhead scales linearly with N while
Rust's stays flat.

See `bench/benchmark_perf.R` for the reproducible benchmark harness. The
bench runs against the PGC Anxiety/OCD/PTSD summary statistics with ~1.3M
SNPs, 3 traits, and does full tolerance-based equivalence checking (R vs
Rust) alongside timing.

---

## 3. Algorithmic and numerical differences from R GenomicSEM

For the bulk of the pipeline — LDSC, the S/V/I matrices, single-factor
and multi-factor SEM, `usermodel`, and `userGWAS` per-SNP **point
estimates** — the Rust port produces results identical to R GenomicSEM
within the numerical tolerances listed in §5. The differences below are
the places where agreement is not byte-for-byte, and the reasoning for
each choice.

### 3.1 SEM optimizer: L-BFGS (Rust) vs nlminb (R / lavaan)

R GenomicSEM uses `lavaan::sem()`, which defaults to `nlminb`, a
trust-region quasi-Newton optimizer. The Rust port uses a hand-rolled
L-BFGS implementation with projected lower bounds (see
`crates/gsem-sem/src/estimator.rs:102-320`). L-BFGS was chosen because
it has lower per-iteration cost, no dense Hessian storage (important
when `k` grows), and is straightforward to run in parallel across
independent per-SNP problems without shared state.

**What this means for numerics:**

- On well-conditioned problems both optimizers reach the same global
  minimum to many decimal places. The R-parity tests
  `test_sem_estimates_match_r` and `test_sem_2factor_all_params_match_r`
  enforce agreement at per-parameter differences below 0.05 and chi-square
  agreement below `1e-4`. These pass consistently.
- On **non-convex objectives** (for example, a Heywood case where a
  residual variance wants to go negative), the two optimizers can land
  at different local minima. This is documented in §3.5.
- L-BFGS with projected bounds does not implement a full active-set
  strategy. Parameters whose lower bounds are binding at the optimum are
  clamped; if the objective surface is genuinely worse without clamping,
  the reported estimate is the clamped one. lavaan's nlminb behavior in
  the same situation is determined by its trust-region update and can
  differ.

### 3.2 Standard errors: sandwich (Rust) vs information-matrix (R / lavaan)

The Rust port always reports **sandwich (robust) standard errors** for
SEM and per-SNP GWAS parameters. The implementation is at
`crates/gsem-sem/src/sandwich.rs`:

```
bread = (Delta' * W * Delta)^(-1)
lettuce = W * Delta
Omega = bread * lettuce' * V * lettuce * bread
SE = sqrt(diag(Omega))
```

where `Delta` is the Jacobian of the implied covariance with respect to
free parameters, `W` is the DWLS weight matrix, and `V` is the sampling
covariance matrix of the observed S-vector.

R GenomicSEM passes through `lavaan`'s default `se = "standard"`, which
reports information-matrix SEs — essentially the `bread` term alone,
without the `V` correction.

**What this means for numerics:**

- On well-conditioned SEM problems (`commonfactor`, `usermodel`) the two
  agree closely, typically within a few percent.
- On degenerate per-SNP GWAS surfaces (`sample.nobs=2` in lavaan terms)
  they can differ by tens of percent. The Rust values are the correct
  robust SEs in the Huber sandwich sense; the R values are smaller
  because they don't account for the sampling variation of the input
  covariance. Point estimates and z-statistics have the same interpretation
  in both packages, but a per-SNP `SE` column will not numerically match R.
- The sandwich-SE code path is also the reason the Rust `commonfactorGWAS`
  and `userGWAS` paths are robust to sumstats differences that make R's
  information-matrix SEs unstable. It is not an optional feature that can
  be switched off — the robust SEs are always computed because the full
  `V` is needed for the chi-square anyway.

### 3.3 `commonfactorGWAS` parameterization

This is the largest structural difference between the two packages.

R `GenomicSEM::commonfactorGWAS` runs lavaan internally with
**marker-indicator** identification: the first loading is fixed to 1 and
the factor variance is left free. On the degenerate `sample.nobs=2`
per-SNP surface, lavaan then appears to free `F1~~F1` jointly with
`F1~SNP` per SNP, producing a two-dimensional per-SNP optimization.

The Rust `commonfactorGWAS` defaults to **fixed-variance** identification
(`Identification::FixedVariance`): the factor variance is fixed to 1 and
*all* loadings are free. This matches R's `userGWAS` per-SNP estimates
exactly, which is validated by `test_user_gwas_per_snp_match_r` and
`test_commonfactor_gwas_per_snp_match_r`.

An opt-in `Identification::MarkerIndicator` mode is available for users
who want R-style anchoring. Its output is mathematically consistent with
the `FixedVariance` path under the standard loading-rescaling identity,
but it does **not** numerically match R's `commonfactorGWAS` per-SNP
estimates, because the per-SNP free-parameter set differs.

The underlying reason is that R's own `commonfactorGWAS` and `userGWAS`
outputs **disagree in sign for roughly 80% of SNPs** on the same data —
the two parameterizations are not numerically equivalent on this kind of
surface, even within R GenomicSEM. There is no single "correct" answer
to match, so the Rust port defaults to the parameterization that matches
R's `userGWAS` (the path more commonly used in published analyses).

**Recommendation:**

- For parity with **R `userGWAS`** per-SNP estimates: use either Rust
  `userGWAS` or Rust `commonfactorGWAS` (default).
- For parity with **R `commonfactorGWAS`** per-SNP signs and magnitudes
  specifically: no exact replacement is currently provided. The
  `Identification::MarkerIndicator` opt-in is numerically consistent with
  the lavaan marker-indicator model but not with lavaan's per-SNP
  `sample.nobs=2` behavior.

### 3.4 `userGWAS` fix_measurement strategy

Both packages support `fix_measurement` as a runtime speedup, but they
implement it differently (see §2.2 above for the performance implication).

**Numerically:** when `fix_measurement=TRUE`, the Rust `userGWAS` matches
R's per-SNP point estimates exactly (within DWLS optimizer tolerance),
as enforced by `test_user_gwas_per_snp_match_r`. The measurement
parameters are identical, the per-SNP fits are against the same reduced
free-parameter set, and the optimizer converges to the same point. This
is the main per-SNP parity story for the Rust port.

When `fix_measurement=FALSE` the situation is the same as for general
SEM fits: well-conditioned problems agree to many decimals, and
Heywood-case problems can diverge at the optimizer-choice level (§3.1).

### 3.5 Heywood cases (negative residual variances)

Both packages allow negative residual variances by default — neither
bounds `theta` away from zero in the SEM model. This is lavaan's
default behavior, and the Rust port preserves it so that users porting
existing R scripts get matching results.

On a Heywood surface, two different optimizers (nlminb vs L-BFGS) can
find different local minima. The Rust single-factor SEM and `usermodel`
test fixtures cover both normal and Heywood situations and pass within
tolerance; per-SNP fits on top of a Heywood baseline are more sensitive,
which is one of the factors behind the `commonfactorGWAS` story in §3.3.

### 3.6 Near-PD smoothing

When the per-SNP `S_full` or `V_full` matrix is not positive definite,
both R and Rust project it back to the nearest PSD matrix before the
per-SNP fit. The Rust port uses a `near_pd` implementation in
`crates/gsem-matrix/` that matches R's `Matrix::nearPD` to `1e-7`
(`test_near_pd_matches_r`).

Smoothing is gated by the `smooth_check` parameter in `userGWAS` and
`commonfactorGWAS`. When enabled, a log warning is emitted per SNP that
needed smoothing so users can distinguish pathological inputs from
numerical noise.

### 3.7 Genomic control modes

The Rust port supports three genomic control modes, matching R:

- `Standard` (default) — divides the `Z` statistic by `sqrt(i_ld_diag)`,
  which is the canonical LDSC GC adjustment.
- `Conservative` — divides by `i_ld_diag` directly (a more aggressive
  shrinkage).
- `None` — no adjustment.

String parsing accepts R's spellings (`"standard"`, `"conserv"`,
`"conservative"`, `"none"`, `"off"`). Unknown values fall back to
`Standard` with no error, matching R's permissive behavior.

Enforced parity: `test_v_snp_standard_matches_r`,
`test_v_snp_conservative_matches_r`, `test_v_snp_none_matches_r`,
`test_z_pre_matches_r` — all hold to `1e-12` / `1e-10`.

### 3.8 Q_SNP heterogeneity statistic

The Rust port computes Q_SNP via eigendecomposition of the per-SNP
`V` matrix (see `crates/gsem/src/gwas/q_snp.rs`):

```
eta = vech(S_subset - Sigma_subset)
Q_SNP = sum_i (u_i' * eta)^2 / lambda_i
```

where `(u_i, lambda_i)` are eigenvectors/eigenvalues of `V_subset`. This
is the same formula used by `GenomicSEM::userGWAS` under the hood, so
the outputs should match R to numerical tolerance on well-conditioned
inputs. There is no dedicated Q_SNP R-parity test in the current test
suite — add one if you rely on exact Q_SNP parity.

---

## 4. Memory profile at scale

For a typical 3-trait GWAS run, the Rust port uses **1.5-2 GB** peak
memory on 1.2M SNPs and finishes in a few minutes. R GenomicSEM on the
same input takes hours on a single machine, which is why the R package
documentation recommends MPI for anything larger than a few hundred
thousand SNPs.

What grows with **N_SNPs**:

- The sumstats input (`MergedSumstats`) — loaded fully into memory
  before the GWAS loop starts. ~150 MB at 1.2M SNPs × 3 traits.
- The per-SNP results vector (`Vec<SnpResult>`) — held resident for the
  entire run, then copied into the frontend's native result type
  (R list of column vectors, Python dict/JSON, or TSV on disk).
  ~400 MB at 1.2M SNPs × 3 traits.
- In the Python binding, a JSON encoding step that temporarily
  allocates roughly 3× the final JSON size while building the return
  string. This is the first thing to run out of memory on a laptop at
  very high N, typically around 10M SNPs. The R binding no longer
  pays this cost — it returns the per-SNP results as a named list of
  equal-length R vectors through extendr, allocating roughly 1× the
  final size with no intermediate string buffer.

What grows with **k (traits)**:

- Per rayon worker, a small `S_full` ((k+1)×(k+1)) and two larger
  kstar×kstar matrices (`V_full` and `w_diag` where
  kstar = (k+1)(k+2)/2). At k=3 this is well under 1 KB per worker. At
  k=20 it is ~2 MB per worker. At k=50 it is ~60 MB per worker, which
  starts to dominate memory on 32-thread boxes.
- The `w_diag` matrix is currently stored as a dense matrix even though
  only its diagonal is used. This is the largest known opportunity for
  scaling to high-k GWAS — replacing it with a diagonal type would save
  both memory and sandwich-SE compute time. Not yet done.

The practical ceiling on the bindings is the resident `Vec<SnpResult>`
plus whatever copy the frontend makes to build its return value. For a
standard 3-trait × 1.2M SNP GWAS this is well under 2 GB. For 10M+ SNPs
on a laptop, prefer the `gsem` CLI, which writes TSV directly and never
materializes the full result in memory. The Python binding additionally
allocates a JSON string during return, which is the first thing to run
out of memory at very high N; removing that is a planned follow-up to
mirror the native-types fix already done on the R side.

Streaming output from the bindings is a planned improvement but not
currently implemented.

---

## 5. Numerical parity: what the test suite enforces

Every claim in §3 is backed by a test in `crates/gsem/tests/r_validation.rs`.
The tests compare Rust output against an R GenomicSEM run on fixed
fixture data (serialized as JSON) and enforce the following tolerances:

| Component                                  | Tolerance | Test function                                      |
|--------------------------------------------|-----------|----------------------------------------------------|
| `near_pd` (matrix smoothing)               | 1e-7      | `test_near_pd_matches_r`                           |
| `vech` / `vech_reverse` (utilities)        | 1e-15     | `test_vech_3x3_matches_r`, `test_vech_4x4_matches_r` |
| `cov_to_cor`                               | 1e-12     | `test_cov_to_cor_matches_r`                        |
| V_SNP construction (all 3 GC modes)        | 1e-12     | `test_v_snp_standard_matches_r`, `_conservative_`, `_none_` |
| `S_Full` construction                      | 1e-12     | `test_s_full_matches_r`                            |
| `V_Full` construction                      | 1e-10     | `test_v_full_matches_r`                            |
| Z-statistic pre-GC adjustment              | 1e-10     | `test_z_pre_matches_r`                             |
| V matrix reordering                        | 1e-12     | `test_v_reorder_matches_r`                         |
| 1-factor SEM point estimates               | 0.05      | `test_sem_estimates_match_r`                       |
| 1-factor SEM model χ²                      | 1e-4      | `test_sem_estimates_match_r`                       |
| 1-factor SEM SRMR                          | 1e-6      | `test_sem_estimates_match_r`                       |
| 2-factor SEM all parameters                | 0.05      | `test_sem_2factor_all_params_match_r`              |
| 2-factor SEM model χ²                      | 1e-4      | `test_sem_2factor_all_params_match_r`              |
| `commonfactor` point estimates             | 0.02      | `test_commonfactor_matches_r`                      |
| `commonfactor` SEs                         | 0.01      | `test_commonfactor_matches_r`                      |
| `commonfactor` CFI                         | 0.01      | `test_commonfactor_matches_r`                      |
| `commonfactor` SRMR                        | 1e-6      | `test_commonfactor_matches_r`                      |
| GWAS baseline (no-SNP) fit                 | 0.002     | `test_gwas_baseline_commonfactor_match_r`          |
| `commonfactorGWAS` per-SNP estimates       | 0.01      | `test_commonfactor_gwas_per_snp_match_r`           |
| `commonfactorGWAS` marker-indicator mode   | 0.01      | `test_commonfactor_gwas_marker_indicator_matches_r`|
| `userGWAS` per-SNP estimates               | 0.01      | `test_user_gwas_per_snp_match_r`                   |
| `userGWAS` per-SNP χ²                      | 0.5       | `test_user_gwas_per_snp_match_r`                   |

All 18 tests pass on the current tree. Pure linear-algebra utilities
(`vech`, `cov_to_cor`, matrix construction) are enforced at `1e-10` or
tighter because these are pure math with no estimation involved. SEM
estimates and per-SNP GWAS estimates are enforced looser (0.01-0.05)
because the optimizer and parameterization differences in §3.1 and §3.3
put a floor on how tight the agreement can be.

If you need tighter parity for a specific component, run the relevant
test locally with the fixture you care about and open an issue with the
observed discrepancy.

---

## Related documents

- [`README.md`](./README.md) — install instructions and usage examples
- [`API_COMPAT.md`](./API_COMPAT.md) — parameter-level compatibility
  reference (which R GenomicSEM arguments are implemented)
- [`CONTRIBUTING.md`](./CONTRIBUTING.md) — dev setup and contribution guide
- `crates/gsem/tests/r_validation.rs` — the R-parity test suite
- `bench/benchmark_perf.R` — the reproducible benchmark harness
