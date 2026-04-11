#!/usr/bin/env python3
"""Benchmark the Python `genomicsem` bindings on the same PGC data used by
benchmark_perf.R.

Writes two CSVs that benchmark_perf.R then merges into the master report:

  python_results.csv      func, impl, time_s, peak_mb, error
  python_artifacts.json   {func: <serialised output for equivalence checks>}

Run from the bench/ directory:
  ../.venv/bin/python run_python_bench.py

Environment variables honoured:
  BENCH_CORES    integer; used for userGWAS/parallel_analysis thread budget.
  BENCH_BIG      "1" to include the 25 000-SNP head-to-head size.
"""
from __future__ import annotations

import csv
import gc
import glob
import gzip
import json
import os
import resource
import sys
import time
import traceback
from pathlib import Path

import numpy as np

import genomicsem as gs

BENCH_DIR = Path(__file__).resolve().parent
os.chdir(BENCH_DIR)

# ---------------------------------------------------------------------------
# Configuration (mirrors benchmark_perf.R)
# ---------------------------------------------------------------------------
ptsd_glob = glob.glob("data/*PTSD*")
ptsd_file = ptsd_glob[0] if ptsd_glob else glob.glob("data/pts_*ea*")[0]
RAW_TRAITS = [
    "data/anxiety.meta.full.cc.tbl.gz",
    "data/ocd_aug2017.gz",
    ptsd_file,
]
LD = "data/eur_w_ld_chr/"
HM3 = "data/w_hm3.snplist"
REF_1000G = "data/reference.1000G.maf.0.005.txt.gz"
TRAIT_NAMES = ["ANX", "OCD", "PTSD"]
SAMPLE_PREV = [0.5, 0.5, 0.5]
POP_PREV = [0.16, 0.02, 0.07]
N_OVERRIDES = [None, 9725.0, 5831.0]

BENCH_CORES = int(os.environ.get("BENCH_CORES", "0") or 0) or None
BENCH_BIG = os.environ.get("BENCH_BIG", "0") in ("1", "true", "TRUE", "True")

# userGWAS head-to-head size grid, matched to benchmark_perf.R.
UG_SIZES = [1000, 5000, 10000, 25000] if BENCH_BIG else [1000, 5000, 10000]

MUNGED_TRAITS = [f"{t}.sumstats.gz" for t in TRAIT_NAMES]

if not Path(REF_1000G).exists():
    sys.exit(f"Missing reference file: {REF_1000G}")

# ---------------------------------------------------------------------------
# Bench infrastructure
# ---------------------------------------------------------------------------
results: list[dict] = []
artifacts: dict[str, object] = {}


def _peak_rss_mb() -> float:
    """Return peak RSS for this process in MB (cross-platform via resource)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss units: bytes on macOS, KB on Linux
    return rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024


def run_bench(name: str, impl: str, fn) -> object:
    """Time `fn()`, capture peak RSS, and append to `results`. Returns the
    function's return value (or None on failure)."""
    gc.collect()
    rss_before = _peak_rss_mb()
    err = ""
    val = None
    t0 = time.perf_counter()
    try:
        val = fn()
        elapsed = time.perf_counter() - t0
    except Exception as e:  # noqa: BLE001
        elapsed = time.perf_counter() - t0
        err = f"{type(e).__name__}: {e}"
        traceback.print_exc()
    rss_after = _peak_rss_mb()
    results.append({
        "func":    name,
        "impl":    impl,
        "time_s":  round(elapsed, 6),
        "peak_mb": round(max(rss_before, rss_after), 1),
        "error":   err,
    })
    return val


def store_artifact(key: str, value: object) -> None:
    """Stash a small JSON-serialisable artifact for cross-impl equivalence checks."""
    artifacts[key] = value


# ---------------------------------------------------------------------------
# 0. Munge prerequisite (skip if already done by R)
# ---------------------------------------------------------------------------
print(f"Python bench: starting (BENCH_CORES={BENCH_CORES}, BENCH_BIG={BENCH_BIG})", flush=True)
if not all(Path(p).exists() for p in MUNGED_TRAITS):
    print("[prereq] munging via Python (R hasn't done it yet)", flush=True)
    gs.munge(
        files=RAW_TRAITS, hm3=HM3, trait_names=TRAIT_NAMES,
        n=None, info_filter=0.9, maf_filter=0.01, out=".",
    )

# ---------------------------------------------------------------------------
# 1. munge (write to a separate name to avoid overwriting the shared output)
# ---------------------------------------------------------------------------
print("[1/12] munge", flush=True)
py_munge_names = [f"{t}_benchPy" for t in TRAIT_NAMES]
for name in py_munge_names:
    Path(f"{name}.sumstats.gz").unlink(missing_ok=True)

run_bench("munge", "Python", lambda: gs.munge(
    files=RAW_TRAITS, hm3=HM3, trait_names=py_munge_names,
    info_filter=0.9, maf_filter=0.01, out=".",
))
for name in py_munge_names:
    Path(f"{name}.sumstats.gz").unlink(missing_ok=True)

# ---------------------------------------------------------------------------
# 2. ldsc — keep result for downstream
# ---------------------------------------------------------------------------
print("[2/12] ldsc", flush=True)
ldsc_res = run_bench("ldsc", "Python", lambda: gs.ldsc(
    traits=MUNGED_TRAITS, sample_prev=SAMPLE_PREV, population_prev=POP_PREV,
    ld=LD, wld=LD, trait_names=TRAIT_NAMES, n_blocks=200,
))
if ldsc_res is None:
    ldsc_res = gs.ldsc(
        traits=MUNGED_TRAITS, sample_prev=SAMPLE_PREV, population_prev=POP_PREV,
        ld=LD, wld=LD, trait_names=TRAIT_NAMES, n_blocks=200,
    )

# ldsc_res is a PyLdscResult; commonfactor/usermodel/rgmodel/user_gwas
# accept it directly via pyany_to_ldsc_result, which reads .s/.v/.i_mat/.n/.m_total.
S_arr = np.asarray(ldsc_res.s, dtype=np.float64)
V_arr = np.asarray(ldsc_res.v, dtype=np.float64)
store_artifact("ldsc_S", S_arr.tolist())


def _loadings_from_result(res: dict | None, op: str = "=~") -> list[float]:
    """Extract factor loadings from a columnar parameters dict."""
    if res is None:
        return []
    params = res.get("parameters") if isinstance(res, dict) else None
    if not isinstance(params, dict):
        return []
    ops = params.get("op", [])
    ests = params.get("est", [])
    return [e for e, o in zip(ests, ops) if o == op]


# ---------------------------------------------------------------------------
# 3. commonfactor
# ---------------------------------------------------------------------------
print("[3/12] commonfactor", flush=True)
cf = run_bench("commonfactor", "Python",
               lambda: gs.commonfactor(ldsc_res, "DWLS"))
store_artifact("commonfactor_loadings", _loadings_from_result(cf))

# ---------------------------------------------------------------------------
# 4. usermodel
# ---------------------------------------------------------------------------
print("[4/12] usermodel", flush=True)
model_str = (
    f"F1 =~ NA*{TRAIT_NAMES[0]} + {TRAIT_NAMES[1]} + {TRAIT_NAMES[2]}\n"
    "F1 ~~ 1*F1\n"
    + "\n".join(f"{t} ~~ {t}" for t in TRAIT_NAMES)
)
um = run_bench("usermodel", "Python",
               lambda: gs.usermodel(ldsc_res, model_str, "DWLS"))
store_artifact("usermodel_loadings", _loadings_from_result(um))

# ---------------------------------------------------------------------------
# 5. rgmodel
# ---------------------------------------------------------------------------
print("[5/12] rgmodel", flush=True)
# Empty model → auto common-factor R matrix. Signature:
# rgmodel(covstruc, model, std_lv=True, estimation=True[DWLS], sub=None)
rg = run_bench("rgmodel", "Python", lambda: gs.rgmodel(ldsc_res, ""))
if isinstance(rg, dict):
    store_artifact("rgmodel_R", rg.get("R"))

# ---------------------------------------------------------------------------
# 6. sumstats
# ---------------------------------------------------------------------------
print("[6/12] sumstats", flush=True)
Path("out_bench").mkdir(exist_ok=True)
py_ss_path = "out_bench/py_merged.tsv"

# gs.sumstats fans reference + per-trait reads across a local rayon
# pool; pass `cores` explicitly so this matches the bench-wide budget.
run_bench("sumstats", "Python", lambda: gs.sumstats(
    files=RAW_TRAITS, ref_dir=REF_1000G, trait_names=TRAIT_NAMES,
    info_filter=0.6, maf_filter=0.01, out=py_ss_path,
    parallel=True, cores=BENCH_CORES,
))

# ---------------------------------------------------------------------------
# 7. write.model
# ---------------------------------------------------------------------------
print("[7/12] write.model", flush=True)
loadings = [[0.7], [0.6], [0.5]]
run_bench("write.model", "Python", lambda: gs.write_model(
    loadings=loadings, names=TRAIT_NAMES, cutoff=0.3,
))

# ---------------------------------------------------------------------------
# Prepare per-size GWAS subsets (1000 / 5000 / 10000 [/ 25000], plus full).
# Builds them up front from the Python-side merged sumstats so the
# userGWAS sweep below can reuse them.
# ---------------------------------------------------------------------------
if not Path(py_ss_path).exists():
    gs.sumstats(
        files=RAW_TRAITS, ref_dir=REF_1000G, trait_names=TRAIT_NAMES,
        info_filter=0.6, maf_filter=0.01, out=py_ss_path,
        parallel=True, cores=BENCH_CORES,
    )

with open(py_ss_path) as f:
    py_total_snps = sum(1 for _ in f) - 1  # drop header

py_ug_sizes = sorted({min(n, py_total_snps) for n in UG_SIZES if n > 0})


def _materialise_subset(dst: str, n_rows: int) -> None:
    with open(py_ss_path) as src, open(dst, "w") as out:
        header = next(src)
        out.write(header)
        for i, line in enumerate(src):
            if i >= n_rows:
                break
            out.write(line)


subset_paths: dict[int, str] = {}
for n_snp in py_ug_sizes:
    dst = f"out_bench/py_subset_{n_snp}.tsv"
    _materialise_subset(dst, n_snp)
    subset_paths[n_snp] = dst

# ---------------------------------------------------------------------------
# 8. commonfactorGWAS — DISABLED (matches benchmark_perf.R section 8).
# R upstream crashes on this PGC input via an unguarded solve() AND uses
# a different identification scheme, so a like-for-like comparison is
# apples-to-oranges. userGWAS below exercises the per-SNP GWAS path.
# ---------------------------------------------------------------------------
print("[8/12] commonfactorGWAS (SKIPPED — see comment above)", flush=True)

# ---------------------------------------------------------------------------
# 9. userGWAS — scaling sweep (parallel only, matches R library section).
# ---------------------------------------------------------------------------
print("[9/12] userGWAS (scaling sweep)", flush=True)
gwas_model = (
    f"F1 =~ NA*{TRAIT_NAMES[0]} + {TRAIT_NAMES[1]} + {TRAIT_NAMES[2]}\n"
    "F1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP"
)

print(f"  head-to-head sizes: {py_ug_sizes} (capped at {py_total_snps}, cores={BENCH_CORES})",
      flush=True)
ug_par_output: object | None = None
for n_snp in py_ug_sizes:
    print(f"  -> N={n_snp}", flush=True)
    out = run_bench("userGWAS", f"Python (N={n_snp})",
                    lambda p=subset_paths[n_snp]: gs.user_gwas(
                        ldsc_res, p, gwas_model,
                        parallel=True, cores=BENCH_CORES))
    if n_snp == py_ug_sizes[0]:
        ug_par_output = out

if ug_par_output is not None and isinstance(ug_par_output, dict):
    snps = ug_par_output.get("snps", {}) if isinstance(ug_par_output.get("snps"), dict) else {}
    conv = snps.get("converged", []) if isinstance(snps, dict) else []
    store_artifact("userGWAS_par_converged", sum(1 for x in conv if x))

# Full-scale "real-world deployment" point, matching the R library
# `Rust (N=..., full)` row. Runs on the entire merged sumstats file.
print(f"  -> Python full-scale: N={py_total_snps}", flush=True)
run_bench("userGWAS", f"Python (N={py_total_snps}, full)",
          lambda: gs.user_gwas(ldsc_res, py_ss_path, gwas_model,
                               parallel=True, cores=BENCH_CORES))

# ---------------------------------------------------------------------------
# 10. paLDSC
# ---------------------------------------------------------------------------
print("[10/12] paLDSC", flush=True)
pa = run_bench("paLDSC", "Python",
               lambda: gs.parallel_analysis(S_arr, V_arr, r=100, p=0.95,
                                            parallel=True, cores=BENCH_CORES))
if isinstance(pa, dict):
    store_artifact("paLDSC_observed", pa.get("observed"))

# ---------------------------------------------------------------------------
# 11. summaryGLS
# ---------------------------------------------------------------------------
print("[11/12] summaryGLS", flush=True)
k = S_arr.shape[0]
kstar = k * (k + 1) // 2
tril = np.tril_indices(k)
y_gls = S_arr[tril]
X_gls = np.arange(1, kstar + 1, dtype=np.float64).reshape(-1, 1)

gls = run_bench("summaryGLS", "Python",
                lambda: gs.summary_gls(x=X_gls, y=y_gls, v=V_arr, intercept=True))
if isinstance(gls, dict):
    store_artifact("summaryGLS_beta", list(gls.get("beta", [])))

# ---------------------------------------------------------------------------
# 12. simLDSC — match the R wrapper's work: load LD scores + total M from
# the LD directory, then hand pre-loaded arrays to `gs.sim_ldsc`. The
# wall-clock time therefore includes the same file I/O the R library
# and CLI paths do.
# ---------------------------------------------------------------------------
print("[12/12] simLDSC", flush=True)


def _load_ld_scores(ld_dir: str) -> tuple[np.ndarray, float]:
    """Read chr 1-22 LD score files and M_5_50 totals from `ld_dir`."""
    ld_parts: list[np.ndarray] = []
    m_total = 0.0
    for chr_i in range(1, 23):
        score_f = Path(ld_dir) / f"{chr_i}.l2.ldscore.gz"
        if score_f.exists():
            with gzip.open(score_f, "rt") as fh:
                header = fh.readline().rstrip("\n").split("\t")
                l2_col = header.index("L2") if "L2" in header else -1
                if l2_col < 0:
                    continue
                vals = [float(line.split("\t")[l2_col])
                        for line in fh if line.strip()]
                ld_parts.append(np.asarray(vals, dtype=np.float64))
        m_f = Path(ld_dir) / f"{chr_i}.l2.M_5_50"
        if m_f.exists():
            m_total += sum(float(x) for x in m_f.read_text().split())
    if not ld_parts:
        raise FileNotFoundError(f"No LD scores found in {ld_dir}")
    return np.concatenate(ld_parts), (m_total or float(sum(len(p) for p in ld_parts)))


def _run_sim_ldsc() -> np.ndarray:
    ld_scores, m_total = _load_ld_scores(LD)
    return gs.sim_ldsc(
        s_matrix=S_arr,
        n_per_trait=[5000.0] * k,
        ld_scores=ld_scores.tolist(),
        m=m_total,
    )


sim_out = run_bench("simLDSC", "Python", _run_sim_ldsc)
if sim_out is not None:
    arr = np.asarray(sim_out)
    store_artifact("simLDSC_shape", list(arr.shape))

# ---------------------------------------------------------------------------
# Persist
# ---------------------------------------------------------------------------
import shutil
shutil.rmtree("out_bench", ignore_errors=True)

with open("python_results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["func", "impl", "time_s", "peak_mb", "error"])
    w.writeheader()
    w.writerows(results)
print("python_results.csv written", flush=True)

with open("python_artifacts.json", "w") as f:
    json.dump(artifacts, f, default=str)
print("python_artifacts.json written", flush=True)
