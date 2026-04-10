#!/usr/bin/env python3
"""Benchmark the Python `genomicsem` bindings on the same PGC data used by
benchmark_perf.R.

Writes two CSVs that benchmark_perf.R then merges into the master report:

  python_results.csv      func, impl, time_s, peak_mb, error
  python_artifacts.json   {func: <serialised output for equivalence checks>}

Run from the bench/ directory:
  ../.venv/bin/python run_python_bench.py
"""
from __future__ import annotations

import csv
import gc
import glob
import json
import os
import resource
import sys
import tempfile
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
GWAS_SUBSET_N = 1000

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
print("Python bench: starting", flush=True)
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
for n in py_munge_names:
    Path(f"{n}.sumstats.gz").unlink(missing_ok=True)

run_bench("munge", "Python", lambda: gs.munge(
    files=RAW_TRAITS, hm3=HM3, trait_names=py_munge_names,
    info_filter=0.9, maf_filter=0.01, out=".",
))
for n in py_munge_names:
    Path(f"{n}.sumstats.gz").unlink(missing_ok=True)

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
covstruc_json = ldsc_res.to_json()
store_artifact("ldsc_S", np.asarray(ldsc_res.s).tolist())

# ---------------------------------------------------------------------------
# 3. commonfactor
# ---------------------------------------------------------------------------
print("[3/12] commonfactor", flush=True)
cf_json = run_bench("commonfactor", "Python",
                    lambda: gs.commonfactor(covstruc_json, "DWLS"))
if cf_json is not None:
    cf = json.loads(cf_json)
    loadings = [p["est"] for p in cf["parameters"] if p["op"] == "=~"]
    store_artifact("commonfactor_loadings", loadings)

# ---------------------------------------------------------------------------
# 4. usermodel
# ---------------------------------------------------------------------------
print("[4/12] usermodel", flush=True)
model_str = (
    f"F1 =~ NA*{TRAIT_NAMES[0]} + {TRAIT_NAMES[1]} + {TRAIT_NAMES[2]}\n"
    "F1 ~~ 1*F1\n"
    + "\n".join(f"{t} ~~ {t}" for t in TRAIT_NAMES)
)
um_json = run_bench("usermodel", "Python",
                    lambda: gs.usermodel(covstruc_json, model_str, "DWLS"))
if um_json is not None:
    um = json.loads(um_json)
    loadings = [p["est"] for p in um["parameters"] if p["op"] == "=~"]
    store_artifact("usermodel_loadings", loadings)

# ---------------------------------------------------------------------------
# 5. rgmodel
# ---------------------------------------------------------------------------
print("[5/12] rgmodel", flush=True)
# rgmodel signature: (covstruc_json, model, std_lv=True, estimation=True[DWLS], sub=None)
# Pass empty model to get the auto common-factor R matrix.
rg_json = run_bench("rgmodel", "Python",
                    lambda: gs.rgmodel(covstruc_json, ""))
if rg_json is not None:
    try:
        rg = json.loads(rg_json)
        store_artifact("rgmodel_R", rg.get("R"))
    except Exception:
        pass

# ---------------------------------------------------------------------------
# 6. sumstats
# ---------------------------------------------------------------------------
print("[6/12] sumstats", flush=True)
Path("out_bench").mkdir(exist_ok=True)
py_ss_path = "out_bench/py_merged.tsv"

run_bench("sumstats", "Python", lambda: gs.sumstats(
    files=RAW_TRAITS, ref_dir=REF_1000G, trait_names=TRAIT_NAMES,
    info_filter=0.6, maf_filter=0.01, out=py_ss_path,
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
# Prepare GWAS subset
# ---------------------------------------------------------------------------
if not Path(py_ss_path).exists():
    gs.sumstats(
        files=RAW_TRAITS, ref_dir=REF_1000G, trait_names=TRAIT_NAMES,
        info_filter=0.6, maf_filter=0.01, out=py_ss_path,
    )

py_subset_path = "out_bench/py_subset.tsv"
with open(py_ss_path) as f_in, open(py_subset_path, "w") as f_out:
    header = next(f_in)
    f_out.write(header)
    for i, line in enumerate(f_in):
        if i >= GWAS_SUBSET_N:
            break
        f_out.write(line)

# ---------------------------------------------------------------------------
# 8. commonfactorGWAS — parallel + serial
# ---------------------------------------------------------------------------
print("[8/12] commonfactorGWAS", flush=True)
cfg_outputs: dict[str, object] = {}
for parallel in (True, False):
    tag = "par" if parallel else "seq"
    j = run_bench(f"commonfactorGWAS", f"Python ({tag})",
                  lambda p=parallel: gs.commonfactor_gwas(
                      covstruc_json, py_subset_path, parallel=p))
    if j is not None:
        try:
            cfg_outputs[tag] = json.loads(j)
        except Exception:
            pass

if "par" in cfg_outputs:
    snp_est = []
    for entry in cfg_outputs["par"]:
        snp_row = next((p for p in entry.get("params", [])
                        if p.get("op") == "~" and p.get("rhs") == "SNP"), None)
        snp_est.append({
            "SNP": entry.get("SNP"),
            "est": snp_row["est"] if snp_row else None,
            "se":  snp_row["se"]  if snp_row else None,
        })
    store_artifact("commonfactorGWAS_par", snp_est)

# ---------------------------------------------------------------------------
# 9. userGWAS — parallel + serial
# ---------------------------------------------------------------------------
print("[9/12] userGWAS", flush=True)
gwas_model = (
    f"F1 =~ NA*{TRAIT_NAMES[0]} + {TRAIT_NAMES[1]} + {TRAIT_NAMES[2]}\n"
    "F1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP"
)
ug_outputs: dict[str, object] = {}
for parallel in (True, False):
    tag = "par" if parallel else "seq"
    j = run_bench(f"userGWAS", f"Python ({tag})",
                  lambda p=parallel: gs.user_gwas(
                      covstruc_json, py_subset_path, gwas_model, parallel=p))
    if j is not None:
        try:
            ug_outputs[tag] = json.loads(j)
        except Exception:
            pass

if "par" in ug_outputs:
    n_conv = sum(1 for e in ug_outputs["par"] if e.get("converged"))
    store_artifact("userGWAS_par_converged", n_conv)

# ---------------------------------------------------------------------------
# 10. paLDSC
# ---------------------------------------------------------------------------
print("[10/12] paLDSC", flush=True)
S_arr = np.asarray(ldsc_res.s)
V_arr = np.asarray(ldsc_res.v)
S_json = json.dumps(S_arr.tolist())
V_json = json.dumps(V_arr.tolist())

pa_json = run_bench("paLDSC", "Python",
                    lambda: gs.parallel_analysis(S_json, V_json, r=100))
if pa_json is not None:
    try:
        pa = json.loads(pa_json)
        store_artifact("paLDSC_observed", pa.get("observed"))
    except Exception:
        pass

# ---------------------------------------------------------------------------
# 11. summaryGLS
# ---------------------------------------------------------------------------
print("[11/12] summaryGLS", flush=True)
k = S_arr.shape[0]
kstar = k * (k + 1) // 2
tril = np.tril_indices(k)
y = S_arr[tril].tolist()
X = [[float(i + 1)] for i in range(kstar)]
V_list = V_arr.tolist()

gls_json = run_bench("summaryGLS", "Python",
                     lambda: gs.summary_gls(x=X, y=y, v=V_list, intercept=True))
if gls_json is not None:
    try:
        gls = json.loads(gls_json)
        store_artifact("summaryGLS_beta",
                       [r.get("beta") for r in gls] if isinstance(gls, list) else gls)
    except Exception:
        pass

# ---------------------------------------------------------------------------
# 12. simLDSC — skipped: the Python `sim_ldsc` binding takes raw LD scores +
# total M, while the R/CLI version reads the LD directory itself. The two
# APIs aren't directly comparable, so we omit Python from this row rather
# than report a misleading number.
# ---------------------------------------------------------------------------
print("[12/12] simLDSC (skipped for Python — different API)", flush=True)

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
