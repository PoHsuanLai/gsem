"""Parity tests for the genomicsem Python binding.

Assert the Rust-backed Python wrappers reproduce the SAME values that R
GenomicSEM produced for the committed synthetic fixtures (the Python counterpart
of crates/gsem/tests/r_validation_*.rs). These catch binding-layer bugs —
numpy<->faer marshalling, dict-covstruc handling, default drift — that the
Rust-core tests can't see.

Fixtures live at the repo-root tests/fixtures/. If that tree isn't reachable
(e.g. wheel installed standalone), the tests skip rather than fail.

Cross-binding triangulation: both this suite and the R testthat suite
(bindings/r/tests/testthat/) assert against the SAME committed R-reference
fixtures. Python == fixture and R == fixture together imply Python == R, so the
two bindings are pinned to each other through the shared reference without a
fragile cross-process harness.
"""

import json
import os
from pathlib import Path

import numpy as np
import pytest

import genomicsem as g


def _fixtures_dir():
    env = os.environ.get("GSEM_FIXTURES_DIR")
    if env and Path(env, "ldsc_synth.json").exists():
        return Path(env)
    here = Path(__file__).resolve()
    for parent in here.parents:
        cand = parent / "tests" / "fixtures"
        if (cand / "ldsc_synth.json").exists():
            return cand
    return None


FIX = _fixtures_dir()
needs_fixtures = pytest.mark.skipif(
    FIX is None, reason="repo-root tests/fixtures/ not found (set GSEM_FIXTURES_DIR)"
)


def load(name):
    with open(FIX / f"{name}.json") as fh:
        return json.load(fh)


@needs_fixtures
def test_ldsc_matches_r_reference():
    fx = load("ldsc_synth")
    munged = [str(FIX / f) for f in fx["munged_files"]]
    ld = str(FIX / fx["ld_dir"])

    res = g.ldsc(
        munged,
        [float("nan")] * len(munged),
        [float("nan")] * len(munged),
        ld,
        ld,
        trait_names=fx["trait_names"],
        chr=fx["chr"],
        n_blocks=fx["n_blocks"],
    )

    # S/V entries are O(1e-3); the binding reproduces R's values to ~1e-7
    # (the same jackknife-level agreement the Rust-core test accepts), so a
    # combined rtol/atol is the honest parity check, not a loose fudge.
    np.testing.assert_allclose(np.asarray(res.s), np.asarray(fx["s"]), rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(np.asarray(res.i_mat), np.asarray(fx["i"]), rtol=1e-3, atol=1e-5)
    np.testing.assert_allclose(np.asarray(res.v), np.asarray(fx["v"]), rtol=1e-4, atol=1e-6)


@needs_fixtures
def test_summary_gls_matches_r_reference():
    fx = load("summary_gls")
    x = np.asarray(fx["x"], dtype=float)
    y = np.asarray(fx["y"], dtype=float)
    v = np.asarray(fx["v"], dtype=float)

    # x already carries the intercept column -> intercept=False (matches fixture).
    # summary_gls returns a dict with keys beta / se / z / p.
    res = g.summary_gls(x, y, v, intercept=False)

    np.testing.assert_allclose(np.asarray(res["beta"]), np.asarray(fx["betas"]), rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(np.asarray(res["se"]), np.asarray(fx["se"]), rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(np.asarray(res["z"]), np.asarray(fx["z"]), rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(np.asarray(res["p"]), np.asarray(fx["pvals"]), rtol=1e-7, atol=1e-10)


@needs_fixtures
def test_summary_gls_intercept_prepend():
    # The fixture's x already carries an all-ones intercept in column 0. Passing
    # x WITHOUT that column and intercept=True must reproduce the same fit as
    # passing the full x with intercept=False. Exercises the intercept-prepend
    # branch of the binding (summary_gls_rust), untested elsewhere.
    fx = load("summary_gls")
    x_full = np.asarray(fx["x"], dtype=float)
    assert np.allclose(x_full[:, 0], 1.0), "fixture col 0 should be the intercept"
    x_no_int = x_full[:, 1:]
    y = np.asarray(fx["y"], dtype=float)
    v = np.asarray(fx["v"], dtype=float)

    res = g.summary_gls(x_no_int, y, v, intercept=True)
    np.testing.assert_allclose(np.asarray(res["beta"]), np.asarray(fx["betas"]), rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(np.asarray(res["se"]), np.asarray(fx["se"]), rtol=1e-8, atol=1e-10)


@needs_fixtures
def test_commonfactor_matches_r_reference():
    fx = load("commonfactor")
    k = np.asarray(fx["s"]).shape[0]
    covstruc = {
        "s": np.asarray(fx["s"], dtype=float),
        "v": np.asarray(fx["v"], dtype=float),
        "i_mat": np.eye(k),  # unused by the fit; binding requires the slot
    }

    # commonfactor returns a dict; res["parameters"] is columnar
    # (parallel lhs/op/rhs/est lists).
    res = g.commonfactor(covstruc, "DWLS")
    params = res["parameters"]
    got = {
        (params["lhs"][j], params["op"][j], params["rhs"][j]): params["est"][j]
        for j in range(len(params["lhs"]))
    }

    # The fixture stores parameters as a JSON array of {lhs,op,rhs,est} rows.
    for row in fx["parameters"]:
        key = (row["lhs"], row["op"], row["rhs"])
        assert key in got, f"missing parameter {key}"
        # Common factor orientation unidentified -> compare |est|.
        assert abs(abs(got[key]) - abs(row["est"])) < 1e-4, f"{key}: {got[key]} vs {row['est']}"
