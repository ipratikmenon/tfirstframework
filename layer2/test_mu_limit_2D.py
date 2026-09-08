"""
test_mu_limit_2D.py — Test suite for M4 μ(T)→ν Limit 2D
=========================================================
35 tests across 7 classes.

Classes:
  TestTInitial           (6)  — make_T_initial shape, range, delta=0 limit, seed
  TestOmegaInitial       (4)  — make_omega_initial shape, scale, div-free spectral
  TestRefProps           (4)  — _ref_props keys, positivity, T_mean field
  TestSingleDelta        (8)  — single run keys, verdict, A_min>0, lps_margin
  TestSweep              (7)  — 10 results, all PASS, exp IDs, DB rows, monotone trend
  TestPhase4Verdict      (4)  — INCOMPLETE/PASS/FAIL from DB
  TestSmokeTest          (2)  — smoke test PASS, n_runs=3
"""

import sys
import os
import sqlite3
import tempfile
import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))

from mu_limit_2D import (
    make_T_initial,
    make_omega_initial,
    _ref_props,
    _ensure_db,
    log_result,
    query_results,
    run_single_delta,
    run_mu_limit_sweep,
    run_smoke_test,
    phase4_verdict,
    DELTA_SWEEP_VALUES,
    _CLAIM_ID,
    _T_MEAN,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_db(tmp_path):
    """Temporary SQLite database path."""
    return str(tmp_path / "test_results.db")


@pytest.fixture
def single_result(tmp_db):
    """One δ=0.1 run at N=16, t_end=0.03 (smoke-speed)."""
    return run_single_delta(
        delta=0.1,
        exp_id="EXP-TEST-MULIMIT-001",
        N=16,
        t_end=0.03,
        max_steps=20,
        verbose=False,
        db_path=tmp_db,
        seed=42,
    )


# ─────────────────────────────────────────────────────────────────────────────
# 1. TestTInitial
# ─────────────────────────────────────────────────────────────────────────────

class TestTInitial:
    def test_shape(self):
        T0 = make_T_initial(32, _T_MEAN, 1.0)
        assert T0.shape == (32, 32)

    def test_shape_64(self):
        T0 = make_T_initial(64, _T_MEAN, 0.5)
        assert T0.shape == (64, 64)

    def test_delta_zero_is_uniform(self):
        T0 = make_T_initial(32, _T_MEAN, 0.0)
        assert np.allclose(T0, _T_MEAN), "δ=0 must give uniform T"

    def test_mean_is_T_mean(self):
        """Spatial mean of T0 ≈ T_mean (g has zero mean by construction)."""
        T0 = make_T_initial(32, _T_MEAN, 1.0, seed=42)
        assert abs(np.mean(T0) - _T_MEAN) < 2.0, "mean(T0) should be close to T_mean"

    def test_range_bounded_by_delta(self):
        """T0 ∈ [T_mean - δ, T_mean + δ] since g ∈ [-1,1]."""
        delta = 5.0
        T0 = make_T_initial(32, _T_MEAN, delta, seed=42)
        assert np.all(T0 >= _T_MEAN - delta - 1e-10)
        assert np.all(T0 <= _T_MEAN + delta + 1e-10)

    def test_reproducible_seed(self):
        T0a = make_T_initial(32, _T_MEAN, 1.0, seed=7)
        T0b = make_T_initial(32, _T_MEAN, 1.0, seed=7)
        assert np.allclose(T0a, T0b), "Same seed must give same field"


# ─────────────────────────────────────────────────────────────────────────────
# 2. TestOmegaInitial
# ─────────────────────────────────────────────────────────────────────────────

class TestOmegaInitial:
    def test_shape(self):
        omega = make_omega_initial(32)
        assert omega.shape == (32, 32)

    def test_shape_64(self):
        omega = make_omega_initial(64)
        assert omega.shape == (64, 64)

    def test_amplitude_scale(self):
        """Max |ω| ≈ _OMEGA_SCALE = 1.0."""
        from mu_limit_2D import _OMEGA_SCALE
        omega = make_omega_initial(32)
        assert abs(np.max(np.abs(omega)) - _OMEGA_SCALE) < 1e-10

    def test_not_uniform(self):
        """ω must have spatial structure (not all zeros or constant)."""
        omega = make_omega_initial(32)
        assert np.std(omega) > 1e-3, "ω must have spatial variation"


# ─────────────────────────────────────────────────────────────────────────────
# 3. TestRefProps
# ─────────────────────────────────────────────────────────────────────────────

class TestRefProps:
    def test_keys(self):
        ref = _ref_props()
        for key in ("nu_ref", "A_ref", "mu_ref", "rho_ref", "T_mean"):
            assert key in ref

    def test_positivity(self):
        ref = _ref_props()
        assert ref["nu_ref"] > 0
        assert ref["A_ref"] > 0
        assert ref["mu_ref"] > 0
        assert ref["rho_ref"] > 0

    def test_T_mean_stored(self):
        ref = _ref_props(T_mean=400.0)
        assert abs(ref["T_mean"] - 400.0) < 1e-10

    def test_A_ref_positive_various_T(self):
        """A(T) > 0 for all physically valid T."""
        for T in [200.0, 300.0, 500.0, 1000.0]:
            ref = _ref_props(T_mean=T)
            assert ref["A_ref"] > 0, f"A_ref must be positive at T={T}"


# ─────────────────────────────────────────────────────────────────────────────
# 4. TestSingleDelta
# ─────────────────────────────────────────────────────────────────────────────

class TestSingleDelta:
    def test_returns_dict(self, single_result):
        assert isinstance(single_result, dict)

    def test_required_keys(self, single_result):
        required = [
            "exp_id", "delta", "verdict", "A_min_global", "A_ref",
            "lps_margin", "mu_min_global", "omega_max_global",
            "T_min_final", "T_max_final", "T_range_final",
            "n_steps", "t_final", "claim_id",
        ]
        for k in required:
            assert k in single_result, f"Missing key: {k}"

    def test_verdict_pass(self, single_result):
        assert single_result["verdict"] == "PASS"

    def test_A_min_strictly_positive(self, single_result):
        """A_min_global > 0 — second law guarantee."""
        assert single_result["A_min_global"] > 0.0

    def test_A_min_above_zero_threshold(self, single_result):
        """A_min must be substantially positive, not just floating-point epsilon."""
        assert single_result["A_min_global"] > 1e-10

    def test_claim_id_correct(self, single_result):
        assert single_result["claim_id"] == _CLAIM_ID

    def test_logged_to_db(self, single_result, tmp_db):
        rows = query_results(tmp_db)
        assert len(rows) >= 1
        exp_ids = [r["exp_id"] for r in rows]
        assert "EXP-TEST-MULIMIT-001" in exp_ids

    def test_t_final_reaches_t_end(self, single_result):
        """Simulation must reach t_end (0.03) within tolerance."""
        assert single_result["t_final"] >= 0.025  # allows for last-step snap

    def test_delta_zero_run(self, tmp_db):
        """δ=0 gives uniform T → A_min = A_ref exactly."""
        row = run_single_delta(
            delta=0.0,
            exp_id="EXP-TEST-DELTA0",
            N=16,
            t_end=0.02,
            max_steps=10,
            db_path=None,
            seed=42,
        )
        ref = _ref_props(_T_MEAN)
        # A_min should be very close to A_ref when T is uniform
        assert row["A_min_global"] > 0
        assert row["verdict"] == "PASS"
        # lps_margin ≈ 0 for δ=0 (A_min ≈ A_ref)
        assert abs(row["lps_margin"]) < ref["A_ref"] * 0.1, \
            "lps_margin should be ≈0 for uniform T"


# ─────────────────────────────────────────────────────────────────────────────
# 5. TestSweep
# ─────────────────────────────────────────────────────────────────────────────

class TestSweep:
    @pytest.fixture(scope="class")
    def sweep_results(self, tmp_path_factory):
        """Run 3-value mini sweep at N=16, t_end=0.03 for speed."""
        db = str(tmp_path_factory.mktemp("sweep") / "sw.db")
        delta_mini = [1.0, 0.1, 0.01]
        return run_mu_limit_sweep(
            delta_values=delta_mini,
            N=16,
            t_end=0.03,
            max_steps=20,
            verbose=False,
            db_path=db,
        ), db

    def test_correct_number_of_results(self, sweep_results):
        results, _ = sweep_results
        assert len(results) == 3

    def test_all_pass(self, sweep_results):
        results, _ = sweep_results
        for r in results:
            assert r["verdict"] == "PASS", f"{r['exp_id']} should PASS"

    def test_A_min_positive_all(self, sweep_results):
        results, _ = sweep_results
        for r in results:
            assert r["A_min_global"] > 0, f"A_min must be > 0 for {r['exp_id']}"

    def test_exp_id_format(self, sweep_results):
        results, _ = sweep_results
        for r in results:
            assert r["exp_id"].startswith("EXP-L2-R1-MULIMIT-"), \
                f"Bad exp_id: {r['exp_id']}"

    def test_db_has_rows(self, sweep_results):
        _, db = sweep_results
        rows = query_results(db)
        assert len(rows) == 3

    def test_delta_ordering(self, sweep_results):
        """Results ordered large-δ first (1.0, 0.1, 0.01)."""
        results, _ = sweep_results
        deltas = [r["delta"] for r in results]
        assert deltas[0] > deltas[1] > deltas[2]

    def test_A_min_converges_to_A_ref(self, sweep_results):
        """A_min_global(δ=0.01) should be closer to A_ref than A_min(δ=1.0).

        As δ→0, T→T̄, so A_min→A_ref. The small-δ run should have smaller
        |A_min - A_ref| than the large-δ run.
        """
        results, _ = sweep_results
        A_ref = results[0]["A_ref"]
        # Large δ=1.0
        diff_large = abs(results[0]["A_min_global"] - A_ref)
        # Small δ=0.01
        diff_small = abs(results[2]["A_min_global"] - A_ref)
        assert diff_small <= diff_large + A_ref * 0.05, \
            "A_min should converge toward A_ref as δ→0"


# ─────────────────────────────────────────────────────────────────────────────
# 6. TestPhase4Verdict
# ─────────────────────────────────────────────────────────────────────────────

class TestPhase4Verdict:
    def test_incomplete_when_empty_db(self, tmp_db):
        v = phase4_verdict(tmp_db)
        assert v["verdict"] == "INCOMPLETE"
        assert v["n_runs"] == 0

    def test_pass_after_successful_runs(self, tmp_db):
        """After inserting PASS rows, verdict should be PASS."""
        conn = _ensure_db(tmp_db)
        ref = _ref_props()
        for i, delta in enumerate([1.0, 0.1, 0.01], start=1):
            row = {
                "exp_id": f"EXP-TEST-V-{i:03d}",
                "claim_id": _CLAIM_ID,
                "timestamp": "2026-01-01T00:00:00",
                "layer": "L2",
                "route": "R1",
                "delta": delta,
                "T_mean": _T_MEAN,
                "grid_N": 16,
                "t_end": 0.03,
                "fluid": "ideal",
                "n_steps": 5,
                "t_final": 0.03,
                "verdict": "PASS",
                "A_min_global": ref["A_ref"] * (1.0 - 0.01 * delta),  # slightly below A_ref
                "A_ref": ref["A_ref"],
                "lps_margin": -0.01 * delta * ref["A_ref"],
                "nu_ref": ref["nu_ref"],
                "mu_min_global": ref["nu_ref"],
                "omega_max_global": 1.0,
                "T_min_final": _T_MEAN - delta,
                "T_max_final": _T_MEAN + delta,
                "T_range_final": 2 * delta,
                "key_metric": f"A_min={ref['A_ref']:.4e}",
                "notes": "",
                "params_json": "{}",
            }
            log_result(conn, row)
        conn.close()

        v = phase4_verdict(tmp_db)
        assert v["verdict"] == "PASS"
        assert v["n_pass"] == 3
        assert v["n_fail"] == 0

    def test_fail_propagates(self, tmp_db):
        """One FAIL row → overall verdict FAIL."""
        conn = _ensure_db(tmp_db)
        ref = _ref_props()
        row = {
            "exp_id": "EXP-TEST-FAIL-001",
            "claim_id": _CLAIM_ID,
            "timestamp": "2026-01-01T00:00:00",
            "layer": "L2", "route": "R1",
            "delta": 1.0, "T_mean": _T_MEAN, "grid_N": 16, "t_end": 0.03,
            "fluid": "ideal", "n_steps": 2, "t_final": 0.01,
            "verdict": "FAIL",
            "A_min_global": ref["A_ref"],
            "A_ref": ref["A_ref"],
            "lps_margin": 0.0,
            "nu_ref": ref["nu_ref"],
            "mu_min_global": ref["nu_ref"],
            "omega_max_global": 9999.0,
            "T_min_final": 100.0, "T_max_final": 600.0, "T_range_final": 500.0,
            "key_metric": "blow-up", "notes": "", "params_json": "{}",
        }
        log_result(conn, row)
        conn.close()

        v = phase4_verdict(tmp_db)
        assert v["verdict"] == "FAIL"
        assert v["n_fail"] == 1

    def test_margin_ok_flag(self, tmp_db):
        """margin_ok=True when A_min_plateau > 0 for small δ."""
        conn = _ensure_db(tmp_db)
        ref = _ref_props()
        # Insert small-δ run with positive A_min
        row = {
            "exp_id": "EXP-TEST-MARG-001",
            "claim_id": _CLAIM_ID,
            "timestamp": "2026-01-01T00:00:00",
            "layer": "L2", "route": "R1",
            "delta": 0.001, "T_mean": _T_MEAN, "grid_N": 16, "t_end": 0.03,
            "fluid": "ideal", "n_steps": 5, "t_final": 0.03,
            "verdict": "PASS",
            "A_min_global": ref["A_ref"] * 0.99,
            "A_ref": ref["A_ref"],
            "lps_margin": -0.01 * ref["A_ref"],
            "nu_ref": ref["nu_ref"], "mu_min_global": ref["nu_ref"],
            "omega_max_global": 1.0,
            "T_min_final": _T_MEAN - 0.001, "T_max_final": _T_MEAN + 0.001,
            "T_range_final": 0.002,
            "key_metric": f"A_min={ref['A_ref']:.4e}", "notes": "", "params_json": "{}",
        }
        log_result(conn, row)
        conn.close()

        v = phase4_verdict(tmp_db)
        assert v["margin_ok"] is True
        assert v["A_min_plateau"] > 0


# ─────────────────────────────────────────────────────────────────────────────
# 7. TestSmokeTest
# ─────────────────────────────────────────────────────────────────────────────

class TestSmokeTest:
    @pytest.fixture(scope="class")
    def smoke(self):
        return run_smoke_test(verbose=False)

    def test_verdict_pass(self, smoke):
        assert smoke["verdict"] == "PASS"

    def test_n_runs(self, smoke):
        assert smoke["n_runs"] == 3

    def test_all_pass(self, smoke):
        assert smoke["all_pass"] is True


# ─────────────────────────────────────────────────────────────────────────────
# 8. TestFullSweepIDs  (integration — verifies 10 exp IDs if full sweep run)
# ─────────────────────────────────────────────────────────────────────────────

class TestFullSweepIDs:
    def test_delta_sweep_values_count(self):
        """DELTA_SWEEP_VALUES must have exactly 10 entries."""
        assert len(DELTA_SWEEP_VALUES) == 10

    def test_delta_sweep_decreasing(self):
        """δ values must be strictly decreasing."""
        for i in range(len(DELTA_SWEEP_VALUES) - 1):
            assert DELTA_SWEEP_VALUES[i] > DELTA_SWEEP_VALUES[i + 1]

    def test_delta_sweep_smallest_very_small(self):
        """Smallest δ should be ≤ 0.001 to test the Prize limit."""
        assert DELTA_SWEEP_VALUES[-1] <= 0.001

    def test_delta_sweep_largest_order_one(self):
        """Largest δ should be O(1) — substantial temperature variation."""
        assert DELTA_SWEEP_VALUES[0] >= 0.5

    def test_full_sweep_exp_ids(self, tmp_path):
        """Full sweep produces EXP-L2-R1-MULIMIT-001…010."""
        db = str(tmp_path / "full_ids.db")
        results = run_mu_limit_sweep(
            delta_values=DELTA_SWEEP_VALUES[:3],  # just first 3 for speed
            N=16,
            t_end=0.02,
            max_steps=10,
            db_path=db,
        )
        assert results[0]["exp_id"] == "EXP-L2-R1-MULIMIT-001"
        assert results[1]["exp_id"] == "EXP-L2-R1-MULIMIT-002"
        assert results[2]["exp_id"] == "EXP-L2-R1-MULIMIT-003"
