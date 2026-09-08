"""
test_lambda_sweep_2D.py — Tests for lambda_sweep_2D.py (M3)

All tests use N=32, t_end=0.05, max_steps=20 (fast) and a temporary DB.
Full 64² runs are done by running lambda_sweep_2D.py directly.

Run with: pytest layer2/test_lambda_sweep_2D.py -v
"""

import os
import sys
import sqlite3
import tempfile

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from lambda_sweep_2D import (
    _ensure_db,
    log_result,
    query_results,
    conjecture_35_verdict,
    run_ccf_sweep,
    run_boussinesq_sweep,
    run_adversarial,
    run_smoke_test,
    BOUS_LAMBDA_VALUES,
)
from T_solver_2D import LAMBDA_SWEEP_VALUES


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture()
def tmp_db(tmp_path):
    """A temporary SQLite database path for each test."""
    return str(tmp_path / "test_results.db")


_FAST = dict(N=32, t_end=0.05, max_steps=20)


# ─────────────────────────────────────────────────────────────────────────────
# Database layer
# ─────────────────────────────────────────────────────────────────────────────

class TestDatabase:

    def test_ensure_db_creates_file(self, tmp_db):
        conn = _ensure_db(tmp_db)
        conn.close()
        assert os.path.exists(tmp_db)

    def test_ensure_db_creates_table(self, tmp_db):
        conn = _ensure_db(tmp_db)
        tables = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()
        conn.close()
        table_names = [t[0] for t in tables]
        assert "experiments" in table_names

    def test_log_result_inserts_row(self, tmp_db):
        conn = _ensure_db(tmp_db)
        row = {
            "exp_id": "EXP-TEST-001", "claim_id": "test_claim",
            "timestamp": "2026-01-01T00:00:00", "layer": 2, "route": 1,
            "profile_type": "CCF", "lambda_val": 1.2, "grid_N": 32,
            "t_end": 0.05, "fluid": "ideal", "n_steps": 10,
            "t_final": 0.05, "verdict": "PASS", "A_min_global": 1e-5,
            "omega_max_global": 0.5, "suppression_ratio_final": 2.0,
            "D_S_ratio_min": None, "T_min_final": 290.0, "T_max_final": 310.0,
            "key_metric": "A_min=1e-5", "notes": "test", "params_json": "{}",
        }
        log_result(conn, row)
        conn.close()
        rows = query_results(db_path=tmp_db)
        assert len(rows) == 1
        assert rows[0]["exp_id"] == "EXP-TEST-001"

    def test_log_result_upserts_on_duplicate_exp_id(self, tmp_db):
        """INSERT OR REPLACE: same exp_id twice should not duplicate."""
        conn = _ensure_db(tmp_db)
        row = {
            "exp_id": "EXP-DUP-001", "claim_id": "x",
            "timestamp": "2026-01-01T00:00:00", "layer": 2, "route": 1,
            "profile_type": "CCF", "lambda_val": 1.2, "grid_N": 32,
            "t_end": 0.05, "fluid": "ideal", "n_steps": 10,
            "t_final": 0.05, "verdict": "PASS", "A_min_global": 1e-5,
            "omega_max_global": 0.5, "suppression_ratio_final": 2.0,
            "D_S_ratio_min": None, "T_min_final": 290.0, "T_max_final": 310.0,
            "key_metric": "v1", "notes": "", "params_json": "{}",
        }
        log_result(conn, row)
        row["key_metric"] = "v2"
        log_result(conn, row)
        conn.close()
        rows = query_results(db_path=tmp_db)
        assert len(rows) == 1, f"Expected 1 row after upsert, got {len(rows)}"

    def test_query_filters_by_layer(self, tmp_db):
        conn = _ensure_db(tmp_db)
        for layer, label in [(2, "L2"), (3, "L3")]:
            log_result(conn, {
                "exp_id": f"EXP-{label}-001", "claim_id": "x",
                "timestamp": "2026-01-01T00:00:00", "layer": layer, "route": 1,
                "profile_type": "CCF", "lambda_val": 1.2, "grid_N": 32,
                "t_end": 0.05, "fluid": "ideal", "n_steps": 5,
                "t_final": 0.05, "verdict": "PASS", "A_min_global": 1e-5,
                "omega_max_global": 0.5, "suppression_ratio_final": 2.0,
                "D_S_ratio_min": None, "T_min_final": 290.0, "T_max_final": 310.0,
                "key_metric": "", "notes": "", "params_json": "{}",
            })
        conn.close()
        rows_l2 = query_results(db_path=tmp_db, layer=2)
        assert len(rows_l2) == 1
        assert rows_l2[0]["exp_id"] == "EXP-L2-001"

    def test_empty_db_returns_empty_list(self, tmp_db):
        rows = query_results(db_path=tmp_db)
        assert rows == []

    def test_nonexistent_db_returns_empty_list(self, tmp_path):
        rows = query_results(db_path=str(tmp_path / "nonexistent.db"))
        assert rows == []


# ─────────────────────────────────────────────────────────────────────────────
# Conjecture 3.5 verdict
# ─────────────────────────────────────────────────────────────────────────────

class TestConjecture35Verdict:

    def test_empty_db_returns_incomplete(self, tmp_db):
        v = conjecture_35_verdict(db_path=tmp_db)
        assert v["verdict"] == "INCOMPLETE"
        assert v["n_runs"] == 0

    def test_all_pass_gives_pass_verdict(self, tmp_db):
        conn = _ensure_db(tmp_db)
        for i in range(3):
            log_result(conn, {
                "exp_id": f"EXP-L2-R1-CCF-{i+1:03d}", "claim_id": "x",
                "timestamp": "2026-01-01T00:00:00", "layer": 2, "route": 1,
                "profile_type": "CCF", "lambda_val": float(i + 1), "grid_N": 32,
                "t_end": 0.05, "fluid": "ideal", "n_steps": 5,
                "t_final": 0.05, "verdict": "PASS", "A_min_global": 1e-5,
                "omega_max_global": 0.5, "suppression_ratio_final": 2.0,
                "D_S_ratio_min": None, "T_min_final": 290.0, "T_max_final": 310.0,
                "key_metric": "", "notes": "", "params_json": "{}",
            })
        conn.close()
        v = conjecture_35_verdict(db_path=tmp_db)
        assert v["verdict"] == "PASS"
        assert v["n_pass"] == 3
        assert v["lambda_c_found"] is False

    def test_one_fail_gives_fail_verdict(self, tmp_db):
        conn = _ensure_db(tmp_db)
        for i, verdict in enumerate(["PASS", "FAIL", "PASS"]):
            log_result(conn, {
                "exp_id": f"EXP-L2-TEST-{i+1:03d}", "claim_id": "x",
                "timestamp": "2026-01-01T00:00:00", "layer": 2, "route": 1,
                "profile_type": "CCF", "lambda_val": float(i + 1), "grid_N": 32,
                "t_end": 0.05, "fluid": "ideal", "n_steps": 5,
                "t_final": 0.05, "verdict": verdict, "A_min_global": 1e-5,
                "omega_max_global": 0.5, "suppression_ratio_final": 2.0,
                "D_S_ratio_min": None, "T_min_final": 290.0, "T_max_final": 310.0,
                "key_metric": "", "notes": "", "params_json": "{}",
            })
        conn.close()
        v = conjecture_35_verdict(db_path=tmp_db)
        assert v["verdict"] == "FAIL"
        assert v["lambda_c_found"] is True


# ─────────────────────────────────────────────────────────────────────────────
# CCF sweep (fast)
# ─────────────────────────────────────────────────────────────────────────────

class TestCCFSweep:

    def test_returns_correct_number_of_results(self, tmp_db):
        lambdas = [1.2, 0.6057]
        results = run_ccf_sweep(lambda_values=lambdas, db_path=tmp_db, **_FAST)
        assert len(results) == 2

    def test_all_ccf_pass_for_canonical_lambdas(self, tmp_db):
        """PRD §7.3: all 6 canonical CCF λ values must PASS."""
        results = run_ccf_sweep(
            lambda_values=LAMBDA_SWEEP_VALUES, db_path=tmp_db, **_FAST
        )
        failures = [r for r in results if r["verdict"] != "PASS"]
        assert len(failures) == 0, (
            f"CCF sweep FAILED for λ: {[r['lambda_val'] for r in failures]}"
        )

    def test_A_min_positive_all_runs(self, tmp_db):
        """A_min_global > 0 for every CCF run — second law."""
        results = run_ccf_sweep(
            lambda_values=LAMBDA_SWEEP_VALUES, db_path=tmp_db, **_FAST
        )
        for r in results:
            assert r["A_min_global"] > 0.0, (
                f"λ={r['lambda_val']}: A_min={r['A_min_global']:.4e} ≤ 0"
            )

    def test_results_logged_to_db(self, tmp_db):
        """Results written to DB; row count matches runs."""
        lambdas = [1.2, 0.1]
        run_ccf_sweep(lambda_values=lambdas, db_path=tmp_db, **_FAST)
        rows = query_results(db_path=tmp_db)
        assert len(rows) == 2

    def test_exp_ids_match_prd(self, tmp_db):
        """Experiment IDs must be EXP-L2-R1-CCF-001…006 for canonical λ values."""
        results = run_ccf_sweep(
            lambda_values=LAMBDA_SWEEP_VALUES, db_path=tmp_db, **_FAST
        )
        for i, r in enumerate(results, start=1):
            assert r["exp_id"] == f"EXP-L2-R1-CCF-{i:03d}", (
                f"Wrong exp_id: {r['exp_id']}, expected EXP-L2-R1-CCF-{i:03d}"
            )

    def test_n_steps_positive(self, tmp_db):
        """Solver must have taken at least 1 step."""
        results = run_ccf_sweep(lambda_values=[0.6057], db_path=tmp_db, **_FAST)
        assert results[0]["n_steps"] >= 1

    def test_critical_lambda_0_4703_passes(self, tmp_db):
        """λ=0.4703 (CCF 2nd unstable) — the hardest canonical case — must PASS."""
        results = run_ccf_sweep(lambda_values=[0.4703], db_path=tmp_db, **_FAST)
        assert results[0]["verdict"] == "PASS", (
            f"Critical λ=0.4703 FAILED: A_min={results[0]['A_min_global']:.4e}"
        )

    def test_small_lambda_0_1_passes(self, tmp_db):
        """λ=0.1 (near adversarial min) must also PASS."""
        results = run_ccf_sweep(lambda_values=[0.1], db_path=tmp_db, **_FAST)
        assert results[0]["verdict"] == "PASS"

    def test_omega_finite(self, tmp_db):
        """Vorticity must remain finite throughout."""
        results = run_ccf_sweep(
            lambda_values=LAMBDA_SWEEP_VALUES, db_path=tmp_db, **_FAST
        )
        for r in results:
            assert np.isfinite(r["omega_max_global"]), (
                f"λ={r['lambda_val']}: omega_max_global is not finite"
            )

    def test_t_final_matches_t_end(self, tmp_db):
        """Solver should reach t_end (within one step dt)."""
        t_end = 0.05
        results = run_ccf_sweep(lambda_values=[1.2], t_end=t_end,
                                max_steps=50, N=32, db_path=tmp_db)
        assert results[0]["t_final"] >= 0.9 * t_end, (
            f"t_final={results[0]['t_final']:.4f} < 0.9*t_end={0.9*t_end:.4f}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Boussinesq sweep (fast)
# ─────────────────────────────────────────────────────────────────────────────

class TestBoussinesqSweep:

    def test_returns_3_results(self, tmp_db):
        results = run_boussinesq_sweep(db_path=tmp_db, **_FAST)
        assert len(results) == 3

    def test_all_boussinesq_pass(self, tmp_db):
        """Conjecture 3.5 must hold for Boussinesq profiles too."""
        results = run_boussinesq_sweep(db_path=tmp_db, **_FAST)
        failures = [r for r in results if r["verdict"] != "PASS"]
        assert len(failures) == 0, (
            f"Boussinesq sweep FAILED for λ: {[r['lambda_val'] for r in failures]}"
        )

    def test_A_min_positive(self, tmp_db):
        results = run_boussinesq_sweep(db_path=tmp_db, **_FAST)
        for r in results:
            assert r["A_min_global"] > 0.0

    def test_exp_ids_bous(self, tmp_db):
        results = run_boussinesq_sweep(db_path=tmp_db, **_FAST)
        for i, r in enumerate(results, start=1):
            assert r["exp_id"] == f"EXP-L2-R1-BOUS-{i:03d}"

    def test_bous_results_in_db(self, tmp_db):
        run_boussinesq_sweep(db_path=tmp_db, **_FAST)
        rows = query_results(db_path=tmp_db)
        bous_rows = [r for r in rows if r["profile_type"] == "Boussinesq"]
        assert len(bous_rows) == 3


# ─────────────────────────────────────────────────────────────────────────────
# Adversarial minimum
# ─────────────────────────────────────────────────────────────────────────────

class TestAdversarial:

    def test_adversarial_passes(self, tmp_db):
        """λ=0.01 adversarial IC must PASS — T-first should still suppress."""
        result = run_adversarial(lambda_val=0.01, db_path=tmp_db, **_FAST)
        assert result["verdict"] == "PASS", (
            f"Adversarial FAIL: A_min={result['A_min_global']:.4e}"
        )

    def test_adversarial_A_min_positive(self, tmp_db):
        result = run_adversarial(lambda_val=0.01, db_path=tmp_db, **_FAST)
        assert result["A_min_global"] > 0.0

    def test_adversarial_exp_id(self, tmp_db):
        result = run_adversarial(lambda_val=0.01, db_path=tmp_db, **_FAST)
        assert result["exp_id"] == "EXP-L2-R1-ADV-001"

    def test_adversarial_logged_to_db(self, tmp_db):
        run_adversarial(lambda_val=0.01, db_path=tmp_db, **_FAST)
        rows = query_results(db_path=tmp_db)
        adv_rows = [r for r in rows if r["profile_type"] == "adversarial_min"]
        assert len(adv_rows) == 1


# ─────────────────────────────────────────────────────────────────────────────
# Smoke test integration
# ─────────────────────────────────────────────────────────────────────────────

class TestSmokeTest:

    def test_smoke_test_passes(self, tmp_db):
        """Built-in smoke test must PASS (3 CCF λ at N=32)."""
        result = run_smoke_test(db_path=tmp_db)
        assert result["verdict"] == "PASS", (
            f"Smoke test FAILED: {result['n_pass']}/{result['n_runs']} pass"
        )

    def test_smoke_test_n_runs(self, tmp_db):
        result = run_smoke_test(db_path=tmp_db)
        assert result["n_runs"] == 3

    def test_smoke_test_all_pass(self, tmp_db):
        result = run_smoke_test(db_path=tmp_db)
        assert result["n_pass"] == result["n_runs"]


# ─────────────────────────────────────────────────────────────────────────────
# End-to-end: full layer 2 M3 verdict from DB
# ─────────────────────────────────────────────────────────────────────────────

class TestM3EndToEnd:

    def test_ccf_plus_bous_gives_pass_verdict(self, tmp_db):
        """After running CCF + Boussinesq, Conjecture 3.5 verdict should be PASS."""
        run_ccf_sweep(lambda_values=LAMBDA_SWEEP_VALUES, db_path=tmp_db, **_FAST)
        run_boussinesq_sweep(db_path=tmp_db, **_FAST)
        v = conjecture_35_verdict(db_path=tmp_db)
        assert v["verdict"] == "PASS", (
            f"Expected PASS, got {v['verdict']}; {v['n_fail']} failures"
        )
        assert v["lambda_c_found"] is False

    def test_total_rows_count(self, tmp_db):
        """CCF (6) + Boussinesq (3) = 9 rows in DB."""
        run_ccf_sweep(lambda_values=LAMBDA_SWEEP_VALUES, db_path=tmp_db, **_FAST)
        run_boussinesq_sweep(db_path=tmp_db, **_FAST)
        rows = query_results(db_path=tmp_db, layer=2)
        assert len(rows) == 9

    def test_A_min_always_positive_across_all_experiments(self, tmp_db):
        """Global A_min over all layer 2 experiments must be > 0."""
        run_ccf_sweep(lambda_values=LAMBDA_SWEEP_VALUES, db_path=tmp_db, **_FAST)
        run_boussinesq_sweep(db_path=tmp_db, **_FAST)
        v = conjecture_35_verdict(db_path=tmp_db)
        assert v["A_min_global"] is not None
        assert v["A_min_global"] > 0.0, (
            f"Global A_min = {v['A_min_global']:.4e} ≤ 0 — second law violation!"
        )
