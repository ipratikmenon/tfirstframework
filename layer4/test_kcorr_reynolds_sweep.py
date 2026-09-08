"""
layer4/test_kcorr_reynolds_sweep.py
Tests for the Reynolds sweep driver (layer4/kcorr_reynolds_sweep.py).

Run: python -m pytest layer4/test_kcorr_reynolds_sweep.py -v
"""

from __future__ import annotations

import os
import sqlite3
import tempfile

import numpy as np
import pytest

from layer4.kcorr_reynolds_sweep import (
    build_run_matrix,
    linregress_slope,
    log_sweep_result,
    run_reynolds_sweep,
    run_single_kcorr_sweep,
    spectral_tail_fraction,
    trend_label,
)

N = 32
_X1D = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
_XX, _YY, _ZZ = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")


# ─────────────────────────────────────────────────────────────────────────────
# build_run_matrix
# ─────────────────────────────────────────────────────────────────────────────

class TestBuildRunMatrix:
    def test_size_and_field_values(self):
        matrix = build_run_matrix()
        assert len(matrix) == 6
        assert {row["ic"] for row in matrix} == {"tg", "adv"}
        tg_nus = [row["nu"] for row in matrix if row["ic"] == "tg"]
        assert tg_nus == [1.0e-3, 5.0e-4, 2.5e-4]
        for row in matrix:
            assert row["N"] == 64
            assert row["n_steps"] == 400
            assert row["seed"] == 42
            assert row["record_every"] == 10

    def test_exp_ids_and_claim_ids(self):
        matrix = build_run_matrix()
        tg_ids = [row["exp_id"] for row in matrix if row["ic"] == "tg"]
        adv_ids = [row["exp_id"] for row in matrix if row["ic"] == "adv"]
        assert tg_ids == [
            "EXP-L4-KCORR-NU-TG-001",
            "EXP-L4-KCORR-NU-TG-002",
            "EXP-L4-KCORR-NU-TG-003",
        ]
        assert adv_ids == [
            "EXP-L4-KCORR-NU-ADV-001",
            "EXP-L4-KCORR-NU-ADV-002",
            "EXP-L4-KCORR-NU-ADV-003",
        ]
        # one shared claim_id per ic across its 3 viscosities
        tg_claims = {row["claim_id"] for row in matrix if row["ic"] == "tg"}
        adv_claims = {row["claim_id"] for row in matrix if row["ic"] == "adv"}
        assert tg_claims == {"kcorr-reynolds-tg"}
        assert adv_claims == {"kcorr-reynolds-adv"}

    def test_invalid_ic_raises(self):
        with pytest.raises(ValueError):
            build_run_matrix(ics=("bogus",))

    def test_default_exp_id_start_and_claim_suffix_unchanged(self):
        # Backward-compat: no args -> identical to the original N=64 sweep.
        matrix = build_run_matrix()
        tg_ids = [row["exp_id"] for row in matrix if row["ic"] == "tg"]
        assert tg_ids == [
            "EXP-L4-KCORR-NU-TG-001",
            "EXP-L4-KCORR-NU-TG-002",
            "EXP-L4-KCORR-NU-TG-003",
        ]
        assert {row["claim_id"] for row in matrix if row["ic"] == "tg"} == {"kcorr-reynolds-tg"}

    def test_custom_exp_id_start_and_claim_suffix(self):
        # Used for the N=128 resolved TG re-run: fresh exp_ids/claim_id so
        # it doesn't collide with or overwrite the original N=64 sweep.
        matrix = build_run_matrix(
            ics=("tg",), N=128, exp_id_start=101, claim_id_suffix="-resolved",
        )
        assert len(matrix) == 3
        assert [row["exp_id"] for row in matrix] == [
            "EXP-L4-KCORR-NU-TG-101",
            "EXP-L4-KCORR-NU-TG-102",
            "EXP-L4-KCORR-NU-TG-103",
        ]
        assert {row["claim_id"] for row in matrix} == {"kcorr-reynolds-tg-resolved"}
        assert all(row["N"] == 128 for row in matrix)


# ─────────────────────────────────────────────────────────────────────────────
# spectral_tail_fraction
# ─────────────────────────────────────────────────────────────────────────────

class TestSpectralTailFraction:
    def test_low_mode_field_has_near_zero_tail(self):
        # k=1 only: far below the (2/3)*(N//3) = 6.67 tail threshold at N=32.
        u = np.sin(_XX)
        v = np.zeros((N, N, N))
        w = np.zeros((N, N, N))
        frac = spectral_tail_fraction(u, v, w, N)
        assert 0.0 <= frac < 1e-8

    def test_high_mode_field_has_near_one_tail(self):
        # k=9 only: inside the retained box (N//3=10) and above the tail
        # threshold (2/3 * 10 = 6.67) -> essentially all energy is "tail".
        u = np.cos(9.0 * _XX)
        v = np.zeros((N, N, N))
        w = np.zeros((N, N, N))
        frac = spectral_tail_fraction(u, v, w, N)
        assert frac > 0.99

    def test_zero_field_returns_zero(self):
        z = np.zeros((N, N, N))
        frac = spectral_tail_fraction(z, z, z, N)
        assert frac == 0.0

    def test_fraction_is_bounded(self):
        rng = np.random.default_rng(0)
        u, v, w = rng.standard_normal((3, N, N, N))
        frac = spectral_tail_fraction(u, v, w, N)
        assert 0.0 <= frac <= 1.0 + 1e-12


# ─────────────────────────────────────────────────────────────────────────────
# linregress_slope / trend_label
# ─────────────────────────────────────────────────────────────────────────────

class TestLinregressSlope:
    def test_recovers_known_positive_slope(self):
        x = np.array([0.0, 1.0, 2.0, 3.0])
        y = 2.0 * x + 1.0
        slope, intercept = linregress_slope(x, y)
        assert slope == pytest.approx(2.0, abs=1e-10)
        assert intercept == pytest.approx(1.0, abs=1e-10)

    def test_recovers_known_negative_slope(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = -3.5 * x + 7.0
        slope, intercept = linregress_slope(x, y)
        assert slope == pytest.approx(-3.5, abs=1e-10)
        assert intercept == pytest.approx(7.0, abs=1e-10)

    def test_flat_data_gives_near_zero_slope(self):
        x = np.array([0.0, 1.0, 2.0, 3.0])
        y = np.array([5.0, 5.0, 5.0, 5.0])
        slope, _ = linregress_slope(x, y)
        assert slope == pytest.approx(0.0, abs=1e-10)

    def test_too_few_points_raises(self):
        with pytest.raises(ValueError):
            linregress_slope(np.array([1.0]), np.array([2.0]))

    def test_identical_x_raises(self):
        with pytest.raises(ValueError):
            linregress_slope(np.array([1.0, 1.0, 1.0]), np.array([1.0, 2.0, 3.0]))

    def test_trend_label_classification(self):
        assert trend_label(-1.0) == "decreasing"
        assert trend_label(1.0) == "increasing"
        assert trend_label(0.0) == "flat"
        assert trend_label(float("nan")) == "undefined"


# ─────────────────────────────────────────────────────────────────────────────
# DB logging roundtrip
# ─────────────────────────────────────────────────────────────────────────────

class TestDbRoundtrip:
    def test_log_sweep_result_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_path = os.path.join(tmp, "results.db")
            row = {
                "exp_id": "EXP-L4-KCORR-NU-TG-001",
                "claim_id": "kcorr-reynolds-tg",
                "timestamp": "2026-07-19T00:00:00",
                "ic": "tg",
                "nu": 1.0e-3,
                "N": 64,
                "n_steps": 400,
                "seed": 42,
                "cK_median": 1.05,
                "cK_max": 2.3,
                "cK_final": 1.1,
                "f_ang_median": 0.4,
                "split_residual_median": 0.01,
                "tail_fraction": 5.0e-5,
                "resolution_warning": False,
                "verdict": "PASS",
                "wall_time_s": 12.3,
                "key_metric": "cK_median=1.0500",
            }
            log_sweep_result(db_path, row)

            conn = sqlite3.connect(db_path)
            cur = conn.execute(
                "SELECT exp_id, claim_id, ic, nu, cK_median, resolution_warning, verdict "
                "FROM layer4_kcorr_sweep WHERE exp_id = ?",
                ("EXP-L4-KCORR-NU-TG-001",),
            )
            got = cur.fetchone()
            conn.close()

            assert got is not None
            assert got[0] == "EXP-L4-KCORR-NU-TG-001"
            assert got[1] == "kcorr-reynolds-tg"
            assert got[2] == "tg"
            assert got[3] == pytest.approx(1.0e-3)
            assert got[4] == pytest.approx(1.05)
            assert got[5] == 0  # False stored as 0
            assert got[6] == "PASS"

    def test_log_sweep_result_upsert_replaces_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_path = os.path.join(tmp, "results.db")
            base_row = {
                "exp_id": "EXP-L4-KCORR-NU-ADV-002",
                "claim_id": "kcorr-reynolds-adv",
                "timestamp": "t0",
                "ic": "adv",
                "nu": 5.0e-4,
                "N": 64,
                "n_steps": 400,
                "seed": 42,
                "cK_median": 1.0,
                "cK_max": 1.0,
                "cK_final": 1.0,
                "f_ang_median": 0.3,
                "split_residual_median": 0.01,
                "tail_fraction": 1e-6,
                "resolution_warning": True,
                "verdict": "PASS",
                "wall_time_s": 1.0,
                "key_metric": "cK_median=1.0000",
            }
            log_sweep_result(db_path, base_row)
            updated_row = {**base_row, "cK_median": 3.0, "resolution_warning": False}
            log_sweep_result(db_path, updated_row)

            conn = sqlite3.connect(db_path)
            cur = conn.execute(
                "SELECT COUNT(*), cK_median, resolution_warning FROM layer4_kcorr_sweep "
                "WHERE exp_id = ?",
                ("EXP-L4-KCORR-NU-ADV-002",),
            )
            count, ck_median, res_warn = cur.fetchone()
            conn.close()

            assert count == 1
            assert ck_median == pytest.approx(3.0)
            assert res_warn == 0

    def test_log_sweep_result_retries_then_succeeds_on_lock(self, monkeypatch):
        # Simulate another agent/session holding a write lock on results.db:
        # the first N-1 attempts raise "database is locked", then it clears.
        with tempfile.TemporaryDirectory() as tmp:
            db_path = os.path.join(tmp, "results.db")
            row = {
                "exp_id": "EXP-L4-KCORR-NU-TG-101",
                "claim_id": "kcorr-reynolds-tg-resolved",
                "timestamp": "t0",
                "ic": "tg",
                "nu": 1.0e-3,
                "N": 128,
                "n_steps": 400,
                "seed": 42,
                "cK_median": 1.0,
                "cK_max": 1.0,
                "cK_final": 1.0,
                "f_ang_median": 0.3,
                "split_residual_median": 0.01,
                "tail_fraction": 1e-6,
                "resolution_warning": False,
                "verdict": "PASS",
                "wall_time_s": 1.0,
                "key_metric": "cK_median=1.0000",
            }

            import layer4.kcorr_reynolds_sweep as mod

            calls = {"n": 0}
            real_ensure_db = mod._ensure_db

            def flaky_ensure_db(path):
                calls["n"] += 1
                if calls["n"] < 3:
                    raise sqlite3.OperationalError("database is locked")
                return real_ensure_db(path)

            monkeypatch.setattr(mod, "_ensure_db", flaky_ensure_db)
            monkeypatch.setattr(mod.time, "sleep", lambda _s: None)  # no real waiting in tests

            log_sweep_result(db_path, row, max_retries=6, initial_backoff_s=0.01)
            assert calls["n"] == 3

            conn = sqlite3.connect(db_path)
            cur = conn.execute(
                "SELECT exp_id FROM layer4_kcorr_sweep WHERE exp_id = ?",
                ("EXP-L4-KCORR-NU-TG-101",),
            )
            assert cur.fetchone() is not None
            conn.close()

    def test_log_sweep_result_gives_up_after_max_retries(self, monkeypatch):
        with tempfile.TemporaryDirectory() as tmp:
            db_path = os.path.join(tmp, "results.db")
            row = {
                "exp_id": "EXP-L4-KCORR-NU-TG-102", "claim_id": "kcorr-reynolds-tg-resolved",
                "timestamp": "t0", "ic": "tg", "nu": 5.0e-4, "N": 128, "n_steps": 400,
                "seed": 42, "cK_median": 1.0, "cK_max": 1.0, "cK_final": 1.0,
                "f_ang_median": 0.3, "split_residual_median": 0.01, "tail_fraction": 1e-6,
                "resolution_warning": False, "verdict": "PASS", "wall_time_s": 1.0,
                "key_metric": "cK_median=1.0000",
            }

            import layer4.kcorr_reynolds_sweep as mod

            def always_locked(_path):
                raise sqlite3.OperationalError("database is locked")

            monkeypatch.setattr(mod, "_ensure_db", always_locked)
            monkeypatch.setattr(mod.time, "sleep", lambda _s: None)

            with pytest.raises(sqlite3.OperationalError):
                log_sweep_result(db_path, row, max_retries=2, initial_backoff_s=0.01)


# ─────────────────────────────────────────────────────────────────────────────
# run_single_kcorr_sweep — smoke test + validation
# ─────────────────────────────────────────────────────────────────────────────

class TestRunSingleKcorrSweep:
    def test_smoke_run_N16_tg(self):
        result = run_single_kcorr_sweep(
            ic="tg", nu=1.0e-3, N=16, n_steps=20, seed=42, record_every=5,
        )
        for key in (
            "cK_median", "cK_max", "cK_final", "f_ang_median",
            "split_residual_median", "tail_fraction", "wall_time_s",
        ):
            assert np.isfinite(result[key]), f"{key} not finite: {result[key]}"
        assert isinstance(result["resolution_warning"], bool)
        assert result["verdict"] in ("PASS", "FAIL")
        assert result["N"] == 16
        assert result["n_steps"] == 20
        assert result["nu"] == pytest.approx(1.0e-3)
        assert "cK_median=" in result["key_metric"]

    def test_invalid_ic_raises(self):
        with pytest.raises(ValueError):
            run_single_kcorr_sweep(ic="bogus", nu=1.0e-3, N=16, n_steps=5)

    def test_invalid_n_steps_raises(self):
        with pytest.raises(ValueError):
            run_single_kcorr_sweep(ic="tg", nu=1.0e-3, N=16, n_steps=0)


# ─────────────────────────────────────────────────────────────────────────────
# run_reynolds_sweep — exp_id_start / claim_id_suffix pass-through (end to end)
# ─────────────────────────────────────────────────────────────────────────────

class TestRunReynoldsSweepIdOverrides:
    def test_custom_ids_reach_the_db(self):
        # Tiny smoke sweep (N=16, few steps) exercising the same
        # exp_id_start/claim_id_suffix path used for the N=128 resolved
        # TG re-run, to confirm the override reaches results.db intact.
        with tempfile.TemporaryDirectory() as tmp:
            db_path = os.path.join(tmp, "results.db")
            result = run_reynolds_sweep(
                nus=(1.0e-3, 5.0e-4), ics=("tg",), N=16, n_steps=10, record_every=5,
                verbose=False, db_path=db_path, exp_id_start=101,
                claim_id_suffix="-resolved",
            )
            assert len(result["runs"]) == 2
            row = result["runs"][0]
            assert row["exp_id"] == "EXP-L4-KCORR-NU-TG-101"
            assert row["claim_id"] == "kcorr-reynolds-tg-resolved"
            assert row["N"] == 16
            assert result["runs"][1]["exp_id"] == "EXP-L4-KCORR-NU-TG-102"

            conn = sqlite3.connect(db_path)
            cur = conn.execute(
                "SELECT exp_id, claim_id, N FROM layer4_kcorr_sweep WHERE exp_id = ?",
                ("EXP-L4-KCORR-NU-TG-101",),
            )
            got = cur.fetchone()
            conn.close()
            assert got == ("EXP-L4-KCORR-NU-TG-101", "kcorr-reynolds-tg-resolved", 16)
