"""
layer4/test_k_correlation.py
Tests for the K-absorption correlation diagnostic (layer4/k_correlation.py).

Run: python -m pytest layer4/test_k_correlation.py -v
"""

from __future__ import annotations

import os
import sqlite3
import tempfile

import numpy as np
import pytest

from layer4.k_correlation import (
    grad_scalar_sq,
    k_correlation_diagnostics,
    kcorr_summary_from_result,
    log_kcorr_result,
    run_and_log_kcorr_experiment,
)
from layer4.geometric_disorder import spectral_wavenumbers
from layer4.alignment_bridge import run_bridge_experiment

N = 32
KX, KY, KZ = spectral_wavenumbers(N)
_X1D = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
_XX, _YY, _ZZ = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")


# ─────────────────────────────────────────────────────────────────────────────
# Analytic field constructors
# ─────────────────────────────────────────────────────────────────────────────

def _twist_field():
    """omega = (cos z, sin z, 0): |omega| = 1 everywhere -> |grad m|^2 = 0."""
    wx = np.cos(_ZZ)
    wy = np.sin(_ZZ)
    wz = np.zeros((N, N, N))
    return wx, wy, wz


def _pure_radial_field(A=2.0, B=1.0):
    """omega = (0, 0, f(x)) with f(x) = A + B*cos(x) > 0: direction constant
    (0,0,1) everywhere -> w = |grad e_hat|^2 = 0 exactly, but magnitude varies.
    """
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = A + B * np.cos(_XX)
    return wx, wy, wz


def _mixed_field(A=2.0, B=1.0):
    """omega = m(x) * (cos z, sin z, 0), m(x) = A + B*cos(x) > 0.
    Magnitude varies with x only, direction varies with z only -> an exact
    orthogonal split (radial and angular gradients live on disjoint axes).
    """
    m = A + B * np.cos(_XX)
    wx = m * np.cos(_ZZ)
    wy = m * np.sin(_ZZ)
    wz = np.zeros((N, N, N))
    return wx, wy, wz


# ─────────────────────────────────────────────────────────────────────────────
# grad_scalar_sq
# ─────────────────────────────────────────────────────────────────────────────

def test_grad_scalar_sq_zero_for_constant_field():
    field = np.full((N, N, N), 3.0)
    g = grad_scalar_sq(field, KX, KY, KZ)
    assert np.max(np.abs(g)) < 1e-20


def test_grad_scalar_sq_shape_matches_input():
    field = np.cos(_XX)
    g = grad_scalar_sq(field, KX, KY, KZ)
    assert g.shape == (N, N, N)


def test_grad_scalar_sq_immutability():
    field = np.cos(_XX)
    field0 = field.copy()
    grad_scalar_sq(field, KX, KY, KZ)
    assert np.array_equal(field, field0)


def test_grad_scalar_sq_matches_analytic_single_mode():
    # f = cos(x) -> |grad f|^2 = sin(x)^2 exactly (spectrally exact, single mode)
    field = np.cos(_XX)
    g = grad_scalar_sq(field, KX, KY, KZ)
    expected = np.sin(_XX) ** 2
    assert np.max(np.abs(g - expected)) < 1e-10


# ─────────────────────────────────────────────────────────────────────────────
# k_correlation_diagnostics — degenerate-covariance conventions
# ─────────────────────────────────────────────────────────────────────────────

def test_k_correlation_twist_field_zero_c_K():
    """m = |omega| is constant -> |grad m|^2 = 0 -> c_K := 0 by convention."""
    wx, wy, wz = _twist_field()
    d = k_correlation_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert d["c_K"] == pytest.approx(0.0, abs=1e-8)
    assert d["mean_grad_m2"] < 1e-15
    assert np.isfinite(d["c_K"])


def test_k_correlation_pure_radial_field_zero_c_K():
    """Direction is constant -> w = |grad e_hat|^2 = 0 -> c_K := 0 by convention."""
    wx, wy, wz = _pure_radial_field()
    d = k_correlation_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert d["c_K"] == pytest.approx(0.0, abs=1e-8)
    assert d["mean_w"] == pytest.approx(0.0, abs=1e-12)
    # magnitude genuinely varies -> mean_grad_m2 is not degenerate
    assert d["mean_grad_m2"] > 0.0
    assert d["f_ang"] == pytest.approx(0.0, abs=1e-8)


def test_k_correlation_mixed_field_finite_positive_c_K_and_small_split_residual():
    wx, wy, wz = _mixed_field()
    d = k_correlation_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert np.isfinite(d["c_K"])
    assert d["c_K"] > 0.0
    assert d["split_residual"] < 0.05
    assert 0.0 <= d["f_ang"] <= 1.0


def test_k_correlation_returns_all_required_keys():
    wx, wy, wz = _mixed_field()
    d = k_correlation_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    for key in (
        "c_K", "f_ang", "mean_grad_m2", "mean_w",
        "E_volume_fraction", "split_residual",
    ):
        assert key in d


def test_k_correlation_f_ang_in_unit_interval_for_various_fields():
    for wx, wy, wz in (_twist_field(), _pure_radial_field(), _mixed_field()):
        d = k_correlation_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
        assert 0.0 <= d["f_ang"] <= 1.0 + 1e-9


def test_k_correlation_e_volume_fraction_full_box_for_these_fields():
    # All three analytic fields have magnitude bounded away from M/4 over the
    # whole box (m in [1,1] or [1,3]), so E should cover the full grid.
    for wx, wy, wz in (_twist_field(), _pure_radial_field(), _mixed_field()):
        d = k_correlation_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
        assert d["E_volume_fraction"] == pytest.approx(1.0, abs=1e-9)


def test_k_correlation_immutability():
    wx, wy, wz = _mixed_field()
    wx0, wy0, wz0 = wx.copy(), wy.copy(), wz.copy()
    k_correlation_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert np.array_equal(wx, wx0)
    assert np.array_equal(wy, wy0)
    assert np.array_equal(wz, wz0)


# ─────────────────────────────────────────────────────────────────────────────
# k_correlation_diagnostics — error handling (shares alignment_bridge's
# _validate_omega_grid contract)
# ─────────────────────────────────────────────────────────────────────────────

def test_k_correlation_all_zero_field_raises_valueerror():
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = np.zeros((N, N, N))
    with pytest.raises(ValueError):
        k_correlation_diagnostics(wx, wy, wz, 2.0 * np.pi / N)


def test_k_correlation_shape_mismatch_raises():
    wx, wy, wz = _mixed_field()
    with pytest.raises(ValueError):
        k_correlation_diagnostics(wx, wy[:-1], wz, 2.0 * np.pi / N)


def test_k_correlation_bad_dx_raises():
    wx, wy, wz = _mixed_field()
    with pytest.raises(ValueError, match="dx"):
        k_correlation_diagnostics(wx, wy, wz, dx=1.0)


# ─────────────────────────────────────────────────────────────────────────────
# Wiring into the experiment runner (alignment_bridge.run_bridge_experiment)
# ─────────────────────────────────────────────────────────────────────────────

def test_run_bridge_experiment_records_c_K_and_f_ang_per_snapshot():
    r = run_bridge_experiment("tg", N=16, n_steps=20, record_every=5, seed=1)
    assert len(r["timeseries"]) > 0
    for rec in r["timeseries"]:
        assert "c_K" in rec and "f_ang" in rec
        assert np.isfinite(rec["c_K"])
        assert np.isfinite(rec["f_ang"])
    assert np.isfinite(r["c_K_max"])
    assert np.isfinite(r["c_K_median"])
    assert np.isfinite(r["c_K_final"])
    assert np.isfinite(r["f_ang_median"])


# ─────────────────────────────────────────────────────────────────────────────
# kcorr_summary_from_result / logging
# ─────────────────────────────────────────────────────────────────────────────

def test_kcorr_summary_from_result_pass_verdict_shape():
    r = run_bridge_experiment("tg", N=16, n_steps=15, record_every=5, seed=2)
    summary = kcorr_summary_from_result(r)
    assert summary["verdict"] in ("PASS", "FAIL")
    assert summary["exp_id"] == "EXP-L4-KCORR-TG-001"
    assert summary["claim_id"] == "k-absorption-corr-tg"
    assert np.isfinite(summary["c_K_median"])
    assert "median_c_K=" in summary["key_metric"]


def test_kcorr_summary_from_result_override_ids():
    r = run_bridge_experiment("shear", N=16, n_steps=10, record_every=5, seed=3)
    summary = kcorr_summary_from_result(
        r, exp_id_override="EXP-L4-KCORR-SH-002", claim_id_override="k-absorption-corr-shear-hi"
    )
    assert summary["exp_id"] == "EXP-L4-KCORR-SH-002"
    assert summary["claim_id"] == "k-absorption-corr-shear-hi"


def test_log_kcorr_result_writes_experiment_and_timeseries_rows():
    r = run_bridge_experiment("tg", N=16, n_steps=10, record_every=5, seed=4)
    summary = kcorr_summary_from_result(r)
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_results.db")
        log_kcorr_result(db_path, summary)

        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT * FROM layer4_kcorr_experiments WHERE exp_id = ?",
            (summary["exp_id"],),
        ).fetchone()
        count = conn.execute(
            "SELECT COUNT(*) FROM layer4_kcorr_timeseries WHERE exp_id = ?",
            (summary["exp_id"],),
        ).fetchone()[0]
        conn.close()

        assert row is not None
        assert row["verdict"] == summary["verdict"]
        assert count == len(summary["timeseries"])


# ─────────────────────────────────────────────────────────────────────────────
# run_and_log_kcorr_experiment — integration smoke test (spec requirement:
# "runner smoke test at N=16 returns c_K/f_ang finite")
# ─────────────────────────────────────────────────────────────────────────────

def test_run_and_log_kcorr_experiment_smoke_n16_finite():
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_results.db")
        summary = run_and_log_kcorr_experiment(
            ic="tg", N=16, n_steps=15, record_every=5, seed=5, db_path=db_path
        )
        assert np.isfinite(summary["c_K_max"])
        assert np.isfinite(summary["c_K_median"])
        assert np.isfinite(summary["c_K_final"])
        assert np.isfinite(summary["f_ang_median"])
        assert summary["verdict"] in ("PASS", "FAIL")

        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT exp_id FROM layer4_kcorr_experiments WHERE exp_id = ?",
            (summary["exp_id"],),
        ).fetchone()
        # Should NOT have written to the bridge tables (kcorr-only runner).
        bridge_tables = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='layer4_bridge_experiments'"
        ).fetchone()
        conn.close()
        assert row is not None
        assert bridge_tables is None
