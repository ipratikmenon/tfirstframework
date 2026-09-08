"""
layer4/test_enstrophy_exponent.py
Tests for the H2 regime-discriminator diagnostic
(layer4/enstrophy_exponent.py): Omega(t) ~ M(t)^p exponent measurement.

Run: python -m pytest layer4/test_enstrophy_exponent.py -v
"""

from __future__ import annotations

import os
import sqlite3
import tempfile

import numpy as np
import pytest

from layer4.enstrophy_exponent import (
    enstrophy_and_M,
    loglog_power_fit,
    select_high_M_mask,
    select_growth_phase_mask,
    run_enstrophy_exponent_experiment,
    log_result,
    _PRIMARY_M_THRESHOLD,
    _FALLBACK_M_THRESHOLD,
)

N = 16
_X1D = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
_DX = 2.0 * np.pi / N


# ─────────────────────────────────────────────────────────────────────────────
# Field constructors
# ─────────────────────────────────────────────────────────────────────────────

def _twist_omega(N, lam=2.0):
    """omega = lam^2*(cos(lam z), sin(lam z), 0); |omega|=lam^2 everywhere
    (constant-magnitude "twist" field: direction rotates with z, magnitude
    is exactly constant across the whole grid — a clean analytic check)."""
    _, _, zz = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")
    wx = (lam ** 2) * np.cos(lam * zz)
    wy = (lam ** 2) * np.sin(lam * zz)
    wz = np.zeros((N, N, N))
    return wx, wy, wz


def _localized_peak_omega(N, base=1.0, peak=9.0):
    """Constant background magnitude `base` everywhere except one grid
    point set to (peak,0,0) — an easy hand-checkable M and Omega."""
    wx = np.full((N, N, N), base)
    wy = np.zeros((N, N, N))
    wz = np.zeros((N, N, N))
    wx[0, 0, 0] = peak
    return wx, wy, wz


# ─────────────────────────────────────────────────────────────────────────────
# enstrophy_and_M — analytic checks
# ─────────────────────────────────────────────────────────────────────────────

def test_enstrophy_and_M_twist_field_matches_hand_calc():
    lam = 2.0
    wx, wy, wz = _twist_omega(N, lam=lam)
    d = enstrophy_and_M(wx, wy, wz, _DX)
    volume = (2.0 * np.pi) ** 3
    expected_Omega = (lam ** 4) * volume   # |omega|^2 = lam^4 everywhere
    assert d["Omega"] == pytest.approx(expected_Omega, rel=1e-10)
    assert d["M"] == pytest.approx(lam ** 2, rel=1e-10)


def test_enstrophy_and_M_localized_peak_matches_hand_calc():
    base, peak = 1.0, 9.0
    wx, wy, wz = _localized_peak_omega(N, base=base, peak=peak)
    d = enstrophy_and_M(wx, wy, wz, _DX)
    volume = (2.0 * np.pi) ** 3
    n_cells = N ** 3
    mean_sq = (base ** 2 * (n_cells - 1) + peak ** 2) / n_cells
    assert d["Omega"] == pytest.approx(mean_sq * volume, rel=1e-10)
    assert d["M"] == pytest.approx(peak, rel=1e-10)


def test_enstrophy_and_M_scales_correctly_under_amplitude_doubling():
    """Omega ~ |omega|^2, so doubling the field quadruples Omega and
    doubles M — a cheap independent sanity check on the formula."""
    wx, wy, wz = _twist_omega(N, lam=1.5)
    d1 = enstrophy_and_M(wx, wy, wz, _DX)
    d2 = enstrophy_and_M(2.0 * wx, 2.0 * wy, 2.0 * wz, _DX)
    assert d2["Omega"] == pytest.approx(4.0 * d1["Omega"], rel=1e-10)
    assert d2["M"] == pytest.approx(2.0 * d1["M"], rel=1e-10)


def test_enstrophy_and_M_bad_dx_raises_valueerror():
    wx, wy, wz = _twist_omega(N)
    with pytest.raises(ValueError):
        enstrophy_and_M(wx, wy, wz, dx=_DX * 1.5)


def test_enstrophy_and_M_shape_mismatch_raises_valueerror():
    wx, wy, wz = _twist_omega(N)
    with pytest.raises(ValueError):
        enstrophy_and_M(wx[:8], wy, wz, _DX)


def test_enstrophy_and_M_immutability():
    wx, wy, wz = _twist_omega(N)
    wx0, wy0, wz0 = wx.copy(), wy.copy(), wz.copy()
    enstrophy_and_M(wx, wy, wz, _DX)
    assert np.array_equal(wx, wx0) and np.array_equal(wy, wy0) and np.array_equal(wz, wz0)


# ─────────────────────────────────────────────────────────────────────────────
# loglog_power_fit — recovers a known synthetic power law
# ─────────────────────────────────────────────────────────────────────────────

def test_loglog_power_fit_recovers_exact_synthetic_exponent():
    """Generate Omega = M^0.5 exactly on a grid of M values; the fit must
    recover p=0.5 (up to floating-point precision) with r_squared ~= 1."""
    M_vals = np.linspace(1.0, 50.0, 30)
    Omega_vals = M_vals ** 0.5
    fit = loglog_power_fit(M_vals, Omega_vals)
    assert fit["p"] == pytest.approx(0.5, abs=1e-8)
    assert fit["r_squared"] == pytest.approx(1.0, abs=1e-6)
    assert fit["n_points"] == 30


def test_loglog_power_fit_recovers_synthetic_exponent_with_noise():
    """Same construction with small log-space noise: p should land close
    to 0.5 and r_squared should stay high but need not be exactly 1."""
    rng = np.random.default_rng(7)
    M_vals = np.linspace(1.0, 50.0, 40)
    log_Omega = 0.5 * np.log(M_vals) + rng.normal(scale=0.01, size=M_vals.shape)
    Omega_vals = np.exp(log_Omega)
    fit = loglog_power_fit(M_vals, Omega_vals)
    assert fit["p"] == pytest.approx(0.5, abs=0.05)
    assert fit["r_squared"] > 0.9


def test_loglog_power_fit_too_few_points_raises():
    with pytest.raises(ValueError):
        loglog_power_fit(np.array([1.0]), np.array([1.0]))


def test_loglog_power_fit_immutability():
    M_vals = np.linspace(1.0, 10.0, 10)
    Omega_vals = M_vals ** 0.3
    M0, O0 = M_vals.copy(), Omega_vals.copy()
    loglog_power_fit(M_vals, Omega_vals)
    assert np.array_equal(M_vals, M0) and np.array_equal(Omega_vals, O0)


# ─────────────────────────────────────────────────────────────────────────────
# select_high_M_mask / select_growth_phase_mask
# ─────────────────────────────────────────────────────────────────────────────

def test_select_high_M_mask_uses_primary_threshold_when_enough_points():
    M_vals = np.array([1.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    threshold, mask, used_fallback = select_high_M_mask(M_vals)
    assert threshold == _PRIMARY_M_THRESHOLD
    assert not used_fallback
    assert np.array_equal(mask, M_vals > _PRIMARY_M_THRESHOLD)


def test_select_high_M_mask_falls_back_when_too_few_points():
    M_vals = np.array([1.0, 2.0, 3.5, 4.0, 6.0])  # only one point > 5
    threshold, mask, used_fallback = select_high_M_mask(M_vals)
    assert threshold == _FALLBACK_M_THRESHOLD
    assert used_fallback
    assert np.array_equal(mask, M_vals > _FALLBACK_M_THRESHOLD)


def test_select_growth_phase_mask_selects_strictly_before_peak():
    M_vals = np.array([1.0, 2.0, 5.0, 9.0, 4.0, 3.0])  # peak at index 3
    mask = select_growth_phase_mask(M_vals)
    assert np.array_equal(mask, np.array([True, True, True, False, False, False]))


def test_select_growth_phase_mask_empty_array():
    assert select_growth_phase_mask(np.array([])).size == 0


# ─────────────────────────────────────────────────────────────────────────────
# run_enstrophy_exponent_experiment — smoke test + error paths
# ─────────────────────────────────────────────────────────────────────────────

def test_run_enstrophy_exponent_experiment_invalid_ic_raises():
    with pytest.raises(ValueError):
        run_enstrophy_exponent_experiment(ic="shear", N=16, n_steps=5)


def test_run_enstrophy_exponent_experiment_smoke_N16():
    result = run_enstrophy_exponent_experiment(
        ic="tg", N=16, n_steps=20, record_every=5, seed=42, verbose=False,
    )
    assert result["verdict"] in ("PASS", "FAIL")
    assert len(result["timeseries"]) > 0
    for rec in result["timeseries"]:
        assert np.isfinite(rec["M"])
        assert np.isfinite(rec["Omega"])
    assert "p_fit_growth_phase" in result and "r_squared" in result


# ─────────────────────────────────────────────────────────────────────────────
# DB logging + retry-on-locked behavior
# ─────────────────────────────────────────────────────────────────────────────

def test_log_result_creates_tables_and_rows():
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "results.db")
        result = run_enstrophy_exponent_experiment(ic="tg", N=16, n_steps=10, record_every=5)
        log_result(db_path, result)

        conn = sqlite3.connect(db_path)
        tables = {r[0] for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )}
        assert "layer4_enstrophy_exponent" in tables
        assert "layer4_enstrophy_exponent_timeseries" in tables

        row = conn.execute(
            "SELECT exp_id, verdict FROM layer4_enstrophy_exponent WHERE exp_id=?",
            (result["exp_id"],),
        ).fetchone()
        assert row is not None
        assert row[1] == result["verdict"]

        ts_count = conn.execute(
            "SELECT COUNT(*) FROM layer4_enstrophy_exponent_timeseries WHERE exp_id=?",
            (result["exp_id"],),
        ).fetchone()[0]
        assert ts_count == len(result["timeseries"])
        conn.close()


def test_log_result_retries_on_locked_then_succeeds(monkeypatch):
    """Simulate a transient 'database is locked' error on the first write
    attempt, then let it succeed — verifies the retry/backoff path runs
    without raising and without sleeping for long."""
    import layer4.enstrophy_exponent as ee_mod

    calls = {"n": 0}
    real_write = ee_mod._write_result

    def _flaky_write(db_path, result):
        calls["n"] += 1
        if calls["n"] == 1:
            raise sqlite3.OperationalError("database is locked")
        return real_write(db_path, result)

    monkeypatch.setattr(ee_mod, "_write_result", _flaky_write)

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "results.db")
        result = run_enstrophy_exponent_experiment(ic="tg", N=16, n_steps=10, record_every=5)
        ee_mod.log_result(db_path, result, backoff_base_s=0.001)
        assert calls["n"] == 2  # one failure, one success
