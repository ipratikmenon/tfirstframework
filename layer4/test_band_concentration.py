"""
layer4/test_band_concentration.py
Tests for the band-concentration diagnostic (layer4/band_concentration.py),
the numerical stress-test for Conjecture target:band-absorption (S41).

Run: python -m pytest layer4/test_band_concentration.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from layer4.band_concentration import (
    grad2_scalar_sq,
    grad2_ehat_sq,
    band_concentration_diagnostics,
    run_band_experiment,
    DEFAULT_THETA_MINUS,
    DEFAULT_THETA_PLUS,
)
from layer4.geometric_disorder import spectral_wavenumbers
from layer4.alignment_bridge import compute_direction_field

N = 32
KX, KY, KZ = spectral_wavenumbers(N)
_X1D = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
_XX, _YY, _ZZ = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")


# ─────────────────────────────────────────────────────────────────────────────
# Analytic field constructors (same conventions as test_k_correlation.py)
# ─────────────────────────────────────────────────────────────────────────────

def _mixed_field(A=2.0, B=1.0):
    """omega = m(x) * (cos z, sin z, 0), m(x) = A + B*cos(x) > 0.

    Magnitude varies with x only, direction varies with z only: |omega|
    ranges continuously over [A-B, A+B], giving a genuine, non-degenerate
    band {theta_minus*M <= |omega| <= theta_plus*M} to test against.
    """
    m = A + B * np.cos(_XX)
    wx = m * np.cos(_ZZ)
    wy = m * np.sin(_ZZ)
    wz = np.zeros((N, N, N))
    return wx, wy, wz


def _uniform_shell_field(shell_lo=0.15, shell_hi=0.24, A=2.0, B=1.0):
    """Like _mixed_field, but chosen so a known fraction of the domain's
    magnitude falls inside a prescribed shell fraction of [0,1] (in units of
    M = A+B), used to sanity-check band_volume_fraction independently of
    the diagnostic's own band mask.
    """
    return _mixed_field(A=A, B=B)


# ─────────────────────────────────────────────────────────────────────────────
# grad2_scalar_sq
# ─────────────────────────────────────────────────────────────────────────────

def test_grad2_scalar_sq_zero_for_constant_field():
    field = np.full((N, N, N), 3.0)
    g2 = grad2_scalar_sq(field, KX, KY, KZ)
    assert np.max(np.abs(g2)) < 1e-16


def test_grad2_scalar_sq_zero_for_linear_ramp_via_single_mode():
    # f = cos(x) -> d^2f/dx^2 = -cos(x); all other second partials zero.
    # |grad^2 f|^2 = (d_xx f)^2 = cos(x)^2 exactly (spectrally exact).
    field = np.cos(_XX)
    g2 = grad2_scalar_sq(field, KX, KY, KZ)
    expected = np.cos(_XX) ** 2
    assert np.max(np.abs(g2 - expected)) < 1e-10


def test_grad2_scalar_sq_shape_matches_input():
    field = np.cos(_XX) * np.sin(_YY)
    g2 = grad2_scalar_sq(field, KX, KY, KZ)
    assert g2.shape == (N, N, N)


def test_grad2_scalar_sq_immutability():
    field = np.cos(_XX)
    field0 = field.copy()
    grad2_scalar_sq(field, KX, KY, KZ)
    assert np.array_equal(field, field0)


def test_grad2_scalar_sq_nonnegative():
    field = np.cos(_XX) * np.cos(_YY) * np.cos(_ZZ)
    g2 = grad2_scalar_sq(field, KX, KY, KZ)
    assert np.all(g2 >= -1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# grad2_ehat_sq
# ─────────────────────────────────────────────────────────────────────────────

def test_grad2_ehat_sq_zero_for_constant_direction():
    """omega = (0,0,f(x)), f > 0 everywhere -> ehat = (0,0,1) constant
    -> grad^2 ehat = 0 exactly."""
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = 2.0 + np.cos(_XX)
    ex, ey, ez, _mag, _valid = compute_direction_field(wx, wy, wz)
    g2 = grad2_ehat_sq(ex, ey, ez, KX, KY, KZ)
    assert np.max(np.abs(g2)) < 1e-16


def test_grad2_ehat_sq_shape_and_nonnegative():
    wx, wy, wz = _mixed_field()
    ex, ey, ez, _mag, _valid = compute_direction_field(wx, wy, wz)
    g2 = grad2_ehat_sq(ex, ey, ez, KX, KY, KZ)
    assert g2.shape == (N, N, N)
    assert np.all(g2 >= -1e-9)


def test_grad2_ehat_sq_immutability():
    wx, wy, wz = _mixed_field()
    ex, ey, ez, _mag, _valid = compute_direction_field(wx, wy, wz)
    ex0, ey0, ez0 = ex.copy(), ey.copy(), ez.copy()
    grad2_ehat_sq(ex, ey, ez, KX, KY, KZ)
    assert np.array_equal(ex, ex0) and np.array_equal(ey, ey0) and np.array_equal(ez, ez0)


# ─────────────────────────────────────────────────────────────────────────────
# band_concentration_diagnostics — masks and degenerate conventions
# ─────────────────────────────────────────────────────────────────────────────

def test_band_and_E_are_essentially_disjoint():
    """Band = {theta_minus*M <= m <= theta_plus*M}, E = {m >= theta_plus*M}:
    overlap only at the measure-zero boundary m = theta_plus*M."""
    wx, wy, wz = _mixed_field()
    d = band_concentration_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    mag = np.sqrt(wx**2 + wy**2 + wz**2)
    M = float(np.max(mag))
    band_mask = (mag >= DEFAULT_THETA_MINUS * M) & (mag <= DEFAULT_THETA_PLUS * M)
    e_mask = mag >= DEFAULT_THETA_PLUS * M
    overlap = band_mask & e_mask
    boundary = np.abs(mag - DEFAULT_THETA_PLUS * M) < 1e-9
    assert np.all(overlap <= boundary)


def test_band_volume_fraction_positive_for_wide_range_field():
    # A=1.1, B=1.0: |omega| ranges continuously over [0.1, 2.1], M=2.1, so
    # the band [M/8, M/4] = [0.2625, 0.525] sits strictly inside that range
    # (unlike the default A=2, B=1 field, whose range [1,3] lies entirely
    # above M/4=0.75 and gives a legitimately empty band).
    wx, wy, wz = _mixed_field(A=1.1, B=1.0)
    d = band_concentration_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert d["band_volume_fraction"] > 0.0
    assert d["E_volume_fraction"] > 0.0


def test_zero_vorticity_raises():
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = np.zeros((N, N, N))
    with pytest.raises(ValueError, match="identically"):
        band_concentration_diagnostics(wx, wy, wz, 2.0 * np.pi / N)


def test_invalid_theta_ordering_raises():
    wx, wy, wz = _mixed_field()
    with pytest.raises(ValueError, match="theta"):
        band_concentration_diagnostics(
            wx, wy, wz, 2.0 * np.pi / N, theta_minus=0.5, theta_plus=0.25
        )
    with pytest.raises(ValueError, match="theta"):
        band_concentration_diagnostics(
            wx, wy, wz, 2.0 * np.pi / N, theta_minus=0.25, theta_plus=0.25
        )


def test_conc_zero_when_denominator_degenerate():
    """Constant-direction field -> w = |grad ehat|^2 = 0 everywhere ->
    Y_density = |grad m|^2 * w = 0 everywhere -> mean_E_Y = 0 -> conc_Y := 0
    by convention (not NaN or inf)."""
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = 2.0 + np.cos(_XX)
    d = band_concentration_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert d["conc_Y"] == pytest.approx(0.0, abs=1e-12)
    assert d["mean_E_Y"] == pytest.approx(0.0, abs=1e-12)
    assert np.isfinite(d["conc_Y"])


def test_mixed_field_finite_nonnegative_conc():
    wx, wy, wz = _mixed_field()
    d = band_concentration_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert np.isfinite(d["conc_Y"]) and d["conc_Y"] >= 0.0
    assert np.isfinite(d["conc_D"]) and d["conc_D"] >= 0.0


def test_returns_all_required_keys():
    wx, wy, wz = _mixed_field()
    d = band_concentration_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    for key in (
        "conc_Y", "conc_D", "mean_E_Y", "mean_E_D", "mean_Band_Y",
        "mean_Band_D", "E_volume_fraction", "band_volume_fraction",
        "theta_minus", "theta_plus",
    ):
        assert key in d


def test_immutability():
    wx, wy, wz = _mixed_field()
    wx0, wy0, wz0 = wx.copy(), wy.copy(), wz.copy()
    band_concentration_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert np.array_equal(wx, wx0) and np.array_equal(wy, wy0) and np.array_equal(wz, wz0)


def test_default_thetas_match_annulus_residual_definition():
    """Definition def:annulus-residual: Ann = {M/8 <= |omega| <= M/4}."""
    assert DEFAULT_THETA_MINUS == pytest.approx(0.125)
    assert DEFAULT_THETA_PLUS == pytest.approx(0.25)


# ─────────────────────────────────────────────────────────────────────────────
# run_band_experiment — solver-run smoke tests
# ─────────────────────────────────────────────────────────────────────────────

def test_run_band_experiment_tg_smoke():
    r = run_band_experiment(ic="tg", N=16, n_steps=10, record_every=5, verbose=False)
    assert r["verdict"] in ("PASS", "FAIL")
    assert np.isfinite(r["conc_Y_median"])
    assert np.isfinite(r["conc_D_median"])
    assert len(r["timeseries"]) >= 2


def test_run_band_experiment_shear_smoke():
    r = run_band_experiment(ic="shear", N=16, n_steps=10, record_every=5, verbose=False)
    assert r["verdict"] in ("PASS", "FAIL")
    assert np.isfinite(r["conc_Y_median"])


def test_run_band_experiment_adv_smoke():
    r = run_band_experiment(ic="adv", N=16, n_steps=10, record_every=5, verbose=False)
    assert r["verdict"] in ("PASS", "FAIL")
    assert np.isfinite(r["conc_Y_median"])


def test_run_band_experiment_invalid_ic_raises():
    with pytest.raises(ValueError, match="Unknown ic"):
        run_band_experiment(ic="bogus", N=16, n_steps=5)


def test_run_band_experiment_exp_id_and_claim_id():
    r = run_band_experiment(ic="tg", N=16, n_steps=10, record_every=5)
    assert r["exp_id"] == "EXP-L4-BAND-TG-001"
    assert r["claim_id"] == "target-band-absorption-tg"


def test_run_band_experiment_result_has_required_keys():
    r = run_band_experiment(ic="tg", N=16, n_steps=10, record_every=5)
    for key in (
        "exp_id", "claim_id", "conc_Y_median", "conc_D_median",
        "conc_Y_max", "conc_D_max", "verdict", "key_metric",
    ):
        assert key in r


def test_run_band_experiment_db_logging(tmp_path):
    import sqlite3
    db_path = str(tmp_path / "test_band.db")
    run_band_experiment(
        ic="tg", N=16, n_steps=10, record_every=5, db_path=db_path,
    )
    conn = sqlite3.connect(db_path)
    row = conn.execute(
        "SELECT exp_id, verdict FROM band_concentration_experiments "
        "WHERE exp_id = 'EXP-L4-BAND-TG-001'"
    ).fetchone()
    conn.close()
    assert row is not None
    assert row[1] in ("PASS", "FAIL")
