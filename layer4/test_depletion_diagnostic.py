"""
layer4/test_depletion_diagnostic.py
Tests for the H2 vortex-stretching depletion diagnostic
(layer4/depletion_diagnostic.py).

Run: python -m pytest layer4/test_depletion_diagnostic.py -v
"""

from __future__ import annotations

import os
import sqlite3
import tempfile

import numpy as np
import pytest

from layer4.depletion_diagnostic import (
    compute_strain_rate_field,
    alpha_field_from_strain,
    depletion_diagnostics,
    run_depletion_experiment,
    log_result,
    _THETA_EFF_M_FLOOR,
)
from layer4.alignment_bridge import compute_direction_field, _vorticity_from_velocity
from layer4.geometric_disorder import spectral_wavenumbers

N = 16
KX, KY, KZ = spectral_wavenumbers(N)
_X1D = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)


# ─────────────────────────────────────────────────────────────────────────────
# Field constructors
# ─────────────────────────────────────────────────────────────────────────────

def _taylor_green_velocity(N, V0=1.0):
    """TG velocity — div-free by construction (exact incompressible field)."""
    xx, yy, zz = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")
    u = V0 * np.sin(xx) * np.cos(yy) * np.cos(zz)
    v = -V0 * np.cos(xx) * np.sin(yy) * np.cos(zz)
    w = np.zeros((N, N, N))
    return u, v, w


def _abc_velocity(N):
    """u=(sin y, sin z, sin x): an ABC-flow-like periodic field (A=B=C=1).

    At (x,y,z)=(0,0,0) (grid index (0,0,0) under linspace(...,endpoint=False)):
      du/dy=cos y=1, dv/dz=cos z=1, dw/dx=cos x=1, all other first
      derivatives vanish there.
      S_xx=S_yy=S_zz=0; S_xy=S_xz=S_yz=0.5*1=0.5.
      omega = curl u = (dw/dy-dv/dz, du/dz-dw/dx, dv/dx-du/dy)
            = (0-1, 0-1, 0-1) = (-1,-1,-1) at the origin.
      ehat = (-1,-1,-1)/sqrt(3).
      alpha = ehat . S . ehat = 2*(ex*ey*Sxy + ex*ez*Sxz + ey*ez*Syz)
            (diagonal terms vanish) = 2*(1/3*0.5)*3 = 1.0 exactly.
    """
    xx, yy, zz = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")
    u = np.sin(yy)
    v = np.sin(zz)
    w = np.sin(xx)
    return u, v, w


def _abc_vorticity(u, v, w):
    return _vorticity_from_velocity(u, v, w, KX, KY, KZ)


def _twist_omega(N, lam=1.0):
    """omega = lam^2*(cos(lam z), sin(lam z), 0); |omega|=lam^2 everywhere."""
    _, _, zz = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")
    wx = (lam ** 2) * np.cos(lam * zz)
    wy = (lam ** 2) * np.sin(lam * zz)
    wz = np.zeros((N, N, N))
    return wx, wy, wz


# ─────────────────────────────────────────────────────────────────────────────
# compute_strain_rate_field
# ─────────────────────────────────────────────────────────────────────────────

def test_strain_rate_traceless_taylor_green():
    u, v, w = _taylor_green_velocity(N)
    S = compute_strain_rate_field(u, v, w, KX, KY, KZ)
    trace = S["Sxx"] + S["Syy"] + S["Szz"]
    assert np.max(np.abs(trace)) < 1e-10


def test_strain_rate_symmetric_by_independent_recomputation():
    """S_xy as returned must equal an independently re-derived 0.5*(dudy+dvdx)."""
    u, v, w = _abc_velocity(N)
    S = compute_strain_rate_field(u, v, w, KX, KY, KZ)

    def _d(f, k):
        return np.real(np.fft.ifftn(1j * k * np.fft.fftn(f)))

    dudy_indep = _d(u, KY)
    dvdx_indep = _d(v, KX)
    expected_Sxy = 0.5 * (dudy_indep + dvdx_indep)
    assert np.allclose(S["Sxy"], expected_Sxy, atol=1e-10)


def test_strain_rate_field_immutability():
    u, v, w = _taylor_green_velocity(N)
    u0, v0, w0 = u.copy(), v.copy(), w.copy()
    compute_strain_rate_field(u, v, w, KX, KY, KZ)
    assert np.array_equal(u, u0)
    assert np.array_equal(v, v0)
    assert np.array_equal(w, w0)


# ─────────────────────────────────────────────────────────────────────────────
# alpha_field_from_strain — analytic hand-derived value
# ─────────────────────────────────────────────────────────────────────────────

def test_alpha_field_abc_flow_analytic_value_at_origin():
    """See _abc_velocity docstring for the hand derivation: alpha(0,0,0)=1.0."""
    u, v, w = _abc_velocity(N)
    wx, wy, wz = _abc_vorticity(u, v, w)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    S = compute_strain_rate_field(u, v, w, KX, KY, KZ)
    alpha = alpha_field_from_strain(ex, ey, ez, S)
    assert abs(alpha[0, 0, 0] - 1.0) < 1e-8


def test_alpha_field_finite_and_shape_abc_flow():
    """ABC-flow-like periodic field: finiteness/shape sanity over the whole grid."""
    u, v, w = _abc_velocity(N)
    wx, wy, wz = _abc_vorticity(u, v, w)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    S = compute_strain_rate_field(u, v, w, KX, KY, KZ)
    alpha = alpha_field_from_strain(ex, ey, ez, S)
    assert alpha.shape == (N, N, N)
    assert np.all(np.isfinite(alpha))


def test_alpha_field_bounded_by_frobenius_norm_of_strain():
    """|ehat . S . ehat| <= ||S||_2 <= ||S||_F since ehat is a unit vector
    (Cauchy-Schwarz / operator-norm <= Frobenius-norm). Sanity bound on
    alpha independent of the depletion phenomenon itself."""
    u, v, w = _abc_velocity(N)
    wx, wy, wz = _abc_vorticity(u, v, w)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    S = compute_strain_rate_field(u, v, w, KX, KY, KZ)
    alpha = alpha_field_from_strain(ex, ey, ez, S)
    frob = np.sqrt(
        S["Sxx"] ** 2 + S["Syy"] ** 2 + S["Szz"] ** 2
        + 2.0 * S["Sxy"] ** 2 + 2.0 * S["Sxz"] ** 2 + 2.0 * S["Syz"] ** 2
    )
    assert np.all(alpha <= frob + 1e-9)
    assert np.all(-alpha <= frob + 1e-9)


# ─────────────────────────────────────────────────────────────────────────────
# depletion_diagnostics — full dict
# ─────────────────────────────────────────────────────────────────────────────

def test_depletion_diagnostics_keys_and_types():
    u, v, w = _abc_velocity(N)
    wx, wy, wz = _abc_vorticity(u, v, w)
    dx = 2.0 * np.pi / N
    d = depletion_diagnostics(u, v, w, wx, wy, wz, dx)
    for key in ("M", "x_star", "alpha_at_max", "alpha_sup_high",
                "depletion_ratio", "theta_eff", "enstrophy"):
        assert key in d
    assert isinstance(d["M"], float)
    assert isinstance(d["x_star"], tuple) and len(d["x_star"]) == 3
    assert all(isinstance(i, int) for i in d["x_star"])
    assert isinstance(d["depletion_ratio"], float)


def test_depletion_ratio_matches_definition():
    u, v, w = _abc_velocity(N)
    wx, wy, wz = _abc_vorticity(u, v, w)
    dx = 2.0 * np.pi / N
    d = depletion_diagnostics(u, v, w, wx, wy, wz, dx)
    assert d["depletion_ratio"] == pytest.approx(d["alpha_sup_high"] / d["M"], rel=1e-12)


def test_alpha_sup_high_ge_alpha_at_max_or_equal_when_singleton():
    """x* always satisfies |omega(x*)|=M >= M/2, so alpha_sup_high (a max
    over that region) must be >= alpha_at_max."""
    u, v, w = _abc_velocity(N)
    wx, wy, wz = _abc_vorticity(u, v, w)
    dx = 2.0 * np.pi / N
    d = depletion_diagnostics(u, v, w, wx, wy, wz, dx)
    assert d["alpha_sup_high"] >= d["alpha_at_max"] - 1e-12


def test_theta_eff_nan_when_M_below_floor():
    """Scale the twist field so M = lam^2 <= _THETA_EFF_M_FLOOR."""
    lam = 1.0  # M = 1.0 << 10
    wx, wy, wz = _twist_omega(N, lam=lam)
    u, v, w = _taylor_green_velocity(N)  # any smooth div-free velocity works
    dx = 2.0 * np.pi / N
    d = depletion_diagnostics(u, v, w, wx, wy, wz, dx)
    assert d["M"] <= _THETA_EFF_M_FLOOR
    assert np.isnan(d["theta_eff"])


def test_theta_eff_defined_and_matches_formula_when_M_above_floor():
    lam = 4.0  # M = lam^2 = 16 > 10
    wx, wy, wz = _twist_omega(N, lam=lam)
    u, v, w = _taylor_green_velocity(N)
    dx = 2.0 * np.pi / N
    d = depletion_diagnostics(u, v, w, wx, wy, wz, dx)
    assert d["M"] > _THETA_EFF_M_FLOOR
    assert np.isfinite(d["theta_eff"])
    expected = np.log(max(d["alpha_sup_high"], 1e-300)) / np.log(d["M"]) - 0.5
    assert d["theta_eff"] == pytest.approx(expected, rel=1e-10)


def test_depletion_diagnostics_all_zero_raises_valueerror():
    zeros = np.zeros((N, N, N))
    dx = 2.0 * np.pi / N
    with pytest.raises(ValueError):
        depletion_diagnostics(zeros, zeros, zeros, zeros, zeros, zeros, dx)


def test_depletion_diagnostics_immutability():
    u, v, w = _abc_velocity(N)
    wx, wy, wz = _abc_vorticity(u, v, w)
    u0, v0, w0 = u.copy(), v.copy(), w.copy()
    wx0, wy0, wz0 = wx.copy(), wy.copy(), wz.copy()
    dx = 2.0 * np.pi / N
    depletion_diagnostics(u, v, w, wx, wy, wz, dx)
    assert np.array_equal(u, u0) and np.array_equal(v, v0) and np.array_equal(w, w0)
    assert np.array_equal(wx, wx0) and np.array_equal(wy, wy0) and np.array_equal(wz, wz0)


def test_depletion_diagnostics_shape_mismatch_raises_valueerror():
    u, v, w = _abc_velocity(N)
    wx, wy, wz = _abc_vorticity(u, v, w)
    dx = 2.0 * np.pi / N
    with pytest.raises(ValueError):
        depletion_diagnostics(u[:8], v, w, wx, wy, wz, dx)


# ─────────────────────────────────────────────────────────────────────────────
# run_depletion_experiment — smoke test + error paths
# ─────────────────────────────────────────────────────────────────────────────

def test_run_depletion_experiment_invalid_ic_raises():
    with pytest.raises(ValueError):
        run_depletion_experiment(ic="shear", N=16, n_steps=5)


def test_run_depletion_experiment_smoke_N16():
    result = run_depletion_experiment(
        ic="tg", N=16, n_steps=20, record_every=5, seed=42, verbose=False,
    )
    assert result["verdict"] in ("PASS", "FAIL")
    assert len(result["timeseries"]) > 0
    for rec in result["timeseries"]:
        assert np.isfinite(rec["M"])
        assert np.isfinite(rec["depletion_ratio"])
    assert np.isfinite(result["depletion_ratio_median"])


# ─────────────────────────────────────────────────────────────────────────────
# DB logging + retry-on-locked behavior
# ─────────────────────────────────────────────────────────────────────────────

def test_log_result_creates_tables_and_rows():
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "results.db")
        result = run_depletion_experiment(ic="tg", N=16, n_steps=10, record_every=5)
        log_result(db_path, result)

        conn = sqlite3.connect(db_path)
        tables = {r[0] for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )}
        assert "layer4_depletion_experiments" in tables
        assert "layer4_depletion_timeseries" in tables

        row = conn.execute(
            "SELECT exp_id, verdict FROM layer4_depletion_experiments WHERE exp_id=?",
            (result["exp_id"],),
        ).fetchone()
        assert row is not None
        assert row[1] == result["verdict"]

        ts_count = conn.execute(
            "SELECT COUNT(*) FROM layer4_depletion_timeseries WHERE exp_id=?",
            (result["exp_id"],),
        ).fetchone()[0]
        assert ts_count == len(result["timeseries"])
        conn.close()


def test_log_result_retries_on_locked_then_succeeds(monkeypatch):
    """Simulate a transient 'database is locked' error on the first write
    attempt, then let it succeed — verifies the retry/backoff path runs
    without raising and without sleeping for long."""
    import layer4.depletion_diagnostic as dd_mod

    calls = {"n": 0}
    real_write = dd_mod._write_result

    def _flaky_write(db_path, result):
        calls["n"] += 1
        if calls["n"] == 1:
            raise sqlite3.OperationalError("database is locked")
        return real_write(db_path, result)

    monkeypatch.setattr(dd_mod, "_write_result", _flaky_write)

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "results.db")
        result = run_depletion_experiment(ic="tg", N=16, n_steps=10, record_every=5)
        dd_mod.log_result(db_path, result, backoff_base_s=0.001)
        assert calls["n"] == 2  # one failure, one success
