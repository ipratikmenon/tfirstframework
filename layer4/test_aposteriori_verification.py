"""
layer4/test_aposteriori_verification.py
Tests for the a posteriori regularity VERIFICATION DIAGNOSTIC
(layer4/aposteriori_verification.py) — NOT a proof of regularity; see the
module docstring's SCOPE section, restated in these tests where relevant.

Run: python -m pytest layer4/test_aposteriori_verification.py -v
"""

from __future__ import annotations

import os
import sqlite3
import tempfile

import numpy as np
import pytest

from layer4.aposteriori_verification import (
    ns_residual,
    estimate_dudt_time_diff,
    residual_norms,
    amplification_rate,
    default_threshold,
    verification_margin,
    windowed_verification,
    windowed_summary,
    run_verification_experiment,
    log_result,
)
import route2_3D as r2  # noqa: E402 — layer3 solver (path already set up by module under test)

N = 16
_X1D = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
_KX, _KY, _KZ, _K2 = r2.make_grid_3D(N)
_DX = 2.0 * np.pi / N


# ─────────────────────────────────────────────────────────────────────────────
# Exact NS solution: 2D Taylor-Green vortex decay embedded in 3D
# u = e^{-2νt}(sin x cos y, -cos x sin y, 0) — an EXACT solution of the full
# (nonlinear) incompressible NS equations, not merely of the linear Stokes
# problem (the matching pressure field cancels the nonlinear term's gradient
# part under Leray projection).
# ─────────────────────────────────────────────────────────────────────────────

def _tg_stokes_velocity(t: float, nu: float):
    xx, yy, zz = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")
    decay = np.exp(-2.0 * nu * t)
    u = decay * np.sin(xx) * np.cos(yy)
    v = -decay * np.cos(xx) * np.sin(yy)
    w = np.zeros((N, N, N))
    return u, v, w


def _tg_stokes_dudt_exact(t: float, nu: float):
    """Analytic ∂_t u_exact = -2ν · u_exact (since u = e^{-2νt}·f(x,y))."""
    u, v, w = _tg_stokes_velocity(t, nu)
    return -2.0 * nu * u, -2.0 * nu * v, -2.0 * nu * w


# ─────────────────────────────────────────────────────────────────────────────
# ns_residual — exact solution vs. perturbed field
# ─────────────────────────────────────────────────────────────────────────────

def test_ns_residual_near_zero_on_exact_solution():
    nu = 0.05
    t = 0.3
    u, v, w = _tg_stokes_velocity(t, nu)
    dudt = _tg_stokes_dudt_exact(t, nu)

    rx, ry, rz = ns_residual(u, v, w, None, nu, _KX, _KY, _KZ, dudt)
    l2 = float(np.sqrt(np.mean(rx**2 + ry**2 + rz**2)))
    # Spectral accuracy on a smooth exact solution: residual should be at
    # machine-precision scale (documented tolerance, not exactly zero due
    # to floating point).
    assert l2 < 1e-9


def test_ns_residual_grows_under_perturbation():
    nu = 0.05
    t = 0.3
    u, v, w = _tg_stokes_velocity(t, nu)
    dudt = _tg_stokes_dudt_exact(t, nu)
    rx0, ry0, rz0 = ns_residual(u, v, w, None, nu, _KX, _KY, _KZ, dudt)
    l2_exact = float(np.sqrt(np.mean(rx0**2 + ry0**2 + rz0**2)))

    rng = np.random.default_rng(0)
    bump = 0.5 * rng.standard_normal(u.shape)
    u_pert = u + bump
    rx1, ry1, rz1 = ns_residual(u_pert, v, w, None, nu, _KX, _KY, _KZ, dudt)
    l2_pert = float(np.sqrt(np.mean(rx1**2 + ry1**2 + rz1**2)))

    assert l2_pert > 100.0 * max(l2_exact, 1e-14)


def test_ns_residual_shape_mismatch_raises():
    u, v, w = _tg_stokes_velocity(0.0, 0.05)
    with pytest.raises(ValueError):
        ns_residual(u, v[:-1], w, None, 0.05, _KX, _KY, _KZ, (u, v, w))


def test_ns_residual_immutable():
    u, v, w = _tg_stokes_velocity(0.2, 0.05)
    u_copy, v_copy, w_copy = u.copy(), v.copy(), w.copy()
    dudt = _tg_stokes_dudt_exact(0.2, 0.05)
    dudt_copy = tuple(a.copy() for a in dudt)
    ns_residual(u, v, w, None, 0.05, _KX, _KY, _KZ, dudt)
    assert np.array_equal(u, u_copy) and np.array_equal(v, v_copy) and np.array_equal(w, w_copy)
    for a, ac in zip(dudt, dudt_copy):
        assert np.array_equal(a, ac)


# ─────────────────────────────────────────────────────────────────────────────
# estimate_dudt_time_diff
# ─────────────────────────────────────────────────────────────────────────────

def test_estimate_dudt_time_diff_basic():
    u0, v0, w0 = _tg_stokes_velocity(0.0, 0.05)
    u1, v1, w1 = _tg_stokes_velocity(0.01, 0.05)
    dudt = estimate_dudt_time_diff(u0, v0, w0, u1, v1, w1, 0.01)
    expected_u = (u1 - u0) / 0.01
    assert np.allclose(dudt[0], expected_u)


def test_estimate_dudt_time_diff_rejects_nonpositive_dt():
    u0, v0, w0 = _tg_stokes_velocity(0.0, 0.05)
    with pytest.raises(ValueError):
        estimate_dudt_time_diff(u0, v0, w0, u0, v0, w0, 0.0)


# ─────────────────────────────────────────────────────────────────────────────
# residual_norms
# ─────────────────────────────────────────────────────────────────────────────

def test_residual_norms_zero_on_zero_field():
    z = np.zeros((N, N, N))
    norms = residual_norms(z, z, z, _DX)
    assert norms["L2"] == 0.0
    assert norms["Hm1"] == 0.0


def test_residual_norms_ordering_hm1_le_l2():
    rng = np.random.default_rng(1)
    # Mid-frequency field: a handful of non-DC Fourier modes.
    f_hat = np.zeros((N, N, N), dtype=complex)
    mid = (np.abs(_KX) > 0) & (np.abs(_KX) <= 4) & (np.abs(_KY) <= 4) & (np.abs(_KZ) <= 4)
    f_hat[mid] = rng.standard_normal(int(np.sum(mid))) + 1j * rng.standard_normal(int(np.sum(mid)))
    rx = np.real(np.fft.ifftn(f_hat))
    ry = np.real(np.fft.ifftn(f_hat)) * 0.5
    rz = np.zeros((N, N, N))
    norms = residual_norms(rx, ry, rz, _DX)
    assert norms["Hm1"] <= norms["L2"] + 1e-12
    assert norms["L2"] > 0.0


def test_residual_norms_bad_dx_raises():
    z = np.zeros((N, N, N))
    with pytest.raises(ValueError):
        residual_norms(z, z, z, dx=1.0)  # inconsistent with N=16 -> 2pi/16


def test_residual_norms_immutable():
    rx = np.random.default_rng(2).standard_normal((N, N, N))
    ry, rz = rx.copy(), rx.copy()
    rx_c, ry_c, rz_c = rx.copy(), ry.copy(), rz.copy()
    residual_norms(rx, ry, rz, _DX)
    assert np.array_equal(rx, rx_c) and np.array_equal(ry, ry_c) and np.array_equal(rz, rz_c)


# ─────────────────────────────────────────────────────────────────────────────
# amplification_rate / default_threshold
# ─────────────────────────────────────────────────────────────────────────────

def test_amplification_rate_zero_on_uniform_field():
    ones = np.ones((N, N, N))
    z = np.zeros((N, N, N))
    a = amplification_rate(ones, z, z, _KX, _KY, _KZ, C=1.0)
    assert a < 1e-10


def test_amplification_rate_scales_with_C():
    u, v, w = _tg_stokes_velocity(0.1, 0.05)
    a1 = amplification_rate(u, v, w, _KX, _KY, _KZ, C=1.0)
    a2 = amplification_rate(u, v, w, _KX, _KY, _KZ, C=2.0)
    assert abs(a2 - 2.0 * a1) < 1e-10


def test_default_threshold_positive_and_scales_with_nu():
    assert default_threshold(1e-3) == pytest.approx(1e-3)
    assert default_threshold(1e-3, C_thresh=2.0) == pytest.approx(2e-3)
    with pytest.raises(ValueError):
        default_threshold(0.0)


# ─────────────────────────────────────────────────────────────────────────────
# verification_margin — the correctness anchors
# ─────────────────────────────────────────────────────────────────────────────

def test_verification_margin_A_zero_exact_quadrature():
    T = 2.0
    n = 21
    times = np.linspace(0.0, T, n)
    r0 = 0.37
    residual_series = np.full(n, r0)
    A_series = np.zeros(n)
    out = verification_margin(times, residual_series, A_series, threshold=1.0)
    assert out["E_final"] == pytest.approx(r0 * T, rel=1e-12)


def test_verification_margin_A_constant_exact_gronwall():
    T = 1.5
    n = 31
    times = np.linspace(0.0, T, n)
    r0 = 0.2
    a = 0.8
    residual_series = np.full(n, r0)
    A_series = np.full(n, a)
    out = verification_margin(times, residual_series, A_series, threshold=1.0)
    expected = r0 * (np.exp(a * T) - 1.0) / a
    assert out["E_final"] == pytest.approx(expected, rel=1e-10)


def test_verification_margin_monotonic_doubling_residual():
    T = 1.0
    n = 11
    times = np.linspace(0.0, T, n)
    A_series = np.zeros(n)
    out1 = verification_margin(times, np.full(n, 0.1), A_series, threshold=1.0)
    out2 = verification_margin(times, np.full(n, 0.2), A_series, threshold=1.0)
    assert out2["E_final"] == pytest.approx(2.0 * out1["E_final"], rel=1e-12)


def test_verification_margin_ratio_reflects_threshold():
    times = np.linspace(0.0, 1.0, 5)
    A_series = np.zeros(5)
    r_series = np.full(5, 1.0)
    out = verification_margin(times, r_series, A_series, threshold=2.0)
    assert out["margin_final"] == pytest.approx(out["E_final"] / 2.0)


def test_verification_margin_degenerate_single_point_raises():
    with pytest.raises(ValueError):
        verification_margin(np.array([0.0]), np.array([0.0]), np.array([0.0]), threshold=1.0)


def test_verification_margin_mismatched_lengths_raises():
    with pytest.raises(ValueError):
        verification_margin(np.array([0.0, 1.0]), np.array([0.0]), np.array([0.0, 1.0]), threshold=1.0)


def test_verification_margin_nonpositive_threshold_raises():
    times = np.array([0.0, 1.0])
    with pytest.raises(ValueError):
        verification_margin(times, np.zeros(2), np.zeros(2), threshold=0.0)


def test_verification_margin_nonincreasing_times_raises():
    with pytest.raises(ValueError):
        verification_margin(np.array([0.0, 0.0]), np.zeros(2), np.zeros(2), threshold=1.0)


def test_verification_margin_immutable():
    times = np.linspace(0.0, 1.0, 5)
    r_series = np.full(5, 0.3)
    A_series = np.full(5, 0.1)
    times_c, r_c, a_c = times.copy(), r_series.copy(), A_series.copy()
    verification_margin(times, r_series, A_series, threshold=1.0)
    assert np.array_equal(times, times_c)
    assert np.array_equal(r_series, r_c)
    assert np.array_equal(A_series, a_c)


# ─────────────────────────────────────────────────────────────────────────────
# windowed_verification / windowed_summary — the fix for the long-window
# margin's known vacuity (see module SCOPE note). Each window independently
# resets E=0 at its own start; reuses verification_margin's own recursion.
# ─────────────────────────────────────────────────────────────────────────────

def test_windowed_verification_A_zero_exact_per_window():
    # T=4.0, dt=0.1 (n=41 points), W=1.0 divides T exactly -> 4 windows,
    # each of true duration 1.0, so E_window == r0 * W exactly (A=0 anchor,
    # same closed form as the long-window test, applied per window).
    T = 4.0
    n = 41
    times = np.linspace(0.0, T, n)
    r0 = 0.37
    residual_series = np.full(n, r0)
    A_series = np.zeros(n)
    windows = windowed_verification(times, residual_series, A_series, window_length=1.0, threshold=1.0)
    assert len(windows) == 4
    for w in windows:
        assert w["t_end"] - w["t_start"] == pytest.approx(1.0, rel=1e-12)
        assert w["E_window"] == pytest.approx(r0 * 1.0, rel=1e-10)
        assert w["gronwall_factor"] == pytest.approx(1.0, rel=1e-12)


def test_windowed_verification_A_constant_exact_gronwall_per_window():
    # Same T/dt/W setup but A ≡ a > 0: each window's E_window should match
    # the closed-form Gronwall solution r0*(exp(a*W)-1)/a EXACTLY (to
    # floating tolerance) — the constant-coefficient anchor, per window.
    T = 4.0
    n = 41
    times = np.linspace(0.0, T, n)
    r0 = 0.2
    a = 0.8
    W = 1.0
    residual_series = np.full(n, r0)
    A_series = np.full(n, a)
    windows = windowed_verification(times, residual_series, A_series, window_length=W, threshold=1.0)
    expected_E = r0 * (np.exp(a * W) - 1.0) / a
    expected_gronwall = np.exp(a * W)
    assert len(windows) == 4
    for w in windows:
        assert w["E_window"] == pytest.approx(expected_E, rel=1e-10)
        assert w["gronwall_factor"] == pytest.approx(expected_gronwall, rel=1e-10)


def test_windowed_verification_partial_trailing_window():
    # T=2.2, W=0.5 -> 4 full windows + one partial trailing window of
    # duration 0.2 (documented behavior: returned as-is, not dropped/padded).
    T = 2.2
    n = 23  # dt = 0.1, so windows land on exact sample points
    times = np.linspace(0.0, T, n)
    r0 = 0.5
    A_series = np.zeros(n)
    windows = windowed_verification(times, np.full(n, r0), A_series, window_length=0.5, threshold=1.0)
    assert len(windows) == 5
    for w in windows[:4]:
        assert w["t_end"] - w["t_start"] == pytest.approx(0.5, rel=1e-9)
    last = windows[-1]
    assert last["t_end"] - last["t_start"] == pytest.approx(0.2, rel=1e-9)
    assert last["t_end"] == pytest.approx(T, rel=1e-12)
    # partial window's E is still r0 * its own (shorter) duration under A=0
    assert last["E_window"] == pytest.approx(r0 * 0.2, rel=1e-9)


def test_windowed_verification_satisfied_flips_around_threshold():
    times = np.linspace(0.0, 1.0, 11)
    A_series = np.zeros(11)
    r0 = 0.5
    # E_window over the single full window [0,1] with A=0 is r0*1.0 = 0.5
    windows_ok = windowed_verification(times, np.full(11, r0), A_series, window_length=1.0, threshold=1.0)
    assert len(windows_ok) == 1
    assert windows_ok[0]["margin_window"] == pytest.approx(0.5)
    assert windows_ok[0]["satisfied"] is True

    windows_bad = windowed_verification(times, np.full(11, r0), A_series, window_length=1.0, threshold=0.4)
    assert windows_bad[0]["margin_window"] == pytest.approx(0.5 / 0.4)
    assert windows_bad[0]["satisfied"] is False


def test_windowed_verification_window_count_matches_ceil_partition():
    T = 3.0
    n = 31  # dt=0.1
    times = np.linspace(0.0, T, n)
    windows = windowed_verification(times, np.full(n, 0.1), np.zeros(n), window_length=0.7, threshold=1.0)
    # ceil(3.0/0.7) = 5 windows (last one partial)
    assert len(windows) == 5
    assert windows[-1]["t_end"] == pytest.approx(T, rel=1e-9)


def test_windowed_verification_too_few_points_per_window_raises():
    # window_length far smaller than the sampling cadence -> every window
    # would hold at most 1 sample point -> degenerate, ValueError.
    times = np.linspace(0.0, 1.0, 3)  # dt = 0.5
    with pytest.raises(ValueError):
        windowed_verification(times, np.zeros(3), np.zeros(3), window_length=0.01, threshold=1.0)


def test_windowed_verification_nonpositive_window_length_raises():
    times = np.linspace(0.0, 1.0, 5)
    with pytest.raises(ValueError):
        windowed_verification(times, np.zeros(5), np.zeros(5), window_length=0.0, threshold=1.0)


def test_windowed_verification_mismatched_lengths_raises():
    with pytest.raises(ValueError):
        windowed_verification(
            np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0]), np.array([0.0, 1.0, 2.0]),
            window_length=0.5, threshold=1.0,
        )


def test_windowed_summary_aggregates_correctly():
    windows = [
        {"t_start": 0.0, "t_end": 1.0, "gronwall_factor": 1.5, "E_window": 0.3, "margin_window": 0.3, "satisfied": True},
        {"t_start": 1.0, "t_end": 2.0, "gronwall_factor": 2.5, "E_window": 1.2, "margin_window": 1.2, "satisfied": False},
        {"t_start": 2.0, "t_end": 3.0, "gronwall_factor": 1.1, "E_window": 0.05, "margin_window": 0.05, "satisfied": True},
    ]
    summary = windowed_summary(windows)
    assert summary["n_windows"] == 3
    assert summary["n_satisfied"] == 2
    assert summary["fraction_satisfied"] == pytest.approx(2.0 / 3.0)
    assert summary["max_gronwall_factor"] == pytest.approx(2.5)
    assert summary["max_margin"] == pytest.approx(1.2)
    assert summary["median_margin"] == pytest.approx(0.3)
    assert summary["worst_window_t_start"] == pytest.approx(1.0)


def test_windowed_summary_empty_raises():
    with pytest.raises(ValueError):
        windowed_summary([])


def test_run_verification_experiment_with_window_lengths_populates_windowed():
    result = run_verification_experiment(
        ic="tg", N=16, n_steps=20, record_every=1, seed=42, verbose=False,
        window_lengths=(0.05, 0.1),
    )
    assert "windowed" in result
    assert set(result["windowed"].keys()) == {0.05, 0.1}
    for w_len, w_data in result["windowed"].items():
        assert "windows" in w_data and "summary" in w_data
        assert w_data["summary"]["n_windows"] == len(w_data["windows"])


# ─────────────────────────────────────────────────────────────────────────────
# run_verification_experiment — smoke test + error paths
# ─────────────────────────────────────────────────────────────────────────────

def test_run_verification_experiment_invalid_ic_raises():
    with pytest.raises(ValueError):
        run_verification_experiment(ic="bogus", N=16, n_steps=5)


def test_run_verification_experiment_invalid_n_steps_raises():
    with pytest.raises(ValueError):
        run_verification_experiment(ic="tg", N=16, n_steps=0)


def test_run_verification_experiment_smoke_N16():
    result = run_verification_experiment(
        ic="tg", N=16, n_steps=20, record_every=5, seed=42, verbose=False,
    )
    assert result["verdict"] in ("PASS", "FAIL", "PARTIAL")
    assert result["exp_id"] == "EXP-L4-VERIFY-TG-001"
    assert result["claim_id"] == "aposteriori-tg"
    assert len(result["timeseries"]) > 0
    for rec in result["timeseries"]:
        assert np.isfinite(rec["residual_L2"])
        assert np.isfinite(rec["residual_Hm1"])
        assert np.isfinite(rec["A_t"])
        assert np.isfinite(rec["E_t"])
        assert np.isfinite(rec["margin_t"])
    for key in ("residual_L2_min", "residual_L2_max", "residual_Hm1_min",
                "residual_Hm1_max", "A_max", "E_final", "margin_final",
                "margin_max", "diagnostic_satisfied", "wall_time_s", "key_metric"):
        assert key in result


def test_run_verification_experiment_shear_smoke_N16():
    result = run_verification_experiment(
        ic="shear", N=16, n_steps=20, record_every=5, seed=42, verbose=False,
    )
    assert result["verdict"] in ("PASS", "FAIL", "PARTIAL")
    assert len(result["timeseries"]) > 0


# ─────────────────────────────────────────────────────────────────────────────
# DB logging + retry-on-locked behavior
# ─────────────────────────────────────────────────────────────────────────────

def test_log_result_creates_tables_and_rows():
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "results.db")
        result = run_verification_experiment(ic="tg", N=16, n_steps=10, record_every=5)
        log_result(db_path, result)

        conn = sqlite3.connect(db_path)
        tables = {r[0] for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )}
        assert "layer4_aposteriori_experiments" in tables
        assert "layer4_aposteriori_timeseries" in tables

        row = conn.execute(
            "SELECT exp_id, verdict FROM layer4_aposteriori_experiments WHERE exp_id=?",
            (result["exp_id"],),
        ).fetchone()
        assert row is not None
        assert row[1] == result["verdict"]

        ts_count = conn.execute(
            "SELECT COUNT(*) FROM layer4_aposteriori_timeseries WHERE exp_id=?",
            (result["exp_id"],),
        ).fetchone()[0]
        assert ts_count == len(result["timeseries"])
        conn.close()


def test_log_result_retries_on_locked_then_succeeds(monkeypatch):
    """Simulate a transient 'database is locked' error on the first write
    attempt, then let it succeed — verifies the retry/backoff path runs
    without raising and without sleeping for long."""
    import layer4.aposteriori_verification as av_mod

    calls = {"n": 0}
    real_write = av_mod._write_result

    def _flaky_write(db_path, result):
        calls["n"] += 1
        if calls["n"] == 1:
            raise sqlite3.OperationalError("database is locked")
        return real_write(db_path, result)

    monkeypatch.setattr(av_mod, "_write_result", _flaky_write)

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "results.db")
        result = run_verification_experiment(ic="tg", N=16, n_steps=10, record_every=5)
        av_mod.log_result(db_path, result, backoff_base_s=0.001)
        assert calls["n"] == 2  # one failure, one success
