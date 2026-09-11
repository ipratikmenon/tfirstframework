"""Tests for layer4/helicity_diagnostics.py (WP4b / WP4c).

Calibrated on planted fields with closed-form answers before any evolved flow
is read, and the PARITY structure of the measurement (signed correlation odd
in iota, absolute correlation even) is asserted here rather than assumed --
it is the reason corr(|tau|, alpha) and not corr(tau, alpha) is the number
reported as a mechanism.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "layer3"))

from layer4.helicity_diagnostics import (  # noqa: E402
    MIN_CORR_POINTS, MIN_FIT_POINTS, MIN_WINDOW_SAMPLES, OMEGA_THRESHOLD,
    MIN_TWIST_SIGNAL, N_SURROGATES, phase_randomise, surrogate_null_correlation,
    BOX_VOL, spectral_wavenumbers_3d, vorticity, velocity_gradient,
    strain_from_gradient, stretching_rate, twist_density, helicity_density,
    helicity_total, current_helicity_density, twist_alpha_correlation,
    helicity_snapshot,
    periodic_label, component_topology, detect_reconnection,
    post_transient_time, window_correlations, sweep_stability,
    _pearson, _residualise,
)
from layer4.axisym_diagnostics import abc_flow  # noqa: E402
from twisted_ring_IC import twisted_ring_ic  # noqa: E402

N = 32
N_RING = 48


def _mirror_x_scalar(f):
    return f[(-np.arange(f.shape[0])) % f.shape[0], :, :]


# ── basic operators ──────────────────────────────────────────────────────────

def test_vorticity_of_abc_flow_equals_itself():
    u, v, w = abc_flow(N)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    wx, wy, wz = vorticity(u, v, w, kx, ky, kz)
    assert np.max(np.abs(wx - u)) < 1e-10


def test_strain_tensor_is_symmetric_and_traceless_for_incompressible_flow():
    u, v, w = abc_flow(N)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    S = strain_from_gradient(velocity_gradient(u, v, w, kx, ky, kz))
    for i in range(3):
        for j in range(3):
            assert np.max(np.abs(S[i][j] - S[j][i])) < 1e-12
    tr = S[0][0] + S[1][1] + S[2][2]
    assert np.max(np.abs(tr)) < 1e-10


def test_stretching_rate_of_a_pure_uniaxial_strain_is_exact():
    """u = (ax, -a y/2, -a z/2) has S = diag(a, -a/2, -a/2); with omega along
    x, alpha must equal a exactly."""
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    kx, ky, kz = spectral_wavenumbers_3d(N)
    a = 0.7
    u = a * np.sin(X)
    v = -0.5 * a * np.sin(Y)
    w = -0.5 * a * np.sin(Z)
    wx = np.ones((N, N, N))
    wy = np.zeros((N, N, N))
    wz = np.zeros((N, N, N))
    alpha, valid, mag = stretching_rate(u, v, w, wx, wy, wz, kx, ky, kz)
    assert np.max(np.abs(alpha - a * np.cos(X))) < 1e-9


def test_twist_density_on_a_uniform_helical_direction_field():
    """e = (-sin th sin qz, sin th cos qz, cos th) has e.(curl e) = -q sin^2 th
    exactly, for every point.  The closed-form calibration of tau."""
    q, th = 3.0, 0.6
    z = np.linspace(0, 2 * np.pi, N, endpoint=False)
    Z = np.broadcast_to(z[None, None, :], (N, N, N))
    st, ct = np.sin(th), np.cos(th)
    wx = -st * np.sin(q * Z)
    wy = st * np.cos(q * Z)
    wz = ct * np.ones((N, N, N))
    kx, ky, kz = spectral_wavenumbers_3d(N)
    tau, valid = twist_density(wx, wy, wz, kx, ky, kz)
    assert np.max(np.abs(tau + q * st ** 2)) < 1e-10


def test_twist_density_vanishes_for_a_uniform_direction_field():
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = np.ones((N, N, N))
    kx, ky, kz = spectral_wavenumbers_3d(N)
    tau, _ = twist_density(wx, wy, wz, kx, ky, kz)
    assert np.max(np.abs(tau)) < 1e-12


def test_twist_density_is_odd_under_the_sign_of_q():
    q, th = 2.0, 0.5
    z = np.linspace(0, 2 * np.pi, N, endpoint=False)
    Z = np.broadcast_to(z[None, None, :], (N, N, N))
    kx, ky, kz = spectral_wavenumbers_3d(N)
    out = []
    for qq in (q, -q):
        wx = -np.sin(th) * np.sin(qq * Z)
        wy = np.sin(th) * np.cos(qq * Z)
        wz = np.cos(th) * np.ones((N, N, N))
        out.append(twist_density(wx, wy, wz, kx, ky, kz)[0])
    assert np.max(np.abs(out[0] + out[1])) < 1e-10


def test_twist_density_is_invariant_under_omega_sign_flip():
    """e -> -e leaves e.(curl e) unchanged, so tau does not depend on the
    arbitrary orientation of the vortex line."""
    ic = twisted_ring_ic(N_RING, iota=1.0)
    kx, ky, kz = spectral_wavenumbers_3d(N_RING)
    wx, wy, wz = vorticity(ic["u"], ic["v"], ic["w"], kx, ky, kz)
    a, _ = twist_density(wx, wy, wz, kx, ky, kz)
    b, _ = twist_density(-wx, -wy, -wz, kx, ky, kz)
    assert np.max(np.abs(a - b)) < 1e-10 * np.max(np.abs(a))


# ── helicity ─────────────────────────────────────────────────────────────────

def test_helicity_of_abc_flow_is_twice_its_energy_integral():
    u, v, w = abc_flow(N)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    wx, wy, wz = vorticity(u, v, w, kx, ky, kz)
    H = helicity_total(u, v, w, wx, wy, wz)
    assert H == pytest.approx(float(np.mean(u**2 + v**2 + w**2) * BOX_VOL), rel=1e-9)


def test_helicity_density_is_the_pointwise_dot_product():
    rng = np.random.default_rng(3)
    f = [rng.standard_normal((8, 8, 8)) for _ in range(6)]
    h = helicity_density(*f)
    assert np.allclose(h, f[0]*f[3] + f[1]*f[4] + f[2]*f[5])


# ── PARITY: the reason |tau| and not tau is the measurement ─────────────────

def test_tau_is_a_pseudoscalar_and_alpha_is_a_scalar_under_reflection():
    """Under x -> -x: u is polar (u_x flips), omega is axial (omega_y, omega_z
    flip).  tau must flip sign; alpha must not.  This identity is why
    corr(tau, alpha) is forced to be ODD in iota and therefore cannot itself
    be read as a mechanism."""
    ic = twisted_ring_ic(N_RING, iota=1.2)
    kx, ky, kz = spectral_wavenumbers_3d(N_RING)
    u, v, w = ic["u"], ic["v"], ic["w"]
    wx, wy, wz = vorticity(u, v, w, kx, ky, kz)
    tau, _ = twist_density(wx, wy, wz, kx, ky, kz)
    alpha, _, _ = stretching_rate(u, v, w, wx, wy, wz, kx, ky, kz)

    M = _mirror_x_scalar
    um, vm, wm = -M(u), M(v), M(w)
    wxm, wym, wzm = vorticity(um, vm, wm, kx, ky, kz)
    taum, _ = twist_density(wxm, wym, wzm, kx, ky, kz)
    alpham, _, _ = stretching_rate(um, vm, wm, wxm, wym, wzm, kx, ky, kz)

    s_tau = float(np.max(np.abs(tau)))
    s_al = float(np.max(np.abs(alpha)))
    assert np.max(np.abs(taum + M(tau))) < 1e-8 * s_tau
    assert np.max(np.abs(alpham - M(alpha))) < 1e-8 * s_al


def test_signed_correlation_is_odd_and_absolute_correlation_even_in_iota():
    """The parity identity, measured through the actual correlation routine on
    the +iota / -iota IC pair."""
    kx, ky, kz = spectral_wavenumbers_3d(N_RING)
    res = {}
    for iota in (1.2, -1.2):
        ic = twisted_ring_ic(N_RING, iota=iota)
        res[iota] = helicity_snapshot(ic["u"], ic["v"], ic["w"], kx, ky, kz)
    a, b = res[1.2], res[-1.2]
    assert a["corr_signed"] == pytest.approx(-b["corr_signed"], abs=1e-8)
    assert a["corr_abs"] == pytest.approx(b["corr_abs"], abs=1e-8)
    assert a["H"] == pytest.approx(-b["H"], rel=1e-8)


# ── correlation machinery ────────────────────────────────────────────────────

def test_pearson_is_exact_on_a_linear_relation():
    x = np.linspace(0, 1, 50)
    assert _pearson(x, 3.0 * x + 1.0) == pytest.approx(1.0, abs=1e-12)
    assert _pearson(x, -3.0 * x + 1.0) == pytest.approx(-1.0, abs=1e-12)


def test_pearson_is_nan_on_a_constant():
    assert np.isnan(_pearson(np.ones(10), np.arange(10)))


def test_residualise_removes_a_quadratic_exactly():
    c = np.linspace(-1, 1, 60)
    y = 2.0 + 3.0 * c - 1.5 * c ** 2
    assert np.max(np.abs(_residualise(y, c))) < 1e-10


def test_partial_correlation_kills_a_pure_confound():
    """If |tau| and alpha are both functions of |omega| alone, the raw
    correlation is 1 but the partial correlation must be ~0.  This is the
    core-ness confound the primary number is designed to remove."""
    rng = np.random.default_rng(7)
    n = 4000
    om = rng.uniform(1.0, 2.0, n)
    tau = om ** 2
    alpha = -3.0 * om ** 2
    alpha = alpha + 0.02 * rng.standard_normal(n)   # keep residuals non-degenerate
    assert _pearson(np.abs(tau), alpha) < -0.98
    assert abs(_pearson(_residualise(np.abs(tau), om),
                        _residualise(alpha, om))) < 0.1


def test_correlation_refuses_to_report_below_min_corr_points():
    """S49 discipline: a correlation on too few points is not a result."""
    mag = np.zeros((10, 10, 10))
    mag[0, 0, 0] = 1.0
    valid = np.ones_like(mag, dtype=bool)
    out = twist_alpha_correlation(np.zeros_like(mag), np.zeros_like(mag),
                                  mag, valid)
    assert out["insufficient"]
    assert "corr_abs" not in out
    assert out["n_points"] < MIN_CORR_POINTS


def test_correlation_mask_is_the_quarter_max_vorticity_set():
    rng = np.random.default_rng(11)
    mag = rng.uniform(0.0, 1.0, (24, 24, 24))
    valid = np.ones_like(mag, dtype=bool)
    out = twist_alpha_correlation(rng.standard_normal(mag.shape),
                                  rng.standard_normal(mag.shape), mag, valid)
    assert out["threshold"] == OMEGA_THRESHOLD
    assert out["n_points"] == int(np.count_nonzero(mag >= 0.25 * mag.max()))


def test_correlation_recovers_a_planted_negative_relation():
    rng = np.random.default_rng(13)
    sh = (24, 24, 24)
    mag = rng.uniform(0.5, 1.0, sh)
    valid = np.ones(sh, dtype=bool)
    tau = rng.standard_normal(sh)
    alpha = -2.0 * np.abs(tau) + 0.01 * rng.standard_normal(sh)
    out = twist_alpha_correlation(tau, alpha, mag, valid)
    assert out["corr_abs"] < -0.95


# ── WP4c: topology / reconnection ────────────────────────────────────────────

def test_periodic_label_merges_components_across_the_boundary():
    m = np.zeros((16, 16, 16), dtype=bool)
    m[0, 4:8, 4:8] = True
    m[-1, 4:8, 4:8] = True
    lab_np, n_np = __import__("scipy.ndimage", fromlist=["label"]).label(m)
    lab, n = periodic_label(m)
    assert n_np == 2
    assert n == 1


def test_periodic_label_keeps_genuinely_separate_components_separate():
    m = np.zeros((16, 16, 16), dtype=bool)
    m[2:4, 2:4, 2:4] = True
    m[9:11, 9:11, 9:11] = True
    _, n = periodic_label(m)
    assert n == 2


def test_component_topology_counts_two_blobs_and_measures_their_distance():
    mag = np.zeros((24, 24, 24))
    mag[3, 3, 3] = 1.0
    mag[9, 3, 3] = 1.0
    out = component_topology(mag, level=0.5)
    assert out["n_components"] == 2
    assert out["min_distance"] == pytest.approx(6 * 2 * np.pi / 24, rel=1e-9)


def test_component_topology_single_blob_has_nan_distance():
    mag = np.zeros((16, 16, 16))
    mag[5, 5, 5] = 1.0
    out = component_topology(mag, level=0.5)
    assert out["n_components"] == 1
    assert np.isnan(out["min_distance"])


def test_detect_reconnection_flags_a_component_count_change():
    series = [{"t": float(i), "n_components": 2 if i < 3 else 1,
               "min_component_distance": np.nan, "H": 1.0} for i in range(6)]
    out = detect_reconnection(series)
    assert out["reconnection_detected"]
    assert out["t_reconnect"] == 3.0
    assert "component_count_change" in out["flags"]


def test_detect_reconnection_flags_a_core_distance_minimum():
    d = [1.0, 0.6, 0.3, 0.5, 0.9]
    series = [{"t": float(i), "n_components": 2,
               "min_component_distance": d[i], "H": 1.0} for i in range(5)]
    out = detect_reconnection(series)
    assert out["flags"].get("core_distance_minimum") == 2.0


def test_detect_reconnection_flags_a_helicity_jump():
    H = [1.0, 1.0, 1.0, 0.5, 0.5]
    series = [{"t": float(i), "n_components": 1,
               "min_component_distance": np.nan, "H": H[i]} for i in range(5)]
    out = detect_reconnection(series)
    assert out["flags"].get("helicity_jump") == 3.0


def test_detect_reconnection_is_silent_on_a_quiet_run():
    series = [{"t": float(i), "n_components": 1,
               "min_component_distance": np.nan, "H": 1.0 - 0.001 * i}
              for i in range(10)]
    out = detect_reconnection(series)
    assert not out["reconnection_detected"]
    assert out["t_reconnect"] is None


# ── time splitting ───────────────────────────────────────────────────────────

def test_post_transient_time_skips_an_initial_spike():
    """HANDOVER failure mode 2 (S43): do not read the correlation inside the
    adjustment transient."""
    ts = np.linspace(0, 10, 41)
    M = np.where(ts < 2.0, np.exp(3.0 * ts), np.exp(6.0) * (1 + 0.001 * ts))
    series = [{"t": float(t), "M": float(m)} for t, m in zip(ts, M)]
    assert 0.0 < post_transient_time(series) <= 3.0


def test_window_correlations_refuses_below_min_window_samples():
    series = [{"t": float(i), "corr_signed": 0.1, "corr_abs": -0.1,
               "partial_corr_abs": -0.1, "n_points": 5000,
               "insufficient": False} for i in range(2)]
    out = window_correlations(series, 0.0, 5.0, "w")
    assert out["insufficient"]
    assert out["n_samples"] < MIN_WINDOW_SAMPLES


def test_window_correlations_reports_median_and_sign_consistency():
    series = [{"t": float(i), "corr_signed": 0.1, "corr_abs": -0.2 - 0.01 * i,
               "partial_corr_abs": -0.3, "n_points": 5000,
               "insufficient": False} for i in range(8)]
    out = window_correlations(series, 0.0, 10.0, "w")
    assert not out["insufficient"]
    assert out["corr_abs_median"] < 0
    assert out["corr_abs_sign_consistency"] == 1.0
    assert out["n_samples"] == 8


def test_window_correlations_excludes_insufficient_snapshots():
    series = ([{"t": float(i), "insufficient": True} for i in range(6)]
              + [{"t": 10.0 + i, "corr_signed": 0.1, "corr_abs": -0.2,
                  "partial_corr_abs": -0.2, "n_points": 5000,
                  "insufficient": False} for i in range(5)])
    out = window_correlations(series, 0.0, 30.0, "w")
    assert out["n_samples"] == 5


def test_sweep_stability_detects_a_sign_flip():
    def mk(iota, c):
        return {"iota": iota, "primary": {"insufficient": False,
                                          "corr_abs_median": c,
                                          "partial_corr_abs_median": c,
                                          "corr_signed_median": c}}
    st = sweep_stability([mk(1.0, -0.2), mk(2.0, 0.3)])
    assert st["corr_abs_sign_stable"] is False
    st2 = sweep_stability([mk(1.0, -0.2), mk(2.0, -0.3)])
    assert st2["corr_abs_sign_stable"] is True


def test_sweep_stability_reports_insufficient_below_two_runs():
    assert sweep_stability([])["insufficient"]


# ── snapshot integration ─────────────────────────────────────────────────────

def test_snapshot_on_a_twisted_ring_has_all_fields():
    ic = twisted_ring_ic(N_RING, iota=1.0)
    kx, ky, kz = spectral_wavenumbers_3d(N_RING)
    s = helicity_snapshot(ic["u"], ic["v"], ic["w"], kx, ky, kz)
    for key in ("H", "M", "enstrophy", "E", "tail_fraction", "n_components",
                "n_points", "corr_abs", "partial_corr_abs"):
        assert key in s
    assert s["n_points"] >= MIN_CORR_POINTS
    assert np.isfinite(s["corr_abs"])


def test_snapshot_helicity_matches_the_ic_construction_value():
    ic = twisted_ring_ic(N_RING, iota=1.5)
    kx, ky, kz = spectral_wavenumbers_3d(N_RING)
    s = helicity_snapshot(ic["u"], ic["v"], ic["w"], kx, ky, kz)
    assert s["H"] == pytest.approx(ic["helicity_ic"], rel=1e-9)


def test_untwisted_ring_has_near_zero_twist_density():
    """The mechanism control, and the test that caught the first (wrong)
    implementation of twist_density: an untwisted ring has omega = f(rho)e_phi,
    whose current helicity omega.(curl omega) vanishes identically, so tau must
    be ~0.  Differentiating a truncated e-hat instead gave +/-0.74 here."""
    kx, ky, kz = spectral_wavenumbers_3d(N_RING)
    ic0 = twisted_ring_ic(N_RING, iota=0.0)
    ic1 = twisted_ring_ic(N_RING, iota=1.0)
    w0 = vorticity(ic0["u"], ic0["v"], ic0["w"], kx, ky, kz)
    w1 = vorticity(ic1["u"], ic1["v"], ic1["w"], kx, ky, kz)
    t0 = twist_density(*w0, kx, ky, kz)[0]
    t1 = twist_density(*w1, kx, ky, kz)[0]
    m0 = np.sqrt(sum(c ** 2 for c in w0))
    mask = m0 >= 0.25 * m0.max()
    assert np.mean(np.abs(t0[mask])) < 0.1 * np.mean(np.abs(t1[mask]))


def test_constants_are_the_documented_thresholds():
    assert MIN_CORR_POINTS == 200
    assert MIN_FIT_POINTS == 8
    assert MIN_WINDOW_SAMPLES == 4
    assert OMEGA_THRESHOLD == 0.25


# ── S50: the degenerate-twist guard and the valid null control ──────────────

def test_degenerate_twist_is_refused_not_reported():
    """The iota = 0 failure, as a test.  A Pearson correlation is
    scale-invariant and so cannot say 'no signal'; handed a pure-noise tau it
    returns an arbitrary number (+0.15 at N=64, -0.42 at N=96 in the S50
    sweep).  The guard refuses instead."""
    rng = np.random.default_rng(5)
    sh = (20, 20, 20)
    mag = np.ones(sh)
    valid = np.ones(sh, dtype=bool)
    tiny_tau = 1e-4 * rng.standard_normal(sh)
    alpha = rng.standard_normal(sh)
    out = twist_alpha_correlation(tiny_tau, alpha, mag, valid)
    assert out["degenerate_twist"]
    assert out["insufficient"]
    assert "corr_abs" not in out


def test_real_twist_signal_is_not_refused():
    rng = np.random.default_rng(6)
    sh = (20, 20, 20)
    mag = np.ones(sh)
    valid = np.ones(sh, dtype=bool)
    tau = 1.0 + 0.2 * rng.standard_normal(sh)
    out = twist_alpha_correlation(tau, rng.standard_normal(sh), mag, valid)
    assert not out["degenerate_twist"]
    assert not out["insufficient"]
    assert out["tau_abs_mean_on_mask"] > MIN_TWIST_SIGNAL


def test_untwisted_ring_snapshot_is_flagged_degenerate():
    ic = twisted_ring_ic(N_RING, iota=0.0)
    kx, ky, kz = spectral_wavenumbers_3d(N_RING)
    s = helicity_snapshot(ic["u"], ic["v"], ic["w"], kx, ky, kz)
    assert s["degenerate_twist"]
    assert s["insufficient"]


def test_twisted_ring_snapshot_is_not_flagged_degenerate():
    ic = twisted_ring_ic(N_RING, iota=1.0)
    kx, ky, kz = spectral_wavenumbers_3d(N_RING)
    s = helicity_snapshot(ic["u"], ic["v"], ic["w"], kx, ky, kz)
    assert not s["degenerate_twist"]


def test_phase_randomise_preserves_the_power_spectrum():
    rng = np.random.default_rng(9)
    f = np.sin(np.linspace(0, 4 * np.pi, 16))[:, None, None] * np.ones((16, 16, 16))
    g = phase_randomise(f, rng)
    assert np.max(np.abs(np.abs(np.fft.fftn(g)) - np.abs(np.fft.fftn(f)))) < 1e-8
    assert np.max(np.abs(g.imag if np.iscomplexobj(g) else 0.0)) == 0.0


def test_phase_randomise_destroys_the_pointwise_relation():
    rng = np.random.default_rng(10)
    f = rng.standard_normal((16, 16, 16))
    g = phase_randomise(f, rng)
    assert abs(_pearson(f, g)) < 0.5


def test_surrogate_null_flags_a_planted_strong_correlation_as_significant():
    rng = np.random.default_rng(12)
    sh = (16, 16, 16)
    tau = np.abs(rng.standard_normal(sh)) + 1.0
    alpha = -3.0 * tau + 0.1 * rng.standard_normal(sh)
    mag = np.ones(sh)
    valid = np.ones(sh, dtype=bool)
    out = surrogate_null_correlation(tau, alpha, mag, valid, n_surrogates=60,
                                     seed=1)
    assert not out["insufficient"]
    assert out["observed"] < -0.9
    assert out["outside_null_95"]
    assert abs(out["z_score"]) > 3.0


def test_surrogate_null_calls_an_unrelated_pair_insignificant():
    rng = np.random.default_rng(14)
    sh = (16, 16, 16)
    tau = np.abs(phase_randomise(rng.standard_normal(sh), rng)) + 1.0
    alpha = phase_randomise(rng.standard_normal(sh), rng)
    mag = np.ones(sh)
    valid = np.ones(sh, dtype=bool)
    out = surrogate_null_correlation(tau, alpha, mag, valid, n_surrogates=60,
                                     seed=2)
    assert not out["outside_null_95"]


def test_surrogate_null_is_wider_than_a_naive_permutation_for_smooth_fields():
    """Why phase randomisation and not a permutation: a permutation destroys
    the spatial autocorrelation and so reports a null that is far too narrow
    for smooth fields, making everything look significant."""
    rng = np.random.default_rng(16)
    sh = (24, 24, 24)

    def smooth(x):
        """Low-pass to k <= 3: a genuinely SMOOTH field, which is the regime
        where the two nulls differ.  (A white-spectrum field has no
        autocorrelation for a permutation to destroy, so the two agree.)"""
        F = np.fft.fftn(x)
        k1 = np.fft.fftfreq(sh[0], d=1.0 / sh[0])
        kk = np.sqrt(k1[:, None, None] ** 2 + k1[None, :, None] ** 2
                     + k1[None, None, :] ** 2)
        return np.real(np.fft.ifftn(np.where(kk <= 3.0, F, 0.0)))

    tau = np.abs(smooth(rng.standard_normal(sh))) + 1.0
    alpha = smooth(rng.standard_normal(sh))
    mag = np.ones(sh)
    valid = np.ones(sh, dtype=bool)
    out = surrogate_null_correlation(tau, alpha, mag, valid, n_surrogates=120,
                                     seed=3)
    perm = np.array([_pearson(rng.permutation(tau.ravel()), alpha.ravel())
                     for _ in range(120)])
    assert out["null_std"] > 3.0 * float(np.std(perm))
