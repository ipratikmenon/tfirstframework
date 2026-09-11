"""Tests for layer4/axisym_diagnostics.py (WP4a).

Calibration first, per the house rule: the diagnostics are checked against
closed forms (exact Burgers vortex; a planted field with a KNOWN NON-ZERO
source) before anything reads them off an evolved flow.
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

from layer4.axisym_diagnostics import (  # noqa: E402
    AXIS_EPS_CELLS, MIN_FIT_POINTS, MP_TOL_REL,
    spectral_wavenumbers_3d, cylindrical_grid, to_cylindrical,
    vorticity_spectral, ddz_spectral, gamma_field, eta_field,
    omega_theta_over_r, swirl_source_eta, swirl_source_naive,
    source_near_axis_report, localise_rz, axisymmetry_defect, polar_resample,
    spectral_tail_fraction, burgers_vortex_fields, planted_swirl_fields,
    abc_flow, axisym_snapshot, gamma_max_principle_check, axisym_swirl_ic,
    calibrate_burgers, calibrate_planted_source,
)

N = 48


# ── grid and cylindrical algebra ─────────────────────────────────────────────

def test_cylindrical_grid_shapes_and_axis_mask():
    r, cp, sp, am = cylindrical_grid(N)
    assert r.shape == (N, N, N)
    assert am.dtype == bool
    # exactly one (i,j) column sits on the axis when the axis is a grid point
    assert np.count_nonzero(am[:, :, 0]) == 1
    assert np.max(r[~am]) <= np.sqrt(2) * np.pi + 1e-12


def test_cos_sin_phi_are_a_unit_vector_off_axis():
    r, cp, sp, am = cylindrical_grid(N)
    n = cp ** 2 + sp ** 2
    assert np.max(np.abs(n[~am] - 1.0)) < 1e-12


def test_to_cylindrical_preserves_norm():
    rng = np.random.default_rng(0)
    r, cp, sp, am = cylindrical_grid(16)
    u, v, w = (rng.standard_normal((16, 16, 16)) for _ in range(3))
    ur, ut, uz = to_cylindrical(u, v, w, cp, sp)
    lhs = u ** 2 + v ** 2 + w ** 2
    rhs = ur ** 2 + ut ** 2 + uz ** 2
    assert np.max(np.abs(lhs - rhs)[~am]) < 1e-10


def test_to_cylindrical_on_pure_swirl():
    r, cp, sp, am = cylindrical_grid(N)
    u, v = -sp * r, cp * r          # u_theta = r, u_r = 0
    ur, ut, _ = to_cylindrical(u, v, np.zeros_like(u), cp, sp)
    assert np.max(np.abs(ur[~am])) < 1e-12
    assert np.max(np.abs((ut - r)[~am])) < 1e-12


def test_ddz_spectral_on_a_known_mode():
    z = np.linspace(0, 2 * np.pi, N, endpoint=False)
    f = np.broadcast_to(np.sin(3 * z)[None, None, :], (N, N, N)).copy()
    _, _, kz = spectral_wavenumbers_3d(N)
    got = ddz_spectral(f, kz)
    want = np.broadcast_to(3 * np.cos(3 * z)[None, None, :], (N, N, N))
    assert np.max(np.abs(got - want)) < 1e-11


def test_ddz_of_z_independent_field_is_exactly_zero():
    r, _, _, _ = cylindrical_grid(N)
    _, _, kz = spectral_wavenumbers_3d(N)
    assert np.max(np.abs(ddz_spectral(r ** 2, kz))) < 1e-12


# ── the three WP4a quantities ────────────────────────────────────────────────

def test_gamma_is_r_times_u_theta():
    r, cp, sp, am = cylindrical_grid(N)
    u, v = -sp * r, cp * r
    assert np.max(np.abs((gamma_field(u, v, r, cp, sp) - r ** 2)[~am])) < 1e-12


def test_gamma_vanishes_for_purely_radial_flow():
    r, cp, sp, am = cylindrical_grid(N)
    u, v = cp * r, sp * r
    assert np.max(np.abs(gamma_field(u, v, r, cp, sp))) < 1e-12


def test_eta_is_gamma_over_r_squared():
    r, cp, sp, am = cylindrical_grid(N)
    u, v = -sp * r ** 2, cp * r ** 2
    G = gamma_field(u, v, r, cp, sp)
    e = eta_field(u, v, r, cp, sp, am)
    rs = np.where(am, 1.0, r)
    assert np.max(np.abs((e - G / rs ** 2)[~am])) < 1e-10


def test_the_two_source_forms_agree_identically_on_smooth_data():
    """(1/r^4) d_z(Gamma^2) == d_z(eta^2) is an exact identity because r does
    not depend on z.  Both branches are implemented; their agreement is the
    near-axis evaluability measurement, so the identity itself is tested."""
    (u, v, w), exact, g = planted_swirl_fields(N)
    _, _, kz = spectral_wavenumbers_3d(N)
    a, _ = swirl_source_eta(u, v, g["r"], g["cos_phi"], g["sin_phi"],
                            g["axis_mask"], kz)
    b = swirl_source_naive(u, v, g["r"], g["cos_phi"], g["sin_phi"],
                           g["axis_mask"], kz)
    sel = (~g["axis_mask"]) & (g["r"] <= 2.0)
    assert np.max(np.abs(a[sel] - b[sel])) < 1e-10 * np.max(np.abs(a[sel]))


def test_omega_theta_is_zero_for_a_purely_axial_vortex():
    r, cp, sp, am = cylindrical_grid(N)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    u, v = -sp * r * np.exp(-r ** 2), cp * r * np.exp(-r ** 2)
    w = np.zeros_like(u)
    wx, wy, wz = vorticity_spectral(u, v, w, kx, ky, kz)
    _, wth = omega_theta_over_r(wx, wy, r, cp, sp, am)
    assert np.max(np.abs(wth)) < 1e-10 * np.max(np.abs(wz))


# ── calibration against closed forms (the house rule) ───────────────────────

def test_burgers_calibration_passes():
    res = calibrate_burgers(N=N)
    assert res["verdict"] == "PASS"
    assert res["gamma_rel_err"] < 1e-12


def test_burgers_gamma_matches_closed_form():
    (u, v, w), exact, g = burgers_vortex_fields(N)
    G = gamma_field(u, v, g["r"], g["cos_phi"], g["sin_phi"])
    sel = ~g["axis_mask"]
    assert np.max(np.abs(G[sel] - exact["Gamma"][sel])) < 1e-12


def test_burgers_source_is_identically_zero():
    """Gamma_Burgers is z-independent, so the source vanishes identically.
    This is exactly why Burgers alone cannot calibrate the source routine --
    a routine returning 0 would also pass.  Hence the planted test below."""
    (u, v, w), exact, g = burgers_vortex_fields(N)
    _, _, kz = spectral_wavenumbers_3d(N)
    src, _ = swirl_source_eta(u, v, g["r"], g["cos_phi"], g["sin_phi"],
                              g["axis_mask"], kz)
    assert np.max(np.abs(src)) < 1e-14


def test_planted_source_calibration_is_nondegenerate_and_passes():
    res = calibrate_planted_source(N=N)
    assert res["verdict"] == "PASS"
    assert res["source_eta_rel_err"] < 1e-10
    (_, _, _), exact, _ = planted_swirl_fields(N)
    assert np.max(np.abs(exact["source"])) > 0.1     # not the trivial zero


def test_planted_source_matches_analytic_dz_eta_squared():
    (u, v, w), exact, g = planted_swirl_fields(N)
    _, _, kz = spectral_wavenumbers_3d(N)
    src, _ = swirl_source_eta(u, v, g["r"], g["cos_phi"], g["sin_phi"],
                              g["axis_mask"], kz)
    sel = (~g["axis_mask"]) & (g["r"] <= 2.0)
    assert np.max(np.abs(src[sel] - exact["source"][sel])) < 1e-10


# ── axisymmetry defect ───────────────────────────────────────────────────────

def test_defect_is_small_on_an_axisymmetric_field():
    ic = axisym_swirl_ic(N)
    g = ic["grid"]
    D = axisymmetry_defect(ic["u"], ic["v"], ic["w"], g["r"], g["cos_phi"],
                           g["sin_phi"], N)
    assert D < 5e-3


def test_defect_is_large_on_a_deliberately_non_axisymmetric_field():
    """The control the earlier (r,z)-binning version of this function failed:
    it reported D ~ 0.19 on an axisymmetric field, i.e. it could not tell the
    two apart."""
    ic = axisym_swirl_ic(N)
    g = ic["grid"]
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    u2 = ic["u"] * (1.0 + 0.5 * np.cos(x)[:, None, None])
    D_ax = axisymmetry_defect(ic["u"], ic["v"], ic["w"], g["r"], g["cos_phi"],
                              g["sin_phi"], N)
    D_na = axisymmetry_defect(u2, ic["v"], ic["w"], g["r"], g["cos_phi"],
                              g["sin_phi"], N)
    assert D_na > 20.0 * D_ax


def test_polar_resample_reproduces_a_known_azimuthal_mode():
    r, cp, sp, am = cylindrical_grid(N)
    f = cp                      # cos(phi): pure m = 1
    rr, tab = polar_resample(f, N, n_r=12, n_theta=32, r_max=2.0)
    fm = np.abs(np.fft.fft(tab, axis=1)) / 32
    assert np.max(fm[:, 1, :]) > 0.4
    assert np.max(fm[:, 2, :]) < 1e-3


# ── maximum principle check ──────────────────────────────────────────────────

def test_mp_passes_on_a_monotone_decreasing_series():
    res = gamma_max_principle_check([1.0, 0.9, 0.8, 0.7])
    assert res["verdict"] == "PASS"
    assert res["max_relative_rise"] == 0.0


def test_mp_fails_on_a_rising_series():
    res = gamma_max_principle_check([1.0, 0.9, 1.2])
    assert res["verdict"] == "FAIL"
    assert res["violated"]
    assert res["max_relative_rise"] == pytest.approx(0.3, rel=1e-9)


def test_mp_tolerance_is_widened_by_the_axisymmetry_defect():
    """A rise smaller than the flow's departure from the axisymmetric class is
    not evidence of a solver bug: Gamma has no maximum principle off that
    class."""
    assert gamma_max_principle_check([1.0, 1.05], [0.0, 0.2])["verdict"] == "PASS"
    assert gamma_max_principle_check([1.0, 1.05], [0.0, 0.0])["verdict"] == "FAIL"


def test_mp_reports_insufficient_on_a_single_point():
    assert gamma_max_principle_check([1.0])["verdict"] == "INSUFFICIENT"


def test_mp_uses_running_minimum_not_just_the_first_value():
    """A dip then a recovery above the dip is a violation even if the series
    never exceeds its initial value."""
    res = gamma_max_principle_check([1.0, 0.5, 0.9], [0.0, 0.0, 0.0])
    assert res["violated"]
    assert res["net_relative_rise"] <= 0.0


# ── axisymmetric IC ──────────────────────────────────────────────────────────

def test_axisym_swirl_ic_is_divergence_free_to_roundoff():
    ic = axisym_swirl_ic(N)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    d = lambda f, k: np.real(np.fft.ifftn(1j * k * np.fft.fftn(f)))
    div = d(ic["u"], kx) + d(ic["v"], ky) + d(ic["w"], kz)
    scale = float(np.sqrt(np.mean(ic["u"] ** 2 + ic["v"] ** 2 + ic["w"] ** 2)))
    assert float(np.sqrt(np.mean(div ** 2))) < 1e-3 * scale


def test_axisym_swirl_ic_has_a_nontrivial_source():
    """If Gamma were z-independent the source would vanish and WP4a would be
    localising nothing."""
    ic = axisym_swirl_ic(N)
    g = ic["grid"]
    _, _, kz = spectral_wavenumbers_3d(N)
    src, _ = swirl_source_eta(ic["u"], ic["v"], g["r"], g["cos_phi"],
                              g["sin_phi"], g["axis_mask"], kz)
    assert np.max(np.abs(src)) > 0.1


def test_axisym_swirl_ic_has_nonzero_swirl():
    ic = axisym_swirl_ic(N)
    g = ic["grid"]
    G = gamma_field(ic["u"], ic["v"], g["r"], g["cos_phi"], g["sin_phi"])
    assert np.max(np.abs(G)) > 0.1


# ── snapshot and utilities ───────────────────────────────────────────────────

def test_snapshot_has_all_wp4a_quantities_and_is_finite():
    ic = axisym_swirl_ic(N)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    s = axisym_snapshot(ic["u"], ic["v"], ic["w"], ic["grid"], kx, ky, kz, N)
    for key in ("Gamma_absmax", "omega_theta_over_r_absmax", "source_absmax",
                "source_naive_absmax", "axisym_defect", "tail_fraction",
                "r_at_Gamma_max", "Gamma_absmax_core", "M"):
        assert key in s and np.isfinite(s[key])


def test_tail_fraction_is_tiny_for_a_smooth_field_and_large_for_noise():
    ic = axisym_swirl_ic(N)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    smooth = spectral_tail_fraction(ic["u"], ic["v"], ic["w"], kx, ky, kz, N)
    rng = np.random.default_rng(1)
    noise = [rng.standard_normal((N, N, N)) for _ in range(3)]
    rough = spectral_tail_fraction(*noise, kx, ky, kz, N)
    assert smooth < 1e-6
    assert rough > 0.3


def test_localise_rz_returns_the_requested_table_shape():
    ic = axisym_swirl_ic(N)
    rc, zc, tab = localise_rz(ic["w"], ic["grid"]["r"], N, n_r=10, n_z=8)
    assert tab.shape == (10, 8)
    assert len(rc) == 10 and len(zc) == 8
    assert np.all(np.isfinite(tab))


def test_near_axis_report_is_zero_discrepancy_on_analytic_data():
    (u, v, w), exact, g = planted_swirl_fields(N)
    _, _, kz = spectral_wavenumbers_3d(N)
    a, _ = swirl_source_eta(u, v, g["r"], g["cos_phi"], g["sin_phi"],
                            g["axis_mask"], kz)
    b = swirl_source_naive(u, v, g["r"], g["cos_phi"], g["sin_phi"],
                           g["axis_mask"], kz)
    rep = source_near_axis_report(a, b, g["r"], g["axis_mask"])
    assert rep["rel_discrepancy_max"] < 1e-8
    assert rep["r_loss_radius"] == 0.0


def test_near_axis_report_detects_a_corrupted_naive_branch():
    """The report must be able to SEE a near-axis failure, or its passing is
    worthless."""
    (u, v, w), exact, g = planted_swirl_fields(N)
    _, _, kz = spectral_wavenumbers_3d(N)
    a, _ = swirl_source_eta(u, v, g["r"], g["cos_phi"], g["sin_phi"],
                            g["axis_mask"], kz)
    b = a.copy()
    inner = (g["r"] < 0.6) & (~g["axis_mask"])
    b[inner] = b[inner] + 50.0 * np.max(np.abs(a))
    rep = source_near_axis_report(a, b, g["r"], g["axis_mask"])
    assert rep["rel_discrepancy_innermost"] > 1.0
    assert rep["r_loss_radius"] > 0.0


def test_abc_flow_is_beltrami():
    """omega = u for the ABC flow -- the property that makes it an exact
    periodic NS solution and hence the solver's calibration target."""
    u, v, w = abc_flow(N)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    wx, wy, wz = vorticity_spectral(u, v, w, kx, ky, kz)
    scale = float(np.max(np.abs(u)))
    assert np.max(np.abs(wx - u)) < 1e-10 * scale
    assert np.max(np.abs(wy - v)) < 1e-10 * scale
    assert np.max(np.abs(wz - w)) < 1e-10 * scale


def test_min_fit_points_constant_is_at_least_three():
    """S49: results.db carries a PASS on a 2-point regression with r^2 = 1.
    This module refuses any fit below MIN_FIT_POINTS."""
    assert MIN_FIT_POINTS >= 3


def test_constants_are_documented_values():
    assert AXIS_EPS_CELLS == 0.5
    assert MP_TOL_REL > 0.0
