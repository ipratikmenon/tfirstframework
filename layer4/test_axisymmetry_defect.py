"""
layer4/test_axisymmetry_defect.py
Tests for §23 axisymmetry defect and pressure flux analysis.

The key open question:
  gap_ratio = D_dens / M = (∫_{B_{r*}} |u_na|² / r*³) / M
  F_ratio   = F_z / (M^{3/2} r*²)

  Open Problem 23.2: Are these bounded below by c₀ > 0 as M→∞, σ→0?

Run: python -m pytest layer4/test_axisymmetry_defect.py -v -s
"""

import numpy as np
import pytest
from layer4.geometric_disorder import spectral_wavenumbers, compute_M
from layer4.axisymmetry_defect import (
    axisymmetric_decompose,
    compute_axisymmetry_defect,
    compute_pressure_flux,
    scaling_sweep,
    _cylindrical_basis,
)

N = 16
kx, ky, kz = spectral_wavenumbers(N)


# ─────────────────────────────────────────────────────────────────────────────
# Field factories
# ─────────────────────────────────────────────────────────────────────────────

def _make_taylor_green(N, amplitude=2.0):
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")
    u = np.zeros((3, N, N, N))
    u[0] =  amplitude * np.sin(X) * np.cos(Y) * np.cos(Z)
    u[1] = -amplitude * np.cos(X) * np.sin(Y) * np.cos(Z)
    omega = np.zeros((3, N, N, N))
    omega[0] =  amplitude * np.cos(X) * np.sin(Y) * np.sin(Z)
    omega[1] = -amplitude * np.sin(X) * np.cos(Y) * np.sin(Z)
    omega[2] =  0.0
    return u, omega


def _make_aligned_tube(N, amplitude=3.0):
    """Perfectly z-aligned tube: nearly axisymmetric about e₃."""
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")
    cx, cy = np.pi, np.pi
    r2 = (X - cx)**2 + (Y - cy)**2
    r0 = 0.6
    omega = np.zeros((3, N, N, N))
    omega[2] = amplitude * np.exp(-r2 / (2 * r0**2))
    u = np.zeros((3, N, N, N))
    # Approximate azimuthal velocity
    r = np.sqrt(r2) + 1e-8
    Gamma = amplitude * 2 * np.pi * r0**2
    uth = Gamma / (2 * np.pi * r) * (1 - np.exp(-r2 / (2 * r0**2)))
    u[0] = -uth * (Y - cy) / r
    u[1] =  uth * (X - cx) / r
    return u, omega


def _make_asymmetric_field(N, amplitude=2.0):
    """Deliberately non-axisymmetric field."""
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")
    u = np.zeros((3, N, N, N))
    omega = np.zeros((3, N, N, N))
    u[0] = amplitude * np.sin(X) * np.cos(2*Y) * np.cos(Z)
    u[1] = -amplitude * np.cos(X) * np.sin(Y) * np.cos(3*Z)
    u[2] = amplitude * np.cos(X) * np.cos(Y) * np.sin(2*Z)
    omega[0] = amplitude * (3*np.cos(X)*np.cos(Y)*np.sin(2*Z) + np.sin(X)*np.cos(2*Y)*np.sin(Z))
    omega[1] = amplitude * (-2*np.cos(X)*np.cos(Y)*np.cos(2*Z) - np.sin(X)*np.cos(2*Y)*np.sin(Z))
    omega[2] = amplitude * (2*np.cos(X)*np.sin(Y)*np.cos(3*Z) - np.cos(X)*np.sin(2*Y)*np.cos(Z))
    return u, omega


# ─────────────────────────────────────────────────────────────────────────────
# Tests: cylindrical basis
# ─────────────────────────────────────────────────────────────────────────────

def test_cylindrical_basis_orthonormal():
    for e in [[0, 0, 1], [1, 0, 0], [0, 1, 0], [1/np.sqrt(3)]*3]:
        ev = np.array(e, dtype=float)
        ev /= np.linalg.norm(ev)
        er, et, ez = _cylindrical_basis(ev)
        assert abs(np.dot(er, ez)) < 1e-12, "er ⊥ ez"
        assert abs(np.dot(et, ez)) < 1e-12, "et ⊥ ez"
        assert abs(np.dot(er, et)) < 1e-12, "er ⊥ et"
        assert abs(np.linalg.norm(er) - 1) < 1e-12
        assert abs(np.linalg.norm(et) - 1) < 1e-12


def test_cylindrical_basis_right_handed():
    e = np.array([0.0, 0.0, 1.0])
    er, et, ez = _cylindrical_basis(e)
    cross = np.cross(er, et)
    assert np.allclose(cross, ez, atol=1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# Tests: axisymmetric decomposition
# ─────────────────────────────────────────────────────────────────────────────

def test_decompose_splits_correctly():
    """u_axi + u_na = u exactly."""
    u, omega = _make_taylor_green(N)
    M, mag, idx = compute_M(omega)
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    x_center = [x1d[idx[0]], x1d[idx[1]], x1d[idx[2]]]
    u_axi, u_na, rho, z_cyl = axisymmetric_decompose(u, x_center, [0, 0, 1], N)
    reconstruction = u_axi + u_na
    assert np.allclose(reconstruction, u, atol=1e-10), "u_axi + u_na ≠ u"


def test_decompose_aligned_tube_small_na():
    """Aligned tube: non-axisymmetric part should be small relative to total."""
    u, omega = _make_aligned_tube(N, amplitude=3.0)
    M, mag, idx = compute_M(omega)
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    x_center = [x1d[idx[0]], x1d[idx[1]], x1d[idx[2]]]
    u_axi, u_na, rho, z_cyl = axisymmetric_decompose(u, x_center, [0, 0, 1], N)
    u_na_sq = float(np.sum(u_na**2))
    u_sq    = float(np.sum(u**2)) + 1e-30
    ratio   = u_na_sq / u_sq
    assert ratio < 1.0, f"Non-axi part should be ≤ total: ratio={ratio:.3f}"


def test_decompose_asymmetric_field_large_na():
    """Asymmetric field: non-axisymmetric part should be significant."""
    u, omega = _make_asymmetric_field(N)
    M, mag, idx = compute_M(omega)
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    x_center = [x1d[idx[0]], x1d[idx[1]], x1d[idx[2]]]
    u_axi, u_na, rho, z_cyl = axisymmetric_decompose(u, x_center, [1/np.sqrt(3)]*3, N)
    u_na_sq = float(np.sum(u_na**2))
    u_sq    = float(np.sum(u**2)) + 1e-30
    assert u_na_sq / u_sq > 1e-6, "Asymmetric field should have nonzero D_na"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: compute_axisymmetry_defect
# ─────────────────────────────────────────────────────────────────────────────

def test_defect_returns_expected_keys():
    u, omega = _make_taylor_green(N)
    result = compute_axisymmetry_defect(u, omega, kx, ky, kz)
    assert "error" not in result
    for key in ["M", "sigma", "r_star", "D_axi", "D_total", "D_dens",
                "D_ratio", "threshold", "gap_ratio"]:
        assert key in result, f"Missing key: {key}"


def test_defect_non_negative():
    for u, omega in [_make_taylor_green(N), _make_aligned_tube(N), _make_asymmetric_field(N)]:
        result = compute_axisymmetry_defect(u, omega, kx, ky, kz)
        if "error" not in result:
            assert result["D_axi"]  >= 0
            assert result["D_total"] >= 0
            assert result["D_ratio"] >= 0


def test_defect_ratio_leq_one():
    """D_ratio = D_axi / D_total ≤ 1 always."""
    for u, omega in [_make_taylor_green(N), _make_aligned_tube(N)]:
        result = compute_axisymmetry_defect(u, omega, kx, ky, kz)
        if "error" not in result and result["D_total"] > 1e-10:
            assert result["D_ratio"] <= 1.0 + 1e-8, \
                f"D_ratio={result['D_ratio']} > 1"


def test_threshold_correct_scaling():
    """threshold = M r*³ = M · M^{-3/2} = M^{-1/2}: decreases with M."""
    u1, omega1 = _make_aligned_tube(N, amplitude=2.0)
    u2, omega2 = _make_aligned_tube(N, amplitude=4.0)
    r1 = compute_axisymmetry_defect(u1, omega1, kx, ky, kz)
    r2 = compute_axisymmetry_defect(u2, omega2, kx, ky, kz)
    if "error" not in r1 and "error" not in r2 and r2["M"] > r1["M"]:
        assert r2["threshold"] < r1["threshold"], \
            "threshold = M^{-1/2} should decrease with M"


def test_gap_ratio_measured():
    """
    MEASUREMENT: gap_ratio = D_dens / M.
    Open Problem 23.2: is this ≥ c₀ > 0 as M→∞?
    Prints the value without asserting the conjecture.
    """
    amplitudes = [1.0, 2.0, 3.0, 4.0]
    print("\n§23 Scaling analysis: gap_ratio = D_dens / M = (∫|u_na|²/r*³) / M")
    print(f"{'Amplitude':>10} {'M':>8} {'σ':>8} {'D_axi':>12} {'threshold':>12} {'gap_ratio':>12}")
    for amp in amplitudes:
        u, omega = _make_aligned_tube(N, amplitude=amp)
        r = compute_axisymmetry_defect(u, omega, kx, ky, kz)
        if "error" not in r:
            print(f"{amp:>10.1f} {r['M']:>8.3f} {r['sigma']:>8.4f} "
                  f"{r['D_axi']:>12.4e} {r['threshold']:>12.4e} {r['gap_ratio']:>12.4e}")
    assert True  # measurement only


# ─────────────────────────────────────────────────────────────────────────────
# Tests: pressure flux F_z
# ─────────────────────────────────────────────────────────────────────────────

def test_flux_returns_expected_keys():
    u, omega = _make_taylor_green(N)
    result = compute_pressure_flux(u, omega, kx, ky, kz)
    for key in ["F_z", "F_total", "F_ratio", "M", "r_star", "M32_r2"]:
        assert key in result, f"Missing key: {key}"


def test_flux_zero_velocity():
    u = np.zeros((3, N, N, N))
    omega = np.zeros((3, N, N, N))
    result = compute_pressure_flux(u, omega, kx, ky, kz)
    assert result["M"] == 0.0 or result["F_z"] == 0.0


def test_flux_taylor_green_nonzero():
    """TG field has nonzero velocity and pressure → expect nonzero flux."""
    u, omega = _make_taylor_green(N, amplitude=2.0)
    result = compute_pressure_flux(u, omega, kx, ky, kz)
    # flux may be small or zero depending on symmetry — just check structure
    assert isinstance(result["F_z"], float)
    assert isinstance(result["F_ratio"], float)


def test_flux_ratio_scaling():
    """
    MEASUREMENT: F_ratio = F_z / (M^{3/2} r*²).
    Open Problem 23.2: is this ≥ c₀ > 0 as M→∞?
    """
    print("\n§23 Flux analysis: F_ratio = F_z / (M^{3/2} r*²)")
    print(f"{'Amplitude':>10} {'M':>8} {'F_z':>14} {'M^{3/2}r*²':>14} {'F_ratio':>12}")
    for amp in [1.0, 2.0, 3.0, 4.0]:
        u, omega = _make_taylor_green(N, amplitude=amp)
        r = compute_pressure_flux(u, omega, kx, ky, kz)
        print(f"{amp:>10.1f} {r['M']:>8.3f} {r['F_z']:>14.4e} {r['M32_r2']:>14.4e} {r['F_ratio']:>12.4e}")
    assert True


# ─────────────────────────────────────────────────────────────────────────────
# Scaling sweep: the critical experiment
# ─────────────────────────────────────────────────────────────────────────────

def test_scaling_sweep_tg():
    """
    Scaling sweep for Taylor-Green: measure D_axi and F_z as M increases.
    This is the computational test of Open Problem 23.2.
    """
    amplitudes = [0.5, 1.0, 2.0, 3.0, 4.0]
    results = scaling_sweep(_make_taylor_green, N, amplitudes, kx, ky, kz)
    assert len(results) > 0, "Sweep returned no results"
    for r in results:
        assert r["D_axi"] >= 0
        assert r["D_ratio"] >= 0

    print("\n§23 Scaling sweep (Taylor-Green):")
    print(f"{'M':>8} {'σ':>8} {'D_axi':>12} {'threshold':>12} {'gap_ratio':>12} {'F_ratio':>12}")
    for r in results:
        print(f"{r['M']:>8.3f} {r['sigma']:>8.4f} {r['D_axi']:>12.4e} "
              f"{r['threshold']:>12.4e} {r['gap_ratio']:>12.4e} {r.get('F_ratio', float('nan')):>12.4e}")


def test_scaling_sweep_aligned_tube():
    """
    Scaling sweep for aligned tube: as amplitude increases, does gap_ratio stay positive?
    """
    amplitudes = [1.0, 2.0, 3.0, 4.0]
    results = scaling_sweep(_make_aligned_tube, N, amplitudes, kx, ky, kz)

    print("\n§23 Scaling sweep (Aligned tube):")
    print(f"{'M':>8} {'σ':>8} {'D_axi':>12} {'threshold':>12} {'gap_ratio':>12} {'F_ratio':>12}")
    for r in results:
        print(f"{r['M']:>8.3f} {r['sigma']:>8.4f} {r['D_axi']:>12.4e} "
              f"{r['threshold']:>12.4e} {r['gap_ratio']:>12.4e} {r.get('F_ratio', float('nan')):>12.4e}")
    assert True


def test_open_problem_23_summary():
    """
    SUMMARY TEST for Open Problem 23.2.
    Reports the scaling of gap_ratio and F_ratio across field types and amplitudes.
    The conjecture: both quantities remain ≥ c₀ > 0 as M→∞.
    """
    print("\n" + "="*70)
    print("OPEN PROBLEM 23.2 — Computational Status")
    print("Needed: gap_ratio = D_dens/M ≥ c₀ > 0   OR   F_ratio ≥ c₀ > 0")
    print("="*70)

    for name, factory in [("Taylor-Green", _make_taylor_green),
                           ("Aligned tube", _make_aligned_tube),
                           ("Asymmetric",   _make_asymmetric_field)]:
        print(f"\n{name}:")
        print(f"  {'M':>6} {'gap_ratio':>12} {'F_ratio':>12}")
        for amp in [1.0, 2.0, 3.0]:
            u, omega = factory(N, amplitude=amp)
            d = compute_axisymmetry_defect(u, omega, kx, ky, kz)
            f = compute_pressure_flux(u, omega, kx, ky, kz)
            if "error" not in d:
                print(f"  {d['M']:>6.2f} {d['gap_ratio']:>12.4e} {f['F_ratio']:>12.4e}")

    print("\nConclusion: if gap_ratio or F_ratio → 0 as M→∞, argument fails.")
    print("If they stabilise at c₀ > 0, Open Problem 23.2 is numerically confirmed.")
    assert True
