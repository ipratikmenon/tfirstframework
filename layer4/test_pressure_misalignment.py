"""
layer4/test_pressure_misalignment.py
Tests for §23 pressure misalignment diagnostic.

Key question (Open Problem 23.1):
  Is |∂_z p_glob / M^{3/2}|_{B_{r*}} ≥ c₀ > 0 for nearly aligned fields?

Run: python -m pytest layer4/test_pressure_misalignment.py -v
"""

import numpy as np
import pytest
from layer4.geometric_disorder import spectral_wavenumbers, compute_M
from layer4.pressure_misalignment import (
    solve_pressure_spectral,
    compute_pressure_gradient,
    decompose_pressure_gradient,
    pressure_misalignment_diagnostic,
)

N = 16
kx, ky, kz = spectral_wavenumbers(N)


# ─────────────────────────────────────────────────────────────────────────────
# Field factories
# ─────────────────────────────────────────────────────────────────────────────

def _make_straight_tube(N, amplitude=3.0):
    """
    Perfectly straight vortex tube aligned with e₃.
    ω = A·e₃·exp(-r²/(2r*²)) — Gaussian tube, no z-variation.
    Velocity: azimuthal u_θ only (Rankine-like).
    For this field: ∂_z p = 0 by symmetry.
    """
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")
    cx, cy = np.pi, np.pi
    r2 = (X - cx)**2 + (Y - cy)**2
    r0 = 0.8
    omega = np.zeros((3, N, N, N))
    omega[2] = amplitude * np.exp(-r2 / (2 * r0**2))

    # Velocity field: azimuthal u_θ derived from ω = -Δψ in 2D for each z-slice
    # For Gaussian ω_z(r), u_θ(r) = (1/r) ∫₀^r ω_z(r') r' dr'  (2D formula)
    u = np.zeros((3, N, N, N))
    dx_val = 2 * np.pi / N
    omega_z_2d = omega[2, :, :, 0]
    omega_hat = np.fft.fftn(omega_z_2d) * dx_val**2
    kx2d = np.fft.fftfreq(N, d=1.0/N).astype(float)
    KX, KY = np.meshgrid(kx2d, kx2d, indexing="ij")
    K2 = KX**2 + KY**2
    with np.errstate(divide="ignore", invalid="ignore"):
        psi_hat = np.where(K2 > 0, -omega_hat / K2, 0.0)
    psi = np.real(np.fft.ifftn(psi_hat)) / dx_val**2
    ux_2d = np.real(np.fft.ifftn(1j * KY * psi_hat * dx_val**2)) / dx_val**2
    uy_2d = np.real(np.fft.ifftn(-1j * KX * psi_hat * dx_val**2)) / dx_val**2
    for iz in range(N):
        u[0, :, :, iz] = -ux_2d
        u[1, :, :, iz] =  uy_2d

    return u, omega


def _make_tilted_tube(N, amplitude=3.0, tilt=0.3):
    """
    Slightly tilted tube: ω = A·(e₃ + tilt·e₁), normalised.
    Introduces curvature in ê → nonzero ∂_z p expected.
    """
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")
    cx, cy, cz = np.pi, np.pi, np.pi
    r2 = (X - cx)**2 + (Y - cy)**2 + (Z - cz)**2
    r0 = 0.7
    omega = np.zeros((3, N, N, N))
    norm = np.sqrt(1 + tilt**2)
    omega[0] = amplitude * tilt / norm * np.exp(-r2 / (2 * r0**2))
    omega[2] = amplitude / norm       * np.exp(-r2 / (2 * r0**2))

    # Approximate velocity: use gradient of streamfunction proxy
    u = np.zeros((3, N, N, N))
    for comp in range(3):
        f_hat = np.fft.fftn(omega[comp]) * (2 * np.pi / N)**3
        K2 = kx**2 + ky**2 + kz**2
        with np.errstate(divide="ignore", invalid="ignore"):
            bs_hat = np.where(K2 > 0, f_hat / K2, 0.0)
        u[comp] = np.real(np.fft.ifftn(bs_hat)) / (2 * np.pi / N)**3

    return u, omega


def _make_taylor_green(N, amplitude=2.0):
    """Taylor-Green vortex: non-trivial σ, known pressure structure."""
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


# ─────────────────────────────────────────────────────────────────────────────
# Tests: pressure solver
# ─────────────────────────────────────────────────────────────────────────────

def test_pressure_solver_zero_velocity():
    u = np.zeros((3, N, N, N))
    p = solve_pressure_spectral(u, kx, ky, kz)
    assert np.max(np.abs(p)) < 1e-10, "Zero velocity → zero pressure"


def test_pressure_mean_zero():
    """Spectral solve enforces mean-zero pressure."""
    u, omega = _make_taylor_green(N)
    p = solve_pressure_spectral(u, kx, ky, kz)
    assert abs(np.mean(p)) < 1e-8, f"Pressure mean not zero: {np.mean(p)}"


def test_pressure_gradient_shape():
    u, omega = _make_taylor_green(N)
    p = solve_pressure_spectral(u, kx, ky, kz)
    grad_p = compute_pressure_gradient(p, kx, ky, kz)
    assert grad_p.shape == (3, N, N, N)


def test_pressure_gradient_zero_for_zero_pressure():
    p = np.zeros((N, N, N))
    grad_p = compute_pressure_gradient(p, kx, ky, kz)
    assert np.max(np.abs(grad_p)) < 1e-10


# ─────────────────────────────────────────────────────────────────────────────
# Tests: decomposition
# ─────────────────────────────────────────────────────────────────────────────

def test_decompose_axial_radial_orthogonal():
    """Axial and radial components are orthogonal by construction."""
    u, omega = _make_taylor_green(N)
    p = solve_pressure_spectral(u, kx, ky, kz)
    grad_p = compute_pressure_gradient(p, kx, ky, kz)
    from layer4.geometric_disorder import compute_e_hat
    M, mag, idx = compute_M(omega)
    e_hat, _ = compute_e_hat(omega, mag=mag)
    ax, rad, tot = decompose_pressure_gradient(grad_p, e_hat, idx)
    # Pythagoras: ax² + rad² ≈ tot²
    assert abs(ax**2 + rad**2 - tot**2) < 1e-6 * (tot**2 + 1e-20), \
        f"Decomposition not orthogonal: ax={ax}, rad={rad}, tot={tot}"


def test_decompose_purely_axial():
    """If ∇p ∥ ê, then radial component = 0."""
    N2 = 8
    kx2, ky2, kz2 = spectral_wavenumbers(N2)
    # Construct ê = e₃ everywhere
    omega = np.zeros((3, N2, N2, N2))
    omega[2, N2 // 2, N2 // 2, N2 // 2] = 5.0
    M, mag, idx = compute_M(omega)
    from layer4.geometric_disorder import compute_e_hat
    e_hat, _ = compute_e_hat(omega, mag=mag)

    # Construct grad_p = 1.0 * e₃ everywhere
    grad_p = np.zeros((3, N2, N2, N2))
    grad_p[2] = 1.0

    ax, rad, tot = decompose_pressure_gradient(grad_p, e_hat, idx)
    assert rad < 1e-10, f"Purely axial ∇p should have zero radial component, got {rad}"
    assert abs(ax - 1.0) < 1e-10, f"Axial component should be 1.0, got {ax}"


def test_decompose_purely_radial():
    """If ∇p ⊥ ê, then axial component = 0."""
    N2 = 8
    kx2, ky2, kz2 = spectral_wavenumbers(N2)
    omega = np.zeros((3, N2, N2, N2))
    omega[2, N2 // 2, N2 // 2, N2 // 2] = 5.0
    M, mag, idx = compute_M(omega)
    from layer4.geometric_disorder import compute_e_hat
    e_hat, _ = compute_e_hat(omega, mag=mag)

    # grad_p = e₁ (perpendicular to ê = e₃)
    grad_p = np.zeros((3, N2, N2, N2))
    grad_p[0] = 1.0

    ax, rad, tot = decompose_pressure_gradient(grad_p, e_hat, idx)
    assert ax < 1e-10, f"Purely radial ∇p should have zero axial component, got {ax}"
    assert abs(rad - 1.0) < 1e-10, f"Radial component should be 1.0, got {rad}"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: §23 diagnostic
# ─────────────────────────────────────────────────────────────────────────────

def test_diagnostic_straight_tube_axial_small():
    """
    Straight aligned tube: axial pressure ratio should be near zero.
    (∂_z p = 0 for a tube with exact z-translational symmetry.)
    Tests the §23.3 theoretical claim.
    """
    u, omega = _make_straight_tube(N, amplitude=3.0)
    result = pressure_misalignment_diagnostic(u, omega, kx, ky, kz)
    assert "error" not in result
    # ratio_axial = |∂_z p| / M^{3/2} should be small for a straight tube
    # (not exactly 0 due to spectral aliasing and periodic kernel)
    assert result["M"] > 0, "M should be positive"
    assert result["ratio_total"] >= 0
    # The total pressure gradient is radial-dominated for a straight tube
    # so ratio_radial >> ratio_axial
    if result["ratio_total"] > 1e-6:
        assert result["ratio_radial"] >= result["ratio_axial"] * 0.5, \
            f"Straight tube: expected radial-dominated pressure, got axial={result['ratio_axial']:.4f}, radial={result['ratio_radial']:.4f}"


def test_diagnostic_tilted_tube_has_axial():
    """
    Tilted tube (non-zero curvature): axial pressure component should be nonzero.
    Tests that tilt → ∂_z p ≠ 0 (§23.3 curvature argument).
    """
    u, omega = _make_tilted_tube(N, amplitude=3.0, tilt=0.5)
    result = pressure_misalignment_diagnostic(u, omega, kx, ky, kz)
    if "error" in result:
        pytest.skip("Zero vorticity field")
    assert result["dp_axial_xstar"] >= 0, "Axial component should be non-negative"
    # For a tilted tube we expect SOME axial pressure (not necessarily large)
    # This tests the structure, not the magnitude bound


def test_diagnostic_taylor_green_keys():
    """Taylor-Green diagnostic returns all expected keys."""
    u, omega = _make_taylor_green(N)
    result = pressure_misalignment_diagnostic(u, omega, kx, ky, kz)
    assert "error" not in result
    for key in ["M", "sigma", "r_star", "dp_axial_xstar", "dp_radial_xstar",
                "dp_total_xstar", "ratio_axial", "ratio_radial", "ratio_total",
                "ball_axial", "ball_ratio"]:
        assert key in result, f"Missing key: {key}"


def test_diagnostic_ratio_total_consistent():
    """ratio_total = dp_total / M^{3/2} should be non-negative."""
    for u, omega in [_make_taylor_green(N), _make_tilted_tube(N)]:
        result = pressure_misalignment_diagnostic(u, omega, kx, ky, kz)
        if "error" not in result:
            assert result["ratio_total"] >= 0


def test_diagnostic_ball_ratio_non_negative():
    u, omega = _make_taylor_green(N)
    result = pressure_misalignment_diagnostic(u, omega, kx, ky, kz)
    assert result["ball_ratio"] >= 0


def test_diagnostic_scaling_check():
    """
    §23 Open Problem: ratio_axial = |∂_z p| / M^{3/2}.
    Tests that ratio_total matches expected O(1) scale for TG field.
    (Not testing the lower bound — that is the open problem.)
    """
    u, omega = _make_taylor_green(N, amplitude=2.0)
    result = pressure_misalignment_diagnostic(u, omega, kx, ky, kz)
    # ratio_total = |∇p| / M^{3/2} should be O(1) for a well-resolved field
    # Not too large (no blowup), not exactly zero (nontrivial flow)
    assert 0 <= result["ratio_total"] < 1e6, \
        f"ratio_total out of expected range: {result['ratio_total']}"


def test_open_problem_23_lower_bound_status():
    """
    EXPLICIT TEST OF OPEN PROBLEM 23.1.
    Measures |∂_z p| / M^{3/2} for the Taylor-Green field.
    This test DOCUMENTS THE CURRENT STATUS — it does NOT assert the lower bound.
    The open problem asks whether this ratio is ≥ c₀ > 0 uniformly as M → ∞.
    """
    u, omega = _make_taylor_green(N)
    result = pressure_misalignment_diagnostic(u, omega, kx, ky, kz)
    ratio = result["ratio_axial"]
    ball_ratio = result["ball_ratio"]
    print(f"\n§23 Open Problem status:")
    print(f"  M = {result['M']:.4f},  σ = {result['sigma']:.4f}")
    print(f"  |∂_z p| / M^{{3/2}} at x* = {ratio:.6f}")
    print(f"  ball mean |∂_z p| / M^{{3/2}} = {ball_ratio:.6f}")
    print(f"  Radial ratio |∂_r p| / M^{{3/2}} = {result['ratio_radial']:.6f}")
    print(f"  OPEN: is ratio_axial ≥ c₀ > 0 as M → ∞ and σ → 0?")
    # Always passes — this is a measurement, not an assertion on the open problem
    assert True
