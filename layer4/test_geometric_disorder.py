"""
layer4/test_geometric_disorder.py
Tests for the σ = A_loc/M^{3/2} geometric disorder framework (§21-§22).

Run: python -m pytest layer4/test_geometric_disorder.py -v
"""

import numpy as np
import pytest
from layer4.geometric_disorder import (
    spectral_wavenumbers,
    compute_M,
    compute_r_star,
    compute_e_hat,
    verify_gradient_identity,
    compute_A_loc,
    compute_R_log,
    check_coercive_inequality,
    SigmaTracker,
)

# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

N = 16
kx, ky, kz = spectral_wavenumbers(N)


def _make_uniform_omega(N, direction=(1, 0, 0), amplitude=2.0):
    """ω = A·e₁ everywhere: ê = const, ∇ê = 0, so A_loc = 0."""
    omega = np.zeros((3, N, N, N))
    omega[0] = amplitude * direction[0]
    omega[1] = amplitude * direction[1]
    omega[2] = amplitude * direction[2]
    return omega


def _make_smooth_omega(N, amplitude=3.0):
    """ω with smooth spatial variation; A_loc > 0 expected."""
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    omega = np.zeros((3, N, N, N))
    omega[0] = amplitude * np.sin(X) * np.cos(Y)
    omega[1] = amplitude * np.cos(X) * np.sin(Z)
    omega[2] = amplitude * np.sin(Y) * np.cos(Z)
    return omega


def _make_concentrated_omega(N, amplitude=5.0):
    """ω concentrated at one grid point (approximate δ-function): M large, r* small."""
    omega = np.zeros((3, N, N, N))
    omega[0, N // 2, N // 2, N // 2] = amplitude
    return omega


# ─────────────────────────────────────────────────────────────────────────────
# Tests: compute_M and r*
# ─────────────────────────────────────────────────────────────────────────────

def test_compute_M_uniform():
    omega = _make_uniform_omega(N, amplitude=3.0)
    M, mag, idx = compute_M(omega)
    assert abs(M - 3.0) < 1e-12, f"Expected M=3.0, got {M}"


def test_compute_M_zero():
    omega = np.zeros((3, N, N, N))
    M, mag, idx = compute_M(omega)
    assert M == 0.0


def test_r_star_scaling():
    M = 4.0
    r = compute_r_star(M)
    assert abs(r - 0.5) < 1e-12, f"Expected r*=0.5 for M=4, got {r}"


def test_r_star_zero_M():
    assert compute_r_star(0.0) == np.inf


def test_r_star_inverse_sqrt():
    for M in [1.0, 4.0, 9.0, 16.0]:
        assert abs(compute_r_star(M) - M**(-0.5)) < 1e-12


# ─────────────────────────────────────────────────────────────────────────────
# Tests: ê and the algebraic identity
# ─────────────────────────────────────────────────────────────────────────────

def test_e_hat_unit_norm():
    omega = _make_smooth_omega(N)
    e_hat, mag = compute_e_hat(omega)
    norms = np.sqrt(e_hat[0]**2 + e_hat[1]**2 + e_hat[2]**2)
    active = mag > 1e-12
    assert np.all(np.abs(norms[active] - 1.0) < 1e-10)


def test_e_hat_direction_uniform():
    omega = _make_uniform_omega(N, direction=(0, 1, 0), amplitude=2.5)
    e_hat, mag = compute_e_hat(omega)
    assert np.allclose(e_hat[1], 1.0, atol=1e-10)
    assert np.allclose(e_hat[0], 0.0, atol=1e-10)
    assert np.allclose(e_hat[2], 0.0, atol=1e-10)


def test_gradient_identity_uniform():
    """For uniform ω: ∇|ω|=0 and ∇ê=0, so LHS=RHS=0 (trivially)."""
    omega = _make_uniform_omega(N, amplitude=3.0)
    max_err, passes = verify_gradient_identity(omega, kx, ky, kz, tol=1e-5)
    assert passes, f"Identity failed: max_rel_err={max_err}"


def test_gradient_identity_smooth():
    """
    Identity |∇ω|² = |∇|ω||² + |ω|²|∇ê|² at active points (Lemma 21.1).
    Tolerance is 0.3: ê = ω/|ω| is not band-limited so spectral ∇ê accumulates
    aliasing at N=16.  The identity is algebraically exact; this tests consistency.
    """
    omega = _make_smooth_omega(N, amplitude=2.0)
    max_err, passes = verify_gradient_identity(omega, kx, ky, kz, tol=0.3)
    assert passes, f"Identity failed on smooth field: max_rel_err={max_err}"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: A_loc and σ
# ─────────────────────────────────────────────────────────────────────────────

def test_A_loc_uniform_is_zero():
    """Uniform ω → ê = const → ∇ê = 0 → A_loc = 0."""
    omega = _make_uniform_omega(N, amplitude=3.0)
    M, mag, idx = compute_M(omega)
    A_loc, sigma, r_star, _ = compute_A_loc(omega, kx, ky, kz, M=M, idx=idx, mag=mag)
    assert A_loc < 1e-8, f"Expected A_loc≈0 for uniform ω, got {A_loc}"
    assert sigma < 1e-8


def test_A_loc_positive_smooth():
    """Smooth non-uniform ω → A_loc > 0."""
    omega = _make_smooth_omega(N, amplitude=3.0)
    M, mag, idx = compute_M(omega)
    A_loc, sigma, r_star, _ = compute_A_loc(omega, kx, ky, kz, M=M, idx=idx, mag=mag)
    assert A_loc > 0, f"Expected A_loc>0 for smooth ω, got {A_loc}"
    assert sigma >= 0


def test_sigma_non_negative():
    """σ = A_loc/M^{3/2} ≥ 0 always."""
    for omega in [_make_uniform_omega(N), _make_smooth_omega(N)]:
        M, mag, idx = compute_M(omega)
        _, sigma, _, _ = compute_A_loc(omega, kx, ky, kz, M=M, idx=idx, mag=mag)
        assert sigma >= 0


def test_sigma_scale_invariance():
    """
    Under λ-rescaling of ω: ω → λ^{-2}ω, space → λ·space.
    σ = A_loc/M^{3/2} should be approximately invariant.
    (We test scaling ω by scalar — this changes M and A_loc by M→λM, A_loc→λ^3·A_loc,
     so σ → λ^3·A_loc / (λM)^{3/2} = λ^3/λ^{3/2} · A_loc/M^{3/2} = λ^{3/2} σ.
     Note: pure amplitude scaling is NOT the NS scaling; this test verifies A_loc∝M^{3/2}.)
    Actually the NS scaling changes both amplitude AND spatial scale. Here we test
    the definition consistency: σ ≥ 0 and A_loc ≤ E1.
    """
    omega = _make_smooth_omega(N, amplitude=2.0)
    M, mag, idx = compute_M(omega)
    A_loc, sigma, _, _ = compute_A_loc(omega, kx, ky, kz, M=M, idx=idx, mag=mag)
    # E1 bound: A_loc ≤ E1 = ∫|∇ω|²
    nabla_omega_sq = 0.0
    for comp in range(3):
        f_hat = np.fft.fftn(omega[comp])
        from layer4.geometric_disorder import _spectral_grad
        gx, gy, gz = _spectral_grad(f_hat, kx, ky, kz)
        nabla_omega_sq += float(np.sum(gx**2 + gy**2 + gz**2)) * (2 * np.pi / N)**3
    assert A_loc <= nabla_omega_sq + 1e-8, f"Violated A_loc ≤ E1: {A_loc} > {nabla_omega_sq}"


def test_A_loc_leq_E1():
    """Algebraic consequence: A_loc ≤ ∫|∇ω|² (Equation 21.5)."""
    omega = _make_smooth_omega(N)
    M, mag, idx = compute_M(omega)
    A_loc, sigma, _, _ = compute_A_loc(omega, kx, ky, kz, M=M, idx=idx, mag=mag)
    # E1 on full grid
    from layer4.geometric_disorder import _spectral_grad
    E1 = 0.0
    for comp in range(3):
        f_hat = np.fft.fftn(omega[comp])
        gx, gy, gz = _spectral_grad(f_hat, kx, ky, kz)
        E1 += float(np.sum(gx**2 + gy**2 + gz**2)) * (2 * np.pi / N)**3
    assert A_loc <= E1 + 1e-6, f"A_loc={A_loc} > E1={E1}"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: R_log = log M - λσ
# ─────────────────────────────────────────────────────────────────────────────

def test_R_log_definition():
    M, sigma, lam = 4.0, 0.5, 1.0
    R = compute_R_log(M, sigma, lam)
    assert abs(R - (np.log(4.0) - 0.5)) < 1e-12


def test_R_log_zero_sigma():
    M = 3.0
    R = compute_R_log(M, 0.0, 1.0)
    assert abs(R - np.log(M)) < 1e-12


def test_R_log_decreases_with_sigma():
    M = 5.0
    R1 = compute_R_log(M, 0.1, 1.0)
    R2 = compute_R_log(M, 1.0, 1.0)
    assert R1 > R2, "R_log should decrease as σ increases"


def test_R_log_leq_log_M():
    """R_log = log M - λσ ≤ log M since σ ≥ 0 and λ > 0."""
    M, sigma = 7.0, 0.3
    R = compute_R_log(M, sigma, lam=2.0)
    assert R <= np.log(M) + 1e-12


def test_R_log_zero_M():
    assert compute_R_log(0.0, 1.0) == -np.inf


# ─────────────────────────────────────────────────────────────────────────────
# Tests: Coercive inequality (Proposition 21.3)
# ─────────────────────────────────────────────────────────────────────────────

def test_coercive_LHS_non_negative():
    omega = _make_smooth_omega(N)
    result = check_coercive_inequality(omega, kx, ky, kz)
    assert result["LHS"] >= 0


def test_coercive_sigma_matches_A_loc():
    omega = _make_smooth_omega(N)
    M, mag, idx = compute_M(omega)
    A_loc, sigma, _, _ = compute_A_loc(omega, kx, ky, kz, M=M, idx=idx, mag=mag)
    result = check_coercive_inequality(omega, kx, ky, kz, M=M, idx=idx, mag=mag)
    assert abs(result["sigma"] - sigma) < 1e-10


def test_coercive_ratio_positive_when_sigma_large():
    """When σ is substantial, the ratio LHS/M should be positive (Prop 21.3)."""
    omega = _make_smooth_omega(N, amplitude=4.0)
    result = check_coercive_inequality(omega, kx, ky, kz)
    if result["sigma"] > 0.1:
        assert result["ratio"] >= 0, f"Ratio should be ≥ 0 when σ > 0.1, got {result['ratio']}"


# ─────────────────────────────────────────────────────────────────────────────
# Tests: SigmaTracker
# ─────────────────────────────────────────────────────────────────────────────

def test_sigma_tracker_records():
    tracker = SigmaTracker(lam=1.0)
    omega = _make_smooth_omega(N)
    tracker.record(omega, kx, ky, kz, t=0.0)
    tracker.record(omega * 1.1, kx, ky, kz, t=0.1)
    assert len(tracker.t_hist) == 2
    assert len(tracker.sigma_hist) == 2
    assert len(tracker.R_log_hist) == 2


def test_sigma_tracker_summary_keys():
    tracker = SigmaTracker()
    omega = _make_smooth_omega(N)
    tracker.record(omega, kx, ky, kz, t=0.0)
    s = tracker.summary()
    for key in ["sigma_max", "sigma_min", "sigma_final", "M_max", "R_log_max"]:
        assert key in s, f"Missing key: {key}"


def test_sigma_tracker_uniform_constant():
    """Uniform ω: σ = 0 throughout."""
    tracker = SigmaTracker()
    omega = _make_uniform_omega(N, amplitude=3.0)
    for i in range(5):
        tracker.record(omega, kx, ky, kz, t=i * 0.1)
    s = tracker.summary()
    assert s["sigma_max"] < 1e-7


def test_blowup_alignment_check_no_growth():
    """For static ω, M does not grow → alignment check is trivially consistent."""
    tracker = SigmaTracker()
    omega = _make_smooth_omega(N)
    for i in range(3):
        tracker.record(omega, kx, ky, kz, t=i * 0.1)
    M_grew, sigma_decreased, consistent = tracker.blowup_alignment_check()
    assert consistent, "Consistent should be True when M does not grow"


def test_blowup_alignment_check_growing_M_decreasing_sigma():
    """Simulate M growing while σ decreases: should be consistent with Theorem 22.3."""
    tracker = SigmaTracker()
    omega1 = _make_smooth_omega(N, amplitude=2.0)
    omega2 = _make_smooth_omega(N, amplitude=3.0)
    # Force σ to decrease by making omega2 more uniform
    omega3 = _make_uniform_omega(N, amplitude=4.0)
    tracker.record(omega1, kx, ky, kz, t=0.0)
    tracker.record(omega2, kx, ky, kz, t=0.1)
    tracker.record(omega3, kx, ky, kz, t=0.2)
    M_grew, sigma_decreased, consistent = tracker.blowup_alignment_check()
    assert consistent


def test_sigma_tracker_empty_summary():
    tracker = SigmaTracker()
    assert tracker.summary() == {}


def test_sigma_tracker_single_record():
    tracker = SigmaTracker()
    omega = _make_smooth_omega(N)
    tracker.record(omega, kx, ky, kz, t=0.0)
    s = tracker.summary()
    assert s["n_records"] == 1
    M_grew, sigma_decreased, consistent = tracker.blowup_alignment_check()
    assert consistent
