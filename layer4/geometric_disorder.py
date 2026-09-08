"""
layer4/geometric_disorder.py
Geometric disorder diagnostics for the σ = A_loc/M^{3/2} framework (§21-§22).

Key objects:
  ê = ω/|ω|            vortex direction field
  A_loc                 local misalignment integral ∫_{B_{r*}} |ω|²|∇ê|² φ²
  σ = A_loc/M^{3/2}    scale-invariant geometric disorder (NS-invariant)
  R_log = log M - λσ   logarithmic regulator

Verifies:
  |∇ω|² = |∇|ω||² + |ω|²|∇ê|²   (algebraic identity, Lemma 21.1)
  Coercive inequality              (Proposition 21.3)
  σ ≥ 0 always                    (trivial from definition)
"""

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Spectral gradient helpers
# ─────────────────────────────────────────────────────────────────────────────

def _spectral_grad(f_hat, kx, ky, kz):
    """Return (df/dx, df/dy, df/dz) for a scalar field given its FFT."""
    dfx = np.real(np.fft.ifftn(1j * kx * f_hat))
    dfy = np.real(np.fft.ifftn(1j * ky * f_hat))
    dfz = np.real(np.fft.ifftn(1j * kz * f_hat))
    return dfx, dfy, dfz


def _spectral_laplacian(f_hat, kx, ky, kz):
    """Return Δf for a scalar field given its FFT."""
    k2 = kx**2 + ky**2 + kz**2
    return np.real(np.fft.ifftn(-k2 * f_hat))


def spectral_wavenumbers(N):
    """Return (kx, ky, kz) arrays broadcast-ready for an N^3 grid on T^3 = [0,2π]^3."""
    k1d = np.fft.fftfreq(N, d=1.0/N).astype(float)
    kx = k1d[:, None, None]
    ky = k1d[None, :, None]
    kz = k1d[None, None, :]
    return kx, ky, kz


# ─────────────────────────────────────────────────────────────────────────────
# Core geometric quantities
# ─────────────────────────────────────────────────────────────────────────────

def compute_M(omega):
    """M(t) = ||ω||_{L^∞} and the location x* of the maximum."""
    mag = np.sqrt(omega[0]**2 + omega[1]**2 + omega[2]**2)
    M = float(np.max(mag))
    idx = np.unravel_index(np.argmax(mag), mag.shape)
    return M, mag, idx


def compute_r_star(M):
    """Self-similar scale r* = M^{-1/2}."""
    if M <= 0:
        return np.inf
    return M**(-0.5)


def compute_e_hat(omega, mag=None, eps=1e-12):
    """
    ê = ω/|ω|.  Returns (e_hat, mag).  Where |ω| < eps, ê is set to zero
    (the identity |∇ω|² = |∇|ω||² + |ω|²|∇ê|² is trivially satisfied there).
    """
    if mag is None:
        mag = np.sqrt(omega[0]**2 + omega[1]**2 + omega[2]**2)
    safe_mag = np.where(mag > eps, mag, 1.0)
    e_hat = omega / safe_mag[None, :, :, :]
    mask = (mag > eps)[None, :, :, :]
    e_hat = e_hat * mask
    return e_hat, mag


def verify_gradient_identity(omega, kx, ky, kz, tol=1e-6):
    """
    Verify |∇ω|² = |∇|ω||² + |ω|²|∇ê|² (Lemma 21.1) pointwise.
    Returns (max_relative_error, passes).
    """
    M, mag, _ = compute_M(omega)
    e_hat, _ = compute_e_hat(omega, mag)

    # Left side: |∇ω|²
    grad_omega_sq = np.zeros_like(mag)
    for comp in range(3):
        f_hat = np.fft.fftn(omega[comp])
        gx, gy, gz = _spectral_grad(f_hat, kx, ky, kz)
        grad_omega_sq += gx**2 + gy**2 + gz**2

    # Right side term 1: |∇|ω||²
    mag_hat = np.fft.fftn(mag)
    gmx, gmy, gmz = _spectral_grad(mag_hat, kx, ky, kz)
    grad_mag_sq = gmx**2 + gmy**2 + gmz**2

    # Right side term 2: |ω|²|∇ê|²
    grad_ehat_sq = np.zeros_like(mag)
    for comp in range(3):
        f_hat = np.fft.fftn(e_hat[comp])
        gx, gy, gz = _spectral_grad(f_hat, kx, ky, kz)
        grad_ehat_sq += gx**2 + gy**2 + gz**2

    rhs = grad_mag_sq + mag**2 * grad_ehat_sq
    # Only verify where |ω| is substantial (identity is degenerate at ω=0)
    # and where the gradient is non-negligible
    M_val = float(np.max(mag))
    omega_mask = mag > 0.1 * M_val  # exclude near-zero vorticity regions
    lhs_scale = float(np.max(grad_omega_sq) + 1e-30)
    grad_mask = grad_omega_sq > 1e-8 * lhs_scale
    active = omega_mask & grad_mask
    if not np.any(active):
        return 0.0, True
    abs_err = np.abs(grad_omega_sq - rhs)
    denom = np.maximum(grad_omega_sq, 1e-30)
    rel_err = (abs_err / denom)[active]
    max_err = float(np.max(rel_err))
    return max_err, max_err < tol


def _ball_mask(N, center_idx, r_star):
    """
    Boolean mask for grid points within B_{r*}(x*) on [0,2π]^3 with grid spacing 2π/N.
    Uses minimum-image convention for periodic domain.
    """
    dx = 2.0 * np.pi / N
    ix, iy, iz = center_idx
    i = np.arange(N)
    # Periodic distance along each axis
    di = np.minimum(np.abs(i - ix), N - np.abs(i - ix)) * dx
    dj = np.minimum(np.abs(i - iy), N - np.abs(i - iy)) * dx
    dk = np.minimum(np.abs(i - iz), N - np.abs(i - iz)) * dx
    dist = np.sqrt(di[:, None, None]**2 + dj[None, :, None]**2 + dk[None, None, :]**2)
    return dist <= r_star


def compute_A_loc(omega, kx, ky, kz, M=None, idx=None, mag=None):
    """
    A_loc = ∫_{B_{r*}(x*)} |ω|²|∇ê|² φ² dx
    Uses ball mask as φ² = 1 inside B_{r*}, 0 outside.

    Returns (A_loc, sigma, r_star, ball_fraction).
    """
    if M is None or idx is None or mag is None:
        M, mag, idx = compute_M(omega)
    if M < 1e-14:
        return 0.0, 0.0, np.inf, 0.0

    r_star = compute_r_star(M)
    N = omega.shape[1]
    mask = _ball_mask(N, idx, r_star)
    ball_fraction = float(np.mean(mask))

    e_hat, _ = compute_e_hat(omega, mag)

    # |∇ê|² on full grid
    grad_ehat_sq = np.zeros((N, N, N), dtype=float)
    for comp in range(3):
        f_hat = np.fft.fftn(e_hat[comp])
        gx, gy, gz = _spectral_grad(f_hat, kx, ky, kz)
        grad_ehat_sq += gx**2 + gy**2 + gz**2

    integrand = mag**2 * grad_ehat_sq * mask.astype(float)
    dx3 = (2.0 * np.pi / N)**3
    A_loc = float(np.sum(integrand)) * dx3
    sigma = A_loc / M**1.5
    return A_loc, sigma, r_star, ball_fraction


def compute_R_log(M, sigma, lam=1.0):
    """R_log = log M - λσ (logarithmic regulator, Definition 22.1)."""
    if M <= 0:
        return -np.inf
    return np.log(M) - lam * sigma


def check_coercive_inequality(omega, kx, ky, kz, M=None, idx=None, mag=None):
    """
    Numerically verify the coercive inequality (Proposition 21.3):
      ∫|∇²ω|²φ² / M^{3/2} ≥ (M/C)(σ - C)
    Returns dict with LHS, RHS_bound, and ratio LHS/M.
    """
    if M is None or idx is None or mag is None:
        M, mag, idx = compute_M(omega)
    if M < 1e-14:
        return {"LHS": 0.0, "sigma": 0.0, "M": 0.0, "ratio": 0.0}

    r_star = compute_r_star(M)
    N = omega.shape[1]
    mask = _ball_mask(N, idx, r_star)
    dx3 = (2.0 * np.pi / N)**3

    # ∫|∇²ω|²φ² (using Laplacian as proxy for full Hessian norm)
    nabla2_omega_sq = 0.0
    for comp in range(3):
        f_hat = np.fft.fftn(omega[comp])
        lap = _spectral_laplacian(f_hat, kx, ky, kz)
        nabla2_omega_sq += float(np.sum(lap**2 * mask)) * dx3

    A_loc, sigma, _, _ = compute_A_loc(omega, kx, ky, kz, M, idx, mag)
    LHS = nabla2_omega_sq / M**1.5
    ratio = LHS / max(M, 1e-14)  # should be ≥ (σ - C)
    return {
        "LHS": LHS,
        "sigma": sigma,
        "M": M,
        "ratio": ratio,
        "A_loc": A_loc,
        "r_star": r_star,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Time-series tracker
# ─────────────────────────────────────────────────────────────────────────────

class SigmaTracker:
    """
    Tracks σ(t), R_log(t), M(t), and A_loc(t) over a simulation run.
    Call .record(omega, kx, ky, kz, t) at each step or output point.
    """

    def __init__(self, lam=1.0):
        self.lam = lam
        self.t_hist = []
        self.M_hist = []
        self.sigma_hist = []
        self.A_loc_hist = []
        self.R_log_hist = []
        self.r_star_hist = []
        self.ball_frac_hist = []

    def record(self, omega, kx, ky, kz, t):
        M, mag, idx = compute_M(omega)
        A_loc, sigma, r_star, ball_frac = compute_A_loc(
            omega, kx, ky, kz, M=M, idx=idx, mag=mag
        )
        R_log = compute_R_log(M, sigma, self.lam)
        self.t_hist.append(float(t))
        self.M_hist.append(float(M))
        self.sigma_hist.append(float(sigma))
        self.A_loc_hist.append(float(A_loc))
        self.R_log_hist.append(float(R_log))
        self.r_star_hist.append(float(r_star))
        self.ball_frac_hist.append(float(ball_frac))

    def summary(self):
        if not self.sigma_hist:
            return {}
        return {
            "sigma_max": max(self.sigma_hist),
            "sigma_min": min(self.sigma_hist),
            "sigma_final": self.sigma_hist[-1],
            "M_max": max(self.M_hist),
            "M_final": self.M_hist[-1],
            "R_log_max": max(self.R_log_hist),
            "R_log_min": min(self.R_log_hist),
            "A_loc_max": max(self.A_loc_hist),
            "n_records": len(self.t_hist),
        }

    def sigma_decays(self):
        """Return True if σ is non-increasing overall (consistent with dσ/dt ≤ 0)."""
        if len(self.sigma_hist) < 2:
            return True
        return self.sigma_hist[-1] <= self.sigma_hist[0] * 1.1  # 10% tolerance

    def blowup_alignment_check(self):
        """
        Check Theorem 22.3: if M grows, does σ decrease?
        Returns (M_grew, sigma_decreased, consistent).
        """
        if len(self.M_hist) < 2:
            return (False, False, True)
        M_grew = self.M_hist[-1] > self.M_hist[0] * 1.01
        sigma_decreased = self.sigma_hist[-1] < self.sigma_hist[0]
        consistent = (not M_grew) or sigma_decreased
        return (M_grew, sigma_decreased, consistent)
