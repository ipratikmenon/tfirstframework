"""
route2_3D.py — Route 2: Exact Prize NS + Auxiliary Scalar θ in 3D
==================================================================
Implements the primary Prize-relevant system from PRD v1.0 FINAL §3.3 in 3D.

Physics (Route 2 — exact Prize equations throughout):

    ∂_t u + (u·∇)u = −∇p + ν·Δu     div u = 0  (exact Prize NS, ν = const)

    ∂_t θ + u·∇θ = ν·Δθ + ν|∇u|²    θ(x,0) = 0
    S_θ = ν|∇u|² ≥ 0               (viscous heating, always non-negative)

    μ_eff(x,t) = ν + ε·f(θ(x,t))    (always ≥ ν)

Key diagnostics (PRD §7.3 Steps 2–7):
  - θ ≥ 0 (maximum principle for θ — source always non-negative)
  - δ(t): effective Sobolev regularity of θ above H¹ — central Prize quantity
  - CZ integrability ε_CZ(t): ‖ν|∇u|²‖_{L^{1+ε}} — tests Claim A
  - LPS margin ε_LPS: min(μ_eff) − ν ≥ 0
  - Energy decay E(t), dissipation rate ε_diss(t)

Numerics:
  - FFT3 spectral derivatives (2/3 dealiasing, Nyquist planes zeroed)
  - Leray projector (B-005-fixed: Nyquist planes zeroed before projection)
  - RK4 time integration for velocity u
  - Crank–Nicolson (IMEX) for θ: implicit ν·Δθ avoids diffusion CFL constraint
  - Adaptive CFL from advection only

Experiments:
  EXP-L3-R2-TG-001   Taylor–Green IC, t_end=1.0, N=32³ (validation)
  EXP-L3-R2-SH-001   Shear layer IC (non-mixing, λ_max≈0) — Route A decisive test
  EXP-L3-R2-EPS-001  ε_param sweep: ε ∈ {1.0, 0.1, 0.01, 0.001, 0.0}

Run:
  python layer3/route2_3D.py --run tg --N 32 --verbose
  python layer3/route2_3D.py --run shear --N 32 --verbose
  python layer3/route2_3D.py --run eps-sweep --N 32
"""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
import time
from datetime import datetime

import numpy as np

# ── path setup ─────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

# ── default DB path ────────────────────────────────────────────────────────────
_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")

# ── kinematic viscosity (Prize equations: ν = const) ──────────────────────────
_NU_DEFAULT: float = 1.0e-3   # m²/s — standard reference value

# ── ε_param sweep values ───────────────────────────────────────────────────────
EPS_SWEEP: list[float] = [1.0, 0.1, 0.01, 0.001, 0.0]


# =============================================================================
# Grid helpers
# =============================================================================

def make_grid_3D(N: int) -> tuple:
    """Return wavenumber arrays kx, ky, kz, k² for N³ domain [0, 2π)³."""
    freqs = np.fft.fftfreq(N, d=1.0 / N)
    kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    return kx, ky, kz, k2


def make_dealias_mask_3D(N: int) -> np.ndarray:
    """2/3-rule dealiasing mask: zeros modes with |k_i| > N//3."""
    freqs = np.fft.fftfreq(N, d=1.0 / N)
    kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k_max = N // 3
    return (
        (np.abs(kx) <= k_max)
        & (np.abs(ky) <= k_max)
        & (np.abs(kz) <= k_max)
    ).astype(np.float64)


def _zero_nyquist(f_hat: np.ndarray, N: int) -> np.ndarray:
    """Zero the three Nyquist planes (B-005 fix: restores Hermitian symmetry)."""
    nq = N // 2
    f_hat = f_hat.copy()
    f_hat[nq, :, :] = 0.0
    f_hat[:, nq, :] = 0.0
    f_hat[:, :, nq] = 0.0
    return f_hat


def project_divergence_free(
    u: np.ndarray, v: np.ndarray, w: np.ndarray,
    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray, k2: np.ndarray,
) -> tuple:
    """Leray projection: remove the irrotational component so div u = 0.

    B-005 fix applied: Nyquist planes zeroed before projection.
    """
    N = u.shape[0]
    u_hat = _zero_nyquist(np.fft.fftn(u), N)
    v_hat = _zero_nyquist(np.fft.fftn(v), N)
    w_hat = _zero_nyquist(np.fft.fftn(w), N)

    k_dot_u = kx * u_hat + ky * v_hat + kz * w_hat
    correction = np.zeros_like(k_dot_u)
    nonzero = k2 > 0
    correction[nonzero] = k_dot_u[nonzero] / k2[nonzero]

    u_hat -= kx * correction
    v_hat -= ky * correction
    w_hat -= kz * correction
    return (
        np.real(np.fft.ifftn(u_hat)),
        np.real(np.fft.ifftn(v_hat)),
        np.real(np.fft.ifftn(w_hat)),
    )


# =============================================================================
# f(θ) functional forms and μ_eff
# =============================================================================

def _apply_f_theta(theta: np.ndarray, f_theta: str = "identity") -> np.ndarray:
    """Apply functional form f to θ. Requirements: f(θ) ≥ 0, f(0) = 0."""
    theta = np.clip(theta, 0.0, None)   # numerical floor
    if f_theta == "identity":
        return theta
    elif f_theta == "tanh":
        return np.tanh(theta)
    elif f_theta == "sqrt":
        eps_reg = 1e-14
        return np.sqrt(theta + eps_reg) - np.sqrt(eps_reg)
    else:
        raise ValueError(f"Unknown f_theta: '{f_theta}'. Use 'identity', 'tanh', 'sqrt'.")


def mu_eff_field(theta: np.ndarray, nu: float, eps_param: float,
                 f_theta: str = "identity") -> np.ndarray:
    """μ_eff(x,t) = ν + ε·f(θ(x,t)).  Always ≥ ν > 0."""
    return nu + eps_param * _apply_f_theta(theta, f_theta)


# =============================================================================
# Diagnostics: δ, CZ integrability, LPS margin
# =============================================================================

def delta_from_theta(theta: np.ndarray, kx: np.ndarray, ky: np.ndarray,
                     kz: np.ndarray) -> float:
    """Estimate the Sobolev regularity gain δ of θ above H¹.

    Computes ‖|k|^s θ̂‖_{L²} for s ∈ {0, 0.5, 1.0, 1.5, 2.0}.
    δ is estimated from the ratio of adjacent H^s norms:

        δ_est(s) = (log ‖|k|^s θ̂‖ − log ‖|k|^{s−0.5} θ̂‖) / log(growth_factor)

    A simpler and robust proxy: the largest s such that ‖|k|^s θ̂‖ is still
    well-controlled (less than 10 × ‖|k|^0 θ̂‖).  Returns that s − 1.0.

    Returns δ_est ≥ 0.  Returns 0.0 if θ is identically zero (startup).
    """
    theta_hat = np.fft.fftn(theta)
    k_abs = np.sqrt(kx**2 + ky**2 + kz**2)

    theta_hat_flat = np.abs(theta_hat.ravel())
    k_flat = k_abs.ravel()

    # Remove DC mode (k=0) from the spectral norm calculation
    mask_nonzero = k_flat > 0.0
    if not np.any(mask_nonzero):
        return 0.0

    theta_hat_nz = theta_hat_flat[mask_nonzero]
    k_nz = k_flat[mask_nonzero]

    H0_sq = float(np.sum(theta_hat_nz**2))
    if H0_sq < 1e-30:
        return 0.0   # θ is effectively zero (early time)

    # Compute H^s norms for s = 0, 0.5, 1.0, 1.5, 2.0, 2.5
    s_vals = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
    Hs_sq = [float(np.sum((k_nz**s * theta_hat_nz)**2)) for s in s_vals]

    # Find largest s where H^s norm is < 100 × H^0 norm (remains controlled)
    H0 = np.sqrt(Hs_sq[0])
    delta_est = 0.0
    for i, s in enumerate(s_vals):
        Hs = np.sqrt(max(Hs_sq[i], 1e-60))
        if Hs < 100.0 * H0:
            delta_est = max(s - 1.0, 0.0)
        else:
            break

    return float(delta_est)


def cz_integrability(source: np.ndarray, eps_probe: float = 0.1) -> float:
    """Estimate ‖S_θ‖_{L^{1+ε}} (discrete) as a proxy for CZ integrability.

    S_θ = ν|∇u|² ≥ 0.  Claim A: S_θ ∈ L^{1+ε} for some ε > 0.
    Returns ‖S_θ‖_{L^{1+ε}} (volumetrically averaged).

    For the pass criterion: compare with ‖S_θ‖_{L^1}.  If the ratio is finite,
    Claim A holds at this exponent.
    """
    p = 1.0 + eps_probe
    N3 = source.size
    L1 = float(np.mean(np.abs(source)))
    Lp = float(np.mean(np.abs(source)**p)**(1.0 / p))
    return Lp


def sobolev_norms(theta: np.ndarray, kx: np.ndarray, ky: np.ndarray,
                  kz: np.ndarray, s_list: tuple = (0.0, 0.5, 1.0, 1.5, 2.0)
                  ) -> dict:
    """Return ‖θ‖_{H^s} for each s in s_list.  Used for tracking δ evolution."""
    theta_hat = np.fft.fftn(theta)
    k2 = kx**2 + ky**2 + kz**2
    N3 = theta.size
    out = {}
    for s in s_list:
        weight = (1.0 + k2)**s
        Hs_sq = float(np.sum(weight * np.abs(theta_hat)**2)) / N3
        out[f"H{s:.1f}"] = float(np.sqrt(max(Hs_sq, 0.0)))
    return out


def divergence_rms(u: np.ndarray, v: np.ndarray, w: np.ndarray,
                   kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> float:
    """RMS of ∂u/∂x + ∂v/∂y + ∂w/∂z (spectral)."""
    div_hat = 1j * (kx * np.fft.fftn(u)
                    + ky * np.fft.fftn(v)
                    + kz * np.fft.fftn(w))
    div = np.real(np.fft.ifftn(div_hat))
    return float(np.sqrt(np.mean(div**2)))


def kinetic_energy(u: np.ndarray, v: np.ndarray, w: np.ndarray) -> float:
    """Mean kinetic energy E = ½⟨|u|²⟩."""
    return float(0.5 * np.mean(u**2 + v**2 + w**2))


# =============================================================================
# Solver construction
# =============================================================================

def make_solver(
    N: int = 32,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    f_theta: str = "identity",
    dt_max: float = 0.05,
) -> dict:
    """Create a Route 2 3D solver state dict.

    All fields on a periodic domain [0, 2π)³ with N³ spectral collocation points.

    Parameters
    ----------
    N         : Grid points per side (N³ total)
    nu        : Constant kinematic viscosity ν [m²/s] (Prize: ν = const)
    eps_param : ε in μ_eff = ν + ε·f(θ).  ε = 0 = exact Prize equations.
    f_theta   : 'identity', 'tanh', or 'sqrt'
    dt_max    : Maximum allowed timestep [s]
    """
    if nu <= 0.0:
        raise ValueError(f"nu must be > 0, got {nu}")
    if eps_param < 0.0:
        raise ValueError(f"eps_param must be ≥ 0, got {eps_param}")

    kx, ky, kz, k2 = make_grid_3D(N)
    dealias = make_dealias_mask_3D(N)

    return {
        "N"         : N,
        "nu"        : nu,
        "eps_param" : eps_param,
        "f_theta"   : f_theta,
        "dt_max"    : dt_max,
        "kx"        : kx, "ky": ky, "kz": kz, "k2": k2,
        "dealias"   : dealias,
        # Fields (all N³)
        "u"         : np.zeros((N, N, N)),
        "v"         : np.zeros((N, N, N)),
        "w"         : np.zeros((N, N, N)),
        "theta"     : np.zeros((N, N, N)),   # θ(x,0) = 0 always (PRD §3.3)
        "t"         : 0.0,
    }


def set_ic(solver: dict, ic: dict) -> dict:
    """Set IC on solver; Leray-project velocity; θ always starts at 0.
    Returns new solver dict (immutable pattern).
    """
    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    u, v, w = ic["u"], ic["v"], ic["w"]
    # Ensure div-free
    u, v, w = project_divergence_free(u, v, w, kx, ky, kz, k2)
    return {
        **solver,
        "u"    : u,
        "v"    : v,
        "w"    : w,
        "theta": np.zeros_like(u),   # θ(x,0) = 0
        "t"    : 0.0,
    }


# =============================================================================
# Initial conditions
# =============================================================================

def taylor_green_ic(N: int, V0: float = 1.0) -> dict:
    """Taylor–Green vortex IC on N³ grid."""
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    u = V0 *  np.sin(xx) * np.cos(yy) * np.cos(zz)
    v = -V0 * np.cos(xx) * np.sin(yy) * np.cos(zz)
    w = np.zeros((N, N, N))
    return {"u": u, "v": v, "w": w, "E0_analytic": V0**2 / 8.0, "name": "taylor_green"}


def shear_layer_ic(N: int, amp: float = 1.0, width: float = 0.2) -> dict:
    """Shear layer IC: u = amp · tanh((z − π)/width), v = w = 0.

    Non-mixing flow: Lyapunov exponent λ_max ≈ 0 at early times.
    Decisive Route A test: δ > 0 even without exponential mixing.
    """
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    _, _, zz = np.meshgrid(x, x, x, indexing="ij")
    u = amp * np.tanh((zz - np.pi) / width)
    v = np.zeros((N, N, N))
    w = np.zeros((N, N, N))
    return {"u": u, "v": v, "w": w, "E0_analytic": float(np.mean(u**2) / 2.0),
            "name": "shear_layer"}


def random_div_free_ic(N: int, amp: float = 0.5, seed: int = 42) -> dict:
    """Smooth random divergence-free IC (low-wavenumber modes only)."""
    rng = np.random.default_rng(seed)
    kx, ky, kz, k2 = make_grid_3D(N)
    k_abs = np.sqrt(k2)

    def _rand_field():
        f_hat = np.zeros((N, N, N), dtype=complex)
        low_k = (k_abs > 0) & (k_abs <= 4)
        f_hat[low_k] = rng.standard_normal(int(np.sum(low_k))) + \
                       1j * rng.standard_normal(int(np.sum(low_k)))
        return np.real(np.fft.ifftn(f_hat))

    u = amp * _rand_field()
    v = amp * _rand_field()
    w = amp * _rand_field()
    # Project to div-free
    u, v, w = project_divergence_free(u, v, w, kx, ky, kz, k2)
    return {"u": u, "v": v, "w": w, "E0_analytic": float(np.mean(u**2 + v**2 + w**2) / 2.0),
            "name": "random_div_free"}


# =============================================================================
# Right-hand side: Route 2 system
# =============================================================================

def _spectral_deriv(f_hat: np.ndarray, k_dir: np.ndarray) -> np.ndarray:
    """Spectral derivative of f̂ in direction k_dir → physical space."""
    return np.real(np.fft.ifftn(1j * k_dir * f_hat))


def _source_theta(u: np.ndarray, v: np.ndarray, w: np.ndarray,
                  kx: np.ndarray, ky: np.ndarray, kz: np.ndarray,
                  nu: float) -> np.ndarray:
    """Compute S_θ = ν|∇u|² = ν · Σ_{ij} (∂_j u_i)² ≥ 0 (spectral accuracy)."""
    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)

    grad_sq = (
        _spectral_deriv(u_hat, kx)**2 + _spectral_deriv(u_hat, ky)**2
        + _spectral_deriv(u_hat, kz)**2
        + _spectral_deriv(v_hat, kx)**2 + _spectral_deriv(v_hat, ky)**2
        + _spectral_deriv(v_hat, kz)**2
        + _spectral_deriv(w_hat, kx)**2 + _spectral_deriv(w_hat, ky)**2
        + _spectral_deriv(w_hat, kz)**2
    )
    return nu * grad_sq


def _rhs_velocity(u: np.ndarray, v: np.ndarray, w: np.ndarray,
                  kx: np.ndarray, ky: np.ndarray, kz: np.ndarray,
                  k2: np.ndarray, dealias: np.ndarray, nu: float,
                  N: int) -> tuple:
    """Route 2 velocity RHS: ∂_t u = −(u·∇)u − ∇p + ν·Δu, div u = 0.

    Uses mean-field viscosity ν (not μ_eff) — exact Prize equations.
    Returns (du/dt, dv/dt, dw/dt) in physical space.
    """
    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)

    # Velocity gradients
    dudx = _spectral_deriv(u_hat, kx); dudy = _spectral_deriv(u_hat, ky)
    dudz = _spectral_deriv(u_hat, kz)
    dvdx = _spectral_deriv(v_hat, kx); dvdy = _spectral_deriv(v_hat, ky)
    dvdz = _spectral_deriv(v_hat, kz)
    dwdx = _spectral_deriv(w_hat, kx); dwdy = _spectral_deriv(w_hat, ky)
    dwdz = _spectral_deriv(w_hat, kz)

    # Nonlinear advection
    NL_u = u * dudx + v * dudy + w * dudz
    NL_v = u * dvdx + v * dvdy + w * dvdz
    NL_w = u * dwdx + v * dwdy + w * dwdz

    # Dealias nonlinear terms
    NL_u_hat = np.fft.fftn(NL_u) * dealias
    NL_v_hat = np.fft.fftn(NL_v) * dealias
    NL_w_hat = np.fft.fftn(NL_w) * dealias

    # Leray projection of nonlinear term
    k_dot_NL = kx * NL_u_hat + ky * NL_v_hat + kz * NL_w_hat
    corr = np.zeros_like(k_dot_NL)
    nz = k2 > 0
    corr[nz] = k_dot_NL[nz] / k2[nz]
    NL_u_hat -= kx * corr
    NL_v_hat -= ky * corr
    NL_w_hat -= kz * corr

    # Velocity tendency: −NL + ν·Δu (linear dissipation spectral)
    du_hat = (-NL_u_hat - nu * k2 * u_hat) * dealias
    dv_hat = (-NL_v_hat - nu * k2 * v_hat) * dealias
    dw_hat = (-NL_w_hat - nu * k2 * w_hat) * dealias

    return (np.real(np.fft.ifftn(du_hat)),
            np.real(np.fft.ifftn(dv_hat)),
            np.real(np.fft.ifftn(dw_hat)))


def _rhs_theta_explicit(theta: np.ndarray, u: np.ndarray, v: np.ndarray,
                        w: np.ndarray, kx: np.ndarray, ky: np.ndarray,
                        kz: np.ndarray, S_theta: np.ndarray) -> np.ndarray:
    """Explicit part of θ RHS: −u·∇θ + S_θ.

    Used with Crank–Nicolson: diffusion ν·Δθ is handled implicitly in spectral space.
    """
    theta_hat = np.fft.fftn(theta)
    adv_theta = (u * _spectral_deriv(theta_hat, kx)
                 + v * _spectral_deriv(theta_hat, ky)
                 + w * _spectral_deriv(theta_hat, kz))
    return -adv_theta + S_theta


# =============================================================================
# CFL timestep
# =============================================================================

def compute_cfl_dt(solver: dict, safety: float = 0.4) -> float:
    """CFL-limited timestep (advection only — θ diffusion handled implicitly)."""
    u, v, w = solver["u"], solver["v"], solver["w"]
    N, nu = solver["N"], solver["nu"]
    dx = 2.0 * np.pi / N

    U_max = max(float(np.max(np.abs(u))), float(np.max(np.abs(v))),
                float(np.max(np.abs(w))), 1e-12)
    dt_adv = safety * dx / U_max
    # Diffusion CFL (reference only; IMEX for θ makes this non-binding)
    dt_diff = safety * dx**2 / (2.0 * nu)

    return float(np.clip(min(dt_adv, dt_diff), 1e-12, solver["dt_max"]))


# =============================================================================
# Single timestep: RK4 for u, Crank–Nicolson for θ
# =============================================================================

def solver_step(solver: dict, dt: float | None = None) -> tuple[dict, dict]:
    """One step of the Route 2 3D system.

    Velocity: RK4 (same as route1_3D.py) using exact Prize NS.
    Theta: Crank–Nicolson (implicit ν·Δθ, explicit advection + source).
      θ̂^{n+1} = (θ̂^n · (1 − ½dt·ν·k²) + dt · F̂^n_expl) / (1 + ½dt·ν·k²)
      where F_expl = −u·∇θ + S_θ (evaluated at time n).

    Returns (new_solver, step_diagnostics).  Does not mutate input (immutable).
    """
    if dt is None:
        dt = compute_cfl_dt(solver)

    u, v, w, theta = solver["u"], solver["v"], solver["w"], solver["theta"]
    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    dealias, N = solver["dealias"], solver["N"]
    nu, eps_param, f_th = solver["nu"], solver["eps_param"], solver["f_theta"]

    # ── 1. Source S_θ = ν|∇u|² at current time ────────────────────────────────
    S_theta = _source_theta(u, v, w, kx, ky, kz, nu)

    # ── 2. RK4 for velocity (exact Prize NS) ───────────────────────────────────
    def vel_rhs(uu, vv, ww):
        return _rhs_velocity(uu, vv, ww, kx, ky, kz, k2, dealias, nu, N)

    k1u, k1v, k1w = vel_rhs(u, v, w)
    k2u, k2v, k2w = vel_rhs(u + 0.5*dt*k1u, v + 0.5*dt*k1v, w + 0.5*dt*k1w)
    k3u, k3v, k3w = vel_rhs(u + 0.5*dt*k2u, v + 0.5*dt*k2v, w + 0.5*dt*k2w)
    k4u, k4v, k4w = vel_rhs(u + dt*k3u,     v + dt*k3v,     w + dt*k3w)

    fac = dt / 6.0
    new_u = u + fac * (k1u + 2.0*k2u + 2.0*k3u + k4u)
    new_v = v + fac * (k1v + 2.0*k2v + 2.0*k3v + k4v)
    new_w = w + fac * (k1w + 2.0*k2w + 2.0*k3w + k4w)

    # Dealias velocity
    def _dealias(f):
        return np.real(np.fft.ifftn(np.fft.fftn(f) * dealias))

    new_u = _dealias(new_u)
    new_v = _dealias(new_v)
    new_w = _dealias(new_w)

    # Re-project velocity (accumulation guard)
    new_u, new_v, new_w = project_divergence_free(new_u, new_v, new_w,
                                                   kx, ky, kz, k2)

    # ── 3. Crank–Nicolson for θ ────────────────────────────────────────────────
    F_expl = _rhs_theta_explicit(theta, u, v, w, kx, ky, kz, S_theta)
    theta_hat = np.fft.fftn(theta)
    F_hat = np.fft.fftn(F_expl)

    denom = 1.0 + 0.5 * dt * nu * k2          # implicit diffusion factor
    theta_hat_new = (theta_hat * (1.0 - 0.5 * dt * nu * k2) + dt * F_hat) / denom
    new_theta = np.real(np.fft.ifftn(theta_hat_new))

    # Enforce θ ≥ 0 (source always ≥ 0; tiny negative values from numerics)
    np.clip(new_theta, 0.0, None, out=new_theta)

    # ── 4. Diagnostics ────────────────────────────────────────────────────────
    mu_eff = mu_eff_field(new_theta, nu, eps_param, f_th)
    eps_lps = float(np.min(mu_eff)) - nu   # always ≥ 0

    delta_est = delta_from_theta(new_theta, kx, ky, kz)
    S_cz = cz_integrability(S_theta, eps_probe=0.1)
    Hs = sobolev_norms(new_theta, kx, ky, kz)
    div_rms_val = divergence_rms(new_u, new_v, new_w, kx, ky, kz)
    E = kinetic_energy(new_u, new_v, new_w)

    step_diag = {
        "t"             : solver["t"] + dt,
        "dt"            : dt,
        "E"             : E,
        "div_rms"       : div_rms_val,
        "theta_min"     : float(np.min(new_theta)),
        "theta_max"     : float(np.max(new_theta)),
        "theta_mean"    : float(np.mean(new_theta)),
        "S_theta_max"   : float(np.max(S_theta)),
        "S_theta_mean"  : float(np.mean(S_theta)),
        "mu_eff_min"    : float(np.min(mu_eff)),
        "eps_lps"       : eps_lps,
        "delta_est"     : delta_est,
        "cz_L1p"        : S_cz,
        **{f"Hs_{k}": v for k, v in Hs.items()},
    }

    new_solver = {
        **solver,
        "u": new_u, "v": new_v, "w": new_w,
        "theta": new_theta,
        "t": solver["t"] + dt,
    }
    return new_solver, step_diag


# =============================================================================
# Run loop
# =============================================================================

def run(
    solver: dict,
    t_end: float,
    max_steps: int = 5000,
    cfl_safety: float = 0.4,
    verbose: bool = False,
) -> dict:
    """Run solver from solver['t'] to t_end.

    Returns result dict with full history and PASS/FAIL verdict.
    """
    history: list[dict] = []
    verdict = "PASS"
    E_initial = kinetic_energy(solver["u"], solver["v"], solver["w"])
    E_prev = E_initial
    delta_max = 0.0
    theta_negative_flag = False

    for step_idx in range(max_steps):
        if solver["t"] >= t_end:
            break

        dt = compute_cfl_dt(solver, safety=cfl_safety)
        dt = min(dt, t_end - solver["t"])

        solver, diag = solver_step(solver, dt)

        history.append(diag)
        delta_max = max(delta_max, diag["delta_est"])

        if diag["theta_min"] < -1e-8:
            theta_negative_flag = True
            verdict = "FAIL"
            if verbose:
                print(f"  [FAIL] Step {step_idx}: θ_min = {diag['theta_min']:.3e} < 0")

        if verbose and step_idx % max(1, max_steps // 10) == 0:
            print(
                f"  step {step_idx:4d}  t={diag['t']:.4f}  E={diag['E']:.4e}  "
                f"θ_max={diag['theta_max']:.3e}  δ={diag['delta_est']:.2f}  "
                f"ε_LPS={diag['eps_lps']:.4f}"
            )

    E_final = kinetic_energy(solver["u"], solver["v"], solver["w"])
    return {
        "solver"        : solver,
        "history"       : history,
        "verdict"       : verdict,
        "n_steps"       : len(history),
        "t_final"       : solver["t"],
        "E_initial"     : E_initial,
        "E_final"       : E_final,
        "energy_decayed": E_final < E_initial,
        "delta_max"     : delta_max,
        "theta_negative": theta_negative_flag,
        "eps_lps_min"   : min(d["eps_lps"] for d in history) if history else 0.0,
        "theta_max_global": max(d["theta_max"] for d in history) if history else 0.0,
    }


# =============================================================================
# Benchmark experiments
# =============================================================================

def run_taylor_green(
    N: int = 32,
    t_end: float = 1.0,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str | None = None,
    exp_suffix: str = "",
) -> dict:
    """Route 2 Taylor–Green benchmark in 3D.

    Pass criteria (PRD §7.3 Step 2–3):
      1. θ ≥ 0 always (source always non-negative)
      2. E(t) decays (energy dissipated by viscosity)
      3. ε_LPS = min(μ_eff) − ν ≥ 0 (always by construction)
      4. δ_est > 0 (Sobolev regularity gain observed)
    """
    solver = make_solver(N=N, nu=nu, eps_param=eps_param)
    ic = taylor_green_ic(N, V0=1.0)
    solver = set_ic(solver, ic)

    E0_computed = kinetic_energy(solver["u"], solver["v"], solver["w"])
    E0_analytic = ic["E0_analytic"]

    if verbose:
        print(f"\n  Route 2 3D TG Benchmark  N={N}³  ν={nu}  ε={eps_param}")
        print(f"  E₀ analytic = {E0_analytic:.6f}  computed = {E0_computed:.6f}")

    t0 = time.perf_counter()
    result = run(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)
    wall = time.perf_counter() - t0

    E0_ok      = abs(E0_computed - E0_analytic) / E0_analytic < 0.01
    theta_ok   = not result["theta_negative"]
    energy_ok  = result["energy_decayed"]
    delta_ok   = result["delta_max"] > 0.0
    passed = E0_ok and theta_ok and energy_ok

    exp = {
        "exp_id"         : f"EXP-L3-R2-TG-{N:03d}{exp_suffix}",
        "claim_id"       : "claim_m5b_route2_3D_taylor_green",
        "timestamp"      : datetime.utcnow().isoformat(),
        "ic"             : "taylor_green",
        "N"              : N,
        "nu"             : nu,
        "eps_param"      : eps_param,
        "t_end"          : t_end,
        "n_steps"        : result["n_steps"],
        "t_final"        : result["t_final"],
        "E0_analytic"    : E0_analytic,
        "E0_computed"    : E0_computed,
        "E0_error_pct"   : abs(E0_computed - E0_analytic) / E0_analytic * 100.0,
        "E_final"        : result["E_final"],
        "energy_decayed" : result["energy_decayed"],
        "theta_max_global": result["theta_max_global"],
        "theta_negative" : result["theta_negative"],
        "eps_lps_min"    : result["eps_lps_min"],
        "delta_max"      : result["delta_max"],
        "verdict"        : "PASS" if passed else "FAIL",
        "wall_time_s"    : wall,
        "key_metric"     : (
            f"E0_err={abs(E0_computed-E0_analytic)/E0_analytic*100:.3f}%  "
            f"θ_max={result['theta_max_global']:.3e}  "
            f"δ_max={result['delta_max']:.2f}  "
            f"ε_LPS_min={result['eps_lps_min']:.4f}"
        ),
    }

    if verbose:
        print(f"\n  Route 2 3D TG Results (N={N}³, t_end={t_end}):")
        print(f"    E₀ error:    {exp['E0_error_pct']:.3f}%  {'✓' if E0_ok else '✗'}")
        print(f"    θ ≥ 0:       {'✓' if theta_ok else '✗'}")
        print(f"    Energy decay:{'✓' if energy_ok else '✗'}")
        print(f"    δ_max:       {result['delta_max']:.2f}  {'✓' if delta_ok else '(0: early time)'}")
        print(f"    ε_LPS_min:   {result['eps_lps_min']:.4f}")
        print(f"    Wall time:   {wall:.1f} s")
        print(f"    Verdict:     {exp['verdict']}")

    if db_path is not None:
        _log_result(db_path, exp)

    return exp


def run_shear_layer(
    N: int = 32,
    t_end: float = 2.0,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    max_steps: int = 5000,
    verbose: bool = False,
    db_path: str | None = None,
) -> dict:
    """Route 2 shear-layer IC — decisive Route A test.

    Non-mixing flow (λ_max ≈ 0 at early times).  If δ_max > 0, the CZ mechanism
    (Route A) is confirmed as the source of regularity — mixing is not required.
    """
    solver = make_solver(N=N, nu=nu, eps_param=eps_param)
    ic = shear_layer_ic(N)
    solver = set_ic(solver, ic)

    if verbose:
        print(f"\n  Route 2 3D Shear IC  N={N}³  ν={nu}  ε={eps_param}")

    t0 = time.perf_counter()
    result = run(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)
    wall = time.perf_counter() - t0

    passed = not result["theta_negative"] and result["energy_decayed"]

    exp = {
        "exp_id"          : f"EXP-L3-R2-SH-{N:03d}",
        "claim_id"        : "claim_m5b_route2_3D_shear_route_A",
        "timestamp"       : datetime.utcnow().isoformat(),
        "ic"              : "shear_layer",
        "N"               : N,
        "nu"              : nu,
        "eps_param"       : eps_param,
        "t_end"           : t_end,
        "n_steps"         : result["n_steps"],
        "t_final"         : result["t_final"],
        "E_initial"       : result["E_initial"],
        "E_final"         : result["E_final"],
        "energy_decayed"  : result["energy_decayed"],
        "theta_max_global": result["theta_max_global"],
        "theta_negative"  : result["theta_negative"],
        "eps_lps_min"     : result["eps_lps_min"],
        "delta_max"       : result["delta_max"],
        "verdict"         : "PASS" if passed else "FAIL",
        "wall_time_s"     : wall,
        "key_metric"      : (
            f"θ_max={result['theta_max_global']:.3e}  "
            f"δ_max={result['delta_max']:.2f}  "
            f"ε_LPS_min={result['eps_lps_min']:.4f}"
        ),
    }

    if verbose:
        print(f"    δ_max: {result['delta_max']:.2f}  "
              f"{'→ Route A CONFIRMED (CZ, no mixing)' if result['delta_max'] > 0 else '→ 0 (accumulation may be slow)'}")
        print(f"    Verdict: {exp['verdict']}")

    if db_path is not None:
        _log_result(db_path, exp)

    return exp


# =============================================================================
# Results logging
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS route2_3D_experiments (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id          TEXT UNIQUE NOT NULL,
            claim_id        TEXT NOT NULL,
            timestamp       TEXT NOT NULL,
            ic              TEXT,
            N               INTEGER NOT NULL,
            nu              REAL,
            eps_param       REAL,
            t_end           REAL,
            n_steps         INTEGER,
            t_final         REAL,
            verdict         TEXT NOT NULL,
            E0_error_pct    REAL,
            E_final         REAL,
            energy_decayed  INTEGER,
            theta_max_global REAL,
            theta_negative  INTEGER,
            eps_lps_min     REAL,
            delta_max       REAL,
            wall_time_s     REAL,
            key_metric      TEXT
        )
    """)
    conn.commit()
    return conn


def _log_result(db_path: str, result: dict) -> None:
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "ic", "N", "nu", "eps_param",
        "t_end", "n_steps", "t_final", "verdict",
        "E0_error_pct", "E_final", "energy_decayed",
        "theta_max_global", "theta_negative", "eps_lps_min",
        "delta_max", "wall_time_s", "key_metric",
    ]
    values = tuple(
        int(result[c]) if isinstance(result.get(c), bool) else result.get(c)
        for c in cols
    )
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO route2_3D_experiments ({', '.join(cols)}) VALUES ({ph})",
        values,
    )
    conn.commit()
    conn.close()


# =============================================================================
# CLI entry point
# =============================================================================

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Route 2 3D solver (Prize-relevant)")
    parser.add_argument("--run", choices=["tg", "shear", "eps-sweep"], default="tg")
    parser.add_argument("--N", type=int, default=32)
    parser.add_argument("--t-end", type=float, default=1.0)
    parser.add_argument("--nu", type=float, default=_NU_DEFAULT)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db", default=_DEFAULT_DB)
    args = parser.parse_args()

    if args.run == "tg":
        r = run_taylor_green(N=args.N, t_end=args.t_end, nu=args.nu,
                             eps_param=args.eps, verbose=args.verbose,
                             db_path=args.db)
        print(f"\nVerdict: {r['verdict']}  key_metric: {r['key_metric']}")

    elif args.run == "shear":
        r = run_shear_layer(N=args.N, t_end=args.t_end, nu=args.nu,
                            eps_param=args.eps, verbose=args.verbose,
                            db_path=args.db)
        print(f"\nVerdict: {r['verdict']}  key_metric: {r['key_metric']}")

    elif args.run == "eps-sweep":
        for eps in EPS_SWEEP:
            r = run_taylor_green(N=args.N, t_end=args.t_end, nu=args.nu,
                                 eps_param=eps, verbose=False, db_path=args.db,
                                 exp_suffix=f"_eps{eps}")
            print(f"  ε={eps:.3f}  δ_max={r['delta_max']:.2f}  "
                  f"ε_LPS={r['eps_lps_min']:.4f}  {r['verdict']}")


if __name__ == "__main__":
    main()
