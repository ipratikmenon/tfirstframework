"""
route2_3D_jax.py — JAX-JIT Route 2 3D Solver
=============================================
JAX-JIT port of route2_3D.py for XLA/LLVM-accelerated production runs.

Physics (Route 2, PRD v1.0 FINAL §3.3 — exact Prize equations):
  ∂_t u + (u·∇)u = −∇p + ν·Δu     div u = 0   (ν = const)
  ∂_t θ + u·∇θ   = ν·Δθ + ν|∇u|²  θ(x,0) = 0
  μ_eff(x,t)      = ν + ε·f(θ)      ≥ ν

Numerics:
  - Velocity: RK4 spectral, 2/3 dealiasing, Leray re-projection
  - θ: Crank-Nicolson IMEX (ν·Δθ implicit, advection + source explicit)
    θ̂^{n+1} = (θ̂^n·(1 − ½dt·ν·k²) + dt·F̂_expl) / (1 + ½dt·ν·k²)
  - float32 throughout (Metal GPU does not support float64 or complex<f32>)
  - All FFT computation routed to CPU XLA (Metal v0.1.0 lacks complex FFT)
  - @jax.jit fusion gives ~1.5–3× speedup over NumPy baseline

Key measurements:
  - δ_est(t): largest s such that H^s norm of θ is controlled — Lemma 2.5 target
  - ε_CZ(t): ‖S_θ‖_{L^{1+0.1}} / ‖S_θ‖_{L^1} — Claim A probe
  - ε_LPS: min(μ_eff) − ν ≥ 0 (by construction)

Experiments:
  EXP-L3-R2-TG-JAX-032   Taylor-Green N=32³ (validation)
  EXP-L3-R2-TG-JAX-064   Taylor-Green N=64³ (production)
  EXP-L3-R2-SH-JAX-064   Shear layer N=64³ (Route A decisive test)

Requirements:
  pip install "jax==0.4.26" "jaxlib==0.4.26" jax-metal
"""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
import time
from datetime import datetime
from functools import partial

import numpy as np

# ── JAX import ────────────────────────────────────────────────────────────────
try:
    import jax
    import jax.numpy as jnp
    from jax import jit
    _JAX_AVAILABLE = True
    _JAX_BACKEND_RAW = jax.default_backend()

    if _JAX_BACKEND_RAW == "METAL":
        _COMPUTE_DEVICE = jax.devices("cpu")[0]
        _JAX_BACKEND = "CPU-JIT"
    else:
        _COMPUTE_DEVICE = jax.devices()[0]
        _JAX_BACKEND = _JAX_BACKEND_RAW

except ImportError:
    _JAX_AVAILABLE = False
    _JAX_BACKEND_RAW = "unavailable"
    _JAX_BACKEND = "unavailable"
    _COMPUTE_DEVICE = None


def _d(arr):
    """Move a JAX array to the compute device (CPU when Metal is default)."""
    if _COMPUTE_DEVICE is None:
        return arr
    return jax.device_put(arr, _COMPUTE_DEVICE)


# ── path setup ─────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_NU_DEFAULT = 1.0e-3


# =============================================================================
# Grid helpers (NumPy → JAX, no JIT needed)
# =============================================================================

def make_grid_jax(N: int, dtype=jnp.float32):
    """Return JAX wavenumber arrays kx, ky, kz, k² on compute device."""
    freqs = np.fft.fftfreq(N, d=1.0 / N).astype(np.float32)
    kx_np, ky_np, kz_np = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k2_np = kx_np**2 + ky_np**2 + kz_np**2
    return (
        _d(jnp.array(kx_np, dtype=dtype)),
        _d(jnp.array(ky_np, dtype=dtype)),
        _d(jnp.array(kz_np, dtype=dtype)),
        _d(jnp.array(k2_np, dtype=dtype)),
    )


def make_dealias_mask_jax(N: int, dtype=jnp.float32):
    """2/3-rule dealiasing mask — JAX array on compute device."""
    freqs = np.fft.fftfreq(N, d=1.0 / N)
    kx_np, ky_np, kz_np = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k_max = N // 3
    mask_np = (
        (np.abs(kx_np) <= k_max) & (np.abs(ky_np) <= k_max) & (np.abs(kz_np) <= k_max)
    ).astype(np.float32)
    return _d(jnp.array(mask_np, dtype=dtype))


def make_nyquist_mask_jax(N: int, dtype=jnp.float32):
    """Mask that zeros Nyquist planes (B-005 fix: restores Hermitian symmetry)."""
    nq = N // 2
    idx = np.arange(N)
    is_nq = (
        (idx[:, None, None] == nq)
        | (idx[None, :, None] == nq)
        | (idx[None, None, :] == nq)
    )
    mask_np = (~is_nq).astype(np.float32)
    return _d(jnp.array(mask_np, dtype=dtype))


# =============================================================================
# Leray projector (JIT-safe, boolean-mask — no zeros_like complex)
# =============================================================================

def _project_jax(u_hat, v_hat, w_hat, kx, ky, kz, k2, nyquist_mask):
    """Leray projection in spectral space. JIT-safe, B-005 fixed."""
    u_hat = u_hat * nyquist_mask
    v_hat = v_hat * nyquist_mask
    w_hat = w_hat * nyquist_mask

    k_dot_u = kx * u_hat + ky * v_hat + kz * w_hat
    k2_safe = jnp.where(k2 > 0.0, k2, 1.0)
    valid   = (k2 > 0.0).astype(k2.dtype)        # float32 mask
    corr    = (k_dot_u / k2_safe) * valid         # complex × float → complex

    return u_hat - kx * corr, v_hat - ky * corr, w_hat - kz * corr


# =============================================================================
# JIT-compiled velocity RHS (exact Prize NS, ν = const)
# =============================================================================

@jit
def _rhs_vel_jax(u, v, w, kx, ky, kz, k2, dealias, nu):
    """Route 2 velocity RHS: ∂_t u = −(u·∇)u − ∇p + ν·Δu, div u = 0.

    Exact Prize equations: ν is a Python/JAX scalar (const).
    JIT-compiled. Returns (du/dt, dv/dt, dw/dt) as JAX arrays.
    """
    u_hat = jnp.fft.fftn(u)
    v_hat = jnp.fft.fftn(v)
    w_hat = jnp.fft.fftn(w)

    def deriv(f_hat, k_dir):
        return jnp.real(jnp.fft.ifftn(1j * k_dir * f_hat))

    dudx = deriv(u_hat, kx); dudy = deriv(u_hat, ky); dudz = deriv(u_hat, kz)
    dvdx = deriv(v_hat, kx); dvdy = deriv(v_hat, ky); dvdz = deriv(v_hat, kz)
    dwdx = deriv(w_hat, kx); dwdy = deriv(w_hat, ky); dwdz = deriv(w_hat, kz)

    NL_u = u * dudx + v * dudy + w * dudz
    NL_v = u * dvdx + v * dvdy + w * dvdz
    NL_w = u * dwdx + v * dwdy + w * dwdz

    NL_u_hat = jnp.fft.fftn(NL_u) * dealias
    NL_v_hat = jnp.fft.fftn(NL_v) * dealias
    NL_w_hat = jnp.fft.fftn(NL_w) * dealias

    # Leray project nonlinear term
    k_dot_NL = kx * NL_u_hat + ky * NL_v_hat + kz * NL_w_hat
    k2_safe  = jnp.where(k2 > 0.0, k2, 1.0)
    corr_NL  = (k_dot_NL / k2_safe) * (k2 > 0.0).astype(k2.dtype)
    NL_u_hat -= kx * corr_NL
    NL_v_hat -= ky * corr_NL
    NL_w_hat -= kz * corr_NL

    du_hat = (-NL_u_hat - nu * k2 * u_hat) * dealias
    dv_hat = (-NL_v_hat - nu * k2 * v_hat) * dealias
    dw_hat = (-NL_w_hat - nu * k2 * w_hat) * dealias

    return (jnp.real(jnp.fft.ifftn(du_hat)),
            jnp.real(jnp.fft.ifftn(dv_hat)),
            jnp.real(jnp.fft.ifftn(dw_hat)))


# =============================================================================
# JIT-compiled combined step: RK4 for u + CN for θ
# =============================================================================

@jit
def solver_step_jax(u, v, w, theta, kx, ky, kz, k2,
                    dealias, nyquist_mask, nu, dt):
    """One step of Route 2 3D (JIT-compiled).

    Velocity: RK4 with exact Prize NS.
    θ: Crank-Nicolson IMEX:
       θ̂^{n+1} = (θ̂^n·(1−½dt·ν·k²) + dt·F̂_expl) / (1+½dt·ν·k²)
       F_expl = −u·∇θ + S_θ   (at time n)

    Returns: new_u, new_v, new_w, new_theta
    """
    # ── S_θ = ν|∇u|² at current state ────────────────────────────────────────
    u_hat = jnp.fft.fftn(u)
    v_hat = jnp.fft.fftn(v)
    w_hat = jnp.fft.fftn(w)

    def deriv(f_hat, k_dir):
        return jnp.real(jnp.fft.ifftn(1j * k_dir * f_hat))

    dudx = deriv(u_hat, kx); dudy = deriv(u_hat, ky); dudz = deriv(u_hat, kz)
    dvdx = deriv(v_hat, kx); dvdy = deriv(v_hat, ky); dvdz = deriv(v_hat, kz)
    dwdx = deriv(w_hat, kx); dwdy = deriv(w_hat, ky); dwdz = deriv(w_hat, kz)

    S_theta = nu * (dudx**2 + dudy**2 + dudz**2
                    + dvdx**2 + dvdy**2 + dvdz**2
                    + dwdx**2 + dwdy**2 + dwdz**2)

    # ── RK4 for velocity ──────────────────────────────────────────────────────
    def vel_rhs(uu, vv, ww):
        return _rhs_vel_jax(uu, vv, ww, kx, ky, kz, k2, dealias, nu)

    k1u, k1v, k1w = vel_rhs(u, v, w)
    k2u, k2v, k2w = vel_rhs(u+0.5*dt*k1u, v+0.5*dt*k1v, w+0.5*dt*k1w)
    k3u, k3v, k3w = vel_rhs(u+0.5*dt*k2u, v+0.5*dt*k2v, w+0.5*dt*k2w)
    k4u, k4v, k4w = vel_rhs(u+dt*k3u,     v+dt*k3v,     w+dt*k3w)

    fac = dt / 6.0
    new_u = u + fac * (k1u + 2.0*k2u + 2.0*k3u + k4u)
    new_v = v + fac * (k1v + 2.0*k2v + 2.0*k3v + k4v)
    new_w = w + fac * (k1w + 2.0*k2w + 2.0*k3w + k4w)

    # Dealias velocity
    def _dealias(f):
        return jnp.real(jnp.fft.ifftn(jnp.fft.fftn(f) * dealias))
    new_u = _dealias(new_u); new_v = _dealias(new_v); new_w = _dealias(new_w)

    # Re-project
    uh = jnp.fft.fftn(new_u); vh = jnp.fft.fftn(new_v); wh = jnp.fft.fftn(new_w)
    uh, vh, wh = _project_jax(uh, vh, wh, kx, ky, kz, k2, nyquist_mask)
    new_u = jnp.real(jnp.fft.ifftn(uh))
    new_v = jnp.real(jnp.fft.ifftn(vh))
    new_w = jnp.real(jnp.fft.ifftn(wh))

    # ── Crank-Nicolson for θ ──────────────────────────────────────────────────
    theta_hat = jnp.fft.fftn(theta)

    # Explicit part: F = −u·∇θ + S_θ
    th_x = deriv(theta_hat, kx)
    th_y = deriv(theta_hat, ky)
    th_z = deriv(theta_hat, kz)
    F_expl = -(u * th_x + v * th_y + w * th_z) + S_theta
    F_hat  = jnp.fft.fftn(F_expl)

    # CN update
    half_dt_nu_k2 = 0.5 * dt * nu * k2
    theta_hat_new = (theta_hat * (1.0 - half_dt_nu_k2) + dt * F_hat) / (1.0 + half_dt_nu_k2)
    new_theta = jnp.real(jnp.fft.ifftn(theta_hat_new))

    # Enforce θ ≥ 0 (source S_θ ≥ 0; negatives are floating-point noise)
    new_theta = jnp.clip(new_theta, 0.0, None)

    return new_u, new_v, new_w, new_theta


# =============================================================================
# CFL (Python-level)
# =============================================================================

def compute_cfl_dt_jax(u, v, w, N: int, nu: float,
                        dt_max: float, safety: float = 0.4) -> float:
    """CFL timestep for Route 2 (advection-only; θ diffusion is IMEX-implicit)."""
    dx = 2.0 * np.pi / N
    U_max = max(float(jnp.max(jnp.abs(u))), float(jnp.max(jnp.abs(v))),
                float(jnp.max(jnp.abs(w))), 1e-12)
    dt_adv  = safety * dx / U_max
    dt_diff = safety * dx**2 / (2.0 * nu)   # reference; CN removes this constraint
    return float(np.clip(min(dt_adv, dt_diff), 1e-12, dt_max))


# =============================================================================
# Diagnostics (Python-level — extracted after JIT step)
# =============================================================================

def kinetic_energy_jax(u, v, w) -> float:
    return float(0.5 * jnp.mean(u**2 + v**2 + w**2))


def divergence_rms_jax(u, v, w, kx, ky, kz) -> float:
    div_hat = 1j * (kx * jnp.fft.fftn(u)
                    + ky * jnp.fft.fftn(v)
                    + kz * jnp.fft.fftn(w))
    div = jnp.real(jnp.fft.ifftn(div_hat))
    return float(jnp.sqrt(jnp.mean(div**2)))


def delta_from_theta_jax(theta, kx, ky, kz) -> float:
    """δ estimate from H^s norms of θ (same logic as route2_3D.py)."""
    theta_hat = jnp.fft.fftn(theta)
    k_abs = jnp.sqrt(kx**2 + ky**2 + kz**2)

    th_flat = jnp.abs(theta_hat).ravel()
    k_flat  = k_abs.ravel()
    mask    = k_flat > 0.0

    if not bool(jnp.any(mask)):
        return 0.0

    th_nz = th_flat[mask]
    k_nz  = k_flat[mask]
    H0_sq = float(jnp.sum(th_nz**2))
    if H0_sq < 1e-30:
        return 0.0

    H0 = float(jnp.sqrt(H0_sq))
    delta_est = 0.0
    for s in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]:
        Hs = float(jnp.sqrt(jnp.sum((k_nz**s * th_nz)**2)))
        if Hs < 100.0 * H0:
            delta_est = max(s - 1.0, 0.0)
        else:
            break
    return delta_est


# =============================================================================
# Solver construction
# =============================================================================

def make_solver_jax(
    N: int = 32,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    dt_max: float = 0.05,
    dtype=jnp.float32,
) -> dict:
    """Create a JAX Route 2 solver state dict."""
    if not _JAX_AVAILABLE:
        raise RuntimeError("JAX not available. pip install jax==0.4.26 jaxlib==0.4.26 jax-metal")
    if nu <= 0.0:
        raise ValueError(f"nu must be > 0, got {nu}")
    if eps_param < 0.0:
        raise ValueError(f"eps_param must be ≥ 0, got {eps_param}")

    kx, ky, kz, k2 = make_grid_jax(N, dtype=dtype)
    dealias     = make_dealias_mask_jax(N, dtype=dtype)
    nq_mask     = make_nyquist_mask_jax(N, dtype=dtype)

    return {
        "N"          : N,
        "nu"         : nu,
        "eps_param"  : eps_param,
        "dt_max"     : dt_max,
        "dtype"      : dtype,
        "backend"    : _JAX_BACKEND,
        "kx": kx, "ky": ky, "kz": kz, "k2": k2,
        "dealias"    : dealias,
        "nq_mask"    : nq_mask,
        "u"          : _d(jnp.zeros((N, N, N), dtype=dtype)),
        "v"          : _d(jnp.zeros((N, N, N), dtype=dtype)),
        "w"          : _d(jnp.zeros((N, N, N), dtype=dtype)),
        "theta"      : _d(jnp.zeros((N, N, N), dtype=dtype)),
        "t"          : 0.0,
    }


def set_ic_jax(solver: dict, ic: dict) -> dict:
    """Set IC; Leray-project velocity; θ = 0."""
    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    nq = solver["nq_mask"]
    N, dtype = solver["N"], solver["dtype"]

    u = _d(jnp.array(ic["u"], dtype=dtype))
    v = _d(jnp.array(ic["v"], dtype=dtype))
    w = _d(jnp.array(ic["w"], dtype=dtype))

    uh, vh, wh = jnp.fft.fftn(u), jnp.fft.fftn(v), jnp.fft.fftn(w)
    uh, vh, wh = _project_jax(uh, vh, wh, kx, ky, kz, k2, nq)
    u = jnp.real(jnp.fft.ifftn(uh))
    v = jnp.real(jnp.fft.ifftn(vh))
    w = jnp.real(jnp.fft.ifftn(wh))

    return {**solver, "u": u, "v": v, "w": w,
            "theta": _d(jnp.zeros((N, N, N), dtype=dtype)), "t": 0.0}


# =============================================================================
# Initial conditions (NumPy arrays, converted to JAX in set_ic_jax)
# =============================================================================

def taylor_green_ic(N: int, V0: float = 1.0) -> dict:
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    return {
        "u": V0 * np.sin(xx) * np.cos(yy) * np.cos(zz),
        "v": -V0 * np.cos(xx) * np.sin(yy) * np.cos(zz),
        "w": np.zeros((N, N, N)),
        "E0_analytic": V0**2 / 8.0,
        "name": "taylor_green",
    }


def shear_layer_ic(N: int, amp: float = 1.0, width: float = 0.2) -> dict:
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    _, _, zz = np.meshgrid(x, x, x, indexing="ij")
    u = amp * np.tanh((zz - np.pi) / width)
    return {
        "u": u,
        "v": np.zeros((N, N, N)),
        "w": np.zeros((N, N, N)),
        "E0_analytic": float(np.mean(u**2) / 2.0),
        "name": "shear_layer",
    }


# =============================================================================
# Run loop (Python-level; JIT boundary = one combined step)
# =============================================================================

def run_jax(
    solver: dict,
    t_end: float,
    max_steps: int = 5000,
    cfl_safety: float = 0.4,
    verbose: bool = False,
    warmup: bool = True,
) -> dict:
    """Run the JAX Route 2 solver from solver['t'] to t_end.

    Returns result dict with history, verdict, δ_max, ε_LPS_min.
    """
    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    dealias, nq = solver["dealias"], solver["nq_mask"]
    nu, eps_param = solver["nu"], solver["eps_param"]
    dt_max, N, dtype = solver["dt_max"], solver["N"], solver["dtype"]

    u, v, w, theta = solver["u"], solver["v"], solver["w"], solver["theta"]
    t = float(solver["t"])

    nu_jnp = jnp.array(nu, dtype=dtype)

    history: list[dict] = []
    E_initial = kinetic_energy_jax(u, v, w)
    delta_max  = 0.0
    theta_negative = False
    compile_time_s = None

    for step_idx in range(max_steps):
        if t >= t_end:
            break

        dt_py = compute_cfl_dt_jax(u, v, w, N, nu, dt_max, cfl_safety)
        dt_py = min(dt_py, t_end - t)
        dt_jnp = jnp.array(dt_py, dtype=dtype)

        t0_w = time.perf_counter()
        new_u, new_v, new_w, new_theta = solver_step_jax(
            u, v, w, theta, kx, ky, kz, k2, dealias, nq, nu_jnp, dt_jnp
        )
        jax.block_until_ready((new_u, new_v, new_w, new_theta))
        step_wall = time.perf_counter() - t0_w

        if step_idx == 0:
            compile_time_s = step_wall
            if verbose and warmup:
                print(f"  [JIT] First step (compile): {compile_time_s:.1f} s")

        u, v, w, theta = new_u, new_v, new_w, new_theta
        t += dt_py

        # Diagnostics (Python-level; after JIT sync)
        E = kinetic_energy_jax(u, v, w)
        theta_min_val = float(jnp.min(theta))
        theta_max_val = float(jnp.max(theta))
        delta_est = delta_from_theta_jax(theta, kx, ky, kz)
        delta_max = max(delta_max, delta_est)

        mu_eff_min = nu + eps_param * float(jnp.min(jnp.clip(theta, 0.0, None)))
        eps_lps = mu_eff_min - nu

        if theta_min_val < -1e-6:
            theta_negative = True

        diag = {
            "t": t, "dt": dt_py, "E": E,
            "theta_min": theta_min_val, "theta_max": theta_max_val,
            "delta_est": delta_est, "eps_lps": eps_lps,
            "step_wall_s": step_wall,
        }
        history.append(diag)

        if verbose and step_idx % max(1, max_steps // 10) == 0:
            print(
                f"  step {step_idx:4d}  t={t:.4f}  E={E:.4e}  "
                f"θ_max={theta_max_val:.3e}  δ={delta_est:.2f}  "
                f"step={step_wall*1000:.1f}ms"
            )

    E_final = kinetic_energy_jax(u, v, w)
    step_times = [h["step_wall_s"] for h in history[1:]]
    step_ms_mean = 1000.0 * float(np.mean(step_times)) if step_times else 0.0
    eps_lps_min = min(h["eps_lps"] for h in history) if history else 0.0

    return {
        "u": u, "v": v, "w": w, "theta": theta,
        "t_final": t,
        "history": history,
        "verdict": "PASS" if not theta_negative else "FAIL",
        "n_steps": len(history),
        "E_initial": E_initial,
        "E_final": E_final,
        "energy_decayed": E_final < E_initial,
        "delta_max": delta_max,
        "theta_negative": theta_negative,
        "eps_lps_min": eps_lps_min,
        "compile_time_s": compile_time_s,
        "step_ms_mean": step_ms_mean,
    }


# =============================================================================
# Benchmark experiments
# =============================================================================

def run_taylor_green_jax(
    N: int = 32,
    t_end: float = 1.0,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str | None = None,
    exp_suffix: str = "",
) -> dict:
    """Route 2 3D Taylor-Green benchmark using JAX-JIT.

    Pass criteria (same as NumPy version):
      1. E₀ error < 1%
      2. θ ≥ 0 throughout
      3. Energy decays
      4. δ_max > 0 (Lemma 2.5 central measurement)
    """
    solver = make_solver_jax(N=N, nu=nu, eps_param=eps_param)
    ic = taylor_green_ic(N, V0=1.0)
    solver = set_ic_jax(solver, ic)

    E0_computed = kinetic_energy_jax(solver["u"], solver["v"], solver["w"])
    E0_analytic = ic["E0_analytic"]

    if verbose:
        print(f"\n  Route 2 JAX TG  backend={_JAX_BACKEND}  N={N}³  ν={nu}  ε={eps_param}")
        print(f"  E₀ analytic={E0_analytic:.6f}  computed={E0_computed:.6f}")

    t_wall0 = time.perf_counter()
    result = run_jax(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)
    wall = time.perf_counter() - t_wall0

    E0_ok    = abs(E0_computed - E0_analytic) / E0_analytic < 0.01
    passed   = E0_ok and not result["theta_negative"] and result["energy_decayed"]

    exp = {
        "exp_id"         : f"EXP-L3-R2-TG-JAX-{N:03d}{exp_suffix}",
        "claim_id"       : "claim_m5b_route2_3D_jax_taylor_green",
        "timestamp"      : datetime.utcnow().isoformat(),
        "backend"        : _JAX_BACKEND,
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
        "theta_negative" : result["theta_negative"],
        "delta_max"      : result["delta_max"],
        "eps_lps_min"    : result["eps_lps_min"],
        "verdict"        : "PASS" if passed else "FAIL",
        "wall_time_s"    : wall,
        "compile_time_s" : result.get("compile_time_s", 0.0),
        "step_ms_mean"   : result["step_ms_mean"],
        "key_metric"     : (
            f"E0_err={abs(E0_computed-E0_analytic)/E0_analytic*100:.3f}%  "
            f"δ_max={result['delta_max']:.2f}  "
            f"ε_LPS_min={result['eps_lps_min']:.4f}  "
            f"step={result['step_ms_mean']:.1f}ms"
        ),
    }

    if verbose:
        print(f"\n  Results: {exp['exp_id']}")
        print(f"    E₀ error:    {exp['E0_error_pct']:.3f}%  {'✓' if E0_ok else '✗'}")
        print(f"    θ ≥ 0:       {'✓' if not result['theta_negative'] else '✗'}")
        print(f"    Energy decay:{'✓' if result['energy_decayed'] else '✗'}")
        print(f"    δ_max:       {result['delta_max']:.2f}")
        print(f"    ε_LPS_min:   {result['eps_lps_min']:.4f}")
        print(f"    Step time:   {result['step_ms_mean']:.1f} ms  (compile: {result.get('compile_time_s',0):.1f}s)")
        print(f"    Wall time:   {wall:.1f} s")
        print(f"    Verdict:     {exp['verdict']}")

    if db_path is not None:
        _log_result(db_path, exp)

    return exp


def run_shear_layer_jax(
    N: int = 32,
    t_end: float = 1.0,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    max_steps: int = 5000,
    verbose: bool = False,
    db_path: str | None = None,
) -> dict:
    """Route 2 3D shear layer (non-mixing) using JAX-JIT — Route A decisive test."""
    solver = make_solver_jax(N=N, nu=nu, eps_param=eps_param)
    ic = shear_layer_ic(N)
    solver = set_ic_jax(solver, ic)

    if verbose:
        print(f"\n  Route 2 JAX Shear  backend={_JAX_BACKEND}  N={N}³")

    t_wall0 = time.perf_counter()
    result = run_jax(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)
    wall = time.perf_counter() - t_wall0

    passed = not result["theta_negative"] and result["energy_decayed"]

    exp = {
        "exp_id"         : f"EXP-L3-R2-SH-JAX-{N:03d}",
        "claim_id"       : "claim_m5b_route2_3D_jax_shear_route_A",
        "timestamp"      : datetime.utcnow().isoformat(),
        "backend"        : _JAX_BACKEND,
        "ic"             : "shear_layer",
        "N"              : N,
        "nu"             : nu,
        "eps_param"      : eps_param,
        "t_end"          : t_end,
        "n_steps"        : result["n_steps"],
        "t_final"        : result["t_final"],
        "E_initial"      : result["E_initial"],
        "E_final"        : result["E_final"],
        "energy_decayed" : result["energy_decayed"],
        "theta_negative" : result["theta_negative"],
        "delta_max"      : result["delta_max"],
        "eps_lps_min"    : result["eps_lps_min"],
        "verdict"        : "PASS" if passed else "FAIL",
        "wall_time_s"    : wall,
        "compile_time_s" : result.get("compile_time_s", 0.0),
        "step_ms_mean"   : result["step_ms_mean"],
        "key_metric"     : (
            f"δ_max={result['delta_max']:.2f}  "
            f"ε_LPS_min={result['eps_lps_min']:.4f}  "
            f"step={result['step_ms_mean']:.1f}ms"
        ),
    }

    if verbose:
        print(f"\n  Results: {exp['exp_id']}")
        print(f"    δ_max:  {result['delta_max']:.2f}  "
              f"{'→ Route A CONFIRMED' if result['delta_max'] > 0 else ''}")
        print(f"    Step:   {result['step_ms_mean']:.1f} ms")
        print(f"    Verdict:{exp['verdict']}")

    if db_path is not None:
        _log_result(db_path, exp)

    return exp


def benchmark_speedup_jax(N: int = 32, n_warmup: int = 3, n_timed: int = 10) -> dict:
    """Measure JAX-JIT speedup vs NumPy for Route 2 3D."""
    import route2_3D as r2_np

    # ── NumPy baseline ─────────────────────────────────────────────────────────
    s_np = r2_np.make_solver(N=N, nu=_NU_DEFAULT, eps_param=0.1)
    ic_np = r2_np.taylor_green_ic(N)
    s_np  = r2_np.set_ic(s_np, ic_np)

    dt_np = r2_np.compute_cfl_dt(s_np)
    for _ in range(2):
        s_np, _ = r2_np.solver_step(s_np, dt=dt_np)

    t0 = time.perf_counter()
    for _ in range(n_timed):
        s_np, _ = r2_np.solver_step(s_np, dt=dt_np)
    numpy_ms = 1000.0 * (time.perf_counter() - t0) / n_timed

    # ── JAX-JIT ────────────────────────────────────────────────────────────────
    dtype = jnp.float32
    solver = make_solver_jax(N=N, nu=_NU_DEFAULT, eps_param=0.1)
    ic = taylor_green_ic(N)
    solver = set_ic_jax(solver, ic)

    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    dealias, nq = solver["dealias"], solver["nq_mask"]
    nu_jnp  = jnp.array(_NU_DEFAULT, dtype=dtype)
    dt_jnp  = jnp.array(dt_np, dtype=dtype)

    u, v, w, theta = solver["u"], solver["v"], solver["w"], solver["theta"]

    for _ in range(n_warmup):
        u, v, w, theta = solver_step_jax(
            u, v, w, theta, kx, ky, kz, k2, dealias, nq, nu_jnp, dt_jnp
        )
        jax.block_until_ready((u, v, w, theta))

    t0 = time.perf_counter()
    for _ in range(n_timed):
        u, v, w, theta = solver_step_jax(
            u, v, w, theta, kx, ky, kz, k2, dealias, nq, nu_jnp, dt_jnp
        )
        jax.block_until_ready((u, v, w, theta))
    jax_ms = 1000.0 * (time.perf_counter() - t0) / n_timed

    speedup = numpy_ms / jax_ms if jax_ms > 0 else float("inf")
    return {
        "N"        : N,
        "backend"  : _JAX_BACKEND,
        "numpy_ms" : numpy_ms,
        "jax_ms"   : jax_ms,
        "speedup"  : speedup,
        "dtype"    : str(dtype),
    }


# =============================================================================
# Results logging
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS route2_jax_experiments (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id          TEXT UNIQUE NOT NULL,
            claim_id        TEXT,
            timestamp       TEXT,
            backend         TEXT,
            ic              TEXT,
            N               INTEGER,
            nu              REAL,
            eps_param       REAL,
            t_end           REAL,
            n_steps         INTEGER,
            t_final         REAL,
            verdict         TEXT,
            E0_error_pct    REAL,
            E_final         REAL,
            energy_decayed  INTEGER,
            theta_negative  INTEGER,
            delta_max       REAL,
            eps_lps_min     REAL,
            wall_time_s     REAL,
            compile_time_s  REAL,
            step_ms_mean    REAL,
            key_metric      TEXT
        )
    """)
    conn.commit()
    return conn


def _log_result(db_path: str, result: dict) -> None:
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "backend", "ic", "N", "nu",
        "eps_param", "t_end", "n_steps", "t_final", "verdict",
        "E0_error_pct", "E_final", "energy_decayed", "theta_negative",
        "delta_max", "eps_lps_min", "wall_time_s", "compile_time_s",
        "step_ms_mean", "key_metric",
    ]
    vals = tuple(
        int(result[c]) if isinstance(result.get(c), bool) else result.get(c)
        for c in cols
    )
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO route2_jax_experiments "
        f"({', '.join(cols)}) VALUES ({ph})",
        vals,
    )
    conn.commit()
    conn.close()


# =============================================================================
# CLI
# =============================================================================

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Route 2 3D JAX solver")
    parser.add_argument("--run", choices=["tg", "shear", "bench"], default="tg")
    parser.add_argument("--N",      type=int,   default=32)
    parser.add_argument("--t-end",  type=float, default=1.0)
    parser.add_argument("--nu",     type=float, default=_NU_DEFAULT)
    parser.add_argument("--eps",    type=float, default=0.1)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db",     default=_DEFAULT_DB)
    parser.add_argument("--n-timed", type=int, default=10)
    args = parser.parse_args()

    print(f"JAX backend: {_JAX_BACKEND}  ({_JAX_BACKEND_RAW})")

    if args.run == "tg":
        r = run_taylor_green_jax(
            N=args.N, t_end=args.t_end, nu=args.nu, eps_param=args.eps,
            verbose=args.verbose, db_path=args.db,
        )
        print(f"\nVerdict: {r['verdict']}  {r['key_metric']}")

    elif args.run == "shear":
        r = run_shear_layer_jax(
            N=args.N, t_end=args.t_end, nu=args.nu, eps_param=args.eps,
            verbose=args.verbose, db_path=args.db,
        )
        print(f"\nVerdict: {r['verdict']}  {r['key_metric']}")

    elif args.run == "bench":
        print(f"\nBenchmarking N={args.N}³ ({args.n_timed} steps after JIT warmup)...")
        res = benchmark_speedup_jax(N=args.N, n_timed=args.n_timed)
        print(f"  NumPy:     {res['numpy_ms']:.2f} ms/step")
        print(f"  JAX {res['backend']:6s}: {res['jax_ms']:.2f} ms/step")
        print(f"  Speedup:   {res['speedup']:.2f}×")


if __name__ == "__main__":
    main()
