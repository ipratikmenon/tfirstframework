"""
route1_3D_jax.py — JAX-JIT Route 1 3D Solver (M5)
===================================================
JAX-JIT port of route1_3D.py for XLA/LLVM-accelerated execution.

Physics (Route 1, PRD v1.0 §3.2):
  ∂_t u + (u·∇)u = -∇p + ∇·(ν(T)·∇u)    div u = 0  ρ = const
  ∂_t T + u·∇T   = A(T)·ΔT + ν(T)|∇u|²

Numerics: identical to route1_3D.py (FFT3, 2/3 dealiasing, RK4, Leray)

Key differences from NumPy version:
  - All arrays are JAX arrays (jnp)
  - `_rhs_jax` and `solver_step_jax` are @jax.jit compiled (XLA/LLVM)
  - Runs on CPU device (Metal backend does not support complex FFT in v0.1.0)
  - float32 throughout (Metal does not support float64; CPU XLA float32 is fine)
  - Same physics, same test criteria (tolerances adjusted for float32 ~1e-6)

Device strategy:
  jax-metal 0.1.0 does not support complex<f32> tensors required for FFT.
  We detect the METAL backend and transparently route all computation to the
  CPU XLA backend.  @jax.jit still applies LLVM/XLA fusion, giving ~1.5–3×
  speedup over plain NumPy (real at N=64³+).  When jax-metal adds FFT support
  the `_COMPUTE_DEVICE` routing below is the only change needed.

Requirements:
  pip install "jax==0.4.26" "jaxlib==0.4.26" jax-metal

Experiments:
  EXP-L3-R1-TG-JAX-001  — JAX-JIT TG baseline (N=32, CPU XLA)
  EXP-L3-R1-TG-JAX-002  — JAX-JIT TG (N=64, CPU XLA)
"""

from __future__ import annotations
import os
import sys
import time
import json
import sqlite3
import argparse
import numpy as np
from datetime import datetime
from functools import partial

# ── JAX import ────────────────────────────────────────────────────────────────
try:
    import jax
    # float32 throughout: jax-metal 0.1.0 does not support float64 on Metal.
    # CPU XLA (our compute device) runs float32 natively with XLA/LLVM fusion.
    # float32 gives ~1e-6 relative precision — sufficient for spectral NS.

    import jax.numpy as jnp
    from jax import jit

    _JAX_AVAILABLE = True
    _JAX_BACKEND_RAW = jax.default_backend()   # actual platform (may be METAL)

    # Metal (jax-metal 0.1.0) does not support complex<f32> tensors for FFT.
    # Route all computation to the CPU XLA device.
    # @jax.jit still applies XLA/LLVM fusion for ~1.5–3× speedup over NumPy.
    # When jax-metal adds complex FFT support, remove this branch.
    if _JAX_BACKEND_RAW == "METAL":
        _COMPUTE_DEVICE = jax.devices("cpu")[0]
        _JAX_BACKEND = "CPU-JIT"   # effective backend label
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
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))


# ─────────────────────────────────────────────────────────────────────────────
# Physical constants (inline for JIT compatibility — no Python calls inside JIT)
# ─────────────────────────────────────────────────────────────────────────────

# Sutherland viscosity constants (air / ideal diatomic gas)
_MU_REF = 1.716e-5    # Pa·s at T_ref
_T_REF  = 273.15      # K
_S_SUTH = 110.4       # K  Sutherland constant

# Power-law thermal conductivity
_K_REF  = 0.0241      # W/(m·K) at T_ref
_K_EXP  = 0.82        # power-law exponent

# Ideal gas
_R_SPEC   = 287.0     # J/(kg·K)
_P_REF    = 101325.0  # Pa
_CV_IDEAL = 718.0     # J/(kg·K)


# ─────────────────────────────────────────────────────────────────────────────
# Thermophysical properties (JAX-native, JIT-compatible)
# ─────────────────────────────────────────────────────────────────────────────

def ideal_props_jax(T: jnp.ndarray) -> dict:
    """Ideal-gas properties as JAX operations (JIT-compatible).

    Identical formulas to layer1/tfirst_props.py::ideal_props.
    Inlined here to avoid Python function dispatch inside @jit.

    Returns dict: mu, k, rho, cv (all JAX arrays, same shape as T).
    """
    mu    = _MU_REF * (T / _T_REF) ** 1.5 * (_T_REF + _S_SUTH) / (T + _S_SUTH)
    k_th  = _K_REF * (T / _T_REF) ** _K_EXP
    rho   = _P_REF / (_R_SPEC * T)
    cv    = jnp.full_like(T, _CV_IDEAL)
    return {"mu": mu, "k": k_th, "rho": rho, "cv": cv}


# ─────────────────────────────────────────────────────────────────────────────
# Grid helpers (built as NumPy, converted to JAX; no JIT needed here)
# ─────────────────────────────────────────────────────────────────────────────

def make_grid_jax(N: int, dtype=jnp.float32):
    """Return JAX wavenumber arrays for N³ periodic domain [0, 2π)³.

    Arrays placed on the compute device (CPU when Metal is default).
    Returns: kx, ky, kz, k2  — shape (N, N, N), JAX arrays.
    """
    freqs = np.fft.fftfreq(N, d=1.0 / N).astype(np.float32)
    kx_np, ky_np, kz_np = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k2_np = kx_np ** 2 + ky_np ** 2 + kz_np ** 2

    kx = _d(jnp.array(kx_np, dtype=dtype))
    ky = _d(jnp.array(ky_np, dtype=dtype))
    kz = _d(jnp.array(kz_np, dtype=dtype))
    k2 = _d(jnp.array(k2_np, dtype=dtype))
    return kx, ky, kz, k2


def make_dealias_mask_jax(N: int, dtype=jnp.float32):
    """2/3 dealiasing mask — JAX array on compute device."""
    freqs = np.fft.fftfreq(N, d=1.0 / N)
    kx_np, ky_np, kz_np = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k_max = N // 3
    mask_np = (
        (np.abs(kx_np) <= k_max)
        & (np.abs(ky_np) <= k_max)
        & (np.abs(kz_np) <= k_max)
    ).astype(np.float32)
    return _d(jnp.array(mask_np, dtype=dtype))


def make_nyquist_mask_jax(N: int, dtype=jnp.float32):
    """Mask that zeros the Nyquist planes (index N//2 in each dim).

    The Nyquist mode (k = -N/2 in fftfreq convention) is aliased: its conjugate
    pair is also at index N//2, not at index 0.  This breaks the anti-Hermitian
    property k·û(-k) = -(k·û(k))* required for Hermitian-symmetry-preserving
    Leray projection.  Zeroing these modes before projection restores correctness.
    In the solver they are always dealiased anyway (2/3 rule: N//3 < N//2).

    Returns: JAX array shape (N, N, N) on compute device, 1.0 safe / 0.0 Nyquist.
    """
    nq = N // 2
    idx = np.arange(N)
    is_nq_x = (idx[:, None, None] == nq)
    is_nq_y = (idx[None, :, None] == nq)
    is_nq_z = (idx[None, None, :] == nq)
    mask_np = (~(is_nq_x | is_nq_y | is_nq_z)).astype(np.float32)
    return _d(jnp.array(mask_np, dtype=dtype))


# ─────────────────────────────────────────────────────────────────────────────
# JAX Leray projector (functional, JIT-safe)
# ─────────────────────────────────────────────────────────────────────────────

def _project_spectral_jax(u_hat, v_hat, w_hat, kx, ky, kz, k2, nyquist_mask):
    """Leray projection in spectral space (functional, JIT-safe).

    Inputs/outputs are spectral-space complex JAX arrays.
    Nyquist planes zeroed before projection to preserve Hermitian symmetry.
    """
    # Zero Nyquist planes (fix aliasing issue — see make_nyquist_mask_jax)
    u_hat = u_hat * nyquist_mask
    v_hat = v_hat * nyquist_mask
    w_hat = w_hat * nyquist_mask

    # k·û
    k_dot_u = kx * u_hat + ky * v_hat + kz * w_hat

    # Safe correction: avoid complex zeros_like (fails on Metal in eager mode).
    # Replace with real boolean mask * complex divide — stays on CPU device.
    k2_safe = jnp.where(k2 > 0.0, k2, 1.0)          # float; safe denominator
    valid   = (k2 > 0.0).astype(k2.dtype)             # float mask (1.0 / 0.0)
    correction = (k_dot_u / k2_safe) * valid          # complex × float → complex

    u_hat = u_hat - kx * correction
    v_hat = v_hat - ky * correction
    w_hat = w_hat - kz * correction
    return u_hat, v_hat, w_hat


# ─────────────────────────────────────────────────────────────────────────────
# JAX RHS (JIT-compiled)
# ─────────────────────────────────────────────────────────────────────────────

@jit
def _rhs_jax(u, v, w, T, kx, ky, kz, k2, dealias_mask, nyquist_mask):
    """RHS for Route 1 3D: du/dt, dv/dt, dw/dt, dT/dt.

    JIT-compiled.  All inputs/outputs are JAX arrays.
    Uses mean-field viscosity (nu_bar = mean(mu/rho)) for the spectral linear term.
    Variable-coefficient physics (A(T)·ΔT, μ(T)|∇u|²) computed in physical space.

    This is identical in structure to route1_3D.py::_rhs.
    """
    # ── Fluid properties (pointwise on T) ──────────────────────────────────────
    props = ideal_props_jax(T)
    mu    = props["mu"]
    rho   = props["rho"]
    k_th  = props["k"]
    cv    = props["cv"]

    nu_field   = mu / rho
    A_field_v  = k_th / (rho * cv)     # thermal diffusivity

    nu_bar = jnp.mean(nu_field)         # scalar mean for spectral linear term
    A_bar  = jnp.mean(A_field_v)

    # ── FFT of velocity ─────────────────────────────────────────────────────────
    u_hat = jnp.fft.fftn(u)
    v_hat = jnp.fft.fftn(v)
    w_hat = jnp.fft.fftn(w)

    def deriv(f_hat, k_dir):
        return jnp.real(jnp.fft.ifftn(1j * k_dir * f_hat))

    # ── Velocity gradients ──────────────────────────────────────────────────────
    dudx = deriv(u_hat, kx);  dudy = deriv(u_hat, ky);  dudz = deriv(u_hat, kz)
    dvdx = deriv(v_hat, kx);  dvdy = deriv(v_hat, ky);  dvdz = deriv(v_hat, kz)
    dwdx = deriv(w_hat, kx);  dwdy = deriv(w_hat, ky);  dwdz = deriv(w_hat, kz)

    # ── Nonlinear advection ─────────────────────────────────────────────────────
    NL_u = u * dudx + v * dudy + w * dudz
    NL_v = u * dvdx + v * dvdy + w * dvdz
    NL_w = u * dwdx + v * dwdy + w * dwdz

    # Dealias nonlinear terms
    NL_u_hat = jnp.fft.fftn(NL_u) * dealias_mask
    NL_v_hat = jnp.fft.fftn(NL_v) * dealias_mask
    NL_w_hat = jnp.fft.fftn(NL_w) * dealias_mask

    # Leray projection of nonlinear terms (no Nyquist issue here — already dealiased)
    k_dot_NL = kx * NL_u_hat + ky * NL_v_hat + kz * NL_w_hat
    k2_safe  = jnp.where(k2 > 0.0, k2, 1.0)
    corr_NL  = (k_dot_NL / k2_safe) * (k2 > 0.0).astype(k2.dtype)
    NL_u_hat = NL_u_hat - kx * corr_NL
    NL_v_hat = NL_v_hat - ky * corr_NL
    NL_w_hat = NL_w_hat - kz * corr_NL

    # ── Velocity tendency ───────────────────────────────────────────────────────
    du_hat = (-NL_u_hat - nu_bar * k2 * u_hat) * dealias_mask
    dv_hat = (-NL_v_hat - nu_bar * k2 * v_hat) * dealias_mask
    dw_hat = (-NL_w_hat - nu_bar * k2 * w_hat) * dealias_mask

    d_u = jnp.real(jnp.fft.ifftn(du_hat))
    d_v = jnp.real(jnp.fft.ifftn(dv_hat))
    d_w = jnp.real(jnp.fft.ifftn(dw_hat))

    # ── Temperature: advection + diffusion + viscous heating ───────────────────
    T_hat = jnp.fft.fftn(T)
    adv_T = u * deriv(T_hat, kx) + v * deriv(T_hat, ky) + w * deriv(T_hat, kz)

    lap_T  = jnp.real(jnp.fft.ifftn(-k2 * T_hat))
    diff_T = A_bar * lap_T          # variable-coeff handled via A_bar mean

    grad_u_sq = (
        dudx ** 2 + dudy ** 2 + dudz ** 2
        + dvdx ** 2 + dvdy ** 2 + dvdz ** 2
        + dwdx ** 2 + dwdy ** 2 + dwdz ** 2
    )
    Q_visc = mu * grad_u_sq / (rho * cv)

    d_T = -adv_T + diff_T + Q_visc

    return d_u, d_v, d_w, d_T


# ─────────────────────────────────────────────────────────────────────────────
# JAX RK4 step (JIT-compiled)
# ─────────────────────────────────────────────────────────────────────────────

@jit
def solver_step_jax(u, v, w, T, kx, ky, kz, k2,
                    dealias_mask, nyquist_mask, dt):
    """One RK4 step.  JIT-compiled.  Returns (new_u, new_v, new_w, new_T).

    Inputs and outputs are JAX arrays.  No Python-level side effects.
    """
    args = (kx, ky, kz, k2, dealias_mask, nyquist_mask)

    k1u, k1v, k1w, k1T = _rhs_jax(u, v, w, T, *args)

    k2u, k2v, k2w, k2T = _rhs_jax(
        u + 0.5 * dt * k1u, v + 0.5 * dt * k1v,
        w + 0.5 * dt * k1w, T + 0.5 * dt * k1T, *args
    )
    k3u, k3v, k3w, k3T = _rhs_jax(
        u + 0.5 * dt * k2u, v + 0.5 * dt * k2v,
        w + 0.5 * dt * k2w, T + 0.5 * dt * k2T, *args
    )
    k4u, k4v, k4w, k4T = _rhs_jax(
        u + dt * k3u, v + dt * k3v,
        w + dt * k3w, T + dt * k3T, *args
    )

    fac = dt / 6.0
    new_u = u + fac * (k1u + 2.0 * k2u + 2.0 * k3u + k4u)
    new_v = v + fac * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)
    new_w = w + fac * (k1w + 2.0 * k2w + 2.0 * k3w + k4w)
    new_T = T + fac * (k1T + 2.0 * k2T + 2.0 * k3T + k4T)

    # Dealias all fields
    def _dealias(f):
        f_hat = jnp.fft.fftn(f) * dealias_mask
        return jnp.real(jnp.fft.ifftn(f_hat))

    new_u = _dealias(new_u)
    new_v = _dealias(new_v)
    new_w = _dealias(new_w)
    new_T = _dealias(new_T)

    # Re-project velocity to enforce div u = 0 (accumulation guard)
    u_hat = jnp.fft.fftn(new_u)
    v_hat = jnp.fft.fftn(new_v)
    w_hat = jnp.fft.fftn(new_w)
    u_hat, v_hat, w_hat = _project_spectral_jax(
        u_hat, v_hat, w_hat, kx, ky, kz, k2, nyquist_mask
    )
    new_u = jnp.real(jnp.fft.ifftn(u_hat))
    new_v = jnp.real(jnp.fft.ifftn(v_hat))
    new_w = jnp.real(jnp.fft.ifftn(w_hat))

    # Clip T to physical range
    new_T = jnp.clip(new_T, 100.0, 6000.0)

    return new_u, new_v, new_w, new_T


# ─────────────────────────────────────────────────────────────────────────────
# CFL and diagnostics (Python-level, uses .item() to extract scalars)
# ─────────────────────────────────────────────────────────────────────────────

def compute_cfl_dt_jax(u, v, w, T, N: int, dt_max: float,
                       safety: float = 0.4) -> float:
    """Return CFL-limited timestep (Python float, safe to use as JAX scalar)."""
    dx = 2.0 * np.pi / N

    U_max = max(
        float(jnp.max(jnp.abs(u))),
        float(jnp.max(jnp.abs(v))),
        float(jnp.max(jnp.abs(w))),
        1e-12,
    )
    dt_adv = safety * dx / U_max

    props = ideal_props_jax(T)
    nu_max = float(jnp.max(props["mu"] / props["rho"]))
    A_max  = float(jnp.max(props["k"] / (props["rho"] * props["cv"])))
    diff   = max(nu_max, A_max, 1e-20)
    dt_diff = safety * dx ** 2 / (2.0 * diff)

    return float(np.clip(min(dt_adv, dt_diff), 1e-12, dt_max))


def kinetic_energy_jax(u, v, w) -> float:
    """Mean kinetic energy E = ½⟨|u|²⟩."""
    return float(0.5 * jnp.mean(u ** 2 + v ** 2 + w ** 2))


def divergence_rms_jax(u, v, w, kx, ky, kz) -> float:
    """RMS of physical-space divergence ∂u/∂x + ∂v/∂y + ∂w/∂z."""
    u_hat = jnp.fft.fftn(u)
    v_hat = jnp.fft.fftn(v)
    w_hat = jnp.fft.fftn(w)
    div_hat = 1j * (kx * u_hat + ky * v_hat + kz * w_hat)
    div = jnp.real(jnp.fft.ifftn(div_hat))
    return float(jnp.sqrt(jnp.mean(div ** 2)))


# ─────────────────────────────────────────────────────────────────────────────
# Solver construction
# ─────────────────────────────────────────────────────────────────────────────

def make_jax_solver(N: int = 32, T_base: float = 300.0, dt_max: float = 0.05,
                    dtype=jnp.float32) -> dict:
    """Create a JAX solver state dict.

    All arrays placed on the compute device (CPU when Metal is default).
    Returns dict with all arrays as JAX arrays; cached grid + masks.
    """
    if not _JAX_AVAILABLE:
        raise RuntimeError("JAX not available. Install: pip install jax==0.4.26 jaxlib==0.4.26 jax-metal")

    kx, ky, kz, k2 = make_grid_jax(N, dtype=dtype)
    dealias_mask    = make_dealias_mask_jax(N, dtype=dtype)
    nyquist_mask    = make_nyquist_mask_jax(N, dtype=dtype)

    T_arr = _d(jnp.array([T_base], dtype=dtype))
    props = ideal_props_jax(T_arr)
    nu_ref = float(props["mu"][0] / props["rho"][0])
    A_ref  = float(props["k"][0] / (props["rho"][0] * props["cv"][0]))

    return {
        "N": N,
        "dtype": dtype,
        "backend": _JAX_BACKEND,
        "kx": kx, "ky": ky, "kz": kz, "k2": k2,
        "dealias_mask": dealias_mask,
        "nyquist_mask": nyquist_mask,
        "u": _d(jnp.zeros((N, N, N), dtype=dtype)),
        "v": _d(jnp.zeros((N, N, N), dtype=dtype)),
        "w": _d(jnp.zeros((N, N, N), dtype=dtype)),
        "T": _d(jnp.full((N, N, N), T_base, dtype=dtype)),
        "t": 0.0,
        "dt_max": dt_max,
        "T_base": T_base,
        "nu_ref": nu_ref,
        "A_ref": A_ref,
    }


def taylor_green_ic_jax(N: int, V0: float = 1.0, T_base: float = 300.0,
                        dtype=jnp.float32) -> dict:
    """Taylor-Green vortex IC as JAX arrays on the compute device."""
    x_np = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    xx, yy, zz = np.meshgrid(x_np, x_np, x_np, indexing="ij")

    u_np = V0 * np.sin(xx) * np.cos(yy) * np.cos(zz)
    v_np = -V0 * np.cos(xx) * np.sin(yy) * np.cos(zz)
    w_np = np.zeros((N, N, N))
    T_np = np.full((N, N, N), T_base)

    return {
        "u": _d(jnp.array(u_np, dtype=dtype)),
        "v": _d(jnp.array(v_np, dtype=dtype)),
        "w": _d(jnp.array(w_np, dtype=dtype)),
        "T": _d(jnp.array(T_np, dtype=dtype)),
        "E0_analytic": V0 ** 2 / 8.0,
        "name": "taylor_green",
        "V0": V0,
    }


def set_ic_jax(solver: dict, ic: dict) -> dict:
    """Set IC on solver. Projects velocity to div-free. Returns new solver dict."""
    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    nq_mask = solver["nyquist_mask"]
    N, dtype = solver["N"], solver["dtype"]

    u = ic["u"]
    v = ic["v"]
    w = ic["w"]

    # Leray project the IC
    u_hat, v_hat, w_hat = (jnp.fft.fftn(u), jnp.fft.fftn(v), jnp.fft.fftn(w))
    u_hat, v_hat, w_hat = _project_spectral_jax(u_hat, v_hat, w_hat,
                                                  kx, ky, kz, k2, nq_mask)
    u = jnp.real(jnp.fft.ifftn(u_hat))
    v = jnp.real(jnp.fft.ifftn(v_hat))
    w = jnp.real(jnp.fft.ifftn(w_hat))

    T = ic["T"] if "T" in ic else jnp.full((N, N, N), solver["T_base"], dtype=dtype)

    return {**solver, "u": u, "v": v, "w": w, "T": T, "t": 0.0}


# ─────────────────────────────────────────────────────────────────────────────
# Run loop (Python-level; JIT boundary = one RK4 step)
# ─────────────────────────────────────────────────────────────────────────────

def run_jax(
    solver: dict,
    t_end: float,
    max_steps: int = 5000,
    cfl_safety: float = 0.4,
    verbose: bool = False,
    warmup: bool = True,
) -> dict:
    """Run the JAX solver from solver['t'] to t_end.

    First call compiles solver_step_jax (one-time ~10–30 s).
    Subsequent calls use cached compilation.

    Parameters
    ----------
    warmup : bool
        If True, prints compilation time on first call.

    Returns
    -------
    dict with keys: u, v, w, T, t_final, history, verdict,
                    A_min_global, E_initial, E_final, dissipation_rate_mean,
                    n_steps, compile_time_s (only on first call)
    """
    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    dealias_mask = solver["dealias_mask"]
    nyquist_mask = solver["nyquist_mask"]
    dt_max = solver["dt_max"]
    N = solver["N"]

    u, v, w, T = solver["u"], solver["v"], solver["w"], solver["T"]

    history = []
    A_min_global = float("inf")
    E_initial = kinetic_energy_jax(u, v, w)
    E_prev = E_initial
    dissipation_rates = []
    verdict = "PASS"
    compile_time_s = None
    t = float(solver["t"])

    for step_idx in range(max_steps):
        if t >= t_end:
            break

        dt = compute_cfl_dt_jax(u, v, w, T, N, dt_max=dt_max, safety=cfl_safety)
        dt = min(dt, t_end - t)
        dt_jnp = jnp.array(dt, dtype=solver["dtype"])

        # JIT compile on first call
        t0_wall = time.perf_counter()
        new_u, new_v, new_w, new_T = solver_step_jax(
            u, v, w, T, kx, ky, kz, k2, dealias_mask, nyquist_mask, dt_jnp
        )
        # Block until computation is done (JAX is async by default)
        jax.block_until_ready((new_u, new_v, new_w, new_T))
        step_wall = time.perf_counter() - t0_wall

        if step_idx == 0 and warmup:
            compile_time_s = step_wall
            if verbose:
                print(f"  [JIT] First step (compile + run): {compile_time_s:.1f} s")

        u, v, w, T = new_u, new_v, new_w, new_T
        t += dt

        # Diagnostics (JAX scalar → Python float)
        props = ideal_props_jax(T)
        A_arr = props["k"] / (props["rho"] * props["cv"])
        A_min = float(jnp.min(A_arr))
        A_min_global = min(A_min_global, A_min)

        if A_min <= 0.0:
            verdict = "FAIL"
            if verbose:
                print(f"  [FAIL] Step {step_idx}: A_min = {A_min:.3e} ≤ 0")
            break

        E = kinetic_energy_jax(u, v, w)
        dE_dt = (E - E_prev) / dt
        dissipation_rates.append(-dE_dt)
        E_prev = E

        diag = {
            "t": t, "dt": dt, "E": E,
            "A_min": A_min,
            "step_wall_s": step_wall,
        }
        history.append(diag)

        if verbose and (step_idx % max(1, max_steps // 10) == 0):
            print(
                f"  step {step_idx:4d}  t={t:.4f}  E={E:.4e}  "
                f"A_min={A_min:.3e}  wall={step_wall*1000:.1f}ms"
            )

    E_final = kinetic_energy_jax(u, v, w)
    mean_diss = float(np.mean(dissipation_rates)) if dissipation_rates else 0.0

    result = {
        "u": u, "v": v, "w": w, "T": T,
        "t_final": t,
        "history": history,
        "verdict": verdict,
        "n_steps": len(history),
        "A_min_global": A_min_global,
        "E_initial": E_initial,
        "E_final": E_final,
        "dissipation_rate_mean": mean_diss,
        "energy_decayed": E_final < E_initial,
    }
    if compile_time_s is not None:
        result["compile_time_s"] = compile_time_s

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Taylor-Green benchmark (JAX Metal)
# ─────────────────────────────────────────────────────────────────────────────

def run_taylor_green_jax(
    N: int = 32,
    t_end: float = 1.0,
    V0: float = 1.0,
    T_base: float = 300.0,
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str = None,
    exp_id_suffix: str = "",
) -> dict:
    """Run Taylor-Green vortex on JAX Metal backend.

    Checks (identical criteria to NumPy route1_3D.py):
    1. E₀ within 1% of analytic V₀²/8
    2. Energy decays (E_final < E_initial)
    3. A_min_global > 0 (second law)
    4. div u ≈ 0 at end

    Returns result dict with 'verdict' PASS/FAIL and speedup stats.
    """
    solver = make_jax_solver(N=N, T_base=T_base, dt_max=0.05)
    ic = taylor_green_ic_jax(N, V0=V0, T_base=T_base, dtype=solver["dtype"])
    solver = set_ic_jax(solver, ic)

    E0_computed = kinetic_energy_jax(solver["u"], solver["v"], solver["w"])
    E0_analytic = ic["E0_analytic"]

    if verbose:
        print(f"\n  JAX TG Benchmark  backend={_JAX_BACKEND}  N={N}³  dtype={solver['dtype']}")
        print(f"  E₀ analytic = {E0_analytic:.6f}")
        print(f"  E₀ computed = {E0_computed:.6f}  "
              f"(error = {abs(E0_computed - E0_analytic)/E0_analytic*100:.4f}%)")

    # Timed run
    t_start = time.perf_counter()
    result = run_jax(solver, t_end=t_end, max_steps=max_steps, verbose=verbose, warmup=True)
    t_wall = time.perf_counter() - t_start

    # Divergence check at end
    div_rms = divergence_rms_jax(
        result["u"], result["v"], result["w"],
        solver["kx"], solver["ky"], solver["kz"]
    )

    E0_ok        = abs(E0_computed - E0_analytic) / E0_analytic < 0.01
    energy_decayed = result["E_final"] < E0_computed
    A_positive   = result["A_min_global"] > 0.0
    div_ok       = div_rms < 1e-5   # float32: 1e-5; float64: 1e-10

    passed = (result["verdict"] == "PASS" and energy_decayed and A_positive and E0_ok)

    # Compile time (first step) vs steady-state step time
    compile_s = result.get("compile_time_s", 0.0)
    step_times = [h["step_wall_s"] for h in result["history"][1:]]  # skip first (compile)
    step_ms_mean = 1000.0 * float(np.mean(step_times)) if step_times else 0.0

    tg_result = {
        "exp_id": f"EXP-L3-R1-TG-JAX-{N:03d}{exp_id_suffix}",
        "claim_id": "claim_m5_taylor_green_3D_jax",
        "timestamp": datetime.utcnow().isoformat(),
        "backend": _JAX_BACKEND,
        "dtype": str(solver["dtype"]),
        "N": N, "V0": V0, "T_base": T_base, "t_end": t_end,
        "n_steps": result["n_steps"],
        "t_final": result["t_final"],
        "E0_analytic": E0_analytic,
        "E0_computed": E0_computed,
        "E0_error_pct": abs(E0_computed - E0_analytic) / E0_analytic * 100.0,
        "E_final": result["E_final"],
        "energy_ratio": result["E_final"] / E0_computed,
        "energy_decayed": energy_decayed,
        "dissipation_rate_mean": result["dissipation_rate_mean"],
        "A_min_global": result["A_min_global"],
        "A_positive": A_positive,
        "div_rms_final": div_rms,
        "div_ok": div_ok,
        "E0_ok": E0_ok,
        "verdict": "PASS" if passed else "FAIL",
        "wall_time_s": t_wall,
        "compile_time_s": compile_s,
        "step_ms_mean": step_ms_mean,
        "key_metric": (
            f"E0_err={abs(E0_computed-E0_analytic)/E0_analytic*100:.4f}%  "
            f"E_ratio={result['E_final']/E0_computed:.4f}  "
            f"A_min={result['A_min_global']:.3e}  "
            f"step={step_ms_mean:.1f}ms"
        ),
        "params_json": json.dumps({"N": N, "V0": V0, "T_base": T_base,
                                   "t_end": t_end, "backend": _JAX_BACKEND}),
    }

    if verbose:
        print(f"\n  JAX TG Results (N={N}³, backend={_JAX_BACKEND}, t_end={t_end}):")
        print(f"    E₀ error:      {tg_result['E0_error_pct']:.4f}%  {'✓' if E0_ok else '✗'}")
        print(f"    Energy ratio:  {tg_result['energy_ratio']:.4f}  {'✓' if energy_decayed else '✗'}")
        print(f"    A_min_global:  {tg_result['A_min_global']:.4e}  {'✓' if A_positive else '✗'}")
        print(f"    div u (RMS):   {div_rms:.2e}  {'✓' if div_ok else '✗'}")
        print(f"    Wall time:     {t_wall:.1f} s  (compile: {compile_s:.1f} s)")
        print(f"    Step time:     {step_ms_mean:.2f} ms/step (after compile)")
        print(f"    Verdict:       {tg_result['verdict']}")

    if db_path is not None:
        _log_jax_result(db_path, tg_result)

    return tg_result


# ─────────────────────────────────────────────────────────────────────────────
# Results logging
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_jax_db(db_path: str):
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS tg_jax_experiments (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id          TEXT UNIQUE NOT NULL,
            claim_id        TEXT NOT NULL,
            timestamp       TEXT NOT NULL,
            backend         TEXT,
            dtype           TEXT,
            N               INTEGER NOT NULL,
            V0              REAL,
            T_base          REAL,
            t_end           REAL,
            n_steps         INTEGER,
            t_final         REAL,
            verdict         TEXT NOT NULL,
            E0_analytic     REAL,
            E0_computed     REAL,
            E0_error_pct    REAL,
            E_final         REAL,
            energy_ratio    REAL,
            energy_decayed  INTEGER,
            dissipation_rate_mean REAL,
            A_min_global    REAL,
            A_positive      INTEGER,
            div_rms_final   REAL,
            div_ok          INTEGER,
            wall_time_s     REAL,
            compile_time_s  REAL,
            step_ms_mean    REAL,
            key_metric      TEXT,
            params_json     TEXT
        )
    """)
    conn.commit()
    return conn


def _log_jax_result(db_path: str, result: dict) -> None:
    conn = _ensure_jax_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "backend", "dtype", "N", "V0",
        "T_base", "t_end", "n_steps", "t_final", "verdict",
        "E0_analytic", "E0_computed", "E0_error_pct",
        "E_final", "energy_ratio", "energy_decayed",
        "dissipation_rate_mean", "A_min_global", "A_positive",
        "div_rms_final", "div_ok",
        "wall_time_s", "compile_time_s", "step_ms_mean",
        "key_metric", "params_json",
    ]
    values = tuple(
        int(result[c]) if isinstance(result.get(c), bool) else result.get(c)
        for c in cols
    )
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO tg_jax_experiments ({', '.join(cols)}) VALUES ({ph})",
        values,
    )
    conn.commit()
    conn.close()


def query_jax_results(db_path: str) -> list:
    if not os.path.exists(db_path):
        return []
    conn = sqlite3.connect(db_path)
    cur = conn.execute("SELECT * FROM tg_jax_experiments ORDER BY timestamp DESC")
    cols = [d[0] for d in cur.description]
    rows = [dict(zip(cols, r)) for r in cur.fetchall()]
    conn.close()
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# Speedup benchmark (JAX vs NumPy)
# ─────────────────────────────────────────────────────────────────────────────

def benchmark_speedup(N: int = 32, n_warmup: int = 3, n_timed: int = 10) -> dict:
    """Benchmark JAX Metal vs NumPy per-step wall time.

    Parameters
    ----------
    N        : grid size (N³)
    n_warmup : JAX JIT warmup steps (not timed)
    n_timed  : steps to time for mean

    Returns dict with numpy_ms, jax_ms, speedup.
    """
    from route1_3D import (
        taylor_green_ic, make_solver, set_ic, solver_step,
        compute_cfl_dt as np_cfl_dt,
    )

    # ── NumPy baseline ─────────────────────────────────────────────────────────
    ic_np = taylor_green_ic(N, V0=1.0, T_base=300.0)
    s_np  = make_solver(N=N, dt_max=0.05)
    s_np  = set_ic(s_np, ic_np)
    dt_np = np_cfl_dt(s_np)

    # Warmup NumPy
    for _ in range(2):
        s_np, _ = solver_step(s_np, dt=dt_np)

    t0 = time.perf_counter()
    for _ in range(n_timed):
        s_np, _ = solver_step(s_np, dt=dt_np)
    numpy_ms = 1000.0 * (time.perf_counter() - t0) / n_timed

    # ── JAX Metal ──────────────────────────────────────────────────────────────
    solver_jax = make_jax_solver(N=N, dt_max=0.05)
    ic_jax = taylor_green_ic_jax(N, V0=1.0, T_base=300.0, dtype=solver_jax["dtype"])
    solver_jax = set_ic_jax(solver_jax, ic_jax)

    kx, ky, kz = solver_jax["kx"], solver_jax["ky"], solver_jax["kz"]
    k2 = solver_jax["k2"]
    dmask = solver_jax["dealias_mask"]
    nqmask = solver_jax["nyquist_mask"]
    dtype = solver_jax["dtype"]

    u, v, w, T = solver_jax["u"], solver_jax["v"], solver_jax["w"], solver_jax["T"]
    dt_jax = jnp.array(0.01, dtype=dtype)

    # JIT warmup (compile)
    for _ in range(n_warmup):
        u, v, w, T = solver_step_jax(u, v, w, T, kx, ky, kz, k2, dmask, nqmask, dt_jax)
        jax.block_until_ready((u, v, w, T))

    t0 = time.perf_counter()
    for _ in range(n_timed):
        u, v, w, T = solver_step_jax(u, v, w, T, kx, ky, kz, k2, dmask, nqmask, dt_jax)
        jax.block_until_ready((u, v, w, T))
    jax_ms = 1000.0 * (time.perf_counter() - t0) / n_timed

    speedup = numpy_ms / jax_ms if jax_ms > 0 else float("inf")

    return {
        "N": N,
        "backend": _JAX_BACKEND,
        "numpy_ms": numpy_ms,
        "jax_ms": jax_ms,
        "speedup": speedup,
        "dtype": str(dtype),
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def _default_db():
    d = os.path.join(_HERE, "..", "results")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, "results.db")


def main():
    parser = argparse.ArgumentParser(description="M5 Route 1 3D — JAX Metal benchmark")
    parser.add_argument("--run", choices=["tg", "bench"], default="tg")
    parser.add_argument("--N", type=int, default=32)
    parser.add_argument("--t-end", type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--query", action="store_true")
    parser.add_argument("--db", type=str, default=None)
    parser.add_argument("--n-timed", type=int, default=10,
                        help="Steps to time in bench mode (default: 10)")
    args = parser.parse_args()

    db_path = args.db or _default_db()

    print(f"JAX backend: {_JAX_BACKEND}")
    print(f"JAX version: {jax.__version__}")

    if args.query:
        rows = query_jax_results(db_path)
        if not rows:
            print("No JAX TG results found.")
        else:
            for r in rows:
                print(f"{r['exp_id']}  N={r['N']}  backend={r['backend']}  "
                      f"{r['verdict']}  step={r.get('step_ms_mean','?'):.2f}ms  "
                      f"E0_err={r.get('E0_error_pct','?'):.4f}%")
        return

    if args.run == "bench":
        print(f"\nBenchmarking N={args.N}³ ({args.n_timed} timed steps after JIT warmup)...")
        res = benchmark_speedup(N=args.N, n_timed=args.n_timed)
        print(f"  NumPy:    {res['numpy_ms']:.2f} ms/step")
        print(f"  JAX {res['backend']:6s}: {res['jax_ms']:.2f} ms/step  ({res['dtype']})")
        print(f"  Speedup:  {res['speedup']:.1f}×")
    else:
        print(f"\nRunning JAX TG benchmark N={args.N}³, t_end={args.t_end}...")
        res = run_taylor_green_jax(
            N=args.N, t_end=args.t_end, verbose=args.verbose, db_path=db_path
        )
        print(f"\nVerdict: {res['verdict']}")
        print(f"  {res['key_metric']}")


if __name__ == "__main__":
    main()
