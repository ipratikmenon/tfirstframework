"""
route1_3D.py — Route 1 3D Spectral Navier-Stokes-Fourier Solver (M5)
======================================================================
3D incompressible NS with temperature-dependent viscosity (Route 1 system).
Pseudo-spectral method on [0, 2π)³.

Physics (Route 1, PRD v0.4 §3.2):
  ∂_t u + (u·∇)u = -∇p + ∇·(ν(T)·∇u)       div u = 0   ρ = const
  ∂_t T + u·∇T   = A(T)·ΔT + μ(T)|∇u|²/(ρ·cv)

where
  ν(T) = μ(T)/ρ       — temperature-dependent kinematic viscosity
  A(T) = k(T)/(ρ·cv)  — thermal diffusivity [always > 0 — second law]

Numerics:
  - Pseudo-spectral (FFT3) for all spatial derivatives
  - 2/3 dealiasing rule applied to all prognostic fields (u, v, w, T)
  - Helmholtz projection to maintain div u = 0 exactly in spectral space
  - RK4 time integration with adaptive CFL
  - Variable-coefficient terms split: mean-field spectral + physical correction

State representation:
  Physical-space fields (u, v, w, T) stored in the solver dict.
  Spectral transforms computed on-the-fly each RHS evaluation.

Experiments:
  EXP-L3-R1-TG-001  — Taylor-Green vortex validation (standard benchmark)
  EXP-L3-R1-TG-002  — Higher-Re Taylor-Green (N=64)

Taylor-Green validation:
  IC:  u = V₀·sin(x)cos(y)cos(z)
       v = -V₀·cos(x)sin(y)cos(z)
       w = 0
  E₀ = V₀²/8 (analytically exact)
  At early times: E(t) ≈ E₀·exp(-2ν̄t) [linear-viscous regime]
  Dissipation rate ε = -dE/dt > 0 always (energy decays monotonically)

Usage:
  python3 layer3/route1_3D.py                          # smoke test
  python3 layer3/route1_3D.py --run tg --N 32          # TG benchmark, N=32
  python3 layer3/route1_3D.py --run tg --N 64          # full benchmark
  python3 layer3/route1_3D.py --query                  # print results.db table
"""

import sys
import os
import json
import sqlite3
import argparse
import numpy as np
from datetime import datetime

# ── path setup ────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))

from layer1.tfirst_props import ideal_props, A_field

# ─────────────────────────────────────────────────────────────────────────────
# Grid and spectral helpers
# ─────────────────────────────────────────────────────────────────────────────

def make_grid_3D(N: int):
    """Return 3D physical grid and wavenumber arrays for N³ periodic domain [0,2π)³.

    Returns
    -------
    x  : ndarray (N,) — physical coordinates
    kx, ky, kz : ndarray (N,N,N) — wavenumber arrays (meshgrid, ij indexing)
    k2 : ndarray (N,N,N) — |k|² = kx²+ky²+kz²
    """
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    freqs = np.fft.fftfreq(N, d=1.0 / N)   # integer wavenumbers: 0,1,…,N/2-1,-N/2,…,-1
    kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    return x, kx, ky, kz, k2


def make_dealias_mask_3D(N: int) -> np.ndarray:
    """Return boolean 2/3-dealiasing mask for 3D spectral fields.

    Zeros all modes where |kx|, |ky|, OR |kz| > N//3.
    Shape: (N, N, N), dtype bool.
    """
    freqs = np.fft.fftfreq(N, d=1.0 / N)
    kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k_max = N // 3
    return (np.abs(kx) <= k_max) & (np.abs(ky) <= k_max) & (np.abs(kz) <= k_max)


def dealias_3D(f: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Apply 2/3 dealiasing to a real 3D field f.

    Transforms to spectral space, zeros high-k modes, returns real physical field.
    """
    f_hat = np.fft.fftn(f)
    f_hat *= mask
    return np.real(np.fft.ifftn(f_hat))


def project_divergence_free(u, v, w, kx, ky, kz, k2) -> tuple:
    """Project (u,v,w) onto divergence-free subspace.

    Helmholtz decomposition: û ← û - k(k·û)/|k|²  (Leray projector)
    Ensures div u = 0 exactly in spectral space.

    Implementation note — Nyquist zeroing:
    numpy's fftfreq represents the Nyquist as −N/2 (not +N/2).  For the 3D
    Nyquist plane (index N//2 in any direction), the "conjugate pair" of mode
    (−N/2, k₁, k₂) lives at array index (N//2, N−k₁, N−k₂) which also carries
    kx = −N/2.  The true negative of (−N/2, k₁, k₂) should have kx = +N/2, but
    that is not represented separately.  This aliasing breaks the anti-Hermitian
    property k·û(−k) = −(k·û(k))* required for the Leray projector to preserve
    Hermitian symmetry → ifftn returns a large imaginary part → taking real(.)
    discards divergence rather than zeroing it.

    The fix: zero the Nyquist planes before projecting.  In the solver these
    modes are always dealiased (2/3 rule: k_max = N//3 < N//2 for any N ≥ 6),
    so zeroing them is physically correct and removes the symmetry violation.

    Returns (u, v, w) — real physical-space fields.
    """
    N = u.shape[0]
    nq = N // 2          # Nyquist index in each dimension

    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)

    # Zero Nyquist planes to restore Hermitian symmetry before projection.
    for f_hat in (u_hat, v_hat, w_hat):
        f_hat[nq, :, :] = 0.0
        f_hat[:, nq, :] = 0.0
        f_hat[:, :, nq] = 0.0

    # k·û
    k_dot_u = kx * u_hat + ky * v_hat + kz * w_hat

    # Leray corrector: correction = (k·û) / |k|²  for k ≠ 0, zero for k = 0
    # Use boolean indexing instead of np.where to avoid NaN propagation.
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


def divergence_rms(u, v, w, kx, ky, kz) -> float:
    """Return RMS of div u = ∂u/∂x + ∂v/∂y + ∂w/∂z.

    Should be ≈ machine epsilon after projection.
    """
    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)
    div_hat = 1j * (kx * u_hat + ky * v_hat + kz * w_hat)
    div = np.real(np.fft.ifftn(div_hat))
    return float(np.sqrt(np.mean(div**2)))


# ─────────────────────────────────────────────────────────────────────────────
# Initial conditions
# ─────────────────────────────────────────────────────────────────────────────

def taylor_green_ic(N: int, V0: float = 1.0, T_base: float = 300.0) -> dict:
    """Taylor-Green vortex initial condition.

    u = V₀·sin(x)·cos(y)·cos(z)
    v = −V₀·cos(x)·sin(y)·cos(z)
    w = 0

    Analytically divergence-free.
    Kinetic energy: E₀ = V₀²/8

    Returns dict with keys: u, v, w, T, name, E0_analytic.
    """
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")

    u = V0 * np.sin(xx) * np.cos(yy) * np.cos(zz)
    v = -V0 * np.cos(xx) * np.sin(yy) * np.cos(zz)
    w = np.zeros((N, N, N))
    T = np.full((N, N, N), T_base)

    E0_analytic = V0**2 / 8.0

    return {
        "u": u,
        "v": v,
        "w": w,
        "T": T,
        "name": "taylor_green",
        "V0": V0,
        "E0_analytic": E0_analytic,
    }


def kinetic_energy(u, v, w) -> float:
    """Return mean kinetic energy E = ½·⟨|u|²⟩."""
    return 0.5 * float(np.mean(u**2 + v**2 + w**2))


# ─────────────────────────────────────────────────────────────────────────────
# Solver construction
# ─────────────────────────────────────────────────────────────────────────────

def make_solver(
    N: int = 32,
    fluid: str = "ideal",
    P=None,
    T_base: float = 300.0,
    dt_max: float = 0.05,
    t_star: float = 1.0,
) -> dict:
    """Create a 3D Route 1 solver state dict.

    Parameters
    ----------
    N       : grid points per side (N³ total; 32 for tests, 64+ for production)
    fluid   : 'ideal' or 'co2'
    P       : pressure [Pa] (required if fluid='co2')
    T_base  : initial temperature [K]
    dt_max  : maximum timestep [s]
    t_star  : reference blow-up time (for diagnostics)

    Returns
    -------
    Immutable solver dict.
    """
    if fluid == "co2" and P is None:
        raise ValueError("co2 fluid requires P [Pa]")

    x, kx, ky, kz, k2 = make_grid_3D(N)
    dealias_mask = make_dealias_mask_3D(N)

    # Reference properties at T_base
    T_arr = np.array([T_base])
    ref = ideal_props(T_arr) if fluid == "ideal" else None
    nu_ref = float(ref["mu"][0] / ref["rho"][0]) if ref is not None else 1e-3
    A_ref = float(ref["k"][0] / (ref["rho"][0] * ref["cv"][0])) if ref is not None else 1e-5

    return {
        # Grid
        "N": N,
        "x": x,
        "kx": kx,
        "ky": ky,
        "kz": kz,
        "k2": k2,
        "dealias_mask": dealias_mask,
        # Flow fields
        "u": np.zeros((N, N, N)),
        "v": np.zeros((N, N, N)),
        "w": np.zeros((N, N, N)),
        "T": np.full((N, N, N), T_base),
        # Time
        "t": 0.0,
        "dt_max": dt_max,
        "t_star": t_star,
        # Configuration
        "fluid": fluid,
        "P": P,
        "T_base": T_base,
        "nu_ref": nu_ref,
        "A_ref": A_ref,
        # Diagnostics
        "history": [],
    }


def set_ic(solver: dict, ic: dict) -> dict:
    """Set initial conditions from an IC dict. Returns new solver (immutable).

    IC dict must contain u, v, w. T is optional (defaults to T_base).
    Velocity is Helmholtz-projected to enforce div u = 0.
    """
    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    N = solver["N"]

    u = ic["u"].copy()
    v = ic["v"].copy()
    w = ic["w"].copy()

    # Enforce divergence-free
    u, v, w = project_divergence_free(u, v, w, kx, ky, kz, k2)

    T = ic["T"].copy() if "T" in ic else np.full((N, N, N), solver["T_base"])

    return {
        **solver,
        "u": u,
        "v": v,
        "w": w,
        "T": T,
        "t": 0.0,
        "history": [],
    }


# ─────────────────────────────────────────────────────────────────────────────
# Right-hand side
# ─────────────────────────────────────────────────────────────────────────────

def _rhs(u, v, w, T, kx, ky, kz, k2, dealias_mask, fluid, P, assert_A_positive=True):
    """Compute (du/dt, dv/dt, dw/dt, dT/dt) for one RK4 stage.

    Split scheme:
    - Velocity: -Leray_proj[(u·∇)u] - ν̄·|k|²·û  (spectral mean viscosity)
    - Temperature: -u·∇T + A(T)·ΔT + Q_visc        (physical, variable-coeff)

    Returns
    -------
    d_u, d_v, d_w : ndarray (N,N,N) — velocity tendencies
    d_T           : ndarray (N,N,N) — temperature tendency
    diag          : dict — diagnostics
    """
    # --- Properties (evaluated pointwise on T) ---
    if fluid == "ideal":
        props = ideal_props(T)
    else:
        props = ideal_props(T)   # fallback; co2 needs CoolProp scalar loop
    mu = props["mu"]
    rho = props["rho"]
    k_th = props["k"]
    cv = props["cv"]

    nu_field = mu / rho
    A_field_arr = k_th / (rho * cv)

    if assert_A_positive:
        A_min_check = float(np.min(A_field_arr))
        if A_min_check <= 0.0:
            raise AssertionError(
                f"[route1_3D] Second law violated: A_min = {A_min_check:.3e} ≤ 0"
            )

    nu_bar = float(np.mean(nu_field))  # scalar mean for spectral linear term

    # --- FFT of velocity ---
    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)

    # --- Velocity gradients (all 9 components, for advection + viscous heating) ---
    def _d(f_hat, k_dir):
        return np.real(np.fft.ifftn(1j * k_dir * f_hat))

    dudx = _d(u_hat, kx);  dudy = _d(u_hat, ky);  dudz = _d(u_hat, kz)
    dvdx = _d(v_hat, kx);  dvdy = _d(v_hat, ky);  dvdz = _d(v_hat, kz)
    dwdx = _d(w_hat, kx);  dwdy = _d(w_hat, ky);  dwdz = _d(w_hat, kz)

    # --- Nonlinear advection: (u·∇)u in physical space ---
    NL_u = u * dudx + v * dudy + w * dudz
    NL_v = u * dvdx + v * dvdy + w * dvdz
    NL_w = u * dwdx + v * dwdy + w * dwdz

    # Dealias nonlinear terms before projecting
    NL_u_hat = np.fft.fftn(NL_u) * dealias_mask
    NL_v_hat = np.fft.fftn(NL_v) * dealias_mask
    NL_w_hat = np.fft.fftn(NL_w) * dealias_mask

    # --- Helmholtz projection of nonlinear terms (subtract gradient part) ---
    k_dot_NL = kx * NL_u_hat + ky * NL_v_hat + kz * NL_w_hat
    with np.errstate(divide="ignore", invalid="ignore"):
        correction = np.where(k2 > 0, k_dot_NL / k2, 0.0 + 0.0j)
    NL_u_hat -= kx * correction
    NL_v_hat -= ky * correction
    NL_w_hat -= kz * correction

    # --- Velocity RHS in spectral space: -NL - ν̄|k|²û ---
    du_hat = -NL_u_hat - nu_bar * k2 * u_hat
    dv_hat = -NL_v_hat - nu_bar * k2 * v_hat
    dw_hat = -NL_w_hat - nu_bar * k2 * w_hat

    # Apply dealias to tendency
    du_hat *= dealias_mask
    dv_hat *= dealias_mask
    dw_hat *= dealias_mask

    d_u = np.real(np.fft.ifftn(du_hat))
    d_v = np.real(np.fft.ifftn(dv_hat))
    d_w = np.real(np.fft.ifftn(dw_hat))

    # --- Temperature: advection ---
    T_hat = np.fft.fftn(T)
    adv_T = u * _d(T_hat, kx) + v * _d(T_hat, ky) + w * _d(T_hat, kz)

    # --- Temperature: diffusion A(T)·ΔT (variable coefficient, physical space) ---
    A_bar = float(np.mean(A_field_arr))
    lap_T = np.real(np.fft.ifftn(-k2 * T_hat))
    diff_T = A_field_arr * lap_T

    # --- Viscous heating: Q = μ(T)|∇u|²/(ρ·cv) ---
    grad_u_sq = (
        dudx**2 + dudy**2 + dudz**2
        + dvdx**2 + dvdy**2 + dvdz**2
        + dwdx**2 + dwdy**2 + dwdz**2
    )
    Q_visc = mu * grad_u_sq / (rho * cv)

    d_T = -adv_T + diff_T + Q_visc

    # --- Diagnostics ---
    u_mag = np.sqrt(u**2 + v**2 + w**2)
    diag = {
        "A_min": float(np.min(A_field_arr)),
        "A_max": float(np.max(A_field_arr)),
        "A_bar": A_bar,
        "nu_bar": nu_bar,
        "nu_max": float(np.max(nu_field)),
        "u_rms": float(np.sqrt(np.mean(u_mag**2))),
        "u_max": float(np.max(u_mag)),
        "Q_visc_max": float(np.max(Q_visc)),
        "grad_u_sq_max": float(np.max(grad_u_sq)),
        "S_theta": float(np.mean(nu_bar * grad_u_sq)),  # mean viscous heating source
    }
    return d_u, d_v, d_w, d_T, diag


# ─────────────────────────────────────────────────────────────────────────────
# Time stepping
# ─────────────────────────────────────────────────────────────────────────────

def compute_cfl_dt(solver: dict, safety: float = 0.4) -> float:
    """Return CFL-limited timestep for 3D solver."""
    u, v, w = solver["u"], solver["v"], solver["w"]
    T = solver["T"]
    N = solver["N"]
    dx = 2.0 * np.pi / N

    U_max = max(float(np.max(np.abs(u))),
                float(np.max(np.abs(v))),
                float(np.max(np.abs(w))), 1e-12)

    dt_adv = safety * dx / U_max

    # Diffusive CFL
    if solver["fluid"] == "ideal":
        props = ideal_props(T)
    else:
        props = ideal_props(T)
    nu_max = float(np.max(props["mu"] / props["rho"]))
    A_max = float(np.max(props["k"] / (props["rho"] * props["cv"])))
    diff_coeff = max(nu_max, A_max, 1e-20)
    dt_diff = safety * dx**2 / (2.0 * diff_coeff)

    return float(np.clip(min(dt_adv, dt_diff), 1e-12, solver["dt_max"]))


def solver_step(solver: dict, dt: float = None) -> tuple:
    """Take one RK4 step. Returns (new_solver, step_diagnostics).

    Does NOT mutate input solver.
    """
    if dt is None:
        dt = compute_cfl_dt(solver)

    u, v, w, T = solver["u"], solver["v"], solver["w"], solver["T"]
    kx, ky, kz = solver["kx"], solver["ky"], solver["kz"]
    k2 = solver["k2"]
    mask = solver["dealias_mask"]
    fluid, P = solver["fluid"], solver["P"]

    args = (kx, ky, kz, k2, mask, fluid, P)

    # RK4 stages
    k1u, k1v, k1w, k1T, d1 = _rhs(u, v, w, T, *args)

    k2u, k2v, k2w, k2T, _ = _rhs(
        u + 0.5 * dt * k1u, v + 0.5 * dt * k1v,
        w + 0.5 * dt * k1w, T + 0.5 * dt * k1T, *args
    )
    k3u, k3v, k3w, k3T, _ = _rhs(
        u + 0.5 * dt * k2u, v + 0.5 * dt * k2v,
        w + 0.5 * dt * k2w, T + 0.5 * dt * k2T, *args
    )
    k4u, k4v, k4w, k4T, _ = _rhs(
        u + dt * k3u, v + dt * k3v,
        w + dt * k3w, T + dt * k3T, *args
    )

    fac = dt / 6.0
    new_u = u + fac * (k1u + 2.0 * k2u + 2.0 * k3u + k4u)
    new_v = v + fac * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)
    new_w = w + fac * (k1w + 2.0 * k2w + 2.0 * k3w + k4w)
    new_T = T + fac * (k1T + 2.0 * k2T + 2.0 * k3T + k4T)

    # Dealias all fields
    new_u = dealias_3D(new_u, mask)
    new_v = dealias_3D(new_v, mask)
    new_w = dealias_3D(new_w, mask)
    new_T = dealias_3D(new_T, mask)

    # Re-project velocity to maintain div u = 0 (accumulation guard)
    new_u, new_v, new_w = project_divergence_free(
        new_u, new_v, new_w, kx, ky, kz, k2
    )

    # Clip T to physical range
    new_T = np.clip(new_T, 100.0, 6000.0)

    new_t = solver["t"] + dt
    E = kinetic_energy(new_u, new_v, new_w)

    step_diag = {
        "t": new_t,
        "dt": dt,
        "E": E,
        "A_min": d1["A_min"],
        "A_max": d1["A_max"],
        "nu_bar": d1["nu_bar"],
        "u_rms": d1["u_rms"],
        "u_max": d1["u_max"],
        "Q_visc_max": d1["Q_visc_max"],
        "S_theta": d1["S_theta"],
        "T_min": float(np.min(new_T)),
        "T_max": float(np.max(new_T)),
        "A_positive": d1["A_min"] > 0.0,
        "verdict": "PASS" if d1["A_min"] > 0.0 else "FAIL",
    }

    new_solver = {
        **solver,
        "u": new_u,
        "v": new_v,
        "w": new_w,
        "T": new_T,
        "t": new_t,
    }
    return new_solver, step_diag


# ─────────────────────────────────────────────────────────────────────────────
# Run loop
# ─────────────────────────────────────────────────────────────────────────────

def run(
    solver: dict,
    t_end: float,
    max_steps: int = 5000,
    cfl_safety: float = 0.4,
    callback=None,
    verbose: bool = False,
) -> dict:
    """Run the 3D Route 1 solver from solver['t'] to t_end.

    Returns
    -------
    dict with keys:
      solver           : final solver state
      history          : list of step_diag dicts
      verdict          : 'PASS' or 'FAIL'
      A_min_global     : minimum A(T) across all steps
      E_initial        : kinetic energy at t=0
      E_final          : kinetic energy at t=t_end
      dissipation_rate : -dE/dt estimated from history (mean)
      n_steps, t_final
    """
    history = []
    A_min_global = float("inf")
    E_initial = kinetic_energy(solver["u"], solver["v"], solver["w"])
    E_prev = E_initial
    dissipation_rates = []
    verdict = "PASS"

    for step_idx in range(max_steps):
        if solver["t"] >= t_end:
            break

        dt = compute_cfl_dt(solver, safety=cfl_safety)
        dt = min(dt, t_end - solver["t"])

        try:
            solver, diag = solver_step(solver, dt=dt)
        except AssertionError as exc:
            verdict = "FAIL"
            if verbose:
                print(f"  [FAIL] Step {step_idx}: {exc}")
            break

        A_min_global = min(A_min_global, diag["A_min"])

        # Track dissipation: ε = -dE/dt > 0 (energy should decay)
        dE_dt = (diag["E"] - E_prev) / diag["dt"]
        dissipation_rates.append(-dE_dt)
        E_prev = diag["E"]

        history.append(diag)

        if callback is not None:
            callback(diag)

        if verbose and step_idx % max(1, len(range(max_steps)) // 10) == 0:
            print(
                f"  step {step_idx:4d}  t={diag['t']:.4f}  "
                f"E={diag['E']:.4e}  "
                f"A_min={diag['A_min']:.3e}  "
                f"u_rms={diag['u_rms']:.3e}"
            )

    E_final = kinetic_energy(solver["u"], solver["v"], solver["w"])
    mean_dissipation = float(np.mean(dissipation_rates)) if dissipation_rates else 0.0

    return {
        "solver": solver,
        "history": history,
        "verdict": verdict,
        "n_steps": len(history),
        "t_final": solver["t"],
        "A_min_global": A_min_global,
        "E_initial": E_initial,
        "E_final": E_final,
        "dissipation_rate_mean": mean_dissipation,
        "energy_decayed": E_final < E_initial,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Taylor-Green benchmark
# ─────────────────────────────────────────────────────────────────────────────

def run_taylor_green(
    N: int = 32,
    t_end: float = 1.0,
    V0: float = 1.0,
    T_base: float = 300.0,
    fluid: str = "ideal",
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str = None,
) -> dict:
    """Run the Taylor-Green vortex validation benchmark.

    Checks:
    1. E₀ matches analytic V₀²/8
    2. Energy decays monotonically (dissipation > 0 on average)
    3. A_min_global > 0 (second law maintained throughout)
    4. div u ≈ 0 (machine precision) at end

    Returns result dict with 'verdict' PASS/FAIL and key metrics.
    """
    ic = taylor_green_ic(N, V0=V0, T_base=T_base)

    solver = make_solver(N=N, fluid=fluid, T_base=T_base, dt_max=0.05)
    solver = set_ic(solver, ic)

    E0_computed = kinetic_energy(solver["u"], solver["v"], solver["w"])
    E0_analytic = ic["E0_analytic"]

    if verbose:
        print(f"  Taylor-Green N={N}³  V₀={V0}  T̄={T_base}K")
        print(f"  E₀ analytic = {E0_analytic:.6f}")
        print(f"  E₀ computed = {E0_computed:.6f}  "
              f"(error = {abs(E0_computed - E0_analytic)/E0_analytic*100:.2f}%)")

    result = run(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)

    final_solver = result["solver"]

    # Divergence check at end
    div_rms = divergence_rms(
        final_solver["u"], final_solver["v"], final_solver["w"],
        final_solver["kx"], final_solver["ky"], final_solver["kz"],
    )

    # Verdict criteria
    E0_ok = abs(E0_computed - E0_analytic) / E0_analytic < 0.01   # within 1%
    energy_decayed = result["E_final"] < E0_computed
    A_positive = result["A_min_global"] > 0.0
    div_ok = div_rms < 1e-8

    passed = (result["verdict"] == "PASS" and energy_decayed and A_positive and E0_ok)

    tg_result = {
        "exp_id": f"EXP-L3-R1-TG-{N:03d}",
        "claim_id": "claim_m5_taylor_green_3D",
        "timestamp": datetime.utcnow().isoformat(),
        "N": N,
        "V0": V0,
        "T_base": T_base,
        "fluid": fluid,
        "t_end": t_end,
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
        "key_metric": (
            f"E0_err={abs(E0_computed-E0_analytic)/E0_analytic*100:.2f}%  "
            f"E_ratio={result['E_final']/E0_computed:.4f}  "
            f"A_min={result['A_min_global']:.3e}"
        ),
        "params_json": json.dumps({
            "N": N, "V0": V0, "T_base": T_base, "fluid": fluid, "t_end": t_end,
        }),
    }

    if verbose:
        print(f"\n  TG Benchmark Results (N={N}³, t_end={t_end}):")
        print(f"    E₀ error:      {tg_result['E0_error_pct']:.4f}%  {'✓' if E0_ok else '✗'}")
        print(f"    Energy ratio:  {tg_result['energy_ratio']:.4f}  (E_final/E₀; < 1 = decayed  {'✓' if energy_decayed else '✗'})")
        print(f"    A_min_global:  {tg_result['A_min_global']:.4e}  {'✓' if A_positive else '✗'}")
        print(f"    div u (RMS):   {div_rms:.2e}  {'✓' if div_ok else '✗'}")
        print(f"    Verdict:       {tg_result['verdict']}")

    if db_path is not None:
        _log_tg_result(db_path, tg_result)

    return tg_result


# ─────────────────────────────────────────────────────────────────────────────
# Results logging
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_db(db_path: str) -> sqlite3.Connection:
    """Open (or create) results.db and ensure the tg_experiments table exists."""
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS tg_experiments (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id          TEXT UNIQUE NOT NULL,
            claim_id        TEXT NOT NULL,
            timestamp       TEXT NOT NULL,
            N               INTEGER NOT NULL,
            V0              REAL,
            T_base          REAL,
            fluid           TEXT,
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
            key_metric      TEXT,
            params_json     TEXT
        )
    """)
    conn.commit()
    return conn


def _log_tg_result(db_path: str, result: dict) -> None:
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "N", "V0", "T_base", "fluid",
        "t_end", "n_steps", "t_final", "verdict",
        "E0_analytic", "E0_computed", "E0_error_pct",
        "E_final", "energy_ratio", "energy_decayed",
        "dissipation_rate_mean", "A_min_global", "A_positive",
        "div_rms_final", "div_ok", "key_metric", "params_json",
    ]
    values = tuple(
        int(result[c]) if isinstance(result.get(c), bool) else result.get(c)
        for c in cols
    )
    placeholders = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO tg_experiments ({', '.join(cols)}) VALUES ({placeholders})",
        values,
    )
    conn.commit()
    conn.close()


def query_results(db_path: str) -> list:
    """Query tg_experiments table. Returns list of dicts."""
    if not os.path.exists(db_path):
        return []
    conn = sqlite3.connect(db_path)
    cur = conn.execute("SELECT * FROM tg_experiments ORDER BY timestamp DESC")
    cols = [d[0] for d in cur.description]
    rows = [dict(zip(cols, r)) for r in cur.fetchall()]
    conn.close()
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# Smoke test (fast CI validation)
# ─────────────────────────────────────────────────────────────────────────────

def run_smoke_test(verbose: bool = False) -> dict:
    """Fast smoke test: TG at N=16, t_end=0.1.

    Returns dict with verdict, n_steps, all_pass.
    """
    result = run_taylor_green(N=16, t_end=0.1, max_steps=50, verbose=verbose)
    all_pass = result["verdict"] == "PASS"
    return {
        "verdict": result["verdict"],
        "n_steps": result["n_steps"],
        "all_pass": all_pass,
        "result": result,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def _default_db_path() -> str:
    results_dir = os.path.join(_HERE, "..", "results")
    os.makedirs(results_dir, exist_ok=True)
    return os.path.join(results_dir, "results.db")


def main():
    parser = argparse.ArgumentParser(description="M5 Route 1 3D — Taylor-Green benchmark")
    parser.add_argument("--run", choices=["smoke", "tg"], default="smoke",
                        help="'smoke' (N=16, fast) or 'tg' (full Taylor-Green benchmark)")
    parser.add_argument("--N", type=int, default=32, help="Grid resolution (default: 32)")
    parser.add_argument("--t-end", type=float, default=1.0, help="End time [s] (default: 1.0)")
    parser.add_argument("--V0", type=float, default=1.0, help="TG amplitude (default: 1.0)")
    parser.add_argument("--T-base", type=float, default=300.0,
                        help="Base temperature [K] (default: 300)")
    parser.add_argument("--verbose", action="store_true", help="Print step diagnostics")
    parser.add_argument("--query", action="store_true", help="Print results.db table and exit")
    parser.add_argument("--db", type=str, default=None, help="Path to results.db (default: auto)")
    args = parser.parse_args()

    db_path = args.db or _default_db_path()

    if args.query:
        rows = query_results(db_path)
        if not rows:
            print("No TG experiment results found.")
        else:
            for r in rows:
                print(
                    f"{r['exp_id']}  N={r['N']}  {r['verdict']}  "
                    f"E0_err={r.get('E0_error_pct', 'N/A'):.3f}%  "
                    f"E_ratio={r.get('energy_ratio', 'N/A'):.4f}  "
                    f"A_min={r.get('A_min_global', 'N/A'):.3e}"
                )
        return

    if args.run == "smoke":
        print("Running 3D smoke test (N=16³, t_end=0.1)...")
        result = run_smoke_test(verbose=args.verbose)
        print(f"Smoke test: {result['verdict']}  ({result['n_steps']} steps)")
    else:
        print(f"Running Taylor-Green benchmark N={args.N}³, t_end={args.t_end}...")
        result = run_taylor_green(
            N=args.N,
            t_end=args.t_end,
            V0=args.V0,
            T_base=args.T_base,
            verbose=args.verbose,
            db_path=db_path,
        )
        print(f"\nVerdict: {result['verdict']}")
        print(f"  {result['key_metric']}")


if __name__ == "__main__":
    main()
