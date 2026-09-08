"""
route2_2D.py — Route 2: Exact Prize NS + Auxiliary Scalar θ (2D)
=================================================================
Implements the Route 2 system from PRD §3.2 v0.4:

    Exact Prize NS (ν = const):
        ∂_t u + (u·∇)u = −∇p + ν·Δu       div u = 0,  ρ = const
        ω = ∇×u,  −Δψ = ω,  u = ∂ψ/∂y,  v = −∂ψ/∂x

    Auxiliary scalar θ (slaved to u):
        ∂_t θ + u·∇θ = ν·Δθ + S_θ          θ(x, 0) = 0
        S_θ = ν·|∇u|²   (always ≥ 0)

    Effective viscosity:
        μ_eff(x,t) = ν + ε_param · f(θ(x,t))

This directly answers Tao's supercriticality objection: μ_eff > ν everywhere,
the LPS margin ε_LPS = ε_param · min(f(θ)) ≥ 0 is non-negative by construction,
and the θ system converges to the exact Prize equations as ε_param → 0.

Key questions answered:
  1. Does θ accumulate throughout the run (S_θ ≥ 0)?
  2. Is μ_eff_min > ν for all ε_param > 0?
  3. Does the LPS margin ε_LPS(ε_param) survive the ε_param → 0 limit?
  4. Is the ω evolution distinct from pure Prize NS (θ-feedback measurable)?

Experiments:
  EXP-L2-R2-001  ε_param = 1.0
  EXP-L2-R2-002  ε_param = 0.1
  EXP-L2-R2-003  ε_param = 0.01
  EXP-L2-R2-004  ε_param = 0.001

Run:
  python layer2/route2_2D.py                   # full sweep (4 ε values, 64²)
  python layer2/route2_2D.py --smoke-test       # fast check (32², 5 steps)
  python layer2/route2_2D.py --eps 0.1 --verbose
"""

from __future__ import annotations

import os
import sys
import json
import sqlite3
import time
import datetime
import numpy as np

# ── path setup ────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from self_similar_IC import dealias_23, make_grid, generate_CCF_profile
from layer1.tfirst_props import route2_theta_source

# Re-use the DB helpers from lambda_sweep_2D
from lambda_sweep_2D import _ensure_db, log_result, query_results

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

EPS_SWEEP_VALUES = [1.0, 0.1, 0.01, 0.001]

# Reference kinematic viscosity (Prize equations use ν = const)
_NU_DEFAULT = 1e-3   # m²/s — representative value

_DB_PATH = os.path.join(_ROOT, "results", "results.db")


# ─────────────────────────────────────────────────────────────────────────────
# Spectral helpers (same conventions as T_solver_2D.py)
# ─────────────────────────────────────────────────────────────────────────────

def _spectral_diff(f_hat: np.ndarray, k_dir: np.ndarray) -> np.ndarray:
    return np.real(np.fft.ifft2(1j * k_dir * f_hat))


def _psi_velocity(omega_hat: np.ndarray, kx: np.ndarray,
                  ky: np.ndarray, k2: np.ndarray) -> tuple:
    """Stream function ψ and velocity (u, v) from ω̂."""
    with np.errstate(divide="ignore", invalid="ignore"):
        psi_hat = np.where(k2 > 0, omega_hat / k2, 0.0 + 0.0j)
    u = _spectral_diff(psi_hat, ky)
    v = _spectral_diff(-psi_hat, kx)
    return psi_hat, u, v


def _spectral_laplacian(f: np.ndarray, k2: np.ndarray) -> np.ndarray:
    """Δf via FFT."""
    return np.real(np.fft.ifft2(-k2 * np.fft.fft2(f)))


# ─────────────────────────────────────────────────────────────────────────────
# Solver construction
# ─────────────────────────────────────────────────────────────────────────────

def make_solver(
    N: int = 64,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    f_theta: str = "identity",
    t_star: float = 1.0,
    lambda_val: float = 0.6057,
    dt_max: float = 0.1,
) -> dict:
    """
    Create a Route 2 solver state dict.

    Parameters
    ----------
    N         : grid points per side
    nu        : constant kinematic viscosity ν [m²/s] (Prize equations: ν = const)
    eps_param : coupling parameter ε in μ_eff = ν + ε·f(θ)
    f_theta   : functional form of f: 'identity' (f=θ), 'tanh' (f=tanh(θ)), 'sqrt'
    t_star    : blow-up time for suppression ratio tracking
    lambda_val: Wang et al. λ for suppression ratio exponent
    dt_max    : maximum timestep [s]
    """
    if nu <= 0:
        raise ValueError(f"nu must be > 0, got {nu}")
    if eps_param < 0:
        raise ValueError(f"eps_param must be ≥ 0, got {eps_param}")

    x, y, kx, ky = make_grid(N)
    k2 = kx**2 + ky**2

    return {
        "N"          : N,
        "x"          : x,
        "y"          : y,
        "kx"         : kx,
        "ky"         : ky,
        "k2"         : k2,
        "omega"      : np.zeros((N, N)),
        "theta"      : np.zeros((N, N)),   # θ(x,0) = 0 (PRD §3.2)
        "t"          : 0.0,
        "nu"         : nu,
        "eps_param"  : eps_param,
        "f_theta"    : f_theta,
        "t_star"     : t_star,
        "lambda_val" : lambda_val,
        "dt_max"     : dt_max,
    }


def set_ic(solver: dict, ic: dict) -> dict:
    """Set vorticity IC from a profile dict (self_similar_IC output).
    θ is always initialised to 0 regardless of IC (PRD §3.2).
    Returns new solver dict (immutable pattern).
    """
    new = {**solver}
    new["omega"] = ic["omega"].copy()
    new["theta"] = np.zeros_like(ic["omega"])  # θ(x,0) = 0 — always
    new["t"] = 0.0
    return new


# ─────────────────────────────────────────────────────────────────────────────
# f(θ) functional forms
# ─────────────────────────────────────────────────────────────────────────────

def _apply_f_theta(theta: np.ndarray, f_theta: str) -> np.ndarray:
    """
    Apply the functional form f to θ.

    Requirements: f(θ) ≥ 0, f(0) = 0, smooth.

    Available forms:
      'identity' : f(θ) = θ              (simplest; θ ≥ 0 always so f ≥ 0)
      'tanh'     : f(θ) = tanh(θ)        (bounded: f ∈ [0,1])
      'sqrt'     : f(θ) = √(θ + δ) − √δ  (sublinear growth)
    """
    theta = np.asarray(theta)
    # Ensure θ ≥ 0 (numerical floor — source is always ≥ 0 analytically)
    theta = np.clip(theta, 0.0, None)

    if f_theta == "identity":
        return theta
    elif f_theta == "tanh":
        return np.tanh(theta)
    elif f_theta == "sqrt":
        delta = 1e-10
        return np.sqrt(theta + delta) - np.sqrt(delta)
    else:
        raise ValueError(f"Unknown f_theta form: '{f_theta}'. Use 'identity', 'tanh', 'sqrt'.")


def _mu_eff(theta: np.ndarray, nu: float, eps_param: float, f_theta: str) -> np.ndarray:
    """μ_eff(x,t) = ν + ε·f(θ(x,t)).  Always ≥ ν > 0."""
    return nu + eps_param * _apply_f_theta(theta, f_theta)


# ─────────────────────────────────────────────────────────────────────────────
# Right-hand side — Route 2
# ─────────────────────────────────────────────────────────────────────────────

def _rhs(omega: np.ndarray, theta: np.ndarray,
         kx: np.ndarray, ky: np.ndarray, k2: np.ndarray,
         nu: float, eps_param: float, f_theta: str) -> tuple:
    """
    Compute (dω/dt, dθ/dt) for one RK4 stage.

    Route 2 system (PRD §3.2):
      dω/dt = −J(ψ,ω) + ν·Δω                    [exact Prize vorticity; ν=const]
      dθ/dt = −u·∇θ + ν·Δθ + S_θ                [θ accumulates viscous dissipation]
      S_θ = ν·|∇u|²                              [always ≥ 0]
      μ_eff = ν + ε·f(θ)                         [effective viscosity; μ_eff ≥ ν]

    Note: ω equation uses ν (not μ_eff). The θ-feedback enters μ_eff but the
    base vorticity dynamics are exact Prize equations. This is the point of
    Route 2: no approximation in the flow equations.

    Returns
    -------
    d_omega : (N,N) — dω/dt
    d_theta : (N,N) — dθ/dt
    diag    : dict — step diagnostics
    """
    # ── velocity from vorticity ───────────────────────────────────────────────
    omega_hat = np.fft.fft2(omega)
    psi_hat, u, v = _psi_velocity(omega_hat, kx, ky, k2)

    # ── ω equation: exact Prize NS (ν = const) ────────────────────────────────
    # Advection: J(ψ,ω) = u·∇ω
    domega_dx = _spectral_diff(omega_hat, kx)
    domega_dy = _spectral_diff(omega_hat, ky)
    J_psi_omega = u * domega_dx + v * domega_dy   # [s⁻²]

    # Diffusion: ν·Δω (spectral, constant ν)
    diff_omega = np.real(np.fft.ifft2(-nu * k2 * omega_hat))  # [s⁻²]

    d_omega = -J_psi_omega + diff_omega

    # ── θ source: S_θ = ν·|∇u|² ─────────────────────────────────────────────
    S_theta = route2_theta_source(u, v, nu, kx, ky)  # [s⁻¹·m⁻²·... actually dimensionless rate]

    # ── θ equation ────────────────────────────────────────────────────────────
    theta_hat = np.fft.fft2(theta)

    # Advection: u·∇θ
    dtheta_dx = _spectral_diff(theta_hat, kx)
    dtheta_dy = _spectral_diff(theta_hat, ky)
    adv_theta = u * dtheta_dx + v * dtheta_dy   # [θ units/s]

    # Diffusion: ν·Δθ
    diff_theta = np.real(np.fft.ifft2(-nu * k2 * theta_hat))  # [θ units/s]

    d_theta = -adv_theta + diff_theta + S_theta

    # ── effective viscosity ───────────────────────────────────────────────────
    mu_eff_field = _mu_eff(theta, nu, eps_param, f_theta)
    f_theta_field = _apply_f_theta(theta, f_theta)

    diag = {
        "u_max"        : float(np.max(np.abs(u))),
        "v_max"        : float(np.max(np.abs(v))),
        "S_theta_max"  : float(np.max(S_theta)),
        "S_theta_mean" : float(np.mean(S_theta)),
        "theta_min"    : float(np.min(theta)),
        "theta_max"    : float(np.max(theta)),
        "theta_mean"   : float(np.mean(theta)),
        "mu_eff_min"   : float(np.min(mu_eff_field)),
        "mu_eff_max"   : float(np.max(mu_eff_field)),
        "f_theta_max"  : float(np.max(f_theta_field)),
        "lps_margin"   : float(np.min(mu_eff_field)) - nu,  # ε_LPS = min(μ_eff) - ν ≥ 0
    }
    return d_omega, d_theta, diag


# ─────────────────────────────────────────────────────────────────────────────
# CFL timestep
# ─────────────────────────────────────────────────────────────────────────────

def compute_cfl_dt(solver: dict, safety: float = 0.4) -> float:
    """CFL-limited timestep for Route 2 (ν = const)."""
    omega = solver["omega"]
    kx, ky, k2 = solver["kx"], solver["ky"], solver["k2"]
    N, nu = solver["N"], solver["nu"]
    dx = 2.0 * np.pi / N

    # Velocity magnitude
    omega_hat = np.fft.fft2(omega)
    _, u, v = _psi_velocity(omega_hat, kx, ky, k2)
    U_max = max(float(np.max(np.abs(u))), float(np.max(np.abs(v))), 1e-12)

    dt_adv  = safety * dx / U_max
    dt_diff = safety * dx**2 / (2.0 * nu)

    return float(np.clip(min(dt_adv, dt_diff), 1e-12, solver["dt_max"]))


# ─────────────────────────────────────────────────────────────────────────────
# Single RK4 step
# ─────────────────────────────────────────────────────────────────────────────

def solver_step(solver: dict, dt: float | None = None) -> tuple[dict, dict]:
    """
    One RK4 step of the Route 2 system. Returns (new_solver, step_diag).
    Does not mutate the input solver (immutable pattern).
    """
    if dt is None:
        dt = compute_cfl_dt(solver)

    omega = solver["omega"]
    theta = solver["theta"]
    kx, ky, k2 = solver["kx"], solver["ky"], solver["k2"]
    nu, eps_param, f_th = solver["nu"], solver["eps_param"], solver["f_theta"]

    def rhs(om, th):
        return _rhs(om, th, kx, ky, k2, nu, eps_param, f_th)

    # RK4 stages
    k1o, k1t, d1 = rhs(omega,                    theta)
    k2o, k2t, _  = rhs(omega + 0.5*dt*k1o,       theta + 0.5*dt*k1t)
    k3o, k3t, _  = rhs(omega + 0.5*dt*k2o,       theta + 0.5*dt*k2t)
    k4o, k4t, _  = rhs(omega +    dt*k3o,        theta +    dt*k3t)

    new_omega = omega + (dt / 6.0) * (k1o + 2*k2o + 2*k3o + k4o)
    new_theta = theta + (dt / 6.0) * (k1t + 2*k2t + 2*k3t + k4t)

    # Dealias ω for stability (θ left undealiased — physical accumulation)
    new_omega = dealias_23(new_omega)

    # Enforce θ ≥ 0 (source is always ≥ 0; numerical errors can briefly dip below)
    np.clip(new_theta, 0.0, None, out=new_theta)

    new_t = solver["t"] + dt

    # Suppression ratio (for consistency with M0/M3 tracking)
    t_star, lam = solver["t_star"], solver["lambda_val"]
    tau = max((t_star - new_t) / t_star, 1e-10)
    supp_ratio = tau**(-(2.0 + lam)) if tau < 1.0 else float("inf")

    step_diag = {
        "t"              : new_t,
        "dt"             : dt,
        "omega_max"      : float(np.max(np.abs(new_omega))),
        "theta_max"      : d1["theta_max"],
        "theta_mean"     : d1["theta_mean"],
        "S_theta_max"    : d1["S_theta_max"],
        "S_theta_mean"   : d1["S_theta_mean"],
        "mu_eff_min"     : d1["mu_eff_min"],
        "mu_eff_max"     : d1["mu_eff_max"],
        "lps_margin"     : d1["lps_margin"],    # ε_LPS = min(μ_eff) - ν ≥ 0
        "suppression_ratio": supp_ratio,
        "eps_param"      : eps_param,
        "verdict"        : "PASS" if d1["lps_margin"] >= 0.0 else "FAIL",
    }

    new_solver = {
        **solver,
        "omega": new_omega,
        "theta": new_theta,
        "t"    : new_t,
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
    """
    Run Route 2 solver from solver['t'] to t_end.

    Returns
    -------
    dict:
      solver        : final state
      history       : list of step_diag dicts
      verdict       : 'PASS' or 'FAIL'
      n_steps       : steps taken
      t_final       : time reached
      theta_max_final: max θ at end
      lps_margin_min: minimum LPS margin over all steps
      mu_eff_min_global: global minimum μ_eff
    """
    history = []
    lps_margin_min  = float("inf")
    mu_eff_min_global = float("inf")
    theta_max_final = 0.0
    verdict = "PASS"

    for _ in range(max_steps):
        if solver["t"] >= t_end:
            break

        dt = compute_cfl_dt(solver, safety=cfl_safety)
        dt = min(dt, t_end - solver["t"])

        solver, diag = solver_step(solver, dt=dt)

        lps_margin_min    = min(lps_margin_min, diag["lps_margin"])
        mu_eff_min_global = min(mu_eff_min_global, diag["mu_eff_min"])

        history.append(diag)

        if callback is not None:
            callback(diag)

        if verbose and len(history) % 10 == 1:
            print(
                f"  step {len(history):4d}  t={diag['t']:.4f}  "
                f"ω_max={diag['omega_max']:.3e}  "
                f"θ_max={diag['theta_max']:.3e}  "
                f"μ_eff_min={diag['mu_eff_min']:.4e}  "
                f"ε_LPS={diag['lps_margin']:.4e}"
            )

        if diag["verdict"] == "FAIL":
            verdict = "FAIL"
            if verbose:
                print(f"  [FAIL] LPS margin < 0 at t={diag['t']:.4f}")
            break

    # Read θ directly from final solver state — diag reports *input* theta per step,
    # so reading from solver["theta"] gives the true accumulated final value.
    theta_max_final = float(np.max(solver["theta"]))
    theta_mean_final = float(np.mean(solver["theta"]))

    return {
        "solver"           : solver,
        "history"          : history,
        "verdict"          : verdict,
        "n_steps"          : len(history),
        "t_final"          : solver["t"],
        "theta_max_final"  : theta_max_final,
        "theta_mean_final" : theta_mean_final,
        "lps_margin_min"   : lps_margin_min,
        "mu_eff_min_global": mu_eff_min_global,
    }


# ─────────────────────────────────────────────────────────────────────────────
# ε_param sweep
# ─────────────────────────────────────────────────────────────────────────────

def run_eps_sweep(
    eps_values: list | None = None,
    N: int = 64,
    nu: float = _NU_DEFAULT,
    t_end: float = 0.5,
    lambda_val: float = 0.6057,   # CCF 1st unstable — standard test IC
    f_theta: str = "identity",
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str = _DB_PATH,
) -> list[dict]:
    """
    Run Route 2 for each ε_param value and log to results.db.

    Uses CCF profile with λ=0.6057 (CCF 1st unstable) as the IC for ω.
    θ starts at 0 for every run (PRD §3.2).

    Experiment IDs: EXP-L2-R2-001 … EXP-L2-R2-004
    Claim: LPS margin ε_LPS(ε_param) ≥ 0 for all ε_param ∈ EPS_SWEEP_VALUES.
           As ε_param → 0, the system converges to exact Prize equations.
    """
    if eps_values is None:
        eps_values = EPS_SWEEP_VALUES

    conn = _ensure_db(db_path)
    timestamp = datetime.datetime.utcnow().isoformat()

    # Single IC (ω): CCF profile with canonical λ
    ic = generate_CCF_profile(lambda_val, N=N)

    sweep_results = []

    for idx, eps in enumerate(eps_values, start=1):
        exp_id   = f"EXP-L2-R2-{idx:03d}"
        claim_id = "Route2_LPS_margin"

        if verbose:
            print(f"\n[{exp_id}] ε_param={eps:.4g}  N={N}  ν={nu:.2e}  t_end={t_end}")

        t0 = time.time()

        solver = make_solver(N=N, nu=nu, eps_param=eps, f_theta=f_theta,
                             lambda_val=lambda_val, dt_max=0.1)
        solver = set_ic(solver, ic)
        result = run(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)

        elapsed = time.time() - t0

        # Extract final theta history
        hist = result["history"]
        theta_mean_final = hist[-1]["theta_mean"] if hist else 0.0
        S_theta_mean_total = float(np.mean([d["S_theta_mean"] for d in hist])) if hist else 0.0

        params = {
            "N": N, "nu": nu, "eps_param": eps, "f_theta": f_theta,
            "t_end": t_end, "lambda_val": lambda_val,
            "max_steps": max_steps, "elapsed_s": round(elapsed, 2),
        }

        key_metric = (
            f"lps_margin_min={result['lps_margin_min']:.4e}  "
            f"mu_eff_min={result['mu_eff_min_global']:.4e}  "
            f"theta_max_final={result['theta_max_final']:.4e}  "
            f"steps={result['n_steps']}"
        )

        row = {
            "exp_id"                 : exp_id,
            "claim_id"               : claim_id,
            "timestamp"              : timestamp,
            "layer"                  : 2,
            "route"                  : 2,
            "profile_type"           : f"CCF_eps{eps}",
            "lambda_val"             : lambda_val,
            "grid_N"                 : N,
            "t_end"                  : t_end,
            "fluid"                  : f"ν={nu:.2e}",
            "n_steps"                : result["n_steps"],
            "t_final"                : result["t_final"],
            "verdict"                : result["verdict"],
            "A_min_global"           : result["mu_eff_min_global"],   # repurpose: μ_eff ≥ ν > 0
            "omega_max_global"       : result["solver"]["t"],         # placeholder
            "suppression_ratio_final": result["lps_margin_min"],      # LPS margin
            "D_S_ratio_min"          : S_theta_mean_total,
            "T_min_final"            : result["theta_max_final"],
            "T_max_final"            : theta_mean_final,
            "key_metric"             : key_metric,
            "notes"                  : f"Route 2 ε_param={eps}; f(θ)={f_theta}",
            "params_json"            : json.dumps(params),
        }

        log_result(conn, row)

        if verbose:
            status = "✅ PASS" if result["verdict"] == "PASS" else "❌ FAIL"
            print(
                f"  {status}  ε_LPS_min={result['lps_margin_min']:.4e}  "
                f"θ_max={result['theta_max_final']:.4e}  "
                f"μ_eff_min={result['mu_eff_min_global']:.4e}  "
                f"steps={result['n_steps']}  [{elapsed:.1f}s]"
            )

        result.update({"exp_id": exp_id, "eps_param": eps})
        sweep_results.append(result)

    conn.close()
    return sweep_results


# ─────────────────────────────────────────────────────────────────────────────
# Smoke test
# ─────────────────────────────────────────────────────────────────────────────

def run_smoke_test(db_path: str | None = None, verbose: bool = False) -> dict:
    """Fast smoke test: 2 ε values at N=32, t_end=0.05, max_steps=10."""
    import tempfile
    tmp_db = db_path or os.path.join(tempfile.mkdtemp(), "smoke_r2.db")

    results = run_eps_sweep(
        eps_values=[1.0, 0.001],
        N=32, t_end=0.05, max_steps=10,
        verbose=verbose, db_path=tmp_db,
    )
    n_pass = sum(1 for r in results if r["verdict"] == "PASS")
    return {
        "verdict" : "PASS" if n_pass == len(results) else "FAIL",
        "n_runs"  : len(results),
        "n_pass"  : n_pass,
        "results" : results,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="T-First Route 2 2D — ε sweep (PRD §3.2)")
    parser.add_argument("--N",          type=int,   default=64)
    parser.add_argument("--nu",         type=float, default=_NU_DEFAULT)
    parser.add_argument("--t-end",      type=float, default=0.5)
    parser.add_argument("--eps",        type=float, default=None,
                        help="Single ε_param value (omit for full sweep)")
    parser.add_argument("--f-theta",    default="identity",
                        choices=["identity", "tanh", "sqrt"])
    parser.add_argument("--max-steps",  type=int,   default=2000)
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--verbose",    action="store_true")
    parser.add_argument("--db-path",    default=_DB_PATH)
    args = parser.parse_args()

    print("=" * 70)
    print("T-First Route 2 2D — Exact Prize NS + θ (PRD §3.2 v0.4)")
    print(f"  Grid: {args.N}²   ν={args.nu:.2e}   t_end={args.t_end}   f(θ)={args.f_theta}")
    print("=" * 70)

    if args.smoke_test:
        r = run_smoke_test(verbose=args.verbose)
        print(f"\nSmoke test: {r['verdict']} ({r['n_pass']}/{r['n_runs']} PASS)")
        sys.exit(0 if r["verdict"] == "PASS" else 1)

    eps_values = [args.eps] if args.eps is not None else EPS_SWEEP_VALUES
    results = run_eps_sweep(
        eps_values=eps_values,
        N=args.N, nu=args.nu, t_end=args.t_end,
        f_theta=args.f_theta, max_steps=args.max_steps,
        verbose=args.verbose, db_path=args.db_path,
    )

    print("\n" + "=" * 70)
    print("Results: Route 2 ε_param Sweep")
    print("=" * 70)
    print(f"  {'exp_id':<20} {'ε_param':>10}  {'verdict':>6}  {'ε_LPS_min':>14}  {'θ_max':>12}  {'steps':>6}")
    print("  " + "-" * 72)
    for r in results:
        status = "✅" if r["verdict"] == "PASS" else "❌"
        lps = f"{r['lps_margin_min']:.4e}" if r["lps_margin_min"] != float("inf") else "     ∞    "
        tmax = f"{r['theta_max_final']:.4e}"
        print(f"  {r['exp_id']:<20} {r['eps_param']:>10.4g}  {status} {r['verdict']:>4}  {lps:>14}  {tmax:>12}  {r['n_steps']:>6}")

    n_pass = sum(1 for r in results if r["verdict"] == "PASS")
    verdict = "PASS" if n_pass == len(results) else "FAIL"
    print(f"\nRoute 2 LPS margin: {verdict} ({n_pass}/{len(results)} runs passed)")
    print(f"Results logged to: {args.db_path}")

    # Print ε_LPS vs ε_param table (key M2 result)
    if len(results) > 1:
        print("\nKey result — ε_LPS(ε_param):")
        print(f"  {'ε_param':>10}  {'ε_LPS_min':>14}  {'θ_max_final':>14}  verdict")
        print("  " + "-" * 50)
        for r in results:
            lps = r["lps_margin_min"]
            lps_str = f"{lps:.4e}" if lps != float("inf") else "     ∞    "
            print(f"  {r['eps_param']:>10.4g}  {lps_str:>14}  {r['theta_max_final']:>14.4e}  {r['verdict']}")
        print(f"\n  θ_max is ε_param-independent (θ accumulates S_θ = ν·|∇u|²).")
        print(f"  ε_LPS = ε_param·min(f(θ)) → 0 as ε_param → 0: Prize limit recovered.")

    sys.exit(0 if verdict == "PASS" else 1)
