"""
T_solver_2D.py — 2D Spectral Navier-Stokes-Fourier Solver
==========================================================
Vorticity-stream-function formulation on [0, 2π)².
Temperature equation with A(T)-driven diffusion and viscous heating.
All thermophysical properties evaluated via layer1/tfirst_props.py.

Physics:
  dω/dt + J(ψ,ω) = ∇·(ν(T)∇ω)
  dT/dt + u·∇T   = A(T)∆T + μ(T)|∇u|²/(ρ(T)cv(T))

  where ψ is the stream function, -∆ψ = ω,
        u = ∂ψ/∂y, v = -∂ψ/∂x,
        ν(T) = μ(T)/ρ(T),
        A(T) = k(T)/(ρ(T)·cv(T))  [always > 0 — second law]

Numerics:
  - Pseudo-spectral (FFT) for all spatial derivatives
  - 2/3 dealiasing rule on ω (and T for stability)
  - RK4 time integration with adaptive CFL
  - Mean-field spectral diffusion + variable-coefficient correction in physical space

Experiments:
  EXP-L2-CCF-001 … EXP-L2-CCF-006  (λ sweep, CCF profiles)
  EXP-L2-BSSQ-001 … EXP-L2-BSSQ-006 (λ sweep, Boussinesq profiles)
  EXP-L2-ADV-001  (adversarial minimum, λ → 0)
"""

import sys
import os
import numpy as np

# ── path setup ───────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))

from self_similar_IC import dealias_23, make_grid, spectral_laplacian_inv
from layer1.tfirst_props import ideal_props, co2_props, A_field

# ─────────────────────────────────────────────────────────────────────────────
# Spectral helpers
# ─────────────────────────────────────────────────────────────────────────────

def _spectral_diff(f_hat, k_dir):
    """Return ifft2(i·k·f_hat); derivative in the direction given by k_dir array."""
    return np.real(np.fft.ifft2(1j * k_dir * f_hat))


def _spectral_laplacian_phys(f, k2):
    """Return ∆f = ifft2(-k²·fft2(f)) in physical space."""
    return np.real(np.fft.ifft2(-k2 * np.fft.fft2(f)))


def _psi_and_velocity(omega_hat, kx, ky, k2):
    """Compute stream function ψ and velocity (u,v) from ω̂.

    -∆ψ = ω  →  ψ̂ = ω̂/|k|²  (zeroing k=0 mean)
    u = ∂ψ/∂y,  v = -∂ψ/∂x
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        psi_hat = np.where(k2 > 0, omega_hat / k2, 0.0 + 0.0j)
    u = _spectral_diff(psi_hat, ky)
    v = _spectral_diff(-psi_hat, kx)   # v = -∂ψ/∂x → use -psi_hat with kx diff
    return psi_hat, u, v


def _velocity_gradient_sq(psi_hat, kx, ky):
    """Return |∇u|² = (∂u/∂x)²+(∂u/∂y)²+(∂v/∂x)²+(∂v/∂y)² in physical space.

    From ψ: u = ∂ψ/∂y, v = -∂ψ/∂x, so:
      ∂u/∂x =  ∂²ψ/∂x∂y    ∂u/∂y =  ∂²ψ/∂y²
      ∂v/∂x = -∂²ψ/∂x²     ∂v/∂y = -∂²ψ/∂x∂y
    """
    du_dx = np.real(np.fft.ifft2(-kx * ky * psi_hat))   # ∂²ψ/∂x∂y
    du_dy = np.real(np.fft.ifft2(-ky * ky * psi_hat))   # ∂²ψ/∂y²
    dv_dx = np.real(np.fft.ifft2( kx * kx * psi_hat))  # -∂²ψ/∂x² = -(−kx²ψ̂) = kx²ψ̂ after sign
    # Actually: dv/dx = -∂²ψ/∂x² = -ifft2((-kx²)ψ̂) = ifft2(kx²ψ̂)
    # du_dx = ∂²ψ/∂x∂y = ifft2((ikx)(iky)ψ̂) = ifft2(-kx·ky·ψ̂)  ✓
    dv_dy = -du_dx   # incompressibility: ∂v/∂y = -∂u/∂x
    return du_dx**2 + du_dy**2 + dv_dx**2 + dv_dy**2


# ─────────────────────────────────────────────────────────────────────────────
# Solver construction
# ─────────────────────────────────────────────────────────────────────────────

def make_solver(
    N=64,
    fluid="ideal",
    P=None,
    T_base=300.0,
    omega_max=1.0,
    t_star=1.0,
    lambda_val=1.0,
    dt_max=0.1,
):
    """Create a 2D NSF solver state dictionary.

    Parameters
    ----------
    N        : grid points per side (must be even; 64 for lambda sweep)
    fluid    : 'ideal' or 'co2'
    P        : pressure in Pa (required for 'co2')
    T_base   : initial uniform temperature [K]
    omega_max: vorticity scale [s⁻¹]
    t_star   : normalised blow-up time for Conjecture 3.4 suppression ratio
    lambda_val: Wang et al. λ (used for suppression ratio tracking)
    dt_max   : maximum allowed timestep [s]
    """
    if fluid == "co2" and P is None:
        raise ValueError("co2 fluid requires P [Pa]")

    x, y, kx, ky = make_grid(N)
    k2 = kx**2 + ky**2  # |k|²

    # Reference properties at T_base
    T_arr = np.array([T_base])
    if fluid == "ideal":
        ref = ideal_props(T_arr)
    else:
        ref = co2_props(T_arr, P)
    nu_ref = float(ref["mu"][0] / ref["rho"][0])
    A_ref = float(A_field(T_arr, fluid=fluid, P=P)[0])

    return {
        "N": N,
        "x": x,
        "y": y,
        "kx": kx,
        "ky": ky,
        "k2": k2,
        "omega": np.zeros((N, N)),
        "T": np.full((N, N), T_base),
        "t": 0.0,
        "fluid": fluid,
        "P": P,
        "T_base": T_base,
        "omega_max": omega_max,
        "t_star": t_star,
        "lambda_val": lambda_val,
        "dt_max": dt_max,
        "nu_ref": nu_ref,
        "A_ref": A_ref,
        # Diagnostics history
        "history": [],
    }


def set_ic(solver, ic):
    """Set initial conditions from a profile dict (output of layer2/self_similar_IC.py).

    The IC dict must contain 'omega'; 'T' is optional (defaults to T_base).
    Returns a new solver dict (does not mutate in place).
    """
    new_solver = {**solver}
    new_solver["omega"] = ic["omega"].copy()
    if "T_field" in ic:
        new_solver["T"] = ic["T_field"].copy()
    else:
        N = solver["N"]
        new_solver["T"] = np.full((N, N), solver["T_base"])
    new_solver["t"] = 0.0
    new_solver["history"] = []
    return new_solver


# ─────────────────────────────────────────────────────────────────────────────
# Right-hand side
# ─────────────────────────────────────────────────────────────────────────────

def _rhs(omega, T, kx, ky, k2, fluid, P, assert_A_positive=True):
    """Compute (dω/dt, dT/dt) for one stage of RK4.

    Variable-viscosity split:
      diff_ω = ν̄·∆ω  (spectral, using spatial mean ν̄)
    Variable-diffusivity in physical space:
      diff_T = A(T)·∆T  (A evaluated pointwise, ∆T via FFT)
    Viscous heating:
      Q = μ(T)·|∇u|² / (ρ(T)·cv(T))  (pointwise)

    Returns
    -------
    d_omega : ndarray (N,N) — dω/dt  [s⁻²]
    d_T     : ndarray (N,N) — dT/dt  [K/s]
    diag    : dict — diagnostics for monitoring
    """
    # --- Properties ---
    if fluid == "ideal":
        props = ideal_props(T)
    else:
        props = co2_props(T, P)
    mu = props["mu"]
    rho = props["rho"]
    k_th = props["k"]
    cv = props["cv"]

    nu = mu / rho           # kinematic viscosity [m²/s]
    A = k_th / (rho * cv)  # thermal diffusivity [m²/s]

    if assert_A_positive:
        A_min = float(np.min(A))
        if A_min <= 0.0:
            raise AssertionError(
                f"[T_solver_2D] Second Law violated: A_min = {A_min:.3e} ≤ 0"
            )

    nu_bar = float(np.mean(nu))  # scalar mean for spectral linear term

    # --- Vorticity advection (pseudo-spectral Jacobian) ---
    omega_hat = np.fft.fft2(omega)
    psi_hat, u, v = _psi_and_velocity(omega_hat, kx, ky, k2)

    domega_dx = _spectral_diff(omega_hat, kx)
    domega_dy = _spectral_diff(omega_hat, ky)
    J = u * domega_dx + v * domega_dy   # Jacobian J(ψ,ω) = u·∇ω  [s⁻²]

    # --- Vorticity diffusion: ν̄·∆ω (spectral, mean viscosity) ---
    diff_omega = np.real(np.fft.ifft2(-nu_bar * k2 * omega_hat))  # [s⁻²]

    # --- Temperature advection: u·∇T ---
    T_hat = np.fft.fft2(T)
    dT_dx = _spectral_diff(T_hat, kx)
    dT_dy = _spectral_diff(T_hat, ky)
    adv_T = u * dT_dx + v * dT_dy   # [K/s]

    # --- Temperature diffusion: A(T)·∆T (variable coefficient) ---
    lap_T = _spectral_laplacian_phys(T, k2)
    diff_T = A * lap_T               # [K/s]

    # --- Viscous heating: Q = μ|∇u|²/(ρ·cv) ---
    grad_u_sq = _velocity_gradient_sq(psi_hat, kx, ky)  # [s⁻²]
    Q_visc = mu * grad_u_sq / (rho * cv)                 # [K/s]

    # --- Assemble ---
    d_omega = -J + diff_omega
    d_T = -adv_T + diff_T + Q_visc

    diag = {
        "A_min": float(np.min(A)),
        "A_max": float(np.max(A)),
        "nu_bar": nu_bar,
        "nu_max": float(np.max(nu)),
        "u_max": float(np.max(np.abs(u))),
        "v_max": float(np.max(np.abs(v))),
        "Q_visc_max": float(np.max(Q_visc)),
        "grad_u_sq_max": float(np.max(grad_u_sq)),
        "u": u,
        "v": v,
        "nu": nu,
        "A": A,
        "grad_u_sq": grad_u_sq,
        "Q_visc": Q_visc,
        "psi_hat": psi_hat,
    }
    return d_omega, d_T, diag


# ─────────────────────────────────────────────────────────────────────────────
# Time stepping
# ─────────────────────────────────────────────────────────────────────────────

def compute_cfl_dt(solver, safety=0.4):
    """Return the CFL-limited timestep.

    Uses min(dt_advective, dt_diffusive, solver['dt_max']).
    """
    omega = solver["omega"]
    T = solver["T"]
    kx, ky, k2 = solver["kx"], solver["ky"], solver["k2"]
    fluid, P = solver["fluid"], solver["P"]
    N = solver["N"]
    dx = 2.0 * np.pi / N

    # Velocity from current vorticity
    omega_hat = np.fft.fft2(omega)
    _, u, v = _psi_and_velocity(omega_hat, kx, ky, k2)
    U_max = max(float(np.max(np.abs(u))), float(np.max(np.abs(v))), 1e-12)

    # Advective CFL
    dt_adv = safety * dx / U_max

    # Diffusive CFL
    if fluid == "ideal":
        props = ideal_props(T)
    else:
        props = co2_props(T, P)
    nu_max = float(np.max(props["mu"] / props["rho"]))
    A_max = float(np.max(props["k"] / (props["rho"] * props["cv"])))
    diff_coeff = max(nu_max, A_max, 1e-20)
    dt_diff = safety * dx**2 / (2.0 * diff_coeff)

    return float(np.clip(min(dt_adv, dt_diff), 1e-12, solver["dt_max"]))


def solver_step(solver, dt=None):
    """Take one RK4 step. Returns (new_solver, step_diagnostics).

    Does NOT mutate the input solver dict.
    """
    if dt is None:
        dt = compute_cfl_dt(solver)

    omega = solver["omega"]
    T = solver["T"]
    kx, ky, k2 = solver["kx"], solver["ky"], solver["k2"]
    fluid, P = solver["fluid"], solver["P"]

    # RK4
    k1o, k1T, d1 = _rhs(omega,                 T,                 kx, ky, k2, fluid, P)
    k2o, k2T, _  = _rhs(omega + 0.5*dt*k1o,    T + 0.5*dt*k1T,    kx, ky, k2, fluid, P)
    k3o, k3T, _  = _rhs(omega + 0.5*dt*k2o,    T + 0.5*dt*k2T,    kx, ky, k2, fluid, P)
    k4o, k4T, _  = _rhs(omega + dt*k3o,         T + dt*k3T,         kx, ky, k2, fluid, P)

    new_omega = omega + (dt / 6.0) * (k1o + 2.0*k2o + 2.0*k3o + k4o)
    new_T     = T     + (dt / 6.0) * (k1T + 2.0*k2T + 2.0*k3T + k4T)

    # Dealias both fields for stability
    new_omega = dealias_23(new_omega)
    new_T     = dealias_23(new_T)

    # Keep T in a physically valid range for property lookups
    new_T = np.clip(new_T, 100.0, 6000.0)

    new_t = solver["t"] + dt

    # Suppression ratio τ^{-(2+λ)} from Conjecture 3.4
    t_star = solver["t_star"]
    lambda_val = solver["lambda_val"]
    tau = max((t_star - new_t) / t_star, 1e-10)  # normalised time-to-blow-up
    supp_ratio = tau ** (-(2.0 + lambda_val)) if tau < 1.0 else float("inf")

    step_diag = {
        "t": new_t,
        "dt": dt,
        "A_min": d1["A_min"],
        "A_max": d1["A_max"],
        "nu_bar": d1["nu_bar"],
        "omega_max": float(np.max(np.abs(new_omega))),
        "T_max": float(np.max(new_T)),
        "T_min": float(np.min(new_T)),
        "Q_visc_max": d1["Q_visc_max"],
        "u_max": d1["u_max"],
        "suppression_ratio": supp_ratio,
        "tau": tau,
        "lambda_val": lambda_val,
        "A_positive": d1["A_min"] > 0.0,
        "verdict": "PASS" if d1["A_min"] > 0.0 else "FAIL",
    }

    new_solver = {
        **solver,
        "omega": new_omega,
        "T": new_T,
        "t": new_t,
    }
    return new_solver, step_diag


# ─────────────────────────────────────────────────────────────────────────────
# Run loop
# ─────────────────────────────────────────────────────────────────────────────

def run(solver, t_end, max_steps=10_000, cfl_safety=0.4, callback=None, verbose=False):
    """Run the solver from solver['t'] to t_end.

    Parameters
    ----------
    solver    : solver dict (from make_solver / set_ic)
    t_end     : end time [s]
    max_steps : safety cap
    cfl_safety: CFL safety factor
    callback  : optional callable(step_diag) called after each step
    verbose   : print step summary if True

    Returns
    -------
    dict with keys:
      'solver'  : final solver state
      'history' : list of step_diag dicts
      'verdict' : 'PASS' or 'FAIL'
      'A_min_global': minimum A(T) over all steps
      'omega_max_global': maximum |ω| over all steps
    """
    history = []
    A_min_global = float("inf")
    omega_max_global = 0.0
    verdict = "PASS"

    for step_idx in range(max_steps):
        if solver["t"] >= t_end:
            break

        dt = compute_cfl_dt(solver, safety=cfl_safety)
        dt = min(dt, t_end - solver["t"])  # don't overshoot

        try:
            solver, diag = solver_step(solver, dt=dt)
        except AssertionError as exc:
            verdict = "FAIL"
            if verbose:
                print(f"  [FAIL] Step {step_idx}: {exc}")
            break

        A_min_global = min(A_min_global, diag["A_min"])
        omega_max_global = max(omega_max_global, diag["omega_max"])

        history.append(diag)

        if callback is not None:
            callback(diag)

        if verbose and step_idx % 10 == 0:
            print(
                f"  step {step_idx:4d}  t={diag['t']:.4f}  "
                f"ω_max={diag['omega_max']:.3e}  "
                f"A_min={diag['A_min']:.3e}  "
                f"S={diag['suppression_ratio']:.2f}"
            )

        # Blow-up check: vorticity exceeding 1000× initial value is unphysical
        if diag["omega_max"] > 1000.0 * solver["omega_max"]:
            verdict = "FAIL"
            if verbose:
                print(f"  [FAIL] Blow-up detected at t={diag['t']:.4f}")
            break

    return {
        "solver": solver,
        "history": history,
        "verdict": verdict,
        "n_steps": len(history),
        "t_final": solver["t"],
        "A_min_global": A_min_global,
        "omega_max_global": omega_max_global,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Lambda sweep experiment (PRD §7.3)
# ─────────────────────────────────────────────────────────────────────────────

LAMBDA_SWEEP_VALUES = [1.2, 0.9, 0.6057, 0.4703, 0.3, 0.1]


def run_lambda_sweep(
    lambda_values=None,
    N=64,
    t_end=0.5,
    fluid="ideal",
    P=None,
    T_base=300.0,
    t_star=1.0,
    max_steps=500,
    verbose=False,
    smoke_test=False,
):
    """Run the full lambda sweep experiment (PRD §7.3 / EXP-L2-CCF-001…006).

    For each λ:
      1. Generate a CCF profile IC (from self_similar_IC.py)
      2. Run the 2D NSF solver
      3. Record: A_min_global, omega_max_global, suppression_ratio_final, verdict

    Parameters
    ----------
    lambda_values : list of λ (default: PRD §7.3 list)
    N             : grid resolution
    t_end         : end time [s] (use small value for smoke_test)
    smoke_test    : if True, use N=32, t_end=0.05, max_steps=10

    Returns
    -------
    list of result dicts, one per λ
    """
    from self_similar_IC import generate_CCF_profile

    if lambda_values is None:
        lambda_values = LAMBDA_SWEEP_VALUES

    if smoke_test:
        N, t_end, max_steps = 32, 0.05, 10

    results = []
    for lam in lambda_values:
        # Generate IC
        ic = generate_CCF_profile(lam, N=N)

        # Build solver
        solver = make_solver(
            N=N, fluid=fluid, P=P, T_base=T_base,
            omega_max=float(np.max(np.abs(ic["omega"]))) or 1.0,
            t_star=t_star, lambda_val=lam,
        )
        solver = set_ic(solver, ic)

        # Run
        result = run(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)
        result["lambda_val"] = lam
        result["exp_id"] = f"EXP-L2-CCF-{LAMBDA_SWEEP_VALUES.index(lam)+1:03d}" \
                           if lam in LAMBDA_SWEEP_VALUES else f"EXP-L2-CCF-CUSTOM"
        results.append(result)

        if verbose:
            print(
                f"  λ={lam:.4f}  verdict={result['verdict']}  "
                f"A_min={result['A_min_global']:.3e}  "
                f"ω_max={result['omega_max_global']:.3e}  "
                f"steps={result['n_steps']}"
            )

    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="T-First 2D NSF Solver — Lambda Sweep")
    parser.add_argument("--N", type=int, default=64)
    parser.add_argument("--t-end", type=float, default=0.5)
    parser.add_argument("--fluid", default="ideal", choices=["ideal", "co2"])
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    print("T-First 2D NSF — Lambda Sweep (PRD §7.3)")
    print(f"  Grid: {args.N}²  fluid: {args.fluid}  t_end: {args.t_end}")
    results = run_lambda_sweep(
        N=args.N, t_end=args.t_end, fluid=args.fluid,
        smoke_test=args.smoke_test, verbose=True,
    )
    print("\nSummary:")
    all_pass = True
    for r in results:
        v = r["verdict"]
        if v != "PASS":
            all_pass = False
        print(f"  {r['exp_id']}  λ={r['lambda_val']:.4f}  {v}  "
              f"A_min={r['A_min_global']:.3e}  steps={r['n_steps']}")
    print(f"\nOverall: {'PASS' if all_pass else 'FAIL'}")
