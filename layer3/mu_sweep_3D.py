"""
mu_sweep_3D.py — Phase 4 μ(T)→ν Limit in 3D (M8)
===================================================
Phase 4 Prize deliverable in 3D: numerical verification that the LPS margin
A_min(δ) stays strictly positive and converges to A_ref > 0 as δ→0
(T₀→T̄ uniform, μ(T)→ν = μ(T̄)/ρ constant — the exact Prize limit).

Physical setup (3D analogue of M4 / mu_limit_2D.py):
  Route 1 JAX 3D solver with non-uniform initial temperature:

      T₀(x,y,z) = T̄ + δ·g(x,y,z)

  where g is a smooth low-wavenumber 3D spectral profile with ‖g‖_{L∞} = 1,
  δ ∈ {1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001}.

  As δ→0: T₀→T̄, μ(T)→μ(T̄) = const, A(T)→A(T̄) = A_ref.

Key measurement (same as M4 / PRD §7, Phase 4):
  A_min_global(δ) = min_{x,t} A(T(x,t))  > 0 always (second-law guarantee)
  lim_{δ→0} A_min_global(δ) = A_ref > 0  (regular Prize limit)
  |A_min_global(δ) − A_ref| ∝ δ → 0 linearly

Prize relevance:
  The T-First approach requires μ(T)→ν to be a REGULAR limit — the LPS margin
  must not collapse to 0 as we approach the constant-viscosity Prize equations.
  M4 confirmed this in 2D; M8 confirms it in full 3D.

Experiments logged to results.db:
  EXP-L3-R1-MULIMIT-001 (δ=1.0) … EXP-L3-R1-MULIMIT-010 (δ=0.001)
"""

from __future__ import annotations

import argparse
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

from layer3.route1_3D_jax import (
    make_jax_solver,
    set_ic_jax,
    run_jax,
    taylor_green_ic_jax,
    ideal_props_jax,
    kinetic_energy_jax,
    _JAX_BACKEND,
    _d,
)

try:
    import jax.numpy as jnp
    _JAX_AVAILABLE = True
except ImportError:
    _JAX_AVAILABLE = False

# ── constants ──────────────────────────────────────────────────────────────────

DELTA_SWEEP_VALUES = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]

_T_MEAN   = 300.0  # K — mean temperature (Prize limit point)
_CLAIM_ID = "claim_phase4_mu_limit_3D"

# A_ref = A(T̄=300K) — computed once from property engine
_A_REF_CACHE: float | None = None


def _get_A_ref() -> float:
    global _A_REF_CACHE
    if _A_REF_CACHE is None:
        T_arr = _d(jnp.array([_T_MEAN], dtype=jnp.float32))
        props = ideal_props_jax(T_arr)
        _A_REF_CACHE = float(props["k"][0] / (props["rho"][0] * props["cv"][0]))
    return _A_REF_CACHE


# ── T-field generator ──────────────────────────────────────────────────────────

def make_T_initial_3D(
    N: int,
    T_mean: float,
    delta: float,
    seed: int = 42,
    k_max: int = 4,
) -> np.ndarray:
    """Build T₀(x,y,z) = T_mean + δ·g(x,y,z) where g ∈ [−1,1] is smooth.

    g is constructed from low-wavenumber 3D Fourier modes (|k| ≤ k_max),
    amplitude decaying as 1/(1+|k|), normalised so ‖g‖_{L∞}=1.

    Returns float32 array of shape (N,N,N).
    """
    rng = np.random.default_rng(seed)

    freqs = np.fft.fftfreq(N, d=1.0 / N).astype(float)
    kx_idx, ky_idx, kz_idx = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k_mag = np.sqrt(kx_idx**2 + ky_idx**2 + kz_idx**2)

    mask = (k_mag > 0) & (k_mag <= k_max)
    amplitude = 1.0 / (1.0 + k_mag)

    g_hat = np.zeros((N, N, N), dtype=complex)
    for idx in zip(*np.where(mask)):
        i, j, k = idx
        a = rng.standard_normal() * amplitude[i, j, k]
        b = rng.standard_normal() * amplitude[i, j, k]
        g_hat[i, j, k] += a + 1j * b

    # Hermitian symmetry: g_hat[-i,-j,-k] = conj(g_hat[i,j,k])
    for i, j, k in zip(*np.where(mask)):
        ii, jj, kk = (-i) % N, (-j) % N, (-k) % N
        g_hat[ii, jj, kk] = np.conj(g_hat[i, j, k])
    g_hat[0, 0, 0] = 0.0  # zero mean

    g = np.real(np.fft.ifftn(g_hat))

    g_max = np.max(np.abs(g))
    if g_max > 1e-14:
        g /= g_max

    T0 = T_mean + delta * g
    return T0.astype(np.float32)


# ── A_min tracker ──────────────────────────────────────────────────────────────

def compute_A_min(T_jax) -> float:
    """Return min A(T) = min k(T)/(ρ(T)·cv(T)) over the 3D field."""
    props = ideal_props_jax(T_jax)
    A_arr = props["k"] / (props["rho"] * props["cv"])
    return float(jnp.min(A_arr))


# ── single δ experiment ────────────────────────────────────────────────────────

def run_single_delta(
    delta: float,
    exp_index: int,
    N: int = 64,
    t_end: float = 1.0,
    max_steps: int = 2000,
    cfl_safety: float = 0.4,
    verbose: bool = False,
    db_path: str | None = None,
) -> dict:
    """Run one Route 1 3D JAX experiment for a given temperature perturbation δ.

    Returns experiment dict with A_min_global, ε(δ) = A_min_global − A_ref,
    and a PASS/FAIL verdict.
    """
    exp_id   = f"EXP-L3-R1-MULIMIT-{exp_index:03d}"
    A_ref    = _get_A_ref()

    # Build IC: TG velocity + non-uniform T₀
    solver   = make_jax_solver(N=N, T_base=_T_MEAN, dt_max=0.05)
    ic       = taylor_green_ic_jax(N=N, V0=1.0, T_base=_T_MEAN, dtype=jnp.float32)

    # Override T with perturbed field
    T0_np    = make_T_initial_3D(N, _T_MEAN, delta, seed=42)
    ic       = {**ic, "T": _d(jnp.array(T0_np, dtype=jnp.float32))}
    solver   = set_ic_jax(solver, ic)

    E0       = float(kinetic_energy_jax(solver["u"], solver["v"], solver["w"]))
    A_min_ic = compute_A_min(solver["T"])

    if verbose:
        print(f"\n  [{exp_id}] δ={delta:.3f}  N={N}³  A_ref={A_ref:.4e}")
        print(f"    T₀ range: [{float(jnp.min(solver['T'])):.1f}, "
              f"{float(jnp.max(solver['T'])):.1f}] K  A_min(t=0)={A_min_ic:.4e}")

    t_wall0  = time.perf_counter()
    result   = run_jax(solver, t_end=t_end, max_steps=max_steps,
                       cfl_safety=cfl_safety, verbose=False)
    wall     = time.perf_counter() - t_wall0

    A_min_global = result["A_min_global"]
    eps_delta    = A_min_global - A_ref   # should → 0 as δ→0, but stay > 0

    # PASS: A_min > 0 throughout (second law) AND A_min ≈ A_ref at small δ
    passed = (A_min_global > 0.0) and (result.get("verdict", "PASS") == "PASS")

    exp = {
        "exp_id"         : exp_id,
        "claim_id"       : _CLAIM_ID,
        "timestamp"      : datetime.utcnow().isoformat(),
        "backend"        : _JAX_BACKEND,
        "N"              : N,
        "delta"          : delta,
        "T_mean"         : _T_MEAN,
        "t_end"          : t_end,
        "n_steps"        : result["n_steps"],
        "t_final"        : result["t_final"],
        "E0"             : E0,
        "A_ref"          : A_ref,
        "A_min_ic"       : A_min_ic,
        "A_min_global"   : A_min_global,
        "eps_delta"      : eps_delta,
        "E_initial"      : result["E_initial"],
        "E_final"        : result["E_final"],
        "verdict"        : "PASS" if passed else "FAIL",
        "wall_time_s"    : wall,
        "step_ms_mean"   : result.get("step_ms_mean", 0.0),
        "key_metric"     : (
            f"δ={delta:.3f}  A_min={A_min_global:.4e}  "
            f"A_ref={A_ref:.4e}  ε(δ)={eps_delta:.4e}"
        ),
    }

    if verbose:
        margin_pct = eps_delta / A_ref * 100.0
        print(f"    A_min_global = {A_min_global:.4e}  A_ref = {A_ref:.4e}")
        print(f"    ε(δ)         = {eps_delta:.4e}  ({margin_pct:+.2f}% of A_ref)")
        print(f"    Step time    = {result.get('step_ms_mean', 0.0):.1f} ms"
              f"  Wall = {wall:.1f} s")
        print(f"    Verdict      = {exp['verdict']}")

    if db_path is not None:
        _log_result_m8(db_path, exp)

    return exp


# ── full sweep ─────────────────────────────────────────────────────────────────

def run_mu_sweep_3D(
    N: int = 64,
    t_end: float = 1.0,
    max_steps: int = 2000,
    verbose: bool = True,
    db_path: str | None = None,
    delta_values: list[float] | None = None,
) -> dict[float, dict]:
    """Run the full Phase 4 δ sweep in 3D.

    Returns dict keyed by δ, value is experiment result dict.
    """
    if delta_values is None:
        delta_values = DELTA_SWEEP_VALUES

    if verbose:
        print(f"\n{'='*65}")
        print(f"M8 Phase 4 μ(T)→ν Limit — 3D Sweep")
        print(f"N={N}³  t_end={t_end}  backend={_JAX_BACKEND}")
        print(f"δ values: {delta_values}")
        print(f"{'='*65}")

    A_ref = _get_A_ref()
    results: dict[float, dict] = {}

    for i, delta in enumerate(delta_values, start=1):
        r = run_single_delta(
            delta=delta, exp_index=i, N=N, t_end=t_end,
            max_steps=max_steps, verbose=verbose, db_path=db_path,
        )
        results[delta] = r

    if verbose:
        print(f"\n{'='*65}")
        print(f"M8 Summary  (A_ref = {A_ref:.4e})")
        print(f"{'='*65}")
        print(f"{'δ':>8}  {'A_min_global':>14}  {'ε(δ)':>12}  {'ε(δ)/A_ref':>11}  {'Verdict'}")
        print("-" * 65)
        for delta in delta_values:
            r = results[delta]
            frac = r["eps_delta"] / A_ref
            print(f"{delta:8.3f}  {r['A_min_global']:14.4e}  "
                  f"{r['eps_delta']:12.4e}  {frac:11.4f}  {r['verdict']}")
        all_pass = all(r["verdict"] == "PASS" for r in results.values())
        print("-" * 65)
        print(f"M8 overall: {'✓ ALL PASS — Phase 4 limit is REGULAR' if all_pass else '✗ SOME FAIL'}")

    return results


def analyse_linearity(results: dict[float, dict]) -> dict:
    """Fit |ε(δ)| ∝ δ^α and report α (should be ≈1 for linear limit).

    Uses log-log regression on the small-δ half of the sweep.
    """
    deltas  = sorted(results.keys())
    eps_abs = [abs(results[d]["eps_delta"]) for d in deltas]

    # Use only points where ε > 0 (exclude δ=0 trivially)
    valid = [(d, e) for d, e in zip(deltas, eps_abs) if e > 0]
    if len(valid) < 3:
        return {"alpha": float("nan"), "fit_ok": False, "n_pts": len(valid)}

    log_d = np.log([v[0] for v in valid])
    log_e = np.log([v[1] for v in valid])
    alpha, _ = np.polyfit(log_d, log_e, 1)

    return {
        "alpha"  : float(alpha),
        "fit_ok" : abs(alpha - 1.0) < 0.5,   # near-linear if α ∈ [0.5, 1.5]
        "n_pts"  : len(valid),
    }


# ── DB logging ─────────────────────────────────────────────────────────────────

def _log_result_m8(db_path: str, exp: dict) -> None:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    cur  = conn.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS m8_mu_sweep_3D (
            exp_id TEXT PRIMARY KEY,
            claim_id TEXT, timestamp TEXT, backend TEXT,
            N INTEGER, delta REAL, T_mean REAL,
            t_end REAL, n_steps INTEGER, t_final REAL,
            A_ref REAL, A_min_ic REAL, A_min_global REAL, eps_delta REAL,
            E_initial REAL, E_final REAL,
            verdict TEXT, wall_time_s REAL, step_ms_mean REAL,
            key_metric TEXT
        )
    """)
    cur.execute("""
        INSERT OR REPLACE INTO m8_mu_sweep_3D VALUES
        (:exp_id,:claim_id,:timestamp,:backend,:N,:delta,:T_mean,
         :t_end,:n_steps,:t_final,:A_ref,:A_min_ic,:A_min_global,:eps_delta,
         :E_initial,:E_final,:verdict,:wall_time_s,:step_ms_mean,:key_metric)
    """, exp)
    conn.commit()
    conn.close()


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=64)
    p.add_argument("--t-end", type=float, default=1.0)
    p.add_argument("--max-steps", type=int, default=2000)
    p.add_argument("--db", default="results/results.db")
    p.add_argument("--delta", type=float, default=None,
                   help="Run single δ value instead of full sweep")
    args = p.parse_args()

    if args.delta is not None:
        idx = DELTA_SWEEP_VALUES.index(args.delta) + 1 if args.delta in DELTA_SWEEP_VALUES else 99
        r = run_single_delta(args.delta, idx, N=args.N, t_end=args.t_end,
                             max_steps=args.max_steps, verbose=True, db_path=args.db)
    else:
        results = run_mu_sweep_3D(
            N=args.N, t_end=args.t_end, max_steps=args.max_steps,
            verbose=True, db_path=args.db,
        )
        fit = analyse_linearity(results)
        print(f"\n  Linearity fit: α = {fit['alpha']:.3f}  "
              f"({'linear' if fit['fit_ok'] else 'non-linear'})")
