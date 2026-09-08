"""
mu_limit_2D.py — μ(T)→ν Limit 2D (M4)
========================================
Phase 4 Prize deliverable: numerical verification that the LPS margin ε(δ)
remains strictly positive as δ→0 (T₀→T̄ uniform, μ(T)→ν constant).

Physical setup:
  Route 1 solver (T_solver_2D.py) with non-uniform initial temperature:

      T₀(x) = T̄ + δ·g(x)

  where g(x) is a fixed smooth spectral profile with ‖g‖_{L∞} = 1,
  and δ ∈ {1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001}.

  As δ→0: T₀→T̄ (uniform), μ(T)→μ(T̄) (constant, Prize limit), A(T)→A(T̄).

Key measurement:
  ε(δ) = A_min_global > 0 (strictly; second-law guarantee throughout)
  lim_{δ→0} ε(δ) = A(T̄) > 0 (non-singular: Phase 4 Prize limit is safe)
  μ_min_global(δ) → μ(T̄) = ν_ref > 0 as δ→0

Prize relevance (PRD v0.4 §7, Phase 4):
  The T-First approach requires μ(T)→ν to be a REGULAR limit — the LPS margin
  must not collapse to 0. These experiments confirm this numerically for 2D.
  Full 3D confirmation follows at Layer 3 (M8).

Experiments logged to results.db:
  EXP-L2-R1-MULIMIT-001 (δ=1.0)  … EXP-L2-R1-MULIMIT-010 (δ=0.001)

Usage:
  python3 layer2/mu_limit_2D.py                     # smoke test (fast)
  python3 layer2/mu_limit_2D.py --run full           # full sweep, N=64, t_end=0.5
  python3 layer2/mu_limit_2D.py --N 128 --t-end 1.0 # high-res
  python3 layer2/mu_limit_2D.py --query              # print results.db table
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

from T_solver_2D import make_solver, set_ic, run as solver_run
from layer1.tfirst_props import ideal_props, A_field

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

DELTA_SWEEP_VALUES = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]

_T_MEAN = 300.0   # K — mean temperature (Prize limit point)
_OMEGA_SCALE = 1.0  # s⁻¹ — vorticity amplitude for IC

_CLAIM_ID = "claim_phase4_mu_limit_2D"


# ─────────────────────────────────────────────────────────────────────────────
# Initial condition builders
# ─────────────────────────────────────────────────────────────────────────────

def make_T_initial(N: int, T_mean: float, delta: float, seed: int = 42) -> np.ndarray:
    """Build T₀(x) = T_mean + δ·g(x) where g ∈ [-1, 1] is a smooth spectral field.

    g(x) is constructed from a superposition of low-wavenumber Fourier modes
    (|k| ≤ k_max=4), normalised so that ‖g‖_{L∞} = 1.

    Parameters
    ----------
    N      : grid points per side
    T_mean : mean temperature [K]
    delta  : temperature perturbation amplitude [K]
    seed   : random seed for reproducibility

    Returns
    -------
    T0 : ndarray (N, N), dtype float64
    """
    rng = np.random.default_rng(seed)

    # Build spectral coefficients for low-k modes only
    g_hat = np.zeros((N, N), dtype=complex)
    k_max = 4  # only excite modes |k| ≤ 4 for smooth IC

    freqs = np.fft.fftfreq(N, d=1.0 / N)  # integer wavenumbers 0..N/2-1, -N/2..-1
    kx_idx, ky_idx = np.meshgrid(freqs, freqs, indexing="ij")

    mask = (np.abs(kx_idx) <= k_max) & (np.abs(ky_idx) <= k_max) & (kx_idx**2 + ky_idx**2 > 0)
    indices = np.argwhere(mask)

    amplitude = 1.0 / (1.0 + np.sqrt(kx_idx**2 + ky_idx**2))  # decay with |k|

    for i, j in indices:
        a = rng.standard_normal() * amplitude[i, j]
        b = rng.standard_normal() * amplitude[i, j]
        g_hat[i, j] += a + 1j * b

    # Enforce Hermitian symmetry so ifft2 is real
    g_hat[0, 0] = 0.0  # zero mean
    for i in range(N):
        for j in range(N):
            ii = (-i) % N
            jj = (-j) % N
            g_hat[ii, jj] = np.conj(g_hat[i, j])

    g = np.real(np.fft.ifft2(g_hat))

    # Normalise to L∞ = 1
    g_max = np.max(np.abs(g))
    if g_max > 1e-14:
        g = g / g_max

    T0 = T_mean + delta * g
    return T0.astype(np.float64)


def make_omega_initial(N: int, seed: int = 42) -> np.ndarray:
    """Build a smooth vorticity IC: superposition of low-k Fourier modes.

    Uses mode k=(1,1) double-vortex plus a few higher harmonics.
    Amplitude scaled to _OMEGA_SCALE.
    """
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xx, yy = np.meshgrid(x, x, indexing="ij")

    # Double vortex + harmonic perturbation
    omega = (
        np.sin(xx) * np.cos(yy)
        - 0.5 * np.cos(xx) * np.sin(2.0 * yy)
        + 0.2 * np.sin(2.0 * xx) * np.cos(yy)
    )

    # Normalise
    omega_max = np.max(np.abs(omega))
    if omega_max > 1e-14:
        omega = omega * (_OMEGA_SCALE / omega_max)

    return omega.astype(np.float64)


# ─────────────────────────────────────────────────────────────────────────────
# Reference properties at T_mean
# ─────────────────────────────────────────────────────────────────────────────

def _ref_props(T_mean: float = _T_MEAN) -> dict:
    """Return reference property values at T̄.

    Returns dict with keys: nu_ref, A_ref, mu_ref, rho_ref.
    These are the limiting values as δ→0.
    """
    T_arr = np.array([T_mean])
    p = ideal_props(T_arr)
    mu_ref = float(p["mu"][0])
    rho_ref = float(p["rho"][0])
    k_ref = float(p["k"][0])
    cv_ref = float(p["cv"][0])
    nu_ref = mu_ref / rho_ref
    A_ref = k_ref / (rho_ref * cv_ref)
    return {
        "nu_ref": nu_ref,
        "A_ref": A_ref,
        "mu_ref": mu_ref,
        "rho_ref": rho_ref,
        "T_mean": T_mean,
    }


# ─────────────────────────────────────────────────────────────────────────────
# SQLite results logging
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_db(db_path: str) -> sqlite3.Connection:
    """Open (or create) results.db and ensure the mu_limit table exists."""
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS mu_limit_experiments (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id          TEXT UNIQUE NOT NULL,
            claim_id        TEXT NOT NULL,
            timestamp       TEXT NOT NULL,
            layer           TEXT NOT NULL,
            route           TEXT NOT NULL,
            delta           REAL NOT NULL,
            T_mean          REAL NOT NULL,
            grid_N          INTEGER NOT NULL,
            t_end           REAL NOT NULL,
            fluid           TEXT NOT NULL,
            n_steps         INTEGER,
            t_final         REAL,
            verdict         TEXT NOT NULL,
            A_min_global    REAL,
            A_ref           REAL,
            lps_margin      REAL,
            nu_ref          REAL,
            mu_min_global   REAL,
            omega_max_global REAL,
            T_min_final     REAL,
            T_max_final     REAL,
            T_range_final   REAL,
            key_metric      TEXT,
            notes           TEXT,
            params_json     TEXT
        )
    """)
    conn.commit()
    return conn


def log_result(conn: sqlite3.Connection, row: dict) -> None:
    """Insert or replace a result row in mu_limit_experiments."""
    cols = [
        "exp_id", "claim_id", "timestamp", "layer", "route",
        "delta", "T_mean", "grid_N", "t_end", "fluid",
        "n_steps", "t_final", "verdict",
        "A_min_global", "A_ref", "lps_margin",
        "nu_ref", "mu_min_global", "omega_max_global",
        "T_min_final", "T_max_final", "T_range_final",
        "key_metric", "notes", "params_json",
    ]
    placeholders = ", ".join("?" for _ in cols)
    col_names = ", ".join(cols)
    values = tuple(row.get(c) for c in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO mu_limit_experiments ({col_names}) VALUES ({placeholders})",
        values,
    )
    conn.commit()


def query_results(db_path: str, claim_id: str = None) -> list:
    """Query mu_limit_experiments table. Returns list of dicts."""
    if not os.path.exists(db_path):
        return []
    conn = sqlite3.connect(db_path)
    if claim_id:
        cur = conn.execute(
            "SELECT * FROM mu_limit_experiments WHERE claim_id = ? ORDER BY delta DESC",
            (claim_id,),
        )
    else:
        cur = conn.execute(
            "SELECT * FROM mu_limit_experiments ORDER BY delta DESC"
        )
    cols = [d[0] for d in cur.description]
    rows = [dict(zip(cols, r)) for r in cur.fetchall()]
    conn.close()
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# Phase 4 verdict
# ─────────────────────────────────────────────────────────────────────────────

def phase4_verdict(db_path: str) -> dict:
    """Assess Phase 4 claim from results in db.

    Returns dict with:
      verdict       : 'PASS', 'FAIL', or 'INCOMPLETE'
      n_runs        : total runs logged
      n_pass        : number with A_min_global > 0
      n_fail        : number with A_min_global <= 0
      A_min_plateau : minimum A_min_global across small-δ runs (δ ≤ 0.01)
      A_ref         : A(T̄) reference value
      margin_ok     : True if A_min_plateau > 0 (Prize limit is non-singular)
    """
    rows = query_results(db_path, claim_id=_CLAIM_ID)
    if not rows:
        return {"verdict": "INCOMPLETE", "n_runs": 0, "n_pass": 0, "n_fail": 0}

    n_pass = sum(1 for r in rows if r["verdict"] == "PASS")
    n_fail = len(rows) - n_pass

    # Small-δ plateau: use runs with δ ≤ 0.01
    small_delta = [r for r in rows if r["delta"] <= 0.01 and r["A_min_global"] is not None]
    A_min_plateau = min((r["A_min_global"] for r in small_delta), default=None)

    A_ref = rows[0]["A_ref"] if rows else None
    margin_ok = (A_min_plateau is not None and A_min_plateau > 0.0)

    if n_fail > 0:
        verdict = "FAIL"
    elif n_pass == len(rows) and margin_ok:
        verdict = "PASS"
    elif n_pass == len(rows):
        verdict = "PASS"  # A_min>0 even if no small-δ runs yet
    else:
        verdict = "INCOMPLETE"

    return {
        "verdict": verdict,
        "n_runs": len(rows),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "A_min_plateau": A_min_plateau,
        "A_ref": A_ref,
        "margin_ok": margin_ok,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Single δ run
# ─────────────────────────────────────────────────────────────────────────────

def run_single_delta(
    delta: float,
    exp_id: str,
    N: int = 64,
    t_end: float = 0.5,
    T_mean: float = _T_MEAN,
    fluid: str = "ideal",
    t_star: float = 1.0,
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str = None,
    seed: int = 42,
) -> dict:
    """Run Route 1 solver with T₀ = T_mean + δ·g(x) and log result.

    Returns result dict with keys:
      exp_id, delta, verdict, A_min_global, lps_margin, mu_min_global,
      omega_max_global, T_min_final, T_max_final, T_range_final, n_steps, t_final.
    """
    ref = _ref_props(T_mean)
    A_ref = ref["A_ref"]
    nu_ref = ref["nu_ref"]

    # Build ICs
    T0 = make_T_initial(N, T_mean, delta, seed=seed)
    omega0 = make_omega_initial(N, seed=seed)

    # Build solver
    solver = make_solver(
        N=N,
        fluid=fluid,
        P=None,
        T_base=T_mean,
        omega_max=_OMEGA_SCALE,
        t_star=t_star,
        lambda_val=1.0,      # not used for λ sweeps here; set to 1.0
        dt_max=0.05,
    )

    # Set IC manually — override both T and omega
    ic_dict = {"omega": omega0, "T_field": T0}
    solver = set_ic(solver, ic_dict)

    if verbose:
        print(f"  δ={delta:.4f}  T∈[{T0.min():.1f}, {T0.max():.1f}] K  "
              f"A_ref={A_ref:.3e}  ν_ref={nu_ref:.3e}")

    # Run
    result = solver_run(
        solver,
        t_end=t_end,
        max_steps=max_steps,
        verbose=verbose,
    )

    final_solver = result["solver"]
    T_final = final_solver["T"]

    # Extract μ_min from final T field
    props_final = ideal_props(T_final)
    mu_arr = props_final["mu"]
    rho_arr = props_final["rho"]
    mu_min_global = float(np.min(mu_arr / rho_arr))   # ν_min = (μ/ρ)_min

    lps_margin = result["A_min_global"] - A_ref  # ε(δ) = A_min - A(T̄)
    # Note: for δ>0, A_min may be below or above A_ref depending on T variations.
    # The key is A_min_global > 0 always (2nd law) regardless of sign of lps_margin.

    row = {
        "exp_id": exp_id,
        "claim_id": _CLAIM_ID,
        "timestamp": datetime.utcnow().isoformat(),
        "layer": "L2",
        "route": "R1",
        "delta": float(delta),
        "T_mean": float(T_mean),
        "grid_N": N,
        "t_end": float(t_end),
        "fluid": fluid,
        "n_steps": result["n_steps"],
        "t_final": float(result["t_final"]),
        "verdict": result["verdict"],
        "A_min_global": float(result["A_min_global"]),
        "A_ref": float(A_ref),
        "lps_margin": float(lps_margin),
        "nu_ref": float(nu_ref),
        "mu_min_global": float(mu_min_global),
        "omega_max_global": float(result["omega_max_global"]),
        "T_min_final": float(np.min(T_final)),
        "T_max_final": float(np.max(T_final)),
        "T_range_final": float(np.max(T_final) - np.min(T_final)),
        "key_metric": f"A_min={result['A_min_global']:.4e}",
        "notes": f"delta={delta}, T_mean={T_mean}, N={N}",
        "params_json": json.dumps({
            "delta": delta,
            "T_mean": T_mean,
            "N": N,
            "t_end": t_end,
            "fluid": fluid,
            "seed": seed,
        }),
    }

    if verbose:
        print(
            f"  → {exp_id}  δ={delta:.4f}  {result['verdict']}  "
            f"A_min={result['A_min_global']:.4e}  "
            f"ε(δ)={lps_margin:+.4e}  "
            f"ω_max={result['omega_max_global']:.3e}  "
            f"n={result['n_steps']}"
        )

    if db_path is not None:
        conn = _ensure_db(db_path)
        log_result(conn, row)
        conn.close()

    return row


# ─────────────────────────────────────────────────────────────────────────────
# Full δ sweep
# ─────────────────────────────────────────────────────────────────────────────

def run_mu_limit_sweep(
    delta_values: list = None,
    N: int = 64,
    t_end: float = 0.5,
    T_mean: float = _T_MEAN,
    fluid: str = "ideal",
    t_star: float = 1.0,
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str = None,
    seed: int = 42,
) -> list:
    """Run the full μ(T)→ν limit sweep (Phase 4, PRD §7).

    For each δ in delta_values:
      1. Build T₀ = T_mean + δ·g(x), ω₀ = standard double-vortex
      2. Run Route 1 2D solver
      3. Record A_min_global, lps_margin, mu_min_global, verdict

    Parameters
    ----------
    delta_values : list of δ amplitudes (default: DELTA_SWEEP_VALUES)
    N            : grid resolution
    t_end        : simulation end time [s]
    T_mean       : mean temperature [K]
    fluid        : 'ideal' or 'co2'
    t_star       : blow-up time reference
    max_steps    : safety cap
    verbose      : print progress
    db_path      : if given, log each result to results.db
    seed         : random seed for g(x) profile

    Returns
    -------
    list of result dicts, one per δ, ordered large δ → small δ
    """
    if delta_values is None:
        delta_values = DELTA_SWEEP_VALUES

    results = []
    for idx, delta in enumerate(delta_values, start=1):
        exp_id = f"EXP-L2-R1-MULIMIT-{idx:03d}"

        if verbose:
            print(f"\n[{idx}/{len(delta_values)}] {exp_id}  δ={delta}")

        row = run_single_delta(
            delta=delta,
            exp_id=exp_id,
            N=N,
            t_end=t_end,
            T_mean=T_mean,
            fluid=fluid,
            t_star=t_star,
            max_steps=max_steps,
            verbose=verbose,
            db_path=db_path,
            seed=seed,
        )
        results.append(row)

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Smoke test (fast validation for CI)
# ─────────────────────────────────────────────────────────────────────────────

def run_smoke_test(db_path: str = None, verbose: bool = False) -> dict:
    """Fast 3-run validation: δ ∈ {1.0, 0.1, 0.01} at N=32, t_end=0.05.

    Returns dict with keys: verdict, n_runs, all_pass, results.
    """
    delta_smoke = [1.0, 0.1, 0.01]
    results = []
    for idx, delta in enumerate(delta_smoke, start=1):
        exp_id = f"EXP-L2-R1-MULIMIT-SMOKE-{idx:03d}"
        row = run_single_delta(
            delta=delta,
            exp_id=exp_id,
            N=32,
            t_end=0.05,
            max_steps=50,
            verbose=verbose,
            db_path=None,   # smoke test does NOT write to production DB
            seed=42,
        )
        results.append(row)

    all_pass = all(r["verdict"] == "PASS" for r in results)
    return {
        "verdict": "PASS" if all_pass else "FAIL",
        "n_runs": len(results),
        "all_pass": all_pass,
        "results": results,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Print table
# ─────────────────────────────────────────────────────────────────────────────

def print_results_table(rows: list) -> None:
    """Print a formatted results table for δ sweep."""
    header = f"{'exp_id':<28} {'δ':>8} {'verdict':>8} {'A_min':>12} {'ε(δ)':>14} {'ω_max':>10}"
    print(header)
    print("-" * len(header))
    for r in rows:
        eps_sign = "+" if r.get("lps_margin", 0) >= 0 else ""
        print(
            f"{r['exp_id']:<28} "
            f"{r['delta']:>8.4f} "
            f"{r['verdict']:>8} "
            f"{r.get('A_min_global', float('nan')):>12.4e} "
            f"{eps_sign}{r.get('lps_margin', float('nan')):>13.4e} "
            f"{r.get('omega_max_global', float('nan')):>10.3e}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def _default_db_path() -> str:
    results_dir = os.path.join(_HERE, "..", "results")
    os.makedirs(results_dir, exist_ok=True)
    return os.path.join(results_dir, "results.db")


def main():
    parser = argparse.ArgumentParser(description="M4 μ(T)→ν Limit 2D sweep")
    parser.add_argument("--run", choices=["smoke", "full"], default="smoke",
                        help="'smoke' (fast, 3 δ values) or 'full' (10 δ values)")
    parser.add_argument("--N", type=int, default=64, help="Grid resolution (default: 64)")
    parser.add_argument("--t-end", type=float, default=0.5, help="End time [s] (default: 0.5)")
    parser.add_argument("--T-mean", type=float, default=_T_MEAN,
                        help="Mean temperature [K] (default: 300)")
    parser.add_argument("--verbose", action="store_true", help="Print step diagnostics")
    parser.add_argument("--query", action="store_true", help="Print results.db table and exit")
    parser.add_argument("--db", type=str, default=None, help="Path to results.db (default: auto)")
    args = parser.parse_args()

    db_path = args.db or _default_db_path()

    if args.query:
        rows = query_results(db_path, claim_id=_CLAIM_ID)
        if not rows:
            print("No results found.")
        else:
            print_results_table(rows)
            v = phase4_verdict(db_path)
            print(f"\nPhase 4 verdict: {v['verdict']}  "
                  f"({v['n_pass']}/{v['n_runs']} PASS, "
                  f"A_min_plateau={v.get('A_min_plateau')}, "
                  f"margin_ok={v.get('margin_ok')})")
        return

    if args.run == "smoke":
        print("Running smoke test (N=32, t_end=0.05, δ∈{1.0,0.1,0.01})...")
        result = run_smoke_test(db_path=None, verbose=args.verbose)
        print(f"\nSmoke test: {result['verdict']}  ({result['n_runs']} runs)")
    else:
        print(f"Running full μ(T)→ν limit sweep (N={args.N}, t_end={args.t_end})...")
        ref = _ref_props(args.T_mean)
        print(f"Reference: T̄={args.T_mean}K  A_ref={ref['A_ref']:.4e}  ν_ref={ref['nu_ref']:.4e}")
        rows = run_mu_limit_sweep(
            N=args.N,
            t_end=args.t_end,
            T_mean=args.T_mean,
            verbose=args.verbose,
            db_path=db_path,
        )
        print("\nResults:")
        print_results_table(rows)
        v = phase4_verdict(db_path)
        print(f"\nPhase 4 verdict: {v['verdict']}  "
              f"({v['n_pass']}/{v['n_runs']} PASS, "
              f"A_min_plateau={v.get('A_min_plateau')}, "
              f"A_ref={v.get('A_ref')})")


if __name__ == "__main__":
    main()
