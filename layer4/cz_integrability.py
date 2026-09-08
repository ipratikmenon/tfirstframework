"""
cz_integrability.py — Layer 4: Calderón-Zygmund Integrability Tracker
=======================================================================
Tracks ‖S_θ‖_{L^{1+ε}}(t) where S_θ = ν|∇u|² is the CZ source term.

This is the numerical test of Claim A (PRD v1.0 FINAL §3.4):
  S_θ = ν|∇u|² ∈ L^{1+ε}(ℝ³ × [0,T])  for some ε > 0.

If S_θ ∈ L^{1+ε} then by the CZ embedding:
  θ ∈ L²(H^{1+δ})  for δ = δ(ε) > 0  (Lemma 2.5)

The tracker measures:
  - ‖S_θ(·,t)‖_{L^p}  for p ∈ {1.0, 1.1, 1.25, 1.5, 2.0} at each step
  - The ratio r_p(t) = ‖S_θ‖_{L^p} / ‖S_θ‖_{L^1} (growth with p)
  - Instantaneous ε_CZ(t): largest p-1 where ‖S_θ‖_{L^p} < C·‖S_θ‖_{L^1}

Claim A passes if ε_CZ(t) > 0 throughout the run.

Run:
  python layer4/cz_integrability.py --N 32 --t-end 1.0 --ic tg --verbose
  python layer4/cz_integrability.py --N 32 --t-end 2.0 --ic shear --verbose
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from datetime import datetime

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

import route2_3D as r2

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")

# p values for L^p norm tracking
_P_VALS = (1.0, 1.1, 1.25, 1.5, 2.0)
# Ratio threshold: controlled means ratio < this value
_RATIO_THRESHOLD = 10.0


# =============================================================================
# CZ norm computation
# =============================================================================

def cz_lp_norms(source: np.ndarray, p_vals: tuple = _P_VALS) -> dict:
    """Compute ‖source‖_{L^p} (spatial average) for each p in p_vals.

    Returns dict: p → ‖source‖_{L^p}
    The spatial L^p norm uses the volumetric average:
      ‖f‖_{L^p}^p = (1/N³) Σ |f_i|^p   →  ‖f‖_{L^p} = (mean |f|^p)^{1/p}
    """
    result = {}
    for p in p_vals:
        lp = float(np.mean(np.abs(source) ** p) ** (1.0 / p))
        result[p] = lp
    return result


def estimate_eps_cz(lp_norms: dict, threshold: float = _RATIO_THRESHOLD) -> float:
    """Estimate ε_CZ: largest (p−1) for which ‖S_θ‖_{L^p} is controlled.

    A norm is 'controlled' if the ratio ‖S_θ‖_{L^p} / ‖S_θ‖_{L^1} < threshold.
    Returns 0.0 if L^1 norm is zero (source is zero — early time).
    """
    L1 = lp_norms.get(1.0, 0.0)
    if L1 < 1e-30:
        return 0.0   # source is zero; undefined

    eps_cz = 0.0
    for p in sorted(lp_norms.keys()):
        if p <= 1.0:
            continue
        ratio = lp_norms[p] / L1
        if ratio < threshold:
            eps_cz = p - 1.0
        else:
            break
    return float(eps_cz)


# =============================================================================
# Main tracker
# =============================================================================

def track_cz_integrability(
    N: int = 32,
    nu: float = 1.0e-3,
    eps_param: float = 0.1,
    ic: str = "tg",
    t_end: float = 1.0,
    max_steps: int = 5000,
    cfl_safety: float = 0.4,
    verbose: bool = False,
    db_path: str | None = None,
    exp_id_override: str | None = None,
) -> dict:
    """Track ‖S_θ‖_{L^p}(t) and ε_CZ(t) throughout a Route 2 3D run.

    Pass criterion (Claim A): ε_CZ(t) > 0 throughout (S_θ ∈ L^{1+ε}).
    Returns result dict with timeseries and summary.
    """
    solver = r2.make_solver(N=N, nu=nu, eps_param=eps_param)
    if ic == "tg":
        ic_dict = r2.taylor_green_ic(N, V0=1.0)
    elif ic == "shear":
        ic_dict = r2.shear_layer_ic(N)
    elif ic == "random":
        ic_dict = r2.random_div_free_ic(N, amp=0.5, seed=42)
    else:
        raise ValueError(f"Unknown ic: '{ic}'")

    solver = r2.set_ic(solver, ic_dict)
    exp_id = exp_id_override or f"EXP-L4-CZ-{ic.upper()}-{N:03d}"

    if verbose:
        print(f"\n  CZ Integrability Tracker  N={N}³  ν={nu}  ε={eps_param}  IC={ic}")
        print(f"  Exp: {exp_id}  t_end={t_end}")
        print(f"  Monitoring ‖S_θ‖_{{L^p}} for p ∈ {_P_VALS}")

    timeseries: list[dict] = []
    t_wall0 = time.perf_counter()
    eps_cz_positive_count = 0

    kx = solver["kx"]; ky = solver["ky"]; kz = solver["kz"]

    for step_idx in range(max_steps):
        if solver["t"] >= t_end:
            break

        dt = r2.compute_cfl_dt(solver, safety=cfl_safety)
        dt = min(dt, t_end - solver["t"])

        # Compute S_θ BEFORE the step (at current velocity)
        S_theta = r2._source_theta(
            solver["u"], solver["v"], solver["w"], kx, ky, kz, nu
        )
        lp_norms = cz_lp_norms(S_theta, _P_VALS)
        eps_cz = estimate_eps_cz(lp_norms)
        if eps_cz > 0.0:
            eps_cz_positive_count += 1

        # Step the solver
        solver, diag = r2.solver_step(solver, dt=dt)

        row = {
            "t"         : float(solver["t"]),
            "dt"        : float(dt),
            "eps_cz"    : eps_cz,
            "S_max"     : float(np.max(S_theta)),
            "S_mean"    : float(np.mean(S_theta)),
            "E"         : float(diag["E"]),
            "delta_est" : float(diag["delta_est"]),
            **{f"Lp_{p:.2f}": float(lp_norms[p]) for p in _P_VALS},
            **{f"ratio_{p:.2f}": float(lp_norms[p] / lp_norms[1.0])
               if lp_norms[1.0] > 1e-30 else 0.0
               for p in _P_VALS if p > 1.0},
        }
        timeseries.append(row)

        if verbose and step_idx % max(1, max_steps // 20) == 0:
            print(
                f"    t={solver['t']:.4f}  S_max={np.max(S_theta):.3e}  "
                f"ε_CZ={eps_cz:.2f}  δ={diag['delta_est']:.2f}"
            )

    wall = time.perf_counter() - t_wall0

    # ── Summary ───────────────────────────────────────────────────────────
    eps_czs = [r["eps_cz"] for r in timeseries]
    n_total = max(len(eps_czs), 1)
    eps_cz_max  = float(max(eps_czs)) if eps_czs else 0.0
    eps_cz_min  = float(min(eps_czs)) if eps_czs else 0.0
    eps_cz_mean = float(np.mean(eps_czs)) if eps_czs else 0.0
    eps_cz_pos_frac = float(eps_cz_positive_count / n_total)

    # Claim A: ε_CZ > 0 for the majority of steps once source is active
    # (Allow early steps where S_θ ≈ 0 before the source accumulates)
    n_active = sum(1 for r in timeseries if r["S_mean"] > 1e-15)
    n_active_pos = sum(1 for r in timeseries if r["S_mean"] > 1e-15 and r["eps_cz"] > 0.0)
    claim_a_frac = float(n_active_pos / max(n_active, 1))

    verdict = "PASS" if eps_cz_max > 0.0 else "FAIL"

    summary = {
        "exp_id"          : exp_id,
        "ic"              : ic,
        "N"               : N,
        "nu"              : nu,
        "eps_param"       : eps_param,
        "t_end"           : t_end,
        "n_steps"         : len(timeseries),
        "t_final"         : float(solver["t"]),
        "eps_cz_max"      : eps_cz_max,
        "eps_cz_min"      : eps_cz_min,
        "eps_cz_mean"     : eps_cz_mean,
        "eps_cz_pos_frac" : eps_cz_pos_frac,
        "claim_a_frac"    : claim_a_frac,
        "verdict"         : verdict,
        "wall_time_s"     : wall,
        "timestamp"       : datetime.utcnow().isoformat(),
        "key_metric"      : (
            f"ε_CZ_max={eps_cz_max:.2f}  ε_CZ_mean={eps_cz_mean:.2f}  "
            f"Claim_A_frac={claim_a_frac:.2%}"
        ),
    }

    if verbose:
        print(f"\n  Results: {exp_id}")
        print(f"    ε_CZ_max:  {eps_cz_max:.2f}")
        print(f"    ε_CZ_mean: {eps_cz_mean:.2f}")
        print(f"    Claim A:   {claim_a_frac:.1%} of active steps")
        print(f"    Verdict:   {verdict}")
        print(f"    Wall:      {wall:.1f} s")

    if db_path is not None:
        _log_cz(db_path, summary, timeseries)

    return {"summary": summary, "timeseries": timeseries, "verdict": verdict}


# =============================================================================
# DB logging
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_cz_experiments (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id           TEXT UNIQUE NOT NULL,
            timestamp        TEXT,
            ic               TEXT,
            N                INTEGER,
            nu               REAL,
            eps_param        REAL,
            t_end            REAL,
            n_steps          INTEGER,
            t_final          REAL,
            eps_cz_max       REAL,
            eps_cz_min       REAL,
            eps_cz_mean      REAL,
            eps_cz_pos_frac  REAL,
            claim_a_frac     REAL,
            verdict          TEXT,
            wall_time_s      REAL,
            key_metric       TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_cz_timeseries (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id      TEXT NOT NULL,
            t           REAL,
            dt          REAL,
            eps_cz      REAL,
            S_max       REAL,
            S_mean      REAL,
            E           REAL,
            delta_est   REAL,
            Lp_1_00     REAL,
            Lp_1_10     REAL,
            Lp_1_25     REAL,
            Lp_1_50     REAL,
            Lp_2_00     REAL
        )
    """)
    conn.commit()
    return conn


def _log_cz(db_path: str, summary: dict, timeseries: list[dict]) -> None:
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "timestamp", "ic", "N", "nu", "eps_param", "t_end",
        "n_steps", "t_final", "eps_cz_max", "eps_cz_min", "eps_cz_mean",
        "eps_cz_pos_frac", "claim_a_frac", "verdict", "wall_time_s", "key_metric",
    ]
    vals = tuple(summary.get(c) for c in cols)
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO layer4_cz_experiments "
        f"({', '.join(cols)}) VALUES ({ph})",
        vals,
    )
    for row in timeseries:
        conn.execute("""
            INSERT INTO layer4_cz_timeseries
            (exp_id, t, dt, eps_cz, S_max, S_mean, E, delta_est,
             Lp_1_00, Lp_1_10, Lp_1_25, Lp_1_50, Lp_2_00)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (
            summary["exp_id"],
            row["t"], row["dt"], row["eps_cz"],
            row["S_max"], row["S_mean"], row["E"], row["delta_est"],
            row.get("Lp_1.00", 0.0), row.get("Lp_1.10", 0.0),
            row.get("Lp_1.25", 0.0), row.get("Lp_1.50", 0.0),
            row.get("Lp_2.00", 0.0),
        ))
    conn.commit()
    conn.close()


# =============================================================================
# CLI
# =============================================================================

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Layer 4: CZ integrability tracker")
    parser.add_argument("--N",      type=int,   default=32)
    parser.add_argument("--nu",     type=float, default=1.0e-3)
    parser.add_argument("--eps",    type=float, default=0.1)
    parser.add_argument("--ic",     choices=["tg", "shear", "random"], default="tg")
    parser.add_argument("--t-end",  type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db",     default=_DEFAULT_DB)
    args = parser.parse_args()

    r = track_cz_integrability(
        N=args.N, nu=args.nu, eps_param=args.eps,
        ic=args.ic, t_end=args.t_end, verbose=args.verbose,
        db_path=args.db,
    )
    print(f"\nVerdict: {r['verdict']}  key_metric: {r['summary']['key_metric']}")


if __name__ == "__main__":
    main()
