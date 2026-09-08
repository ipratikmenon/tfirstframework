"""
delta_extractor.py — Layer 4: δ(t) Timeseries Extractor
=========================================================
Runs a Route 2 3D solver and extracts the full δ(t) timeseries — the
Sobolev regularity gain of θ above H¹ — along with H^s norms at each step.

δ > 0 is the central quantity in Lemma 2.5 (PRD v1.0 FINAL):
  θ ∈ L²(H^{1+δ}) for any δ > 0  ⟹  u satisfies Prodi-Serrin  ⟹  global smoothness.

Output:
  - JSON timeseries: t, δ(t), H^s norms for s ∈ {0, 0.5, 1.0, 1.5, 2.0, 2.5}
  - Summary: δ_min, δ_max, δ_mean, first time δ > 0
  - SQLite row in layer4_delta_experiments table

Run:
  python layer4/delta_extractor.py --N 32 --t-end 1.0 --ic tg --verbose
  python layer4/delta_extractor.py --N 32 --t-end 2.0 --ic shear --verbose
  python layer4/delta_extractor.py --N 32 --t-end 1.0 --ic random --eps 0.0
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

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

import route2_3D as r2

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_S_LIST = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5)


# =============================================================================
# Core extraction
# =============================================================================

def extract_delta_timeseries(
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
    """Run Route 2 3D solver and extract δ(t) timeseries.

    Parameters
    ----------
    N          : Grid size (N³)
    nu         : Kinematic viscosity
    eps_param  : ε in μ_eff = ν + ε·f(θ)
    ic         : 'tg' | 'shear' | 'random'
    t_end      : Integration end time
    max_steps  : Maximum timesteps
    cfl_safety : CFL safety factor
    verbose    : Print progress
    db_path    : SQLite DB path (None = don't log)

    Returns
    -------
    dict with keys: exp_id, timeseries, summary, verdict
    """
    # ── Build solver and IC ──────────────────────────────────────────────────
    solver = r2.make_solver(N=N, nu=nu, eps_param=eps_param)

    if ic == "tg":
        ic_dict = r2.taylor_green_ic(N, V0=1.0)
    elif ic == "shear":
        ic_dict = r2.shear_layer_ic(N)
    elif ic == "random":
        ic_dict = r2.random_div_free_ic(N, amp=0.5, seed=42)
    else:
        raise ValueError(f"Unknown ic: '{ic}'. Use 'tg', 'shear', 'random'.")

    solver = r2.set_ic(solver, ic_dict)
    exp_id = exp_id_override or f"EXP-L4-DELTA-{ic.upper()}-{N:03d}"

    if verbose:
        print(f"\n  δ(t) Extractor  N={N}³  ν={nu}  ε={eps_param}  IC={ic}")
        print(f"  Exp: {exp_id}  t_end={t_end}")

    # ── Integration loop with full H^s recording ──────────────────────────
    timeseries: list[dict] = []
    t_wall0 = time.perf_counter()
    delta_first_positive: float | None = None

    for step_idx in range(max_steps):
        if solver["t"] >= t_end:
            break

        dt = r2.compute_cfl_dt(solver, safety=cfl_safety)
        dt = min(dt, t_end - solver["t"])
        solver, diag = r2.solver_step(solver, dt=dt)

        # H^s norms for full tracking
        kx, ky, kz = solver["kx"], solver["ky"], solver["kz"]
        Hs = r2.sobolev_norms(solver["theta"], kx, ky, kz, s_list=_S_LIST)

        row = {
            "t"          : float(solver["t"]),
            "dt"         : float(dt),
            "delta_est"  : float(diag["delta_est"]),
            "E"          : float(diag["E"]),
            "theta_max"  : float(diag["theta_max"]),
            "S_theta_max": float(diag["S_theta_max"]),
            "eps_lps"    : float(diag["eps_lps"]),
            **{k: float(v) for k, v in Hs.items()},
        }
        timeseries.append(row)

        if delta_first_positive is None and diag["delta_est"] > 0.0:
            delta_first_positive = float(solver["t"])

        if verbose and step_idx % max(1, max_steps // 20) == 0:
            print(
                f"    t={solver['t']:.4f}  δ={diag['delta_est']:.2f}  "
                f"E={diag['E']:.4e}  θ_max={diag['theta_max']:.3e}"
            )

    wall = time.perf_counter() - t_wall0

    # ── Summary ────────────────────────────────────────────────────────────
    deltas = [r["delta_est"] for r in timeseries]
    delta_min  = float(min(deltas)) if deltas else 0.0
    delta_max  = float(max(deltas)) if deltas else 0.0
    delta_mean = float(np.mean(deltas)) if deltas else 0.0

    # δ sustained: fraction of steps where δ > 0
    delta_positive_frac = float(sum(1 for d in deltas if d > 0.0) / max(len(deltas), 1))

    verdict = "PASS" if delta_max > 0.0 else "FAIL"

    summary = {
        "exp_id"              : exp_id,
        "ic"                  : ic,
        "N"                   : N,
        "nu"                  : nu,
        "eps_param"           : eps_param,
        "t_end"               : t_end,
        "n_steps"             : len(timeseries),
        "t_final"             : float(solver["t"]),
        "delta_min"           : delta_min,
        "delta_max"           : delta_max,
        "delta_mean"          : delta_mean,
        "delta_positive_frac" : delta_positive_frac,
        "delta_first_positive": delta_first_positive,
        "verdict"             : verdict,
        "wall_time_s"         : wall,
        "timestamp"           : datetime.utcnow().isoformat(),
        "key_metric"          : (
            f"δ_max={delta_max:.2f}  δ_mean={delta_mean:.2f}  "
            f"δ_frac={delta_positive_frac:.2%}  "
            f"t_first_δ>0={delta_first_positive}"
        ),
    }

    if verbose:
        print(f"\n  Results: {exp_id}")
        print(f"    δ_max:   {delta_max:.2f}")
        print(f"    δ_mean:  {delta_mean:.2f}")
        print(f"    δ>0:     {delta_positive_frac:.1%} of steps")
        print(f"    t(δ>0):  {delta_first_positive}")
        print(f"    Verdict: {verdict}")
        print(f"    Wall:    {wall:.1f} s")

    if db_path is not None:
        _log_delta(db_path, summary, timeseries)

    return {"summary": summary, "timeseries": timeseries, "verdict": verdict}


# =============================================================================
# ε sweep: δ(ε) at fixed IC and t_end
# =============================================================================

def delta_eps_sweep(
    N: int = 32,
    nu: float = 1.0e-3,
    ic: str = "tg",
    t_end: float = 1.0,
    eps_vals: tuple = (1.0, 0.1, 0.01, 0.001, 0.0),
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str | None = None,
) -> list[dict]:
    """Run δ extractor for each ε in eps_vals.

    The key test: does δ_max > 0 hold at ε=0 (exact Prize equations)?
    If yes, the CZ source structure alone (S_θ = ν|∇u|²) delivers Lemma 2.5.
    """
    results = []
    for eps in eps_vals:
        if verbose:
            print(f"\n  ε = {eps:.3f}  {'(EXACT PRIZE)' if eps == 0.0 else ''}")
        r = extract_delta_timeseries(
            N=N, nu=nu, eps_param=eps, ic=ic, t_end=t_end,
            max_steps=max_steps, verbose=False,
            db_path=db_path,
            exp_id_override=f"EXP-L4-DELTA-EPS-{ic.upper()}-{N:03d}-eps{eps}",
        )
        s = r["summary"]
        if verbose:
            print(
                f"    δ_max={s['delta_max']:.2f}  δ_mean={s['delta_mean']:.2f}  "
                f"{s['verdict']}"
            )
        results.append(s)

    if verbose:
        print("\n  ε sweep summary:")
        print(f"  {'ε':>8}  {'δ_max':>7}  {'δ_mean':>7}  Verdict")
        for s in results:
            print(
                f"  {s['eps_param']:8.3f}  {s['delta_max']:7.2f}  "
                f"{s['delta_mean']:7.2f}  {s['verdict']}"
            )
        eps0 = [s for s in results if s["eps_param"] == 0.0]
        if eps0 and eps0[0]["delta_max"] > 0.0:
            print("\n  ✓ δ > 0 at ε=0: EXACT PRIZE EQUATIONS confirm Lemma 2.5")

    return results


# =============================================================================
# DB logging
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_delta_experiments (
            id                   INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id               TEXT UNIQUE NOT NULL,
            timestamp            TEXT NOT NULL,
            ic                   TEXT,
            N                    INTEGER,
            nu                   REAL,
            eps_param            REAL,
            t_end                REAL,
            n_steps              INTEGER,
            t_final              REAL,
            delta_min            REAL,
            delta_max            REAL,
            delta_mean           REAL,
            delta_positive_frac  REAL,
            delta_first_positive REAL,
            verdict              TEXT,
            wall_time_s          REAL,
            key_metric           TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_delta_timeseries (
            id         INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id     TEXT NOT NULL,
            t          REAL,
            dt         REAL,
            delta_est  REAL,
            E          REAL,
            theta_max  REAL,
            S_theta_max REAL,
            eps_lps    REAL,
            H0_0       REAL,
            H0_5       REAL,
            H1_0       REAL,
            H1_5       REAL,
            H2_0       REAL,
            H2_5       REAL
        )
    """)
    conn.commit()
    return conn


def _log_delta(db_path: str, summary: dict, timeseries: list[dict]) -> None:
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "timestamp", "ic", "N", "nu", "eps_param",
        "t_end", "n_steps", "t_final",
        "delta_min", "delta_max", "delta_mean",
        "delta_positive_frac", "delta_first_positive",
        "verdict", "wall_time_s", "key_metric",
    ]
    vals = tuple(summary.get(c) for c in cols)
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO layer4_delta_experiments "
        f"({', '.join(cols)}) VALUES ({ph})",
        vals,
    )
    # Timeseries rows
    s_keys = [f"H{s:.1f}" for s in _S_LIST]
    ts_cols = [f"H{s:.1f}".replace(".", "_") for s in _S_LIST]
    for row in timeseries:
        conn.execute("""
            INSERT INTO layer4_delta_timeseries
            (exp_id, t, dt, delta_est, E, theta_max, S_theta_max, eps_lps,
             H0_0, H0_5, H1_0, H1_5, H2_0, H2_5)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (
            summary["exp_id"],
            row["t"], row["dt"], row["delta_est"],
            row["E"], row["theta_max"], row["S_theta_max"], row["eps_lps"],
            row.get("H0.0", 0.0), row.get("H0.5", 0.0), row.get("H1.0", 0.0),
            row.get("H1.5", 0.0), row.get("H2.0", 0.0), row.get("H2.5", 0.0),
        ))
    conn.commit()
    conn.close()


# =============================================================================
# Query utilities
# =============================================================================

def query_delta_results(db_path: str, ic: str | None = None) -> list[dict]:
    """Return all δ experiment rows, optionally filtered by IC type."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    q = "SELECT * FROM layer4_delta_experiments"
    args: tuple = ()
    if ic:
        q += " WHERE ic = ?"
        args = (ic,)
    q += " ORDER BY timestamp DESC"
    rows = [dict(r) for r in conn.execute(q, args)]
    conn.close()
    return rows


def print_delta_table(rows: list[dict]) -> None:
    """Pretty-print a δ experiment summary table."""
    if not rows:
        print("  (no rows)")
        return
    print(f"\n  {'exp_id':<45} {'ε':>7}  {'δ_max':>6}  {'δ_mean':>6}  {'Verdict'}")
    print("  " + "-" * 80)
    for r in rows:
        print(
            f"  {r['exp_id']:<45} {r['eps_param']:7.3f}  "
            f"{r['delta_max']:6.2f}  {r['delta_mean']:6.2f}  {r['verdict']}"
        )


# =============================================================================
# CLI
# =============================================================================

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Layer 4: δ(t) timeseries extractor")
    parser.add_argument("--N",      type=int,   default=32)
    parser.add_argument("--nu",     type=float, default=1.0e-3)
    parser.add_argument("--eps",    type=float, default=0.1)
    parser.add_argument("--ic",     choices=["tg", "shear", "random"], default="tg")
    parser.add_argument("--t-end",  type=float, default=1.0)
    parser.add_argument("--mode",   choices=["single", "eps-sweep", "query"], default="single")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db",     default=_DEFAULT_DB)
    parser.add_argument("--save-ts", action="store_true",
                        help="Save timeseries JSON to results/ directory")
    args = parser.parse_args()

    if args.mode == "single":
        r = extract_delta_timeseries(
            N=args.N, nu=args.nu, eps_param=args.eps,
            ic=args.ic, t_end=args.t_end, verbose=args.verbose,
            db_path=args.db,
        )
        print(f"\nVerdict: {r['verdict']}  key_metric: {r['summary']['key_metric']}")
        if args.save_ts:
            ts_path = os.path.join(_ROOT, "results",
                                   f"{r['summary']['exp_id']}_timeseries.json")
            with open(ts_path, "w") as f:
                json.dump(r["timeseries"], f, indent=2)
            print(f"Timeseries saved → {ts_path}")

    elif args.mode == "eps-sweep":
        results = delta_eps_sweep(
            N=args.N, nu=args.nu, ic=args.ic, t_end=args.t_end,
            verbose=args.verbose, db_path=args.db,
        )
        print_delta_table(results)

    elif args.mode == "query":
        rows = query_delta_results(args.db)
        print_delta_table(rows)


if __name__ == "__main__":
    main()
