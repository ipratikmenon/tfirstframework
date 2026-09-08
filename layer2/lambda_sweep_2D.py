"""
lambda_sweep_2D.py — Lambda Sweep Driver with results.db Logging (PRD §7.3)
==============================================================================
Runs the full 2D lambda sweep experiment from PRD §7.3 and logs every result
to `results/results.db` (SQLite).

Experiment IDs:
  EXP-L2-R1-CCF-001 … EXP-L2-R1-CCF-006   (CCF profiles, 6 λ values)
  EXP-L2-R1-BOUS-001 … EXP-L2-R1-BOUS-003 (Boussinesq profiles)
  EXP-L2-R1-ADV-001                         (adversarial min, λ→0)

Claim tested: Conjecture 3.5 — μ(T) scaling defeats Wang et al. self-similar
blow-up for ALL λ > 0. PASS = all λ suppressed, no λ_c threshold found.

Pass criterion (PRD §7.3):
  For every run:
    A_min_global > 0          — thermal diffusivity stays positive (second law)
    suppression_ratio ≥ 1     — viscous heating rate ≥ blow-up rate
    omega_max < blow-up bound — no vorticity blowup detected

Run:
  python layer2/lambda_sweep_2D.py                  # full sweep, 64²
  python layer2/lambda_sweep_2D.py --smoke-test     # fast check, 32²
  python layer2/lambda_sweep_2D.py --profile ccf    # CCF only
  python layer2/lambda_sweep_2D.py --profile bous   # Boussinesq only
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

from T_solver_2D import (
    make_solver, set_ic, run, run_lambda_sweep,
    LAMBDA_SWEEP_VALUES,
)
from self_similar_IC import (
    generate_CCF_profile,
    generate_Boussinesq_profile,
    generate_adversarial_min,
)

# ─────────────────────────────────────────────────────────────────────────────
# Results database
# ─────────────────────────────────────────────────────────────────────────────

_DB_PATH = os.path.join(_ROOT, "results", "results.db")


def _ensure_db(db_path: str = _DB_PATH) -> sqlite3.Connection:
    """Open (or create) the results SQLite database and ensure the table exists."""
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS experiments (
            id            INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id        TEXT    NOT NULL UNIQUE,
            claim_id      TEXT    NOT NULL,
            timestamp     TEXT    NOT NULL,
            layer         INTEGER NOT NULL,
            route         INTEGER NOT NULL,
            profile_type  TEXT    NOT NULL,
            lambda_val    REAL,
            grid_N        INTEGER NOT NULL,
            t_end         REAL    NOT NULL,
            fluid         TEXT    NOT NULL,
            n_steps       INTEGER,
            t_final       REAL,
            verdict       TEXT    NOT NULL,
            A_min_global  REAL,
            omega_max_global REAL,
            suppression_ratio_final REAL,
            D_S_ratio_min REAL,
            T_min_final   REAL,
            T_max_final   REAL,
            key_metric    TEXT,
            notes         TEXT,
            params_json   TEXT
        )
    """)
    conn.commit()
    return conn


def log_result(conn: sqlite3.Connection, row: dict) -> None:
    """Insert one experiment result. If exp_id already exists, update it."""
    cols = [
        "exp_id", "claim_id", "timestamp", "layer", "route",
        "profile_type", "lambda_val", "grid_N", "t_end", "fluid",
        "n_steps", "t_final", "verdict", "A_min_global", "omega_max_global",
        "suppression_ratio_final", "D_S_ratio_min", "T_min_final", "T_max_final",
        "key_metric", "notes", "params_json",
    ]
    placeholders = ", ".join("?" for _ in cols)
    col_names = ", ".join(cols)
    values = tuple(row.get(c) for c in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO experiments ({col_names}) VALUES ({placeholders})",
        values,
    )
    conn.commit()


# ─────────────────────────────────────────────────────────────────────────────
# Sweep helpers
# ─────────────────────────────────────────────────────────────────────────────

def _extract_suppression_final(history: list) -> float:
    """Get the final finite suppression ratio from step history."""
    for d in reversed(history):
        sr = d.get("suppression_ratio", float("inf"))
        if not (sr == float("inf") or sr != sr):  # not inf and not NaN
            return float(sr)
    return float("inf")


def _extract_ds_ratio_min(history: list) -> float | None:
    """Get minimum D/S ratio if recorded in history."""
    ratios = [d["D_S_ratio"] for d in history if "D_S_ratio" in d]
    return float(min(ratios)) if ratios else None


def _extract_T_final(result: dict) -> tuple[float, float]:
    """Get final T_min, T_max from last history entry."""
    hist = result.get("history", [])
    if not hist:
        return float("nan"), float("nan")
    last = hist[-1]
    return last.get("T_min", float("nan")), last.get("T_max", float("nan"))


# ─────────────────────────────────────────────────────────────────────────────
# CCF sweep
# ─────────────────────────────────────────────────────────────────────────────

def run_ccf_sweep(
    lambda_values: list | None = None,
    N: int = 64,
    t_end: float = 0.5,
    fluid: str = "ideal",
    P: float | None = None,
    T_base: float = 300.0,
    t_star: float = 1.0,
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str = _DB_PATH,
) -> list[dict]:
    """
    Run CCF lambda sweep and log each result to results.db.

    Experiment IDs: EXP-L2-R1-CCF-001 … EXP-L2-R1-CCF-006
    Claim: Conjecture 3.5 — A_min_global > 0 and suppression_ratio ≥ 1 for all λ.
    """
    if lambda_values is None:
        lambda_values = LAMBDA_SWEEP_VALUES

    conn = _ensure_db(db_path)
    timestamp = datetime.datetime.utcnow().isoformat()

    sweep_results = []

    for idx, lam in enumerate(lambda_values, start=1):
        exp_id = f"EXP-L2-R1-CCF-{idx:03d}"
        claim_id = "Conjecture_3.5"

        if verbose:
            print(f"\n[{exp_id}] λ={lam:.4f}  N={N}  t_end={t_end}")

        t0 = time.time()

        # Generate CCF IC
        ic = generate_CCF_profile(lam, N=N)
        omega_ic_max = float(np.max(np.abs(ic["omega"]))) or 1.0

        # Build and run solver
        solver = make_solver(
            N=N, fluid=fluid, P=P, T_base=T_base,
            omega_max=omega_ic_max, t_star=t_star, lambda_val=lam,
        )
        solver = set_ic(solver, ic)
        result = run(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)

        elapsed = time.time() - t0

        # Extract metrics
        supp_ratio_final = _extract_suppression_final(result["history"])
        ds_ratio_min     = _extract_ds_ratio_min(result["history"])
        T_min_f, T_max_f = _extract_T_final(result)

        verdict = result["verdict"]
        A_min   = result["A_min_global"]

        key_metric = (
            f"A_min={A_min:.4e}  "
            f"supp_ratio_final={supp_ratio_final:.4e}  "
            f"omega_max={result['omega_max_global']:.4e}  "
            f"steps={result['n_steps']}"
        )

        params = {
            "N": N, "t_end": t_end, "fluid": fluid, "T_base": T_base,
            "t_star": t_star, "lambda_val": lam, "max_steps": max_steps,
            "elapsed_s": round(elapsed, 2),
        }

        row = {
            "exp_id"                 : exp_id,
            "claim_id"               : claim_id,
            "timestamp"              : timestamp,
            "layer"                  : 2,
            "route"                  : 1,
            "profile_type"           : "CCF",
            "lambda_val"             : lam,
            "grid_N"                 : N,
            "t_end"                  : t_end,
            "fluid"                  : fluid,
            "n_steps"                : result["n_steps"],
            "t_final"                : result["t_final"],
            "verdict"                : verdict,
            "A_min_global"           : A_min,
            "omega_max_global"       : result["omega_max_global"],
            "suppression_ratio_final": supp_ratio_final if not np.isinf(supp_ratio_final) else None,
            "D_S_ratio_min"          : ds_ratio_min,
            "T_min_final"            : T_min_f,
            "T_max_final"            : T_max_f,
            "key_metric"             : key_metric,
            "notes"                  : f"PRD §7.3 CCF λ={lam:.4f}",
            "params_json"            : json.dumps(params),
        }

        log_result(conn, row)

        if verbose:
            status = "✅ PASS" if verdict == "PASS" else "❌ FAIL"
            print(f"  {status}  A_min={A_min:.4e}  ω_max={result['omega_max_global']:.4e}"
                  f"  steps={result['n_steps']}  [{elapsed:.1f}s]")

        result.update({"exp_id": exp_id, "lambda_val": lam})
        sweep_results.append(result)

    conn.close()
    return sweep_results


# ─────────────────────────────────────────────────────────────────────────────
# Boussinesq sweep
# ─────────────────────────────────────────────────────────────────────────────

# Boussinesq λ values from PRD §4.3 (3 unstable profiles)
BOUS_LAMBDA_VALUES = [1.9206, 1.3991, 1.1843]


def run_boussinesq_sweep(
    lambda_values: list | None = None,
    N: int = 64,
    t_end: float = 0.5,
    fluid: str = "ideal",
    P: float | None = None,
    T_base: float = 300.0,
    t_star: float = 1.0,
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str = _DB_PATH,
) -> list[dict]:
    """
    Run Boussinesq profile sweep and log to results.db.

    Experiment IDs: EXP-L2-R1-BOUS-001 … EXP-L2-R1-BOUS-003
    Claim: Conjecture 3.5 holds for Boussinesq unstable profiles too.
    """
    if lambda_values is None:
        lambda_values = BOUS_LAMBDA_VALUES

    conn = _ensure_db(db_path)
    timestamp = datetime.datetime.utcnow().isoformat()
    sweep_results = []

    for idx, lam in enumerate(lambda_values, start=1):
        exp_id = f"EXP-L2-R1-BOUS-{idx:03d}"
        claim_id = "Conjecture_3.5_Boussinesq"

        if verbose:
            print(f"\n[{exp_id}] λ={lam:.4f}  N={N}  t_end={t_end}")

        t0 = time.time()

        ic = generate_Boussinesq_profile(lam, N=N)
        omega_ic_max = float(np.max(np.abs(ic["omega"]))) or 1.0

        solver = make_solver(
            N=N, fluid=fluid, P=P, T_base=T_base,
            omega_max=omega_ic_max, t_star=t_star, lambda_val=lam,
        )
        solver = set_ic(solver, ic)
        result = run(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)

        elapsed = time.time() - t0

        supp_ratio_final = _extract_suppression_final(result["history"])
        ds_ratio_min     = _extract_ds_ratio_min(result["history"])
        T_min_f, T_max_f = _extract_T_final(result)
        verdict = result["verdict"]
        A_min   = result["A_min_global"]

        key_metric = (
            f"A_min={A_min:.4e}  "
            f"supp_ratio_final={supp_ratio_final:.4e}  "
            f"omega_max={result['omega_max_global']:.4e}  "
            f"steps={result['n_steps']}"
        )

        params = {
            "N": N, "t_end": t_end, "fluid": fluid, "T_base": T_base,
            "t_star": t_star, "lambda_val": lam, "max_steps": max_steps,
            "elapsed_s": round(elapsed, 2),
        }

        row = {
            "exp_id"                 : exp_id,
            "claim_id"               : claim_id,
            "timestamp"              : timestamp,
            "layer"                  : 2,
            "route"                  : 1,
            "profile_type"           : "Boussinesq",
            "lambda_val"             : lam,
            "grid_N"                 : N,
            "t_end"                  : t_end,
            "fluid"                  : fluid,
            "n_steps"                : result["n_steps"],
            "t_final"                : result["t_final"],
            "verdict"                : verdict,
            "A_min_global"           : A_min,
            "omega_max_global"       : result["omega_max_global"],
            "suppression_ratio_final": supp_ratio_final if not np.isinf(supp_ratio_final) else None,
            "D_S_ratio_min"          : ds_ratio_min,
            "T_min_final"            : T_min_f,
            "T_max_final"            : T_max_f,
            "key_metric"             : key_metric,
            "notes"                  : f"PRD §7.3 Boussinesq λ={lam:.4f}",
            "params_json"            : json.dumps(params),
        }

        log_result(conn, row)

        if verbose:
            status = "✅ PASS" if verdict == "PASS" else "❌ FAIL"
            print(f"  {status}  A_min={A_min:.4e}  ω_max={result['omega_max_global']:.4e}"
                  f"  steps={result['n_steps']}  [{elapsed:.1f}s]")

        result.update({"exp_id": exp_id, "lambda_val": lam})
        sweep_results.append(result)

    conn.close()
    return sweep_results


# ─────────────────────────────────────────────────────────────────────────────
# Adversarial minimum
# ─────────────────────────────────────────────────────────────────────────────

def run_adversarial(
    lambda_val: float = 0.01,
    N: int = 64,
    t_end: float = 0.5,
    fluid: str = "ideal",
    P: float | None = None,
    T_base: float = 300.0,
    t_star: float = 1.0,
    max_steps: int = 2000,
    verbose: bool = False,
    db_path: str = _DB_PATH,
) -> dict:
    """
    Run the adversarial minimum IC (λ→0) and log to results.db.

    Experiment ID: EXP-L2-R1-ADV-001
    This is the hardest case for T-first — weakest thermal feedback.
    """
    exp_id   = "EXP-L2-R1-ADV-001"
    claim_id = "Conjecture_3.5_adversarial"

    if verbose:
        print(f"\n[{exp_id}] λ={lambda_val:.4f} (adversarial min)  N={N}")

    conn = _ensure_db(db_path)
    timestamp = datetime.datetime.utcnow().isoformat()
    t0 = time.time()

    ic = generate_adversarial_min(lambda_val, N=N)
    omega_ic_max = float(np.max(np.abs(ic["omega"]))) or 1.0

    solver = make_solver(
        N=N, fluid=fluid, P=P, T_base=T_base,
        omega_max=omega_ic_max, t_star=t_star, lambda_val=lambda_val,
    )
    solver = set_ic(solver, ic)
    result = run(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)

    elapsed = time.time() - t0

    supp_ratio_final = _extract_suppression_final(result["history"])
    ds_ratio_min     = _extract_ds_ratio_min(result["history"])
    T_min_f, T_max_f = _extract_T_final(result)
    verdict = result["verdict"]
    A_min   = result["A_min_global"]

    key_metric = (
        f"A_min={A_min:.4e}  "
        f"supp_ratio_final={supp_ratio_final:.4e}  "
        f"omega_max={result['omega_max_global']:.4e}  "
        f"steps={result['n_steps']}"
    )

    params = {
        "N": N, "t_end": t_end, "fluid": fluid, "T_base": T_base,
        "t_star": t_star, "lambda_val": lambda_val, "max_steps": max_steps,
        "elapsed_s": round(elapsed, 2),
    }

    row = {
        "exp_id"                 : exp_id,
        "claim_id"               : claim_id,
        "timestamp"              : timestamp,
        "layer"                  : 2,
        "route"                  : 1,
        "profile_type"           : "adversarial_min",
        "lambda_val"             : lambda_val,
        "grid_N"                 : N,
        "t_end"                  : t_end,
        "fluid"                  : fluid,
        "n_steps"                : result["n_steps"],
        "t_final"                : result["t_final"],
        "verdict"                : verdict,
        "A_min_global"           : A_min,
        "omega_max_global"       : result["omega_max_global"],
        "suppression_ratio_final": supp_ratio_final if not np.isinf(supp_ratio_final) else None,
        "D_S_ratio_min"          : ds_ratio_min,
        "T_min_final"            : T_min_f,
        "T_max_final"            : T_max_f,
        "key_metric"             : key_metric,
        "notes"                  : f"PRD §7.3 adversarial min λ→0 (λ={lambda_val})",
        "params_json"            : json.dumps(params),
    }

    log_result(conn, row)
    conn.close()

    if verbose:
        status = "✅ PASS" if verdict == "PASS" else "❌ FAIL"
        print(f"  {status}  A_min={A_min:.4e}  steps={result['n_steps']}  [{elapsed:.1f}s]")

    result.update({"exp_id": exp_id, "lambda_val": lambda_val})
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Results querying
# ─────────────────────────────────────────────────────────────────────────────

def query_results(
    db_path: str = _DB_PATH,
    claim_id: str | None = None,
    layer: int | None = None,
) -> list[dict]:
    """Load experiment rows from results.db. Optionally filter by claim_id or layer."""
    if not os.path.exists(db_path):
        return []
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    q = "SELECT * FROM experiments WHERE 1=1"
    params: list = []
    if claim_id is not None:
        q += " AND claim_id = ?"
        params.append(claim_id)
    if layer is not None:
        q += " AND layer = ?"
        params.append(layer)
    q += " ORDER BY id"
    rows = [dict(r) for r in conn.execute(q, params).fetchall()]
    conn.close()
    return rows


def print_results_table(rows: list[dict]) -> None:
    """Pretty-print experiment results."""
    if not rows:
        print("  (no results)")
        return
    header = f"  {'exp_id':<25} {'λ':>8}  {'N':>4}  {'verdict':>6}  {'A_min':>12}  {'ω_max':>12}  {'steps':>6}"
    print(header)
    print("  " + "-" * 80)
    for r in rows:
        lam = f"{r['lambda_val']:.4f}" if r["lambda_val"] is not None else "  ——  "
        A   = f"{r['A_min_global']:.4e}" if r["A_min_global"] is not None else "    ——    "
        om  = f"{r['omega_max_global']:.4e}" if r["omega_max_global"] is not None else "    ——    "
        status = "✅" if r["verdict"] == "PASS" else "❌"
        print(f"  {r['exp_id']:<25} {lam:>8}  {r['grid_N']:>4}  "
              f"{status} {r['verdict']:>4}  {A:>12}  {om:>12}  {r['n_steps']:>6}")


def conjecture_35_verdict(db_path: str = _DB_PATH) -> dict:
    """
    Compute overall Conjecture 3.5 verdict from logged results.

    Returns dict with keys:
      verdict         : 'PASS' or 'FAIL' or 'INCOMPLETE'
      n_runs          : total runs in DB for layer 2, route 1
      n_pass          : runs with verdict=PASS
      n_fail          : runs with verdict=FAIL
      lambda_c_found  : False if all PASS (Conjecture 3.5 confirmed)
      A_min_global    : minimum A_min across all runs
    """
    rows = query_results(db_path=db_path, layer=2)
    if not rows:
        return {"verdict": "INCOMPLETE", "n_runs": 0, "n_pass": 0, "n_fail": 0}

    n_pass = sum(1 for r in rows if r["verdict"] == "PASS")
    n_fail = sum(1 for r in rows if r["verdict"] == "FAIL")
    A_mins = [r["A_min_global"] for r in rows if r["A_min_global"] is not None]

    verdict = "PASS" if n_fail == 0 else "FAIL"

    return {
        "verdict"       : verdict,
        "n_runs"        : len(rows),
        "n_pass"        : n_pass,
        "n_fail"        : n_fail,
        "lambda_c_found": n_fail > 0,
        "A_min_global"  : float(min(A_mins)) if A_mins else None,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Smoke test (fast validation — used by test_lambda_sweep.py)
# ─────────────────────────────────────────────────────────────────────────────

def run_smoke_test(db_path: str | None = None, verbose: bool = False) -> dict:
    """
    Fast smoke test: 3 CCF λ values at N=32, t_end=0.05, max_steps=20.
    Used by test_lambda_sweep_2D.py. Does NOT write to main results.db unless
    db_path is explicitly provided.

    Returns:
        dict with 'verdict', 'n_runs', 'n_pass', 'results'
    """
    import tempfile

    tmp_db = db_path or os.path.join(tempfile.mkdtemp(), "smoke_test.db")

    smoke_lambdas = [1.2, 0.6057, 0.1]  # representative low/mid/high
    results = run_ccf_sweep(
        lambda_values=smoke_lambdas,
        N=32,
        t_end=0.05,
        max_steps=20,
        verbose=verbose,
        db_path=tmp_db,
    )

    n_pass = sum(1 for r in results if r["verdict"] == "PASS")
    verdict = "PASS" if n_pass == len(results) else "FAIL"

    return {
        "verdict": verdict,
        "n_runs" : len(results),
        "n_pass" : n_pass,
        "results": results,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="T-First 2D Lambda Sweep — PRD §7.3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python layer2/lambda_sweep_2D.py --smoke-test
  python layer2/lambda_sweep_2D.py --profile ccf --N 64 --t-end 0.5
  python layer2/lambda_sweep_2D.py --profile all --N 64 --verbose
  python layer2/lambda_sweep_2D.py --query
"""
    )
    parser.add_argument("--N",         type=int,   default=64,      help="Grid size")
    parser.add_argument("--t-end",     type=float, default=0.5,     help="End time")
    parser.add_argument("--fluid",     default="ideal",              help="Fluid type")
    parser.add_argument("--max-steps", type=int,   default=2000,    help="Max timesteps")
    parser.add_argument("--profile",   default="ccf",
                        choices=["ccf", "bous", "adv", "all"],      help="Profile type(s)")
    parser.add_argument("--smoke-test",action="store_true",          help="Quick smoke test")
    parser.add_argument("--verbose",   action="store_true",          help="Verbose output")
    parser.add_argument("--query",     action="store_true",          help="Print stored results")
    parser.add_argument("--db-path",   default=_DB_PATH,             help="Path to results.db")
    args = parser.parse_args()

    print("=" * 70)
    print("T-First 2D Lambda Sweep — PRD §7.3")
    print(f"  Grid: {args.N}²   fluid: {args.fluid}   t_end: {args.t_end}")
    print("=" * 70)

    if args.query:
        rows = query_results(db_path=args.db_path, layer=2)
        print(f"\nStored results in {args.db_path} (layer 2):")
        print_results_table(rows)
        v = conjecture_35_verdict(db_path=args.db_path)
        print(f"\nConjecture 3.5 verdict: {v['verdict']}")
        print(f"  {v['n_pass']}/{v['n_runs']} PASS, λ_c found: {v['lambda_c_found']}")
        sys.exit(0)

    if args.smoke_test:
        print("\nSmoke test (N=32, t_end=0.05, 3 λ values) ...")
        r = run_smoke_test(verbose=args.verbose)
        print(f"\nSmoke test: {r['verdict']} ({r['n_pass']}/{r['n_runs']} PASS)")
        sys.exit(0 if r["verdict"] == "PASS" else 1)

    kw = dict(N=args.N, t_end=args.t_end, fluid=args.fluid,
              max_steps=args.max_steps, verbose=args.verbose,
              db_path=args.db_path)

    all_results = []

    if args.profile in ("ccf", "all"):
        print("\n--- CCF Lambda Sweep ---")
        ccf = run_ccf_sweep(**kw)
        all_results.extend(ccf)

    if args.profile in ("bous", "all"):
        print("\n--- Boussinesq Lambda Sweep ---")
        bous = run_boussinesq_sweep(**kw)
        all_results.extend(bous)

    if args.profile in ("adv", "all"):
        print("\n--- Adversarial Minimum ---")
        adv = run_adversarial(**kw)
        all_results.append(adv)

    # Print summary table
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    header = f"  {'exp_id':<25} {'λ':>8}  {'verdict':>6}  {'A_min':>12}  {'steps':>6}"
    print(header)
    print("  " + "-" * 65)
    for r in all_results:
        lam = f"{r.get('lambda_val', 0):.4f}"
        A   = f"{r['A_min_global']:.4e}" if r.get('A_min_global') else "    ——    "
        status = "✅" if r["verdict"] == "PASS" else "❌"
        print(f"  {r.get('exp_id','?'):<25} {lam:>8}  {status} {r['verdict']:>4}  {A:>12}  {r['n_steps']:>6}")

    n_pass = sum(1 for r in all_results if r["verdict"] == "PASS")
    n_total = len(all_results)
    overall = "PASS" if n_pass == n_total else "FAIL"
    print(f"\nConjecture 3.5: {overall} ({n_pass}/{n_total} runs passed)")
    print(f"Results logged to: {args.db_path}")
    sys.exit(0 if overall == "PASS" else 1)
