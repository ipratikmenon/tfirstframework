"""
layer4/enstrophy_exponent.py — Regime Discriminator for H2 (Omega vs M scaling)
================================================================================
Context (§cf-s37, in progress in a parallel session; this module does not
touch the .tex file): whether H2 (localized Constantin-Fefferman depletion
at the shrinking self-similar scale r*(t)) holds unconditionally, or reduces
to the same open problem as sigma*-decay, hinges on how the global enstrophy
Omega(t) = ||omega(.,t)||^2_{L2} = integral_{T^3} |omega(x,t)|^2 dx grows
relative to M(t) = max_x |omega(x,t)| near strong vorticity growth:

  Regime A: Omega = O(1)      (bounded)            -> exponent p ~ 0
  Regime B: Omega ~ M^{1/2}   (self-similar Type-I) -> exponent p ~ 0.5

Measures the empirical exponent p in Omega(t) ~ M(t)^p during the existing
TG stretching run (the only run with real dynamic range in M, ~15x growth)
and, for comparison, the adversarial run (M stays ~flat there -> fit
expected to be uninformative).

House rules: spectral-only vorticity/enstrophy computation reusing existing
helpers; no hardcoded fluid properties; functions < 50 lines; immutable
(no input array mutated); every logged result carries claim_id + key_metric
+ verdict; seed always set/logged. Imports FROM layer3/route2_3D.py and
layer4/alignment_bridge.py WITHOUT modifying either.

Run: python layer4/enstrophy_exponent.py --ic tg --N 64 --n-steps 400 --verbose
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from datetime import datetime

import numpy as np
# ── path setup (matches layer4/alignment_bridge.py + depletion_diagnostic.py) ──
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

import route2_3D as r2  # noqa: E402 — Route 2 3D solver (reused, not reimplemented)
from layer4.alignment_bridge import (  # noqa: E402
    _vorticity_from_velocity,
    _validate_omega_grid,
    adv_vortex_tubes_ic,
)

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_SEED = 42
_IC_TAGS = {"tg": "TG", "adv": "ADV"}
_DB_MAX_RETRIES = 6         # sqlite "database is locked" retry policy
_DB_BACKOFF_BASE_S = 0.15
_PRIMARY_M_THRESHOLD = 5.0   # per task spec, mirrors depletion's M>10 idea
_FALLBACK_M_THRESHOLD = 3.0  # used if primary leaves too few points
_MIN_POINTS_FOR_PRIMARY = 6
_MIN_POINTS_FOR_FIT = 2      # a log-log line needs >= 2 points, period

# ── Core primitive: Omega and M from one vorticity snapshot ─────────────────

def enstrophy_and_M(
    wx: np.ndarray, wy: np.ndarray, wz: np.ndarray, dx: float,
) -> dict:
    """Omega = ||omega||_{L2}^2 = integral |omega|^2 dx over the periodic
    box [0,2pi)^3, and M = max|omega|.

    Direct-sum route: Omega = mean(|omega|^2) * Volume, Volume = (N*dx)^3
    (Riemann sum on the periodic grid's own cell volume dx^3=(2pi/N)^3,
    exact for the FFT's own quadrature; equivalent to spectral Parseval).

    Does not mutate wx, wy, wz. Raises ValueError via
    layer4.alignment_bridge._validate_omega_grid on shape/spacing issues.
    """
    N = _validate_omega_grid(wx, wy, wz, dx)
    mag_sq = wx ** 2 + wy ** 2 + wz ** 2
    volume = (N * dx) ** 3
    Omega = float(np.mean(mag_sq) * volume)
    M = float(np.sqrt(np.max(mag_sq)))
    return {"Omega": Omega, "M": M}

# ── Regression helper: fit p in Omega ~ M^p via log-log linear regression ──

def loglog_power_fit(M_vals: np.ndarray, Omega_vals: np.ndarray) -> dict:
    """Fit log(Omega)=p*log(M)+c by OLS. Returns dict{p, intercept,
    r_squared, n_points}. Raises ValueError if fewer than 2 usable
    (M>0, Omega>0, finite) points supplied. Immutable."""
    M_arr = np.asarray(M_vals, dtype=float)
    O_arr = np.asarray(Omega_vals, dtype=float)
    usable = np.isfinite(M_arr) & np.isfinite(O_arr) & (M_arr > 0) & (O_arr > 0)
    n = int(np.sum(usable))
    if n < _MIN_POINTS_FOR_FIT:
        raise ValueError(
            f"loglog_power_fit: need >= {_MIN_POINTS_FOR_FIT} positive-finite "
            f"points, got {n}"
        )
    log_M = np.log(M_arr[usable])
    log_O = np.log(O_arr[usable])
    p, intercept = np.polyfit(log_M, log_O, 1)

    pred = p * log_M + intercept
    ss_res = float(np.sum((log_O - pred) ** 2))
    ss_tot = float(np.sum((log_O - np.mean(log_O)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 1e-300 else float("nan")

    return {
        "p": float(p), "intercept": float(intercept),
        "r_squared": float(r_squared), "n_points": n,
    }

def select_high_M_mask(
    M_vals: np.ndarray,
    primary: float = _PRIMARY_M_THRESHOLD,
    fallback: float = _FALLBACK_M_THRESHOLD,
    min_points: int = _MIN_POINTS_FOR_PRIMARY,
) -> tuple[float, np.ndarray, bool]:
    """Boolean mask for M_vals > threshold; try `primary` first, fall back
    to `fallback` if that leaves fewer than min_points. Returns
    (threshold_used, mask, used_fallback). Immutable."""
    M_arr = np.asarray(M_vals, dtype=float)
    mask_primary = M_arr > primary
    if int(np.sum(mask_primary)) >= min_points:
        return primary, mask_primary, False
    mask_fallback = M_arr > fallback
    return fallback, mask_fallback, True

def select_growth_phase_mask(M_vals: np.ndarray) -> np.ndarray:
    """Boolean mask for the stretching phase: snapshots strictly before M's
    peak index (argmax) — the physically relevant regime near a putative
    singularity. Immutable."""
    M_arr = np.asarray(M_vals, dtype=float)
    if M_arr.size == 0:
        return np.zeros(0, dtype=bool)
    peak_idx = int(np.argmax(M_arr))
    mask = np.zeros(M_arr.shape, dtype=bool)
    mask[:peak_idx] = True
    return mask

def _safe_fit(M_vals: np.ndarray, Omega_vals: np.ndarray, mask: np.ndarray) -> dict:
    """loglog_power_fit restricted to `mask`, degrading to a NaN-filled
    dict when there aren't enough points — report honestly, no forced p."""
    n_selected = int(np.sum(mask))
    if n_selected < _MIN_POINTS_FOR_FIT:
        return {"p": float("nan"), "intercept": float("nan"),
                "r_squared": float("nan"), "n_points": n_selected}
    try:
        return loglog_power_fit(np.asarray(M_vals)[mask], np.asarray(Omega_vals)[mask])
    except ValueError:
        return {"p": float("nan"), "intercept": float("nan"),
                "r_squared": float("nan"), "n_points": n_selected}

# ── Production experiment: run Route 2 3D solver, track Omega vs M ─────────

def run_enstrophy_exponent_experiment(
    ic: str = "tg",
    N: int = 64,
    n_steps: int = 400,
    nu: float | None = 1e-3,
    eps_param: float = 0.1,
    record_every: int = 10,
    seed: int = _DEFAULT_SEED,
    cfl_safety: float = 0.4,
    verbose: bool = False,
    db_path: str | None = None,
    exp_id_override: str | None = None,
    claim_id_override: str | None = None,
    adv_amp: float | None = None,
) -> dict:
    """Run the existing Route 2 3D solver and track (M, Omega) to fit the
    empirical exponent p in Omega(t) ~ M(t)^p.

    Mirrors run_bridge_experiment / run_depletion_experiment's solver setup
    exactly (same make_solver / set_ic / compute_cfl_dt / solver_step
    sequence, same IC constructors, same seed/record cadence) so this run
    is directly comparable to EXP-L4-DEPLETION-TG-001/ADV-001 — no solver
    logic is reimplemented here.
    """
    if ic not in _IC_TAGS:
        raise ValueError(f"Unknown ic='{ic}'; use 'tg' or 'adv'.")
    if n_steps <= 0:
        raise ValueError(f"n_steps must be > 0, got {n_steps}")

    if nu is None:
        nu = r2._NU_DEFAULT

    # Seed set and logged for reproducibility provenance (identical
    # rationale to alignment_bridge.run_bridge_experiment).
    _rng = np.random.default_rng(seed)

    solver = r2.make_solver(N=N, nu=nu, eps_param=eps_param)
    if ic == "tg":
        ic_dict = r2.taylor_green_ic(N)
    else:
        adv_kwargs = {"seed": seed}
        if adv_amp is not None:
            adv_kwargs["amp"] = adv_amp
        ic_dict = adv_vortex_tubes_ic(N, **adv_kwargs)
    solver = r2.set_ic(solver, ic_dict)

    dx = 2.0 * np.pi / N
    kx, ky, kz = solver["kx"], solver["ky"], solver["kz"]

    timeseries: list[dict] = []
    t0 = time.perf_counter()

    for step_idx in range(n_steps):
        dt = r2.compute_cfl_dt(solver, safety=cfl_safety)
        solver, _diag = r2.solver_step(solver, dt)

        is_last = step_idx == n_steps - 1
        if step_idx % record_every == 0 or is_last:
            wx, wy, wz = _vorticity_from_velocity(
                solver["u"], solver["v"], solver["w"], kx, ky, kz
            )
            em = enstrophy_and_M(wx, wy, wz, dx)
            record = {"step": step_idx, "t": solver["t"], **em}
            timeseries.append(record)
            if verbose:
                print(
                    f"  step {step_idx:4d}  t={record['t']:.4f}  "
                    f"M={em['M']:.4e}  Omega={em['Omega']:.4e}"
                )

    wall = time.perf_counter() - t0
    result = _summarize_run(timeseries, ic, N, n_steps, seed, nu, eps_param,
                             solver["t"], wall, exp_id_override, claim_id_override)
    result["timeseries"] = timeseries

    if verbose:
        print(
            f"\n  Enstrophy exponent ({ic}, N={N}, n_steps={n_steps})  "
            f"M={result['M_range']}  Omega={result['Omega_range']}\n"
            f"  p_fit_all={result['p_fit_all']:.4f} (thr={result['threshold_used']}, "
            f"n={result['n_points_used']})  "
            f"p_fit_growth_phase={result['p_fit_growth_phase']:.4f} "
            f"(n={result['n_points_growth_phase']}, r2={result['r_squared']:.4f})  "
            f"verdict={result['verdict']}\n  note: {result['note']}"
        )

    if db_path is not None:
        log_result(db_path, result)

    return result

def _summarize_run(
    timeseries: list[dict], ic: str, N: int, n_steps: int, seed: int,
    nu: float, eps_param: float, t_final: float, wall: float,
    exp_id_override: str | None, claim_id_override: str | None,
) -> dict:
    """Fold a raw (step,t,M,Omega) timeseries into the logged summary dict."""
    M_vals = np.array([rec["M"] for rec in timeseries], dtype=float)
    Omega_vals = np.array([rec["Omega"] for rec in timeseries], dtype=float)

    threshold_used, mask_high, used_fallback = select_high_M_mask(M_vals)
    fit_all = _safe_fit(M_vals, Omega_vals, mask_high)

    mask_growth = select_growth_phase_mask(M_vals)
    fit_growth = _safe_fit(M_vals, Omega_vals, mask_growth)

    note = (
        f"used fallback M>{_FALLBACK_M_THRESHOLD} threshold "
        f"({int(np.sum(M_vals > _PRIMARY_M_THRESHOLD))} points at "
        f"M>{_PRIMARY_M_THRESHOLD} < min {_MIN_POINTS_FOR_PRIMARY})"
        if used_fallback else
        f"used primary M>{_PRIMARY_M_THRESHOLD} threshold"
    )

    verdict = "PASS" if (
        np.isfinite(fit_growth["p"]) and np.isfinite(fit_growth["r_squared"])
    ) else "FAIL"

    tag = _IC_TAGS[ic]
    exp_id = exp_id_override or f"EXP-L4-ENSTROPHY-EXPONENT-{tag}-001"
    claim_id = claim_id_override or f"enstrophy-exponent-{ic}"

    M_min, M_max = (float(np.min(M_vals)), float(np.max(M_vals))) if M_vals.size else (0.0, 0.0)
    O_min, O_max = (float(np.min(Omega_vals)), float(np.max(Omega_vals))) if Omega_vals.size else (0.0, 0.0)

    return {
        "exp_id": exp_id, "claim_id": claim_id,
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic, "N": N, "n_steps": n_steps, "seed": seed, "nu": nu,
        "eps_param": eps_param, "t_final": t_final,
        "p_fit_all": fit_all["p"], "p_fit_growth_phase": fit_growth["p"],
        "threshold_used": threshold_used,
        "n_points_used": fit_all["n_points"],
        "n_points_growth_phase": fit_growth["n_points"],
        "r_squared": fit_growth["r_squared"], "r_squared_all": fit_all["r_squared"],
        "M_range": f"[{M_min:.4e},{M_max:.4e}]",
        "Omega_range": f"[{O_min:.4e},{O_max:.4e}]",
        "verdict": verdict, "wall_time_s": wall, "note": note,
        "key_metric": (
            f"p_fit_growth_phase={fit_growth['p']:.4f} "
            f"r_squared={fit_growth['r_squared']:.4f} "
            f"n_points={fit_growth['n_points']}"
        ),
    }

# ── Results logging (layer4 table pattern) — retry-on-locked backoff ───────

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path, timeout=30.0)
    conn.execute(
        "CREATE TABLE IF NOT EXISTS layer4_enstrophy_exponent ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, exp_id TEXT UNIQUE NOT NULL, "
        "claim_id TEXT NOT NULL, timestamp TEXT NOT NULL, ic TEXT, "
        "N INTEGER, n_steps INTEGER, seed INTEGER, nu REAL, eps_param REAL, "
        "t_final REAL, p_fit_all REAL, p_fit_growth_phase REAL, "
        "threshold_used REAL, n_points_used INTEGER, "
        "n_points_growth_phase INTEGER, r_squared REAL, r_squared_all REAL, "
        "M_range TEXT, Omega_range TEXT, verdict TEXT NOT NULL, "
        "wall_time_s REAL, key_metric TEXT, note TEXT)"
    )
    conn.execute(
        "CREATE TABLE IF NOT EXISTS layer4_enstrophy_exponent_timeseries ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, exp_id TEXT NOT NULL, "
        "step INTEGER, t REAL, M REAL, Omega REAL)"
    )
    conn.commit()
    return conn

def _write_result(db_path: str, result: dict) -> None:
    """Single-attempt DB write (summary UPSERT + append-only timeseries)."""
    conn = _ensure_db(db_path)
    try:
        cols = [
            "exp_id", "claim_id", "timestamp", "ic", "N", "n_steps", "seed",
            "nu", "eps_param", "t_final", "p_fit_all", "p_fit_growth_phase",
            "threshold_used", "n_points_used", "n_points_growth_phase",
            "r_squared", "r_squared_all", "M_range", "Omega_range",
            "verdict", "wall_time_s", "key_metric", "note",
        ]
        values = tuple(result.get(c) for c in cols)
        ph = ", ".join("?" for _ in cols)
        conn.execute(
            f"INSERT OR REPLACE INTO layer4_enstrophy_exponent ({', '.join(cols)}) "
            f"VALUES ({ph})",
            values,
        )
        for rec in result.get("timeseries", []):
            conn.execute(
                "INSERT INTO layer4_enstrophy_exponent_timeseries "
                "(exp_id, step, t, M, Omega) VALUES (?, ?, ?, ?, ?)",
                (result["exp_id"], rec["step"], rec["t"], rec["M"], rec["Omega"]),
            )
        conn.commit()
    finally:
        conn.close()

def log_result(
    db_path: str, result: dict,
    max_retries: int = _DB_MAX_RETRIES,
    backoff_base_s: float = _DB_BACKOFF_BASE_S,
) -> None:
    """Append one enstrophy-exponent experiment (+ timeseries) to results.db.
    Retries with short exponential backoff on sqlite3.OperationalError
    ("database is locked") since other scripts may write concurrently;
    re-raises immediately on any other error or once retries exhaust."""
    last_err: Exception | None = None
    for attempt in range(max_retries):
        try:
            _write_result(db_path, result)
            return
        except sqlite3.OperationalError as exc:
            if "locked" not in str(exc).lower() or attempt == max_retries - 1:
                raise
            last_err = exc
            time.sleep(backoff_base_s * (2 ** attempt))
    if last_err is not None:  # pragma: no cover — defensive, unreachable
        raise last_err

# ── CLI entry point ──────────────────────────────────────────────────────

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="H2 regime discriminator: fit Omega(t) ~ M(t)^p"
    )
    parser.add_argument("--ic", choices=["tg", "adv"], default="tg")
    parser.add_argument("--N", type=int, default=64)
    parser.add_argument("--n-steps", type=int, default=400)
    parser.add_argument("--record-every", type=int, default=10)
    parser.add_argument("--seed", type=int, default=_DEFAULT_SEED)
    parser.add_argument("--nu", type=float, default=None)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db", default=_DEFAULT_DB)
    args = parser.parse_args()

    r = run_enstrophy_exponent_experiment(
        ic=args.ic, N=args.N, n_steps=args.n_steps,
        record_every=args.record_every, seed=args.seed, nu=args.nu,
        eps_param=args.eps, verbose=args.verbose, db_path=args.db,
    )
    print(f"\nVerdict: {r['verdict']}  key_metric: {r['key_metric']}")

if __name__ == "__main__":
    main()
