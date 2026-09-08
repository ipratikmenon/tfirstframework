"""
layer4/kcorr_reynolds_sweep.py — Reynolds Sweep for the K-Absorption
Correlation Diagnostic (c_K)
================================================================
S34 measured c_K ≈ 1 (decorrelated) at a single viscosity, nu=1e-3. That
is consistent with the refined K-absorption conjecture c_K <= C/𝓛, but
cannot distinguish "flat O(1)" from "decreasing like 1/log(Re)" — both
predict c_K = O(1) at one moderate Reynolds number.

This module re-runs the c_K diagnostic across three viscosities
(nu = 1e-3, 5e-4, 2.5e-4 — a 4x Reynolds range) for two ICs (Taylor-Green,
adversarial anti-parallel vortex tubes), and fits median c_K against
log(1/nu): c_K_median(nu) ~ slope*log(1/nu) + intercept.
  slope < 0  -> c_K decreases with Re -> supports c_K <= C/𝓛
  slope ~ 0  -> c_K flat O(1)         -> supports plain O(1) absorbability
  slope > 0  -> c_K increases with Re -> warning sign against K-absorption

Because nu decreases while N, n_steps stay fixed (N=64, 400 steps), each
run also carries a resolution-sanity check: the spectral tail fraction
(energy in the top 1/3 of dealiased-retained wavenumbers / total energy)
at the final step. A large tail fraction flags the run as under-resolved
(resolution_warning=True) rather than silently trusting the c_K trend.

Reuses (does not reimplement): layer3/route2_3D.py (solver, ICs, dealias
mask, grid), layer4/alignment_bridge.py (adv_vortex_tubes_ic,
_vorticity_from_velocity), layer4/k_correlation.py
(k_correlation_diagnostics). Only new orchestration lives here: the
Reynolds run matrix, the tail-fraction diagnostic, the log-linear
regression, and the sweep-level results table.

House rules: spectral (FFT) derivatives only; all functions pure /
immutable; every logged result carries claim_id + key_metric + verdict;
random seed always set and logged.

Run:
  python layer4/kcorr_reynolds_sweep.py --verbose
  python layer4/kcorr_reynolds_sweep.py --N 16 --n-steps 20 --verbose   # smoke
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

import route2_3D as r2  # noqa: E402 — Route 2 3D solver (reused, not reimplemented)
from layer4.alignment_bridge import (  # noqa: E402
    adv_vortex_tubes_ic,
    _vorticity_from_velocity,
)
from layer4.k_correlation import k_correlation_diagnostics  # noqa: E402

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_NUS = (1.0e-3, 5.0e-4, 2.5e-4)
_DEFAULT_ICS = ("tg", "adv")
_DEFAULT_N = 64
_DEFAULT_N_STEPS = 400
_DEFAULT_SEED = 42
_DEFAULT_RECORD_EVERY = 10
_KCORR_PASS_MEDIAN_CK = 10.0
_TAIL_FRACTION_WARN = 1.0e-3
_IC_TAGS = {"tg": "TG", "adv": "ADV"}


# --- Resolution-sanity diagnostic: spectral tail fraction ---

def spectral_tail_fraction(u: np.ndarray, v: np.ndarray, w: np.ndarray, N: int) -> float:
    """Fraction of kinetic energy in the top 1/3 of retained wavenumbers.

    "Retained" = the 2/3-rule dealiasing box the solver already uses
    (route2_3D.make_dealias_mask_3D: |k_i| <= N//3 per axis). Within that
    box, the "top 1/3" shell is the isotropic annulus
    |k| in ((2/3)*k_max_retained, k_max_retained], k_max_retained = N//3.
    fraction = sum(energy in that shell) / sum(energy over the whole domain).
    Large fraction => energy piling up near the grid cutoff => under-resolved
    at this (N, nu). Does not mutate u, v, w.
    """
    kx, ky, kz, k2 = r2.make_grid_3D(N)
    dealias = r2.make_dealias_mask_3D(N).astype(bool)
    k_mag = np.sqrt(k2)
    k_max_retained = N // 3
    tail_mask = dealias & (k_mag > (2.0 / 3.0) * k_max_retained)

    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)
    energy_density = np.abs(u_hat) ** 2 + np.abs(v_hat) ** 2 + np.abs(w_hat) ** 2

    total = float(np.sum(energy_density))
    if total <= 0.0:
        return 0.0
    tail = float(np.sum(energy_density[tail_mask]))
    return tail / total


# --- Log-linear regression helper: median c_K vs log(1/nu) ---

def linregress_slope(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """OLS slope/intercept for y = slope*x + intercept (dependency-free, no scipy).

    Raises ValueError if fewer than 2 points are given or all x are
    identical (slope undefined). Does not mutate x, y.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if x_arr.size != y_arr.size:
        raise ValueError(f"x and y must have the same length; got {x_arr.size}, {y_arr.size}")
    if x_arr.size < 2:
        raise ValueError("linregress_slope requires at least 2 points")

    x_mean = float(np.mean(x_arr))
    y_mean = float(np.mean(y_arr))
    denom = float(np.sum((x_arr - x_mean) ** 2))
    if denom <= 0.0:
        raise ValueError("linregress_slope: all x values are identical; slope undefined")

    slope = float(np.sum((x_arr - x_mean) * (y_arr - y_mean)) / denom)
    intercept = y_mean - slope * x_mean
    return slope, intercept


def trend_label(slope: float, threshold: float = 0.05) -> str:
    """Human-readable classification of a c_K-vs-log(1/nu) regression slope."""
    if not np.isfinite(slope):
        return "undefined"
    if slope < -threshold:
        return "decreasing"
    return "increasing" if slope > threshold else "flat"


# --- Run matrix ---

def build_run_matrix(
    nus: tuple[float, ...] = _DEFAULT_NUS,
    ics: tuple[str, ...] = _DEFAULT_ICS,
    N: int = _DEFAULT_N,
    n_steps: int = _DEFAULT_N_STEPS,
    seed: int = _DEFAULT_SEED,
    record_every: int = _DEFAULT_RECORD_EVERY,
    exp_id_start: int = 1,
    claim_id_suffix: str = "",
) -> list[dict]:
    """Build the ordered list of (ic, nu) run specs for the Reynolds sweep.

    Outer loop over ics, inner loop over nus (decreasing nu = increasing
    Re). exp_id = EXP-L4-KCORR-NU-{TAG}-{exp_id_start, exp_id_start+1, ...}
    (default exp_id_start=1 -> 001,002,003, matching the original N=64
    sweep); one claim_id per ic (kcorr-reynolds-{ic}{claim_id_suffix})
    shared across its viscosities, since the claim is about the trend
    across Re, not any single run. claim_id_suffix lets a separate,
    higher-resolution re-run of the same (ic, nu) points (e.g. the N=128
    resolved TG re-run) log under a distinct claim_id — e.g.
    claim_id_suffix="-resolved" -> "kcorr-reynolds-tg-resolved" — without
    colliding with or overwriting the original N=64 sweep's claim.
    """
    matrix: list[dict] = []
    for ic in ics:
        if ic not in _IC_TAGS:
            raise ValueError(f"Unknown ic='{ic}'; use 'tg' or 'adv'.")
        tag = _IC_TAGS[ic]
        for i, nu in enumerate(nus):
            matrix.append({
                "ic": ic,
                "nu": float(nu),
                "N": N,
                "n_steps": n_steps,
                "seed": seed,
                "record_every": record_every,
                "exp_id": f"EXP-L4-KCORR-NU-{tag}-{exp_id_start + i:03d}",
                "claim_id": f"kcorr-reynolds-{ic}{claim_id_suffix}",
            })
    return matrix


# --- Single production run ---

def run_single_kcorr_sweep(
    ic: str,
    nu: float,
    N: int = _DEFAULT_N,
    n_steps: int = _DEFAULT_N_STEPS,
    seed: int = _DEFAULT_SEED,
    record_every: int = _DEFAULT_RECORD_EVERY,
    verbose: bool = False,
) -> dict:
    """Run the Route 2 3D solver at one (ic, nu) point and summarize c_K.

    Reuses route2_3D (solver stepping, ICs, grid), alignment_bridge
    (adv_vortex_tubes_ic, _vorticity_from_velocity), and
    k_correlation.k_correlation_diagnostics — no solver/diagnostic logic
    is reimplemented here. Returns cK_median/max/final, f_ang_median,
    split_residual_median, tail_fraction, resolution_warning, verdict,
    key_metric. Verdict PASS iff every recorded c_K is finite AND
    cK_median <= 10.
    """
    if ic not in _IC_TAGS:
        raise ValueError(f"Unknown ic='{ic}'; use 'tg' or 'adv'.")
    if n_steps <= 0:
        raise ValueError(f"n_steps must be > 0, got {n_steps}")

    solver = r2.make_solver(N=N, nu=nu, eps_param=0.1)
    if ic == "tg":
        ic_dict = r2.taylor_green_ic(N)
    else:
        ic_dict = adv_vortex_tubes_ic(N, seed=seed)
    solver = r2.set_ic(solver, ic_dict)

    dx = 2.0 * np.pi / N
    kx, ky, kz = solver["kx"], solver["ky"], solver["kz"]

    c_K_vals: list[float] = []
    f_ang_vals: list[float] = []
    split_residual_vals: list[float] = []

    t0 = time.perf_counter()
    for step_idx in range(n_steps):
        dt = r2.compute_cfl_dt(solver, safety=0.4)
        solver, _diag = r2.solver_step(solver, dt)

        is_last = step_idx == n_steps - 1
        if step_idx % record_every == 0 or is_last:
            wx, wy, wz = _vorticity_from_velocity(
                solver["u"], solver["v"], solver["w"], kx, ky, kz
            )
            try:
                kc = k_correlation_diagnostics(wx, wy, wz, dx)
            except ValueError:
                kc = {"c_K": 0.0, "f_ang": 0.0, "split_residual": 0.0}
            c_K_vals.append(kc["c_K"])
            f_ang_vals.append(kc["f_ang"])
            split_residual_vals.append(kc["split_residual"])
            if verbose:
                print(
                    f"    [{ic} nu={nu:.2e}] step {step_idx:4d} t={solver['t']:.4f} "
                    f"c_K={kc['c_K']:.4e} f_ang={kc['f_ang']:.4e}"
                )
    wall = time.perf_counter() - t0

    tail_fraction = spectral_tail_fraction(solver["u"], solver["v"], solver["w"], N)
    resolution_warning = bool(tail_fraction > _TAIL_FRACTION_WARN)

    cK_arr = np.array(c_K_vals, dtype=float)
    all_finite = bool(np.all(np.isfinite(cK_arr))) if cK_arr.size else False
    cK_median = float(np.median(cK_arr)) if cK_arr.size else float("nan")
    cK_max = float(np.max(cK_arr)) if cK_arr.size else float("nan")
    cK_final = float(cK_arr[-1]) if cK_arr.size else float("nan")
    f_ang_median = float(np.median(f_ang_vals)) if f_ang_vals else float("nan")
    split_residual_median = float(np.median(split_residual_vals)) if split_residual_vals else float("nan")

    verdict = "PASS" if (all_finite and np.isfinite(cK_median) and cK_median <= _KCORR_PASS_MEDIAN_CK) else "FAIL"

    return {
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic,
        "nu": float(nu),
        "N": N,
        "n_steps": n_steps,
        "seed": seed,
        "cK_median": cK_median,
        "cK_max": cK_max,
        "cK_final": cK_final,
        "f_ang_median": f_ang_median,
        "split_residual_median": split_residual_median,
        "tail_fraction": tail_fraction,
        "resolution_warning": resolution_warning,
        "verdict": verdict,
        "wall_time_s": wall,
        "key_metric": f"cK_median={cK_median:.4f}",
    }


# --- Results logging (layer4-module table pattern; see k_correlation.py) ---

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path, timeout=10.0)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_kcorr_sweep (
            id                     INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id                 TEXT UNIQUE NOT NULL,
            claim_id               TEXT NOT NULL,
            timestamp              TEXT NOT NULL,
            ic                     TEXT,
            nu                     REAL,
            N                      INTEGER,
            n_steps                INTEGER,
            seed                   INTEGER,
            cK_median              REAL,
            cK_max                 REAL,
            cK_final               REAL,
            f_ang_median           REAL,
            split_residual_median  REAL,
            tail_fraction          REAL,
            resolution_warning     INTEGER,
            verdict                TEXT NOT NULL,
            wall_time_s            REAL,
            key_metric             TEXT
        )
    """)
    conn.commit()
    return conn


def log_sweep_result(
    db_path: str,
    row: dict,
    max_retries: int = 6,
    initial_backoff_s: float = 0.5,
) -> None:
    """Append one Reynolds-sweep run to results.db (UPSERT keyed on exp_id).

    Other layer4 modules may be logging to the same results.db concurrently
    (separate sessions/agents). SQLite raises sqlite3.OperationalError
    ("database is locked") under write contention; retry with exponential
    backoff (0.5s, 1s, 2s, 4s, 8s, 16s by default) before giving up and
    re-raising, rather than silently dropping the result.
    """
    cols = [
        "exp_id", "claim_id", "timestamp", "ic", "nu", "N", "n_steps", "seed",
        "cK_median", "cK_max", "cK_final", "f_ang_median", "split_residual_median",
        "tail_fraction", "resolution_warning", "verdict", "wall_time_s", "key_metric",
    ]
    values = tuple(
        int(row[c]) if c == "resolution_warning" else row.get(c) for c in cols
    )
    ph = ", ".join("?" for _ in cols)

    backoff = initial_backoff_s
    last_err: sqlite3.OperationalError | None = None
    for attempt in range(max_retries + 1):
        try:
            conn = _ensure_db(db_path)
            try:
                conn.execute(
                    f"INSERT OR REPLACE INTO layer4_kcorr_sweep ({', '.join(cols)}) "
                    f"VALUES ({ph})",
                    values,
                )
                conn.commit()
            finally:
                conn.close()
            return
        except sqlite3.OperationalError as err:
            if "locked" not in str(err).lower() or attempt == max_retries:
                raise
            last_err = err
            time.sleep(backoff)
            backoff *= 2.0
    if last_err is not None:  # pragma: no cover — unreachable given the raise above
        raise last_err


# --- Full sweep orchestration ---

def run_reynolds_sweep(
    nus: tuple[float, ...] = _DEFAULT_NUS,
    ics: tuple[str, ...] = _DEFAULT_ICS,
    N: int = _DEFAULT_N,
    n_steps: int = _DEFAULT_N_STEPS,
    seed: int = _DEFAULT_SEED,
    record_every: int = _DEFAULT_RECORD_EVERY,
    verbose: bool = False,
    db_path: str | None = _DEFAULT_DB,
    exp_id_start: int = 1,
    claim_id_suffix: str = "",
) -> dict:
    """Run the full Reynolds sweep, log every run, fit the c_K trend.

    exp_id_start / claim_id_suffix are forwarded to build_run_matrix (see
    its docstring) so a higher-resolution re-run of the same (ic, nu)
    points — e.g. the N=128 resolved TG re-run, exp_id_start=101,
    claim_id_suffix="-resolved" — logs under fresh exp_ids/claim_id
    instead of colliding with the original sweep.

    Returns {'runs': [...], 'regressions': {ic: {'slope','intercept','trend'}}}.
    """
    matrix = build_run_matrix(
        nus=nus, ics=ics, N=N, n_steps=n_steps, seed=seed, record_every=record_every,
        exp_id_start=exp_id_start, claim_id_suffix=claim_id_suffix,
    )

    runs: list[dict] = []
    for spec in matrix:
        if verbose:
            print(f"\n== {spec['exp_id']} ({spec['ic']}, nu={spec['nu']:.2e}) ==")
        result = run_single_kcorr_sweep(
            ic=spec["ic"], nu=spec["nu"], N=spec["N"], n_steps=spec["n_steps"],
            seed=spec["seed"], record_every=spec["record_every"], verbose=verbose,
        )
        row = {**result, "exp_id": spec["exp_id"], "claim_id": spec["claim_id"]}
        runs.append(row)
        if db_path is not None:
            log_sweep_result(db_path, row)
        if verbose:
            print(f"  -> {row['verdict']}  {row['key_metric']}  tail_fraction={row['tail_fraction']:.2e}"
                  f"  resolution_warning={row['resolution_warning']}")

    regressions: dict[str, dict] = {}
    for ic in ics:
        ic_runs = [r for r in runs if r["ic"] == ic]
        x = np.log(1.0 / np.array([r["nu"] for r in ic_runs], dtype=float))
        y = np.array([r["cK_median"] for r in ic_runs], dtype=float)
        slope, intercept = linregress_slope(x, y)
        regressions[ic] = {
            "slope": slope,
            "intercept": intercept,
            "trend": trend_label(slope),
        }
        if verbose:
            print(f"\n  regression [{ic}]: slope={slope:.4f} intercept={intercept:.4f} "
                  f"trend={regressions[ic]['trend']}")

    return {"runs": runs, "regressions": regressions}


# --- CLI entry point ---

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Reynolds sweep for the K-absorption c_K diagnostic")
    parser.add_argument("--N", type=int, default=_DEFAULT_N)
    parser.add_argument("--n-steps", type=int, default=_DEFAULT_N_STEPS)
    parser.add_argument("--record-every", type=int, default=_DEFAULT_RECORD_EVERY)
    parser.add_argument("--seed", type=int, default=_DEFAULT_SEED)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db", default=_DEFAULT_DB)
    args = parser.parse_args()

    result = run_reynolds_sweep(
        N=args.N, n_steps=args.n_steps, record_every=args.record_every,
        seed=args.seed, verbose=args.verbose, db_path=args.db,
    )
    for ic, reg in result["regressions"].items():
        print(f"[{ic}] slope={reg['slope']:.4f} trend={reg['trend']}")


if __name__ == "__main__":
    main()
