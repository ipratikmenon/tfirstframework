"""
layer4/sigma_cycles.py — Alignment/Disorder Cycle Analysis (Cascade Diagnostic)
================================================================================
Motivation (condensed). sigma_star(t) is a scale-invariant "coherence
deficit" of the vorticity direction field near the point of maximal
vorticity (layer4/alignment_bridge.py: bridge_diagnostics -> 'sigma_star').
Low = locally aligned/coherent vortex lines near the candidate blow-up
site; high = geometrically disordered. Existing runs show sigma_star
rising, then partially receding, as max-vorticity M(t) approaches its peak.

HYPOTHESIS UNDER TEST (not assumed true): turbulence near a candidate
singularity may proceed through repeated ALIGN -> STRETCH -> AMPLIFY ->
DESTABILIZE -> DISORDER -> DISSIPATE -> REALIGN cycles. Tao's "fluid
computer" cascade construction for finite-time blow-up requires successive
cycles at GEOMETRICALLY SHRINKING period and scale (sum of periods finite).
This module tests that, on the program's own data: shrinking periods +
minima marching toward 0 => cascade-leaning (disorder barrier threatened);
O(1)/growing periods + bounded-away-from-0 minima => persistent turbulence
without cascade; fewer than 3 complete cycles => genuinely INSUFFICIENT to
tell — the EXPECTED, HONEST outcome on the short runs logged so far.

THIS IS A FALSIFICATION TEST OF THE PROGRAM'S OWN THESIS. A cascade-leaning
result undermines the disorder barrier and is a valid, reportable,
PASS-worthy scientific outcome — "PASS" means "ran correctly and produced
an honest verdict", never "the thesis was confirmed". cascade_verdict has
no value that is inherently good or bad for the program: CASCADE-LEANING,
PERSISTENT, INSUFFICIENT-CYCLES, and AMBIGUOUS are all legitimate outcomes.

Late-time-growth caveat (definitional artifact): LATE-TIME sigma_star
growth — after M(t) peaks and the flow decays — is probably a DEFINITIONAL
ARTIFACT: as M falls, the high-vorticity set {|omega| >= M/4} used inside
bridge_diagnostics becomes unselective and picks up decaying debris rather
than genuine coherent structure. This module therefore always reports
BOTH the full trajectory and the GROWTH PHASE ONLY (steps up to and
including the time of max M), clearly labelled; growth-phase numbers are
the trustworthy ones.

See layer4/sigma_cycles_run.py for the new long, densely-sampled unforced
production runs (Part 2) and its unforced-flow-limitation docstring.

Reused, not reimplemented: alignment_bridge.py's bridge_diagnostics
(sigma_star is computed exclusively there); kcorr_reynolds_sweep.py's
linregress_slope. New here: peak-finding/cycle detection, cycle-ratio and
cascade-verdict logic, the growth-phase split, and this module's tables.

House rules: pure/immutable functions; every logged result carries
claim_id + key_metric + verdict; retry-with-backoff on a locked results.db
(depletion_diagnostic.py pattern).

Run: python layer4/sigma_cycles.py --analyze-existing --verbose
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

from layer4.kcorr_reynolds_sweep import linregress_slope  # noqa: E402

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_ALLOWED_TIMESERIES_TABLES = {
    "layer4_bridge_timeseries",
    "layer4_sigma_cycles_timeseries",
}
_DB_MAX_RETRIES = 6
_DB_BACKOFF_BASE_S = 0.5
_CASCADE_PERIOD_RATIO_THRESHOLD = 0.8
_MIN_COMPLETE_CYCLES = 3

_EXISTING_EXPERIMENTS = [
    "EXP-L4-BRIDGE-TG-001", "EXP-L4-BRIDGE-SH-001",
    "EXP-L4-BRIDGE-TG-002", "EXP-L4-BRIDGE-SH-002",
    "EXP-L4-BRIDGE-ADV-001", "EXP-L4-BRIDGE-ADV-002",
]


# Peak finding: hand-rolled (scipy is unused elsewhere in layer4, so this
# avoids introducing scipy.signal.find_peaks as a new layer4 dependency) —
# light moving-average smoothing + a prominence filter.

def _moving_average(values: np.ndarray, window: int) -> np.ndarray:
    """Edge-padded moving average; same length as input. Does not mutate values."""
    values = np.asarray(values, dtype=float)
    if window <= 1 or values.size == 0:
        return values.copy()
    w = int(window)
    pad_left = w // 2
    pad_right = w - 1 - pad_left
    padded = np.pad(values, (pad_left, pad_right), mode="edge")
    kernel = np.ones(w) / w
    return np.convolve(padded, kernel, mode="valid")


def _find_local_extrema(values: np.ndarray, prominence_ratio: float) -> list:
    """Alternating local minima/maxima of `values`, small wiggles collapsed.
    Candidates are strict local min/max of the raw sequence; adjacent pairs
    whose value difference is below prominence_ratio * (max - min) are
    iteratively deleted (noise, not genuine turning points). Does not
    mutate values.
    """
    values = np.asarray(values, dtype=float)
    n = values.size
    if n < 3:
        return []
    data_range = float(np.max(values) - np.min(values))
    if data_range <= 0.0:
        return []
    min_prom = prominence_ratio * data_range

    extrema: list = []
    for i in range(1, n - 1):
        if values[i] < values[i - 1] and values[i] <= values[i + 1]:
            extrema.append([i, "min", float(values[i])])
        elif values[i] > values[i - 1] and values[i] >= values[i + 1]:
            extrema.append([i, "max", float(values[i])])

    changed = True
    while changed and len(extrema) > 1:
        changed = False
        for j in range(len(extrema) - 1):
            if abs(extrema[j][2] - extrema[j + 1][2]) < min_prom:
                del extrema[j : j + 2]
                changed = True
                break
    return extrema


# --- Cycle detection (Part 1) ---

def detect_cycles(
    times: np.ndarray,
    sigma_star: np.ndarray,
    smoothing_window: int = 3,
    min_prominence_ratio: float = 0.25,
) -> list:
    """Detect min -> max -> min cycles in a sigma_star(t) series.

    Extrema come from a light moving-average smoothing of sigma_star
    (smoothing_window points) with a prominence filter (min_prominence_ratio
    fraction of the smoothed range). Reported values are read from the
    *original* series at the detected indices; consecutive cycles share
    their boundary minimum. Raises ValueError on empty input, length
    mismatch, or < 2 points. Does not mutate times or sigma_star.
    """
    times = np.asarray(times, dtype=float)
    sigma_star = np.asarray(sigma_star, dtype=float)
    if times.size == 0 or sigma_star.size == 0:
        raise ValueError("detect_cycles: times and sigma_star must be non-empty")
    if times.size != sigma_star.size:
        raise ValueError(
            f"detect_cycles: times and sigma_star must have equal length; "
            f"got {times.size}, {sigma_star.size}"
        )
    if times.size < 2:
        raise ValueError("detect_cycles: need at least 2 points")

    smoothed = _moving_average(sigma_star, smoothing_window)
    extrema = _find_local_extrema(smoothed, min_prominence_ratio)

    cycles: list = []
    i = 0
    while i + 2 < len(extrema):
        e0, e1, e2 = extrema[i], extrema[i + 1], extrema[i + 2]
        if e0[1] == "min" and e1[1] == "max" and e2[1] == "min":
            i0, i1, i2 = e0[0], e1[0], e2[0]
            sigma_min_start = float(sigma_star[i0])
            sigma_max = float(sigma_star[i1])
            sigma_min_end = float(sigma_star[i2])
            cycles.append({
                "t_min_start": float(times[i0]),
                "t_max": float(times[i1]),
                "t_min_end": float(times[i2]),
                "period": float(times[i2] - times[i0]),
                "sigma_min_start": sigma_min_start,
                "sigma_max": sigma_max,
                "sigma_min_end": sigma_min_end,
                "amplitude": float(sigma_max - 0.5 * (sigma_min_start + sigma_min_end)),
            })
            i += 2  # shared boundary minimum starts the next cycle
        else:
            i += 1
    return cycles


def cycle_ratios(cycles: list) -> dict:
    """Consecutive-cycle period/amplitude ratios and the sigma_min trend.

    Returns period_ratios/amplitude_ratios (length len(cycles)-1),
    sigma_min_values (one per cycle boundary), geom_mean_period_ratio and
    sigma_min_trend_slope (OLS slope of sigma_min_values vs. boundary
    index; both nan if fewer than 2 cycles), and sigma_min_decreasing
    (True iff the slope is finite and < 0). Raises ValueError if cycles is
    not a list. Does not mutate cycles.
    """
    if not isinstance(cycles, list):
        raise ValueError("cycle_ratios: cycles must be a list")

    n = len(cycles)
    sigma_min_values = [float(c["sigma_min_start"]) for c in cycles]
    if cycles:
        sigma_min_values.append(float(cycles[-1]["sigma_min_end"]))

    if n < 2:
        return {
            "n_cycles": n, "period_ratios": [], "amplitude_ratios": [],
            "sigma_min_values": sigma_min_values,
            "geom_mean_period_ratio": float("nan"),
            "sigma_min_trend_slope": float("nan"), "sigma_min_decreasing": False,
        }

    periods = np.array([c["period"] for c in cycles], dtype=float)
    amplitudes = np.array([c["amplitude"] for c in cycles], dtype=float)
    period_ratios = (periods[1:] / periods[:-1]).tolist()
    amplitude_ratios = [
        float(amplitudes[k + 1] / amplitudes[k]) if amplitudes[k] != 0.0 else float("nan")
        for k in range(len(amplitudes) - 1)
    ]

    pr = np.asarray(period_ratios, dtype=float)
    valid_pr = pr[np.isfinite(pr) & (pr > 0.0)]
    geom_mean_period_ratio = (
        float(np.exp(np.mean(np.log(valid_pr)))) if valid_pr.size else float("nan")
    )

    x = np.arange(len(sigma_min_values), dtype=float)
    y = np.asarray(sigma_min_values, dtype=float)
    try:
        slope, _ = linregress_slope(x, y)
    except ValueError:
        slope = float("nan")

    return {
        "n_cycles": n,
        "period_ratios": [float(r) for r in period_ratios],
        "amplitude_ratios": amplitude_ratios,
        "sigma_min_values": sigma_min_values,
        "geom_mean_period_ratio": geom_mean_period_ratio,
        "sigma_min_trend_slope": float(slope),
        "sigma_min_decreasing": bool(np.isfinite(slope) and slope < 0.0),
    }


def cascade_verdict(cycles: list, ratios: dict) -> str:
    """Classify cycles: CASCADE-LEANING / PERSISTENT / INSUFFICIENT-CYCLES /
    AMBIGUOUS. INSUFFICIENT-CYCLES (< 3 complete cycles) is the EXPECTED,
    honest answer on short runs — never stretch 1-2 cycles into a trend.
    CASCADE-LEANING needs geom_mean_period_ratio < 0.8 AND decreasing
    sigma_min minima (both required by Tao's cascade construction).
    PERSISTENT needs ratio >= 0.8 AND non-decreasing minima. Mixed signal
    => AMBIGUOUS. None of these outcomes is inherently good or bad for
    the program (see module docstring).
    """
    if len(cycles) < _MIN_COMPLETE_CYCLES:
        return "INSUFFICIENT-CYCLES"

    gm = ratios.get("geom_mean_period_ratio", float("nan"))
    decreasing = bool(ratios.get("sigma_min_decreasing", False))
    if not np.isfinite(gm):
        return "AMBIGUOUS"

    if gm < _CASCADE_PERIOD_RATIO_THRESHOLD and decreasing:
        return "CASCADE-LEANING"
    if gm >= _CASCADE_PERIOD_RATIO_THRESHOLD and not decreasing:
        return "PERSISTENT"
    return "AMBIGUOUS"


def _phase_report(
    times: np.ndarray, sigma_vals: np.ndarray,
    smoothing_window: int, min_prominence_ratio: float,
) -> dict:
    """One phase's (full-trajectory or growth-phase-only) cycle report."""
    cycles = (
        detect_cycles(times, sigma_vals, smoothing_window, min_prominence_ratio)
        if times.size >= 3 else []
    )
    ratios = cycle_ratios(cycles)
    verdict = cascade_verdict(cycles, ratios)
    return {"cycles": cycles, "ratios": ratios, "verdict": verdict, "n_points": int(times.size)}


def _analyze_arrays(
    times, M_vals, sigma_vals,
    smoothing_window: int = 3, min_prominence_ratio: float = 0.25,
) -> dict:
    """Shared core: cycle-analyze a (times, M, sigma_star) timeseries, over
    the full trajectory and restricted to the growth phase (steps up to
    and including the time of max M) — see module docstring's "late-time-
    growth caveat". Does not mutate inputs.
    """
    times = np.asarray(times, dtype=float)
    M_vals = np.asarray(M_vals, dtype=float)
    sigma_vals = np.asarray(sigma_vals, dtype=float)
    if times.size == 0:
        raise ValueError("_analyze_arrays: empty timeseries")
    if not (times.size == M_vals.size == sigma_vals.size):
        raise ValueError("_analyze_arrays: times, M_vals, sigma_vals must have equal length")

    full = _phase_report(times, sigma_vals, smoothing_window, min_prominence_ratio)

    peak_idx = int(np.argmax(M_vals))
    growth_times = times[: peak_idx + 1]
    growth_sigma = sigma_vals[: peak_idx + 1]
    growth = _phase_report(growth_times, growth_sigma, smoothing_window, min_prominence_ratio)

    return {
        "full_trajectory": full, "growth_phase": growth,
        "m_peak_index": peak_idx, "m_peak_time": float(times[peak_idx]),
        "n_points": int(times.size),
    }


def analyze_logged_experiment(
    exp_id: str,
    db_path: str = _DEFAULT_DB,
    table: str = "layer4_bridge_timeseries",
    smoothing_window: int = 3,
    min_prominence_ratio: float = 0.25,
) -> dict:
    """Pull (step, t, M, sigma_star) for exp_id from `table` and cycle-analyze
    it. `table` must be a table this program logs M/sigma_star timeseries
    to (layer4_bridge_timeseries for the six existing alignment_bridge
    runs, layer4_sigma_cycles_timeseries for the new long-cycle runs in
    sigma_cycles_run.py). Raises ValueError for an unknown table or if no
    rows are found for exp_id.
    """
    if table not in _ALLOWED_TIMESERIES_TABLES:
        raise ValueError(f"analyze_logged_experiment: unknown table {table!r}")

    conn = sqlite3.connect(db_path)
    try:
        cur = conn.execute(
            f"SELECT step, t, M, sigma_star FROM {table} WHERE exp_id = ? ORDER BY step",
            (exp_id,),
        )
        rows = cur.fetchall()
    finally:
        conn.close()

    if not rows:
        raise ValueError(
            f"analyze_logged_experiment: no rows found for exp_id={exp_id!r} in {table!r}"
        )

    _steps, times, M_vals, sigma_vals = zip(*rows)
    analysis = _analyze_arrays(times, M_vals, sigma_vals, smoothing_window, min_prominence_ratio)
    return {"exp_id": exp_id, "source_table": table, **analysis}


# --- Results logging (retry-on-locked backoff, layer4/depletion_diagnostic.py pattern) ---

def _ensure_cycles_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path, timeout=30.0)
    conn.execute(
        "CREATE TABLE IF NOT EXISTS layer4_sigma_cycles ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, exp_id TEXT NOT NULL, "
        "source_exp_id TEXT NOT NULL, phase TEXT NOT NULL, timestamp TEXT NOT NULL, "
        "n_cycles INTEGER, geom_mean_period_ratio REAL, sigma_min_trend_slope REAL, "
        "sigma_min_first REAL, sigma_min_last REAL, verdict TEXT NOT NULL, "
        "key_metric TEXT, notes TEXT)"
    )
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_sigma_cycles_timeseries (
            id       INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id   TEXT NOT NULL,
            step     INTEGER,
            t        REAL,
            M        REAL,
            sigma_star REAL
        )
    """)
    conn.commit()
    return conn


def build_cycle_log_row(
    exp_id: str, source_exp_id: str, phase: str, phase_report: dict, notes: str = "",
) -> dict:
    """Build one layer4_sigma_cycles row dict from a _phase_report() result."""
    ratios = phase_report["ratios"]
    sigma_vals = ratios.get("sigma_min_values", [])
    gm = ratios.get("geom_mean_period_ratio", float("nan"))
    slope = ratios.get("sigma_min_trend_slope", float("nan"))
    n_cycles = len(phase_report["cycles"])
    return {
        "exp_id": exp_id, "source_exp_id": source_exp_id, "phase": phase,
        "timestamp": datetime.utcnow().isoformat(), "n_cycles": n_cycles,
        "geom_mean_period_ratio": gm, "sigma_min_trend_slope": slope,
        "sigma_min_first": float(sigma_vals[0]) if sigma_vals else float("nan"),
        "sigma_min_last": float(sigma_vals[-1]) if sigma_vals else float("nan"),
        "verdict": phase_report["verdict"],
        "key_metric": (
            f"n_cycles={n_cycles} geom_mean_period_ratio={gm:.4f} "
            f"sigma_min_trend_slope={slope:.4e} verdict={phase_report['verdict']}"
        ),
        "notes": notes,
    }


def _write_cycle_row(db_path: str, row: dict, timeseries: list | None) -> None:
    conn = _ensure_cycles_db(db_path)
    try:
        cols = [
            "exp_id", "source_exp_id", "phase", "timestamp", "n_cycles",
            "geom_mean_period_ratio", "sigma_min_trend_slope",
            "sigma_min_first", "sigma_min_last", "verdict", "key_metric", "notes",
        ]
        values = tuple(row.get(c) for c in cols)
        ph = ", ".join("?" for _ in cols)
        conn.execute(
            f"INSERT INTO layer4_sigma_cycles ({', '.join(cols)}) VALUES ({ph})", values,
        )
        if timeseries:
            for rec in timeseries:
                conn.execute(
                    """
                    INSERT INTO layer4_sigma_cycles_timeseries
                    (exp_id, step, t, M, sigma_star) VALUES (?, ?, ?, ?, ?)
                    """,
                    (row["exp_id"], rec["step"], rec["t"], rec["M"], rec["sigma_star"]),
                )
        conn.commit()
    finally:
        conn.close()


def log_cycle_analysis(
    db_path: str,
    row: dict,
    timeseries: list | None = None,
    max_retries: int = _DB_MAX_RETRIES,
    backoff_base_s: float = _DB_BACKOFF_BASE_S,
) -> None:
    """Append one layer4_sigma_cycles row (+ optional timeseries) to results.db.
    Retries with short exponential backoff on sqlite3.OperationalError
    ("database is locked"); re-raises on any other error or once retries
    are exhausted. `row` must supply the columns used by build_cycle_log_row.
    """
    last_err: Exception | None = None
    for attempt in range(max_retries):
        try:
            _write_cycle_row(db_path, row, timeseries)
            return
        except sqlite3.OperationalError as exc:
            if "locked" not in str(exc).lower() or attempt == max_retries - 1:
                raise
            last_err = exc
            time.sleep(backoff_base_s * (2 ** attempt))
    if last_err is not None:  # pragma: no cover — defensive, unreachable
        raise last_err


def log_existing_experiment_analysis(
    db_path: str, exp_id: str, report: dict, notes: str = "",
) -> None:
    """Log both phases of an analyze_logged_experiment() report. Does not
    duplicate the source timeseries (already in layer4_bridge_timeseries).
    """
    full_row = build_cycle_log_row(exp_id, exp_id, "full", report["full_trajectory"], notes)
    growth_row = build_cycle_log_row(exp_id, exp_id, "growth", report["growth_phase"], notes)
    log_cycle_analysis(db_path, full_row, timeseries=None)
    log_cycle_analysis(db_path, growth_row, timeseries=None)


# --- CLI entry point ---

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Sigma* cycle analysis of existing logged experiments")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db", default=_DEFAULT_DB)
    args = parser.parse_args()

    for exp_id in _EXISTING_EXPERIMENTS:
        report = analyze_logged_experiment(exp_id, db_path=args.db)
        log_existing_experiment_analysis(args.db, exp_id, report)
        if args.verbose:
            ft, gp = report["full_trajectory"], report["growth_phase"]
            print(f"{exp_id}: full n_cycles={len(ft['cycles'])} verdict={ft['verdict']}  |  "
                  f"growth n_cycles={len(gp['cycles'])} verdict={gp['verdict']}")


if __name__ == "__main__":
    main()
