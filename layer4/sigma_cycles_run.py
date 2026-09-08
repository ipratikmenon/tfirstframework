"""
layer4/sigma_cycles_run.py — Long, Densely-Sampled Cycle-Detection Runs (Part 2)
=================================================================================
Companion to layer4/sigma_cycles.py (read that module's docstring first for
the full scientific motivation and falsification framing). The existing
logged sigma_star runs are ~12-20 time units / 41 snapshots — too short to
contain 3+ complete align/disorder cycles, so cycle_ratios/cascade_verdict
there mostly return INSUFFICIENT-CYCLES (the expected, honest answer, not
a failure). This module runs longer, higher-Reynolds, MUCH more densely
sampled (record_every~2 instead of ~10) Route 2 3D runs specifically to
give the cycle detector a real chance at 3+ cycles.

Unforced-flow limitation (read before drawing conclusions from these runs):
this reuses layer3/route2_3D.py's UNFORCED solver exactly as
layer4/alignment_bridge.py does — no external forcing is added, staying
inside the Prize's unforced setting (div u = 0). Unforced flows decay, so
any run here can only ever exhibit a handful of cycles before kinetic
energy dissipates away. Testing whether cycles recur *forever* at
shrinking period (the full cascade claim) would require sustained forcing
to hold the flow near a statistically steady state — that is a mechanism
study, not a Prize-relevant computation, and is explicitly NOT attempted
here or anywhere in this task.

Reused, not reimplemented: layer3/route2_3D.py (make_solver, taylor_green_ic,
shear_layer_ic, set_ic, compute_cfl_dt, solver_step), layer4/alignment_bridge.py
(bridge_diagnostics, _vorticity_from_velocity, adv_vortex_tubes_ic — sigma_star
itself is only ever computed by bridge_diagnostics), layer4/kcorr_reynolds_sweep.py
(spectral_tail_fraction), and this package's own layer4/sigma_cycles.py
(_analyze_arrays, log_cycle_analysis, build_cycle_log_row). New here: the
run loop's wall-clock budget (so a single foreground call stays well under
the ~9-minute limit regardless of n_steps) and its own results logging.

House rules: pure/immutable functions; every logged result carries
claim_id + key_metric + verdict; random seed always set and logged;
retry-with-backoff on a locked results.db.

Run:
  python layer4/sigma_cycles_run.py --ic tg --N 64 --nu 2.5e-4 --n-steps 1600 --verbose
  python layer4/sigma_cycles_run.py --ic shear --N 64 --nu 2.5e-4 --n-steps 1600 --verbose
  python layer4/sigma_cycles_run.py --ic tg --N 16 --n-steps 20 --verbose   # smoke test
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from datetime import datetime

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

import route2_3D as r2  # noqa: E402 — Route 2 3D solver (reused, not reimplemented)
import layer4.alignment_bridge as ab  # noqa: E402 — sigma_star / vorticity helpers (reused)
from layer4.kcorr_reynolds_sweep import spectral_tail_fraction  # noqa: E402
from layer4.sigma_cycles import (  # noqa: E402
    _analyze_arrays,
    build_cycle_log_row,
    log_cycle_analysis,
)

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_SEED = 42
_DEFAULT_MAX_WALL_S = 480.0
_TAIL_FRACTION_WARN = 1.0e-3
_IC_TAGS = {"tg": "TG", "shear": "SH", "adv": "ADV"}
_IC_CLAIM_SUFFIX = {"tg": "tg", "shear": "shear", "adv": "adv"}


def run_long_cycle_experiment(
    ic: str,
    N: int,
    n_steps: int,
    nu: float,
    seed: int,
    record_every: int = 2,
    cfl_safety: float = 0.4,
    max_wall_s: float = _DEFAULT_MAX_WALL_S,
    verbose: bool = False,
    db_path: str | None = None,
    exp_id_override: str | None = None,
    claim_id_override: str | None = None,
) -> dict:
    """Run route2_3D (via alignment_bridge's helpers) tracking M/sigma_star
    at high temporal resolution, to try to catch multiple align/disorder
    cycles within a single, necessarily finite, decaying, unforced run.

    Reuses layer3/route2_3D.py's solver/IC/step functions and
    layer4/alignment_bridge.py's bridge_diagnostics/vorticity helper
    exactly as run_bridge_experiment does — no solver or sigma_star logic
    is reimplemented here. New here: dense recording (record_every~2 vs.
    alignment_bridge's typical 10), a wall-clock budget (max_wall_s) so a
    single call stays well under the foreground-call time limit even for
    a large n_steps, and the resolution-sanity tail-fraction check.

    If the wall-clock budget is hit before n_steps completes, the run
    stops early and reports the actual steps_completed and
    stopped_early_wall_budget=True — never silently truncated without
    saying so. See module docstring for the unforced-flow limitation:
    this can only ever show a handful of cycles before the flow decays.
    """
    if ic not in _IC_TAGS:
        raise ValueError(f"Unknown ic='{ic}'; use 'tg', 'shear', or 'adv'.")
    if n_steps <= 0:
        raise ValueError(f"n_steps must be > 0, got {n_steps}")
    if record_every <= 0:
        raise ValueError(f"record_every must be > 0, got {record_every}")

    solver = r2.make_solver(N=N, nu=nu, eps_param=0.1)
    if ic == "tg":
        ic_dict = r2.taylor_green_ic(N)
    elif ic == "shear":
        ic_dict = r2.shear_layer_ic(N)
    else:
        ic_dict = ab.adv_vortex_tubes_ic(N, seed=seed)
    solver = r2.set_ic(solver, ic_dict)

    dx = 2.0 * np.pi / N
    kx, ky, kz = solver["kx"], solver["ky"], solver["kz"]

    timeseries: list = []
    t0 = time.perf_counter()
    steps_completed = 0
    stopped_early_wall_budget = False

    for step_idx in range(n_steps):
        dt = r2.compute_cfl_dt(solver, safety=cfl_safety)
        solver, _diag = r2.solver_step(solver, dt)
        steps_completed = step_idx + 1

        is_last = step_idx == n_steps - 1
        if step_idx % record_every == 0 or is_last:
            wx, wy, wz = ab._vorticity_from_velocity(
                solver["u"], solver["v"], solver["w"], kx, ky, kz
            )
            try:
                bd = ab.bridge_diagnostics(wx, wy, wz, dx)
                M_val, sigma_val = bd["M"], bd["sigma_star"]
            except ValueError:
                M_val, sigma_val = 0.0, 0.0
            timeseries.append({
                "step": step_idx, "t": solver["t"], "M": M_val, "sigma_star": sigma_val,
            })
            if verbose:
                print(f"  step {step_idx:5d}  t={solver['t']:.4f}  M={M_val:.4e}  "
                      f"sigma*={sigma_val:.4e}")

        if (time.perf_counter() - t0) > max_wall_s:
            stopped_early_wall_budget = True
            break

    wall = time.perf_counter() - t0
    tail_fraction = spectral_tail_fraction(solver["u"], solver["v"], solver["w"], N)
    resolution_warning = bool(tail_fraction > _TAIL_FRACTION_WARN)

    times = [rec["t"] for rec in timeseries]
    M_vals = [rec["M"] for rec in timeseries]
    sigma_vals = [rec["sigma_star"] for rec in timeseries]
    analysis = _analyze_arrays(times, M_vals, sigma_vals)

    tag = _IC_TAGS[ic]
    exp_id = exp_id_override or f"EXP-L4-CYCLE-{tag}-001"
    claim_id = claim_id_override or f"sigma-cycles-{_IC_CLAIM_SUFFIX[ic]}"

    result = {
        "exp_id": exp_id, "claim_id": claim_id,
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic, "N": N, "n_steps_requested": n_steps,
        "steps_completed": steps_completed,
        "stopped_early_wall_budget": stopped_early_wall_budget,
        "seed": seed, "nu": nu, "record_every": record_every,
        "t_final": solver["t"], "wall_time_s": wall,
        "tail_fraction": tail_fraction, "resolution_warning": resolution_warning,
        "timeseries": timeseries, "analysis": analysis,
    }

    if verbose:
        gp, ft = analysis["growth_phase"], analysis["full_trajectory"]
        print(f"\n  Sigma-cycle experiment ({ic}, N={N}, steps={steps_completed}/{n_steps})")
        print(f"    t_final={solver['t']:.4f}  wall={wall:.1f}s  "
              f"tail_fraction={tail_fraction:.3e}  resolution_warning={resolution_warning}")
        print(f"    full-trajectory:  n_cycles={len(ft['cycles'])}  verdict={ft['verdict']}")
        print(f"    growth-phase:     n_cycles={len(gp['cycles'])}  verdict={gp['verdict']}")

    if db_path is not None:
        _log_long_cycle_result(db_path, result)

    return result


def _log_long_cycle_result(db_path: str, result: dict) -> None:
    """Log both phases (full + growth) of a run_long_cycle_experiment result,
    plus the raw timeseries (attached to the "full" row only)."""
    exp_id = result["exp_id"]
    analysis = result["analysis"]
    notes = (
        f"tail_fraction={result['tail_fraction']:.4e} "
        f"resolution_warning={result['resolution_warning']} "
        f"steps_completed={result['steps_completed']}/{result['n_steps_requested']} "
        f"stopped_early_wall_budget={result['stopped_early_wall_budget']}"
    )
    full_row = build_cycle_log_row(exp_id, exp_id, "full", analysis["full_trajectory"], notes)
    growth_row = build_cycle_log_row(exp_id, exp_id, "growth", analysis["growth_phase"], notes)
    log_cycle_analysis(db_path, full_row, timeseries=result["timeseries"])
    log_cycle_analysis(db_path, growth_row, timeseries=None)


def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Long, dense-sampled sigma_star cycle run")
    parser.add_argument("--ic", choices=["tg", "shear", "adv"], default="tg")
    parser.add_argument("--N", type=int, default=64)
    parser.add_argument("--n-steps", type=int, default=1600)
    parser.add_argument("--nu", type=float, default=2.5e-4)
    parser.add_argument("--record-every", type=int, default=2)
    parser.add_argument("--seed", type=int, default=_DEFAULT_SEED)
    parser.add_argument("--max-wall-s", type=float, default=_DEFAULT_MAX_WALL_S)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db", default=_DEFAULT_DB)
    args = parser.parse_args()

    r = run_long_cycle_experiment(
        ic=args.ic, N=args.N, n_steps=args.n_steps, nu=args.nu, seed=args.seed,
        record_every=args.record_every, max_wall_s=args.max_wall_s,
        verbose=args.verbose, db_path=args.db,
    )
    print(f"\nsteps_completed={r['steps_completed']}/{r['n_steps_requested']}  "
          f"resolution_warning={r['resolution_warning']}")


if __name__ == "__main__":
    main()
