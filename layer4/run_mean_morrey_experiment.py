"""
layer4/run_mean_morrey_experiment.py

Run mean Morrey diagnostic on Taylor-Green (TG) and shear velocity fields,
then log results to results/results.db.

Experiments:
  EXP-L4-MM-TG-001  — Taylor-Green vortex IC at N=32³, deterministic centres
  EXP-L4-MM-SH-001  — Shear IC u=(sin(y),0,0) at N=32³, deterministic centres
  EXP-L4-MM-TG-002  — Taylor-Green vortex IC at N=32³, random centres (n=64, seed=42)
  EXP-L4-MM-SH-002  — Shear IC u=(sin(y),0,0) at N=32³, random centres (n=64, seed=42)
"""
from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import numpy as np

# Ensure project root is importable
_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

from layer4.mean_morrey_diagnostic import (  # noqa: E402
    compute_mean_morrey_field,
    compute_mean_morrey_field_random,
    init_db,
    log_result,
    run_mean_morrey_analysis,
)

# ─────────────────────────────────────────────────────────────────────────────
# Grid / domain
# ─────────────────────────────────────────────────────────────────────────────

N = 32
L = 2 * np.pi
x = np.linspace(0, L, N, endpoint=False)
xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")

# ─────────────────────────────────────────────────────────────────────────────
# Taylor-Green IC
# ─────────────────────────────────────────────────────────────────────────────

u_tg = np.zeros((3, N, N, N))
u_tg[0] = np.sin(xx) * np.cos(yy) * np.cos(zz)
u_tg[1] = -np.cos(xx) * np.sin(yy) * np.cos(zz)
u_tg[2] = 0.0  # divergence-free in 2D extension

tg_snaps = [u_tg, u_tg * 0.85, u_tg * 0.72]
t_tg = [0.0, 0.05, 0.10]

# ─────────────────────────────────────────────────────────────────────────────
# Shear IC
# ─────────────────────────────────────────────────────────────────────────────

u_shear = np.zeros((3, N, N, N))
u_shear[0] = np.sin(yy)  # u = (sin(y), 0, 0)

shear_snaps = [u_shear, u_shear * 0.90, u_shear * 0.81]
t_shear = [0.0, 0.05, 0.10]

# ─────────────────────────────────────────────────────────────────────────────
# r and alpha values
# ─────────────────────────────────────────────────────────────────────────────

r_values = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2]
alpha_values = [0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0]

# ─────────────────────────────────────────────────────────────────────────────
# Run diagnostics
# ─────────────────────────────────────────────────────────────────────────────

print("Running mean Morrey diagnostic — EXP-L4-MM-TG-001 (Taylor-Green) …")
result_tg = run_mean_morrey_analysis(
    exp_id="EXP-L4-MM-TG-001",
    u_snapshots=tg_snaps,
    t_values=t_tg,
    grid_N=N,
    ic_type="TG",
    r_values=r_values,
    alpha_values=alpha_values,
)

print("Running mean Morrey diagnostic — EXP-L4-MM-SH-001 (Shear) …")
result_shear = run_mean_morrey_analysis(
    exp_id="EXP-L4-MM-SH-001",
    u_snapshots=shear_snaps,
    t_values=t_shear,
    grid_N=N,
    ic_type="shear",
    r_values=r_values,
    alpha_values=alpha_values,
)

print("Running mean Morrey diagnostic — EXP-L4-MM-TG-002 (Taylor-Green, random centres) …")
result_tg2 = run_mean_morrey_analysis(
    exp_id="EXP-L4-MM-TG-002",
    u_snapshots=tg_snaps,
    t_values=t_tg,
    grid_N=N,
    ic_type="TG",
    r_values=r_values,
    alpha_values=alpha_values,
    n_centres_random=64,
    seed=42,
)

print("Running mean Morrey diagnostic — EXP-L4-MM-SH-002 (Shear, random centres) …")
result_shear2 = run_mean_morrey_analysis(
    exp_id="EXP-L4-MM-SH-002",
    u_snapshots=shear_snaps,
    t_values=t_shear,
    grid_N=N,
    ic_type="shear",
    r_values=r_values,
    alpha_values=alpha_values,
    n_centres_random=64,
    seed=42,
)

# ─────────────────────────────────────────────────────────────────────────────
# Log to database
# ─────────────────────────────────────────────────────────────────────────────

db_path = _ROOT / "results" / "results.db"
db_path.parent.mkdir(parents=True, exist_ok=True)

conn = sqlite3.connect(str(db_path))
try:
    init_db(conn)
    log_result(conn, result_tg)
    log_result(conn, result_shear)
    log_result(conn, result_tg2)
    log_result(conn, result_shear2)

    row_count = conn.execute(
        "SELECT COUNT(*) FROM mean_morrey_results"
    ).fetchone()[0]
    print(f"\nDatabase write succeeded — mean_morrey_results now has {row_count} row(s).")
finally:
    conn.close()

# ─────────────────────────────────────────────────────────────────────────────
# Summary table
# ─────────────────────────────────────────────────────────────────────────────

def _fmt(val: float, fmt: str = ".4e") -> str:
    try:
        if not (val == val):   # NaN check without math import
            return "nan"
        return format(val, fmt)
    except (TypeError, ValueError):
        return str(val)


header = (
    f"{'IC':<8} {'exp_id':<22} {'alpha_thresh':>12} {'morrey_max':>12} "
    f"{'scaling_exp':>12} {'verdict':>7}"
)
sep = "-" * len(header)

print()
print("=" * len(header))
print("MEAN MORREY DIAGNOSTIC — SUMMARY")
print("=" * len(header))
print(header)
print(sep)

for res in [result_tg, result_shear, result_tg2, result_shear2]:
    at = _fmt(res["alpha_threshold"], ".3f")
    mm = _fmt(res["morrey_max"], ".4e")
    se = _fmt(res["scaling_exponent"], ".3f")
    row = (
        f"{res['ic_type']:<8} {res['exp_id']:<22} {at:>12} {mm:>12} "
        f"{se:>12} {res['verdict']:>7}"
    )
    print(row)

print(sep)

# ─────────────────────────────────────────────────────────────────────────────
# M_alpha values at alpha=1.0 vs r  (both ICs)
# ─────────────────────────────────────────────────────────────────────────────

print()
print("M_α(r) at α=1.0 — how seminorm varies with scale")
print("-" * 50)
print(f"{'r':>8}   {'TG M_1.0':>14}   {'Shear M_1.0':>14}")
print("-" * 50)

# Re-compute morrey fields for the final snapshots to retrieve (r, alpha) values
tg_mm, _ = compute_mean_morrey_field(u_tg, r_values, alpha_values)
sh_mm, _ = compute_mean_morrey_field(u_shear, r_values, alpha_values)

for r in r_values:
    tg_val = tg_mm.get((r, 1.0), float("nan"))
    sh_val = sh_mm.get((r, 1.0), float("nan"))
    print(f"{r:>8.3f}   {_fmt(tg_val):>14}   {_fmt(sh_val):>14}")

print("-" * 50)
print()
print("Done.")
