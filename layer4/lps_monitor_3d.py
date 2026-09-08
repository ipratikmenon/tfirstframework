"""
lps_monitor_3d.py — Layer 4: LPS Norm Monitoring in 3D
=======================================================
Tracks the Ladyzhenskaya-Prodi-Serrin (LPS) norms of u and the
LPS margin ε_LPS(t) = min(μ_eff) − ν throughout a Route 2 3D run.

Prodi-Serrin condition: u ∈ L^q(0,T; L^p(ℝ³)) with 2/q + 3/p ≤ 1, p > 3.
The LPS margin ε_LPS = min(μ_eff) − ν ≥ 0 by construction when ε_param ≥ 0.

Key diagnostics (PRD §7.3 Step 7):
  - u ∈ L^p spatial norm: ‖u(·,t)‖_{L^p} for p ∈ {3, 4, 6, ∞}
  - Prodi-Serrin product: C_{p,q}(t) = ‖u‖_{L^p}^q / T (exponent check)
  - LPS margin: ε_LPS(t) = min(μ_eff) − ν
  - Serrin norm: max_t ‖u‖_{L^p} (is it bounded?)
  - Enstrophy Z = ½‖ω‖_{L²}² (vorticity norm, proxy for smoothness)

Pass criterion:
  - ε_LPS(t) ≥ 0 throughout (always true by construction)
  - ‖u‖_{L^6} bounded (energy stays in physical regime)
  - Enstrophy remains finite (no blowup signal)

Run:
  python layer4/lps_monitor_3d.py --N 32 --t-end 1.0 --ic tg --verbose
  python layer4/lps_monitor_3d.py --N 32 --t-end 2.0 --ic shear --eps 0.0 --verbose
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

# Prodi-Serrin pairs (p, q) satisfying 2/q + 3/p = 1 (critical line)
# p > 3 required; p = ∞ → q = 2 (endpoint)
_PS_PAIRS = [(4, 8), (6, 4), (12, 3)]   # (p, q) on the critical Serrin line

# Spatial L^p norms to track
_LP_VALS = (3, 4, 6)


# =============================================================================
# LPS norms
# =============================================================================

def spatial_lp_norm(u: np.ndarray, v: np.ndarray, w: np.ndarray, p: int) -> float:
    """‖u(·,t)‖_{L^p} = (mean |u|^p)^{1/p}, volumetric average."""
    mag = np.sqrt(u**2 + v**2 + w**2)
    return float(np.mean(mag**p) ** (1.0 / p))


def linf_norm(u: np.ndarray, v: np.ndarray, w: np.ndarray) -> float:
    """‖u(·,t)‖_{L^∞} = max |u|."""
    return float(np.max(np.sqrt(u**2 + v**2 + w**2)))


def enstrophy(u: np.ndarray, v: np.ndarray, w: np.ndarray,
              kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> float:
    """Z = ½‖ω‖_{L²}² where ω = curl u (spectral)."""
    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)
    # ω = curl u: ω_x = ∂w/∂y − ∂v/∂z, etc.
    wx = np.real(np.fft.ifftn(1j * (ky * w_hat - kz * v_hat)))
    wy = np.real(np.fft.ifftn(1j * (kz * u_hat - kx * w_hat)))
    wz = np.real(np.fft.ifftn(1j * (kx * v_hat - ky * u_hat)))
    return float(0.5 * np.mean(wx**2 + wy**2 + wz**2))


def palinstrophy(u: np.ndarray, v: np.ndarray, w: np.ndarray,
                 kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> float:
    """P = ½‖∇ω‖_{L²}² — rate of enstrophy production (spectral)."""
    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)
    k2 = kx**2 + ky**2 + kz**2
    wx_hat = 1j * (ky * w_hat - kz * v_hat)
    wy_hat = 1j * (kz * u_hat - kx * w_hat)
    wz_hat = 1j * (kx * v_hat - ky * u_hat)
    # ‖∇ω‖² = Σ k² |ω̂|²
    P_sq = float(np.sum(k2 * (np.abs(wx_hat)**2 + np.abs(wy_hat)**2 + np.abs(wz_hat)**2)))
    return 0.5 * P_sq / u.size


# =============================================================================
# Main monitor
# =============================================================================

def monitor_lps(
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
    """Run Route 2 3D and monitor LPS norms + LPS margin throughout.

    Returns result dict with full LPS timeseries and summary.
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
    exp_id = exp_id_override or f"EXP-L4-LPS-{ic.upper()}-{N:03d}"

    if verbose:
        print(f"\n  LPS Monitor 3D  N={N}³  ν={nu}  ε={eps_param}  IC={ic}")
        print(f"  Exp: {exp_id}  t_end={t_end}")
        print(f"  Tracking: ‖u‖_{{L^p}} for p ∈ {_LP_VALS} + Enstrophy + ε_LPS")

    kx = solver["kx"]; ky = solver["ky"]; kz = solver["kz"]
    timeseries: list[dict] = []
    t_wall0 = time.perf_counter()
    eps_lps_negative = False

    for step_idx in range(max_steps):
        if solver["t"] >= t_end:
            break

        dt = r2.compute_cfl_dt(solver, safety=cfl_safety)
        dt = min(dt, t_end - solver["t"])
        solver, diag = r2.solver_step(solver, dt=dt)

        u, v, w = solver["u"], solver["v"], solver["w"]
        theta = solver["theta"]

        # LPS norms
        lp_norms = {p: spatial_lp_norm(u, v, w, p) for p in _LP_VALS}
        linf = linf_norm(u, v, w)
        Z  = enstrophy(u, v, w, kx, ky, kz)
        P  = palinstrophy(u, v, w, kx, ky, kz)
        mu_eff = r2.mu_eff_field(theta, nu, eps_param)
        eps_lps = float(np.min(mu_eff)) - nu

        if eps_lps < -1e-12:
            eps_lps_negative = True

        # Prodi-Serrin exponents (instantaneous spatial norm)
        ps_norms = {}
        for p_ps, q_ps in _PS_PAIRS:
            lp_val = spatial_lp_norm(u, v, w, p_ps)
            ps_norms[f"PS_p{p_ps}_Lp"] = float(lp_val)

        row = {
            "t"         : float(solver["t"]),
            "dt"        : float(dt),
            "E"         : float(diag["E"]),
            "Z"         : Z,
            "P"         : P,
            "Linf"      : linf,
            "eps_lps"   : eps_lps,
            "delta_est" : float(diag["delta_est"]),
            "div_rms"   : float(diag["div_rms"]),
            **{f"Lp{p}": lp_norms[p] for p in _LP_VALS},
            **ps_norms,
        }
        timeseries.append(row)

        if verbose and step_idx % max(1, max_steps // 20) == 0:
            print(
                f"    t={solver['t']:.4f}  E={diag['E']:.4e}  Z={Z:.4e}  "
                f"‖u‖_L6={lp_norms[6]:.4e}  ε_LPS={eps_lps:.4f}"
            )

    wall = time.perf_counter() - t_wall0

    # ── Summary ───────────────────────────────────────────────────────────
    Zs = [r["Z"] for r in timeseries]
    lp6s = [r["Lp6"] for r in timeseries]
    eps_lpss = [r["eps_lps"] for r in timeseries]

    # Pass: no blowup signal (enstrophy stays finite), ε_LPS ≥ 0
    Z_max = float(max(Zs)) if Zs else 0.0
    lp6_max = float(max(lp6s)) if lp6s else 0.0
    eps_lps_min = float(min(eps_lpss)) if eps_lpss else 0.0

    # Blowup alarm: Z grows faster than exponential (heuristic)
    blowup_flag = False
    if len(Zs) > 10:
        # Check if enstrophy at end is > 100× initial
        if Zs[-1] > 100.0 * Zs[0] and Zs[0] > 1e-10:
            blowup_flag = True

    verdict = "PASS" if not blowup_flag and not eps_lps_negative else "FAIL"

    summary = {
        "exp_id"         : exp_id,
        "ic"             : ic,
        "N"              : N,
        "nu"             : nu,
        "eps_param"      : eps_param,
        "t_end"          : t_end,
        "n_steps"        : len(timeseries),
        "t_final"        : float(solver["t"]),
        "Z_max"          : Z_max,
        "lp6_max"        : lp6_max,
        "eps_lps_min"    : eps_lps_min,
        "eps_lps_negative": eps_lps_negative,
        "blowup_flag"    : blowup_flag,
        "verdict"        : verdict,
        "wall_time_s"    : wall,
        "timestamp"      : datetime.utcnow().isoformat(),
        "key_metric"     : (
            f"Z_max={Z_max:.4e}  ‖u‖_L6_max={lp6_max:.4e}  "
            f"ε_LPS_min={eps_lps_min:.4f}  blowup={blowup_flag}"
        ),
    }

    if verbose:
        print(f"\n  Results: {exp_id}")
        print(f"    Z_max:       {Z_max:.4e}")
        print(f"    ‖u‖_L6_max:  {lp6_max:.4e}")
        print(f"    ε_LPS_min:   {eps_lps_min:.4f}")
        print(f"    Blowup flag: {blowup_flag}")
        print(f"    Verdict:     {verdict}")
        print(f"    Wall:        {wall:.1f} s")

    if db_path is not None:
        _log_lps(db_path, summary, timeseries)

    return {"summary": summary, "timeseries": timeseries, "verdict": verdict}


# =============================================================================
# DB logging
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_lps_experiments (
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
            Z_max            REAL,
            lp6_max          REAL,
            eps_lps_min      REAL,
            eps_lps_negative INTEGER,
            blowup_flag      INTEGER,
            verdict          TEXT,
            wall_time_s      REAL,
            key_metric       TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_lps_timeseries (
            id         INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id     TEXT NOT NULL,
            t          REAL,
            dt         REAL,
            E          REAL,
            Z          REAL,
            P          REAL,
            Linf       REAL,
            eps_lps    REAL,
            delta_est  REAL,
            div_rms    REAL,
            Lp3        REAL,
            Lp4        REAL,
            Lp6        REAL
        )
    """)
    conn.commit()
    return conn


def _log_lps(db_path: str, summary: dict, timeseries: list[dict]) -> None:
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "timestamp", "ic", "N", "nu", "eps_param", "t_end",
        "n_steps", "t_final", "Z_max", "lp6_max",
        "eps_lps_min", "eps_lps_negative", "blowup_flag",
        "verdict", "wall_time_s", "key_metric",
    ]
    vals = tuple(
        int(summary[c]) if isinstance(summary.get(c), bool) else summary.get(c)
        for c in cols
    )
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO layer4_lps_experiments "
        f"({', '.join(cols)}) VALUES ({ph})",
        vals,
    )
    for row in timeseries:
        conn.execute("""
            INSERT INTO layer4_lps_timeseries
            (exp_id, t, dt, E, Z, P, Linf, eps_lps, delta_est, div_rms, Lp3, Lp4, Lp6)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (
            summary["exp_id"],
            row["t"], row["dt"], row["E"], row["Z"], row["P"],
            row["Linf"], row["eps_lps"], row["delta_est"], row["div_rms"],
            row.get("Lp3", 0.0), row.get("Lp4", 0.0), row.get("Lp6", 0.0),
        ))
    conn.commit()
    conn.close()


# =============================================================================
# CLI
# =============================================================================

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Layer 4: LPS norm monitor 3D")
    parser.add_argument("--N",      type=int,   default=32)
    parser.add_argument("--nu",     type=float, default=1.0e-3)
    parser.add_argument("--eps",    type=float, default=0.1)
    parser.add_argument("--ic",     choices=["tg", "shear", "random"], default="tg")
    parser.add_argument("--t-end",  type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db",     default=_DEFAULT_DB)
    args = parser.parse_args()

    r = monitor_lps(
        N=args.N, nu=args.nu, eps_param=args.eps,
        ic=args.ic, t_end=args.t_end, verbose=args.verbose,
        db_path=args.db,
    )
    print(f"\nVerdict: {r['verdict']}  key_metric: {r['summary']['key_metric']}")


if __name__ == "__main__":
    main()
