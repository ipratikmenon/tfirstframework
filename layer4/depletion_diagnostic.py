"""
layer4/depletion_diagnostic.py — H2 Vortex-Stretching Depletion Diagnostic
================================================================
Feeds open target H2 (localized Constantin–Fefferman at the self-similar
scale): is the vortex-stretching rate at the vorticity max DEPLETED
below its naive size?

Definitions (velocity (u,v,w), vorticity ω=(wx,wy,wz), periodic N³ grid,
domain [0,2π)³, spacing dx=2π/N):

  S_ij   = (∂_i u_j + ∂_j u_i)/2     strain-rate tensor, spectral, 6
                                       unique components (S_ij = S_ji)
  ê      = ω/|ω|                     vorticity direction (reused from
                                       layer4.alignment_bridge
                                       .compute_direction_field)
  α(x,t) = ê·S·ê = Σ_ij ê_i S_ij ê_j  stretching rate along ê, pointwise
  M      = max|ω|,  x* = argmax|ω|

  alpha_at_max   = α(x*,t)
  alpha_sup_high = max{ α(x,t) : |ω(x,t)| ≥ M/2 }   (sup over the
                     high-vorticity region — a slightly more robust proxy
                     than the single point x*; always includes x* since
                     |ω(x*)|=M ≥ M/2)
  depletion_ratio = alpha_sup_high / M    naive bound gives O(1); H2
                     depletion means this stays « 1 and/or DECREASES as
                     M grows
  theta_eff = log(max(alpha_sup_high, tiny)) / log(M) − 0.5
                     effective depletion exponent, defined ONLY for M>10
                     (else NaN — a single-snapshot log-log slope is not
                     meaningful until M is well away from O(1))

H2's criterion (task spec): α ≤ C·M^{1/2+θ} with θ<1/2 (naive: θ=1/2,
α~‖∇u‖_∞~M). depletion_ratio=α/M is the direct numerical read of this;
theta_eff is a secondary log-log-slope read of the same trend.

Numerical caveat: alpha_sup_high inherits α's dependence on ê, truncated
to zero below the 1e-12·M floor (Gibbs-ringing caveat documented in
alignment_bridge). {|ω|≥M/2} sits far inside that floor in practice.

House rules: spectral (FFT) derivatives only; no hardcoded fluid
properties; all functions pure/immutable; every logged result carries
claim_id + key_metric + verdict; random seed always set/logged. Reads
ONLY from layer3/route2_3D.py and layer4/alignment_bridge.py — does not
modify either.

Run:
  python layer4/depletion_diagnostic.py --ic tg --N 64 --n-steps 400 --verbose
  python layer4/depletion_diagnostic.py --ic adv --N 64 --n-steps 400 --verbose
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from datetime import datetime

import numpy as np

# ── path setup (matches layer4/alignment_bridge.py convention) ──────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

import route2_3D as r2  # noqa: E402 — Route 2 3D solver (reused, not reimplemented)
from layer4.alignment_bridge import (  # noqa: E402
    compute_direction_field,
    _vorticity_from_velocity,
    _validate_omega_grid,
    adv_vortex_tubes_ic,
)
from layer4.geometric_disorder import spectral_wavenumbers  # noqa: E402

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_SEED = 42
_IC_TAGS = {"tg": "TG", "adv": "ADV"}
_THETA_EFF_M_FLOOR = 10.0   # theta_eff only defined for M > this (per task spec)
_LOG_TINY = 1e-300          # floor for log(alpha_sup_high) when alpha_sup_high <= 0
_DB_MAX_RETRIES = 6         # sqlite "database is locked" retry policy
_DB_BACKOFF_BASE_S = 0.15


# ── Core strain-rate / stretching-rate primitives ───────────────────────────

def compute_strain_rate_field(
    u: np.ndarray, v: np.ndarray, w: np.ndarray,
    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray,
) -> dict:
    """S_ij = (∂_i u_j + ∂_j u_i)/2, via spectral (FFT) derivatives.

    Computes the 9 velocity-gradient components and folds them into the 6
    unique symmetric strain-rate entries (S_ij = S_ji, so S_yx == S_xy
    etc. are not stored separately). Does not mutate u, v, w.

    Returns dict with keys: Sxx, Syy, Szz, Sxy, Sxz, Syz — each shape ==
    u.shape.
    """
    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)

    def _d(f_hat, k_dir):
        return np.real(np.fft.ifftn(1j * k_dir * f_hat))

    dudx, dudy, dudz = _d(u_hat, kx), _d(u_hat, ky), _d(u_hat, kz)
    dvdx, dvdy, dvdz = _d(v_hat, kx), _d(v_hat, ky), _d(v_hat, kz)
    dwdx, dwdy, dwdz = _d(w_hat, kx), _d(w_hat, ky), _d(w_hat, kz)

    return {
        "Sxx": dudx,
        "Syy": dvdy,
        "Szz": dwdz,
        "Sxy": 0.5 * (dudy + dvdx),
        "Sxz": 0.5 * (dudz + dwdx),
        "Syz": 0.5 * (dvdz + dwdy),
    }


def alpha_field_from_strain(
    ex: np.ndarray, ey: np.ndarray, ez: np.ndarray, S: dict,
) -> np.ndarray:
    """α(x) = ê_i S_ij(x) ê_j, pointwise vortex-stretching rate along ê.

    S is the dict returned by compute_strain_rate_field (6 unique
    symmetric components). Does not mutate ex, ey, ez, or S's arrays.
    """
    return (
        ex ** 2 * S["Sxx"] + ey ** 2 * S["Syy"] + ez ** 2 * S["Szz"]
        + 2.0 * ex * ey * S["Sxy"]
        + 2.0 * ex * ez * S["Sxz"]
        + 2.0 * ey * ez * S["Syz"]
    )


def depletion_diagnostics(
    u: np.ndarray, v: np.ndarray, w: np.ndarray,
    wx: np.ndarray, wy: np.ndarray, wz: np.ndarray,
    dx: float,
) -> dict:
    """Compute the full H2 depletion diagnostic dict for one snapshot.

    Parameters
    ----------
    u, v, w    : velocity components, shape (N, N, N)
    wx, wy, wz : vorticity components, shape (N, N, N), = curl(u,v,w)
    dx         : grid spacing, must equal 2π/N (validated)

    Returns dict with keys: M, x_star, alpha_at_max, alpha_sup_high,
    depletion_ratio, theta_eff, enstrophy. Raises ValueError on shape
    mismatch, bad dx, or an identically-zero vorticity field (M=0: the
    direction field and depletion_ratio are undefined).

    Does not mutate any input array.
    """
    N = _validate_omega_grid(wx, wy, wz, dx)
    if not (u.shape == v.shape == w.shape == wx.shape):
        raise ValueError(
            f"velocity and vorticity must share a shape; got u={u.shape}, "
            f"wx={wx.shape}"
        )

    mag = np.sqrt(wx ** 2 + wy ** 2 + wz ** 2)
    M = float(np.max(mag))
    if M < 1e-14:
        raise ValueError(
            "depletion_diagnostics: vorticity field is identically zero "
            "(M ≈ 0); the stretching-rate diagnostics are undefined."
        )

    kx, ky, kz = spectral_wavenumbers(N)
    ex, ey, ez, _mag2, valid = compute_direction_field(wx, wy, wz)
    S = compute_strain_rate_field(u, v, w, kx, ky, kz)
    alpha = alpha_field_from_strain(ex, ey, ez, S)

    idx = np.unravel_index(int(np.argmax(mag)), mag.shape)
    x_star = tuple(int(i) for i in idx)
    alpha_at_max = float(alpha[idx])

    mask_high = (mag >= (M / 2.0)) & valid
    if np.any(mask_high):
        alpha_sup_high = float(np.max(alpha[mask_high]))
    else:
        # x* itself always satisfies |ω|=M >= M/2, so this branch is only
        # reachable if x* were excluded by the (far looser) valid_mask
        # floor — a degenerate fallback, documented rather than silent.
        alpha_sup_high = alpha_at_max

    depletion_ratio = float(alpha_sup_high / M)

    if M > _THETA_EFF_M_FLOOR:
        theta_eff = float(
            np.log(max(alpha_sup_high, _LOG_TINY)) / np.log(M) - 0.5
        )
    else:
        theta_eff = float("nan")

    enstrophy = float(0.5 * np.mean(wx ** 2 + wy ** 2 + wz ** 2))

    return {
        "M": M,
        "x_star": x_star,
        "alpha_at_max": alpha_at_max,
        "alpha_sup_high": alpha_sup_high,
        "depletion_ratio": depletion_ratio,
        "theta_eff": theta_eff,
        "enstrophy": enstrophy,
    }


# ── Production experiment: run Route 2 3D solver, track depletion diagnostics ──

def run_depletion_experiment(
    ic: str = "tg",
    N: int = 64,
    n_steps: int = 400,
    nu: float | None = None,
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
    """Run the existing Route 2 3D solver and track H2 depletion diagnostics.

    Mirrors layer4/alignment_bridge.py's run_bridge_experiment solver setup
    exactly (same make_solver / set_ic / compute_cfl_dt / solver_step
    sequence, same IC constructors) — no solver logic is reimplemented
    here. Vorticity is derived from velocity every `record_every` steps via
    alignment_bridge's own curl helper (reused, not reimplemented).

    Parameters
    ----------
    ic            : 'tg' (Taylor–Green) or 'adv' (anti-parallel Biot–Savart
                    vortex tubes — the vortex-stretching-heavy IC)
    N, n_steps, eps_param, record_every : as in alignment_bridge's version.
    nu       : kinematic viscosity; defaults to route2_3D._NU_DEFAULT
    seed     : random seed, set/logged for provenance ('tg' draws no
               randomness, 'adv' does — see alignment_bridge)
    db_path  : if given, log via log_result() (retry-on-locked backoff)
    adv_amp  : only used when ic=='adv'; overrides adv_vortex_tubes_ic's
               default peak-vorticity amplitude

    Returns a dict with the full timeseries and a PASS/FAIL verdict.
    Verdict PASS iff every recorded depletion_ratio, M, alpha_at_max, and
    alpha_sup_high is finite AND median(depletion_ratio) < 1. The
    scientific content is the depletion_ratio TREND vs M, reported
    regardless of verdict.
    """
    if ic not in _IC_TAGS:
        raise ValueError(f"Unknown ic='{ic}'; use 'tg' or 'adv'.")
    if n_steps <= 0:
        raise ValueError(f"n_steps must be > 0, got {n_steps}")

    if nu is None:
        nu = r2._NU_DEFAULT

    # Seed set and logged for reproducibility provenance (see
    # alignment_bridge.run_bridge_experiment for the identical rationale).
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
            try:
                dd = depletion_diagnostics(
                    solver["u"], solver["v"], solver["w"], wx, wy, wz, dx
                )
            except ValueError:
                # Vorticity has (numerically) vanished; record a degenerate,
                # explicitly-flagged zero point rather than dropping the step.
                dd = {
                    "M": 0.0, "x_star": (0, 0, 0),
                    "alpha_at_max": 0.0, "alpha_sup_high": 0.0,
                    "depletion_ratio": 0.0, "theta_eff": float("nan"),
                    "enstrophy": 0.0,
                }
            record = {"step": step_idx, "t": solver["t"], **dd}
            timeseries.append(record)
            if verbose:
                print(
                    f"  step {step_idx:4d}  t={record['t']:.4f}  M={dd['M']:.4e}  "
                    f"alpha_max={dd['alpha_at_max']:.4e}  alpha_sup={dd['alpha_sup_high']:.4e}  "
                    f"ratio={dd['depletion_ratio']:.4e}  theta_eff={dd['theta_eff']:.4f}"
                )

    wall = time.perf_counter() - t0

    ratios = np.array([rec["depletion_ratio"] for rec in timeseries], dtype=float)
    M_vals = np.array([rec["M"] for rec in timeseries], dtype=float)
    alpha_max_vals = np.array([rec["alpha_at_max"] for rec in timeseries], dtype=float)
    alpha_sup_vals = np.array([rec["alpha_sup_high"] for rec in timeseries], dtype=float)
    theta_vals = np.array([rec["theta_eff"] for rec in timeseries], dtype=float)
    enstrophy_vals = np.array([rec["enstrophy"] for rec in timeseries], dtype=float)

    all_finite = bool(
        np.all(np.isfinite(ratios))
        and np.all(np.isfinite(M_vals))
        and np.all(np.isfinite(alpha_max_vals))
        and np.all(np.isfinite(alpha_sup_vals))
    )
    median_ratio = float(np.median(ratios)) if ratios.size else float("nan")
    max_ratio = float(np.max(ratios)) if ratios.size else float("nan")
    final_ratio = float(ratios[-1]) if ratios.size else float("nan")

    verdict = "PASS" if (
        all_finite and np.isfinite(median_ratio) and median_ratio < 1.0
    ) else "FAIL"

    finite_theta = theta_vals[np.isfinite(theta_vals)]
    theta_eff_median = float(np.median(finite_theta)) if finite_theta.size else float("nan")

    tag = _IC_TAGS[ic]
    exp_id = exp_id_override or f"EXP-L4-DEPLETION-{tag}-001"
    claim_id = claim_id_override or f"depletion-h2-{ic}"

    M_max = float(np.max(M_vals)) if M_vals.size else 0.0
    M_min = float(np.min(M_vals)) if M_vals.size else 0.0

    alpha_max_med = float(np.median(alpha_max_vals)) if alpha_max_vals.size else float("nan")
    alpha_sup_med = float(np.median(alpha_sup_vals)) if alpha_sup_vals.size else float("nan")
    enstrophy_final = float(enstrophy_vals[-1]) if enstrophy_vals.size else float("nan")

    result = {
        "exp_id": exp_id, "claim_id": claim_id,
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic, "N": N, "n_steps": n_steps, "seed": seed, "nu": nu,
        "eps_param": eps_param, "t_final": solver["t"],
        "M_max": M_max, "M_min": M_min,
        "alpha_at_max_median": alpha_max_med, "alpha_sup_high_median": alpha_sup_med,
        "depletion_ratio_max": max_ratio, "depletion_ratio_median": median_ratio,
        "depletion_ratio_final": final_ratio, "theta_eff_median": theta_eff_median,
        "enstrophy_final": enstrophy_final, "verdict": verdict, "wall_time_s": wall,
        "key_metric": (
            f"median_depletion_ratio={median_ratio:.4f} max_depletion_ratio={max_ratio:.4f} "
            f"final_depletion_ratio={final_ratio:.4f} M_range=[{M_min:.3e},{M_max:.3e}]"
        ),
        "timeseries": timeseries,
    }

    if verbose:
        print(
            f"\n  H2 Depletion ({ic}, N={N}, n_steps={n_steps})  "
            f"M=[{M_min:.4e},{M_max:.4e}]  median={median_ratio:.4f}  "
            f"max={max_ratio:.4f}  final={final_ratio:.4f}  "
            f"theta_eff_med={theta_eff_median}  verdict={verdict}"
        )

    if db_path is not None:
        log_result(db_path, result)

    return result


# ── Results logging (layer4 table pattern) — retry-on-locked backoff, since
# other agents may write results.db concurrently ────────────────────────────

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path, timeout=30.0)
    conn.execute(
        "CREATE TABLE IF NOT EXISTS layer4_depletion_experiments ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, exp_id TEXT UNIQUE NOT NULL, "
        "claim_id TEXT NOT NULL, timestamp TEXT NOT NULL, ic TEXT, "
        "N INTEGER, n_steps INTEGER, seed INTEGER, nu REAL, eps_param REAL, "
        "t_final REAL, M_max REAL, M_min REAL, alpha_at_max_median REAL, "
        "alpha_sup_high_median REAL, depletion_ratio_max REAL, "
        "depletion_ratio_median REAL, depletion_ratio_final REAL, "
        "theta_eff_median REAL, enstrophy_final REAL, verdict TEXT NOT NULL, "
        "wall_time_s REAL, key_metric TEXT)"
    )
    conn.execute(
        "CREATE TABLE IF NOT EXISTS layer4_depletion_timeseries ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, exp_id TEXT NOT NULL, "
        "step INTEGER, t REAL, M REAL, alpha_at_max REAL, "
        "alpha_sup_high REAL, depletion_ratio REAL, theta_eff REAL, "
        "enstrophy REAL)"
    )
    conn.commit()
    return conn


def _write_result(db_path: str, result: dict) -> None:
    """Single-attempt DB write (summary UPSERT + append-only timeseries)."""
    conn = _ensure_db(db_path)
    try:
        cols = [
            "exp_id", "claim_id", "timestamp", "ic", "N", "n_steps", "seed",
            "nu", "eps_param", "t_final", "M_max", "M_min",
            "alpha_at_max_median", "alpha_sup_high_median",
            "depletion_ratio_max", "depletion_ratio_median",
            "depletion_ratio_final", "theta_eff_median", "enstrophy_final",
            "verdict", "wall_time_s", "key_metric",
        ]
        values = tuple(result.get(c) for c in cols)
        ph = ", ".join("?" for _ in cols)
        conn.execute(
            f"INSERT OR REPLACE INTO layer4_depletion_experiments ({', '.join(cols)}) "
            f"VALUES ({ph})",
            values,
        )
        for rec in result.get("timeseries", []):
            conn.execute(
                """
                INSERT INTO layer4_depletion_timeseries
                (exp_id, step, t, M, alpha_at_max, alpha_sup_high,
                 depletion_ratio, theta_eff, enstrophy)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    result["exp_id"], rec["step"], rec["t"], rec["M"],
                    rec["alpha_at_max"], rec["alpha_sup_high"],
                    rec["depletion_ratio"], rec["theta_eff"], rec["enstrophy"],
                ),
            )
        conn.commit()
    finally:
        conn.close()


def log_result(
    db_path: str, result: dict,
    max_retries: int = _DB_MAX_RETRIES,
    backoff_base_s: float = _DB_BACKOFF_BASE_S,
) -> None:
    """Append one depletion experiment (+ timeseries) to results.db.

    Retries with short exponential backoff on sqlite3.OperationalError
    ("database is locked") — other layer4 experiment scripts may be
    writing to the same results.db concurrently. Re-raises immediately on
    any other error, or once retries are exhausted.
    """
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
    parser = argparse.ArgumentParser(description="H2 vortex-stretching depletion diagnostic")
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

    r = run_depletion_experiment(
        ic=args.ic, N=args.N, n_steps=args.n_steps,
        record_every=args.record_every, seed=args.seed, nu=args.nu,
        eps_param=args.eps, verbose=args.verbose, db_path=args.db,
    )
    print(f"\nVerdict: {r['verdict']}  key_metric: {r['key_metric']}")


if __name__ == "__main__":
    main()
