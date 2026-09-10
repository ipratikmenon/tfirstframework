"""
layer4/scale_diagnostics.py — Diagnostics AT SCALE r* (S44 follow-up)
======================================================================
S44 (`rem:s44-numerics`, proofs/claim_a_3d_proof_attempt.tex §sec:relation-s44)
names three quantities that bear on Conjecture `conj:lc-s44` / `conj:K-refined`
/ `target:band-absorption` and that this program has NOT measured:

  1. |E|/|Q(r*)| and |Band|/|Q(r*)|  -- AT SCALE r*, not over the whole box.
     S44 flags explicitly that S42-S43's `band_concentration.py` volume
     fractions are taken over the whole periodic box and are therefore
     "motivation and not evidence" for the volume hypotheses (NB)/(NB-E)
     (`eq:s44-NB`, `eq:s44-NBE`).
  2. c_K^band -- the correlation constant on the transition band, which
     `prop:s44-linfty-route` shows must satisfy the DECAYING bound
     c_K^band <= C/L for the sup-norm route to close.
  3. The L-DEPENDENCE of both correlation constants. This is the sharp
     point: S34-S35 measured c_K ~ 1 and explicitly could not distinguish
     "c_K = O(1)" from the stronger "c_K <= C/L" that the closure actually
     needs. Boundedness is not enough; decay is what is required.

This module measures all three, restricted to the ball B_{r*}(x*) around the
vorticity maximum -- the spatial section of the parabolic cylinder
Q(r*) = B_{r*}(x*) x (t0-(r*)^2, t0).

SCOPE CAVEAT (stated up front, not buried): a single snapshot gives the
SPATIAL ball B_{r*}(x*), not the full space-time cylinder Q(r*). The volume
fractions reported here are |E n B_{r*}|/|B_{r*}|, which coincide with the
space-time fractions only if they are roughly steady in t. The per-snapshot
time series is therefore recorded so that steadiness can be checked rather
than assumed.

The L-dependence test:
  L = log(e + ||u||_{H^3}/M)  (definition at proof doc line ~5671)
  varies along a run (M and ||u||_{H^3} both evolve) and across a viscosity
  sweep. Regressing log(c_K) on log(L) over the collected samples gives a
  slope:
      slope ~ -1  ->  consistent with the decaying form c_K <= C/L
      slope ~  0  ->  only boundedness, c_K = O(1); the closure route of
                      `prop:s44-linfty-route` would NOT be supported
  This distinction is the whole point; a PASS on boundedness alone is
  explicitly NOT evidence for the conjecture as stated.

Reuses (does not reimplement): alignment_bridge's compute_direction_field,
grad_ehat_sq, _validate_omega_grid, _vorticity_from_velocity,
adv_vortex_tubes_ic, bridge_diagnostics; k_correlation's grad_scalar_sq;
geometric_disorder's spectral_wavenumbers and _ball_mask; route2_3D's solver.

House rules: spectral (FFT) derivatives only; pure functions; every logged
result carries claim_id + key_metric + verdict.
"""

from __future__ import annotations

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

import route2_3D as r2  # noqa: E402
from layer4.geometric_disorder import spectral_wavenumbers, _ball_mask  # noqa: E402
from layer4.alignment_bridge import (  # noqa: E402
    _validate_omega_grid,
    compute_direction_field,
    grad_ehat_sq,
    _vorticity_from_velocity,
    adv_vortex_tubes_ic,
    bridge_diagnostics,
)
from layer4.k_correlation import grad_scalar_sq  # noqa: E402

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_SEED = 42
_IC_TAGS = {"tg": "TG", "shear": "SH", "adv": "ADV"}

THETA_MINUS = 0.125   # Band lower level (Definition def:annulus-residual)
THETA_PLUS = 0.25     # Band upper level == E lower level
_ZERO_TOL = 1e-12


def sobolev_H3_norm(u: np.ndarray, v: np.ndarray, w: np.ndarray,
                    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> float:
    """||u||_{H^3} for the velocity field, spectrally.

    ||u||^2_{H^3} = sum_k (1+|k|^2)^3 |uhat(k)|^2, summed over the three
    components, with the FFT normalisation that makes this the continuum
    norm on [0,2pi)^3. Does not mutate inputs.
    """
    N = u.shape[0]
    k2 = kx**2 + ky**2 + kz**2
    weight = (1.0 + k2) ** 3
    total = 0.0
    for comp in (u, v, w):
        ch = np.fft.fftn(comp) / (N**3)
        total += float(np.sum(weight * (np.abs(ch) ** 2)))
    return float(np.sqrt(max(total, 0.0)))


def L_factor(u: np.ndarray, v: np.ndarray, w: np.ndarray, M: float,
             kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> float:
    """L = log(e + ||u||_{H^3}/M), the logarithmic factor of the proof document.

    Returns log(e) = 1 in the degenerate case M <= 0.
    """
    if M <= _ZERO_TOL:
        return 1.0
    return float(np.log(np.e + sobolev_H3_norm(u, v, w, kx, ky, kz) / M))


def _corr_ratio(gm2: np.ndarray, w: np.ndarray, mask: np.ndarray) -> float:
    """c_K[1_R] = <gm2*w>_R / (<gm2>_R <w>_R) on the region R = mask.

    Same degenerate-covariance convention as k_correlation.py: returns 0.0
    when either factor averages to (numerically) zero on R, or R is empty.
    """
    if not np.any(mask):
        return 0.0
    a = float(np.mean(gm2[mask]))
    b = float(np.mean(w[mask]))
    if a <= _ZERO_TOL or b <= _ZERO_TOL:
        return 0.0
    return float(np.mean(gm2[mask] * w[mask]) / (a * b))


def scale_diagnostics(u: np.ndarray, v: np.ndarray, w_vel: np.ndarray,
                      dx: float) -> dict:
    """Diagnostics restricted to the ball B_{r*}(x*), from a velocity snapshot.

    Returns dict with:
      M, r_star, L
      frac_E_rstar     = |E n B_{r*}| / |B_{r*}|          (for (NB-E))
      frac_band_rstar  = |Band n B_{r*}| / |B_{r*}|       (for (NB))
      c_K_rstar        = correlation constant on E n B_{r*}
      c_K_band         = correlation constant on Band n B_{r*}   [S44 item 2]
      cK_times_L, cKband_times_L   -- the products that must stay BOUNDED
                                      if the decaying form c_K <= C/L holds
      ball_cells       -- number of grid cells in B_{r*} (resolution guard)

    Raises ValueError on an identically-zero vorticity field (same contract
    as bridge_diagnostics). Does not mutate inputs.
    """
    N = u.shape[0]
    kx, ky, kz = spectral_wavenumbers(N)
    wx, wy, wz = _vorticity_from_velocity(u, v, w_vel, kx, ky, kz)
    _validate_omega_grid(wx, wy, wz, dx)

    mag = np.sqrt(wx**2 + wy**2 + wz**2)
    M = float(np.max(mag))
    if M < 1e-14:
        raise ValueError("scale_diagnostics: vorticity identically zero; "
                         "r* and E are undefined.")

    r_star = M ** -0.5
    idx = np.unravel_index(int(np.argmax(mag)), mag.shape)
    ball = _ball_mask(N, (int(idx[0]), int(idx[1]), int(idx[2])), r_star)
    ball_cells = int(np.count_nonzero(ball))

    E = ball & (mag >= THETA_PLUS * M)
    Band = ball & (mag >= THETA_MINUS * M) & (mag <= THETA_PLUS * M)

    gm2 = grad_scalar_sq(mag, kx, ky, kz)
    ex, ey, ez, _m2, _valid = compute_direction_field(wx, wy, wz)
    wang = grad_ehat_sq(ex, ey, ez, kx, ky, kz)

    c_K_rstar = _corr_ratio(gm2, wang, E)
    c_K_band = _corr_ratio(gm2, wang, Band)
    L = L_factor(u, v, w_vel, M, kx, ky, kz)

    denom = float(max(ball_cells, 1))
    return {
        "M": M,
        "r_star": float(r_star),
        "L": L,
        "ball_cells": ball_cells,
        "frac_E_rstar": float(np.count_nonzero(E)) / denom,
        "frac_band_rstar": float(np.count_nonzero(Band)) / denom,
        "c_K_rstar": c_K_rstar,
        "c_K_band": c_K_band,
        "cK_times_L": c_K_rstar * L,
        "cKband_times_L": c_K_band * L,
    }


def run_scale_experiment(
    ic: str = "tg",
    N: int = 64,
    n_steps: int = 200,
    nu: float | None = None,
    eps_param: float = 0.1,
    record_every: int = 10,
    seed: int = _DEFAULT_SEED,
    cfl_safety: float = 0.4,
    verbose: bool = False,
    db_path: str | None = None,
    exp_id_override: str | None = None,
    claim_id_override: str | None = None,
    ic_dict_override: dict | None = None,
) -> dict:
    """Run Route 2 3D (exact unforced Prize NS, f=0) tracking scale-r* diagnostics.

    `ic_dict_override` lets a caller supply a custom IC dict (used by the
    jet strain sweep) while reusing this runner unchanged.

    Verdict semantics (deliberately weak, and named as such): PASS means only
    that every recorded quantity is finite and the correlation constants stay
    bounded. It is NOT evidence for the decaying form c_K <= C/L -- that is
    what the L-regression in `l_dependence_summary` addresses.
    """
    if ic_dict_override is None and ic not in _IC_TAGS:
        raise ValueError(f"Unknown ic='{ic}'; use 'tg', 'shear', or 'adv'.")
    if n_steps <= 0:
        raise ValueError(f"n_steps must be > 0, got {n_steps}")
    if nu is None:
        nu = r2._NU_DEFAULT

    solver = r2.make_solver(N=N, nu=nu, eps_param=eps_param)
    if ic_dict_override is not None:
        ic_dict = ic_dict_override
    elif ic == "tg":
        ic_dict = r2.taylor_green_ic(N)
    elif ic == "shear":
        ic_dict = r2.shear_layer_ic(N)
    else:
        ic_dict = adv_vortex_tubes_ic(N, seed=seed)
    solver = r2.set_ic(solver, ic_dict)

    dx = 2.0 * np.pi / N
    timeseries: list[dict] = []
    t0 = time.perf_counter()

    for step_idx in range(n_steps):
        dt = r2.compute_cfl_dt(solver, safety=cfl_safety)
        solver, _diag = r2.solver_step(solver, dt)
        if step_idx % record_every == 0 or step_idx == n_steps - 1:
            try:
                sd = scale_diagnostics(solver["u"], solver["v"], solver["w"], dx)
            except ValueError:
                continue
            try:
                bd = bridge_diagnostics(
                    *_vorticity_from_velocity(
                        solver["u"], solver["v"], solver["w"],
                        solver["kx"], solver["ky"], solver["kz"]),
                    dx)
                sd["sigma_star"] = bd["sigma_star"]
                sd["bridge_ratio"] = bd["bridge_ratio"]
            except ValueError:
                sd["sigma_star"] = 0.0
                sd["bridge_ratio"] = 0.0
            sd["step"] = step_idx
            sd["t"] = solver["t"]
            timeseries.append(sd)
            if verbose:
                print(f"  step {step_idx:4d} t={sd['t']:.3f} L={sd['L']:.3f} "
                      f"fE={sd['frac_E_rstar']:.4f} fB={sd['frac_band_rstar']:.4f} "
                      f"cK={sd['c_K_rstar']:.3f} cKb={sd['c_K_band']:.3f} "
                      f"sig*={sd['sigma_star']:.3e}")

    wall = time.perf_counter() - t0
    if not timeseries:
        raise RuntimeError("run_scale_experiment: no samples recorded")

    def med(key):
        return float(np.median([r[key] for r in timeseries]))

    def mx(key):
        return float(np.max([r[key] for r in timeseries]))

    all_finite = all(
        np.isfinite([r[k] for k in ("L", "c_K_rstar", "c_K_band",
                                    "frac_E_rstar", "frac_band_rstar")]).all()
        for r in timeseries
    )
    verdict = "PASS" if (all_finite and med("c_K_rstar") <= 10.0
                         and med("c_K_band") <= 10.0) else "FAIL"

    tag = _IC_TAGS.get(ic, ic.upper())
    exp_id = exp_id_override or f"EXP-L4-SCALE-{tag}-001"
    claim_id = claim_id_override or f"s44-scale-diagnostics-{ic}"

    result = {
        "exp_id": exp_id,
        "claim_id": claim_id,
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic,
        "N": N,
        "n_steps": n_steps,
        "seed": seed,
        "nu": nu,
        "n_samples": len(timeseries),
        "L_median": med("L"),
        "L_max": mx("L"),
        "frac_E_rstar_median": med("frac_E_rstar"),
        "frac_band_rstar_median": med("frac_band_rstar"),
        "c_K_rstar_median": med("c_K_rstar"),
        "c_K_rstar_max": mx("c_K_rstar"),
        "c_K_band_median": med("c_K_band"),
        "c_K_band_max": mx("c_K_band"),
        "cK_times_L_median": med("cK_times_L"),
        "cK_times_L_max": mx("cK_times_L"),
        "cKband_times_L_median": med("cKband_times_L"),
        "cKband_times_L_max": mx("cKband_times_L"),
        "ball_cells_median": med("ball_cells"),
        "sigma_star_median": med("sigma_star"),
        "verdict": verdict,
        "wall_time_s": wall,
        "key_metric": (
            f"fE={med('frac_E_rstar'):.4f} fB={med('frac_band_rstar'):.4f} "
            f"cK={med('c_K_rstar'):.3f} cKband={med('c_K_band'):.3f} "
            f"L={med('L'):.3f}"
        ),
        "timeseries": timeseries,
    }

    if verbose:
        print(f"\n  scale-r* experiment ({ic}, N={N}, {len(timeseries)} samples)")
        print(f"    |E|/|B_r*| median:    {med('frac_E_rstar'):.4f}")
        print(f"    |Band|/|B_r*| median: {med('frac_band_rstar'):.4f}")
        print(f"    c_K at r* median:     {med('c_K_rstar'):.4f}")
        print(f"    c_K^band median:      {med('c_K_band'):.4f}")
        print(f"    L median:             {med('L'):.4f}")
        print(f"    ball cells median:    {med('ball_cells'):.0f}")
        print(f"    verdict:              {verdict}")

    if db_path is not None:
        log_result(db_path, result)
    return result


def l_dependence_summary(results: list[dict]) -> dict:
    """Regress log(c_K) on log(L) across pooled samples from several runs.

    THE point of this module. `prop:s44-linfty-route` needs the DECAYING
    bound c_K <= C/L, not mere boundedness; S34-S35 could not distinguish
    them. Pooling samples across runs (and, ideally, across a viscosity
    sweep, which is what moves L) gives:

        slope ~ -1  -> consistent with c_K ~ C/L  (decay; supports the route)
        slope ~  0  -> boundedness only           (does NOT support it)

    Samples with a degenerate (zero) correlation constant are dropped, and
    the count of dropped samples is reported -- they carry no information
    about the exponent and silently including them would bias the fit toward
    the boundedness reading.
    """
    Ls, cKs, cKbs = [], [], []
    for r in results:
        for s in r["timeseries"]:
            Ls.append(s["L"]); cKs.append(s["c_K_rstar"]); cKbs.append(s["c_K_band"])
    Ls = np.array(Ls); cKs = np.array(cKs); cKbs = np.array(cKbs)

    def fit(y):
        ok = (y > _ZERO_TOL) & (Ls > 1.0 + 1e-9) & np.isfinite(y) & np.isfinite(Ls)
        n_used, n_drop = int(ok.sum()), int((~ok).sum())
        if n_used < 3:
            return {"slope": float("nan"), "r2": float("nan"),
                    "n_used": n_used, "n_dropped": n_drop}
        X = np.log(Ls[ok]); Y = np.log(y[ok])
        A = np.vstack([X, np.ones_like(X)]).T
        coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
        pred = A @ coef
        ss_res = float(np.sum((Y - pred) ** 2))
        ss_tot = float(np.sum((Y - Y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        return {"slope": float(coef[0]), "intercept": float(coef[1]),
                "r2": float(r2), "n_used": n_used, "n_dropped": n_drop}

    return {
        "c_K_vs_L": fit(cKs),
        "c_K_band_vs_L": fit(cKbs),
        "L_range": [float(np.min(Ls)), float(np.max(Ls))] if Ls.size else [np.nan, np.nan],
    }


def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS scale_diagnostics_experiments (
            id                      INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id                  TEXT UNIQUE NOT NULL,
            claim_id                TEXT,
            timestamp               TEXT,
            ic                      TEXT,
            N                       INTEGER,
            n_steps                 INTEGER,
            seed                    INTEGER,
            nu                      REAL,
            n_samples               INTEGER,
            L_median                REAL,
            L_max                   REAL,
            frac_E_rstar_median     REAL,
            frac_band_rstar_median  REAL,
            c_K_rstar_median        REAL,
            c_K_rstar_max           REAL,
            c_K_band_median         REAL,
            c_K_band_max            REAL,
            cK_times_L_median       REAL,
            cK_times_L_max          REAL,
            cKband_times_L_median   REAL,
            cKband_times_L_max      REAL,
            ball_cells_median       REAL,
            sigma_star_median       REAL,
            verdict                 TEXT,
            wall_time_s             REAL,
            key_metric              TEXT
        )
    """)
    conn.commit()
    return conn


def log_result(db_path: str, result: dict) -> None:
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "ic", "N", "n_steps", "seed", "nu",
        "n_samples", "L_median", "L_max", "frac_E_rstar_median",
        "frac_band_rstar_median", "c_K_rstar_median", "c_K_rstar_max",
        "c_K_band_median", "c_K_band_max", "cK_times_L_median",
        "cK_times_L_max", "cKband_times_L_median", "cKband_times_L_max",
        "ball_cells_median", "sigma_star_median", "verdict", "wall_time_s",
        "key_metric",
    ]
    vals = tuple(result.get(c) for c in cols)
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO scale_diagnostics_experiments "
        f"({', '.join(cols)}) VALUES ({ph})", vals)
    conn.commit()
    conn.close()
