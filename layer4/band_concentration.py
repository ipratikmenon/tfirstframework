"""
layer4/band_concentration.py — Band Concentration Diagnostic (S41 "core term")
================================================================================
Numerical stress-test for Conjecture target:band-absorption
(proofs/claim_a_3d_proof_attempt.tex, \\S sec:annulus-s41, session S41) — the
"core term" left open after the annulus level iteration (Remark
rem:level-iteration) was PROVED not to close (Proposition s41-nocontraction,
s41-coarea). See PROGRESS.md's S41 KEY FINDING.

The conjecture asserts that the radial-angular coupling density on the
transition band

    Band = {M/8 <= |omega| <= M/4}

(Definition def:annulus-residual's Ann) is NOT concentrated there relative
to its typical value on the broader high-vorticity region

    E = {|omega| >= M/4}

(the identical E used throughout layer4/k_correlation.py's c_K
diagnostic — the conjecture is explicitly stated in the proof doc to be
"the same analytic type as Conjecture conj:K-refined", i.e. the c_K
decorrelation hypothesis, localized to the band).

Two band-concentration ratios test the two inequalities of the conjecture:

  conc_Y = mean_Band(|grad m|^2 * w) / mean_E(|grad m|^2 * w)
           tests eq:s41-bandY: concentration of the (S1)-family integrand
           |grad m|^2 * w -- the term S41 proved CANNOT be absorbed by any
           choice of Young parameter in the level-iteration recursion.

  conc_D = mean_Band(|grad^2 ehat|^2) / mean_E(|grad^2 ehat|^2)
           tests eq:s41-bandD: concentration of the angular-dissipation
           density that the (S2)-family's Young split reduces to.

Reading convention (identical to c_K, S34-S35): conc_Y, conc_D ~ O(1) with
no growth across an amplitude/viscosity sweep is NUMERICAL SUPPORT for the
conjecture, not a proof -- the program never proved c_K <= C/L either, only
found c_K ~ 1 with no Re-growth and reported that as "supported." A large
or growing conc_Y/conc_D is a genuine, useful negative finding (a concrete
lead toward an adversarial counterexample for future analytical work), not
a failure of this diagnostic.

Reuses (does not reimplement): layer4.alignment_bridge's
compute_direction_field, grad_ehat_sq, _validate_omega_grid,
_vorticity_from_velocity, adv_vortex_tubes_ic; layer4.k_correlation's
grad_scalar_sq; layer4.geometric_disorder's spectral_wavenumbers;
layer3.route2_3D's make_solver, set_ic, taylor_green_ic, shear_layer_ic,
compute_cfl_dt, solver_step (exact unforced Prize NS, f=0 throughout, same
solver every other layer4 diagnostic in this program uses).

House rules: spectral (FFT) derivatives only; pure functions (no mutation
of inputs); every logged result carries claim_id + key_metric + verdict.
"""

from __future__ import annotations

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
from layer4.geometric_disorder import spectral_wavenumbers  # noqa: E402
from layer4.alignment_bridge import (  # noqa: E402
    _validate_omega_grid,
    compute_direction_field,
    grad_ehat_sq,
    _vorticity_from_velocity,
    adv_vortex_tubes_ic,
)
from layer4.k_correlation import grad_scalar_sq  # noqa: E402

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_SEED = 42
_IC_TAGS = {"tg": "TG", "shear": "SH", "adv": "ADV"}

# Definition def:annulus-residual's transition band: Ann = {M/8 <= |omega| <= M/4}.
DEFAULT_THETA_MINUS = 0.125
DEFAULT_THETA_PLUS = 0.25

# Convention (documented, not derived): degenerate covariance (empty Band or
# empty E, or an identically-zero denominator) reports conc = 0.0 -- same
# "trivially satisfied" reading as k_correlation.py's c_K := 0 convention.
_ZERO_TOL = 1e-12

# Same PASS threshold convention as k_correlation.py's _KCORR_PASS_MEDIAN_CK:
# an arbitrary but documented "stays O(1) or below" bar, not a proved bound.
_BAND_PASS_MEDIAN_CONC = 10.0


# =============================================================================
# Second-derivative spectral primitives (new; siblings of
# k_correlation.grad_scalar_sq and alignment_bridge.grad_ehat_sq)
# =============================================================================

def grad2_scalar_sq(
    field: np.ndarray, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray
) -> np.ndarray:
    """|grad^2 f|^2 = sum_{j,k in {x,y,z}} (d_j d_k f)^2 for a scalar field f.

    Spectral: d_j d_k f = ifft(-k_j*k_k * fft(f)). Sums all 9 ordered pairs
    (j,k) — the Frobenius norm squared of the Hessian, same "sum over all
    ordered index pairs" convention as grad_ehat_sq's sum over (i,k).

    Does not mutate field.
    """
    f_hat = np.fft.fftn(field)
    ks = (kx, ky, kz)
    total = np.zeros_like(field, dtype=float)
    for kj in ks:
        for kk in ks:
            d2 = np.real(np.fft.ifftn(-kj * kk * f_hat))
            total = total + d2**2
    return total


def grad2_ehat_sq(
    ex: np.ndarray, ey: np.ndarray, ez: np.ndarray,
    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray,
) -> np.ndarray:
    """|grad^2 ehat|^2 = sum_i |grad^2 ehat_i|^2, spectral second derivatives.

    Same Gibbs-artifact caveat as grad_ehat_sq (ehat is truncated to zero
    outside the valid vorticity set): read only on masks well inside the
    high-vorticity region, as band_concentration_diagnostics does.

    Does not mutate ex, ey, ez.
    """
    total = np.zeros_like(ex, dtype=float)
    for comp in (ex, ey, ez):
        total = total + grad2_scalar_sq(comp, kx, ky, kz)
    return total


# =============================================================================
# Core band-concentration diagnostic
# =============================================================================

def band_concentration_diagnostics(
    wx: np.ndarray, wy: np.ndarray, wz: np.ndarray, dx: float,
    theta_minus: float = DEFAULT_THETA_MINUS,
    theta_plus: float = DEFAULT_THETA_PLUS,
) -> dict:
    """Compute the band-concentration diagnostic dict for one vorticity snapshot.

    Parameters
    ----------
    wx, wy, wz   : vorticity components, shape (N, N, N), periodic domain [0,2π)^3
    dx           : grid spacing, must equal 2π/N (validated)
    theta_minus, theta_plus : band = {theta_minus*M <= |omega| <= theta_plus*M};
                   defaults are Definition def:annulus-residual's (1/8, 1/4).
                   Must satisfy 0 < theta_minus < theta_plus <= 1.

    Returns dict with keys: conc_Y, conc_D, mean_E_Y, mean_E_D, mean_Band_Y,
    mean_Band_D, E_volume_fraction, band_volume_fraction, theta_minus,
    theta_plus. Raises ValueError on shape mismatch, bad dx, an
    identically-zero vorticity field (same contract as bridge_diagnostics /
    k_correlation_diagnostics), or theta_minus >= theta_plus.

    Does not mutate wx, wy, wz.
    """
    if not (0.0 < theta_minus < theta_plus <= 1.0):
        raise ValueError(
            f"require 0 < theta_minus < theta_plus <= 1; got "
            f"theta_minus={theta_minus}, theta_plus={theta_plus}"
        )

    N = _validate_omega_grid(wx, wy, wz, dx)

    mag = np.sqrt(wx**2 + wy**2 + wz**2)
    M = float(np.max(mag))
    if M < 1e-14:
        raise ValueError(
            "band_concentration_diagnostics: vorticity field is identically "
            "zero (M ≈ 0); E and Band are undefined."
        )

    kx, ky, kz = spectral_wavenumbers(N)

    E = mag >= (theta_plus * M)
    Band = (mag >= (theta_minus * M)) & (mag <= (theta_plus * M))

    E_volume_fraction = float(np.count_nonzero(E)) / float(N**3)
    band_volume_fraction = float(np.count_nonzero(Band)) / float(N**3)

    grad_m2 = grad_scalar_sq(mag, kx, ky, kz)  # |∇m|²
    ex, ey, ez, _mag2, _valid = compute_direction_field(wx, wy, wz)
    w = grad_ehat_sq(ex, ey, ez, kx, ky, kz)  # |∇ê|²
    grad2_e = grad2_ehat_sq(ex, ey, ez, kx, ky, kz)  # |∇²ê|²

    Y_density = grad_m2 * w       # (S1)-family integrand, eq:s41-bandY
    D_density = grad2_e           # (S2)/angular-dissipation density, eq:s41-bandD

    mean_E_Y = float(np.mean(Y_density[E])) if np.any(E) else 0.0
    mean_E_D = float(np.mean(D_density[E])) if np.any(E) else 0.0
    mean_Band_Y = float(np.mean(Y_density[Band])) if np.any(Band) else 0.0
    mean_Band_D = float(np.mean(D_density[Band])) if np.any(Band) else 0.0

    conc_Y = mean_Band_Y / mean_E_Y if mean_E_Y > _ZERO_TOL else 0.0
    conc_D = mean_Band_D / mean_E_D if mean_E_D > _ZERO_TOL else 0.0

    return {
        "conc_Y": float(conc_Y),
        "conc_D": float(conc_D),
        "mean_E_Y": mean_E_Y,
        "mean_E_D": mean_E_D,
        "mean_Band_Y": mean_Band_Y,
        "mean_Band_D": mean_Band_D,
        "E_volume_fraction": E_volume_fraction,
        "band_volume_fraction": band_volume_fraction,
        "theta_minus": theta_minus,
        "theta_plus": theta_plus,
    }


# =============================================================================
# Production experiment: run the existing Route 2 3D solver, track the
# band-concentration diagnostic (mirrors alignment_bridge.run_bridge_experiment
# structurally; does not call it, to avoid re-running the solver twice)
# =============================================================================

def run_band_experiment(
    ic: str = "tg",
    N: int = 64,
    n_steps: int = 400,
    nu: float | None = None,
    eps_param: float = 0.1,
    record_every: int = 10,
    seed: int = _DEFAULT_SEED,
    cfl_safety: float = 0.4,
    theta_minus: float = DEFAULT_THETA_MINUS,
    theta_plus: float = DEFAULT_THETA_PLUS,
    verbose: bool = False,
    db_path: str | None = None,
    exp_id_override: str | None = None,
    claim_id_override: str | None = None,
    adv_amp: float | None = None,
) -> dict:
    """Run Route 2 3D (exact unforced Prize NS, f=0) and track the
    band-concentration diagnostic against Conjecture target:band-absorption.

    Verdict PASS iff every recorded conc_Y and conc_D is finite AND
    median(conc_Y) <= 10 AND median(conc_D) <= 10 -- i.e. the band density
    stays O(1) or below relative to the broader high-vorticity region,
    same reading convention as k_correlation.py's c_K threshold.
    """
    if ic not in _IC_TAGS:
        raise ValueError(f"Unknown ic='{ic}'; use 'tg', 'shear', or 'adv'.")
    if n_steps <= 0:
        raise ValueError(f"n_steps must be > 0, got {n_steps}")

    if nu is None:
        nu = r2._NU_DEFAULT

    solver = r2.make_solver(N=N, nu=nu, eps_param=eps_param)
    if ic == "tg":
        ic_dict = r2.taylor_green_ic(N)
    elif ic == "shear":
        ic_dict = r2.shear_layer_ic(N)
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
                bcd = band_concentration_diagnostics(
                    wx, wy, wz, dx, theta_minus, theta_plus
                )
            except ValueError:
                # Vorticity has (numerically) vanished; record a degenerate,
                # explicitly-flagged zero point rather than dropping the step.
                bcd = {
                    "conc_Y": 0.0, "conc_D": 0.0,
                    "mean_E_Y": 0.0, "mean_E_D": 0.0,
                    "mean_Band_Y": 0.0, "mean_Band_D": 0.0,
                    "E_volume_fraction": 0.0, "band_volume_fraction": 0.0,
                    "theta_minus": theta_minus, "theta_plus": theta_plus,
                }
            record = {"step": step_idx, "t": solver["t"], **bcd}
            timeseries.append(record)
            if verbose:
                print(
                    f"  step {step_idx:4d}  t={record['t']:.4f}  "
                    f"conc_Y={bcd['conc_Y']:.4e}  conc_D={bcd['conc_D']:.4e}  "
                    f"band_vol_frac={bcd['band_volume_fraction']:.4e}"
                )

    wall = time.perf_counter() - t0

    conc_Y_vals = np.array([rec["conc_Y"] for rec in timeseries], dtype=float)
    conc_D_vals = np.array([rec["conc_D"] for rec in timeseries], dtype=float)
    band_vol_vals = np.array(
        [rec["band_volume_fraction"] for rec in timeseries], dtype=float
    )

    all_finite = bool(np.all(np.isfinite(conc_Y_vals)) and np.all(np.isfinite(conc_D_vals)))
    conc_Y_max = float(np.max(conc_Y_vals)) if conc_Y_vals.size else float("nan")
    conc_Y_median = float(np.median(conc_Y_vals)) if conc_Y_vals.size else float("nan")
    conc_Y_final = float(conc_Y_vals[-1]) if conc_Y_vals.size else float("nan")
    conc_D_max = float(np.max(conc_D_vals)) if conc_D_vals.size else float("nan")
    conc_D_median = float(np.median(conc_D_vals)) if conc_D_vals.size else float("nan")
    conc_D_final = float(conc_D_vals[-1]) if conc_D_vals.size else float("nan")
    band_vol_median = float(np.median(band_vol_vals)) if band_vol_vals.size else float("nan")

    verdict = "PASS" if (
        all_finite
        and np.isfinite(conc_Y_median) and conc_Y_median <= _BAND_PASS_MEDIAN_CONC
        and np.isfinite(conc_D_median) and conc_D_median <= _BAND_PASS_MEDIAN_CONC
    ) else "FAIL"

    tag = _IC_TAGS[ic]
    exp_id = exp_id_override or f"EXP-L4-BAND-{tag}-001"
    claim_id = claim_id_override or f"target-band-absorption-{ic}"

    result = {
        "exp_id": exp_id,
        "claim_id": claim_id,
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic,
        "N": N,
        "n_steps": n_steps,
        "seed": seed,
        "nu": nu,
        "theta_minus": theta_minus,
        "theta_plus": theta_plus,
        "t_final": solver["t"],
        "conc_Y_max": conc_Y_max,
        "conc_Y_median": conc_Y_median,
        "conc_Y_final": conc_Y_final,
        "conc_D_max": conc_D_max,
        "conc_D_median": conc_D_median,
        "conc_D_final": conc_D_final,
        "band_volume_fraction_median": band_vol_median,
        "verdict": verdict,
        "wall_time_s": wall,
        "key_metric": (
            f"median_conc_Y={conc_Y_median:.4f} median_conc_D={conc_D_median:.4f} "
            f"max_conc_Y={conc_Y_max:.4f} max_conc_D={conc_D_max:.4f}"
        ),
        "timeseries": timeseries,
    }

    if verbose:
        print(f"\n  Band-concentration experiment ({ic}, N={N}, n_steps={n_steps})")
        print(f"    median conc_Y: {conc_Y_median:.4f}")
        print(f"    median conc_D: {conc_D_median:.4f}")
        print(f"    max conc_Y:    {conc_Y_max:.4f}")
        print(f"    max conc_D:    {conc_D_max:.4f}")
        print(f"    verdict:       {verdict}")

    if db_path is not None:
        log_result(db_path, result)

    return result


# =============================================================================
# Results logging (layer4-module table pattern; see alignment_bridge.py,
# k_correlation.py's use of results.db)
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS band_concentration_experiments (
            id                          INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id                      TEXT UNIQUE NOT NULL,
            claim_id                    TEXT,
            timestamp                   TEXT,
            ic                          TEXT,
            N                           INTEGER,
            n_steps                     INTEGER,
            seed                        INTEGER,
            nu                          REAL,
            theta_minus                 REAL,
            theta_plus                  REAL,
            t_final                     REAL,
            conc_Y_max                  REAL,
            conc_Y_median               REAL,
            conc_Y_final                REAL,
            conc_D_max                  REAL,
            conc_D_median               REAL,
            conc_D_final                REAL,
            band_volume_fraction_median REAL,
            verdict                     TEXT,
            wall_time_s                 REAL,
            key_metric                  TEXT
        )
    """)
    conn.commit()
    return conn


def log_result(db_path: str, result: dict) -> None:
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "ic", "N", "n_steps", "seed", "nu",
        "theta_minus", "theta_plus", "t_final", "conc_Y_max", "conc_Y_median",
        "conc_Y_final", "conc_D_max", "conc_D_median", "conc_D_final",
        "band_volume_fraction_median", "verdict", "wall_time_s", "key_metric",
    ]
    vals = tuple(result.get(c) for c in cols)
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO band_concentration_experiments "
        f"({', '.join(cols)}) VALUES ({ph})",
        vals,
    )
    conn.commit()
    conn.close()


# =============================================================================
# CLI
# =============================================================================

def main() -> None:  # pragma: no cover
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--ic", default="tg", choices=["tg", "shear", "adv", "all"])
    p.add_argument("--N", type=int, default=64)
    p.add_argument("--n-steps", type=int, default=400)
    p.add_argument("--nu", type=float, default=None)
    p.add_argument("--db", default=_DEFAULT_DB)
    args = p.parse_args()

    ics = ["tg", "shear", "adv"] if args.ic == "all" else [args.ic]
    for ic in ics:
        run_band_experiment(
            ic=ic, N=args.N, n_steps=args.n_steps, nu=args.nu,
            verbose=True, db_path=args.db,
        )


if __name__ == "__main__":  # pragma: no cover
    main()
