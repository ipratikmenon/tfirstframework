"""
layer4/k_correlation.py — K-Absorption Correlation Diagnostic (c_K)
================================================================
Numerical diagnostic for the K-absorption conjecture (S34): the first
Moser rung of the Bridge Lemma (§moser-s33 in the proof doc) reduces
Obstacle (i) to absorbing the radial-angular coupling term

    K = M^-2 ∬ |∇|ω||² (w + ...)

with a factor θ_K < 1.  The orthogonal split of the full gradient of
vorticity,

    |∇ω|² = |∇m|² + m²·|∇ê|²          (m = |ω|, ê = ω/|ω|)

makes the radial term |∇m|² and the angular term w = |∇ê|² competitors:
if they are positively correlated, K is hard to absorb; if they are not
positively correlated (c_K ≲ 1), the refined K-absorption conjecture is
supported.

Definitions (given ω = (wx, wy, wz) on a periodic N³ grid, domain
[0, 2π)³, spacing dx = 2π/N; reuses layer4.alignment_bridge's grid
validation, direction field, and spectral-gradient primitives):

  E            = {|ω| ≥ M/4}                      (full box, no ball
                                                     restriction — same
                                                     high-vorticity floor
                                                     as alignment_bridge's
                                                     mask4, simpler and
                                                     sufficient here)
  m            = |ω|                              (pointwise magnitude)
  |∇m|²        = spectral gradient-squared of m
  w            = |∇ê|²                            (alignment_bridge.grad_ehat_sq)
  |∇ω|²        = spectral gradient-squared of the full vector field ω
                 (alignment_bridge.grad_ehat_sq applied directly to
                 (wx, wy, wz) — the function is generic in its three
                 input components, so no reimplementation is needed)

  mean_E(X)    = arithmetic mean of X over grid points in E

  c_K          = mean_E(|∇m|²·w) / (mean_E(|∇m|²) · mean_E(w))
                 c_K = 1  -> decorrelated
                 c_K < 1  -> anti-correlated (favorable for K-absorption)
                 c_K >> 1 -> positively correlated (unfavorable)
                 By convention (documented, not a numerical accident):
                 c_K := 0 whenever mean_E(|∇m|²) = 0 OR mean_E(w) = 0 —
                 in both cases the covariance is degenerate (one factor
                 is identically zero on E), and the "no positive
                 correlation" reading (c_K ≤ 1) is trivially satisfied.

  f_ang        = mean_E(m²·w) / mean_E(|∇ω|²)     equipartition fraction,
                 in [0, 1]: how much of the local enstrophy-gradient
                 budget is carried by the angular (alignment) term vs.
                 the radial (magnitude) term.

  split_residual = |mean_E(|∇m|² + m²·w) - mean_E(|∇ω|²)| / mean_E(|∇ω|²)
                 Sanity check on the orthogonal-split identity above.
                 For analytic fields with a clean magnitude/direction
                 separation (see test_k_correlation.py) this is < 1e-6.
                 For solver snapshots we do NOT hard-assert this < 5%
                 inside this function: |∇m|² is computed as a clean
                 global FFT derivative of the (smooth, untruncated)
                 magnitude field m, but w = |∇ê|² inherits the Gibbs-type
                 ringing documented in alignment_bridge.py's module
                 docstring (ê is forced to zero outside the 1e-12·M
                 floor, so its global FFT derivative rings near that
                 truncation boundary).  E = {|ω| ≥ M/4} sits well inside
                 that floor, which keeps split_residual small in
                 practice, but production callers should still read it
                 as a diagnostic quality indicator, not a proof.

House rules: spectral (FFT) derivatives only, no finite differences; all
functions pure / immutable (return new arrays, never mutate inputs);
every logged result carries claim_id + key_metric + verdict.

Run (as a library — see layer4/alignment_bridge.py for the CLI entry
point that drives the underlying Route 2 3D solver runs):
  python -c "from layer4.k_correlation import k_correlation_diagnostics"
"""

from __future__ import annotations

import os
import sqlite3
import sys
from datetime import datetime

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _ROOT)

from layer4.alignment_bridge import (  # noqa: E402
    _validate_omega_grid,
    compute_direction_field,
    grad_ehat_sq,
)
from layer4.geometric_disorder import spectral_wavenumbers  # noqa: E402

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_KCORR_PASS_MEDIAN_CK = 10.0

# Below this, mean_grad_m2 / mean_w are treated as "exactly zero" for the
# c_K degenerate-covariance convention. Genuine floating-point FFT noise on
# an analytically-constant field (e.g. |grad m|^2 for a pure-twist field)
# lands around 1e-28..1e-32; any physically meaningful gradient on the
# vorticity scales used throughout this program (O(1)-O(1e3)) is many
# orders of magnitude above this floor.
_ZERO_TOL = 1e-12


# =============================================================================
# Core K-correlation primitives
# =============================================================================

def grad_scalar_sq(
    field: np.ndarray, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray
) -> np.ndarray:
    """|∇f|² for a scalar field f, via spectral (FFT) derivatives.

    Same FFT convention as alignment_bridge.grad_ehat_sq. Does not mutate
    field.
    """
    f_hat = np.fft.fftn(field)
    dx_ = np.real(np.fft.ifftn(1j * kx * f_hat))
    dy_ = np.real(np.fft.ifftn(1j * ky * f_hat))
    dz_ = np.real(np.fft.ifftn(1j * kz * f_hat))
    return dx_**2 + dy_**2 + dz_**2


def k_correlation_diagnostics(
    wx: np.ndarray, wy: np.ndarray, wz: np.ndarray, dx: float
) -> dict:
    """Compute the K-correlation diagnostic dict for one vorticity snapshot.

    Parameters
    ----------
    wx, wy, wz : vorticity components, shape (N, N, N), periodic domain [0,2π)^3
    dx         : grid spacing, must equal 2π/N (validated)

    Returns dict with keys: c_K, f_ang, mean_grad_m2, mean_w,
    E_volume_fraction, split_residual.  Raises ValueError on shape
    mismatch, bad dx, or an identically-zero vorticity field (E and M are
    undefined at M=0) — same contract as
    alignment_bridge.bridge_diagnostics.

    Does not mutate wx, wy, wz.
    """
    N = _validate_omega_grid(wx, wy, wz, dx)

    mag = np.sqrt(wx**2 + wy**2 + wz**2)
    M = float(np.max(mag))
    if M < 1e-14:
        raise ValueError(
            "k_correlation_diagnostics: vorticity field is identically zero "
            "(M ≈ 0); E = {|ω| ≥ M/4} is undefined."
        )

    kx, ky, kz = spectral_wavenumbers(N)

    E = mag >= (M / 4.0)
    E_volume_fraction = float(np.count_nonzero(E)) / float(N**3)

    grad_m2 = grad_scalar_sq(mag, kx, ky, kz)  # |∇m|²

    ex, ey, ez, _mag2, _valid = compute_direction_field(wx, wy, wz)
    w = grad_ehat_sq(ex, ey, ez, kx, ky, kz)  # |∇ê|²

    grad_omega2 = grad_ehat_sq(wx, wy, wz, kx, ky, kz)  # |∇ω|² (generic reuse)

    grad_m2_E = grad_m2[E]
    w_E = w[E]
    mag_E = mag[E]
    grad_omega2_E = grad_omega2[E]

    mean_grad_m2 = float(np.mean(grad_m2_E))
    mean_w = float(np.mean(w_E))
    mean_grad_m2_w = float(np.mean(grad_m2_E * w_E))
    mean_m2_w = float(np.mean((mag_E**2) * w_E))
    mean_grad_omega2 = float(np.mean(grad_omega2_E))

    if mean_grad_m2 <= _ZERO_TOL or mean_w <= _ZERO_TOL:
        c_K = 0.0
    else:
        c_K = mean_grad_m2_w / (mean_grad_m2 * mean_w)

    if mean_grad_omega2 > _ZERO_TOL:
        f_ang = float(np.clip(mean_m2_w / mean_grad_omega2, 0.0, 1.0))
        split_residual = float(
            abs((mean_grad_m2 + mean_m2_w) - mean_grad_omega2) / mean_grad_omega2
        )
    else:
        f_ang = 0.0
        split_residual = 0.0

    return {
        "c_K": float(c_K),
        "f_ang": f_ang,
        "mean_grad_m2": mean_grad_m2,
        "mean_w": mean_w,
        "E_volume_fraction": E_volume_fraction,
        "split_residual": split_residual,
    }


# =============================================================================
# Experiment-level summarization (reduces a run_bridge_experiment
# timeseries — which already carries per-step c_K/f_ang/split_residual via
# alignment_bridge.run_bridge_experiment's K-correlation wiring — to a
# PASS/FAIL summary, without re-running the solver)
# =============================================================================

def kcorr_summary_from_result(
    result: dict,
    exp_id_override: str | None = None,
    claim_id_override: str | None = None,
) -> dict:
    """Reduce a alignment_bridge.run_bridge_experiment() result dict to a
    K-correlation summary row + PASS/FAIL verdict.

    Verdict PASS iff every recorded c_K and f_ang is finite AND
    median(c_K) <= 10 — i.e. the K-radial/angular coupling stays O(1) or
    below over the run. The scientific c_K value is reported regardless
    of verdict; the threshold only gates the PASS/FAIL label.
    """
    timeseries = result.get("timeseries", [])
    c_K_vals = np.array([rec["c_K"] for rec in timeseries], dtype=float)
    f_ang_vals = np.array([rec["f_ang"] for rec in timeseries], dtype=float)
    split_residual_vals = np.array(
        [rec["split_residual"] for rec in timeseries], dtype=float
    )

    all_finite = bool(np.all(np.isfinite(c_K_vals)) and np.all(np.isfinite(f_ang_vals)))
    c_K_max = float(np.max(c_K_vals)) if c_K_vals.size else float("nan")
    c_K_median = float(np.median(c_K_vals)) if c_K_vals.size else float("nan")
    c_K_final = float(c_K_vals[-1]) if c_K_vals.size else float("nan")
    f_ang_median = float(np.median(f_ang_vals)) if f_ang_vals.size else float("nan")
    split_residual_min = float(np.min(split_residual_vals)) if split_residual_vals.size else float("nan")
    split_residual_max = float(np.max(split_residual_vals)) if split_residual_vals.size else float("nan")

    verdict = "PASS" if (
        all_finite
        and np.isfinite(c_K_median)
        and c_K_median <= _KCORR_PASS_MEDIAN_CK
    ) else "FAIL"

    ic = result.get("ic", "")
    tag = {"tg": "TG", "shear": "SH", "adv": "ADV"}.get(ic, ic.upper() if ic else "")
    exp_id = exp_id_override or f"EXP-L4-KCORR-{tag}-001"
    claim_id = claim_id_override or f"k-absorption-corr-{ic}"

    return {
        "exp_id": exp_id,
        "claim_id": claim_id,
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic,
        "N": result.get("N"),
        "n_steps": result.get("n_steps"),
        "seed": result.get("seed"),
        "nu": result.get("nu"),
        "c_K_max": c_K_max,
        "c_K_median": c_K_median,
        "c_K_final": c_K_final,
        "f_ang_median": f_ang_median,
        "split_residual_min": split_residual_min,
        "split_residual_max": split_residual_max,
        "verdict": verdict,
        "wall_time_s": result.get("wall_time_s"),
        "key_metric": f"median_c_K={c_K_median:.4f}",
        "timeseries": timeseries,
    }


# =============================================================================
# Results logging (layer4-module table pattern; see alignment_bridge.py)
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_kcorr_experiments (
            id                   INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id               TEXT UNIQUE NOT NULL,
            claim_id             TEXT NOT NULL,
            timestamp            TEXT NOT NULL,
            ic                   TEXT,
            N                    INTEGER,
            n_steps              INTEGER,
            seed                 INTEGER,
            nu                   REAL,
            c_K_max              REAL,
            c_K_median           REAL,
            c_K_final            REAL,
            f_ang_median         REAL,
            split_residual_min   REAL,
            split_residual_max   REAL,
            verdict              TEXT NOT NULL,
            wall_time_s          REAL,
            key_metric           TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_kcorr_timeseries (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id           TEXT NOT NULL,
            step             INTEGER,
            t                REAL,
            c_K              REAL,
            f_ang            REAL,
            split_residual   REAL,
            E_volume_fraction REAL
        )
    """)
    conn.commit()
    return conn


def log_kcorr_result(db_path: str, summary: dict) -> None:
    """Append one K-correlation experiment (+ its timeseries) to results.db.

    Follows the layer4 per-module table convention (see
    alignment_bridge.log_result): an UPSERT of the summary row keyed on
    exp_id, plus append-only timeseries rows.
    """
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "ic", "N", "n_steps", "seed",
        "nu", "c_K_max", "c_K_median", "c_K_final", "f_ang_median",
        "split_residual_min", "split_residual_max",
        "verdict", "wall_time_s", "key_metric",
    ]
    values = tuple(summary.get(c) for c in cols)
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO layer4_kcorr_experiments ({', '.join(cols)}) "
        f"VALUES ({ph})",
        values,
    )
    for rec in summary.get("timeseries", []):
        conn.execute(
            """
            INSERT INTO layer4_kcorr_timeseries
            (exp_id, step, t, c_K, f_ang, split_residual, E_volume_fraction)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                summary["exp_id"], rec["step"], rec["t"], rec["c_K"],
                rec["f_ang"], rec["split_residual"], rec["E_volume_fraction"],
            ),
        )
    conn.commit()
    conn.close()


def run_and_log_kcorr_experiment(
    ic: str = "tg",
    N: int = 64,
    n_steps: int = 400,
    nu: float | None = None,
    record_every: int = 10,
    seed: int = 42,
    verbose: bool = False,
    db_path: str | None = _DEFAULT_DB,
    exp_id_override: str | None = None,
    claim_id_override: str | None = None,
) -> dict:
    """Run the Route 2 3D solver via alignment_bridge.run_bridge_experiment
    (which now also records c_K/f_ang/split_residual on every snapshot),
    summarize the K-correlation diagnostic, and log it to
    layer4_kcorr_experiments / layer4_kcorr_timeseries.

    Does not log to layer4_bridge_experiments / layer4_bridge_timeseries —
    callers who want both (e.g. the escalated adversarial run) should call
    alignment_bridge.run_bridge_experiment(..., db_path=...) themselves
    first and pass its result dict into kcorr_summary_from_result +
    log_kcorr_result directly.
    """
    from layer4.alignment_bridge import run_bridge_experiment

    result = run_bridge_experiment(
        ic=ic, N=N, n_steps=n_steps, nu=nu, record_every=record_every,
        seed=seed, verbose=verbose, db_path=None,
    )
    summary = kcorr_summary_from_result(
        result, exp_id_override=exp_id_override, claim_id_override=claim_id_override
    )
    if db_path is not None:
        log_kcorr_result(db_path, summary)
    return summary
