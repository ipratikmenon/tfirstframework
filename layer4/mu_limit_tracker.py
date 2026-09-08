"""
mu_limit_tracker.py — μ(T)→ν Limit Tracker (M10)
==================================================
Standalone tracker for the Phase 4 Prize limit:
  ε(δ) = |A_min(δ) − A_ref| → 0  as  δ → 0

Consolidates M4 (2D) and M8 (3D) results into a single analysis tool.
Works on any sequence of (δ, A_min) pairs — from mu_limit_2D.py,
mu_sweep_3D.py, or synthetic data.

Functions
---------
compute_eps(delta_values, A_min_values, A_ref)
    ε(δ) = |A_min(δ) − A_ref| for each δ.

fit_scaling_exponent(delta_values, eps_values, min_points)
    Log-log fit: ε ~ C · δ^α.  Returns α, C, R².

mu_limit_verdict(alpha, eps_min, alpha_tolerance)
    PASS if α ∈ (1−tol, 1+tol) (linear scaling → regular limit).

run_mu_limit_analysis(delta_values, A_min_values, A_ref)
    Full analysis dict with scaling fit + verdict.

compare_2d_3d(results_2d, results_3d)
    Side-by-side comparison of M4 (2D) and M8 (3D) scaling.
"""

from __future__ import annotations

import numpy as np
from typing import Sequence


# ─────────────────────────────────────────────────────────────────────────────
# Core computations
# ─────────────────────────────────────────────────────────────────────────────

def compute_eps(
    delta_values: Sequence[float],
    A_min_values: Sequence[float],
    A_ref: float,
) -> np.ndarray:
    """ε(δ) = |A_min(δ) − A_ref| for each δ."""
    A_min = np.asarray(A_min_values, dtype=float)
    return np.abs(A_min - A_ref)


def fit_scaling_exponent(
    delta_values: Sequence[float],
    eps_values: Sequence[float],
    min_points: int = 3,
) -> dict:
    """Log-log fit: ε ~ C · δ^α.

    Excludes zero ε values (can't take log of 0).

    Returns
    -------
    dict with alpha, C, R2, fit_ok, n_points.
    """
    delta = np.asarray(delta_values, dtype=float)
    eps   = np.asarray(eps_values,   dtype=float)

    valid = (delta > 0.0) & (eps > 0.0)
    delta = delta[valid]
    eps   = eps[valid]

    if len(eps) < min_points:
        return {
            "alpha": float("nan"), "C": float("nan"),
            "R2": float("nan"), "fit_ok": False, "n_points": int(valid.sum()),
        }

    log_d = np.log(delta)
    log_e = np.log(eps)

    A_mat = np.column_stack([np.ones_like(log_d), log_d])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, log_e, rcond=None)
    except np.linalg.LinAlgError:
        return {
            "alpha": float("nan"), "C": float("nan"),
            "R2": float("nan"), "fit_ok": False, "n_points": len(eps),
        }

    log_C, alpha = coeffs
    C = float(np.exp(log_C))

    log_e_pred = A_mat @ coeffs
    ss_res = float(np.sum((log_e - log_e_pred)**2))
    ss_tot = float(np.sum((log_e - np.mean(log_e))**2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {
        "alpha":    float(alpha),
        "C":        C,
        "R2":       float(R2) if np.isfinite(R2) else float("nan"),
        "fit_ok":   True,
        "n_points": len(eps),
    }


def mu_limit_verdict(
    alpha: float,
    eps_values: Sequence[float],
    alpha_tolerance: float = 0.3,
) -> str:
    """PASS if α ∈ (1−tol, 1+tol): linear ε(δ)→0 as δ→0 (regular Prize limit).

    Also requires max(ε) > 0 (so there's something to fit) and
    all ε ≥ 0 (A_min ≤ A_ref or A_min ≥ A_ref — monotone convergence not required,
    only magnitude → 0 linearly).
    """
    eps = np.asarray(eps_values, dtype=float)
    if np.all(eps == 0.0):
        return "PASS"   # perfect convergence
    if not np.isfinite(alpha):
        return "INSUFFICIENT_DATA"
    if abs(alpha - 1.0) <= alpha_tolerance:
        return "PASS"
    return "FAIL"


# ─────────────────────────────────────────────────────────────────────────────
# Full analysis
# ─────────────────────────────────────────────────────────────────────────────

def run_mu_limit_analysis(
    delta_values: Sequence[float],
    A_min_values: Sequence[float],
    A_ref: float,
    alpha_tolerance: float = 0.3,
    label: str = "",
) -> dict:
    """Full μ(T)→ν limit analysis.

    Parameters
    ----------
    delta_values  : perturbation amplitudes (decreasing sequence)
    A_min_values  : A_min_global at each δ
    A_ref         : A(T̄) — reference value for uniform T
    alpha_tolerance: tolerance on linear scaling exponent
    label         : optional label ('2D', '3D', etc.)

    Returns
    -------
    dict with:
      delta_values, A_min_values, A_ref, eps_values,
      scaling_fit (dict from fit_scaling_exponent),
      verdict, label,
      A_min_plateau (mean of last 3 values),
      A_min_plateau_pct (deviation from A_ref),
      eps_range (min, max).
    """
    delta = np.asarray(delta_values, dtype=float)
    A_min = np.asarray(A_min_values, dtype=float)
    eps   = compute_eps(delta, A_min, A_ref)

    fit = fit_scaling_exponent(delta, eps)
    verdict = mu_limit_verdict(fit["alpha"], eps, alpha_tolerance)

    # Plateau: last 3 δ values (smallest δ → closest to ν limit)
    n_plateau = min(3, len(A_min))
    plateau = float(np.mean(A_min[-n_plateau:])) if n_plateau > 0 else float("nan")
    plateau_pct = abs(plateau - A_ref) / A_ref * 100.0 if A_ref > 0 else float("nan")

    return {
        "label":           label,
        "delta_values":    delta.tolist(),
        "A_min_values":    A_min.tolist(),
        "A_ref":           A_ref,
        "eps_values":      eps.tolist(),
        "scaling_fit":     fit,
        "verdict":         verdict,
        "A_min_plateau":   plateau,
        "A_min_plateau_pct": plateau_pct,
        "eps_range":       (float(np.min(eps)), float(np.max(eps))),
        "alpha":           fit["alpha"],
        "R2":              fit["R2"],
    }


# ─────────────────────────────────────────────────────────────────────────────
# Comparison: M4 (2D) vs M8 (3D)
# ─────────────────────────────────────────────────────────────────────────────

def compare_2d_3d(
    results_2d: dict,
    results_3d: dict,
) -> dict:
    """Side-by-side comparison of M4 (2D) and M8 (3D) mu-limit analyses.

    Both inputs are dicts from run_mu_limit_analysis.

    Returns consistency verdict:
      CONSISTENT if |α_2D − α_3D| < 0.2 and both PASS.
    """
    alpha_2d = results_2d.get("alpha", float("nan"))
    alpha_3d = results_3d.get("alpha", float("nan"))

    both_pass = (results_2d.get("verdict") == "PASS" and
                 results_3d.get("verdict") == "PASS")

    if np.isfinite(alpha_2d) and np.isfinite(alpha_3d):
        alpha_diff = abs(alpha_2d - alpha_3d)
        consistent = both_pass and alpha_diff < 0.2
    else:
        alpha_diff = float("nan")
        consistent = False

    verdict = "CONSISTENT" if consistent else ("PASS (one dim)" if both_pass else "INCONSISTENT")

    return {
        "verdict":      verdict,
        "consistent":   consistent,
        "alpha_2d":     alpha_2d,
        "alpha_3d":     alpha_3d,
        "alpha_diff":   alpha_diff,
        "A_ref_2d":     results_2d.get("A_ref"),
        "A_ref_3d":     results_3d.get("A_ref"),
        "plateau_2d":   results_2d.get("A_min_plateau"),
        "plateau_3d":   results_3d.get("A_min_plateau"),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: load from mu_sweep results
# ─────────────────────────────────────────────────────────────────────────────

def load_from_sweep_results(sweep_results: dict, A_ref: float, label: str = "") -> dict:
    """Build analysis from the dict returned by mu_sweep_3D.run_mu_sweep_3D.

    sweep_results: dict mapping δ → {A_min_global, ...}
    """
    deltas = sorted(sweep_results.keys(), reverse=True)
    A_mins = [sweep_results[d]["A_min_global"] for d in deltas]
    return run_mu_limit_analysis(deltas, A_mins, A_ref, label=label)
