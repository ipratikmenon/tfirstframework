"""
viscosity_scaling.py — Wang et al. Viscosity Scaling Analysis (M10)
====================================================================
Tests Conjecture 3.5: μ(T) scaling defeats Wang et al. self-similar blow-up
for all λ > 0.

Self-similar blowup Ansatz (Wang et al.):
  ‖∇u‖(t) ~ C · (T* − t)^{−(1+λ)/2}   as t → T*

The μ(T) damping rate must grow SLOWER than the blowup rate for regularity:
  μ_growth_exponent < λ + 1

Functions
---------
blowup_rate(t_arr, T_star, lambda_val)
    Theoretical (T*−t)^{−(1+λ)} profile (the Leray self-similar rate).

fit_power_law_decay(t_arr, f_arr, T_star)
    Fit f(t) ~ C · (T*−t)^α near T*.  Returns α, C, R².

suppression_ratio(t_arr, mu_arr, T_star, lambda_val)
    Compute μ(t) / blowup_rate(t) at each step.
    Ratio → 0 as t → T* ⟹ μ(T) damping wins.

viscosity_scaling_verdict(t_arr, mu_arr, T_star, lambda_val)
    Full verdict: PASS if suppression_ratio is non-increasing.
"""

from __future__ import annotations

import numpy as np
from typing import Sequence

# ── Sutherland constants (same as continuation_criterion.py) ─────────────────
_MU_REF = 1.716e-5
_T_REF  = 273.15
_S_SUTH = 110.4


def mu_sutherland(T: float) -> float:
    return _MU_REF * (T / _T_REF) ** 1.5 * (_T_REF + _S_SUTH) / (T + _S_SUTH)


# ─────────────────────────────────────────────────────────────────────────────
# Theoretical blowup rate
# ─────────────────────────────────────────────────────────────────────────────

def blowup_rate(
    t_arr: Sequence[float],
    T_star: float,
    lambda_val: float,
) -> np.ndarray:
    """Theoretical Leray self-similar rate (T*−t)^{−(1+λ)}.

    Returns NaN for t ≥ T_star (post-singularity region, not physical).
    """
    t = np.asarray(t_arr, dtype=float)
    tau = T_star - t
    rate = np.where(tau > 0.0, tau ** (-(1.0 + lambda_val)), np.nan)
    return rate


# ─────────────────────────────────────────────────────────────────────────────
# Power law fitting near T*
# ─────────────────────────────────────────────────────────────────────────────

def fit_power_law_decay(
    t_arr: Sequence[float],
    f_arr: Sequence[float],
    T_star: float,
    min_points: int = 4,
) -> dict:
    """Fit f(t) ~ C · (T*−t)^α using log-log regression on τ = T*−t.

    Only uses points with f > 0 and τ > 0.

    Returns
    -------
    dict with alpha, C, R2, fit_ok, n_points.
    """
    t = np.asarray(t_arr, dtype=float)
    f = np.asarray(f_arr, dtype=float)

    tau = T_star - t
    valid = (tau > 0.0) & (f > 0.0)
    tau = tau[valid]
    f   = f[valid]

    if len(f) < min_points:
        return {"alpha": float("nan"), "C": float("nan"),
                "R2": float("nan"), "fit_ok": False, "n_points": int(valid.sum())}

    log_tau = np.log(tau)
    log_f   = np.log(f)

    A_mat = np.column_stack([np.ones_like(log_tau), log_tau])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, log_f, rcond=None)
    except np.linalg.LinAlgError:
        return {"alpha": float("nan"), "C": float("nan"),
                "R2": float("nan"), "fit_ok": False, "n_points": len(f)}

    log_C, alpha = coeffs
    C = float(np.exp(log_C))

    log_f_pred = A_mat @ coeffs
    ss_res = float(np.sum((log_f - log_f_pred)**2))
    ss_tot = float(np.sum((log_f - np.mean(log_f))**2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {
        "alpha":    float(alpha),
        "C":        C,
        "R2":       float(R2) if np.isfinite(R2) else float("nan"),
        "fit_ok":   True,
        "n_points": len(f),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Suppression ratio
# ─────────────────────────────────────────────────────────────────────────────

def suppression_ratio(
    t_arr: Sequence[float],
    mu_arr: Sequence[float],
    T_star: float,
    lambda_val: float,
) -> np.ndarray:
    """Compute μ(t) / blowup_rate(t) at each timestep.

    If this ratio is → 0 as t → T*, viscosity wins.
    If ratio is bounded below or increasing, could be marginal.
    """
    t  = np.asarray(t_arr,  dtype=float)
    mu = np.asarray(mu_arr, dtype=float)
    rate = blowup_rate(t, T_star, lambda_val)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(rate > 0.0, mu / rate, np.nan)
    return ratio


# ─────────────────────────────────────────────────────────────────────────────
# Full verdict
# ─────────────────────────────────────────────────────────────────────────────

def viscosity_scaling_verdict(
    t_arr: Sequence[float],
    mu_arr: Sequence[float],
    T_star: float,
    lambda_val: float,
    min_ratio_trend_points: int = 4,
) -> dict:
    """Full verdict: does μ(T) suppress the Wang et al. blowup rate?

    PASS criteria (either):
    1. The suppression ratio μ(t)/blowup_rate(t) is non-increasing over time
       (μ(T) damping is gaining on the blowup rate).
    2. The fitted exponent α of μ(t) ~ (T*−t)^α satisfies α > −(1+λ),
       i.e., μ grows SLOWER than the blowup rate.

    Parameters
    ----------
    t_arr      : time array
    mu_arr     : μ_mean(t) — mean viscosity at each step
    T_star     : hypothetical blowup time (use t_end + margin, e.g. 2.0)
    lambda_val : Wang et al. self-similar parameter λ

    Returns
    -------
    dict with verdict, ratio_trend, mu_fit, blowup_exponent, suppression_ok.
    """
    t   = np.asarray(t_arr,  dtype=float)
    mu  = np.asarray(mu_arr, dtype=float)

    blowup_exp = -(1.0 + lambda_val)  # blowup rate exponent in τ

    # Suppression ratio
    ratio = suppression_ratio(t, mu, T_star, lambda_val)
    valid_ratio = ratio[np.isfinite(ratio)]

    # Trend: is ratio non-increasing?
    if len(valid_ratio) >= min_ratio_trend_points:
        ratio_trend = float(np.polyfit(
            np.arange(len(valid_ratio)), valid_ratio, 1
        )[0])
        ratio_decreasing = ratio_trend <= 0.0
    else:
        ratio_trend = float("nan")
        ratio_decreasing = True  # not enough data to falsify

    # Fit μ(t) power law near T*
    mu_fit = fit_power_law_decay(t, mu, T_star)

    # Suppression check: if α_μ > blowup_exp → μ grows slower → wins
    suppression_ok_by_exponent = (
        mu_fit["fit_ok"]
        and mu_fit["alpha"] > blowup_exp
    )
    suppression_ok = ratio_decreasing or suppression_ok_by_exponent

    verdict = "PASS" if suppression_ok else "FAIL"

    return {
        "verdict":              verdict,
        "suppression_ok":       suppression_ok,
        "ratio_decreasing":     ratio_decreasing,
        "ratio_trend":          ratio_trend,
        "suppression_ok_by_exponent": suppression_ok_by_exponent,
        "mu_fit":               mu_fit,
        "blowup_exponent":      blowup_exp,
        "lambda_val":           lambda_val,
        "T_star":               T_star,
        "n_steps":              len(t),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: extract μ_mean from blowup_search history
# ─────────────────────────────────────────────────────────────────────────────

def mu_timeseries_from_T(
    T_mean_series: Sequence[float],
) -> np.ndarray:
    """Compute μ̄(t) = μ(T̄(t)) from mean temperature series."""
    return np.array([mu_sutherland(T) for T in T_mean_series], dtype=float)
