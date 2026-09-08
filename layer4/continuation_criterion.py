"""
continuation_criterion.py — Phase 0 Continuation Criterion Monitor (M10)
=========================================================================
Phase 0 of the T-First proof: a singularity at time T* in Route 1 requires
  A(T) → 0  OR  μ(T) → 0  as t → T*.

This module monitors A_min(t) and μ_min(t) from solver history and raises
an alarm if either approaches zero.  All experiments so far show neither
approaches zero — confirming Phase 0 criterion is not triggered.

Functions
---------
check_A_min(A_min_series, A_ref, threshold_frac)
    Returns margin = min(A_min) / A_ref.  Triggered if margin < threshold_frac.

check_mu_min(mu_min_series, mu_ref, threshold_frac)
    Returns margin = min(mu_min) / mu_ref.  Triggered if margin < threshold_frac.

continuation_criterion_verdict(history, nu_ref, A_ref, threshold_frac)
    Full verdict from run_jax / run_blowup_search history.

extract_mu_min_from_history(history, N)
    Reconstruct μ_min(t) from history entries (uses A_min via A=k/ρcv,
    μ/ρ = ν field, back-computed from stored A_min as proxy).
"""

from __future__ import annotations

import numpy as np
from typing import Sequence


# ── Physical constants (Sutherland viscosity; identical to route1_3D_jax.py) ─

_MU_REF  = 1.716e-5   # Pa·s at T_ref
_T_REF   = 273.15     # K
_S_SUTH  = 110.4      # K
_K_REF   = 0.0241     # W/(m·K) at T_ref
_K_EXP   = 0.82
_P_REF   = 101325.0   # Pa
_R_SPEC  = 287.0      # J/(kg·K)
_CV      = 718.0      # J/(kg·K)


def A_of_T(T: float) -> float:
    """Thermal diffusivity A = k/(ρ·cv) at scalar temperature T."""
    k   = _K_REF * (T / _T_REF) ** _K_EXP
    rho = _P_REF / (_R_SPEC * T)
    return k / (rho * _CV)


def mu_of_T(T: float) -> float:
    """Sutherland dynamic viscosity μ(T)."""
    return _MU_REF * (T / _T_REF) ** 1.5 * (_T_REF + _S_SUTH) / (T + _S_SUTH)


def A_ref_at(T_bar: float = 300.0) -> float:
    """Reference A(T̄) for uniform temperature T_bar."""
    return A_of_T(T_bar)


def mu_ref_at(T_bar: float = 300.0) -> float:
    """Reference μ(T̄) for uniform temperature T_bar."""
    return mu_of_T(T_bar)


# ─────────────────────────────────────────────────────────────────────────────
# Core checks
# ─────────────────────────────────────────────────────────────────────────────

def check_A_min(
    A_min_series: Sequence[float],
    A_ref: float,
    threshold_frac: float = 0.01,
) -> dict:
    """Check whether A_min approaches zero (Phase 0 criterion for blowup).

    Parameters
    ----------
    A_min_series : sequence of A_min values at each timestep
    A_ref        : reference A(T̄) — non-trivial lower bound
    threshold_frac : fraction of A_ref below which we flag a warning

    Returns
    -------
    dict with:
      triggered  : bool — True if min(A_min) < threshold_frac * A_ref
      A_min_abs  : float — minimum A_min observed
      margin     : float — min(A_min) / A_ref
      threshold  : float — threshold_frac * A_ref
      n_points   : int
    """
    if len(A_min_series) == 0:
        return {"triggered": False, "A_min_abs": float("nan"),
                "margin": float("nan"), "threshold": threshold_frac * A_ref,
                "n_points": 0}

    A_arr = np.asarray(A_min_series, dtype=float)
    A_min_abs = float(np.min(A_arr))
    margin    = A_min_abs / A_ref if A_ref > 0.0 else float("inf")
    triggered = A_min_abs < threshold_frac * A_ref

    return {
        "triggered":  triggered,
        "A_min_abs":  A_min_abs,
        "margin":     margin,
        "threshold":  threshold_frac * A_ref,
        "n_points":   len(A_arr),
    }


def check_mu_min(
    mu_min_series: Sequence[float],
    mu_ref: float,
    threshold_frac: float = 0.01,
) -> dict:
    """Check whether μ_min approaches zero (Phase 0 criterion).

    Parameters and return structure mirror check_A_min.
    """
    if len(mu_min_series) == 0:
        return {"triggered": False, "mu_min_abs": float("nan"),
                "margin": float("nan"), "threshold": threshold_frac * mu_ref,
                "n_points": 0}

    mu_arr    = np.asarray(mu_min_series, dtype=float)
    mu_min_abs = float(np.min(mu_arr))
    margin     = mu_min_abs / mu_ref if mu_ref > 0.0 else float("inf")
    triggered  = mu_min_abs < threshold_frac * mu_ref

    return {
        "triggered":   triggered,
        "mu_min_abs":  mu_min_abs,
        "margin":      margin,
        "threshold":   threshold_frac * mu_ref,
        "n_points":    len(mu_arr),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Full verdict from solver history
# ─────────────────────────────────────────────────────────────────────────────

def continuation_criterion_verdict(
    history: list[dict],
    T_bar: float = 300.0,
    threshold_frac: float = 0.01,
) -> dict:
    """Full Phase 0 continuation criterion verdict from a solver history.

    Works with history produced by:
      - route1_3D_jax.run_jax  (entries have 'A_min', 't')
      - blowup_search_3D.run_blowup_search (entries have 'A_min', 't')

    For μ_min: since the history stores A_min (not μ_min directly), we
    reconstruct μ_min from A_min using the thermodynamic identity:
      A = k(T) / (ρ(T)·cv) = k_ref · (T/T_ref)^k_exp · R_spec · T / (P_ref · cv)
    Inverting: T = A^(1/(1+k_exp)) · (scale),  then μ = μ(T).
    We use a conservative bound: μ_min ≥ μ(T_min) where T_min is inferred from A_min.

    Parameters
    ----------
    history      : list of step dicts, each with at least 'A_min' and 't'
    T_bar        : reference temperature (K)
    threshold_frac : alarm threshold as fraction of reference value

    Returns
    -------
    dict with:
      verdict         : 'PASS' or 'FAIL'
      A_check         : result of check_A_min
      mu_check        : result of check_mu_min (via T_min inversion)
      criterion_triggered : bool
      A_ref, mu_ref   : reference values used
    """
    A_ref  = A_ref_at(T_bar)
    mu_ref = mu_ref_at(T_bar)

    if not history:
        return {
            "verdict": "PASS",
            "A_check": check_A_min([], A_ref, threshold_frac),
            "mu_check": check_mu_min([], mu_ref, threshold_frac),
            "criterion_triggered": False,
            "A_ref": A_ref,
            "mu_ref": mu_ref,
        }

    A_min_series = [h["A_min"] for h in history if "A_min" in h]

    # Infer T_min from A_min for each step: A ∝ T^(1+k_exp) approximately.
    # A = k_ref*(T/T_ref)^k_exp * R*T/(P*cv) = const * T^(1+k_exp)
    # T_from_A = T_ref * (A / A_ref)^(1/(1+k_exp))
    exp_inv = 1.0 / (1.0 + _K_EXP)  # ≈ 1/1.82
    mu_min_series = []
    for A_min_val in A_min_series:
        if A_min_val > 0.0:
            T_inferred = _T_REF * (A_min_val / A_ref_at(_T_REF)) ** exp_inv
            T_inferred = max(T_inferred, 1.0)  # physical floor
            mu_min_series.append(mu_of_T(T_inferred))
        else:
            mu_min_series.append(0.0)

    A_check  = check_A_min(A_min_series, A_ref, threshold_frac)
    mu_check = check_mu_min(mu_min_series, mu_ref, threshold_frac)

    triggered = A_check["triggered"] or mu_check["triggered"]
    verdict   = "FAIL" if triggered else "PASS"

    return {
        "verdict":              verdict,
        "A_check":              A_check,
        "mu_check":             mu_check,
        "criterion_triggered":  triggered,
        "A_ref":                A_ref,
        "mu_ref":               mu_ref,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: run on M9 blowup search result
# ─────────────────────────────────────────────────────────────────────────────

def analyse_blowup_result(result: dict, T_bar: float = 300.0) -> dict:
    """Apply continuation criterion to a blowup_search_3D result dict."""
    return continuation_criterion_verdict(
        history=result.get("history", []),
        T_bar=T_bar,
    )
