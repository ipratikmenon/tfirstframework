"""
gronwall_fitter.py — Grönwall Energy Decay Fitter (M10)
========================================================
Fits the H¹ energy E1(t) = ½‖∇u‖²_{L²} to an exponential decay model

    E1(t) = E1(0) · exp(−c·t)

using a least-squares log-linear fit.  The damping constant c > 0 is
consistent with Grönwall-type regularity: dE1/dt ≤ −c·E1 implies
E1 stays bounded and decays — no blowup.

Functions
---------
extract_E1_jax(u, v, w, kx, ky, kz)
    Compute E1 = ½‖∇u‖²_{L²} from JAX arrays (spectral).

fit_gronwall(E1_arr, t_arr)
    Fit E1(t) = E1(0)·exp(−c·t).  Returns c, R², verdict.

fit_energy_timeseries(history, E1_key='E')
    Convenience: extract energy from solver history and fit.

gronwall_verdict(c, R2, c_threshold, R2_threshold)
    PASS if c > c_threshold (positive damping) and R² > R2_threshold.
"""

from __future__ import annotations

import numpy as np
from typing import Sequence


# ─────────────────────────────────────────────────────────────────────────────
# H¹ energy extraction (spectral, works on JAX or NumPy arrays)
# ─────────────────────────────────────────────────────────────────────────────

def extract_E1_spectral(
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    kz: np.ndarray,
) -> float:
    """E1 = ½‖∇u‖²_{L²} computed spectrally.

    ‖∇u‖² = ∑_{i,j} ‖∂_j u_i‖² = Parseval: ∑_k |k|²·(|û|² + |v̂|² + |ŵ|²).
    """
    u_np = np.asarray(u, dtype=np.float32)
    v_np = np.asarray(v, dtype=np.float32)
    w_np = np.asarray(w, dtype=np.float32)
    kx_np = np.asarray(kx, dtype=np.float32)
    ky_np = np.asarray(ky, dtype=np.float32)
    kz_np = np.asarray(kz, dtype=np.float32)

    k2 = kx_np**2 + ky_np**2 + kz_np**2
    N3 = u_np.size

    u_hat = np.fft.fftn(u_np)
    v_hat = np.fft.fftn(v_np)
    w_hat = np.fft.fftn(w_np)

    power = (np.abs(u_hat)**2 + np.abs(v_hat)**2 + np.abs(w_hat)**2) / N3**2
    E1 = float(0.5 * np.sum(k2 * power))
    return E1


# ─────────────────────────────────────────────────────────────────────────────
# Grönwall fit
# ─────────────────────────────────────────────────────────────────────────────

def fit_gronwall(
    E1_arr: Sequence[float],
    t_arr: Sequence[float],
    min_points: int = 5,
) -> dict:
    """Fit E1(t) = E1(0) · exp(−c·t) by log-linear least squares.

    Parameters
    ----------
    E1_arr    : sequence of E1 values
    t_arr     : corresponding time values
    min_points: minimum data points required for a fit

    Returns
    -------
    dict with:
      c          : float — damping constant (> 0 means decay)
      c_se       : float — standard error of c
      E1_0_fit   : float — fitted initial value
      R2         : float — coefficient of determination
      fit_ok     : bool — fit converged with enough points
      n_points   : int
      verdict    : 'PASS' / 'FAIL' / 'INSUFFICIENT_DATA'
    """
    E1 = np.asarray(E1_arr, dtype=float)
    t  = np.asarray(t_arr,  dtype=float)

    positive = (E1 > 0.0)
    E1 = E1[positive]
    t  = t[positive]

    if len(E1) < min_points:
        return {
            "c": float("nan"), "c_se": float("nan"),
            "E1_0_fit": float("nan"), "R2": float("nan"),
            "fit_ok": False, "n_points": int(positive.sum()),
            "verdict": "INSUFFICIENT_DATA",
        }

    # Log-linear: log(E1) = log(E1_0) − c·t
    log_E1 = np.log(E1)

    # Weighted least squares: [1, -t] @ [log_E1_0, c] = log_E1
    A_mat = np.column_stack([np.ones_like(t), t])
    try:
        coeffs, residuals, rank, sv = np.linalg.lstsq(A_mat, log_E1, rcond=None)
    except np.linalg.LinAlgError:
        return {
            "c": float("nan"), "c_se": float("nan"),
            "E1_0_fit": float("nan"), "R2": float("nan"),
            "fit_ok": False, "n_points": len(E1),
            "verdict": "FIT_ERROR",
        }

    log_E1_0_fit, neg_c = coeffs
    c = -neg_c
    E1_0_fit = float(np.exp(log_E1_0_fit))

    # R²
    log_E1_pred = A_mat @ coeffs
    ss_res = float(np.sum((log_E1 - log_E1_pred)**2))
    ss_tot = float(np.sum((log_E1 - np.mean(log_E1))**2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # Standard error of c (from diagonal of covariance matrix)
    n = len(E1)
    p = 2
    if n > p and ss_res > 0:
        sigma2 = ss_res / (n - p)
        try:
            cov = sigma2 * np.linalg.inv(A_mat.T @ A_mat)
            c_se = float(np.sqrt(cov[1, 1]))
        except np.linalg.LinAlgError:
            c_se = float("nan")
    else:
        c_se = 0.0

    fit_ok = np.isfinite(c) and np.isfinite(R2)

    return {
        "c":         float(c),
        "c_se":      c_se,
        "E1_0_fit":  E1_0_fit,
        "R2":        float(R2) if np.isfinite(R2) else float("nan"),
        "fit_ok":    fit_ok,
        "n_points":  len(E1),
        "verdict":   gronwall_verdict(c, R2) if fit_ok else "INSUFFICIENT_DATA",
    }


def gronwall_verdict(
    c: float,
    R2: float,
    c_threshold: float = 0.0,
    R2_threshold: float = 0.5,
) -> str:
    """PASS if c > c_threshold (positive damping) and R² > R2_threshold."""
    if not (np.isfinite(c) and np.isfinite(R2)):
        return "INSUFFICIENT_DATA"
    if c > c_threshold and R2 > R2_threshold:
        return "PASS"
    return "FAIL"


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: fit from solver history
# ─────────────────────────────────────────────────────────────────────────────

def fit_energy_timeseries(
    history: list[dict],
    E_key: str = "E",
    t_key: str = "t",
) -> dict:
    """Extract E(t) from solver history and fit Grönwall decay.

    Works with history produced by run_jax (entries have 'E' and 't').
    Note: the history 'E' = ½‖u‖²_{L²} (L² energy), not H¹.
    Use extract_E1_spectral for the H¹ norm if u, v, w arrays are available.

    Parameters
    ----------
    history : list of step dicts
    E_key   : key for energy in history entries
    t_key   : key for time in history entries

    Returns
    -------
    Grönwall fit dict (same as fit_gronwall).
    """
    E1_vals = [h[E_key] for h in history if E_key in h and t_key in h]
    t_vals  = [h[t_key]  for h in history if E_key in h and t_key in h]
    return fit_gronwall(E1_vals, t_vals)
