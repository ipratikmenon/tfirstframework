"""
Mean Morrey Diagnostic — layer4/mean_morrey_diagnostic.py

Computes the mean Morrey seminorm M_α(r) = r^{-1+α}|⟨u⟩_{B_r}|³
for velocity field snapshots. Used to numerically verify the
gap characterisation from Corollary 14.6 of the proof document.

Claim: M_α(r) ≤ C₀ for some α > 0 (mean Morrey condition)
is the single remaining obstacle to Claim A for all Leray-Hopf solutions.
"""
from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Module-level constants
# ─────────────────────────────────────────────────────────────────────────────

DB_PATH = Path(__file__).parent.parent / "results" / "results.db"
DEFAULT_R_VALUES: List[float] = [0.1, 0.2, 0.314, 0.5, 0.785, 1.0, 1.571]   # in [0, 2π]
DEFAULT_ALPHA_VALUES: List[float] = [0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0]
L_BOX: float = 2 * np.pi       # periodic box side length

# Boundedness ceiling for PASS check; generous upper bound
_C0: float = 1e6


# ─────────────────────────────────────────────────────────────────────────────
# Ball mask
# ─────────────────────────────────────────────────────────────────────────────

def ball_mask(N: int, cx: int, cy: int, cz: int, r_frac: float) -> np.ndarray:
    """Boolean mask for ball of radius r_frac*L_BOX centred at grid index (cx,cy,cz).

    Uses periodic minimum-image distance on an N×N×N grid with box length L_BOX.

    Parameters
    ----------
    N      : grid size
    cx, cy, cz : centre grid indices (integers)
    r_frac : radius as fraction of L_BOX

    Returns
    -------
    mask : shape (N, N, N), dtype bool
    """
    r_phys = r_frac * L_BOX           # physical radius
    dx = L_BOX / N                    # grid spacing

    idx = np.arange(N)
    ii, jj, kk = np.meshgrid(idx, idx, idx, indexing="ij")

    # Periodic minimum-image distance (in grid-cell units)
    di = np.abs(ii - cx)
    dj = np.abs(jj - cy)
    dk = np.abs(kk - cz)

    di = np.minimum(di, N - di).astype(float)
    dj = np.minimum(dj, N - dj).astype(float)
    dk = np.minimum(dk, N - dk).astype(float)

    # Convert to physical distance and compare with r_phys
    dist_sq = (di * dx) ** 2 + (dj * dx) ** 2 + (dk * dx) ** 2
    return dist_sq <= r_phys ** 2


# ─────────────────────────────────────────────────────────────────────────────
# Spatial mean over a ball
# ─────────────────────────────────────────────────────────────────────────────

def compute_spatial_mean(
    u_field: np.ndarray,
    cx: int,
    cy: int,
    cz: int,
    r_frac: float,
) -> Tuple[np.ndarray, float]:
    """Spatial mean of u over B_r(x₀) where x₀ = (cx,cy,cz) in grid indices.

    Parameters
    ----------
    u_field : shape (3, N, N, N)
    cx, cy, cz : centre grid indices
    r_frac  : radius as fraction of L_BOX

    Returns
    -------
    mean_vec : shape (3,)
    mean_mag : float — |mean_vec|
    """
    N = u_field.shape[1]
    mask = ball_mask(N, cx, cy, cz, r_frac)
    n_pts = int(mask.sum())
    if n_pts == 0:
        return np.zeros(3), 0.0

    mean_vec = np.array([float(u_field[c][mask].mean()) for c in range(3)])
    mean_mag = float(np.linalg.norm(mean_vec))
    return mean_vec, mean_mag


# ─────────────────────────────────────────────────────────────────────────────
# Mean Morrey field over radius / alpha grid
# ─────────────────────────────────────────────────────────────────────────────

def _default_centres(N: int, n_centres: int = 8) -> List[Tuple[int, int, int]]:
    """Return up to n_centres grid indices sampling the periodic box.

    Includes box centre and a grid of evenly-spaced octant points.
    """
    pts: List[Tuple[int, int, int]] = [(N // 2, N // 2, N // 2)]  # box centre
    step = max(1, N // max(1, int(round(n_centres ** (1 / 3)))))
    for i in range(0, N, step):
        for j in range(0, N, step):
            for k in range(0, N, step):
                pt = (i, j, k)
                if pt not in pts:
                    pts.append(pt)
                if len(pts) >= n_centres:
                    return pts[:n_centres]
    return pts[:n_centres]


def compute_mean_morrey_field_random(
    u: np.ndarray,
    r_values: Sequence[float],
    alpha_values: Sequence[float],
    n_centres: int = 32,
    seed: Optional[int] = None,
) -> Tuple[Dict[Tuple[float, float], float], np.ndarray]:
    """Compute M_α(r) = r^{-1+α} * max_{random x₀} |⟨u⟩_{B_r}|³.

    Random centres avoid symmetry-node bias of fixed grids on periodic fields.
    Returns (morrey_dict, centres_xyz) where centres_xyz has shape (n_centres, 3).
    """
    N = u.shape[1]
    rng = np.random.default_rng(seed)
    # Sample n_centres random grid indices uniformly in [0, N-1]^3
    centres_xyz = rng.integers(0, N, size=(n_centres, 3))

    # For each r, find the maximum |mean| over all random centres
    max_mean_mag: Dict[float, float] = {}
    for r_phys in r_values:
        r_frac = r_phys / L_BOX
        max_mag = 0.0
        for cx, cy, cz in centres_xyz:
            _, mag = compute_spatial_mean(u, int(cx), int(cy), int(cz), r_frac)
            if mag > max_mag:
                max_mag = mag
        max_mean_mag[r_phys] = max_mag

    # Compute M_alpha for each (r, alpha) pair
    morrey_dict: Dict[Tuple[float, float], float] = {}
    for r_phys in r_values:
        max_mag = max_mean_mag[r_phys]
        max_mag_cube = max_mag ** 3
        for alpha in alpha_values:
            exponent = -1.0 + alpha
            if r_phys <= 0.0:
                m_val = float("inf") if exponent < 0 else 0.0
            elif exponent < 0 and r_phys < 1e-15:
                m_val = float("inf")
            else:
                m_val = (r_phys ** exponent) * max_mag_cube
            morrey_dict[(r_phys, alpha)] = m_val

    return morrey_dict, centres_xyz


def compute_mean_morrey_field(
    u_field: np.ndarray,
    r_values: Sequence[float],
    alpha_values: Sequence[float],
    n_centres: int = 8,
) -> Tuple[Dict[Tuple[float, float], float], Dict[float, float]]:
    """Compute M_α(r) for all (r, α) pairs, maximised over sampled centres.

    M_α(r) = r^{-1+α} * max_{x₀} |⟨u⟩_{B_r(x₀)}|³

    Parameters
    ----------
    u_field      : shape (3, N, N, N)
    r_values     : list of radii in physical units ∈ [0, L_BOX]
    alpha_values : list of α exponents
    n_centres    : number of centre points to sample

    Returns
    -------
    morrey_dict  : dict (r, alpha) → M_alpha value
    max_mean_mag : dict r → max spatial-mean magnitude over centres
    """
    N = u_field.shape[1]
    centres = _default_centres(N, n_centres)

    # For each r, find the maximum |mean| over all sampled centres
    max_mean_mag: Dict[float, float] = {}
    for r_phys in r_values:
        r_frac = r_phys / L_BOX
        max_mag = 0.0
        for (cx, cy, cz) in centres:
            _, mag = compute_spatial_mean(u_field, cx, cy, cz, r_frac)
            if mag > max_mag:
                max_mag = mag
        max_mean_mag[r_phys] = max_mag

    # Compute M_alpha for each (r, alpha) pair
    morrey_dict: Dict[Tuple[float, float], float] = {}
    for r_phys in r_values:
        max_mag = max_mean_mag[r_phys]
        max_mag_cube = max_mag ** 3
        for alpha in alpha_values:
            exponent = -1.0 + alpha
            if r_phys <= 0.0:
                m_val = float("inf") if exponent < 0 else 0.0
            elif exponent < 0 and r_phys < 1e-15:
                m_val = float("inf")
            else:
                m_val = (r_phys ** exponent) * max_mag_cube
            morrey_dict[(r_phys, alpha)] = m_val

    return morrey_dict, max_mean_mag


# ─────────────────────────────────────────────────────────────────────────────
# Scaling fit
# ─────────────────────────────────────────────────────────────────────────────

def fit_morrey_scaling(
    r_values: Sequence[float],
    morrey_dict: Dict[Tuple[float, float], float],
    alpha: float,
) -> Tuple[float, float, float]:
    """Fit M_α(r) ~ C * r^γ using log-log least-squares regression.

    Parameters
    ----------
    r_values   : physical radii used in morrey_dict
    morrey_dict: dict (r, alpha) → M_alpha value
    alpha      : the specific α to fit

    Returns
    -------
    C   : prefactor (amplitude)
    γ   : scaling exponent (γ ≥ 0 → bounded as r→0; γ < 0 → diverges)
    R2  : coefficient of determination of the log-log fit
    """
    r_arr = np.asarray([r for r in r_values if r > 0.0], dtype=float)
    m_arr = np.array([morrey_dict.get((r, alpha), np.nan) for r in r_arr])

    # Keep only finite, positive values for log-log fit
    valid = np.isfinite(m_arr) & (m_arr > 0.0) & (r_arr > 0.0)
    if valid.sum() < 2:
        return 0.0, 0.0, 0.0

    log_r = np.log(r_arr[valid])
    log_m = np.log(m_arr[valid])

    # Linear fit: log_m = log_C + γ * log_r
    A = np.column_stack([np.ones_like(log_r), log_r])
    coeffs, residuals, rank, _ = np.linalg.lstsq(A, log_m, rcond=None)
    log_C, gamma = float(coeffs[0]), float(coeffs[1])
    C = float(np.exp(log_C))

    # R² = 1 - SS_res / SS_tot
    log_m_pred = log_C + gamma * log_r
    ss_res = float(np.sum((log_m - log_m_pred) ** 2))
    ss_tot = float(np.sum((log_m - log_m.mean()) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0

    return C, gamma, R2


# ─────────────────────────────────────────────────────────────────────────────
# Alpha threshold detection
# ─────────────────────────────────────────────────────────────────────────────

def find_alpha_threshold(
    r_values: Sequence[float],
    morrey_dict: Dict[Tuple[float, float], float],
    alpha_values: Sequence[float],
) -> Tuple[float, dict]:
    """Find smallest α such that M_α(r) is bounded (γ ≥ 0) as r → 0.

    Parameters
    ----------
    r_values     : physical radii used in morrey_dict
    morrey_dict  : dict (r, alpha) → M_alpha value
    alpha_values : α values to test (in ascending order)

    Returns
    -------
    alpha_min  : smallest α with γ ≥ 0 (nan if none found)
    fit_details: dict mapping each alpha to (C, gamma, R2)
    """
    fit_details: dict = {}
    alpha_min = float("nan")

    for alpha in sorted(alpha_values):
        C, gamma, R2 = fit_morrey_scaling(r_values, morrey_dict, alpha)
        fit_details[alpha] = {"C": C, "gamma": gamma, "R2": R2}
        if np.isnan(alpha_min) and gamma >= 0.0:
            alpha_min = float(alpha)

    return alpha_min, fit_details


# ─────────────────────────────────────────────────────────────────────────────
# Full analysis pipeline
# ─────────────────────────────────────────────────────────────────────────────

def run_mean_morrey_analysis(
    exp_id: str,
    u_snapshots: Sequence[np.ndarray],
    t_values: Sequence[float],
    grid_N: int,
    r_values: Optional[Sequence[float]] = None,
    alpha_values: Optional[Sequence[float]] = None,
    ic_type: str = "unknown",
    n_centres_random: int = 0,
    seed: Optional[int] = None,
) -> dict:
    """Full mean Morrey analysis over velocity snapshots.

    n_centres_random > 0: use random-centre sampling (avoids symmetry bias);
    seed controls reproducibility.  n_centres_random <= 0: deterministic grid.
    Returns dict with exp_id, ic_type, grid_N, alpha_threshold, morrey_max,
    morrey_at_r_min, scaling_exponent, verdict, notes.
    """
    if r_values is None:
        r_values = DEFAULT_R_VALUES
    if alpha_values is None:
        alpha_values = DEFAULT_ALPHA_VALUES

    r_list = list(r_values)
    a_list = list(alpha_values)

    # Accumulate maximum M_alpha over all snapshots
    morrey_accum: Dict[Tuple[float, float], float] = {
        (r, a): 0.0 for r in r_list for a in a_list
    }

    n_snapshots = 0
    for snap in u_snapshots:
        if n_centres_random > 0:
            mm, _ = compute_mean_morrey_field_random(
                snap, r_list, a_list,
                n_centres=n_centres_random,
                seed=seed,
            )
        else:
            mm, _ = compute_mean_morrey_field(snap, r_list, a_list)
        for key, val in mm.items():
            if np.isfinite(val):
                morrey_accum[key] = max(morrey_accum[key], val)
            else:
                morrey_accum[key] = float("inf")
        n_snapshots += 1

    # Find alpha threshold and scaling exponent
    alpha_threshold, fit_details = find_alpha_threshold(r_list, morrey_accum, a_list)

    # Morrey value at smallest radius (most sensitive probe)
    r_min = min(r_list)
    morrey_at_r_min = float(
        max(morrey_accum.get((r_min, a), 0.0) for a in a_list)
    )

    # Overall maximum Morrey value across all (r, alpha)
    finite_vals = [v for v in morrey_accum.values() if np.isfinite(v)]
    morrey_max = float(max(finite_vals)) if finite_vals else float("inf")

    # Representative scaling exponent: from alpha_threshold fit if available
    if np.isfinite(alpha_threshold) and alpha_threshold in fit_details:
        scaling_exponent = float(fit_details[alpha_threshold]["gamma"])
    else:
        scaling_exponent = float("nan")

    # PASS: found finite alpha_threshold in valid range with bounded seminorm
    verdict = "FAIL"
    if (
        np.isfinite(alpha_threshold)
        and alpha_threshold <= 2.0
        and np.isfinite(morrey_max)
        and morrey_max <= _C0
    ):
        verdict = "PASS"

    at_str = f"{alpha_threshold:.3f}" if np.isfinite(alpha_threshold) else "nan"
    se_str = f"{scaling_exponent:.3f}" if np.isfinite(scaling_exponent) else "nan"
    sampling_str = (
        f"random(n={n_centres_random}, seed={seed})"
        if n_centres_random > 0
        else "deterministic"
    )
    notes = (
        f"alpha_threshold={at_str}, "
        f"morrey_max={morrey_max:.4e}, "
        f"scaling_exponent={se_str}, "
        f"n_snapshots={n_snapshots}, "
        f"sampling={sampling_str}"
    )

    return {
        "exp_id":           exp_id,
        "ic_type":          ic_type,
        "grid_N":           grid_N,
        "alpha_threshold":  alpha_threshold,
        "morrey_max":       morrey_max,
        "morrey_at_r_min":  morrey_at_r_min,
        "scaling_exponent": scaling_exponent,
        "verdict":          verdict,
        "notes":            notes,
        "timestamp":        datetime.now(timezone.utc).isoformat(),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Database helpers
# ─────────────────────────────────────────────────────────────────────────────

_CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS mean_morrey_results (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    exp_id            TEXT    NOT NULL,
    timestamp         TEXT,
    grid_N            INTEGER,
    ic_type           TEXT,
    alpha_threshold   REAL,
    morrey_max        REAL,
    morrey_at_r_min   REAL,
    scaling_exponent  REAL,
    verdict           TEXT,
    notes             TEXT
)
"""

_INSERT_RESULT = """
INSERT INTO mean_morrey_results
    (exp_id, timestamp, grid_N, ic_type,
     alpha_threshold, morrey_max, morrey_at_r_min, scaling_exponent,
     verdict, notes)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


def init_db(conn: sqlite3.Connection) -> None:
    """Create mean_morrey_results table if it does not already exist."""
    conn.execute(_CREATE_TABLE)
    conn.commit()


def log_result(conn: sqlite3.Connection, result: dict) -> None:
    """Insert a result dict into the mean_morrey_results table.

    Calls init_db automatically so the table is guaranteed to exist.

    Parameters
    ----------
    conn   : open sqlite3 connection
    result : output of run_mean_morrey_analysis
    """
    init_db(conn)
    conn.execute(
        _INSERT_RESULT,
        (
            result.get("exp_id", ""),
            result.get("timestamp", datetime.now(timezone.utc).isoformat()),
            result.get("grid_N", 0),
            result.get("ic_type", "unknown"),
            result.get("alpha_threshold", None),
            result.get("morrey_max", None),
            result.get("morrey_at_r_min", None),
            result.get("scaling_exponent", None),
            result.get("verdict", "FAIL"),
            result.get("notes", ""),
        ),
    )
    conn.commit()
