"""
layer4/alignment_bridge.py — Numerical Bridge Lemma Diagnostic
================================================================
Tests the numerical plausibility of the "Bridge Lemma": whether the
POINTWISE coherence of the vorticity direction field ê = ω/|ω| near the
point of maximal vorticity is controlled by its L²-AVERAGED coherence
(σ-type quantity), at the self-similar scale r* = M^{-1/2}.

Physics / definitions (given ω = (wx, wy, wz) on a periodic N³ grid,
domain [0, 2π)³, spacing dx = 2π/N):

  M            = max |ω|                         (pointwise magnitude)
  r*           = M^{-1/2}                        (self-similar scale)
  x*           = argmax |ω|                       (blow-up candidate site)
  ê            = ω / |ω|  where |ω| ≥ floor = 1e-12·M, else ê = 0 (excluded)
  |∇ê|²        = Σ_{i,k} (∂_k ê_i)²               (spectral derivatives)
  mask4        = {|ω| ≥ M/4}, mask2 = {|ω| ≥ M/2}  (high-vorticity sets)
  B(r*)        = periodic (minimum-image) ball of radius r* about x*

  σ*  = (1/r*) · Σ_{B(r*) ∩ mask4} |∇ê|² dx³        (L²-averaged coherence)
  L∞* = max_{B(r*) ∩ mask2} |∇ê|                     (pointwise coherence)

  Bridge ratio  R = (L∞* · r*) / max(√σ*, 1e-10)

The Bridge Lemma is numerically plausible iff R stays O(1)–O(100) and does
not diverge over the course of a run (verdict thresholds documented in
run_bridge_experiment).

Numerical caveat (Gibbs artifacts): ê is discontinuous across the
floor/mask boundary (set identically to zero outside the valid region),
so a *global* FFT derivative of ê exhibits Gibbs-type ringing near that
boundary. This is mitigated — not eliminated — by only ever reading
|∇ê|² on the high-vorticity ball ∩ mask4/mask2 sets, i.e. well inside the
region where ê is smooth and away from the ê=0 truncation boundary
(mask4/mask2 require |ω| ≥ M/4 or M/2, which is far above the 1e-12·M
floor). Users of this module should treat σ* and L∞* as scale-invariant
DIAGNOSTIC indicators, not as spectrally-exact quantities near very
sharp vorticity fronts.

House rules: spectral (FFT) derivatives only, no finite differences;
no hardcoded fluid properties (nu comes from layer3/route2_3D.py's own
default); all functions pure / immutable (return new arrays, never
mutate inputs); every logged result carries claim_id + key_metric +
verdict; random seed always set and logged even when the ICs used here
are deterministic.

Run:
  python layer4/alignment_bridge.py --ic tg --N 32 --n-steps 200 --verbose
  python layer4/alignment_bridge.py --ic shear --N 32 --n-steps 200 --verbose
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from datetime import datetime

import numpy as np

# ── path setup (matches layer4/delta_extractor.py convention) ────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

import route2_3D as r2  # noqa: E402  — Route 2 3D solver (reused, not reimplemented)
from layer4.geometric_disorder import spectral_wavenumbers, _ball_mask  # noqa: E402
from layer3.blowup_search_3D import (  # noqa: E402 — Biot–Savart tube IC, reused not reimplemented
    make_vortex_omega,
    biot_savart_spectral,
    _AMP_DEFAULT as _ADV_AMP_DEFAULT,
    _TUBE_SIGMA as _ADV_SIGMA_DEFAULT,
    _TUBE_SEP as _ADV_SEP_DEFAULT,
)

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_SEED = 42
_IC_TAGS = {"tg": "TG", "shear": "SH", "adv": "ADV"}
_IC_CLAIM_SUFFIX = {"tg": "tg", "shear": "shear", "adv": "adv"}


# =============================================================================
# Core Bridge Lemma primitives
# =============================================================================

def compute_direction_field(
    wx: np.ndarray, wy: np.ndarray, wz: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Direction field ê = ω/|ω|, restricted to the valid (non-degenerate) set.

    floor = 1e-12 · M, M = max|ω|.  Where |ω| < floor (or M == 0 identically),
    ê is set to (0, 0, 0) and marked invalid via valid_mask=False — such
    points must be excluded from all downstream diagnostics.

    Does not mutate wx, wy, wz.

    Returns
    -------
    ex, ey, ez  : components of ê, each shape == wx.shape
    mag         : |ω|, shape == wx.shape
    valid_mask  : bool array, True where ê is well-defined
    """
    mag = np.sqrt(wx**2 + wy**2 + wz**2)
    M = float(np.max(mag))

    if M <= 0.0:
        valid_mask = np.zeros_like(mag, dtype=bool)
    else:
        floor = 1e-12 * M
        valid_mask = mag >= floor

    safe_mag = np.where(valid_mask, mag, 1.0)  # avoid 0/0 outside valid set
    ex = np.where(valid_mask, wx / safe_mag, 0.0)
    ey = np.where(valid_mask, wy / safe_mag, 0.0)
    ez = np.where(valid_mask, wz / safe_mag, 0.0)
    return ex, ey, ez, mag, valid_mask


def grad_ehat_sq(
    ex: np.ndarray, ey: np.ndarray, ez: np.ndarray,
    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray,
) -> np.ndarray:
    """|∇ê|² = Σ_{i,k} (∂_k ê_i)², via spectral (FFT) derivatives.

    Computes all 9 components ∂_k ê_i (i, k ∈ {x, y, z}) and sums their
    squares.  kx, ky, kz are FFT wavenumber arrays consistent with the
    periodic [0, 2π)^N grid convention used throughout this program
    (see layer4.geometric_disorder.spectral_wavenumbers).

    See module docstring for the Gibbs-artifact caveat: this is a global
    FFT derivative of a field that is truncated to zero outside the valid
    vorticity set, so ringing can appear near that truncation boundary.
    Callers should only read this array on masks well inside the
    high-vorticity region (see bridge_diagnostics).

    Does not mutate ex, ey, ez.
    """
    total = np.zeros_like(ex, dtype=float)
    for comp in (ex, ey, ez):
        comp_hat = np.fft.fftn(comp)
        dx_ = np.real(np.fft.ifftn(1j * kx * comp_hat))
        dy_ = np.real(np.fft.ifftn(1j * ky * comp_hat))
        dz_ = np.real(np.fft.ifftn(1j * kz * comp_hat))
        total = total + dx_**2 + dy_**2 + dz_**2
    return total


def _validate_omega_grid(wx: np.ndarray, wy: np.ndarray, wz: np.ndarray, dx: float) -> int:
    """Shared shape/spacing validation for a periodic N³ vorticity snapshot.

    Raises ValueError on shape mismatch, a non-cubic grid, or a dx that is
    inconsistent with the periodic [0, 2π)^N convention (dx must equal
    2π/N).  Returns N.  Factored out of bridge_diagnostics so that
    layer4/k_correlation.py's k_correlation_diagnostics can reuse the exact
    same contract without duplicating it.
    """
    if not (wx.shape == wy.shape == wz.shape):
        raise ValueError(
            f"wx, wy, wz must share a shape; got {wx.shape}, {wy.shape}, {wz.shape}"
        )
    if wx.ndim != 3 or wx.shape[0] != wx.shape[1] or wx.shape[1] != wx.shape[2]:
        raise ValueError(f"vorticity fields must be N x N x N; got shape {wx.shape}")

    N = wx.shape[0]
    expected_dx = 2.0 * np.pi / N
    if abs(dx - expected_dx) > 1e-9 * max(1.0, expected_dx):
        raise ValueError(
            f"dx={dx} inconsistent with N={N} periodic grid (expected dx=2π/N={expected_dx})"
        )
    return N


def bridge_diagnostics(wx: np.ndarray, wy: np.ndarray, wz: np.ndarray, dx: float) -> dict:
    """Compute the full Bridge Lemma diagnostic dict for one vorticity snapshot.

    Parameters
    ----------
    wx, wy, wz : vorticity components, shape (N, N, N), periodic domain [0,2π)^3
    dx         : grid spacing, must equal 2π/N (validated)

    Returns dict with keys: M, r_star, x_star, sigma_star, linf_grad_ehat,
    bridge_ratio.  Raises ValueError on shape mismatch, bad dx, or an
    identically-zero vorticity field (r* = M^{-1/2} undefined at M=0).

    Does not mutate wx, wy, wz.
    """
    N = _validate_omega_grid(wx, wy, wz, dx)

    mag = np.sqrt(wx**2 + wy**2 + wz**2)
    M = float(np.max(mag))
    if M < 1e-14:
        raise ValueError(
            "bridge_diagnostics: vorticity field is identically zero (M ≈ 0); "
            "r* = M^{-1/2} and the Bridge Lemma direction field are undefined."
        )

    ex, ey, ez, mag, valid_mask = compute_direction_field(wx, wy, wz)
    kx, ky, kz = spectral_wavenumbers(N)
    grad_sq = grad_ehat_sq(ex, ey, ez, kx, ky, kz)

    r_star = M ** -0.5
    idx = np.unravel_index(int(np.argmax(mag)), mag.shape)
    x_star = (int(idx[0]), int(idx[1]), int(idx[2]))

    mask4 = mag >= (M / 4.0)
    mask2 = mag >= (M / 2.0)
    ball = _ball_mask(N, idx, r_star)

    region4 = ball & mask4 & valid_mask
    region2 = ball & mask2 & valid_mask

    dx3 = dx ** 3
    if np.any(region4):
        sigma_star = float((1.0 / r_star) * np.sum(grad_sq[region4]) * dx3)
    else:
        sigma_star = 0.0

    if np.any(region2):
        linf_grad_ehat = float(np.sqrt(np.max(grad_sq[region2])))
    else:
        linf_grad_ehat = 0.0

    bridge_ratio = float(
        (linf_grad_ehat * r_star) / max(np.sqrt(max(sigma_star, 0.0)), 1e-10)
    )

    return {
        "M": M,
        "r_star": float(r_star),
        "x_star": x_star,
        "sigma_star": sigma_star,
        "linf_grad_ehat": linf_grad_ehat,
        "bridge_ratio": bridge_ratio,
    }


# =============================================================================
# Vorticity from Route 2 solver velocity (spectral curl; not a re-derivation
# of the solver itself — the solver only carries u, v, w)
# =============================================================================

def _vorticity_from_velocity(
    u: np.ndarray, v: np.ndarray, w: np.ndarray,
    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ω = ∇ × u, computed spectrally.  Does not mutate u, v, w."""
    u_hat = np.fft.fftn(u)
    v_hat = np.fft.fftn(v)
    w_hat = np.fft.fftn(w)

    dwdy = np.real(np.fft.ifftn(1j * ky * w_hat))
    dvdz = np.real(np.fft.ifftn(1j * kz * v_hat))
    dudz = np.real(np.fft.ifftn(1j * kz * u_hat))
    dwdx = np.real(np.fft.ifftn(1j * kx * w_hat))
    dvdx = np.real(np.fft.ifftn(1j * kx * v_hat))
    dudy = np.real(np.fft.ifftn(1j * ky * u_hat))

    wx = dwdy - dvdz
    wy = dudz - dwdx
    wz = dvdx - dudy
    return wx, wy, wz


# =============================================================================
# Adversarial IC: anti-parallel vortex tubes (Biot–Savart, reused from
# layer3/blowup_search_3D.py — velocity construction only; the solver itself
# remains route2_3D, i.e. exact Prize NS with nu=const, not Route 1's mu(T))
# =============================================================================

def adv_vortex_tubes_ic(
    N: int,
    amp: float = _ADV_AMP_DEFAULT,
    sigma: float = _ADV_SIGMA_DEFAULT,
    separation: float = _ADV_SEP_DEFAULT,
    seed: int = _DEFAULT_SEED,
) -> dict:
    """Anti-parallel Gaussian vortex-tube IC, recovered as velocity via
    Biot–Savart, in the same {'u','v','w','name'} dict shape as
    taylor_green_ic / shear_layer_ic / random_div_free_ic (route2_3D.py).

    Reuses layer3.blowup_search_3D.make_vortex_omega (vorticity construction)
    and .biot_savart_spectral (u_hat = i·k×ω̂/|k|²) verbatim — no
    reimplementation. Wavenumber convention (np.fft.fftfreq(N, d=1/N),
    meshgrid indexing='ij') matches route2_3D.make_grid_3D exactly, so the
    Biot–Savart inversion is consistent with the Route 2 solver's own grid.

    The caller (run_bridge_experiment, via r2.set_ic) Leray-projects this
    velocity before the run starts, mirroring the 'tg'/'shear' IC path.
    """
    kx, ky, kz, k2 = r2.make_grid_3D(N)
    omega_x, omega_y, omega_z = make_vortex_omega(
        N, amp=amp, sigma=sigma, separation=separation, seed=seed
    )
    u, v, w = biot_savart_spectral(omega_x, omega_y, omega_z, kx, ky, kz, k2)
    return {
        "u": u.astype(np.float64),
        "v": v.astype(np.float64),
        "w": w.astype(np.float64),
        "name": "adv_vortex_tubes",
    }


def _k_correlation_snapshot(wx: np.ndarray, wy: np.ndarray, wz: np.ndarray, dx: float) -> dict:
    """Lazy-import wrapper around layer4.k_correlation.k_correlation_diagnostics.

    The import is deferred to call time (rather than module load time)
    because layer4/k_correlation.py imports from this module — a
    module-level import here would be circular.  By the time this function
    is first called, alignment_bridge has already finished loading, so the
    import resolves cleanly.  Does not mutate wx, wy, wz.
    """
    from layer4.k_correlation import k_correlation_diagnostics
    return k_correlation_diagnostics(wx, wy, wz, dx)


# =============================================================================
# Production experiment: run Route 2 3D solver, track Bridge diagnostics
# =============================================================================

def run_bridge_experiment(
    ic: str = "tg",
    N: int = 32,
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
    adv_amp: float | None = None,
) -> dict:
    """Run the existing Route 2 3D solver and track Bridge Lemma diagnostics.

    Reuses layer3/route2_3D.py (make_solver, taylor_green_ic / shear_layer_ic,
    set_ic, compute_cfl_dt, solver_step) — no solver logic is reimplemented
    here.  Velocity evolves under exact Prize NS (Route 2's velocity RHS does
    not depend on eps_param); vorticity is derived from velocity every
    `record_every` steps and fed through bridge_diagnostics.

    Parameters
    ----------
    ic            : 'tg' (Taylor–Green), 'shear' (shear layer), or 'adv'
                    (anti-parallel Biot–Savart vortex tubes — adversarial
                    alignment-stress IC, see adv_vortex_tubes_ic)
    N             : grid points per side (N³ total)
    n_steps       : number of solver steps to run
    nu            : kinematic viscosity; defaults to route2_3D._NU_DEFAULT
    eps_param     : Route 2 ε (does not affect velocity/vorticity)
    record_every  : compute Bridge diagnostics every this-many steps
    seed          : random seed, set and logged even though the ICs used
                    here (TG, shear) are deterministic — no randomness is
                    actually drawn, but the seed is recorded for provenance
    verbose       : print per-record progress
    db_path       : if given, log the result via log_result()
    exp_id_override, claim_id_override : if given, override the default
                    'EXP-L4-BRIDGE-{TAG}-001' / 'bridge-lemma-{ic}' naming
                    (used for higher-resolution / repeat production runs,
                    e.g. 'EXP-L4-BRIDGE-TG-002', 'bridge-lemma-tg-64')
    adv_amp       : only used when ic == 'adv'; overrides
                    adv_vortex_tubes_ic's default peak-vorticity amplitude
                    (layer3.blowup_search_3D._AMP_DEFAULT) — used to build
                    escalated-amplitude adversarial ICs (e.g. the S34
                    combined stretching+alignment stress test). Ignored for
                    ic in {'tg', 'shear'}.

    Returns a dict with the full timeseries and a PASS/FAIL verdict.

    Verdict PASS iff:
      - every recorded sigma_star is finite
      - every recorded bridge_ratio is finite
      - max(bridge_ratio) <= 1000
      - final bridge_ratio / median bridge_ratio < 10  (no late-time divergence)
    """
    if ic not in _IC_TAGS:
        raise ValueError(f"Unknown ic='{ic}'; use 'tg', 'shear', or 'adv'.")
    if n_steps <= 0:
        raise ValueError(f"n_steps must be > 0, got {n_steps}")

    if nu is None:
        nu = r2._NU_DEFAULT

    # Seed set and logged for reproducibility provenance; the 'tg'/'shear'
    # ICs are deterministic (RNG unused downstream); 'adv' does draw a small
    # symmetry-breaking perturbation from this seed (see make_vortex_omega).
    _rng = np.random.default_rng(seed)

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
                bd = bridge_diagnostics(wx, wy, wz, dx)
            except ValueError:
                # Vorticity has (numerically) vanished; record a degenerate,
                # explicitly-flagged zero point rather than dropping the step.
                bd = {
                    "M": 0.0, "r_star": float("inf"), "x_star": (0, 0, 0),
                    "sigma_star": 0.0, "linf_grad_ehat": 0.0, "bridge_ratio": 0.0,
                }
            try:
                kc = _k_correlation_snapshot(wx, wy, wz, dx)
            except ValueError:
                # Same degenerate (identically-zero vorticity) fallback as bd.
                kc = {
                    "c_K": 0.0, "f_ang": 0.0, "mean_grad_m2": 0.0, "mean_w": 0.0,
                    "E_volume_fraction": 0.0, "split_residual": 0.0,
                }
            record = {"step": step_idx, "t": solver["t"], **bd, **kc}
            timeseries.append(record)
            if verbose:
                print(
                    f"  step {step_idx:4d}  t={record['t']:.4f}  M={bd['M']:.4e}  "
                    f"r*={bd['r_star']:.4e}  sigma*={bd['sigma_star']:.4e}  "
                    f"R={bd['bridge_ratio']:.4e}  c_K={kc['c_K']:.4e}  "
                    f"f_ang={kc['f_ang']:.4e}"
                )

    wall = time.perf_counter() - t0

    bridge_ratios = np.array([rec["bridge_ratio"] for rec in timeseries], dtype=float)
    sigma_stars = np.array([rec["sigma_star"] for rec in timeseries], dtype=float)
    M_vals = np.array([rec["M"] for rec in timeseries], dtype=float)
    r_star_vals = np.array([rec["r_star"] for rec in timeseries], dtype=float)

    all_sigma_finite = bool(np.all(np.isfinite(sigma_stars)))
    all_ratio_finite = bool(np.all(np.isfinite(bridge_ratios)))
    max_br = float(np.max(bridge_ratios)) if bridge_ratios.size else float("nan")
    median_br = float(np.median(bridge_ratios)) if bridge_ratios.size else float("nan")
    final_br = float(bridge_ratios[-1]) if bridge_ratios.size else float("nan")
    ratio_final_median = final_br / max(median_br, 1e-12) if np.isfinite(final_br) else float("inf")

    verdict = "PASS" if (
        all_sigma_finite
        and all_ratio_finite
        and np.isfinite(max_br) and max_br <= 1000.0
        and np.isfinite(ratio_final_median) and ratio_final_median < 10.0
    ) else "FAIL"

    tag = _IC_TAGS[ic]
    exp_id = exp_id_override or f"EXP-L4-BRIDGE-{tag}-001"
    claim_id = claim_id_override or f"bridge-lemma-{_IC_CLAIM_SUFFIX[ic]}"

    sigma_min = float(np.min(sigma_stars)) if sigma_stars.size else 0.0
    sigma_max = float(np.max(sigma_stars)) if sigma_stars.size else 0.0

    # K-correlation summary stats (see layer4/k_correlation.py) — computed
    # here from the already-recorded per-step c_K/f_ang/split_residual so
    # that callers (e.g. layer4/k_correlation.py's kcorr_summary_from_result)
    # don't need to re-run the solver.
    c_K_vals = np.array([rec["c_K"] for rec in timeseries], dtype=float)
    f_ang_vals = np.array([rec["f_ang"] for rec in timeseries], dtype=float)
    split_residual_vals = np.array([rec["split_residual"] for rec in timeseries], dtype=float)

    result = {
        "exp_id": exp_id,
        "claim_id": claim_id,
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic,
        "N": N,
        "n_steps": n_steps,
        "seed": seed,
        "nu": nu,
        "eps_param": eps_param,
        "t_final": solver["t"],
        "M_max": float(np.max(M_vals)) if M_vals.size else 0.0,
        "r_star_min": float(np.min(r_star_vals)) if r_star_vals.size else float("inf"),
        "sigma_star_min": sigma_min,
        "sigma_star_max": sigma_max,
        "bridge_ratio_max": max_br,
        "bridge_ratio_median": median_br,
        "bridge_ratio_final": final_br,
        "c_K_max": float(np.max(c_K_vals)) if c_K_vals.size else float("nan"),
        "c_K_median": float(np.median(c_K_vals)) if c_K_vals.size else float("nan"),
        "c_K_final": float(c_K_vals[-1]) if c_K_vals.size else float("nan"),
        "f_ang_median": float(np.median(f_ang_vals)) if f_ang_vals.size else float("nan"),
        "split_residual_min": float(np.min(split_residual_vals)) if split_residual_vals.size else float("nan"),
        "split_residual_max": float(np.max(split_residual_vals)) if split_residual_vals.size else float("nan"),
        "verdict": verdict,
        "wall_time_s": wall,
        "key_metric": (
            f"max_bridge_ratio={max_br:.4f} median_bridge_ratio={median_br:.4f} "
            f"final_bridge_ratio={final_br:.4f} "
            f"sigma_star_range=[{sigma_min:.3e},{sigma_max:.3e}]"
        ),
        "timeseries": timeseries,
    }

    if verbose:
        print(f"\n  Bridge Lemma experiment ({ic}, N={N}, n_steps={n_steps})")
        print(f"    max R:    {max_br:.4f}")
        print(f"    median R: {median_br:.4f}")
        print(f"    final R:  {final_br:.4f}")
        print(f"    verdict:  {verdict}")

    if db_path is not None:
        log_result(db_path, result)

    return result


# =============================================================================
# Results logging (layer4-module table pattern; see delta_extractor.py,
# mean_morrey_diagnostic.py)
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_bridge_experiments (
            id                   INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id               TEXT UNIQUE NOT NULL,
            claim_id             TEXT NOT NULL,
            timestamp            TEXT NOT NULL,
            ic                   TEXT,
            N                    INTEGER,
            n_steps              INTEGER,
            seed                 INTEGER,
            nu                   REAL,
            eps_param            REAL,
            t_final              REAL,
            M_max                REAL,
            r_star_min           REAL,
            sigma_star_min       REAL,
            sigma_star_max       REAL,
            bridge_ratio_max     REAL,
            bridge_ratio_median  REAL,
            bridge_ratio_final   REAL,
            verdict              TEXT NOT NULL,
            wall_time_s          REAL,
            key_metric           TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_bridge_timeseries (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id           TEXT NOT NULL,
            step             INTEGER,
            t                REAL,
            M                REAL,
            r_star           REAL,
            sigma_star       REAL,
            linf_grad_ehat   REAL,
            bridge_ratio     REAL
        )
    """)
    conn.commit()
    return conn


def log_result(db_path: str, result: dict) -> None:
    """Append one Bridge Lemma experiment (+ its timeseries) to results.db.

    Follows the layer4 per-module table convention (see delta_extractor.py,
    mean_morrey_diagnostic.py): an UPSERT of the summary row keyed on
    exp_id, plus append-only timeseries rows.
    """
    conn = _ensure_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "ic", "N", "n_steps", "seed",
        "nu", "eps_param", "t_final", "M_max", "r_star_min",
        "sigma_star_min", "sigma_star_max",
        "bridge_ratio_max", "bridge_ratio_median", "bridge_ratio_final",
        "verdict", "wall_time_s", "key_metric",
    ]
    values = tuple(result.get(c) for c in cols)
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO layer4_bridge_experiments ({', '.join(cols)}) "
        f"VALUES ({ph})",
        values,
    )
    for rec in result.get("timeseries", []):
        conn.execute(
            """
            INSERT INTO layer4_bridge_timeseries
            (exp_id, step, t, M, r_star, sigma_star, linf_grad_ehat, bridge_ratio)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                result["exp_id"], rec["step"], rec["t"], rec["M"],
                rec["r_star"], rec["sigma_star"], rec["linf_grad_ehat"],
                rec["bridge_ratio"],
            ),
        )
    conn.commit()
    conn.close()


# =============================================================================
# CLI entry point
# =============================================================================

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(description="Bridge Lemma numerical diagnostic")
    parser.add_argument("--ic", choices=["tg", "shear", "adv"], default="tg")
    parser.add_argument("--N", type=int, default=32)
    parser.add_argument("--n-steps", type=int, default=200)
    parser.add_argument("--record-every", type=int, default=10)
    parser.add_argument("--seed", type=int, default=_DEFAULT_SEED)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db", default=_DEFAULT_DB)
    args = parser.parse_args()

    r = run_bridge_experiment(
        ic=args.ic, N=args.N, n_steps=args.n_steps,
        record_every=args.record_every, seed=args.seed,
        eps_param=args.eps, verbose=args.verbose, db_path=args.db,
    )
    print(f"\nVerdict: {r['verdict']}  key_metric: {r['key_metric']}")


if __name__ == "__main__":
    main()
