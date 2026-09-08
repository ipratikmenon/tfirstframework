"""
layer4/aposteriori_verification.py — A Posteriori Regularity Verification Diagnostic
=====================================================================================

*** SCOPE — READ BEFORE USING ANY OUTPUT OF THIS MODULE ***

Implements a "verification-style diagnostic" in the shape of the a
posteriori regularity verification technique for 3D incompressible NS
(Chernyshenko–Constantin–Robinson–Titi 2007; Morosi–Pizzocchero). Idea:
given an approximate (numerical) trajectory u_a on [0,T], form its
residual against the NS momentum equation; if the accumulated residual,
amplified by a Gronwall factor built from the computed solution's own
norms, stays below a threshold, one can conclude the EXACT solution
sharing u_a's initial datum is regular on [0,T].

**THIS MODULE DOES NOT PRODUCE A PROOF.** A rigorous computer-assisted
proof in the CCRT/Morosi–Pizzocchero sense requires rigorous (interval-
arithmetic / validated-numerics) bounds on the residual and every norm
entering the Gronwall estimate. Everything here is ordinary IEEE-754
floating point. Three choices are NOT rigorously justified, flagged at
every call site:

  1. Amplification rate A(t) = C·‖∇u_a(t)‖_{L²} (`amplification_rate`) — a
     defensible functional form (energy method), but C is NOT rigorously
     derived; C=1.0 by default, treat any change as a modeling choice.
  2. `default_threshold` is a chosen ν-proportional scale, NOT a rigorously
     derived local-existence radius — a modeling decision.
  3. The residual's ∂_t u_a is a TIME-differenced estimate between
     consecutive solver states (`estimate_dudt_time_diff`) — the honest
     residual of the discrete trajectory, including its own time-
     integration error (the "spectral only" house rule concerns SPATIAL
     discretization only).

Consequently `verification_margin` returns a "verification_margin" — never
a "proof_certificate" — and `run_verification_experiment` reports a
"diagnostic satisfied" verdict, NEVER "regularity proved". A PASS is
evidence consistent with regularity on the run window and the correct
precursor to a rigorous a posteriori proof — it is NOT that proof, and
must never be represented as one (results.db, docstrings, reports alike).

*** KNOWN VACUITY OF THE LONG-WINDOW MARGIN — READ BEFORE INTERPRETING ANY
    `margin_final` / `margin_max` OVER A LONG RUN ***

`verification_margin` accumulates the Gronwall factor exp(∫_0^T A(s) ds)
over the WHOLE run window [0, T]. For production runs over T ~ O(10) time
units with A(t) ~ O(1)–O(3.5), that factor is exp(~10–70) ~ 1e4–1e30. The
long-window margin is therefore amplified by an astronomically large,
essentially arbitrary constant, and PASS (margin < 1 throughout) becomes
practically unreachable regardless of how small the residual actually is.
A long-window PARTIAL/FAIL verdict under these conditions carries NO
INFORMATION about the quality of the computed trajectory — it is vacuous
by construction, not a finding about the solution. This is exactly why the
CCRT/Morosi–Pizzocchero a posteriori literature restricts verification to
SHORT time windows where the amplification factor stays O(1)–O(10).

The fix is `windowed_verification` / `windowed_summary` below: a sliding,
RESTARTING verification where the Gronwall energy E is reset to 0 at the
start of every window of length W, so each window's amplification factor
is only exp(∫_{t_k}^{t_k+W} A ds) — the short-window quantity that can
actually be satisfied or violated informatively. **The windowed margin is
the meaningful diagnostic quantity; the long-window margin above is kept
only for backward compatibility / historical comparison and is known to be
vacuous over long runs.** Windowing fixes vacuity, not rigor: the windowed
margin is still floating point, still uses the same non-rigorous C in
`amplification_rate` and the same modeling-choice `default_threshold` — it
remains a VERIFICATION-STYLE DIAGNOSTIC, NOT a proof.

PHYSICS: r = ∂_t u_a + (u_a·∇)u_a + ∇p_a − ν·Δu_a. `ns_residual` evaluates
the spatial operator exactly as layer3/route2_3D.py's own solver does
(`_rhs_velocity`, imported verbatim: nonlinear term, Leray projection
implicitly removing ∇p_a, spectral viscous term — none reimplemented);
∂_t u_a comes from `estimate_dudt_time_diff`. `residual_norms` reports
‖r‖_{L²} and a spectral ‖r‖_{H^{-1}} ((1+|k|²)^{-1} energy weight) — the
energy estimate for e=u−u_a pairs the H^{-1} part of r against ‖e‖_{L²}
(one derivative moved by parts), so `verification_margin` is fed H^{-1}.
`verification_margin` accumulates dE/dt = A(t)E(t) + ‖r(t)‖, E(0)=0, via
the exact frozen-coefficient (exponential-Euler) step per sub-interval:
E_{i+1}=E_i·exp(A_i dt)+r_i·(exp(A_i dt)−1)/A_i (A_i≠0), else E_i+r_i·dt —
which reduces EXACTLY to E(T)=r0·T (A≡0) or r0·(e^{aT}−1)/a (A≡a>0) for
constant series (the correctness anchor, tested in
test_aposteriori_verification.py). margin(t)=E(t)/threshold; "satisfied"
over [0,T] iff margin(t)<1 throughout (never "regularity proved").

House rules: spectral (FFT) derivatives only for spatial operators; no
finite differences in space; no hardcoded fluid properties; all public
functions pure/immutable; every logged result carries claim_id +
key_metric + verdict; seed always set/logged. Reads ONLY from
layer3/route2_3D.py — does not modify it, layer3/*, proofs/*,
alignment_bridge.py, k_correlation.py, depletion_diagnostic.py, or
enstrophy_exponent.py.

Run:
  python layer4/aposteriori_verification.py --ic tg --N 32 --n-steps 200 --verbose
  python layer4/aposteriori_verification.py --ic shear --N 32 --n-steps 200 --verbose
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from datetime import datetime

import numpy as np

# ── path setup (matches layer4/alignment_bridge.py, depletion_diagnostic.py) ──
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

import route2_3D as r2  # noqa: E402 — Route 2 3D solver (reused, not reimplemented)

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_SEED = 42
_IC_TAGS = {"tg": "TG", "shear": "SH"}
_DB_MAX_RETRIES = 6         # sqlite "database is locked" retry policy
_DB_BACKOFF_BASE_S = 0.15
_A_ZERO_TOL = 1e-14         # |A| below this treated as exactly 0 in the Gronwall step


# =============================================================================
# Residual of the momentum equation (spatial operator via route2_3D, ∂_t via
# time-differencing — see module docstring, "residual scope" §3)
# =============================================================================

def ns_residual(
    u: np.ndarray, v: np.ndarray, w: np.ndarray,
    p_hat_or_None,
    nu: float,
    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray,
    dudt_estimate: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """r = dudt_estimate − [−(u·∇)u − ∇p + ν·Δu], evaluated at (u, v, w).

    Spatial operator via layer3.route2_3D._rhs_velocity, IMPORTED verbatim
    (not reimplemented): nonlinear advection, dealiased, Leray-projected
    (implicitly removes ∇p for a divergence-free field — hence
    p_hat_or_None is accepted but not consumed), plus spectral ν·Δu.
    dudt_estimate is (dudt_u, dudt_v, dudt_w), typically from
    `estimate_dudt_time_diff`. Immutable. Returns (rx, ry, rz).
    """
    if not (u.shape == v.shape == w.shape):
        raise ValueError(f"u, v, w must share a shape; got {u.shape}, {v.shape}, {w.shape}")
    if u.ndim != 3 or u.shape[0] != u.shape[1] or u.shape[1] != u.shape[2]:
        raise ValueError(f"u, v, w must be N x N x N; got shape {u.shape}")

    N = u.shape[0]
    k2 = kx**2 + ky**2 + kz**2
    dealias = r2.make_dealias_mask_3D(N)

    rhs_u, rhs_v, rhs_w = r2._rhs_velocity(u, v, w, kx, ky, kz, k2, dealias, nu, N)

    dudt_u, dudt_v, dudt_w = dudt_estimate
    rx = dudt_u - rhs_u
    ry = dudt_v - rhs_v
    rz = dudt_w - rhs_w
    return rx, ry, rz


def estimate_dudt_time_diff(
    u_prev: np.ndarray, v_prev: np.ndarray, w_prev: np.ndarray,
    u_next: np.ndarray, v_next: np.ndarray, w_next: np.ndarray,
    dt: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """∂_t u_a ≈ (u_next − u_prev)/dt — PRIMARY (a) definition (module
    docstring): the discrete trajectory's own time-derivative, so the
    residual also carries the scheme's time-integration error. Immutable.
    """
    if dt <= 0.0:
        raise ValueError(f"dt must be > 0, got {dt}")
    return (
        (u_next - u_prev) / dt,
        (v_next - v_prev) / dt,
        (w_next - w_prev) / dt,
    )


# =============================================================================
# Residual norms: L² and a spectral H^{-1}
# =============================================================================

def residual_norms(rx: np.ndarray, ry: np.ndarray, rz: np.ndarray, dx: float) -> dict:
    """‖r‖_{L²} (volume-mean RMS) and ‖r‖_{H^{-1}} ((1+|k|²)^{-1} spectral
    energy weight) for a residual vector field. H^{-1} feeds
    `verification_margin` (the energy estimate for e=u−u_a pairs the
    residual, integrated by parts once, against ‖e‖_{L²}); weight ≤ 1
    everywhere so H^{-1} ≤ L² always, both zero iff r≡0. Wavenumbers come
    from route2_3D.make_grid_3D(N), N=rx.shape[0]; dx validated == 2π/N.
    Immutable. Raises ValueError on shape mismatch or inconsistent dx.
    """
    if not (rx.shape == ry.shape == rz.shape):
        raise ValueError(f"rx, ry, rz must share a shape; got {rx.shape}, {ry.shape}, {rz.shape}")
    if rx.ndim != 3 or rx.shape[0] != rx.shape[1] or rx.shape[1] != rx.shape[2]:
        raise ValueError(f"residual fields must be N x N x N; got shape {rx.shape}")

    N = rx.shape[0]
    expected_dx = 2.0 * np.pi / N
    if abs(dx - expected_dx) > 1e-9 * max(1.0, expected_dx):
        raise ValueError(
            f"dx={dx} inconsistent with N={N} periodic grid (expected dx=2π/N={expected_dx})"
        )

    l2 = float(np.sqrt(np.mean(rx**2 + ry**2 + rz**2)))

    _, _, _, k2 = r2.make_grid_3D(N)
    weight = 1.0 / (1.0 + k2)

    def _hm1_energy(f: np.ndarray) -> float:
        f_hat = np.fft.fftn(f)
        return float(np.sum(weight * np.abs(f_hat) ** 2))

    n6 = float(N) ** 6   # Parseval normalization for numpy's unnormalized fftn
    hm1_sq = (_hm1_energy(rx) + _hm1_energy(ry) + _hm1_energy(rz)) / n6
    hm1 = float(np.sqrt(max(hm1_sq, 0.0)))

    return {"L2": l2, "Hm1": hm1}


# =============================================================================
# Gronwall amplification rate A(t) — see SCOPE note #1 in module docstring
# =============================================================================

def amplification_rate(
    u: np.ndarray, v: np.ndarray, w: np.ndarray,
    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray,
    C: float = 1.0,
) -> float:
    """A(t) = C·‖∇u_a(t)‖_{L²} (RMS over all 9 gradient components) — the
    chosen amplification rate for dE/dt ≤ A(t)E(t) + ‖r(t)‖. C is NOT
    rigorously determined (SCOPE note #1). Immutable.
    """
    u_hat, v_hat, w_hat = np.fft.fftn(u), np.fft.fftn(v), np.fft.fftn(w)

    def _d(f_hat: np.ndarray, k_dir: np.ndarray) -> np.ndarray:
        return np.real(np.fft.ifftn(1j * k_dir * f_hat))

    grad_sq = (
        _d(u_hat, kx) ** 2 + _d(u_hat, ky) ** 2 + _d(u_hat, kz) ** 2
        + _d(v_hat, kx) ** 2 + _d(v_hat, ky) ** 2 + _d(v_hat, kz) ** 2
        + _d(w_hat, kx) ** 2 + _d(w_hat, ky) ** 2 + _d(w_hat, kz) ** 2
    )
    grad_l2 = float(np.sqrt(np.mean(grad_sq)))
    return C * grad_l2


def default_threshold(nu: float, C_thresh: float = 1.0) -> float:
    """ν-proportional scale for `verification_margin`'s threshold — a
    chosen surrogate for "error no longer controlled", NOT a rigorously
    derived local-existence radius (SCOPE note #2). C_thresh=1.0 default.
    """
    if nu <= 0.0:
        raise ValueError(f"nu must be > 0, got {nu}")
    return C_thresh * nu


# =============================================================================
# Gronwall / energy-inequality accumulation
# =============================================================================

def verification_margin(
    times: np.ndarray,
    residual_series: np.ndarray,
    A_series: np.ndarray,
    threshold: float,
) -> dict:
    """Accumulate E(t) via the exact frozen-coefficient (exponential-Euler)
    Gronwall recursion on dE/dt = A(t)·E(t) + r(t), E(0) = 0:

        E_{i+1} = E_i·exp(A_i·dt) + r_i·(exp(A_i·dt) − 1)/A_i   (A_i ≠ 0)
        E_{i+1} = E_i + r_i·dt                                  (A_i ≈ 0)

    with A_i, r_i the values at t_i, dt = t_{i+1} − t_i. Reduces EXACTLY to
    the closed-form Gronwall solution when A, r are constant over the
    series — E(T)=r0·T (A≡0), E(T)=r0·(e^{aT}−1)/a (A≡a>0) — the
    correctness anchor tested in test_aposteriori_verification.py.

    margin(t) = E(t)/threshold (see `default_threshold`, SCOPE note #2);
    margin<1 means "diagnostic satisfied on [0,t]", NOT "regularity proved".

    times, residual_series, A_series: equal-length arrays, len ≥ 2,
    strictly increasing times. threshold must be > 0. Returns dict with
    E_t, E_final, margin_t, margin_final. Immutable. Raises ValueError on
    any degenerate configuration (see checks below).
    """
    times = np.asarray(times, dtype=float)
    residual_series = np.asarray(residual_series, dtype=float)
    A_series = np.asarray(A_series, dtype=float)

    if times.size < 2:
        raise ValueError(
            f"verification_margin requires at least 2 time points to integrate; got {times.size}"
        )
    if not (times.size == residual_series.size == A_series.size):
        raise ValueError(
            "verification_margin: times, residual_series, A_series must have equal "
            f"length; got {times.size}, {residual_series.size}, {A_series.size}"
        )
    if threshold <= 0.0:
        raise ValueError(f"verification_margin: threshold must be > 0, got {threshold}")
    if np.any(np.diff(times) <= 0.0):
        raise ValueError("verification_margin: times must be strictly increasing")

    n = times.size
    E = np.zeros(n, dtype=float)
    for i in range(1, n):
        dt = times[i] - times[i - 1]
        a_i = A_series[i - 1]
        r_i = residual_series[i - 1]
        if abs(a_i) < _A_ZERO_TOL:
            E[i] = E[i - 1] + r_i * dt
        else:
            growth = np.exp(a_i * dt)
            E[i] = E[i - 1] * growth + r_i * (growth - 1.0) / a_i

    margin_t = E / threshold
    return {
        "E_t": E,
        "E_final": float(E[-1]),
        "margin_t": margin_t,
        "margin_final": float(margin_t[-1]),
    }


# =============================================================================
# Windowed (sliding / restarting) verification — see module SCOPE note on the
# known vacuity of the long-window margin above. Each window is an
# INDEPENDENT verification problem: E is reset to 0 at the window's start,
# so the Gronwall amplification factor is only exp(∫ over the window), not
# over the whole run. Reuses `verification_margin`'s own frozen-coefficient
# recursion for the accumulation on each window's sub-series — no second
# quadrature is introduced.
# =============================================================================

def windowed_verification(
    times: np.ndarray,
    residual_norms: np.ndarray,
    A_values: np.ndarray,
    window_length: float,
    threshold: float,
) -> list[dict]:
    """Partition [times[0], times[-1]] into consecutive windows of length
    `window_length`, starting at times[0]. On EACH window independently,
    reset the Gronwall energy to E=0 at the window's start and accumulate
    via `verification_margin` on that window's own (times, residual, A)
    sub-series — i.e. E_window = verification_margin(sub_times, sub_r,
    sub_A, threshold)["E_final"], the SAME exact frozen-coefficient
    (exponential-Euler) recursion used for the long-window margin, just
    restarted at each window boundary. This is the fix for the long-window
    margin's known vacuity (see module docstring): the amplification factor
    accrued per window is only exp(∫_{t_start}^{t_end} A ds), not
    exp(∫_0^T A ds).

    Partitioning: windows start at t_start = times[0] + k·window_length for
    k = 0, 1, 2, ... while t_start is still strictly less than times[-1].
    Each window's nominal end t_start + window_length is clipped to
    times[-1]. CONSEQUENCE: if (times[-1] - times[0]) is not an exact
    multiple of window_length, the LAST window is a "partial trailing
    window" shorter than window_length — it is returned as-is (its actual,
    shorter duration), never dropped and never padded to full length. This
    matches the intent of "verify each short window on the data actually
    available" rather than silently discarding trailing data.

    Consecutive windows share their boundary sample point (the last point
    of window k is the first point of window k+1); this is intentional —
    each window's accumulation still independently resets E=0 there.

    A window with fewer than 2 recorded sample points (i.e. window_length
    is smaller than the recording cadence) cannot be integrated by the
    frozen-coefficient recursion (which requires >= 2 points) and is
    skipped, not returned. If skipping empties the result entirely,
    ValueError is raised (a window_length that can never be evaluated is a
    degenerate configuration, not a valid diagnostic run).

    gronwall_factor = exp(∫_{t_start}^{t_end} A ds) for that window,
    evaluated with the SAME piecewise-constant, left-endpoint,
    frozen-coefficient rule as the homogeneous part of
    `verification_margin`'s recursion: exp(Σ_i A_i·dt_i) over the window's
    own subintervals (mathematically the product of each subinterval's
    exp(A_i·dt_i) growth factor, since exp(Σ) = Π exp).

    times, residual_norms, A_values: equal-length arrays, len >= 2, strictly
    increasing times (same shape contract as `verification_margin`).
    window_length, threshold: both must be > 0.

    Returns a list of dicts ordered by t_start:
    {t_start, t_end, gronwall_factor, E_window, margin_window, satisfied}
    where satisfied = (margin_window < 1.0). Raises ValueError on any
    degenerate configuration: too few time points, mismatched lengths,
    non-increasing times, non-positive window_length or threshold, or an
    empty resulting window list.
    """
    times = np.asarray(times, dtype=float)
    residual_norms = np.asarray(residual_norms, dtype=float)
    A_values = np.asarray(A_values, dtype=float)

    if times.size < 2:
        raise ValueError(
            f"windowed_verification requires at least 2 time points to integrate; got {times.size}"
        )
    if not (times.size == residual_norms.size == A_values.size):
        raise ValueError(
            "windowed_verification: times, residual_norms, A_values must have equal "
            f"length; got {times.size}, {residual_norms.size}, {A_values.size}"
        )
    if window_length <= 0.0:
        raise ValueError(f"windowed_verification: window_length must be > 0, got {window_length}")
    if threshold <= 0.0:
        raise ValueError(f"windowed_verification: threshold must be > 0, got {threshold}")
    if np.any(np.diff(times) <= 0.0):
        raise ValueError("windowed_verification: times must be strictly increasing")

    t0 = float(times[0])
    t_final = float(times[-1])
    eps = 1e-9 * max(1.0, abs(t_final - t0))

    windows: list[dict] = []
    t_start = t0
    while t_start < t_final - eps:
        t_end = min(t_start + window_length, t_final)
        idx = np.nonzero((times >= t_start - eps) & (times <= t_end + eps))[0]
        if idx.size >= 2:
            sub_times = times[idx]
            sub_r = residual_norms[idx]
            sub_A = A_values[idx]
            margin_out = verification_margin(sub_times, sub_r, sub_A, threshold)
            dt = np.diff(sub_times)
            gronwall_factor = float(np.exp(np.sum(sub_A[:-1] * dt)))
            windows.append({
                "t_start": float(sub_times[0]),
                "t_end": float(sub_times[-1]),
                "gronwall_factor": gronwall_factor,
                "E_window": margin_out["E_final"],
                "margin_window": margin_out["margin_final"],
                "satisfied": bool(margin_out["margin_final"] < 1.0),
            })
        t_start += window_length

    if not windows:
        raise ValueError(
            "windowed_verification: no window contained >= 2 sample points — "
            f"window_length={window_length} is too small relative to the sampling cadence "
            "in `times`; increase window_length or the recording density"
        )
    return windows


def windowed_summary(windows: list[dict]) -> dict:
    """Aggregate a `windowed_verification` window list into one summary dict:
    {n_windows, n_satisfied, fraction_satisfied, max_gronwall_factor,
    max_margin, median_margin, worst_window_t_start}. worst_window_t_start
    is the t_start of the window with the LARGEST margin_window (the
    window furthest from / least comfortably within "satisfied").
    fraction_satisfied == 1.0 means every window independently satisfied
    the (still non-rigorous) diagnostic — see module SCOPE, this is still
    not a proof. Raises ValueError if `windows` is empty.
    """
    if not windows:
        raise ValueError("windowed_summary: windows must be a non-empty list")

    n_windows = len(windows)
    satisfied_flags = [bool(w["satisfied"]) for w in windows]
    n_satisfied = int(sum(satisfied_flags))
    margins = np.array([float(w["margin_window"]) for w in windows], dtype=float)
    gronwalls = np.array([float(w["gronwall_factor"]) for w in windows], dtype=float)
    worst_idx = int(np.argmax(margins))

    return {
        "n_windows": n_windows,
        "n_satisfied": n_satisfied,
        "fraction_satisfied": float(n_satisfied) / float(n_windows),
        "max_gronwall_factor": float(np.max(gronwalls)),
        "max_margin": float(np.max(margins)),
        "median_margin": float(np.median(margins)),
        "worst_window_t_start": float(windows[worst_idx]["t_start"]),
    }


# =============================================================================
# Production experiment: run the existing Route 2 3D solver, track the
# residual / Gronwall diagnostic every `record_every` steps
# =============================================================================

def run_verification_experiment(
    ic: str = "tg",
    N: int = 32,
    n_steps: int = 200,
    nu: float | None = None,
    seed: int = _DEFAULT_SEED,
    record_every: int = 10,
    C_gronwall: float = 1.0,
    C_thresh: float = 1.0,
    cfl_safety: float = 0.4,
    verbose: bool = False,
    db_path: str | None = None,
    exp_id_override: str | None = None,
    claim_id_override: str | None = None,
    window_lengths: tuple[float, ...] = (),
) -> dict:
    """Run layer3/route2_3D.py's own solver (make_solver, taylor_green_ic /
    shear_layer_ic, set_ic, compute_cfl_dt, solver_step — imported, none
    reimplemented) and track the a posteriori verification diagnostic.

    Each recorded step: ∂_t u_a via single-step time-difference (pre→post
    state, `estimate_dudt_time_diff`); spatial operator at the PRE-step
    state (`ns_residual`); residual norms (`residual_norms`); A(t) at the
    POST-step state (`amplification_rate`). At the end,
    `verification_margin` accumulates E(t), margin(t)=E(t)/threshold over
    the recorded H^{-1} residual series.

    Verdict ("diagnostic satisfied" semantics only, never "regularity
    proved" — see module SCOPE): PASS iff every recorded quantity is
    finite AND margin(t)<1 throughout; PARTIAL if finite but margin(t)≥1
    somewhere; FAIL if any quantity is non-finite.

    ic: 'tg' or 'shear'. N: grid points/side. n_steps: solver steps. nu:
    defaults to route2_3D._NU_DEFAULT. seed: logged for provenance (both
    ICs deterministic). record_every: diagnostic cadence. C_gronwall: the
    C in amplification_rate (SCOPE #1). C_thresh: see default_threshold
    (SCOPE #2). db_path: if given, log via log_result(). window_lengths: if
    non-empty, ALSO run `windowed_verification` + `windowed_summary` at each
    given window length on the recorded (times, residual_Hm1, A) series
    (post-processing the same trajectory — no re-simulation), using the same
    `threshold`; results land under result["windowed"][window_length] as
    {"windows": [...], "summary": {...}}. This is the non-vacuous
    counterpart to the long-window margin above (see module SCOPE note on
    vacuity). Empty tuple (default) preserves prior behavior exactly.

    Raises ValueError for unknown ic or non-positive n_steps.
    """
    if ic not in _IC_TAGS:
        raise ValueError(f"Unknown ic='{ic}'; use 'tg' or 'shear'.")
    if n_steps <= 0:
        raise ValueError(f"n_steps must be > 0, got {n_steps}")

    if nu is None:
        nu = r2._NU_DEFAULT
    threshold = default_threshold(nu, C_thresh=C_thresh)

    # Seed set and logged for reproducibility provenance (see alignment_bridge.py
    # convention) — 'tg'/'shear' ICs here are deterministic; no draw is used.
    _rng = np.random.default_rng(seed)

    solver = r2.make_solver(N=N, nu=nu, eps_param=0.0)
    ic_dict = r2.taylor_green_ic(N) if ic == "tg" else r2.shear_layer_ic(N)
    solver = r2.set_ic(solver, ic_dict)

    dx = 2.0 * np.pi / N
    kx, ky, kz = solver["kx"], solver["ky"], solver["kz"]

    record_steps = [0]
    times = [0.0]
    residual_L2 = [0.0]
    residual_Hm1 = [0.0]
    A_vals = [amplification_rate(solver["u"], solver["v"], solver["w"], kx, ky, kz, C=C_gronwall)]

    t0 = time.perf_counter()
    for step_idx in range(n_steps):
        dt = r2.compute_cfl_dt(solver, safety=cfl_safety)
        prev_u, prev_v, prev_w = solver["u"], solver["v"], solver["w"]
        new_solver, _diag = r2.solver_step(solver, dt)

        is_last = step_idx == n_steps - 1
        if step_idx % record_every == 0 or is_last:
            dudt = estimate_dudt_time_diff(
                prev_u, prev_v, prev_w,
                new_solver["u"], new_solver["v"], new_solver["w"], dt,
            )
            rx, ry, rz = ns_residual(prev_u, prev_v, prev_w, None, nu, kx, ky, kz, dudt)
            norms = residual_norms(rx, ry, rz, dx)
            a_val = amplification_rate(
                new_solver["u"], new_solver["v"], new_solver["w"], kx, ky, kz, C=C_gronwall
            )

            record_steps.append(step_idx + 1)
            times.append(new_solver["t"])
            residual_L2.append(norms["L2"])
            residual_Hm1.append(norms["Hm1"])
            A_vals.append(a_val)

            if verbose:
                print(
                    f"  step {step_idx:4d}  t={new_solver['t']:.4f}  "
                    f"L2={norms['L2']:.4e}  Hm1={norms['Hm1']:.4e}  A={a_val:.4e}"
                )

        solver = new_solver

    wall = time.perf_counter() - t0

    times_arr = np.array(times, dtype=float)
    residual_L2_arr = np.array(residual_L2, dtype=float)
    residual_Hm1_arr = np.array(residual_Hm1, dtype=float)
    A_arr = np.array(A_vals, dtype=float)

    margin = verification_margin(times_arr, residual_Hm1_arr, A_arr, threshold)

    all_finite = bool(
        np.all(np.isfinite(residual_L2_arr)) and np.all(np.isfinite(residual_Hm1_arr))
        and np.all(np.isfinite(A_arr)) and np.all(np.isfinite(margin["E_t"]))
        and np.all(np.isfinite(margin["margin_t"]))
    )
    margin_below_1 = bool(np.all(margin["margin_t"] < 1.0)) if all_finite else False
    verdict = "FAIL" if not all_finite else ("PASS" if margin_below_1 else "PARTIAL")

    tag = _IC_TAGS[ic]
    exp_id = exp_id_override or f"EXP-L4-VERIFY-{tag}-001"
    claim_id = claim_id_override or f"aposteriori-{ic}"

    timeseries = [
        {"step": s, "t": float(t), "residual_L2": float(l2v), "residual_Hm1": float(hm1v),
         "A_t": float(a), "E_t": float(e), "margin_t": float(m)}
        for s, t, l2v, hm1v, a, e, m in zip(
            record_steps, times, residual_L2, residual_Hm1, A_vals,
            margin["E_t"], margin["margin_t"],
        )
    ]

    margin_max = float(np.max(margin["margin_t"])) if margin["margin_t"].size else float("nan")

    windowed: dict = {}
    for w_len in window_lengths:
        w_windows = windowed_verification(times_arr, residual_Hm1_arr, A_arr, w_len, threshold)
        w_summary = windowed_summary(w_windows)
        windowed[w_len] = {"windows": w_windows, "summary": w_summary}
        if verbose:
            print(
                f"    windowed W={w_len:.3f}: n_windows={w_summary['n_windows']}  "
                f"fraction_satisfied={w_summary['fraction_satisfied']:.3f}  "
                f"max_gronwall_factor={w_summary['max_gronwall_factor']:.4e}  "
                f"max_margin={w_summary['max_margin']:.4e}  "
                f"median_margin={w_summary['median_margin']:.4e}"
            )

    result = {
        "exp_id": exp_id, "claim_id": claim_id,
        "timestamp": datetime.utcnow().isoformat(),
        "ic": ic, "N": N, "n_steps": n_steps, "seed": seed, "nu": nu,
        "C_gronwall": C_gronwall, "threshold": threshold, "t_final": solver["t"],
        "residual_L2_min": float(np.min(residual_L2_arr)),
        "residual_L2_max": float(np.max(residual_L2_arr)),
        "residual_Hm1_min": float(np.min(residual_Hm1_arr)),
        "residual_Hm1_max": float(np.max(residual_Hm1_arr)),
        "A_max": float(np.max(A_arr)) if A_arr.size else float("nan"),
        "E_final": float(margin["E_final"]),
        "margin_final": float(margin["margin_final"]),
        "margin_max": margin_max,
        "diagnostic_satisfied": margin_below_1,
        "verdict": verdict, "wall_time_s": wall,
        "key_metric": (
            f"E_final={margin['E_final']:.4e}  margin_final={margin['margin_final']:.4e}  "
            f"margin_max={margin_max:.4e}  threshold={threshold:.4e}  "
            f"(diagnostic satisfied, NOT a proof of regularity)"
        ),
        "timeseries": timeseries,
    }
    if windowed:
        result["windowed"] = windowed

    if verbose:
        print(f"\n  A posteriori verification diagnostic ({ic}, N={N}, n_steps={n_steps})")
        print(f"    threshold:     {threshold:.4e}")
        print(f"    E_final:       {margin['E_final']:.4e}")
        print(f"    margin_final:  {margin['margin_final']:.4e}")
        print(f"    margin_max:    {margin_max:.4e}")
        print(f"    verdict:       {verdict}  (diagnostic satisfied, NOT regularity proved)")

    if db_path is not None:
        log_result(db_path, result)

    return result


# =============================================================================
# Results logging (layer4 per-module table pattern; retry-on-locked backoff,
# see layer4/depletion_diagnostic.py — other agents may write results.db
# concurrently)
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path, timeout=30.0)
    conn.execute(
        "CREATE TABLE IF NOT EXISTS layer4_aposteriori_experiments ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, exp_id TEXT UNIQUE NOT NULL, "
        "claim_id TEXT NOT NULL, timestamp TEXT NOT NULL, ic TEXT, "
        "N INTEGER, n_steps INTEGER, seed INTEGER, nu REAL, "
        "C_gronwall REAL, threshold REAL, t_final REAL, "
        "residual_L2_min REAL, residual_L2_max REAL, "
        "residual_Hm1_min REAL, residual_Hm1_max REAL, A_max REAL, "
        "E_final REAL, margin_final REAL, margin_max REAL, "
        "diagnostic_satisfied INTEGER, verdict TEXT NOT NULL, "
        "wall_time_s REAL, key_metric TEXT)"
    )
    conn.execute(
        "CREATE TABLE IF NOT EXISTS layer4_aposteriori_timeseries ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, exp_id TEXT NOT NULL, "
        "step INTEGER, t REAL, residual_L2 REAL, residual_Hm1 REAL, "
        "A_t REAL, E_t REAL, margin_t REAL)"
    )
    # windowed (sliding/restarting) verification — the non-vacuous fix; see
    # module SCOPE note. window_length is an ADDITIONAL column beyond the
    # minimal {exp_id, window_index, t_start, t_end, gronwall_factor,
    # E_window, margin_window, satisfied} spec: without it, rows from the
    # same exp_id logged at multiple window lengths (as the production
    # reruns do — W in {0.5, 1.0, 2.0}) would collide on window_index.
    conn.execute(
        "CREATE TABLE IF NOT EXISTS layer4_aposteriori_windows ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, exp_id TEXT NOT NULL, "
        "window_length REAL, window_index INTEGER, "
        "t_start REAL, t_end REAL, gronwall_factor REAL, "
        "E_window REAL, margin_window REAL, satisfied INTEGER)"
    )
    conn.commit()
    return conn


def _write_result(db_path: str, result: dict) -> None:
    """Single-attempt DB write (summary UPSERT + append-only timeseries)."""
    conn = _ensure_db(db_path)
    try:
        cols = [
            "exp_id", "claim_id", "timestamp", "ic", "N", "n_steps", "seed",
            "nu", "C_gronwall", "threshold", "t_final",
            "residual_L2_min", "residual_L2_max",
            "residual_Hm1_min", "residual_Hm1_max", "A_max",
            "E_final", "margin_final", "margin_max",
            "diagnostic_satisfied", "verdict", "wall_time_s", "key_metric",
        ]
        values = tuple(
            int(result[c]) if c == "diagnostic_satisfied" else result.get(c)
            for c in cols
        )
        ph = ", ".join("?" for _ in cols)
        conn.execute(
            f"INSERT OR REPLACE INTO layer4_aposteriori_experiments ({', '.join(cols)}) "
            f"VALUES ({ph})",
            values,
        )
        for rec in result.get("timeseries", []):
            conn.execute(
                """
                INSERT INTO layer4_aposteriori_timeseries
                (exp_id, step, t, residual_L2, residual_Hm1, A_t, E_t, margin_t)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    result["exp_id"], rec["step"], rec["t"],
                    rec["residual_L2"], rec["residual_Hm1"],
                    rec["A_t"], rec["E_t"], rec["margin_t"],
                ),
            )
        for w_len, w_data in result.get("windowed", {}).items():
            for w_idx, w in enumerate(w_data["windows"]):
                conn.execute(
                    """
                    INSERT INTO layer4_aposteriori_windows
                    (exp_id, window_length, window_index, t_start, t_end,
                     gronwall_factor, E_window, margin_window, satisfied)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        result["exp_id"], float(w_len), w_idx,
                        w["t_start"], w["t_end"], w["gronwall_factor"],
                        w["E_window"], w["margin_window"], int(w["satisfied"]),
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
    """Append one a posteriori verification experiment (+ timeseries) to
    results.db. Retries with short exponential backoff on
    sqlite3.OperationalError ("database is locked") — other layer4
    experiment scripts may be writing to the same results.db concurrently.
    Re-raises immediately on any other error, or once retries are exhausted.
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


# =============================================================================
# CLI entry point
# =============================================================================

def main() -> None:  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="A posteriori regularity VERIFICATION DIAGNOSTIC "
                    "(Morosi-Pizzocchero shape, floating point — NOT a proof)"
    )
    parser.add_argument("--ic", choices=["tg", "shear"], default="tg")
    parser.add_argument("--N", type=int, default=32)
    parser.add_argument("--n-steps", type=int, default=200)
    parser.add_argument("--record-every", type=int, default=10)
    parser.add_argument("--seed", type=int, default=_DEFAULT_SEED)
    parser.add_argument("--nu", type=float, default=None)
    parser.add_argument("--c-gronwall", type=float, default=1.0)
    parser.add_argument("--c-thresh", type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--db", default=_DEFAULT_DB)
    args = parser.parse_args()

    r = run_verification_experiment(
        ic=args.ic, N=args.N, n_steps=args.n_steps,
        record_every=args.record_every, seed=args.seed, nu=args.nu,
        C_gronwall=args.c_gronwall, C_thresh=args.c_thresh,
        verbose=args.verbose, db_path=args.db,
    )
    print(f"\nVerdict (diagnostic satisfied, NOT a proof): {r['verdict']}  "
          f"key_metric: {r['key_metric']}")


if __name__ == "__main__":
    main()
