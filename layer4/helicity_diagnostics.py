"""
layer4/helicity_diagnostics.py — Does twist locally deplete stretching? (WP4b/4c)
=================================================================================
HANDOVER-S46-OPUS §WP4b.  Diagnostics for the twisted vortex ring of
`layer3/twisted_ring_IC.py`:

  H = int u.omega          total helicity
  h = u.omega              local helicity density
  tau = e.(curl e)         local TWIST DENSITY of the vortex-line field,
                           e = omega/|omega|
  alpha = e.S.e            local STRETCHING RATE, S = sym grad u

and the pointwise correlation between twist density and alpha on the set
{|omega| >= M/4}.

THE QUESTION, and why it is worth asking (the Wall-6 test).
S46 recorded the warning that Kerr's anti-parallel tubes are helicity-neutral
by symmetry, so GLOBAL helicity conservation is not by itself a defence against
blowup.  The claim worth chasing is LOCAL: twist density depleting the
stretching rate pointwise.  That is a LEMMA CANDIDATE, not a blowup proxy --
which is exactly why it is information-bearing where c_K was not (W6).  A
regular flow can exhibit local depletion or fail to, and either answer means
something.

  DISCRIMINATING STATEMENT (written before any run):
    On a regular flow both |tau| and alpha are finite, O(1)-measurable fields
    on {|omega| >= M/4}; if twist locally depletes stretching then
    corr(|tau|, alpha) is NEGATIVE there, post-transient and pre-reconnection,
    with a sign that persists as |iota| grows and across two resolutions --
    and a regular flow is equally capable of showing a positive or a null
    correlation, so either outcome carries information.

THE PARITY TRAP -- why the SIGNED correlation cannot be the measurement.
tau is a PSEUDOSCALAR: under a mirror reflection tau -> -tau, while alpha is a
true scalar and is unchanged.  If the IC at -iota is the mirror image of the IC
at +iota (it is, for this construction), then

    corr(tau, alpha)  is ODD in iota      -- it MUST flip sign with iota
    corr(|tau|, alpha) is EVEN in iota    -- it may be stable

So a signed correlation that flips sign across iota is not evidence against a
mechanism; it is a symmetry identity, and reporting it as "a correlation that
flips sign with iota is not a mechanism" would be reading the mirror symmetry
of the initial data.  The PARITY-CORRECT measurement -- and the only one whose
sign-stability across iota is meaningful -- is corr(|tau|, alpha).  Both are
computed; the odd/even structure of the pair is used as a correctness check on
the IC family (test_helicity_diagnostics.py asserts it).

THE CORE-NESS CONFOUND.  |tau| and alpha are both elevated in the tube core,
so a raw correlation on {|omega| >= M/4} can report "depletion" that is really
"the core is different from the sheath".  `partial_corr_abs` residualises both
fields on a quadratic in |omega| before correlating, which removes the
leading form of that confound.  It, not the raw correlation, is the primary
number.

TIME-SPLITTING (HANDOVER §WP4b failure mode 2, S43).  The IC is a Biot-Savart
inversion, not an NS solution, so the first turnover carries an adjustment
transient; earlier sessions were misled by reading early-transient spikes as
convergence behaviour.  Every correlation is reported per time window, and a
post-transient window is defined explicitly by |d ln M / dt| falling below half
its initial value.

RECONNECTION (WP4c).  A correlation measured across a reconnection event is
measuring two different flows.  Topology change is flagged by (i) the number of
connected components of {|omega| >= theta_rec * M} under PERIODIC connectivity,
(ii) the minimum distance between distinct components passing through a local
minimum, and (iii) a jump in H relative to its viscous baseline.  The
correlation windows are then cut at the first flagged time.

MINIMUM POINT COUNTS, stated rather than assumed (S49 rule; results.db carries
a legacy PASS on a 2-point regression with r^2 = 1.0000):
    MIN_CORR_POINTS = 200   -- no correlation reported on fewer masked points
    MIN_FIT_POINTS  = 8     -- no regression or slope reported on fewer points
    MIN_WINDOW_SAMPLES = 4  -- no per-window correlation on fewer snapshots

House rules: spectral (FFT) derivatives only; float64; every logged result
carries exp_id + claim_id + key_metric + verdict.  Route 2 is the exact Prize
system (nu constant), so there is no layer1/tfirst_props.py lookup to make.
"""

from __future__ import annotations

import json
import os
import sqlite3
import sys
import time
from datetime import datetime

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

from layer4.alignment_bridge import compute_direction_field  # noqa: E402

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")

MIN_CORR_POINTS = 200
MIN_FIT_POINTS = 8
MIN_WINDOW_SAMPLES = 4
OMEGA_THRESHOLD = 0.25          # the {|omega| >= M/4} set of the handover
RECONNECT_OMEGA_LEVEL = 0.5     # level set used for component counting
MIN_TWIST_SIGNAL = 1.0e-2       # mean|tau| on the mask below which the twist
                                # field carries no signal and the correlation
                                # is NOT reported.  See `twist_alpha_correlation`.
N_SURROGATES = 200              # phase-randomised surrogates per null test
BOX_VOL = (2.0 * np.pi) ** 3


# =============================================================================
# Spectral field operators
# =============================================================================

def spectral_wavenumbers_3d(N: int):
    k1 = np.fft.fftfreq(N, d=1.0 / N).astype(np.float64)
    return k1[:, None, None], k1[None, :, None], k1[None, None, :]


def _grad(f, kx, ky, kz):
    fh = np.fft.fftn(f)
    return (np.real(np.fft.ifftn(1j * kx * fh)),
            np.real(np.fft.ifftn(1j * ky * fh)),
            np.real(np.fft.ifftn(1j * kz * fh)))


def vorticity(u, v, w, kx, ky, kz):
    """omega = curl u, spectrally."""
    ux, uy, uz = _grad(u, kx, ky, kz)
    vx, vy, vz = _grad(v, kx, ky, kz)
    wx_, wy_, wz_ = _grad(w, kx, ky, kz)
    return wy_ - vz, uz - wx_, vx - uy


def velocity_gradient(u, v, w, kx, ky, kz):
    """The 9 components of grad u, as a nested tuple J[i][j] = d_j u_i."""
    return (_grad(u, kx, ky, kz), _grad(v, kx, ky, kz), _grad(w, kx, ky, kz))


def strain_from_gradient(J):
    """S_ij = (d_j u_i + d_i u_j)/2 from J[i][j].

    S_ij is DOWNSTREAM of u (CLAUDE.md naming table) -- never a primary
    variable.
    """
    return [[0.5 * (J[i][j] + J[j][i]) for j in range(3)] for i in range(3)]


def stretching_rate(u, v, w, wx, wy, wz, kx, ky, kz):
    """alpha = e.S.e, e = omega/|omega|.  Returns (alpha, valid_mask)."""
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    J = velocity_gradient(u, v, w, kx, ky, kz)
    S = strain_from_gradient(J)
    e = (ex, ey, ez)
    alpha = np.zeros_like(mag)
    for i in range(3):
        for j in range(3):
            alpha = alpha + e[i] * S[i][j] * e[j]
    return np.where(valid, alpha, 0.0), valid, mag


def twist_density(wx, wy, wz, kx, ky, kz, floor_frac: float = 1e-6):
    """tau = e.(curl e), the local twist rate of the vortex-line field,
    computed through the EXACT identity

        e.(curl e)  ==  omega.(curl omega) / |omega|^2 ,      e = omega/|omega|

    which follows because curl(omega/|omega|) = curl(omega)/|omega|
    + grad(1/|omega|) x omega, and the second term is killed by the dot with
    e (a triple product carrying omega twice).

    WHY THE IDENTITY AND NOT THE DEFINITION.  Computing curl(e) by FFT
    directly, from an e set to zero outside the valid vorticity set, was the
    first implementation here and it is WRONG: e is a UNIT vector right up to
    the edge of the vorticity support and 0 outside it, so it carries a jump
    of size 1, and the spectral derivative of that jump rings across the whole
    box with O(1) amplitude.  On an UNTWISTED ring -- where tau must vanish
    identically -- that version returned tau = +/-0.74 against a true-signal
    scale of about 1.5, i.e. the contamination was the same order as the
    measurement.  It was caught by
    `test_untwisted_ring_has_near_zero_twist_density`, which is why that
    control exists.  The identity differentiates omega, which is smooth, so
    there is no truncation and no ringing.

    tau is the local twist rate of vortex lines about each other -- the
    gauge-free, pointwise form of the twist term in the Calugareanu-White
    decomposition, and the local field whose prescribed value along the tube
    core is the rotational transform iota.  It is a PSEUDOSCALAR.

    Returns (tau, valid) with tau set to 0 where |omega| < floor_frac * M.
    """
    mag2 = wx ** 2 + wy ** 2 + wz ** 2
    M2 = float(np.max(mag2))
    if M2 <= 0:
        return np.zeros_like(wx), np.zeros_like(wx, dtype=bool)
    valid = mag2 >= (floor_frac ** 2) * M2
    cx, cy, cz = _curl_of_vector(wx, wy, wz, kx, ky, kz)
    num = wx * cx + wy * cy + wz * cz
    den = np.where(valid, mag2, 1.0)
    return np.where(valid, num / den, 0.0), valid


def _curl_of_vector(ax, ay, az, kx, ky, kz):
    axx, axy, axz = _grad(ax, kx, ky, kz)
    ayx, ayy, ayz = _grad(ay, kx, ky, kz)
    azx, azy, azz = _grad(az, kx, ky, kz)
    return azy - ayz, axz - azx, ayx - axy


def current_helicity_density(wx, wy, wz, kx, ky, kz):
    """omega.(curl omega) -- the numerator of tau, recorded for context."""
    cx, cy, cz = _curl_of_vector(wx, wy, wz, kx, ky, kz)
    return wx * cx + wy * cy + wz * cz


def helicity_density(u, v, w, wx, wy, wz):
    return u * wx + v * wy + w * wz


def helicity_total(u, v, w, wx, wy, wz):
    """H = int_box u.omega dx, on [0,2pi)^3."""
    return float(np.mean(helicity_density(u, v, w, wx, wy, wz)) * BOX_VOL)


# =============================================================================
# Correlations — the measurement
# =============================================================================

def _pearson(a, b):
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    if a.size < 2:
        return float("nan")
    sa, sb = a.std(), b.std()
    if sa <= 0 or sb <= 0:
        return float("nan")
    return float(np.mean((a - a.mean()) * (b - b.mean())) / (sa * sb))


def _residualise(y, ctrl):
    """Residual of y after linear least squares on [1, c, c^2]."""
    c = np.asarray(ctrl, float).ravel()
    X = np.column_stack([np.ones_like(c), c, c ** 2])
    beta, *_ = np.linalg.lstsq(X, np.asarray(y, float).ravel(), rcond=None)
    return np.asarray(y, float).ravel() - X @ beta


def twist_alpha_correlation(tau, alpha, mag, valid,
                            threshold: float = OMEGA_THRESHOLD) -> dict:
    """Pointwise correlation of twist density with stretching on {|omega|>=tM}.

    Returns the SIGNED correlation (odd in iota by parity -- a symmetry
    identity, not a mechanism), the ABSOLUTE correlation corr(|tau|, alpha)
    (parity-even, the meaningful one), and the PARTIAL absolute correlation
    with |omega| residualised out (the primary number: it removes the
    core-ness confound).
    """
    M = float(np.max(mag)) if mag.size else 0.0
    if M <= 0:
        return {"n_points": 0, "insufficient": True, "M": 0.0}
    mask = valid & (mag >= threshold * M)
    n = int(np.count_nonzero(mask))
    tau_abs_mean = float(np.mean(np.abs(tau[mask]))) if n else 0.0
    # DEGENERATE-TWIST GUARD (S50).  A Pearson correlation is scale-invariant,
    # so it CANNOT report "there is no signal": handed a tau field that is pure
    # numerical residue it returns an arbitrary number.  That is exactly what
    # the untwisted (iota = 0) ring does -- mean|tau| there is ~2e-3 at N=64 and
    # ~1e-3 at N=96, i.e. converging to zero, while the iota >= 0.5 signal is
    # 0.72-1.95 and agrees to four significant figures between those two grids.
    # Read as a null control, iota = 0 gave +0.15 at N=64 and -0.42 at N=96 --
    # a sign flip under refinement that is a property of correlating noise, not
    # of the flow.  So iota = 0 is NOT a valid null control for this statistic
    # and is refused here rather than reported.  The valid null is
    # `surrogate_null_correlation` below.
    degenerate = tau_abs_mean < MIN_TWIST_SIGNAL
    out = {"n_points": n, "M": M, "threshold": threshold,
           "mask_volume_fraction": n / mag.size,
           "tau_abs_mean_on_mask": tau_abs_mean,
           "degenerate_twist": bool(degenerate),
           "insufficient": (n < MIN_CORR_POINTS) or bool(degenerate)}
    if out["insufficient"]:
        return out
    t = tau[mask]
    al = alpha[mask]
    om = mag[mask]
    out.update({
        "corr_signed": _pearson(t, al),
        "corr_abs": _pearson(np.abs(t), al),
        "partial_corr_abs": _pearson(_residualise(np.abs(t), om),
                                     _residualise(al, om)),
        "tau_abs_mean": float(np.mean(np.abs(t))),
        "alpha_mean": float(np.mean(al)),
        "alpha_std": float(np.std(al)),
        "alpha_positive_fraction": float(np.mean(al > 0)),
    })
    return out


def phase_randomise(field, rng):
    """A surrogate with the SAME power spectrum (hence the same spatial
    autocorrelation) as `field`, but randomised phases.

    Phases are taken from the FFT of a real Gaussian field, so they carry the
    Hermitian symmetry a real inverse transform needs.
    """
    F = np.fft.fftn(field)
    G = np.fft.fftn(rng.standard_normal(field.shape))
    ph = G / np.maximum(np.abs(G), 1e-300)
    return np.real(np.fft.ifftn(np.abs(F) * ph))


def surrogate_null_correlation(tau, alpha, mag, valid,
                               threshold: float = OMEGA_THRESHOLD,
                               n_surrogates: int = N_SURROGATES,
                               seed: int = 42) -> dict:
    """The valid null control for corr(|tau|, alpha).

    WHY A SURROGATE AND NOT THE iota = 0 RUN.  An untwisted ring has no twist
    signal at all, so correlating its |tau| against alpha correlates noise and
    returns an arbitrary value (see the guard in `twist_alpha_correlation`).
    The question a null control must answer is different: given fields with
    THESE marginals and THIS spatial autocorrelation, how large a correlation
    arises by chance?  Phase randomisation answers exactly that -- it preserves
    the power spectrum of |tau| (so the surrogate is as smooth and as
    spatially correlated as the real field) while destroying its pointwise
    alignment with alpha.

    A plain random permutation would NOT do: it destroys the autocorrelation
    too and so reports a null far narrower than the truth for smooth fields.

    Returns the null mean and standard deviation, the two-sided 95% band, and
    the z-score of the observed correlation against that null.
    """
    M = float(np.max(mag)) if mag.size else 0.0
    if M <= 0:
        return {"insufficient": True}
    mask = valid & (mag >= threshold * M)
    n = int(np.count_nonzero(mask))
    if n < MIN_CORR_POINTS:
        return {"insufficient": True, "n_points": n}
    at = np.abs(tau)
    obs = _pearson(at[mask], alpha[mask])
    rng = np.random.default_rng(seed)
    null = np.empty(n_surrogates)
    for i in range(n_surrogates):
        null[i] = _pearson(np.abs(phase_randomise(at, rng))[mask], alpha[mask])
    null = null[np.isfinite(null)]
    if null.size < MIN_FIT_POINTS:
        return {"insufficient": True, "n_surrogates_valid": int(null.size)}
    mu, sd = float(np.mean(null)), float(np.std(null))
    return {
        "insufficient": False,
        "n_points": n,
        "n_surrogates": int(null.size),
        "observed": obs,
        "null_mean": mu,
        "null_std": sd,
        "null_q025": float(np.percentile(null, 2.5)),
        "null_q975": float(np.percentile(null, 97.5)),
        "z_score": float((obs - mu) / sd) if sd > 0 else float("nan"),
        "outside_null_95": bool(obs < np.percentile(null, 2.5)
                                or obs > np.percentile(null, 97.5)),
    }


def helicity_snapshot(u, v, w, kx, ky, kz,
                      threshold: float = OMEGA_THRESHOLD) -> dict:
    """All WP4b quantities for one velocity snapshot."""
    wx, wy, wz = vorticity(u, v, w, kx, ky, kz)
    alpha, valid, mag = stretching_rate(u, v, w, wx, wy, wz, kx, ky, kz)
    tau, _ = twist_density(wx, wy, wz, kx, ky, kz)
    h = helicity_density(u, v, w, wx, wy, wz)
    corr = twist_alpha_correlation(tau, alpha, mag, valid, threshold)
    N = u.shape[0]
    kk = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
    E = sum(np.abs(np.fft.fftn(c) / N ** 3) ** 2 for c in (u, v, w))
    tot = float(E.sum())
    rec = component_topology(mag, RECONNECT_OMEGA_LEVEL)
    return {
        "H": float(np.mean(h) * BOX_VOL),
        "h_absmax": float(np.max(np.abs(h))),
        "M": float(np.max(mag)),
        "enstrophy": float(np.mean(mag ** 2) * BOX_VOL),
        "E": float(0.5 * np.mean(u ** 2 + v ** 2 + w ** 2) * BOX_VOL),
        "tail_fraction": float(E[kk > (2.0 / 3.0) * (N / 2)].sum() / tot) if tot > 0 else 0.0,
        "n_components": rec["n_components"],
        "min_component_distance": rec["min_distance"],
        **corr,
    }


# =============================================================================
# WP4c — reconnection / topology-change detection
# =============================================================================

def periodic_label(mask):
    """Connected components of a boolean mask with PERIODIC connectivity.

    scipy.ndimage.label does not wrap; components touching opposite faces are
    merged here by union-find over the wrapped face pairs.
    """
    from scipy.ndimage import label
    lab, n = label(mask)
    if n <= 1:
        return lab, int(n)
    parent = list(range(n + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)

    for axis in range(3):
        lo = np.take(lab, 0, axis=axis)
        hi = np.take(lab, -1, axis=axis)
        both = (lo > 0) & (hi > 0)
        for A, B in zip(lo[both].ravel(), hi[both].ravel()):
            union(int(A), int(B))

    roots = {}
    out = np.zeros_like(lab)
    nxt = 0
    for i in range(1, n + 1):
        r = find(i)
        if r not in roots:
            nxt += 1
            roots[r] = nxt
        out[lab == i] = roots[r]
    return out, nxt


def component_topology(mag, level: float = RECONNECT_OMEGA_LEVEL):
    """Number of components of {|omega| >= level*M} and their minimum
    separation (periodic minimum-image), as the topology observable."""
    M = float(np.max(mag)) if mag.size else 0.0
    if M <= 0:
        return {"n_components": 0, "min_distance": float("nan")}
    lab, n = periodic_label(mag >= level * M)
    if n <= 1:
        return {"n_components": int(n), "min_distance": float("nan")}
    N = mag.shape[0]
    dxg = 2.0 * np.pi / N
    cents = []
    for i in range(1, n + 1):
        idx = np.argwhere(lab == i)
        # sub-sample large components; centroid suffices for a distance TREND
        if idx.shape[0] > 4000:
            idx = idx[:: idx.shape[0] // 4000]
        cents.append(idx * dxg)
    best = float("inf")
    for i in range(n):
        for j in range(i + 1, n):
            A, B = cents[i], cents[j]
            if A.shape[0] * B.shape[0] > 4_000_000:
                A = A[::4]
                B = B[::4]
            d = A[:, None, :] - B[None, :, :]
            d = np.abs(d)
            d = np.minimum(d, 2.0 * np.pi - d)
            best = min(best, float(np.sqrt((d ** 2).sum(-1)).min()))
    return {"n_components": int(n), "min_distance": best}


def detect_reconnection(series, h_jump_rel: float = 0.05,
                        h_scale_floor: float = 1e-3) -> dict:
    """Flag the first topology change in a run.

    Three independent flags:
      (a) n_components changes from its initial value;
      (b) min_component_distance passes through a local minimum (approach then
          separation -- the signature of a reconnection event);
      (c) |H(t) - H(0)| / |H(0)| exceeds h_jump_rel over one sampling interval
          (a JUMP, not the slow viscous drift).
    Returns the earliest flagged time, or None.
    """
    ts = np.array([s["t"] for s in series], float)
    nc = np.array([s["n_components"] for s in series], float)
    md = np.array([s["min_component_distance"] for s in series], float)
    H = np.array([s["H"] for s in series], float)
    flags = {}

    if nc.size and not np.all(nc == nc[0]):
        i = int(np.argmax(nc != nc[0]))
        flags["component_count_change"] = float(ts[i])

    fin = np.isfinite(md)
    if fin.sum() >= 3:
        tt, mm = ts[fin], md[fin]
        for i in range(1, len(mm) - 1):
            if mm[i] < mm[i - 1] and mm[i] < mm[i + 1]:
                flags["core_distance_minimum"] = float(tt[i])
                break

    # The relative-jump test is meaningless when H(0) is ~0 by symmetry -- the
    # UNTWISTED ring is exactly that case, and an earlier version of this
    # function flagged every iota = 0 run as a reconnection because a 1e-18
    # helicity moved by orders of magnitude in relative terms.  Scale H by the
    # field's own helicity capacity 2*E*M rather than by H(0) when H(0) is
    # negligible against it.
    if H.size >= 2:
        cap = np.array([2.0 * s.get("E", 0.0) * s.get("M", 0.0) for s in series],
                       float)
        ref = max(abs(H[0]), h_scale_floor * float(np.max(cap)) if cap.size else 0.0)
        if ref > 0:
            d = np.abs(np.diff(H)) / ref
            j = np.where(d > h_jump_rel)[0]
            if j.size:
                flags["helicity_jump"] = float(ts[j[0] + 1])

    t_rec = min(flags.values()) if flags else None
    return {"flags": flags, "t_reconnect": t_rec,
            "reconnection_detected": t_rec is not None}


# =============================================================================
# Time-splitting
# =============================================================================

def post_transient_time(series, frac: float = 0.5) -> float:
    """First t at which |d ln M / dt| falls below `frac` times its initial value.

    The IC is a Biot-Savart inversion, not an NS solution, so the first
    turnover is an adjustment transient.  HANDOVER §WP4b failure mode 2 (S43):
    do not read the correlation inside it.
    """
    ts = np.array([s["t"] for s in series], float)
    M = np.array([s["M"] for s in series], float)
    if ts.size < 3 or np.any(M <= 0):
        return float(ts[0]) if ts.size else 0.0
    rate = np.abs(np.gradient(np.log(M), ts))
    r0 = float(np.max(rate[: max(2, ts.size // 5)]))
    if r0 <= 0:
        return float(ts[0])
    below = np.where(rate < frac * r0)[0]
    return float(ts[below[0]]) if below.size else float(ts[0])


def window_correlations(series, t_lo, t_hi, label) -> dict:
    """Aggregate the per-snapshot correlations over a time window."""
    sel = [s for s in series if t_lo <= s["t"] <= t_hi
           and not s.get("insufficient", True)]
    if len(sel) < MIN_WINDOW_SAMPLES:
        return {"window": label, "t_lo": t_lo, "t_hi": t_hi,
                "n_samples": len(sel), "insufficient": True,
                "reason": f"fewer than MIN_WINDOW_SAMPLES={MIN_WINDOW_SAMPLES}"}
    out = {"window": label, "t_lo": t_lo, "t_hi": t_hi,
           "n_samples": len(sel), "insufficient": False}
    for key in ("corr_signed", "corr_abs", "partial_corr_abs"):
        vals = np.array([s[key] for s in sel], float)
        vals = vals[np.isfinite(vals)]
        if vals.size < MIN_WINDOW_SAMPLES:
            out[key + "_median"] = float("nan")
            out[key + "_iqr"] = float("nan")
            continue
        out[key + "_median"] = float(np.median(vals))
        out[key + "_iqr"] = float(np.percentile(vals, 75) - np.percentile(vals, 25))
        out[key + "_sign_consistency"] = float(
            max(np.mean(vals > 0), np.mean(vals < 0)))
    out["n_points_median"] = float(np.median([s["n_points"] for s in sel]))
    return out


# =============================================================================
# Experiment runner — drives layer3/route2_3D_jax.py, sweeping iota
# =============================================================================

def run_twisted_ring(N: int = 64, iota: float = 1.0, nu: float = 1.0e-3,
                     t_end: float = 4.0, max_steps: int = 800,
                     sample_every: int = 8, R0: float | None = None,
                     a: float | None = None, circulation: float = 1.0,
                     seed: int = 42, cfl_safety: float = 0.4,
                     verbose: bool = False) -> dict:
    """Evolve one twisted ring with the existing Route-2 JAX solver (float64)
    and record the WP4b/4c diagnostics on a time series of snapshots."""
    os.environ.setdefault("JAX_ENABLE_X64", "1")
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    import route2_3D_jax as r2j  # noqa: E402
    from twisted_ring_IC import twisted_ring_ic, DEFAULT_R0, DEFAULT_A  # noqa: E402

    np.random.seed(seed)
    t0 = time.perf_counter()
    R0 = DEFAULT_R0 if R0 is None else R0
    a = DEFAULT_A if a is None else a

    ic = twisted_ring_ic(N, iota=iota, R0=R0, a=a, circulation=circulation,
                         seed=seed)
    kx, ky, kz = spectral_wavenumbers_3d(N)

    solver = r2j.make_solver_jax(N=N, nu=nu, eps_param=0.1, dt_max=0.05,
                                 dtype=jnp.float64)
    solver = r2j.set_ic_jax(solver, ic)
    u, v, w, th = solver["u"], solver["v"], solver["w"], solver["theta"]
    nuj = jnp.array(nu, dtype=jnp.float64)
    t = 0.0
    series = []

    snap = helicity_snapshot(np.asarray(u), np.asarray(v), np.asarray(w),
                             kx, ky, kz)
    snap.update({"t": 0.0, "step": 0})
    series.append(snap)

    for step in range(1, max_steps + 1):
        if t >= t_end:
            break
        dt = r2j.compute_cfl_dt_jax(u, v, w, N, nu, solver["dt_max"], cfl_safety)
        dt = min(dt, t_end - t)
        u, v, w, th = r2j.solver_step_jax(
            u, v, w, th, solver["kx"], solver["ky"], solver["kz"],
            solver["k2"], solver["dealias"], solver["nq_mask"], nuj,
            jnp.array(dt, dtype=jnp.float64))
        t += dt
        if step % sample_every == 0 or t >= t_end:
            snap = helicity_snapshot(np.asarray(u), np.asarray(v),
                                     np.asarray(w), kx, ky, kz)
            snap.update({"t": t, "step": step})
            series.append(snap)
            if verbose:
                print(f"   step {step:4d} t={t:.3f} M={snap['M']:.3e} "
                      f"H={snap['H']:+.4e} "
                      f"c_abs={snap.get('corr_abs', float('nan')):+.4f} "
                      f"c_part={snap.get('partial_corr_abs', float('nan')):+.4f} "
                      f"tail={snap['tail_fraction']:.2e}")

    rec = detect_reconnection(series)
    t_trans = post_transient_time(series)
    t_fin = series[-1]["t"]
    t_cut = rec["t_reconnect"] if rec["t_reconnect"] is not None else t_fin

    thirds = [
        window_correlations(series, 0.0, t_fin / 3.0, "early"),
        window_correlations(series, t_fin / 3.0, 2.0 * t_fin / 3.0, "middle"),
        window_correlations(series, 2.0 * t_fin / 3.0, t_fin, "late"),
    ]
    primary = window_correlations(series, t_trans, t_cut,
                                  "post_transient_pre_reconnection")

    return {
        "N": N, "iota": float(iota), "nu": nu, "seed": seed,
        "R0": R0, "a": a, "circulation": circulation,
        "t_final": t_fin, "n_steps": series[-1]["step"],
        "n_samples": len(series),
        "helicity_ic": ic["helicity_ic"],
        "omega_div_defect": ic["omega_div_defect"],
        "core_cells": ic["resolution"]["core_cells"],
        "H0": series[0]["H"], "H_final": series[-1]["H"],
        "H_rel_change": (abs(series[-1]["H"] - series[0]["H"]) / abs(series[0]["H"])
                         if abs(series[0]["H"]) > 0 else float("nan")),
        "M_max": float(max(s["M"] for s in series)),
        "tail_fraction_max": float(max(s["tail_fraction"] for s in series)),
        "t_transient": t_trans, "reconnection": rec,
        "windows": thirds, "primary": primary,
        "series": series,
        "wall_time_s": time.perf_counter() - t0,
    }


def iota_sweep(iotas=(-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0), N: int = 64,
               nu: float = 1.0e-3, t_end: float = 4.0, verbose: bool = True,
               **kw) -> dict:
    """Sweep the rotational transform.  Includes +/- pairs deliberately: the
    parity identity (signed correlation odd, absolute correlation even) is the
    correctness check on the IC family.

    COMMON WINDOW.  Each run's own `primary` window starts at its own
    `post_transient_time`, and those differ across iota (they came out 0.0 for
    |iota| = 2 and 0.8 elsewhere in the first S50 sweep).  Comparing
    correlations taken over DIFFERENT time ranges is not a comparison, so after
    the runs finish every run's `primary` is RECOMPUTED on one window common to
    the whole sweep: from the largest post-transient time over the sweep to the
    earliest reconnection (or the end).  The per-run adaptive window is kept as
    `primary_own_window` for reference.
    """
    runs = []
    for io in iotas:
        r = run_twisted_ring(N=N, iota=io, nu=nu, t_end=t_end, verbose=False, **kw)
        runs.append(r)
        if verbose:
            p = r["primary"]
            print(f"  iota={io:+5.2f} N={N}  "
                  f"c_signed={p.get('corr_signed_median', float('nan')):+.4f}  "
                  f"c_abs={p.get('corr_abs_median', float('nan')):+.4f}  "
                  f"c_partial={p.get('partial_corr_abs_median', float('nan')):+.4f}  "
                  f"nsamp={p.get('n_samples')}  "
                  f"recon={r['reconnection']['reconnection_detected']}  "
                  f"tail={r['tail_fraction_max']:.2e}")
    t_lo = max(r["t_transient"] for r in runs)
    rec_times = [r["reconnection"]["t_reconnect"] for r in runs
                 if r["reconnection"]["t_reconnect"] is not None]
    t_hi = min([min(r["t_final"] for r in runs)] + rec_times)
    for r in runs:
        r["primary_own_window"] = r["primary"]
        r["primary"] = window_correlations(r["series"], t_lo, t_hi,
                                           "sweep_common_window")
    if verbose:
        print(f"  [common window] t in [{t_lo:.3f}, {t_hi:.3f}]", flush=True)
        for r in runs:
            p = r["primary"]
            print(f"   iota={r['iota']:+5.2f}  "
                  f"c_abs={p.get('corr_abs_median', float('nan')):+.4f}  "
                  f"c_partial={p.get('partial_corr_abs_median', float('nan')):+.4f}  "
                  f"n={p.get('n_samples')}", flush=True)
    return {"N": N, "nu": nu, "t_end": t_end, "runs": runs,
            "common_window": [t_lo, t_hi],
            "stability": sweep_stability(runs)}


def sweep_stability(runs) -> dict:
    """Is the correlation's sign stable across iota?

    Reports, for each of the three correlation measures, the median over the
    non-degenerate |iota| > 0 runs and the fraction sharing the modal sign.
    The ABSOLUTE and PARTIAL measures are the ones for which sign-stability is
    meaningful; the SIGNED measure is odd in iota by parity and its sign-flip
    is recorded as a symmetry check, not as instability.
    """
    out = {}
    sel = [r for r in runs if abs(r["iota"]) > 1e-12
           and not r["primary"].get("insufficient", True)]
    out["n_runs_used"] = len(sel)
    if len(sel) < 2:
        out["insufficient"] = True
        return out
    out["insufficient"] = False
    for key in ("corr_signed", "corr_abs", "partial_corr_abs"):
        vals = np.array([r["primary"].get(key + "_median", np.nan) for r in sel])
        good = np.isfinite(vals)
        v = vals[good]
        if v.size < 2:
            continue
        out[key + "_median"] = float(np.median(v))
        out[key + "_min"] = float(np.min(v))
        out[key + "_max"] = float(np.max(v))
        out[key + "_sign_stable"] = bool(np.all(v > 0) or np.all(v < 0))
        out[key + "_modal_sign_fraction"] = float(max(np.mean(v > 0), np.mean(v < 0)))
    # parity check: is the signed correlation odd in iota?
    by_iota = {round(r["iota"], 6): r["primary"].get("corr_signed_median", np.nan)
               for r in sel}
    pairs = [(k, by_iota[k], by_iota[-k]) for k in by_iota
             if -k in by_iota and k > 0
             and np.isfinite(by_iota[k]) and np.isfinite(by_iota[-k])]
    if pairs:
        odd = [abs(a + b) / max(abs(a), abs(b), 1e-30) for _, a, b in pairs]
        out["signed_parity_oddness_max_rel_residual"] = float(np.max(odd))
        by_abs = {round(r["iota"], 6): r["primary"].get("corr_abs_median", np.nan)
                  for r in sel}
        ev = [abs(by_abs[k] - by_abs[-k]) / max(abs(by_abs[k]), abs(by_abs[-k]), 1e-30)
              for k, _, _ in pairs if np.isfinite(by_abs.get(k, np.nan))
              and np.isfinite(by_abs.get(-k, np.nan))]
        if ev:
            out["abs_parity_evenness_max_rel_residual"] = float(np.max(ev))
    return out


# =============================================================================
# results.db
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_helicity_experiments (
            id                         INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id                     TEXT UNIQUE NOT NULL,
            claim_id                   TEXT NOT NULL,
            timestamp                  TEXT NOT NULL,
            ic                         TEXT,
            iota                       REAL,
            N                          INTEGER,
            nu                         REAL,
            seed                       INTEGER,
            R0                         REAL,
            core_a                     REAL,
            core_cells                 REAL,
            omega_div_defect           REAL,
            n_steps                    INTEGER,
            n_samples                  INTEGER,
            t_final                    REAL,
            t_transient                REAL,
            reconnection_detected      INTEGER,
            t_reconnect                REAL,
            H0                         REAL,
            H_final                    REAL,
            H_rel_change               REAL,
            M_max                      REAL,
            tail_fraction_max          REAL,
            corr_signed_median         REAL,
            corr_abs_median            REAL,
            partial_corr_abs_median    REAL,
            corr_abs_sign_consistency  REAL,
            primary_window_samples     INTEGER,
            primary_n_points_median    REAL,
            verdict                    TEXT NOT NULL,
            wall_time_s                REAL,
            key_metric                 TEXT NOT NULL,
            notes                      TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_helicity_timeseries (
            id                     INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id                 TEXT NOT NULL,
            claim_id               TEXT NOT NULL,
            step                   INTEGER,
            t                      REAL,
            H                      REAL,
            M                      REAL,
            enstrophy              REAL,
            E                      REAL,
            corr_signed            REAL,
            corr_abs               REAL,
            partial_corr_abs       REAL,
            n_points               INTEGER,
            n_components           INTEGER,
            min_component_distance REAL,
            tail_fraction          REAL
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_helicity_sweep (
            id                       INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id                   TEXT UNIQUE NOT NULL,
            claim_id                 TEXT NOT NULL,
            timestamp                TEXT NOT NULL,
            N                        INTEGER,
            nu                       REAL,
            t_end                    REAL,
            iotas                    TEXT,
            n_runs_used              INTEGER,
            corr_abs_median          REAL,
            corr_abs_min             REAL,
            corr_abs_max             REAL,
            corr_abs_sign_stable     INTEGER,
            partial_corr_abs_median  REAL,
            partial_corr_abs_min     REAL,
            partial_corr_abs_max     REAL,
            partial_sign_stable      INTEGER,
            signed_parity_residual   REAL,
            abs_parity_residual      REAL,
            verdict                  TEXT NOT NULL,
            key_metric               TEXT NOT NULL,
            notes                    TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_helicity_null_control (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id           TEXT UNIQUE NOT NULL,
            claim_id         TEXT NOT NULL,
            timestamp        TEXT NOT NULL,
            N                INTEGER,
            iota             REAL,
            t_snapshot       REAL,
            tau_abs_mean     REAL,
            degenerate_twist INTEGER,
            n_points         INTEGER,
            observed         REAL,
            null_mean        REAL,
            null_std         REAL,
            null_q025        REAL,
            null_q975        REAL,
            z_score          REAL,
            outside_null_95  INTEGER,
            verdict          TEXT NOT NULL,
            key_metric       TEXT NOT NULL,
            notes            TEXT
        )
    """)
    conn.commit()
    return conn


def log_run(db_path: str, exp_id: str, claim_id: str, res: dict,
            verdict: str, notes: str = "") -> None:
    if not claim_id:
        raise ValueError("claim_id is mandatory — a result without one is not a result")
    conn = _ensure_db(db_path)
    p = res["primary"]
    key_metric = (
        f"partial_corr(|tau|,alpha)_median={p.get('partial_corr_abs_median', float('nan')):+.4f} "
        f"corr(|tau|,alpha)_median={p.get('corr_abs_median', float('nan')):+.4f} "
        f"n_win_samples={p.get('n_samples')} n_pts_median={p.get('n_points_median')}")
    conn.execute("""
        INSERT OR REPLACE INTO layer4_helicity_experiments
        (exp_id, claim_id, timestamp, ic, iota, N, nu, seed, R0, core_a,
         core_cells, omega_div_defect, n_steps, n_samples, t_final,
         t_transient, reconnection_detected, t_reconnect, H0, H_final,
         H_rel_change, M_max, tail_fraction_max, corr_signed_median,
         corr_abs_median, partial_corr_abs_median, corr_abs_sign_consistency,
         primary_window_samples, primary_n_points_median, verdict,
         wall_time_s, key_metric, notes)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, (exp_id, claim_id, datetime.now().isoformat(),
          f"twisted_ring_iota{res['iota']:+.3f}", res["iota"], res["N"],
          res["nu"], res["seed"], res["R0"], res["a"], res["core_cells"],
          res["omega_div_defect"], res["n_steps"], res["n_samples"],
          res["t_final"], res["t_transient"],
          int(res["reconnection"]["reconnection_detected"]),
          res["reconnection"]["t_reconnect"], res["H0"], res["H_final"],
          res["H_rel_change"], res["M_max"], res["tail_fraction_max"],
          p.get("corr_signed_median"), p.get("corr_abs_median"),
          p.get("partial_corr_abs_median"), p.get("corr_abs_sign_consistency"),
          p.get("n_samples"), p.get("n_points_median"), verdict,
          res["wall_time_s"], key_metric, notes))
    conn.execute("DELETE FROM layer4_helicity_timeseries WHERE exp_id = ?", (exp_id,))
    for s in res["series"]:
        conn.execute("""
            INSERT INTO layer4_helicity_timeseries
            (exp_id, claim_id, step, t, H, M, enstrophy, E, corr_signed,
             corr_abs, partial_corr_abs, n_points, n_components,
             min_component_distance, tail_fraction)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (exp_id, claim_id, s["step"], s["t"], s["H"], s["M"], s["enstrophy"], s["E"],
              s.get("corr_signed"), s.get("corr_abs"), s.get("partial_corr_abs"),
              s.get("n_points"), s["n_components"], s["min_component_distance"],
              s["tail_fraction"]))
    conn.commit()
    conn.close()


def log_null_control(db_path: str, exp_id: str, claim_id: str, row: dict,
                     notes: str = "") -> None:
    """Log one surrogate-null control point (WP4b, S50)."""
    if not claim_id:
        raise ValueError("claim_id is mandatory")
    conn = _ensure_db(db_path)
    km = (f"observed={row['observed']:+.4f} null={row['null_mean']:+.4f}"
          f"+/-{row['null_std']:.4f} z={row['z']:+.2f} "
          f"mean|tau|={row['tau_abs_mean']:.5f} n_pts={row['n_points']}")
    verdict = ("DEGENERATE" if row["degenerate"]
               else ("SIGNIFICANT" if row["significant"] else "NULL"))
    conn.execute("""
        INSERT OR REPLACE INTO layer4_helicity_null_control
        (exp_id, claim_id, timestamp, N, iota, t_snapshot, tau_abs_mean,
         degenerate_twist, n_points, observed, null_mean, null_std, null_q025,
         null_q975, z_score, outside_null_95, verdict, key_metric, notes)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, (exp_id, claim_id, datetime.now().isoformat(), row["N"], row["iota"],
          row["t"], row["tau_abs_mean"], int(bool(row["degenerate"])),
          row["n_points"], row["observed"], row["null_mean"], row["null_std"],
          row["q025"], row["q975"], row["z"], int(bool(row["significant"])),
          verdict, km, notes))
    conn.commit(); conn.close()


def log_sweep(db_path: str, exp_id: str, claim_id: str, sweep: dict,
              verdict: str, notes: str = "") -> None:
    if not claim_id:
        raise ValueError("claim_id is mandatory")
    conn = _ensure_db(db_path)
    st = sweep["stability"]
    key_metric = (
        f"partial_corr(|tau|,alpha) median={st.get('partial_corr_abs_median', float('nan')):+.4f} "
        f"range=[{st.get('partial_corr_abs_min', float('nan')):+.4f},"
        f"{st.get('partial_corr_abs_max', float('nan')):+.4f}] "
        f"sign_stable={st.get('partial_corr_abs_sign_stable')} "
        f"n_runs={st.get('n_runs_used')}")
    conn.execute("""
        INSERT OR REPLACE INTO layer4_helicity_sweep
        (exp_id, claim_id, timestamp, N, nu, t_end, iotas, n_runs_used,
         corr_abs_median, corr_abs_min, corr_abs_max, corr_abs_sign_stable,
         partial_corr_abs_median, partial_corr_abs_min, partial_corr_abs_max,
         partial_sign_stable, signed_parity_residual, abs_parity_residual,
         verdict, key_metric, notes)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, (exp_id, claim_id, datetime.now().isoformat(), sweep["N"], sweep["nu"],
          sweep["t_end"], json.dumps([r["iota"] for r in sweep["runs"]]),
          st.get("n_runs_used"), st.get("corr_abs_median"),
          st.get("corr_abs_min"), st.get("corr_abs_max"),
          int(bool(st.get("corr_abs_sign_stable", False))),
          st.get("partial_corr_abs_median"), st.get("partial_corr_abs_min"),
          st.get("partial_corr_abs_max"),
          int(bool(st.get("partial_corr_abs_sign_stable", False))),
          st.get("signed_parity_oddness_max_rel_residual"),
          st.get("abs_parity_evenness_max_rel_residual"),
          verdict, key_metric, notes))
    conn.commit()
    conn.close()


__all__ = [
    "MIN_CORR_POINTS", "MIN_FIT_POINTS", "MIN_WINDOW_SAMPLES",
    "OMEGA_THRESHOLD", "RECONNECT_OMEGA_LEVEL",
    "spectral_wavenumbers_3d", "vorticity", "velocity_gradient",
    "strain_from_gradient", "stretching_rate", "twist_density",
    "helicity_density", "helicity_total", "current_helicity_density", "twist_alpha_correlation",
    "helicity_snapshot", "periodic_label", "component_topology",
    "detect_reconnection", "post_transient_time", "window_correlations",
    "phase_randomise", "surrogate_null_correlation",
    "MIN_TWIST_SIGNAL", "N_SURROGATES", "run_twisted_ring", "iota_sweep", "sweep_stability", "log_run", "log_sweep", "log_null_control",
    "make_helicity_figure",
]


# =============================================================================
# Figures
# =============================================================================

def make_helicity_figure(sweep: dict, out_path: str,
                         extra_runs: dict | None = None) -> str:
    """Correlation time series per iota, the iota dependence, and the
    resolution check."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:      # a missing plotting library must never lose a run
        return ""

    runs = sweep["runs"]
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.3))

    for r in runs:
        s = [x for x in r["series"] if not x.get("insufficient", True)]
        if not s:
            continue
        ax[0].plot([x["t"] for x in s], [x["partial_corr_abs"] for x in s],
                   marker=".", ms=3, label=f"$\\iota$={r['iota']:+.1f}")
        ax[0].axvline(r["t_transient"], color="grey", lw=0.4, alpha=0.4)
    ax[0].axhline(0.0, color="k", lw=0.8, ls="--")
    ax[0].set_xlabel("t")
    ax[0].set_ylabel(r"partial corr($|\tau|,\alpha$ | $|\omega|$)")
    ax[0].set_title("grey lines: end of adjustment transient")
    ax[0].legend(fontsize=7, ncol=2)

    io = [r["iota"] for r in runs]
    for key, mark, lab in (("corr_signed_median", "s", r"corr($\tau,\alpha$)  [odd by parity]"),
                           ("corr_abs_median", "o", r"corr($|\tau|,\alpha$)"),
                           ("partial_corr_abs_median", "^",
                            r"partial corr($|\tau|,\alpha$ | $|\omega|$)")):
        y = [r["primary"].get(key, np.nan) for r in runs]
        ax[1].plot(io, y, mark + "-", label=lab)
    ax[1].axhline(0.0, color="k", lw=0.8, ls="--")
    ax[1].set_xlabel(r"rotational transform $\iota$")
    ax[1].set_title("post-transient, pre-reconnection medians")
    ax[1].legend(fontsize=7)
    ax[1].grid(alpha=0.3)

    if extra_runs:
        for N, rl in sorted(extra_runs.items()):
            ax[2].plot([r["iota"] for r in rl],
                       [r["primary"].get("partial_corr_abs_median", np.nan)
                        for r in rl], "o-", label=f"N={N}")
    ax[2].plot(io, [r["primary"].get("partial_corr_abs_median", np.nan)
                    for r in runs], "o-", label=f"N={sweep['N']}")
    ax[2].axhline(0.0, color="k", lw=0.8, ls="--")
    ax[2].set_xlabel(r"$\iota$")
    ax[2].set_ylabel("partial corr")
    ax[2].set_title("resolution check")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path
