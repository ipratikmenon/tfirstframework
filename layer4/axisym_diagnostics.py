"""
layer4/axisym_diagnostics.py — Axisymmetric diagnostics (WP4a, HANDOVER-S46-OPUS §WP4)
=======================================================================================
WP4a asks for axisymmetric diagnostics on the EXISTING 3D periodic Route-2 JAX
solver (`layer3/route2_3D_jax.py`) — no new solver.  Three quantities:

  1. Gamma = r * u_theta     — obeys a maximum principle for axisymmetric NS.
  2. omega_theta / r         — the quantity whose transport equation carries the
                               whole difficulty of axisymmetric-with-swirl.
  3. Source  d_z(Gamma^2)/r^4 — localised in (r, z).  This single term is the
                               entire obstruction; everything else is benign.

The axisymmetric-with-swirl system (Hou-Li form), about the z axis:

    d_t Gamma  + u_r d_r Gamma  + u_z d_z Gamma
        = nu ( d_rr - (1/r) d_r + d_zz ) Gamma                          (G)

    d_t (w/r) + u_r d_r (w/r) + u_z d_z (w/r)
        = nu ( d_rr + (3/r) d_r + d_zz ) (w/r)  +  (1/r^4) d_z(Gamma^2) (W)

with w = omega_theta.  Equation (G) has NO zeroth-order term, so by the
parabolic maximum principle

        t -> max_x |Gamma(x,t)|   is NON-INCREASING.                    (MP)

WHY THIS MODULE IS WORTH RUNNING (the Wall-6 test, HANDOVER §WP4 preamble).
Before any experiment this program now requires a one-sentence statement of what
a REGULAR flow could show that would distinguish the competing hypotheses.  For
the three measurements here:

  (MP) Gamma maximum principle -- DISCRIMINATING, and not about physics at all:
       a regular axisymmetric flow MUST show max|Gamma| non-increasing, so any
       increase above the measured axisymmetry defect and the time-truncation
       floor is a SOLVER BUG.  The outcome is falsifiable on exactly the flows
       this solver can produce, which is what c_K never was.

  Source near-axis evaluability -- DISCRIMINATING as a gate on WP4d: for smooth
       axisymmetric flow u_theta ~ c*r near the axis, so
           Gamma^2/r^4 = (u_theta/r)^2 = eta^2
       is FINITE on the axis.  A regular flow therefore tells us whether the
       r^-4 form can be evaluated at achievable resolution or is dominated by
       1/r cancellation noise -- i.e. whether a true (r,z) solver (WP4d) would
       be measuring the term or measuring its own roundoff.  Either answer is
       information.

  omega_theta/r localisation -- DESCRIPTIVE, not a hypothesis test.  Recorded
       as context, and labelled as such; it is not offered as evidence.

NUMERICS NOTE -- the two forms of the source term.
Because r is independent of z, the following is an EXACT identity, not an
approximation:
        (1/r^4) d_z(Gamma^2)  ==  d_z( Gamma^2 / r^4 )  ==  d_z( eta^2 ),
        eta := u_theta / r.
The `eta` form is regular on the axis; the `r^-4` form divides two quantities
that both vanish to fourth order there.  Computing BOTH and reporting their
relative discrepancy as a function of r is the near-axis evaluability
measurement described above.

CALIBRATION (house rule: calibrate on an exact or planted solution first).
`burgers_vortex_fields` supplies the exact Burgers vortex, for which all three
diagnostics have closed forms:
        Gamma_B(r) = (Circ/2pi)(1 - exp(-a r^2 / 4nu))     (z-independent)
        omega_theta == 0,      source == 0 identically.
Because Burgers gives only the TRIVIAL value of the source, a second planted
field with a known NON-ZERO source is supplied by `planted_swirl_fields`:
        eta(r,z) = c exp(-r^2/s^2) sin(z)
          =>  source = d_z(eta^2) = c^2 exp(-2 r^2/s^2) sin(2z).

HONEST CAVEAT, STATED UP FRONT.  The Burgers vortex is NOT periodic: its strain
field u_r = -(a/2) r, u_z = a z is unbounded and cannot live on T^3.  It is
therefore NOT a steady state of this periodic solver and the solver cannot be
expected to hold it.  It is used here to calibrate the DIAGNOSTICS against
closed forms, which is what the house rule actually requires.  The solver's own
"holds an exact solution" test is done separately on a Beltrami (ABC) flow,
which IS an exact periodic NS solution decaying as exp(-nu t).

PROPERTY ENGINE.  Route 2 is the exact Prize system: nu is a constant numerical
parameter, not a thermophysical property, so there is no Layer-1 lookup to make
here.  No property formula is inlined -- there is none to inline.

House rules: spectral (FFT) derivatives only; float64 throughout (the Gamma
maximum-principle test is a cancellation test and is meaningless in float32);
every logged result carries exp_id + claim_id + key_metric + verdict.
"""

from __future__ import annotations

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

_DEFAULT_DB = os.path.join(_ROOT, "results", "results.db")
_DEFAULT_SEED = 42

# ── constants, stated rather than buried ─────────────────────────────────────
AXIS_EPS_CELLS = 0.5     # grid points with r < AXIS_EPS_CELLS*dx are masked
                         # from every r-divided quantity (the axis column).
MIN_FIT_POINTS = 8       # S49 rule: no regression is reported on fewer points
                         # than this.  (results.db carries a legacy PASS on a
                         # 2-point fit with r^2 = 1.0000; see S37 log.)
MIN_CORR_POINTS = 200    # minimum masked grid points for a reported correlation
MP_TOL_REL = 1e-9        # relative growth of max|Gamma| attributable to float64
                         # roundoff in a single step


# =============================================================================
# Grid and cylindrical decomposition
# =============================================================================

def spectral_wavenumbers_3d(N: int):
    """(kx, ky, kz) broadcast-ready for an N^3 grid on T^3 = [0,2pi)^3.

    Same convention as route2_3D.make_grid_3D / route2_3D_jax.make_grid_jax:
    np.fft.fftfreq(N, d=1/N), meshgrid indexing='ij'.
    """
    k1d = np.fft.fftfreq(N, d=1.0 / N).astype(np.float64)
    return k1d[:, None, None], k1d[None, :, None], k1d[None, None, :]


def cylindrical_grid(N: int, x0: float | None = None, y0: float | None = None):
    """Cylindrical radius r and azimuth phi about an axis parallel to z.

    The axis passes through (x0, y0); both default to the box centre pi.
    Returns (r, cos_phi, sin_phi, axis_mask) where axis_mask is True on the
    points EXCLUDED from r-divided quantities (r < AXIS_EPS_CELLS*dx).

    cos/sin phi are returned rather than phi itself: every use below is
    through cos/sin, and forming them directly from (dx, dy)/r avoids an
    arctan2 round trip.
    """
    if x0 is None:
        x0 = np.pi
    if y0 is None:
        y0 = np.pi
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    dxg = 2.0 * np.pi / N
    # minimum-image displacement from the axis
    ddx = (x - x0 + np.pi) % (2.0 * np.pi) - np.pi
    ddy = (x - y0 + np.pi) % (2.0 * np.pi) - np.pi
    DX = ddx[:, None, None] * np.ones((1, N, 1))
    DY = ddy[None, :, None] * np.ones((N, 1, 1))
    DX = np.broadcast_to(ddx[:, None, None], (N, N, N))
    DY = np.broadcast_to(ddy[None, :, None], (N, N, N))
    r = np.sqrt(DX ** 2 + DY ** 2)
    axis_mask = r < AXIS_EPS_CELLS * dxg
    r_safe = np.where(axis_mask, 1.0, r)
    cos_phi = np.where(axis_mask, 1.0, DX / r_safe)
    sin_phi = np.where(axis_mask, 0.0, DY / r_safe)
    return (np.ascontiguousarray(r), np.ascontiguousarray(cos_phi),
            np.ascontiguousarray(sin_phi), np.ascontiguousarray(axis_mask))


def to_cylindrical(u, v, w, cos_phi, sin_phi):
    """Cartesian (u,v,w) -> cylindrical (u_r, u_theta, u_z).  Pure."""
    u_r = u * cos_phi + v * sin_phi
    u_th = -u * sin_phi + v * cos_phi
    return u_r, u_th, w


def vorticity_spectral(u, v, w, kx, ky, kz):
    """omega = curl u, spectrally.  Does not mutate inputs."""
    uh, vh, wh = np.fft.fftn(u), np.fft.fftn(v), np.fft.fftn(w)
    d = lambda fh, k: np.real(np.fft.ifftn(1j * k * fh))
    wx = d(wh, ky) - d(vh, kz)
    wy = d(uh, kz) - d(wh, kx)
    wz = d(vh, kx) - d(uh, ky)
    return wx, wy, wz


def ddz_spectral(f, kz):
    """d_z f by FFT along the last axis (z is periodic; r, phi are not)."""
    return np.real(np.fft.ifft(1j * kz * np.fft.fft(f, axis=2), axis=2))


# =============================================================================
# The three WP4a quantities
# =============================================================================

def gamma_field(u, v, r, cos_phi, sin_phi):
    """Gamma = r * u_theta.  Regular everywhere including the axis."""
    _, u_th, _ = to_cylindrical(u, v, np.zeros_like(u), cos_phi, sin_phi)
    return r * u_th


def eta_field(u, v, r, cos_phi, sin_phi, axis_mask):
    """eta = u_theta / r = Gamma / r^2.  Finite on the axis for smooth flow
    (u_theta ~ c r), but 0/0 on the grid: the axis column is masked to 0 and
    reported separately via `axis_mask`."""
    _, u_th, _ = to_cylindrical(u, v, np.zeros_like(u), cos_phi, sin_phi)
    r_safe = np.where(axis_mask, 1.0, r)
    return np.where(axis_mask, 0.0, u_th / r_safe)


def omega_theta_over_r(wx, wy, r, cos_phi, sin_phi, axis_mask):
    """(omega_theta)/r, with the axis column masked to 0."""
    w_th = -wx * sin_phi + wy * cos_phi
    r_safe = np.where(axis_mask, 1.0, r)
    return np.where(axis_mask, 0.0, w_th / r_safe), w_th


def swirl_source_eta(u, v, r, cos_phi, sin_phi, axis_mask, kz):
    """Source d_z(Gamma^2)/r^4 in the axis-REGULAR form d_z(eta^2).

    Exact identity (r is z-independent), and the numerically stable branch.
    """
    eta = eta_field(u, v, r, cos_phi, sin_phi, axis_mask)
    return ddz_spectral(eta ** 2, kz), eta


def swirl_source_naive(u, v, r, cos_phi, sin_phi, axis_mask, kz):
    """Source in the literal form (1/r^4) d_z(Gamma^2).

    Kept deliberately naive: this is the branch whose near-axis failure is the
    thing being measured.  Axis column masked to 0 (it is 0/0 there).
    """
    G = gamma_field(u, v, r, cos_phi, sin_phi)
    dzG2 = ddz_spectral(G ** 2, kz)
    r_safe = np.where(axis_mask, 1.0, r)
    return np.where(axis_mask, 0.0, dzG2 / r_safe ** 4)


def source_near_axis_report(src_eta, src_naive, r, axis_mask, n_bins: int = 12):
    """Relative discrepancy between the two source forms, binned in r.

    This is the WP4d gate: if the naive r^-4 branch departs from the regular
    branch by O(1) inside the first few cells, a true (r,z) solver evaluating
    the term in that form would be measuring its own roundoff.

    Returns dict with r bin centres, per-bin relative L2 discrepancy, and the
    innermost bin at which the discrepancy first exceeds 1 (== total loss).
    """
    valid = ~axis_mask
    rv = r[valid]
    a = src_eta[valid]
    b = src_naive[valid]
    r_max = float(np.max(rv))
    edges = np.linspace(0.0, r_max, n_bins + 1)
    centres, rel = [], []
    for i in range(n_bins):
        m = (rv >= edges[i]) & (rv < edges[i + 1])
        if np.count_nonzero(m) < 1:
            continue
        num = float(np.sqrt(np.mean((a[m] - b[m]) ** 2)))
        den = float(np.sqrt(np.mean(a[m] ** 2)))
        centres.append(0.5 * (edges[i] + edges[i + 1]))
        rel.append(num / den if den > 0 else (0.0 if num == 0 else np.inf))
    centres = np.asarray(centres)
    rel = np.asarray(rel)
    bad = np.where(rel > 1.0)[0]
    return {
        "r_bin_centres": centres.tolist(),
        "rel_discrepancy": rel.tolist(),
        "r_loss_radius": float(centres[bad[-1]]) if bad.size else 0.0,
        "rel_discrepancy_max": float(np.max(rel)) if rel.size else 0.0,
        "rel_discrepancy_innermost": float(rel[0]) if rel.size else 0.0,
    }


def localise_rz(field, r, N, n_r: int = 16, n_z: int = 16, reduce: str = "absmax"):
    """Localise a 3D field on an (r, z) grid by azimuthal reduction.

    reduce='absmax' -> max |field| over the azimuth in each (r,z) bin
    reduce='mean'   -> azimuthal mean
    Returns (r_centres, z_centres, table[n_r, n_z]).
    """
    r_max = float(np.max(r))
    edges = np.linspace(0.0, r_max, n_r + 1)
    r_c = 0.5 * (edges[:-1] + edges[1:])
    z = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    z_idx = np.array_split(np.arange(N), n_z)
    z_c = np.array([float(np.mean(z[ii])) for ii in z_idx])
    out = np.zeros((n_r, n_z))
    for i in range(n_r):
        m_r = (r >= edges[i]) & (r < edges[i + 1])
        for j, ii in enumerate(z_idx):
            sub = field[:, :, ii][m_r[:, :, ii]]
            if sub.size == 0:
                continue
            out[i, j] = float(np.max(np.abs(sub))) if reduce == "absmax" else float(np.mean(sub))
    return r_c, z_c, out


# =============================================================================
# Axisymmetry defect — the caveat, measured rather than assumed
# =============================================================================

def polar_resample(field, N, n_r=24, n_theta=32, x0=None, y0=None, r_max=None):
    """Resample an N^3 Cartesian field onto a polar grid (r, theta, z).

    Cubic-spline interpolation with periodic wrap in x and y (the box is
    periodic, so 'grid-wrap' is the correct boundary rule).  Returns
    (r_centres, table[n_r, n_theta, N]).
    """
    from scipy.ndimage import map_coordinates
    if x0 is None:
        x0 = np.pi
    if y0 is None:
        y0 = np.pi
    if r_max is None:
        r_max = 0.9 * np.pi
    dxg = 2.0 * np.pi / N
    rr = np.linspace(r_max / n_r, r_max, n_r)
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    R, TH = np.meshgrid(rr, th, indexing="ij")
    xs = (x0 + R * np.cos(TH)) / dxg
    ys = (y0 + R * np.sin(TH)) / dxg
    out = np.empty((n_r, n_theta, N))
    for k in range(N):
        out[:, :, k] = map_coordinates(field[:, :, k], [xs, ys], order=3,
                                       mode="grid-wrap")
    return rr, out


def axisymmetry_defect(u, v, w, r, cos_phi, sin_phi, N,
                       n_r: int = 24, n_theta: int = 32, r_max=None):
    """Relative departure from axisymmetry, by azimuthal Fourier decomposition.

    An axisymmetric IC in a PERIODIC box is not exactly axisymmetric: the
    periodic images break the rotational symmetry.  The Gamma maximum
    principle holds only for genuinely axisymmetric flow, so the size of this
    defect is the tolerance against which any max|Gamma| growth must be judged.

        D = sqrt( sum_{m != 0} |f_m|^2 / sum_m |f_m|^2 )

    with f_m the azimuthal Fourier coefficients of the three CYLINDRICAL
    components on a polar resampling.  An earlier (r,z)-binning version of this
    function was discarded: binning conflates within-bin RADIAL variation with
    azimuthal variation and reported D ~ 0.19 on an analytically axisymmetric
    field.  The Fourier form has no such leak.
    """
    u_r, u_th, u_z = to_cylindrical(u, v, w, cos_phi, sin_phi)
    num = 0.0
    den = 0.0
    for comp in (u_r, u_th, u_z):
        _, tab = polar_resample(comp, N, n_r=n_r, n_theta=n_theta, r_max=r_max)
        fm = np.fft.fft(tab, axis=1) / n_theta
        p = np.abs(fm) ** 2
        num += float(np.sum(p[:, 1:, :]))
        den += float(np.sum(p))
    return float(np.sqrt(num / den)) if den > 0 else 0.0


# =============================================================================
# Exact / planted calibration fields
# =============================================================================

def burgers_vortex_fields(N: int, a: float = 1.0, nu: float = 1.0e-2,
                          circulation: float = 1.0,
                          x0: float | None = None, y0: float | None = None):
    """Exact Burgers vortex on the N^3 grid, plus its closed-form diagnostics.

        u_r = -(a/2) r,   u_z = a (z - z0),
        u_theta = (Circ / 2 pi r) (1 - exp(-a r^2 / 4 nu))

    Closed forms:
        Gamma_B(r)   = (Circ/2pi)(1 - exp(-a r^2/4nu))        (z-independent)
        omega_theta  == 0
        source       == 0     (Gamma is z-independent)
        omega_z(r)   = (Circ a / 4 pi nu) exp(-a r^2 / 4 nu)

    NOT periodic (the strain is unbounded): this is a DIAGNOSTIC calibration,
    not a solver steady state.  See the module docstring.
    """
    r, cos_phi, sin_phi, axis_mask = cylindrical_grid(N, x0, y0)
    z = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    zz = np.broadcast_to(z[None, None, :], (N, N, N)) - np.pi

    r_safe = np.where(axis_mask, 1.0, r)
    core = 1.0 - np.exp(-a * r ** 2 / (4.0 * nu))
    u_th = np.where(axis_mask, 0.0, circulation * core / (2.0 * np.pi * r_safe))
    u_r = -0.5 * a * r
    u_z = a * zz

    u = u_r * cos_phi - u_th * sin_phi
    v = u_r * sin_phi + u_th * cos_phi
    w = u_z

    exact = {
        "Gamma": circulation * core / (2.0 * np.pi),
        "omega_theta": np.zeros((N, N, N)),
        "source": np.zeros((N, N, N)),
        "omega_z": circulation * a / (4.0 * np.pi * nu) * np.exp(-a * r ** 2 / (4.0 * nu)),
    }
    grid = {"r": r, "cos_phi": cos_phi, "sin_phi": sin_phi, "axis_mask": axis_mask}
    return (np.ascontiguousarray(u), np.ascontiguousarray(v),
            np.ascontiguousarray(w)), exact, grid


def planted_swirl_fields(N: int, c: float = 1.0, s: float = 1.0,
                         x0: float | None = None, y0: float | None = None):
    """Planted field with a known NON-ZERO source (Burgers gives only zero).

        eta(r,z) = c exp(-r^2/s^2) sin(z),   u_theta = r eta
      => source = d_z(eta^2) = c^2 exp(-2 r^2/s^2) sin(2 z)

    Purely azimuthal velocity, so div u = 0 exactly and the field is periodic
    in z.  The radial profile is Gaussian, hence effectively periodic in (x,y)
    for s well inside the box.
    """
    r, cos_phi, sin_phi, axis_mask = cylindrical_grid(N, x0, y0)
    z = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    zz = np.broadcast_to(z[None, None, :], (N, N, N))
    eta = c * np.exp(-r ** 2 / s ** 2) * np.sin(zz)
    u_th = r * eta
    u = -u_th * sin_phi
    v = u_th * cos_phi
    w = np.zeros((N, N, N))
    exact = {
        "eta": eta,
        "Gamma": r ** 2 * eta,
        "source": c ** 2 * np.exp(-2.0 * r ** 2 / s ** 2) * np.sin(2.0 * zz),
    }
    grid = {"r": r, "cos_phi": cos_phi, "sin_phi": sin_phi, "axis_mask": axis_mask}
    return (np.ascontiguousarray(u), np.ascontiguousarray(v),
            np.ascontiguousarray(w)), exact, grid


def abc_flow(N: int, A: float = 1.0, B: float = 1.0, C: float = 1.0):
    """ABC (Arnold-Beltrami-Childress) flow: an EXACT periodic NS solution.

    It is Beltrami (omega = u), so (u.grad)u = grad(|u|^2/2) is absorbed into
    the pressure and the flow decays exactly as u(t) = u(0) exp(-nu t),
    E(t) = E(0) exp(-2 nu t).  This is the solver's "holds an exact solution"
    calibration in the geometry the periodic solver actually has -- the Burgers
    vortex cannot play that role because it is not periodic.
    """
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return (A * np.sin(Z) + C * np.cos(Y),
            B * np.sin(X) + A * np.cos(Z),
            C * np.sin(Y) + B * np.cos(X))


# =============================================================================
# Snapshot diagnostic
# =============================================================================

def spectral_tail_fraction(u, v, w, kx, ky, kz, N: int, frac: float = 2.0 / 3.0) -> float:
    """Fraction of kinetic energy above frac*k_max — the resolution monitor.

    2/3 is the dealiasing cutoff: energy there is what the 2/3 rule is about to
    truncate, so a rising tail fraction is loss of resolution.  Recorded on
    every snapshot because the Gamma maximum-principle test turned out to be
    sensitive to exactly this (WP4a, S50).
    """
    kk = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
    E = np.zeros_like(kk)
    for c in (u, v, w):
        E = E + np.abs(np.fft.fftn(c) / N ** 3) ** 2
    tot = float(E.sum())
    if tot <= 0.0:
        return 0.0
    return float(E[kk > frac * (N / 2)].sum() / tot)


def axisym_snapshot(u, v, w, grid, kx, ky, kz, N: int) -> dict:
    """All WP4a quantities for one velocity snapshot."""
    r, cp, sp, am = grid["r"], grid["cos_phi"], grid["sin_phi"], grid["axis_mask"]
    G = gamma_field(u, v, r, cp, sp)
    wx, wy, wz = vorticity_spectral(u, v, w, kx, ky, kz)
    wtr, w_th = omega_theta_over_r(wx, wy, r, cp, sp, am)
    src_eta, eta = swirl_source_eta(u, v, r, cp, sp, am, kz)
    src_naive = swirl_source_naive(u, v, r, cp, sp, am, kz)
    near = source_near_axis_report(src_eta, src_naive, r, am)
    omega_mag = np.sqrt(wx ** 2 + wy ** 2 + wz ** 2)
    return {
        "Gamma_absmax": float(np.max(np.abs(G))),
        "Gamma_min": float(np.min(G)),
        "Gamma_max": float(np.max(G)),
        "omega_theta_over_r_absmax": float(np.max(np.abs(wtr))),
        "omega_theta_absmax": float(np.max(np.abs(w_th))),
        "source_absmax": float(np.max(np.abs(src_eta))),
        "source_naive_absmax": float(np.max(np.abs(src_naive))),
        "source_rel_discrepancy_innermost": near["rel_discrepancy_innermost"],
        "source_rel_discrepancy_max": near["rel_discrepancy_max"],
        "source_r_loss_radius": near["r_loss_radius"],
        "eta_absmax": float(np.max(np.abs(eta))),
        "M": float(np.max(omega_mag)),
        "axisym_defect": axisymmetry_defect(u, v, w, r, cp, sp, N),
        "tail_fraction": spectral_tail_fraction(u, v, w, kx, ky, kz, N),
        "r_at_Gamma_max": float(r[np.unravel_index(np.argmax(np.abs(G)), G.shape)]),
        "Gamma_absmax_core": float(np.max(np.abs(G[r <= 2.0]))),
        "_near_axis": near,
    }


def gamma_max_principle_check(gamma_series, defect_series=None,
                              tol_rel: float = MP_TOL_REL) -> dict:
    """Test (MP): is t -> max|Gamma| non-increasing?

    Returns the worst RELATIVE increase over the series,
        growth = max_n ( G_n - min_{m<=n} G_m ) / G_0,
    together with a verdict.  A violation larger than the axisymmetry defect
    (which measures how far the periodic-box flow is from the axisymmetric
    class where (MP) is a theorem) and larger than `tol_rel` is a SOLVER BUG,
    not physics -- HANDOVER §WP4a.
    """
    g = np.asarray(gamma_series, dtype=float)
    if g.size < 2:
        return {"verdict": "INSUFFICIENT", "n_points": int(g.size)}
    g0 = float(g[0])
    running_min = np.minimum.accumulate(g)
    rise = g - running_min
    worst_i = int(np.argmax(rise))
    growth_rel = float(rise[worst_i] / g0) if g0 > 0 else float("inf")
    monotone_rel = float((np.max(g) - g0) / g0) if g0 > 0 else float("inf")
    defect = float(np.max(defect_series)) if defect_series is not None and len(defect_series) else 0.0
    tol = max(tol_rel, defect)
    return {
        "n_points": int(g.size),
        "Gamma0": g0,
        "Gamma_absmax_over_run": float(np.max(g)),
        "Gamma_final": float(g[-1]),
        "max_relative_rise": growth_rel,
        "net_relative_rise": monotone_rel,
        "rise_index": worst_i,
        "axisym_defect_max": defect,
        "tolerance_used": tol,
        "violated": bool(growth_rel > tol),
        "verdict": "FAIL" if growth_rel > tol else "PASS",
    }


# =============================================================================
# Axisymmetric IC for the periodic solver
# =============================================================================

def axisym_swirl_ic(N: int, psi0: float = 1.0, gamma0: float = 1.0,
                    s: float = 1.0, z_mod: float = 0.5,
                    x0: float | None = None, y0: float | None = None) -> dict:
    """Axisymmetric IC with swirl, EXACTLY divergence-free and periodic in z.

    Poloidal part from a Stokes stream function; swirl from Gamma directly:

        psi(r,z)   = psi0  r^2 exp(-r^2/s^2) sin z
        u_r = -(1/r) d_z psi = -psi0 r exp(-r^2/s^2) cos z
        u_z =  (1/r) d_r psi =  psi0 (2 - 2 r^2/s^2) exp(-r^2/s^2) sin z
        Gamma(r,z) = gamma0 r^2 exp(-r^2/s^2) (1 + z_mod sin z)
        u_theta    = Gamma / r

    div u = 0 identically for any stream function in this form, so the solver's
    Leray projection is the identity on it to roundoff -- itself a check.
    Gamma depends on z, so the source d_z(Gamma^2)/r^4 is NON-trivial: the
    thing WP4a is meant to localise is actually present.

    Every component carries a factor r or r^2 near the axis, as smoothness
    requires, so eta = u_theta/r is finite there.
    """
    r, cos_phi, sin_phi, axis_mask = cylindrical_grid(N, x0, y0)
    z = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    zz = np.broadcast_to(z[None, None, :], (N, N, N))
    g = np.exp(-r ** 2 / s ** 2)

    u_r = -psi0 * r * g * np.cos(zz)
    u_z = psi0 * (2.0 - 2.0 * r ** 2 / s ** 2) * g * np.sin(zz)
    u_th = gamma0 * r * g * (1.0 + z_mod * np.sin(zz))

    u = u_r * cos_phi - u_th * sin_phi
    v = u_r * sin_phi + u_th * cos_phi
    return {
        "u": np.ascontiguousarray(u),
        "v": np.ascontiguousarray(v),
        "w": np.ascontiguousarray(u_z),
        "name": "axisym_swirl",
        "grid": {"r": r, "cos_phi": cos_phi, "sin_phi": sin_phi,
                 "axis_mask": axis_mask},
    }


# =============================================================================
# Experiment runner — drives layer3/route2_3D_jax.py (no new solver)
# =============================================================================

def run_axisym_experiment(
    N: int = 64,
    nu: float = 1.0e-3,
    t_end: float = 1.0,
    max_steps: int = 400,
    sample_every: int = 5,
    psi0: float = 1.0,
    gamma0: float = 1.0,
    s: float = 1.0,
    seed: int = _DEFAULT_SEED,
    cfl_safety: float = 0.4,
    verbose: bool = False,
) -> dict:
    """Evolve an axisymmetric swirling IC with the Route-2 JAX solver and track
    the WP4a diagnostics, in float64.

    float64 is mandatory here: (MP) is a cancellation test, and in float32 the
    roundoff floor (~1e-7 relative) sits above any physically meaningful
    Gamma drift, which would make the maximum-principle check vacuous.
    """
    os.environ.setdefault("JAX_ENABLE_X64", "1")
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    import route2_3D_jax as r2j  # noqa: E402

    np.random.seed(seed)
    t_wall0 = time.perf_counter()

    ic = axisym_swirl_ic(N, psi0=psi0, gamma0=gamma0, s=s)
    grid = ic["grid"]
    kx, ky, kz = spectral_wavenumbers_3d(N)

    solver = r2j.make_solver_jax(N=N, nu=nu, eps_param=0.1,
                                 dt_max=0.05, dtype=jnp.float64)
    solver = r2j.set_ic_jax(solver, ic)

    # Leray-projection residual: the IC is analytically divergence-free, so
    # this measures how much the projection had to change it.
    proj_res = float(np.sqrt(np.mean(
        (np.asarray(solver["u"]) - ic["u"]) ** 2 +
        (np.asarray(solver["v"]) - ic["v"]) ** 2 +
        (np.asarray(solver["w"]) - ic["w"]) ** 2)))
    ic_rms = float(np.sqrt(np.mean(ic["u"] ** 2 + ic["v"] ** 2 + ic["w"] ** 2)))
    proj_res_rel = proj_res / ic_rms if ic_rms > 0 else 0.0

    jkx, jky, jkz = solver["kx"], solver["ky"], solver["kz"]
    k2, dealias, nq = solver["k2"], solver["dealias"], solver["nq_mask"]
    nu_j = jnp.array(nu, dtype=jnp.float64)

    u, v, w, theta = solver["u"], solver["v"], solver["w"], solver["theta"]
    t = 0.0
    series: list[dict] = []

    snap0 = axisym_snapshot(np.asarray(u), np.asarray(v), np.asarray(w),
                            grid, kx, ky, kz, N)
    snap0.update({"t": 0.0, "step": 0, "dt": 0.0})
    series.append(snap0)

    for step in range(1, max_steps + 1):
        if t >= t_end:
            break
        dt = r2j.compute_cfl_dt_jax(u, v, w, N, nu, solver["dt_max"], cfl_safety)
        dt = min(dt, t_end - t)
        u, v, w, theta = r2j.solver_step_jax(
            u, v, w, theta, jkx, jky, jkz, k2, dealias, nq,
            nu_j, jnp.array(dt, dtype=jnp.float64))
        t += dt
        if step % sample_every == 0 or t >= t_end:
            snap = axisym_snapshot(np.asarray(u), np.asarray(v), np.asarray(w),
                                   grid, kx, ky, kz, N)
            snap.update({"t": t, "step": step, "dt": dt})
            series.append(snap)
            if verbose:
                print(f"  step {step:4d} t={t:.4f} |G|max={snap['Gamma_absmax']:.6e} "
                      f"D_axi={snap['axisym_defect']:.3e} "
                      f"src={snap['source_absmax']:.3e}")

    mp = gamma_max_principle_check([sg["Gamma_absmax"] for sg in series],
                                   [sg["axisym_defect"] for sg in series])
    mp_core = gamma_max_principle_check([sg["Gamma_absmax_core"] for sg in series],
                                        [sg["axisym_defect"] for sg in series])

    last = series[-1]
    r_c, z_c, src_table = localise_rz(
        swirl_source_eta(np.asarray(u), np.asarray(v), grid["r"],
                         grid["cos_phi"], grid["sin_phi"],
                         grid["axis_mask"], kz)[0], grid["r"], N)

    return {
        "exp_id": None,
        "N": N, "nu": nu, "t_end": t_end, "seed": seed,
        "ic": "axisym_swirl",
        "n_steps": series[-1]["step"],
        "n_samples": len(series),
        "t_final": t,
        "projection_residual_rel": proj_res_rel,
        "mp": mp,
        "mp_core": mp_core,
        "tail_fraction_final": series[-1]["tail_fraction"],
        "tail_fraction_max": float(max(sg["tail_fraction"] for sg in series)),
        "series": series,
        "source_table_r": r_c.tolist(),
        "source_table_z": z_c.tolist(),
        "source_table": src_table.tolist(),
        "source_rel_discrepancy_innermost_max": float(
            max(sg["source_rel_discrepancy_innermost"] for sg in series)),
        "source_r_loss_radius_max": float(
            max(sg["source_r_loss_radius"] for sg in series)),
        "axisym_defect_final": last["axisym_defect"],
        "verdict": mp["verdict"],
        "wall_time_s": time.perf_counter() - t_wall0,
    }


# =============================================================================
# Calibration drivers (run before any experiment is trusted)
# =============================================================================

def calibrate_burgers(N: int = 64, a: float = 1.0, nu: float = 1.0e-2,
                      circulation: float = 1.0, r_fit_max: float = 1.5) -> dict:
    """Diagnostics vs the exact Burgers vortex closed forms.

    Checked on r <= r_fit_max only: outside that radius the non-periodic
    Burgers field is being wrapped by the box and the comparison is against a
    field the grid does not actually represent.  Stated, not hidden.
    """
    (u, v, w), exact, grid = burgers_vortex_fields(N, a=a, nu=nu,
                                                   circulation=circulation)
    kx, ky, kz = spectral_wavenumbers_3d(N)
    r, cp, sp, am = grid["r"], grid["cos_phi"], grid["sin_phi"], grid["axis_mask"]
    sel = (r <= r_fit_max) & (~am)

    G = gamma_field(u, v, r, cp, sp)
    err_G = float(np.max(np.abs(G[sel] - exact["Gamma"][sel])) /
                  max(np.max(np.abs(exact["Gamma"][sel])), 1e-300))

    wx, wy, wz = vorticity_spectral(u, v, w, kx, ky, kz)
    wtr, w_th = omega_theta_over_r(wx, wy, r, cp, sp, am)
    scale_w = float(np.max(np.abs(wz[sel])))
    err_wth = float(np.max(np.abs(w_th[sel])) / max(scale_w, 1e-300))

    src, _ = swirl_source_eta(u, v, r, cp, sp, am, kz)
    eta = eta_field(u, v, r, cp, sp, am)
    scale_src = float(np.max(eta[sel] ** 2))
    err_src = float(np.max(np.abs(src[sel])) / max(scale_src, 1e-300))

    return {
        "N": N, "r_fit_max": r_fit_max,
        "gamma_rel_err": err_G,
        "omega_theta_rel_err": err_wth,
        "source_rel_err": err_src,
        "verdict": "PASS" if (err_G < 1e-12 and err_wth < 1e-8 and err_src < 1e-8) else "FAIL",
    }


def calibrate_planted_source(N: int = 64, c: float = 1.0, s: float = 1.0,
                             r_fit_max: float = 2.0) -> dict:
    """Source diagnostic vs a planted field with a KNOWN NON-ZERO source.

    Burgers gives source == 0, which cannot distinguish a correct source
    routine from one that returns zero.  This is the non-degenerate half of
    the calibration.
    """
    (u, v, w), exact, grid = planted_swirl_fields(N, c=c, s=s)
    _, _, kz = spectral_wavenumbers_3d(N)
    r, cp, sp, am = grid["r"], grid["cos_phi"], grid["sin_phi"], grid["axis_mask"]
    sel = (r <= r_fit_max) & (~am)

    src_eta, eta = swirl_source_eta(u, v, r, cp, sp, am, kz)
    src_naive = swirl_source_naive(u, v, r, cp, sp, am, kz)
    ref = exact["source"]
    scale = float(np.max(np.abs(ref[sel])))
    err_eta = float(np.max(np.abs(src_eta[sel] - ref[sel])) / scale)
    err_naive = float(np.max(np.abs(src_naive[sel] - ref[sel])) / scale)
    near = source_near_axis_report(src_eta, src_naive, r, am)
    return {
        "N": N, "r_fit_max": r_fit_max,
        "source_eta_rel_err": err_eta,
        "source_naive_rel_err": err_naive,
        "near_axis": near,
        "verdict": "PASS" if err_eta < 1e-10 else "FAIL",
    }


def calibrate_solver_beltrami(N: int = 32, nu: float = 1.0e-2,
                              t_end: float = 0.5, max_steps: int = 200) -> dict:
    """The solver's 'holds an exact solution' test, on a field that IS periodic.

    ABC flow is Beltrami, so u(t) = u(0) exp(-nu t) exactly and
    E(t)/E(0) = exp(-2 nu t).  Reported as the relative error in that ratio.
    The Burgers vortex cannot serve this purpose (not periodic).
    """
    os.environ.setdefault("JAX_ENABLE_X64", "1")
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    import route2_3D_jax as r2j  # noqa: E402

    u0, v0, w0 = abc_flow(N)
    solver = r2j.make_solver_jax(N=N, nu=nu, eps_param=0.1, dt_max=0.02,
                                 dtype=jnp.float64)
    solver = r2j.set_ic_jax(solver, {"u": u0, "v": v0, "w": w0})
    E0 = r2j.kinetic_energy_jax(solver["u"], solver["v"], solver["w"])
    out = r2j.run_jax(solver, t_end=t_end, max_steps=max_steps, verbose=False)
    E1 = out["E_final"]
    expected = float(np.exp(-2.0 * nu * out["t_final"]))
    got = E1 / E0
    rel = abs(got - expected) / expected
    # shape error: u(t) should be a pure rescaling of u(0)
    scale = np.exp(-nu * out["t_final"])
    shape_err = float(np.sqrt(np.mean(
        (np.asarray(out["u"]) - scale * np.asarray(solver["u"])) ** 2 +
        (np.asarray(out["v"]) - scale * np.asarray(solver["v"])) ** 2 +
        (np.asarray(out["w"]) - scale * np.asarray(solver["w"])) ** 2))
        / np.sqrt(np.mean(np.asarray(solver["u"]) ** 2 +
                          np.asarray(solver["v"]) ** 2 +
                          np.asarray(solver["w"]) ** 2)))
    return {
        "N": N, "nu": nu, "t_final": out["t_final"], "n_steps": out["n_steps"],
        "energy_ratio": got, "energy_ratio_expected": expected,
        "energy_rel_err": rel, "shape_rel_err": shape_err,
        "verdict": "PASS" if (rel < 1e-6 and shape_err < 1e-5) else "FAIL",
    }


# =============================================================================
# results.db
# =============================================================================

def _ensure_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_axisym_experiments (
            id                          INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id                      TEXT UNIQUE NOT NULL,
            claim_id                    TEXT NOT NULL,
            timestamp                   TEXT NOT NULL,
            ic                          TEXT,
            N                           INTEGER,
            nu                          REAL,
            seed                        INTEGER,
            n_steps                     INTEGER,
            n_samples                   INTEGER,
            t_final                     REAL,
            projection_residual_rel     REAL,
            Gamma0                      REAL,
            Gamma_absmax_over_run       REAL,
            Gamma_final                 REAL,
            mp_max_relative_rise        REAL,
            mp_tolerance_used           REAL,
            mp_violated                 INTEGER,
            axisym_defect_max           REAL,
            tail_fraction_max           REAL,
            source_absmax_final         REAL,
            source_rel_disc_innermost   REAL,
            source_r_loss_radius        REAL,
            verdict                     TEXT NOT NULL,
            wall_time_s                 REAL,
            key_metric                  TEXT NOT NULL,
            notes                       TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_axisym_timeseries (
            id                     INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id                 TEXT NOT NULL,
            claim_id               TEXT NOT NULL,
            step                   INTEGER,
            t                      REAL,
            Gamma_absmax           REAL,
            omega_theta_over_r_absmax REAL,
            source_absmax          REAL,
            source_naive_absmax    REAL,
            axisym_defect          REAL,
            M                      REAL,
            tail_fraction          REAL,
            r_at_Gamma_max         REAL
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS layer4_axisym_calibration (
            id           INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id       TEXT UNIQUE NOT NULL,
            claim_id     TEXT NOT NULL,
            timestamp    TEXT NOT NULL,
            kind         TEXT NOT NULL,
            N            INTEGER,
            verdict      TEXT NOT NULL,
            key_metric   TEXT NOT NULL,
            detail_json  TEXT
        )
    """)
    conn.commit()
    return conn


def log_axisym_result(db_path: str, exp_id: str, claim_id: str, res: dict,
                      notes: str = "") -> None:
    """Append one WP4a experiment (+ timeseries).  claim_id and key_metric are
    both NOT NULL by schema: this program does not record a result without
    them (16 legacy rows lack claim_id; no seventeenth)."""
    if not claim_id:
        raise ValueError("claim_id is mandatory — a result without one is not a result")
    conn = _ensure_db(db_path)
    mp = res["mp"]
    last = res["series"][-1]
    key_metric = (f"max_rel_rise(|Gamma|_inf)={mp['max_relative_rise']:.3e} "
                  f"tol={mp['tolerance_used']:.3e} "
                  f"D_axi_max={mp['axisym_defect_max']:.3e} "
                  f"tail_frac_max={res['tail_fraction_max']:.3e}")
    conn.execute("""
        INSERT OR REPLACE INTO layer4_axisym_experiments
        (exp_id, claim_id, timestamp, ic, N, nu, seed, n_steps, n_samples,
         t_final, projection_residual_rel, Gamma0, Gamma_absmax_over_run,
         Gamma_final, mp_max_relative_rise, mp_tolerance_used, mp_violated,
         axisym_defect_max, tail_fraction_max, source_absmax_final,
         source_rel_disc_innermost, source_r_loss_radius, verdict,
         wall_time_s, key_metric, notes)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, (exp_id, claim_id, datetime.now().isoformat(), res["ic"], res["N"],
          res["nu"], res["seed"], res["n_steps"], res["n_samples"],
          res["t_final"], res["projection_residual_rel"], mp["Gamma0"],
          mp["Gamma_absmax_over_run"], mp["Gamma_final"],
          mp["max_relative_rise"], mp["tolerance_used"], int(mp["violated"]),
          mp["axisym_defect_max"], res["tail_fraction_max"], last["source_absmax"],
          res["source_rel_discrepancy_innermost_max"],
          res["source_r_loss_radius_max"], res["verdict"],
          res["wall_time_s"], key_metric, notes))
    conn.execute("DELETE FROM layer4_axisym_timeseries WHERE exp_id = ?", (exp_id,))
    for sg in res["series"]:
        conn.execute("""
            INSERT INTO layer4_axisym_timeseries
            (exp_id, claim_id, step, t, Gamma_absmax, omega_theta_over_r_absmax,
             source_absmax, source_naive_absmax, axisym_defect, M,
             tail_fraction, r_at_Gamma_max)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
        """, (exp_id, claim_id, sg["step"], sg["t"], sg["Gamma_absmax"],
              sg["omega_theta_over_r_absmax"], sg["source_absmax"],
              sg["source_naive_absmax"], sg["axisym_defect"], sg["M"],
              sg["tail_fraction"], sg["r_at_Gamma_max"]))
    conn.commit()
    conn.close()


def log_calibration(db_path: str, exp_id: str, claim_id: str, kind: str,
                    res: dict, key_metric: str) -> None:
    if not claim_id:
        raise ValueError("claim_id is mandatory")
    import json
    conn = _ensure_db(db_path)
    safe = {k: v for k, v in res.items() if not isinstance(v, np.ndarray)}
    conn.execute("""
        INSERT OR REPLACE INTO layer4_axisym_calibration
        (exp_id, claim_id, timestamp, kind, N, verdict, key_metric, detail_json)
        VALUES (?,?,?,?,?,?,?,?)
    """, (exp_id, claim_id, datetime.now().isoformat(), kind, res.get("N"),
          res["verdict"], key_metric, json.dumps(safe, default=float)))
    conn.commit()
    conn.close()


__all__ = [
    "spectral_wavenumbers_3d", "cylindrical_grid", "to_cylindrical",
    "vorticity_spectral", "ddz_spectral", "gamma_field", "eta_field",
    "omega_theta_over_r", "swirl_source_eta", "swirl_source_naive",
    "source_near_axis_report", "localise_rz", "axisymmetry_defect",
    "burgers_vortex_fields", "planted_swirl_fields", "abc_flow",
    "axisym_snapshot", "gamma_max_principle_check", "axisym_swirl_ic",
    "polar_resample", "spectral_tail_fraction",
    "run_axisym_experiment", "calibrate_burgers", "calibrate_planted_source",
    "calibrate_solver_beltrami", "log_axisym_result", "log_calibration",
    "make_axisym_figure",
    "AXIS_EPS_CELLS", "MIN_FIT_POINTS", "MIN_CORR_POINTS", "MP_TOL_REL",
]


# =============================================================================
# Figures
# =============================================================================

def make_axisym_figure(results_by_N: dict, out_path: str) -> str:
    """Gamma maximum-principle convergence + the (r,z) source localisation.

    `results_by_N` maps N -> the dict returned by run_axisym_experiment.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:      # a missing plotting library must never lose a run
        return ""

    Ns = sorted(results_by_N)
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.3))

    for N in Ns:
        r = results_by_N[N]
        t = [s["t"] for s in r["series"]]
        g = np.array([s["Gamma_absmax"] for s in r["series"]])
        ax[0].plot(t, g / g[0], marker=".", label=f"N={N}")
    ax[0].axhline(1.0, color="k", lw=0.8, ls="--")
    ax[0].set_xlabel("t")
    ax[0].set_ylabel(r"$\max|\Gamma|(t)\,/\,\max|\Gamma|(0)$")
    ax[0].set_title("(MP): must be non-increasing")
    ax[0].legend(fontsize=8)

    rise = [results_by_N[N]["mp"]["max_relative_rise"] for N in Ns]
    tail = [results_by_N[N]["tail_fraction_max"] for N in Ns]
    dax = [results_by_N[N]["mp"]["axisym_defect_max"] for N in Ns]
    ax[1].loglog(Ns, rise, "o-", label="max relative rise of $|\\Gamma|_\\infty$")
    ax[1].loglog(Ns, tail, "s-", label="spectral tail fraction")
    ax[1].loglog(Ns, dax, "^-", label="axisymmetry defect $D$")
    ax[1].set_xlabel("N")
    ax[1].set_title("the violation converges away")
    ax[1].legend(fontsize=8)
    ax[1].grid(True, which="both", alpha=0.3)

    N = Ns[-1]
    tab = np.array(results_by_N[N]["source_table"])
    rc = results_by_N[N]["source_table_r"]
    zc = results_by_N[N]["source_table_z"]
    im = ax[2].pcolormesh(zc, rc, np.log10(tab + 1e-30), shading="auto")
    ax[2].set_xlabel("z")
    ax[2].set_ylabel("r")
    ax[2].set_title(r"$\log_{10}\max_\phi|\partial_z(\Gamma^2)/r^4|$"
                    f"  (N={N}, final t)")
    fig.colorbar(im, ax=ax[2])

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path
