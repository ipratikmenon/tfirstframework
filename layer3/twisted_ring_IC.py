"""
layer3/twisted_ring_IC.py — Twisted vortex ring, parameterised by rotational
transform iota (WP4b, HANDOVER-S46-OPUS §WP4b)
=============================================================================
The stellarator geometry of S46, read as a vortex tube: a torus whose field
lines wind poloidally as they go round toroidally.  `iota` is the rotational
transform -- poloidal turns per toroidal circuit.  This is the geometric
embodiment of helicity (linking = twist + writhe, Calugareanu-White; the
rotational transform is the twist term).

CONSTRUCTION.  Local toroidal frame about the core circle of major radius R0
lying in the plane z = z0, centred on the box axis:

    s   = cylindrical radius from the axis,    phi = azimuth
    rho = sqrt((s - R0)^2 + (z - z0)^2)        (poloidal distance from core)
    chi = atan2(z - z0, s - R0)                (poloidal angle)

    e_phi = (-sin phi,  cos phi, 0)                        (toroidal)
    e_chi = (-sin chi cos phi, -sin chi sin phi, cos chi)  (poloidal)

A field line at poloidal radius rho with dchi/dphi = iota has arc-length ratio
    (rho dchi) / (R0 dphi) = iota rho / R0,
so the vorticity of a tube with rotational transform iota is

    omega = f(rho) [ e_phi + (iota rho / R0) e_chi ],
    f(rho) = (Circ / pi a^2) exp(-rho^2/a^2) * cutoff(rho)

Velocity is recovered by Biot-Savart in Fourier space, u_hat = i k x omega_hat
/ |k|^2, which returns a divergence-free u whose curl is the divergence-free
projection P(omega) of the constructed field.  The toroidal metric factor
h_phi = R0 + rho cos chi makes the constructed omega divergence-free only to
O(a/R0); `omega_div_defect` measures ||omega - curl u|| / ||omega|| so that
this known O(a/R0) correction is reported rather than assumed small.

NAMED FAILURE MODES GUARDED HERE (HANDOVER §WP4b):
  (1) "A twisted ring with zero net helicity by symmetry."  H != 0 for
      iota != 0 is asserted in layer3/test_twisted_ring_IC.py before any run.
  (3) "Ring core under 4 cells (S40)."  `check_resolution` raises unless
      a/dx >= MIN_CORE_CELLS = 4.

FLOAT64.  `layer3.blowup_search_3D.biot_savart_spectral` exists and is the
natural thing to reuse, but it returns float32 (it was written for the Metal
backend).  WP4a's Gamma maximum-principle test and WP4b's correlations are
cancellation-sensitive, so this module keeps a float64 Biot-Savart.  It is
NOT an independent reimplementation: test_twisted_ring_IC.py asserts the two
agree to float32 tolerance on the same input.

PROPERTY ENGINE.  Route 2 is the exact Prize system (nu constant, a numerical
parameter).  There is no thermophysical property here to route through
layer1/tfirst_props.py.
"""

from __future__ import annotations

import numpy as np

MIN_CORE_CELLS = 4.0      # S40 rule: a vortex core thinner than this is noise
DEFAULT_R0 = 1.35
DEFAULT_A = 0.55
DEFAULT_RHO_CUT_A = 2.0   # Gaussian truncated at 2a (exp(-4) ~ 1.8e-2), with
                          # a cosine taper over the outer 20% so the truncation
                          # does not ring in Fourier space.
# Why these values and not a wider cut: the constraints
#     a/dx >= 4                (S40, HANDOVER failure mode 3)
#     rho_cut_a * a  <  R0     (tube must not pass through the axis)
#     R0 + rho_cut_a * a < pi  (ring inside the box, clear of its images)
# are jointly tight.  With rho_cut_a = 3 they admit no (R0, a) that is also
# resolved at N = 48-64, which is how the earlier R0=1.0, a=0.5, cut=3a choice
# came to be self-intersecting; the third guard below now rejects it.


def check_resolution(N: int, a: float = DEFAULT_A, R0: float = DEFAULT_R0,
                     rho_cut_a: float = DEFAULT_RHO_CUT_A) -> dict:
    """Grid-adequacy guard.  Raises on failure mode (3): core under 4 cells.

    Also checks that the ring plus its cut-off core fits inside the box with
    margin, so the periodic images are not part of the physics.
    """
    dx = 2.0 * np.pi / N
    cells = a / dx
    outer = R0 + rho_cut_a * a
    margin = np.pi - outer
    info = {"N": N, "dx": dx, "core_cells": float(cells), "outer_radius": float(outer),
            "box_margin": float(margin), "a": a, "R0": R0}
    if cells < MIN_CORE_CELLS:
        raise ValueError(
            f"ring core under-resolved: a/dx = {cells:.2f} < {MIN_CORE_CELLS} "
            f"(N={N}, a={a}).  HANDOVER §WP4b named failure mode 3 (S40).")
    if margin <= 0.0:
        raise ValueError(f"ring does not fit in the box: outer radius {outer:.3f} "
                         f"exceeds half-width pi")
    if rho_cut_a * a >= R0:
        raise ValueError(
            f"ring tube self-intersects through the axis: rho_cut = "
            f"{rho_cut_a * a:.3f} >= R0 = {R0:.3f}.  The toroidal frame is "
            f"singular there and the constructed omega acquires a spurious "
            f"divergence; pick a smaller core or a larger major radius.")
    info["rho_cut"] = float(rho_cut_a * a)
    return info


def toroidal_frame(N: int, R0: float = DEFAULT_R0,
                   x0: float | None = None, y0: float | None = None,
                   z0: float | None = None):
    """Local toroidal coordinates (rho, chi) and unit vectors e_phi, e_chi.

    Returns dict with rho, and the Cartesian components of e_phi and e_chi.
    """
    if x0 is None:
        x0 = np.pi
    if y0 is None:
        y0 = np.pi
    if z0 is None:
        z0 = np.pi
    g = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    dxv = (g - x0 + np.pi) % (2.0 * np.pi) - np.pi
    dyv = (g - y0 + np.pi) % (2.0 * np.pi) - np.pi
    dzv = (g - z0 + np.pi) % (2.0 * np.pi) - np.pi
    DX = np.broadcast_to(dxv[:, None, None], (N, N, N))
    DY = np.broadcast_to(dyv[None, :, None], (N, N, N))
    DZ = np.broadcast_to(dzv[None, None, :], (N, N, N))

    s = np.sqrt(DX ** 2 + DY ** 2)
    s_safe = np.where(s > 1e-12, s, 1.0)
    cphi = np.where(s > 1e-12, DX / s_safe, 1.0)
    sphi = np.where(s > 1e-12, DY / s_safe, 0.0)

    ds = s - R0
    rho = np.sqrt(ds ** 2 + DZ ** 2)
    rho_safe = np.where(rho > 1e-12, rho, 1.0)
    cchi = np.where(rho > 1e-12, ds / rho_safe, 1.0)
    schi = np.where(rho > 1e-12, DZ / rho_safe, 0.0)

    return {
        "rho": np.ascontiguousarray(rho),
        "s": np.ascontiguousarray(s),
        "e_phi": (np.ascontiguousarray(-sphi), np.ascontiguousarray(cphi),
                  np.zeros((N, N, N))),
        "e_chi": (np.ascontiguousarray(-schi * cphi),
                  np.ascontiguousarray(-schi * sphi),
                  np.ascontiguousarray(cchi)),
        "cos_phi": np.ascontiguousarray(cphi),
        "sin_phi": np.ascontiguousarray(sphi),
    }


def twisted_ring_omega(N: int, iota: float = 0.0, R0: float = DEFAULT_R0,
                       a: float = DEFAULT_A, circulation: float = 1.0,
                       rho_cut_a: float = DEFAULT_RHO_CUT_A,
                       x0=None, y0=None, z0=None):
    """Vorticity of a twisted vortex ring with rotational transform iota."""
    fr = toroidal_frame(N, R0=R0, x0=x0, y0=y0, z0=z0)
    rho = fr["rho"]
    rho_cut = rho_cut_a * a
    # C^1 cosine taper over the outer 20% of the cut radius, so the truncated
    # Gaussian does not ring in Fourier space.
    taper_lo = 0.8 * rho_cut
    tap = np.ones_like(rho)
    band = (rho > taper_lo) & (rho <= rho_cut)
    tap = np.where(band, 0.5 * (1.0 + np.cos(np.pi * (rho - taper_lo) /
                                             (rho_cut - taper_lo))), tap)
    tap = np.where(rho > rho_cut, 0.0, tap)

    f = circulation / (np.pi * a ** 2) * np.exp(-rho ** 2 / a ** 2) * tap
    tw = iota * rho / R0
    ex = f * (fr["e_phi"][0] + tw * fr["e_chi"][0])
    ey = f * (fr["e_phi"][1] + tw * fr["e_chi"][1])
    ez = f * (fr["e_phi"][2] + tw * fr["e_chi"][2])
    return (np.ascontiguousarray(ex), np.ascontiguousarray(ey),
            np.ascontiguousarray(ez)), fr


def biot_savart_f64(ox, oy, oz, kx, ky, kz):
    """u_hat = i k x omega_hat / |k|^2, in float64.

    Divergence-free by construction; curl u = P(omega).  float64 mirror of
    layer3.blowup_search_3D.biot_savart_spectral (which returns float32);
    agreement between the two is asserted in the test suite.
    """
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    k2s = np.where(k2 > 0.0, k2, 1.0)
    valid = (k2 > 0.0)
    oxh, oyh, ozh = np.fft.fftn(ox), np.fft.fftn(oy), np.fft.fftn(oz)
    cx = ky * ozh - kz * oyh
    cy = kz * oxh - kx * ozh
    cz = kx * oyh - ky * oxh
    u = np.real(np.fft.ifftn(np.where(valid, 1j * cx / k2s, 0.0)))
    v = np.real(np.fft.ifftn(np.where(valid, 1j * cy / k2s, 0.0)))
    w = np.real(np.fft.ifftn(np.where(valid, 1j * cz / k2s, 0.0)))
    return (np.ascontiguousarray(u, dtype=np.float64),
            np.ascontiguousarray(v, dtype=np.float64),
            np.ascontiguousarray(w, dtype=np.float64))


def twisted_ring_ic(N: int, iota: float = 0.0, R0: float = DEFAULT_R0,
                    a: float = DEFAULT_A, circulation: float = 1.0,
                    rho_cut_a: float = DEFAULT_RHO_CUT_A,
                    check_res: bool = True, seed: int = 42) -> dict:
    """Twisted vortex ring IC in the {'u','v','w','name'} dict shape used by
    route2_3D / route2_3D_jax set_ic.

    Reported alongside: `omega_div_defect` = ||omega - curl u||/||omega||, the
    size of the O(a/R0) toroidal divergence of the constructed field that the
    Biot-Savart projection removed.
    """
    res = check_resolution(N, a=a, R0=R0, rho_cut_a=rho_cut_a) if check_res else {}
    k1 = np.fft.fftfreq(N, d=1.0 / N).astype(np.float64)
    kx, ky, kz = k1[:, None, None], k1[None, :, None], k1[None, None, :]

    (ox, oy, oz), fr = twisted_ring_omega(N, iota=iota, R0=R0, a=a,
                                          circulation=circulation,
                                          rho_cut_a=rho_cut_a)
    u, v, w = biot_savart_f64(ox, oy, oz, kx, ky, kz)

    # curl of the recovered velocity == P(omega); the residual is the
    # divergence part the construction carried.
    d = lambda fh, k: np.real(np.fft.ifftn(1j * k * fh))
    uh, vh, wh = np.fft.fftn(u), np.fft.fftn(v), np.fft.fftn(w)
    cx = d(wh, ky) - d(vh, kz)
    cy = d(uh, kz) - d(wh, kx)
    cz = d(vh, kx) - d(uh, ky)
    num = float(np.sqrt(np.mean((cx - ox) ** 2 + (cy - oy) ** 2 + (cz - oz) ** 2)))
    den = float(np.sqrt(np.mean(ox ** 2 + oy ** 2 + oz ** 2)))

    H = float(np.mean(u * cx + v * cy + w * cz) * (2.0 * np.pi) ** 3)
    return {
        "u": u, "v": v, "w": w,
        "name": f"twisted_ring_iota{iota:+.3f}",
        "iota": float(iota), "R0": float(R0), "a": float(a),
        "circulation": float(circulation), "seed": int(seed),
        "omega_div_defect": num / den if den > 0 else 0.0,
        "helicity_ic": H,
        "resolution": res,
        "rho": fr["rho"],
    }


__all__ = ["MIN_CORE_CELLS", "DEFAULT_R0", "DEFAULT_A", "DEFAULT_RHO_CUT_A",
           "check_resolution", "toroidal_frame", "twisted_ring_omega",
           "biot_savart_f64", "twisted_ring_ic"]
