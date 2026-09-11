"""Tests for layer3/twisted_ring_IC.py (WP4b).

The named failure modes of HANDOVER-S46-OPUS §WP4b are tests here, not hopes:
  (1) zero net helicity by symmetry  -> test_helicity_nonzero_for_nonzero_iota
  (3) ring core under 4 cells (S40)  -> test_core_resolution_guard_raises
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from twisted_ring_IC import (  # noqa: E402
    MIN_CORE_CELLS, DEFAULT_R0, DEFAULT_A, check_resolution, toroidal_frame,
    twisted_ring_omega, biot_savart_f64, twisted_ring_ic,
)

# N=32 cannot resolve a ring at all under the S40 rule (a/dx >= 4 forces
# a >= 0.79, which closes the ring hole for any R0 that fits the box), so 48
# is the smallest grid used here.  That is the guard working, not a test
# inconvenience.
N_SMALL = 48
N_MED = 48
A_OK = 0.55
R0_OK = 1.35


def _k(N):
    k1 = np.fft.fftfreq(N, d=1.0 / N).astype(np.float64)
    return k1[:, None, None], k1[None, :, None], k1[None, None, :]


def _curl(u, v, w, kx, ky, kz):
    d = lambda f, k: np.real(np.fft.ifftn(1j * k * np.fft.fftn(f)))
    return (d(w, ky) - d(v, kz), d(u, kz) - d(w, kx), d(v, kx) - d(u, ky))


# ── resolution guard (failure mode 3, S40) ───────────────────────────────────

def test_core_resolution_guard_raises():
    """a/dx < 4 must raise, not warn."""
    with pytest.raises(ValueError, match="under-resolved"):
        check_resolution(N=16, a=0.1)


def test_core_resolution_guard_passes_at_default():
    info = check_resolution(N=64, a=DEFAULT_A, R0=DEFAULT_R0)
    assert info["core_cells"] >= MIN_CORE_CELLS
    assert info["box_margin"] > 0.0


def test_self_intersecting_tube_is_rejected():
    """The geometry bug found at S50: rho_cut >= R0 puts tube material on the
    axis, where the toroidal frame is singular, and the constructed vorticity
    picks up a spurious divergence."""
    with pytest.raises(ValueError, match="self-intersects"):
        check_resolution(N=64, a=0.50, R0=1.0, rho_cut_a=3.0)


def test_default_geometry_does_not_self_intersect():
    info = check_resolution(N=64, a=DEFAULT_A, R0=DEFAULT_R0)
    assert info["rho_cut"] < info["R0"]


def test_ring_must_fit_in_box():
    with pytest.raises(ValueError, match="does not fit"):
        check_resolution(N=256, a=0.5, R0=3.0)


def test_ic_respects_guard_by_default():
    with pytest.raises(ValueError):
        twisted_ring_ic(16, iota=1.0, a=0.1)


# ── toroidal frame ───────────────────────────────────────────────────────────

def test_frame_vectors_orthonormal():
    fr = toroidal_frame(N_SMALL)
    ep, ec = fr["e_phi"], fr["e_chi"]
    dot = ep[0] * ec[0] + ep[1] * ec[1] + ep[2] * ec[2]
    assert np.max(np.abs(dot)) < 1e-12
    for e in (ep, ec):
        nrm = np.sqrt(e[0] ** 2 + e[1] ** 2 + e[2] ** 2)
        assert np.max(np.abs(nrm - 1.0)) < 1e-12


def test_rho_vanishes_on_core_circle():
    fr = toroidal_frame(N_MED, R0=DEFAULT_R0)
    assert float(np.min(fr["rho"])) < 2.0 * np.pi / N_MED


# ── Biot-Savart ──────────────────────────────────────────────────────────────

def test_biot_savart_gives_divergence_free_velocity():
    kx, ky, kz = _k(N_SMALL)
    (ox, oy, oz), _ = twisted_ring_omega(N_SMALL, iota=1.0, a=A_OK, R0=R0_OK)
    u, v, w = biot_savart_f64(ox, oy, oz, kx, ky, kz)
    d = lambda f, k: np.real(np.fft.ifftn(1j * k * np.fft.fftn(f)))
    div = d(u, kx) + d(v, ky) + d(w, kz)
    scale = float(np.sqrt(np.mean(u ** 2 + v ** 2 + w ** 2)))
    assert float(np.sqrt(np.mean(div ** 2))) < 1e-13 * max(scale, 1e-30)


def test_biot_savart_f64_matches_existing_float32_routine():
    """Not an independent reimplementation: must agree with
    layer3.blowup_search_3D.biot_savart_spectral to float32 tolerance."""
    from blowup_search_3D import biot_savart_spectral
    kx, ky, kz = _k(N_SMALL)
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    (ox, oy, oz), _ = twisted_ring_omega(N_SMALL, iota=1.0, a=A_OK, R0=R0_OK)
    a = biot_savart_f64(ox, oy, oz, kx, ky, kz)
    b = biot_savart_spectral(ox, oy, oz, kx, ky, kz, k2)
    scale = float(np.max(np.abs(a[0])))
    for p, q in zip(a, b):
        assert np.max(np.abs(p - np.asarray(q, float))) < 1e-5 * max(scale, 1e-30)


def test_biot_savart_output_is_float64():
    kx, ky, kz = _k(N_SMALL)
    (ox, oy, oz), _ = twisted_ring_omega(N_SMALL, iota=0.5, a=A_OK, R0=R0_OK)
    for c in biot_savart_f64(ox, oy, oz, kx, ky, kz):
        assert c.dtype == np.float64


# ── HANDOVER failure mode (1): H != 0 for iota != 0 ──────────────────────────

def test_helicity_nonzero_for_nonzero_iota():
    """The named failure mode: a twist that is helicity-neutral by symmetry.

    Checked BEFORE any run, as the handover requires.
    """
    h0 = twisted_ring_ic(N_MED, iota=0.0, a=A_OK, R0=R0_OK)["helicity_ic"]
    for iota in (0.5, 1.0, 2.0):
        h = twisted_ring_ic(N_MED, iota=iota, a=A_OK, R0=R0_OK)["helicity_ic"]
        assert abs(h) > 1e-3, f"iota={iota} gave helicity-neutral ring H={h}"
        assert abs(h) > 1e6 * abs(h0)


def test_untwisted_ring_has_zero_helicity():
    h = twisted_ring_ic(N_MED, iota=0.0, a=A_OK, R0=R0_OK)["helicity_ic"]
    assert abs(h) < 1e-12


def test_helicity_is_odd_and_linear_in_iota():
    hs = {i: twisted_ring_ic(N_MED, iota=i, a=A_OK, R0=R0_OK)["helicity_ic"]
          for i in (-2.0, -1.0, 1.0, 2.0)}
    assert hs[1.0] * hs[-1.0] < 0
    assert abs(hs[1.0] + hs[-1.0]) < 1e-9 * abs(hs[1.0])
    assert abs(hs[2.0] / hs[1.0] - 2.0) < 1e-6


# ── constructed field quality ────────────────────────────────────────────────

def test_omega_div_defect_is_small_and_grows_with_iota():
    """The O(a/R0) toroidal divergence is reported, not assumed away."""
    d0 = twisted_ring_ic(N_MED, iota=0.0, a=A_OK, R0=R0_OK)["omega_div_defect"]
    d1 = twisted_ring_ic(N_MED, iota=1.0, a=A_OK, R0=R0_OK)["omega_div_defect"]
    d2 = twisted_ring_ic(N_MED, iota=2.0, a=A_OK, R0=R0_OK)["omega_div_defect"]
    # iota = 0: analytically divergence-free (f(rho) e_phi has no divergence),
    # so what is left is pure grid truncation and it CONVERGES.  iota != 0
    # carries the genuine O(a/R0) toroidal divergence, which does not.
    assert d0 < 2e-3
    # the toroidal correction is O(a/R0) per unit iota; a/R0 = 0.41 here
    assert 0.0 < d1 < 0.20
    assert d2 > d1


def test_untwisted_divergence_defect_is_discretisation_and_converges():
    """iota=0 is analytically divergence-free, so its defect must fall with N;
    the iota != 0 defect is the physical O(a/R0) toroidal term and must not."""
    d_untwisted = [twisted_ring_ic(N, iota=0.0)["omega_div_defect"]
                   for N in (48, 96)]
    assert d_untwisted[1] < 0.25 * d_untwisted[0]
    d_twisted = [twisted_ring_ic(N, iota=1.0)["omega_div_defect"]
                 for N in (48, 96)]
    assert abs(d_twisted[1] - d_twisted[0]) < 0.02 * d_twisted[0]


def test_curl_of_velocity_recovers_projected_vorticity():
    kx, ky, kz = _k(N_MED)
    ic = twisted_ring_ic(N_MED, iota=1.0, a=A_OK, R0=R0_OK)
    (ox, oy, oz), _ = twisted_ring_omega(N_MED, iota=1.0, a=A_OK, R0=R0_OK)
    cx, cy, cz = _curl(ic["u"], ic["v"], ic["w"], kx, ky, kz)
    num = float(np.sqrt(np.mean((cx - ox) ** 2 + (cy - oy) ** 2 + (cz - oz) ** 2)))
    den = float(np.sqrt(np.mean(ox ** 2 + oy ** 2 + oz ** 2)))
    assert num / den == pytest.approx(ic["omega_div_defect"], rel=1e-9)


def test_vorticity_is_compactly_supported_away_from_box_edge():
    (ox, oy, oz), fr = twisted_ring_omega(N_MED, iota=1.0, a=A_OK, R0=R0_OK)
    mag = np.sqrt(ox ** 2 + oy ** 2 + oz ** 2)
    assert np.max(mag[0, :, :]) == 0.0
    assert np.max(mag[:, :, 0]) == 0.0


def test_taper_makes_omega_exactly_zero_outside_cut():
    (ox, oy, oz), fr = twisted_ring_omega(N_MED, iota=1.0, a=A_OK, R0=R0_OK, rho_cut_a=3.0)
    mag = np.sqrt(ox ** 2 + oy ** 2 + oz ** 2)
    assert np.max(mag[fr["rho"] > 3.0 * A_OK + 1e-9]) == 0.0


def test_ic_dict_shape_matches_solver_contract():
    ic = twisted_ring_ic(N_SMALL, iota=1.0, a=A_OK, R0=R0_OK)
    for key in ("u", "v", "w", "name"):
        assert key in ic
    for key in ("u", "v", "w"):
        assert ic[key].shape == (N_SMALL,) * 3
        assert np.all(np.isfinite(ic[key]))


def test_ic_is_deterministic():
    a = twisted_ring_ic(N_SMALL, iota=1.0, a=A_OK, R0=R0_OK)
    b = twisted_ring_ic(N_SMALL, iota=1.0, a=A_OK, R0=R0_OK)
    assert np.array_equal(a["u"], b["u"])


def test_mirror_pair_ic_is_related_by_parity():
    """IC(-iota) is the mirror image of IC(+iota): the premise of the parity
    argument in layer4/helicity_diagnostics.py.  Checked via helicity, the
    pseudoscalar invariant."""
    hp = twisted_ring_ic(N_MED, iota=1.5, a=A_OK, R0=R0_OK)["helicity_ic"]
    hm = twisted_ring_ic(N_MED, iota=-1.5, a=A_OK, R0=R0_OK)["helicity_ic"]
    assert hp == pytest.approx(-hm, rel=1e-10)


def test_circulation_scales_velocity_linearly():
    a = twisted_ring_ic(N_SMALL, iota=1.0, a=A_OK, R0=R0_OK, circulation=1.0)
    b = twisted_ring_ic(N_SMALL, iota=1.0, a=A_OK, R0=R0_OK, circulation=2.5)
    scale = float(np.max(np.abs(a["u"])))
    assert np.max(np.abs(b["u"] - 2.5 * a["u"])) < 1e-12 * scale


def test_helicity_scales_as_circulation_squared():
    a = twisted_ring_ic(N_MED, iota=1.0, a=A_OK, R0=R0_OK, circulation=1.0)
    b = twisted_ring_ic(N_MED, iota=1.0, a=A_OK, R0=R0_OK, circulation=3.0)
    assert b["helicity_ic"] == pytest.approx(9.0 * a["helicity_ic"], rel=1e-10)
