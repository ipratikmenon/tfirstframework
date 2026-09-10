"""
layer4/test_scale_diagnostics.py
Tests for the scale-r* diagnostics (layer4/scale_diagnostics.py) and the
bipolar jet IC / strain sweep (layer3/quasar_jet_IC.py).

Run: python -m pytest layer4/test_scale_diagnostics.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from layer4.scale_diagnostics import (
    sobolev_H3_norm,
    L_factor,
    _corr_ratio,
    scale_diagnostics,
    run_scale_experiment,
    l_dependence_summary,
    THETA_MINUS,
    THETA_PLUS,
)
from layer4.geometric_disorder import spectral_wavenumbers
from layer3.quasar_jet_IC import quasar_jet_ic, JET_SWEEP_RATIOS

N = 32
DX = 2.0 * np.pi / N
KX, KY, KZ = spectral_wavenumbers(N)
_X1D = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
_XX, _YY, _ZZ = np.meshgrid(_X1D, _X1D, _X1D, indexing="ij")


def _tg(N_=N):
    x = np.linspace(0.0, 2.0 * np.pi, N_, endpoint=False)
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    return (np.sin(xx) * np.cos(yy) * np.cos(zz),
            -np.cos(xx) * np.sin(yy) * np.cos(zz),
            np.zeros((N_, N_, N_)))


# ── sobolev_H3_norm / L_factor ───────────────────────────────────────────────

def test_H3_norm_zero_for_zero_field():
    z = np.zeros((N, N, N))
    assert sobolev_H3_norm(z, z, z, KX, KY, KZ) == pytest.approx(0.0, abs=1e-14)


def test_H3_norm_constant_field_matches_L2():
    """For a constant field only k=0 contributes, weight (1+0)^3 = 1, so the
    H^3 norm equals the L^2 norm = |c| (unit-measure normalisation)."""
    c = np.full((N, N, N), 2.0)
    z = np.zeros((N, N, N))
    assert sobolev_H3_norm(c, z, z, KX, KY, KZ) == pytest.approx(2.0, rel=1e-10)


def test_H3_norm_single_mode_exact():
    """u = sin(x): one conjugate pair at |k|=1, weight (1+1)^3 = 8.
    ||u||_{H^3} = sqrt(8 * ||sin||_{L^2}^2) = sqrt(8 * 1/2) = 2."""
    u = np.sin(_XX)
    z = np.zeros((N, N, N))
    assert sobolev_H3_norm(u, z, z, KX, KY, KZ) == pytest.approx(2.0, rel=1e-10)


def test_H3_norm_grows_with_wavenumber():
    z = np.zeros((N, N, N))
    lo = sobolev_H3_norm(np.sin(_XX), z, z, KX, KY, KZ)
    hi = sobolev_H3_norm(np.sin(4 * _XX), z, z, KX, KY, KZ)
    assert hi > lo


def test_L_factor_at_least_one_and_finite():
    u, v, w = _tg()
    L = L_factor(u, v, w, 1.0, KX, KY, KZ)
    assert np.isfinite(L) and L >= 1.0


def test_L_factor_degenerate_M_returns_one():
    u, v, w = _tg()
    assert L_factor(u, v, w, 0.0, KX, KY, KZ) == pytest.approx(1.0)


def test_L_factor_decreases_as_M_grows():
    u, v, w = _tg()
    assert L_factor(u, v, w, 100.0, KX, KY, KZ) < L_factor(u, v, w, 0.1, KX, KY, KZ)


# ── _corr_ratio conventions ──────────────────────────────────────────────────

def test_corr_ratio_uncorrelated_is_one():
    g = np.full((N, N, N), 3.0)
    w = np.full((N, N, N), 5.0)
    mask = np.zeros((N, N, N), bool); mask[:8] = True
    assert _corr_ratio(g, w, mask) == pytest.approx(1.0, rel=1e-12)


def test_corr_ratio_empty_mask_is_zero():
    g = np.ones((N, N, N)); w = np.ones((N, N, N))
    assert _corr_ratio(g, w, np.zeros((N, N, N), bool)) == 0.0


def test_corr_ratio_degenerate_factor_is_zero():
    g = np.ones((N, N, N)); w = np.zeros((N, N, N))
    mask = np.ones((N, N, N), bool)
    assert _corr_ratio(g, w, mask) == 0.0


def test_corr_ratio_positive_correlation_exceeds_one():
    """g and w both large on the same half of the mask -> c_K > 1."""
    g = np.zeros((N, N, N)); w = np.zeros((N, N, N))
    mask = np.zeros((N, N, N), bool); mask[:8] = True
    g[:4] = 10.0; g[4:8] = 1.0
    w[:4] = 10.0; w[4:8] = 1.0
    assert _corr_ratio(g, w, mask) > 1.0


# ── scale_diagnostics ────────────────────────────────────────────────────────

def test_scale_diagnostics_keys_and_finiteness():
    u, v, w = _tg()
    d = scale_diagnostics(u, v, w, DX)
    for k in ("M", "r_star", "L", "ball_cells", "frac_E_rstar",
              "frac_band_rstar", "c_K_rstar", "c_K_band",
              "cK_times_L", "cKband_times_L"):
        assert k in d and np.isfinite(d[k])


def test_scale_diagnostics_fractions_in_unit_interval():
    u, v, w = _tg()
    d = scale_diagnostics(u, v, w, DX)
    assert 0.0 <= d["frac_E_rstar"] <= 1.0
    assert 0.0 <= d["frac_band_rstar"] <= 1.0


def test_scale_diagnostics_ball_is_nonempty():
    u, v, w = _tg()
    assert scale_diagnostics(u, v, w, DX)["ball_cells"] >= 1


def test_scale_diagnostics_r_star_matches_M():
    u, v, w = _tg()
    d = scale_diagnostics(u, v, w, DX)
    assert d["r_star"] == pytest.approx(d["M"] ** -0.5, rel=1e-12)


def test_scale_diagnostics_products_consistent():
    u, v, w = _tg()
    d = scale_diagnostics(u, v, w, DX)
    assert d["cK_times_L"] == pytest.approx(d["c_K_rstar"] * d["L"], rel=1e-12)
    assert d["cKband_times_L"] == pytest.approx(d["c_K_band"] * d["L"], rel=1e-12)


def test_scale_diagnostics_zero_vorticity_raises():
    z = np.zeros((N, N, N))
    with pytest.raises(ValueError, match="identically zero"):
        scale_diagnostics(z, z, z, DX)


def test_scale_diagnostics_immutability():
    u, v, w = _tg()
    u0, v0, w0 = u.copy(), v.copy(), w.copy()
    scale_diagnostics(u, v, w, DX)
    assert np.array_equal(u, u0) and np.array_equal(v, v0) and np.array_equal(w, w0)


def test_band_and_E_levels_are_the_documented_ones():
    """Definition def:annulus-residual: Ann = {M/8 <= |omega| <= M/4}."""
    assert THETA_MINUS == pytest.approx(0.125)
    assert THETA_PLUS == pytest.approx(0.25)


# ── quasar_jet_ic ────────────────────────────────────────────────────────────

@pytest.mark.parametrize("sr", JET_SWEEP_RATIOS)
def test_jet_ic_shape_and_energy(sr):
    ic = quasar_jet_ic(N=N, strain_ratio=sr)
    for k in ("u", "v", "w"):
        assert ic[k].shape == (N, N, N)
    assert ic["E0_analytic"] > 1e-12


def test_jet_ic_energy_is_fixed_across_the_sweep():
    """The swept parameter must be geometry, not amplitude: E0 must not vary."""
    es = [quasar_jet_ic(N=N, strain_ratio=sr)["E0_analytic"] for sr in JET_SWEEP_RATIOS]
    for e in es:
        assert e == pytest.approx(es[0], rel=1e-9)


def test_jet_ic_axial_component_is_bipolar():
    """u_z must be odd about the mid-plane z = pi: outflow both directions."""
    ic = quasar_jet_ic(N=N, strain_ratio=1.0, tilt=0.0)
    w = ic["w"]
    mid = N // 2
    above = w[:, :, mid + 2].mean()
    below = w[:, :, mid - 2].mean()
    assert above * below < 0 or (abs(above) < 1e-12 and abs(below) < 1e-12)


def test_jet_ic_zero_strain_has_no_axial_flow():
    ic = quasar_jet_ic(N=N, strain_ratio=0.0, tilt=0.0)
    assert np.max(np.abs(ic["w"])) < 1e-12


def test_jet_ic_strain_increases_axial_share():
    """Higher strain_ratio puts a larger share of the (fixed) energy in u_z."""
    def share(sr):
        ic = quasar_jet_ic(N=N, strain_ratio=sr, tilt=0.0)
        tot = np.mean(ic["u"]**2 + ic["v"]**2 + ic["w"]**2)
        return float(np.mean(ic["w"]**2) / tot)
    assert share(4.0) > share(1.0) > share(0.25)


def test_jet_ic_tilt_breaks_alignment():
    """tilt=0 gives a near-degenerate direction field; tilt>0 must change it."""
    a = quasar_jet_ic(N=N, strain_ratio=1.0, tilt=0.0)
    b = quasar_jet_ic(N=N, strain_ratio=1.0, tilt=0.2)
    assert not np.allclose(a["u"], b["u"])


def test_jet_ic_negative_strain_raises():
    with pytest.raises(ValueError, match="strain_ratio"):
        quasar_jet_ic(N=N, strain_ratio=-1.0)


def test_jet_ic_bad_core_frac_raises():
    with pytest.raises(ValueError, match="core_frac"):
        quasar_jet_ic(N=N, strain_ratio=1.0, core_frac=0.0)


def test_jet_ic_metadata():
    ic = quasar_jet_ic(N=N, strain_ratio=2.0)
    assert ic["name"] == "quasar_jet"
    assert ic["strain_ratio"] == 2.0
    assert ic["N"] == N


# ── runner + L-regression ────────────────────────────────────────────────────

def test_run_scale_experiment_smoke():
    r = run_scale_experiment(ic="tg", N=16, n_steps=6, record_every=3)
    assert r["verdict"] in ("PASS", "FAIL")
    assert r["n_samples"] >= 1
    for k in ("frac_E_rstar_median", "c_K_rstar_median", "c_K_band_median",
              "L_median", "key_metric", "exp_id", "claim_id"):
        assert k in r


def test_run_scale_experiment_with_jet_ic_override():
    ic = quasar_jet_ic(N=16, strain_ratio=1.0)
    r = run_scale_experiment(ic="jet", N=16, n_steps=6, record_every=3,
                             ic_dict_override=ic)
    assert r["n_samples"] >= 1
    assert np.isfinite(r["c_K_band_median"])


def test_run_scale_experiment_invalid_ic_raises():
    with pytest.raises(ValueError, match="Unknown ic"):
        run_scale_experiment(ic="bogus", N=16, n_steps=4)


def test_l_dependence_summary_recovers_a_planted_slope():
    """Synthetic samples with c_K = 3 * L^{-1} must regress to slope ~ -1."""
    Ls = np.linspace(1.5, 12.0, 40)
    fake = [{"timeseries": [
        {"L": float(L), "c_K_rstar": 3.0 / L, "c_K_band": 2.0 / (L**0.5)}
        for L in Ls]}]
    s = l_dependence_summary(fake)
    assert s["c_K_vs_L"]["slope"] == pytest.approx(-1.0, abs=1e-6)
    assert s["c_K_band_vs_L"]["slope"] == pytest.approx(-0.5, abs=1e-6)
    assert s["c_K_vs_L"]["r2"] == pytest.approx(1.0, abs=1e-6)


def test_l_dependence_summary_flat_case_gives_zero_slope():
    Ls = np.linspace(1.5, 12.0, 40)
    fake = [{"timeseries": [
        {"L": float(L), "c_K_rstar": 1.7, "c_K_band": 1.1} for L in Ls]}]
    s = l_dependence_summary(fake)
    assert s["c_K_vs_L"]["slope"] == pytest.approx(0.0, abs=1e-9)


def test_l_dependence_summary_drops_degenerate_samples():
    fake = [{"timeseries": [
        {"L": 2.0, "c_K_rstar": 0.0, "c_K_band": 0.0},
        {"L": 3.0, "c_K_rstar": 1.0, "c_K_band": 1.0},
        {"L": 4.0, "c_K_rstar": 1.0, "c_K_band": 1.0},
        {"L": 5.0, "c_K_rstar": 1.0, "c_K_band": 1.0},
    ]}]
    s = l_dependence_summary(fake)
    assert s["c_K_vs_L"]["n_dropped"] == 1
    assert s["c_K_vs_L"]["n_used"] == 3
