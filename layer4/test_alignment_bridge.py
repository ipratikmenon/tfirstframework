"""
layer4/test_alignment_bridge.py
Tests for the Bridge Lemma numerical diagnostic (layer4/alignment_bridge.py).

Run: python -m pytest layer4/test_alignment_bridge.py -v
"""

from __future__ import annotations

import os
import sqlite3
import tempfile

import numpy as np
import pytest

from layer4.alignment_bridge import (
    compute_direction_field,
    grad_ehat_sq,
    bridge_diagnostics,
    run_bridge_experiment,
    log_result,
    _vorticity_from_velocity,
    adv_vortex_tubes_ic,
)
from layer4.geometric_disorder import spectral_wavenumbers

import route2_3D as r2  # noqa: E402 — layer3 is on sys.path via alignment_bridge import above


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures / field constructors
# ─────────────────────────────────────────────────────────────────────────────

N = 16
KX, KY, KZ = spectral_wavenumbers(N)


def _solid_body_omega(N, direction=(0.0, 0.0, 1.0), amplitude=1.0):
    """ω = A·d everywhere (constant direction): ê constant, ∇ê = 0 exactly."""
    d = np.array(direction, dtype=float)
    d /= np.linalg.norm(d)
    wx = np.full((N, N, N), amplitude * d[0])
    wy = np.full((N, N, N), amplitude * d[1])
    wz = np.full((N, N, N), amplitude * d[2])
    return wx, wy, wz


def _twist_field(N, lam=1.0):
    """ω = λ²·(cos(λz), sin(λz), 0); |ω| = λ² everywhere; |∇ê|² = λ² everywhere."""
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    _, _, Z = np.meshgrid(x, x, x, indexing="ij")
    wx = (lam ** 2) * np.cos(lam * Z)
    wy = (lam ** 2) * np.sin(lam * Z)
    wz = np.zeros((N, N, N))
    return wx, wy, wz


def _low_vorticity_patch_omega(N):
    """Mostly ω=(0,0,1); a small corner is exactly zero -> masked / excluded."""
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = np.ones((N, N, N))
    wz[0:2, 0:2, 0:2] = 0.0
    return wx, wy, wz


# ─────────────────────────────────────────────────────────────────────────────
# compute_direction_field
# ─────────────────────────────────────────────────────────────────────────────

def test_direction_field_shapes_and_dtypes():
    wx, wy, wz = _solid_body_omega(N)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    assert ex.shape == wx.shape
    assert ey.shape == wx.shape
    assert ez.shape == wx.shape
    assert mag.shape == wx.shape
    assert valid.shape == wx.shape
    assert valid.dtype == bool


def test_direction_field_unit_norm_where_valid():
    wx, wy, wz = _twist_field(N, lam=1.0)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    norm = np.sqrt(ex ** 2 + ey ** 2 + ez ** 2)
    assert np.all(valid)
    assert np.max(np.abs(norm[valid] - 1.0)) < 1e-10


def test_direction_field_matches_omega_over_magnitude_solid_body():
    wx, wy, wz = _solid_body_omega(N, direction=(0, 0, 1), amplitude=3.0)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    assert np.allclose(ex, 0.0)
    assert np.allclose(ey, 0.0)
    assert np.allclose(ez, 1.0)
    assert np.allclose(mag, 3.0)


def test_direction_field_masks_zero_region():
    wx, wy, wz = _low_vorticity_patch_omega(N)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    # The exactly-zero corner must be excluded
    assert not np.any(valid[0:2, 0:2, 0:2])
    # ê is forced to zero there
    assert np.allclose(ex[0:2, 0:2, 0:2], 0.0)
    assert np.allclose(ey[0:2, 0:2, 0:2], 0.0)
    assert np.allclose(ez[0:2, 0:2, 0:2], 0.0)
    # elsewhere it's valid and normalized
    assert np.all(valid[2:, 2:, 2:])


def test_direction_field_all_zero_field_is_fully_invalid():
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = np.zeros((N, N, N))
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    assert not np.any(valid)
    assert np.allclose(ex, 0.0) and np.allclose(ey, 0.0) and np.allclose(ez, 0.0)


def test_direction_field_immutability():
    wx, wy, wz = _twist_field(N, lam=1.0)
    wx0, wy0, wz0 = wx.copy(), wy.copy(), wz.copy()
    compute_direction_field(wx, wy, wz)
    assert np.array_equal(wx, wx0)
    assert np.array_equal(wy, wy0)
    assert np.array_equal(wz, wz0)


# ─────────────────────────────────────────────────────────────────────────────
# grad_ehat_sq
# ─────────────────────────────────────────────────────────────────────────────

def test_grad_ehat_sq_zero_for_constant_direction():
    wx, wy, wz = _solid_body_omega(N, direction=(1, 0, 0), amplitude=2.0)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    g = grad_ehat_sq(ex, ey, ez, KX, KY, KZ)
    assert np.max(np.abs(g)) < 1e-20


def test_grad_ehat_sq_twist_field_exact_value():
    wx, wy, wz = _twist_field(N, lam=1.0)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    g = grad_ehat_sq(ex, ey, ez, KX, KY, KZ)
    # Analytic: e_hat = (cos z, sin z, 0); d/dz e_hat = (-sin z, cos z, 0)
    # |grad e_hat|^2 = 1 everywhere, exactly (spectrally exact single-mode field)
    assert np.max(np.abs(g - 1.0)) < 1e-10


def test_grad_ehat_sq_shape_matches_input():
    wx, wy, wz = _twist_field(N, lam=1.0)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    g = grad_ehat_sq(ex, ey, ez, KX, KY, KZ)
    assert g.shape == (N, N, N)


def test_grad_ehat_sq_nonnegative():
    wx, wy, wz = _low_vorticity_patch_omega(N)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    g = grad_ehat_sq(ex, ey, ez, KX, KY, KZ)
    assert np.all(g >= -1e-12)


def test_grad_ehat_sq_immutability():
    wx, wy, wz = _twist_field(N, lam=1.0)
    ex, ey, ez, mag, valid = compute_direction_field(wx, wy, wz)
    ex0, ey0, ez0 = ex.copy(), ey.copy(), ez.copy()
    grad_ehat_sq(ex, ey, ez, KX, KY, KZ)
    assert np.array_equal(ex, ex0)
    assert np.array_equal(ey, ey0)
    assert np.array_equal(ez, ez0)


# ─────────────────────────────────────────────────────────────────────────────
# bridge_diagnostics — basic contract
# ─────────────────────────────────────────────────────────────────────────────

def test_bridge_diagnostics_returns_all_keys():
    wx, wy, wz = _twist_field(N, lam=1.0)
    d = bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    for key in ("M", "r_star", "x_star", "sigma_star", "linf_grad_ehat", "bridge_ratio"):
        assert key in d


def test_bridge_diagnostics_all_finite():
    wx, wy, wz = _twist_field(N, lam=1.0)
    d = bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert np.isfinite(d["M"])
    assert np.isfinite(d["r_star"])
    assert np.isfinite(d["sigma_star"])
    assert np.isfinite(d["linf_grad_ehat"])
    assert np.isfinite(d["bridge_ratio"])


def test_bridge_diagnostics_r_star_formula():
    wx, wy, wz = _solid_body_omega(N, amplitude=4.0)
    d = bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert d["M"] == pytest.approx(4.0, abs=1e-10)
    assert d["r_star"] == pytest.approx(4.0 ** -0.5, rel=1e-10)


def test_bridge_diagnostics_x_star_is_argmax():
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = np.full((N, N, N), 0.1)
    wz[5, 6, 7] = 9.0  # a single hot spot
    d = bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert d["x_star"] == (5, 6, 7)
    assert d["M"] == pytest.approx(9.0, abs=1e-10)


def test_bridge_diagnostics_solid_body_sigma_and_linf_are_zero():
    wx, wy, wz = _solid_body_omega(N, direction=(0, 1, 0), amplitude=5.0)
    d = bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert d["sigma_star"] == pytest.approx(0.0, abs=1e-20)
    assert d["linf_grad_ehat"] == pytest.approx(0.0, abs=1e-20)
    # bridge_ratio = 0 / max(0, 1e-10) = 0
    assert d["bridge_ratio"] == pytest.approx(0.0, abs=1e-6)


def test_bridge_diagnostics_twist_field_positive_sigma_and_linf():
    wx, wy, wz = _twist_field(N, lam=1.0)
    d = bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert d["sigma_star"] > 0.0
    # |grad e_hat|^2 = 1 exactly, so L_inf sqrt(|grad e_hat|^2) = 1
    assert d["linf_grad_ehat"] == pytest.approx(1.0, abs=1e-6)


def test_bridge_diagnostics_bridge_ratio_matches_definition():
    wx, wy, wz = _twist_field(N, lam=1.0)
    d = bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    expected = (d["linf_grad_ehat"] * d["r_star"]) / max(np.sqrt(d["sigma_star"]), 1e-10)
    assert d["bridge_ratio"] == pytest.approx(expected, rel=1e-10)


def test_bridge_diagnostics_masking_excludes_low_vorticity_region():
    wx, wy, wz = _low_vorticity_patch_omega(N)
    d = bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert np.isfinite(d["sigma_star"])
    assert np.isfinite(d["linf_grad_ehat"])
    assert np.isfinite(d["bridge_ratio"])


# ─────────────────────────────────────────────────────────────────────────────
# bridge_diagnostics — scale invariance (numerical, with documented grid caveat)
# ─────────────────────────────────────────────────────────────────────────────

def test_bridge_diagnostics_scale_invariance_of_sigma_star():
    """sigma_star = A_loc/r* is a scale-invariant (self-similar) quantity:
    under omega -> lambda^2 * omega(lambda x), sigma_star is unchanged in the
    continuum limit.  On a finite grid with the *same* N for both fields
    (lambda=2 shrinks the physical ball radius r* by a factor of 2, so it is
    resolved by fewer grid points), agreement is only approximate.  We use a
    resolution (N=96) where empirically the discretization error is well
    under 5% (observed ~0.2%), and document this as a grid-dependent caveat
    rather than an exact identity.
    """
    N_hi = 96
    wx1, wy1, wz1 = _twist_field(N_hi, lam=1.0)
    d1 = bridge_diagnostics(wx1, wy1, wz1, 2.0 * np.pi / N_hi)

    wx2, wy2, wz2 = _twist_field(N_hi, lam=2.0)
    d2 = bridge_diagnostics(wx2, wy2, wz2, 2.0 * np.pi / N_hi)

    rel_diff = abs(d2["sigma_star"] - d1["sigma_star"]) / d1["sigma_star"]
    assert rel_diff < 0.05, (
        f"sigma_star not scale-invariant within grid tolerance: "
        f"lambda=1 -> {d1['sigma_star']}, lambda=2 -> {d2['sigma_star']}, "
        f"rel_diff={rel_diff:.4f}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# bridge_diagnostics — immutability
# ─────────────────────────────────────────────────────────────────────────────

def test_bridge_diagnostics_immutability():
    wx, wy, wz = _twist_field(N, lam=1.0)
    wx0, wy0, wz0 = wx.copy(), wy.copy(), wz.copy()
    bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)
    assert np.array_equal(wx, wx0)
    assert np.array_equal(wy, wy0)
    assert np.array_equal(wz, wz0)


# ─────────────────────────────────────────────────────────────────────────────
# bridge_diagnostics — error handling
# ─────────────────────────────────────────────────────────────────────────────

def test_bridge_diagnostics_all_zero_raises_valueerror():
    wx = np.zeros((N, N, N))
    wy = np.zeros((N, N, N))
    wz = np.zeros((N, N, N))
    with pytest.raises(ValueError, match="identically zero"):
        bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)


def test_bridge_diagnostics_shape_mismatch_raises():
    wx, wy, wz = _twist_field(N, lam=1.0)
    with pytest.raises(ValueError):
        bridge_diagnostics(wx, wy[:-1], wz, 2.0 * np.pi / N)


def test_bridge_diagnostics_non_cubic_grid_raises():
    wx = np.ones((N, N, N + 1))
    wy = np.ones((N, N, N + 1))
    wz = np.ones((N, N, N + 1))
    with pytest.raises(ValueError):
        bridge_diagnostics(wx, wy, wz, 2.0 * np.pi / N)


def test_bridge_diagnostics_bad_dx_raises():
    wx, wy, wz = _twist_field(N, lam=1.0)
    with pytest.raises(ValueError, match="dx"):
        bridge_diagnostics(wx, wy, wz, dx=1.0)  # inconsistent with N


# ─────────────────────────────────────────────────────────────────────────────
# _vorticity_from_velocity (curl helper used by run_bridge_experiment)
# ─────────────────────────────────────────────────────────────────────────────

def test_vorticity_from_velocity_immutability():
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    u = np.sin(xx) * np.cos(yy) * np.cos(zz)
    v = -np.cos(xx) * np.sin(yy) * np.cos(zz)
    w = np.zeros((N, N, N))
    u0, v0, w0 = u.copy(), v.copy(), w.copy()
    _vorticity_from_velocity(u, v, w, KX, KY, KZ)
    assert np.array_equal(u, u0)
    assert np.array_equal(v, v0)
    assert np.array_equal(w, w0)


def test_vorticity_from_velocity_finite_and_shaped():
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    u = np.sin(xx) * np.cos(yy) * np.cos(zz)
    v = -np.cos(xx) * np.sin(yy) * np.cos(zz)
    w = np.zeros((N, N, N))
    wx, wy, wz = _vorticity_from_velocity(u, v, w, KX, KY, KZ)
    assert wx.shape == (N, N, N)
    assert np.all(np.isfinite(wx))
    assert np.all(np.isfinite(wy))
    assert np.all(np.isfinite(wz))


# ─────────────────────────────────────────────────────────────────────────────
# run_bridge_experiment — integration smoke tests
# ─────────────────────────────────────────────────────────────────────────────

def test_run_bridge_experiment_tg_smoke():
    r = run_bridge_experiment("tg", N=16, n_steps=20, record_every=5, seed=1)
    assert r["verdict"] in ("PASS", "FAIL")
    assert np.isfinite(r["bridge_ratio_max"])
    assert np.isfinite(r["bridge_ratio_median"])
    assert np.isfinite(r["bridge_ratio_final"])
    assert len(r["timeseries"]) > 0
    assert r["seed"] == 1


def test_run_bridge_experiment_shear_smoke():
    r = run_bridge_experiment("shear", N=16, n_steps=20, record_every=5, seed=1)
    assert r["verdict"] in ("PASS", "FAIL")
    assert np.isfinite(r["bridge_ratio_max"])
    assert len(r["timeseries"]) > 0


def test_run_bridge_experiment_timeseries_entries_have_expected_keys():
    r = run_bridge_experiment("tg", N=16, n_steps=15, record_every=5, seed=2)
    for rec in r["timeseries"]:
        for key in ("step", "t", "M", "r_star", "sigma_star", "linf_grad_ehat", "bridge_ratio"):
            assert key in rec


def test_run_bridge_experiment_exp_id_and_claim_id():
    r_tg = run_bridge_experiment("tg", N=16, n_steps=10, record_every=5, seed=3)
    r_sh = run_bridge_experiment("shear", N=16, n_steps=10, record_every=5, seed=3)
    assert r_tg["exp_id"] == "EXP-L4-BRIDGE-TG-001"
    assert r_tg["claim_id"] == "bridge-lemma-tg"
    assert r_sh["exp_id"] == "EXP-L4-BRIDGE-SH-001"
    assert r_sh["claim_id"] == "bridge-lemma-shear"


def test_run_bridge_experiment_invalid_ic_raises():
    with pytest.raises(ValueError):
        run_bridge_experiment(ic="bogus", N=16, n_steps=10)


def test_run_bridge_experiment_invalid_n_steps_raises():
    with pytest.raises(ValueError):
        run_bridge_experiment(ic="tg", N=16, n_steps=0)


# ─────────────────────────────────────────────────────────────────────────────
# log_result — SQLite logging
# ─────────────────────────────────────────────────────────────────────────────

def test_log_result_writes_experiment_row():
    r = run_bridge_experiment("tg", N=16, n_steps=10, record_every=5, seed=4)
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_results.db")
        log_result(db_path, r)

        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT * FROM layer4_bridge_experiments WHERE exp_id = ?",
            (r["exp_id"],),
        ).fetchone()
        conn.close()

        assert row is not None
        assert row["claim_id"] == r["claim_id"]
        assert row["verdict"] == r["verdict"]
        assert row["N"] == 16


def test_log_result_writes_timeseries_rows():
    r = run_bridge_experiment("tg", N=16, n_steps=10, record_every=5, seed=5)
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_results.db")
        log_result(db_path, r)

        conn = sqlite3.connect(db_path)
        count = conn.execute(
            "SELECT COUNT(*) FROM layer4_bridge_timeseries WHERE exp_id = ?",
            (r["exp_id"],),
        ).fetchone()[0]
        conn.close()

        assert count == len(r["timeseries"])


def test_run_bridge_experiment_db_path_logs_automatically():
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_results.db")
        r = run_bridge_experiment(
            "shear", N=16, n_steps=10, record_every=5, seed=6, db_path=db_path
        )
        assert os.path.exists(db_path)
        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT exp_id FROM layer4_bridge_experiments WHERE exp_id = ?",
            (r["exp_id"],),
        ).fetchone()
        conn.close()
        assert row is not None


# ─────────────────────────────────────────────────────────────────────────────
# adv_vortex_tubes_ic — adversarial anti-parallel vortex-tube IC (S33)
# ─────────────────────────────────────────────────────────────────────────────

_N_ADV = 16


def test_adv_vortex_tubes_ic_returns_expected_keys_and_shapes():
    ic = adv_vortex_tubes_ic(_N_ADV, seed=42)
    for key in ("u", "v", "w", "name"):
        assert key in ic
    assert ic["u"].shape == (_N_ADV, _N_ADV, _N_ADV)
    assert ic["v"].shape == (_N_ADV, _N_ADV, _N_ADV)
    assert ic["w"].shape == (_N_ADV, _N_ADV, _N_ADV)
    assert ic["name"] == "adv_vortex_tubes"


def test_adv_vortex_tubes_ic_finite():
    ic = adv_vortex_tubes_ic(_N_ADV, seed=42)
    assert np.all(np.isfinite(ic["u"]))
    assert np.all(np.isfinite(ic["v"]))
    assert np.all(np.isfinite(ic["w"]))


def test_adv_vortex_tubes_ic_is_divergence_free_after_solver_projection():
    """The IC is fed through r2.set_ic (Leray projection) exactly as the
    'tg'/'shear' ICs are, before any run starts.  After that projection,
    the spectral RMS divergence must be < 1e-10.
    """
    solver = r2.make_solver(N=_N_ADV)
    ic = adv_vortex_tubes_ic(_N_ADV, seed=42)
    solver = r2.set_ic(solver, ic)
    div_rms = r2.divergence_rms(
        solver["u"], solver["v"], solver["w"],
        solver["kx"], solver["ky"], solver["kz"],
    )
    assert div_rms < 1e-10


def test_adv_vortex_tubes_ic_has_nonzero_enstrophy():
    ic = adv_vortex_tubes_ic(_N_ADV, seed=42)
    kx, ky, kz = spectral_wavenumbers(_N_ADV)
    wx, wy, wz = _vorticity_from_velocity(ic["u"], ic["v"], ic["w"], kx, ky, kz)
    enstrophy = float(np.mean(wx**2 + wy**2 + wz**2))
    assert enstrophy > 0.0
    assert np.isfinite(enstrophy)


def test_adv_vortex_tubes_ic_reproducible_for_same_seed():
    ic1 = adv_vortex_tubes_ic(_N_ADV, seed=7)
    ic2 = adv_vortex_tubes_ic(_N_ADV, seed=7)
    assert np.array_equal(ic1["u"], ic2["u"])
    assert np.array_equal(ic1["v"], ic2["v"])
    assert np.array_equal(ic1["w"], ic2["w"])


def test_adv_vortex_tubes_ic_differs_across_seeds():
    # Note: the seeded symmetry-breaking noise is injected only into omega_x
    # (make_vortex_omega), and via Biot-Savart's k-cross-product structure
    # (cross_x = ky*oz_hat - kz*oy_hat) omega_x feeds only v_hat and w_hat,
    # not u_hat — so u is seed-independent by construction; check v instead.
    ic1 = adv_vortex_tubes_ic(_N_ADV, seed=7)
    ic2 = adv_vortex_tubes_ic(_N_ADV, seed=8)
    assert not np.array_equal(ic1["v"], ic2["v"])


def test_run_bridge_experiment_adv_smoke_finite_diagnostics():
    r = run_bridge_experiment("adv", N=16, n_steps=10, record_every=5, seed=42)
    assert r["verdict"] in ("PASS", "FAIL")
    assert np.isfinite(r["bridge_ratio_max"])
    assert np.isfinite(r["bridge_ratio_median"])
    assert np.isfinite(r["bridge_ratio_final"])
    assert len(r["timeseries"]) > 0
    for rec in r["timeseries"]:
        assert np.isfinite(rec["sigma_star"])
        assert np.isfinite(rec["bridge_ratio"])


def test_run_bridge_experiment_adv_exp_id_and_claim_id_default():
    r = run_bridge_experiment("adv", N=16, n_steps=10, record_every=5, seed=42)
    assert r["exp_id"] == "EXP-L4-BRIDGE-ADV-001"
    assert r["claim_id"] == "bridge-lemma-adv"


def test_run_bridge_experiment_exp_id_and_claim_id_override():
    r = run_bridge_experiment(
        "tg", N=16, n_steps=10, record_every=5, seed=42,
        exp_id_override="EXP-L4-BRIDGE-TG-002",
        claim_id_override="bridge-lemma-tg-64",
    )
    assert r["exp_id"] == "EXP-L4-BRIDGE-TG-002"
    assert r["claim_id"] == "bridge-lemma-tg-64"


def test_run_bridge_experiment_adv_ic_still_rejects_bogus():
    with pytest.raises(ValueError):
        run_bridge_experiment(ic="not_a_real_ic", N=16, n_steps=10)
