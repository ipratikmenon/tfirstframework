"""
test_blowup_search_3D.py — Test suite for M9 Blowup Search 3D.

Coverage:
  - make_vortex_omega: shape, anti-parallel property, Gaussian structure
  - biot_savart_spectral: div-free output, energy positive, ω reconstruction
  - make_adversarial_ic_jax: shape, T field, div u ≈ 0, Z0 positive
  - compute_enstrophy_jax / compute_palinstrophy_jax: positive for non-trivial fields
  - check_blowup: correct detection logic
  - run_blowup_search: short N=16 run — A_min>0, verdict PASS, all keys present
  - EXPERIMENT_IDS: completeness and format
"""

from __future__ import annotations

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import jax.numpy as jnp
from layer3.blowup_search_3D import (
    make_vortex_omega,
    biot_savart_spectral,
    make_adversarial_ic_jax,
    compute_enstrophy_jax,
    compute_palinstrophy_jax,
    compute_u_linf_jax,
    check_blowup,
    run_blowup_search,
    EXPERIMENT_IDS,
    _T_ADV_DEFAULT,
    _AMP_DEFAULT,
    _TUBE_SIGMA,
    _TUBE_SEP,
    _d,
)
from layer3.route1_3D_jax import (
    make_grid_jax,
    make_dealias_mask_jax,
    make_nyquist_mask_jax,
    divergence_rms_jax,
)

N_SMALL = 16


# ── helpers ──────────────────────────────────────────────────────────────────

def _make_kxyz_np(N):
    freqs = np.fft.fftfreq(N, d=1.0 / N).astype(np.float32)
    kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    return kx, ky, kz, k2


# ── TestMakeVortexOmega ───────────────────────────────────────────────────────

class TestMakeVortexOmega:
    def test_shapes(self):
        ox, oy, oz = make_vortex_omega(N_SMALL)
        assert ox.shape == (N_SMALL, N_SMALL, N_SMALL)
        assert oy.shape == (N_SMALL, N_SMALL, N_SMALL)
        assert oz.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_dtype_float32(self):
        ox, oy, oz = make_vortex_omega(N_SMALL)
        assert ox.dtype == np.float32
        assert oz.dtype == np.float32

    def test_omega_z_antisymmetric(self):
        """omega_z should be zero-mean (tube1 - tube2 integrates to 0)."""
        _, _, oz = make_vortex_omega(N_SMALL)
        assert abs(float(np.mean(oz))) < 0.5, \
            f"mean(oz) = {np.mean(oz):.4f} should be near 0"

    def test_omega_z_has_both_signs(self):
        """Anti-parallel: oz has both positive and negative values."""
        _, _, oz = make_vortex_omega(N_SMALL)
        assert float(np.max(oz)) > 0.0
        assert float(np.min(oz)) < 0.0

    def test_amplitude_bounded(self):
        amp = 3.0
        ox, oy, oz = make_vortex_omega(N_SMALL, amp=amp)
        assert float(np.max(np.abs(oz))) <= amp + 1e-3, \
            f"max|oz|={np.max(np.abs(oz)):.4f} > amp={amp}"

    def test_reproducible(self):
        ox1, _, oz1 = make_vortex_omega(N_SMALL, seed=7)
        ox2, _, oz2 = make_vortex_omega(N_SMALL, seed=7)
        np.testing.assert_array_equal(oz1, oz2)
        np.testing.assert_array_equal(ox1, ox2)

    def test_different_seeds_differ(self):
        ox1, _, _ = make_vortex_omega(N_SMALL, seed=0)
        ox2, _, _ = make_vortex_omega(N_SMALL, seed=99)
        assert not np.allclose(ox1, ox2), "Different seeds must give different noise"

    def test_omega_y_is_zero(self):
        """ω_y component should be zero (only z-axis vorticity + x-noise)."""
        _, oy, _ = make_vortex_omega(N_SMALL)
        np.testing.assert_array_equal(oy, 0.0)


# ── TestBiotSavart ────────────────────────────────────────────────────────────

class TestBiotSavart:
    @pytest.fixture(scope="class")
    def bs_velocity(self):
        """Biot-Savart velocity for a standard vorticity field."""
        N = N_SMALL
        kx, ky, kz, k2 = _make_kxyz_np(N)
        ox, oy, oz = make_vortex_omega(N)
        u, v, w = biot_savart_spectral(ox, oy, oz, kx, ky, kz, k2)
        return u, v, w, kx, ky, kz, k2

    def test_shapes(self, bs_velocity):
        u, v, w, *_ = bs_velocity
        assert u.shape == (N_SMALL, N_SMALL, N_SMALL)
        assert v.shape == (N_SMALL, N_SMALL, N_SMALL)
        assert w.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_dtype_float32(self, bs_velocity):
        u, v, w, *_ = bs_velocity
        assert u.dtype == np.float32
        assert v.dtype == np.float32

    def test_divergence_free(self, bs_velocity):
        """div u should be near zero after Biot-Savart recovery."""
        u, v, w, kx, ky, kz, k2 = bs_velocity
        u_hat = np.fft.fftn(u)
        v_hat = np.fft.fftn(v)
        w_hat = np.fft.fftn(w)
        div_hat = 1j * (kx * u_hat + ky * v_hat + kz * w_hat)
        div_rms = float(np.sqrt(np.mean(np.real(np.fft.ifftn(div_hat))**2)))
        assert div_rms < 1e-4, f"div u RMS = {div_rms:.3e} too large"

    def test_velocity_not_zero(self, bs_velocity):
        u, v, w, *_ = bs_velocity
        E = 0.5 * np.mean(u**2 + v**2 + w**2)
        assert E > 0.0, "Biot-Savart must produce non-zero velocity from non-trivial vorticity"

    def test_zero_vorticity_gives_zero_velocity(self):
        N = N_SMALL
        kx, ky, kz, k2 = _make_kxyz_np(N)
        ox = np.zeros((N, N, N), dtype=np.float32)
        oy = np.zeros((N, N, N), dtype=np.float32)
        oz = np.zeros((N, N, N), dtype=np.float32)
        u, v, w = biot_savart_spectral(ox, oy, oz, kx, ky, kz, k2)
        np.testing.assert_allclose(u, 0.0, atol=1e-6)
        np.testing.assert_allclose(v, 0.0, atol=1e-6)
        np.testing.assert_allclose(w, 0.0, atol=1e-6)


# ── TestAdversarialIC ─────────────────────────────────────────────────────────

class TestAdversarialIC:
    @pytest.fixture(scope="class")
    def adv_ic(self):
        return make_adversarial_ic_jax(N_SMALL, T_base=_T_ADV_DEFAULT, seed=42)

    def test_shapes(self, adv_ic):
        for key in ["u", "v", "w", "T"]:
            assert adv_ic[key].shape == (N_SMALL, N_SMALL, N_SMALL), \
                f"{key} shape mismatch"

    def test_T_field_uniform(self, adv_ic):
        """T₀ should be uniform at T_base."""
        T = np.array(adv_ic["T"])
        np.testing.assert_allclose(T, _T_ADV_DEFAULT, rtol=1e-4)

    def test_T_field_positive(self, adv_ic):
        assert float(jnp.min(adv_ic["T"])) > 0.0

    def test_Z0_positive(self, adv_ic):
        assert adv_ic["Z0"] > 0.0, "Initial enstrophy must be > 0"

    def test_velocity_div_free(self, adv_ic):
        """set_ic_jax applies Leray projection — div u must be small."""
        kx, ky, kz, k2 = make_grid_jax(N_SMALL)
        div = divergence_rms_jax(adv_ic["u"], adv_ic["v"], adv_ic["w"], kx, ky, kz)
        assert div < 1e-4, f"div u = {div:.3e} after Leray projection"

    def test_energy_positive(self, adv_ic):
        E = float(0.5 * jnp.mean(
            adv_ic["u"]**2 + adv_ic["v"]**2 + adv_ic["w"]**2
        ))
        assert E > 0.0, "Adversarial IC must have non-zero kinetic energy"

    def test_name(self, adv_ic):
        assert adv_ic["name"] == "adversarial_antitube"

    def test_metadata_keys(self, adv_ic):
        for key in ["N", "T_base", "amp", "sigma", "separation", "Z0"]:
            assert key in adv_ic, f"Missing key: {key}"


# ── TestEnstrophy ─────────────────────────────────────────────────────────────

class TestEnstrophy:
    @pytest.fixture(scope="class")
    def grid(self):
        return make_grid_jax(N_SMALL)

    def test_enstrophy_positive(self, grid):
        kx, ky, kz, k2 = grid
        adv = make_adversarial_ic_jax(N_SMALL)
        Z = compute_enstrophy_jax(adv["u"], adv["v"], adv["w"], kx, ky, kz)
        assert Z > 0.0, f"Z = {Z:.4e} must be > 0"

    def test_enstrophy_zero_velocity(self, grid):
        kx, ky, kz, k2 = grid
        u = _d(jnp.zeros((N_SMALL, N_SMALL, N_SMALL)))
        v = _d(jnp.zeros((N_SMALL, N_SMALL, N_SMALL)))
        w = _d(jnp.zeros((N_SMALL, N_SMALL, N_SMALL)))
        Z = compute_enstrophy_jax(u, v, w, kx, ky, kz)
        assert abs(Z) < 1e-10, f"Z={Z:.4e} should be 0 for zero velocity"

    def test_palinstrophy_positive(self, grid):
        kx, ky, kz, k2 = grid
        adv = make_adversarial_ic_jax(N_SMALL)
        P = compute_palinstrophy_jax(adv["u"], adv["v"], adv["w"], kx, ky, kz, k2)
        assert P >= 0.0, f"P = {P:.4e} must be ≥ 0"

    def test_u_linf_positive(self):
        adv = make_adversarial_ic_jax(N_SMALL)
        linf = compute_u_linf_jax(adv["u"], adv["v"], adv["w"])
        assert linf > 0.0


# ── TestBlowupDetector ────────────────────────────────────────────────────────

class TestBlowupDetector:
    def test_no_blowup_at_threshold(self):
        """Z = factor * Z0 exactly: NOT a blowup (strict inequality)."""
        Z0 = 1.0
        factor = 10.0
        assert not check_blowup(Z=factor * Z0, Z0=Z0, factor=factor)

    def test_blowup_detected_above_threshold(self):
        Z0 = 1.0
        factor = 10.0
        assert check_blowup(Z=factor * Z0 + 1e-6, Z0=Z0, factor=factor)

    def test_no_blowup_below_threshold(self):
        Z0 = 1.0
        factor = 20.0
        assert not check_blowup(Z=5.0 * Z0, Z0=Z0, factor=factor)

    def test_zero_Z0_never_triggers(self):
        """Zero initial enstrophy: blowup check should never fire."""
        assert not check_blowup(Z=1e10, Z0=0.0)

    def test_negative_Z0_never_triggers(self):
        assert not check_blowup(Z=1e10, Z0=-1.0)

    def test_custom_factor(self):
        Z0 = 2.0
        assert check_blowup(Z=101.0, Z0=Z0, factor=50.0)
        assert not check_blowup(Z=99.0, Z0=Z0, factor=50.0)


# ── TestRunBlowupSearch ───────────────────────────────────────────────────────

class TestRunBlowupSearch:
    """Short runs at N=16 for CI speed."""

    @pytest.fixture(scope="class")
    def result_n16(self):
        return run_blowup_search(
            N=N_SMALL, T_base=_T_ADV_DEFAULT, t_end=0.05,
            max_steps=10, verbose=False,
        )

    def test_verdict_pass(self, result_n16):
        assert result_n16["verdict"] == "PASS", \
            f"verdict={result_n16['verdict']}, A_min={result_n16['A_min_global']:.3e}"

    def test_A_min_positive(self, result_n16):
        assert result_n16["A_min_global"] > 0.0, \
            f"A_min_global={result_n16['A_min_global']:.3e} must be > 0"

    def test_blowup_not_detected(self, result_n16):
        assert not result_n16["blowup_detected"], \
            "Blowup should not be detected for physical T>0"

    def test_required_keys(self, result_n16):
        required = [
            "exp_id", "verdict", "A_min_global", "Z0", "Z_max",
            "P_max", "u_linf_max", "blowup_detected", "n_steps",
            "t_final", "wall_time_s", "key_metric", "history",
        ]
        for k in required:
            assert k in result_n16, f"Missing key: {k}"

    def test_Z_max_finite(self, result_n16):
        assert np.isfinite(result_n16["Z_max"]), "Z_max must be finite"
        assert result_n16["Z_max"] > 0.0

    def test_history_non_empty(self, result_n16):
        assert len(result_n16["history"]) > 0

    def test_history_keys(self, result_n16):
        for entry in result_n16["history"]:
            for k in ["t", "dt", "A_min", "Z", "P", "u_linf"]:
                assert k in entry, f"history entry missing key: {k}"

    def test_exp_id_format(self, result_n16):
        assert result_n16["exp_id"].startswith("EXP-L3-R1-ADV-")

    def test_t_final_positive(self, result_n16):
        assert result_n16["t_final"] > 0.0

    def test_room_temperature_also_passes(self):
        """T_base=300K (standard room temp) should also pass."""
        r = run_blowup_search(
            N=N_SMALL, T_base=300.0, t_end=0.05,
            max_steps=10, verbose=False,
        )
        assert r["verdict"] == "PASS"
        assert r["A_min_global"] > 0.0


# ── TestExperimentIDs ─────────────────────────────────────────────────────────

class TestExperimentIDs:
    def test_has_three_resolutions(self):
        assert len(EXPERIMENT_IDS) == 3

    def test_expected_keys(self):
        for N in [64, 128, 256]:
            assert N in EXPERIMENT_IDS, f"N={N} missing from EXPERIMENT_IDS"

    def test_id_format(self):
        for N, eid in EXPERIMENT_IDS.items():
            assert eid.startswith("EXP-L3-R1-ADV-"), f"{eid} has wrong prefix"
            # Should contain the N value
            assert str(N) in eid, f"{eid} does not contain N={N}"

    def test_all_unique(self):
        ids = list(EXPERIMENT_IDS.values())
        assert len(set(ids)) == len(ids), "Experiment IDs must be unique"
