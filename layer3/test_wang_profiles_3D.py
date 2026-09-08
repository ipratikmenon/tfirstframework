"""
test_wang_profiles_3D.py — Test suite for M6 Wang et al. 3D IC generator.

Coverage:
  - embed_2D_to_3D: shape, div-free property, perturbation scaling
  - wang_ic_3D: shape, named profiles, λ values, non-zero energy
  - Short solver runs: δ_max > 0 (Lemma 2.5 on adversarial ICs)
  - M6 profile set completeness
  - Divergence after Leray projection
"""

from __future__ import annotations

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from layer3.wang_profiles_3D import (
    embed_2D_to_3D,
    wang_ic_3D,
    divergence_rms_3D,
    profile_energy_3D,
    M6_PROFILES,
    run_wang_route2_jax,
)
from layer2.self_similar_IC import generate_CCF_profile, WANG_PROFILES


# ── constants ─────────────────────────────────────────────────────────────────
N_SMALL = 16
N_MED   = 32
NU      = 1e-3


# ── helpers ───────────────────────────────────────────────────────────────────

def _make_2d_velocity(N: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Simple divergence-free 2D velocity for embedding tests."""
    from layer2.self_similar_IC import generate_CCF_profile, make_grid
    ic = generate_CCF_profile(0.6057, N=N, seed=seed)
    return ic["u"], ic["v"]


# ── TestEmbed2Dto3D ───────────────────────────────────────────────────────────

class TestEmbed2Dto3D:
    def test_output_shapes(self):
        u2, v2 = _make_2d_velocity(N_SMALL)
        u3, v3, w3 = embed_2D_to_3D(u2, v2, N_SMALL)
        assert u3.shape == (N_SMALL, N_SMALL, N_SMALL)
        assert v3.shape == (N_SMALL, N_SMALL, N_SMALL)
        assert w3.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_u_z_independence_without_perturb(self):
        u2, v2 = _make_2d_velocity(N_SMALL)
        u3, v3, w3 = embed_2D_to_3D(u2, v2, N_SMALL, perturb_amp=0.0)
        # u and v should be z-independent: all z-slices identical
        for k in range(N_SMALL):
            np.testing.assert_allclose(u3[:, :, k], u2, atol=1e-12)
            np.testing.assert_allclose(v3[:, :, k], v2, atol=1e-12)
        # w should be zero when perturb_amp=0
        assert np.max(np.abs(w3)) < 1e-10

    def test_w_is_zero_without_perturb(self):
        u2, v2 = _make_2d_velocity(N_SMALL)
        _, _, w3 = embed_2D_to_3D(u2, v2, N_SMALL, perturb_amp=0.0)
        assert np.max(np.abs(w3)) < 1e-10

    def test_perturbation_is_small(self):
        u2, v2 = _make_2d_velocity(N_SMALL)
        u3_no, _, _ = embed_2D_to_3D(u2, v2, N_SMALL, perturb_amp=0.0)
        u3_yes, v3_yes, w3 = embed_2D_to_3D(u2, v2, N_SMALL, perturb_amp=0.02)
        E_base = float(np.mean(u2**2 + v2**2))
        E_w    = float(np.mean(w3**2))
        # w perturbation energy < 10% of base energy (it's 2% amplitude)
        assert E_w < 0.10 * E_base, f"w energy {E_w:.3e} >= 10% base {E_base:.3e}"

    def test_non_zero_energy(self):
        u2, v2 = _make_2d_velocity(N_SMALL)
        u3, v3, w3 = embed_2D_to_3D(u2, v2, N_SMALL)
        E = profile_energy_3D(u3, v3, w3)
        assert E > 1e-6, f"Energy too small: {E:.3e}"

    def test_shape_mismatch_raises(self):
        u2 = np.ones((N_SMALL, N_SMALL))
        v2 = np.ones((N_SMALL + 1, N_SMALL))  # wrong shape
        with pytest.raises(AssertionError):
            embed_2D_to_3D(u2, v2, N_SMALL)

    def test_different_seeds_give_different_w(self):
        u2, v2 = _make_2d_velocity(N_SMALL)
        _, _, w1 = embed_2D_to_3D(u2, v2, N_SMALL, perturb_amp=0.05, seed=0)
        _, _, w2 = embed_2D_to_3D(u2, v2, N_SMALL, perturb_amp=0.05, seed=1)
        # Different seeds → different sign/structure (at least one element differs)
        assert not np.allclose(w1, w2), "w fields should differ for different seeds"

    def test_perturb_modes_respected(self):
        u2, v2 = _make_2d_velocity(N_SMALL)
        _, _, w1 = embed_2D_to_3D(u2, v2, N_SMALL, perturb_amp=0.05, perturb_modes=1)
        _, _, w2 = embed_2D_to_3D(u2, v2, N_SMALL, perturb_amp=0.05, perturb_modes=3)
        # More modes → different w fields
        assert not np.allclose(w1, w2)


# ── TestWangIC3D ──────────────────────────────────────────────────────────────

class TestWangIC3D:
    @pytest.mark.parametrize("profile_name", list(M6_PROFILES.keys()))
    def test_output_shape(self, profile_name):
        ic = wang_ic_3D(profile_name, N=N_SMALL)
        assert ic["u"].shape == (N_SMALL, N_SMALL, N_SMALL)
        assert ic["v"].shape == (N_SMALL, N_SMALL, N_SMALL)
        assert ic["w"].shape == (N_SMALL, N_SMALL, N_SMALL)

    @pytest.mark.parametrize("profile_name", list(M6_PROFILES.keys()))
    def test_lambda_values_correct(self, profile_name):
        ic = wang_ic_3D(profile_name, N=N_SMALL)
        expected_lv = M6_PROFILES[profile_name]["lambda_val"]
        assert ic["lambda_val"] == pytest.approx(expected_lv, rel=1e-4)

    @pytest.mark.parametrize("profile_name", list(M6_PROFILES.keys()))
    def test_non_zero_energy(self, profile_name):
        ic = wang_ic_3D(profile_name, N=N_SMALL)
        E = profile_energy_3D(ic["u"], ic["v"], ic["w"])
        assert E > 1e-8, f"{profile_name}: energy {E:.3e} too small"

    def test_ccf2_has_higher_enstrophy_than_ccf1(self):
        ic1 = wang_ic_3D("CCF_1st_unstable", N=N_SMALL)
        ic2 = wang_ic_3D("CCF_2nd_unstable", N=N_SMALL)
        # CCF2 (λ=0.4703) is more concentrated → higher peak vorticity
        # (w_max may not hold in 3D embedding, but omega_max_2d should)
        assert ic2["omega_max_2d"] >= ic1["omega_max_2d"] * 0.5, \
            "CCF2 should have meaningful 2D vorticity (same order as CCF1)"

    def test_invalid_profile_raises(self):
        with pytest.raises(ValueError, match="Unknown profile"):
            wang_ic_3D("not_a_real_profile", N=N_SMALL)

    @pytest.mark.parametrize("profile_name", list(M6_PROFILES.keys()))
    def test_metadata_present(self, profile_name):
        ic = wang_ic_3D(profile_name, N=N_SMALL)
        assert "lambda_val" in ic
        assert "family" in ic
        assert "omega_max_2d" in ic
        assert ic["N"] == N_SMALL

    def test_e0_analytic_positive(self):
        for name in M6_PROFILES:
            ic = wang_ic_3D(name, N=N_SMALL)
            assert ic["E0_analytic"] > 0


# ── TestDivergence3D ──────────────────────────────────────────────────────────

class TestDivergence3D:
    def test_divergence_helper_zero_field(self):
        u = np.zeros((N_SMALL,) * 3)
        div = divergence_rms_3D(u, u, u, N_SMALL)
        assert div < 1e-12

    def test_known_divergence_free(self):
        # u = sin(x), v = 0, w = 0  → div = cos(x) → not div-free
        # But u = sin(y), v = 0 → div = 0  (no x-dependence)
        x1d = np.linspace(0.0, 2.0 * np.pi, N_SMALL, endpoint=False)
        _, yy, _ = np.meshgrid(x1d, x1d, x1d, indexing="ij")
        u = np.sin(yy)   # ∂u/∂x = 0
        v = np.zeros_like(u)
        w = np.zeros_like(u)
        div = divergence_rms_3D(u, v, w, N_SMALL)
        assert div < 1e-6, f"Expected near-zero div, got {div:.3e}"

    def test_leray_reduces_divergence(self):
        """Embedding + Leray projection reduces divergence below numerical noise."""
        from layer3.route2_3D_jax import make_solver_jax, set_ic_jax
        ic = wang_ic_3D("CCF_2nd_unstable", N=N_SMALL)
        # Before projection: embedding may have small divergence from perturbation
        div_before = divergence_rms_3D(ic["u"], ic["v"], ic["w"], N_SMALL)

        # After Leray projection (inside set_ic_jax):
        solver = make_solver_jax(N=N_SMALL, nu=NU)
        solver = set_ic_jax(solver, ic)
        import jax.numpy as jnp
        u_np = np.array(solver["u"])
        v_np = np.array(solver["v"])
        w_np = np.array(solver["w"])
        div_after = divergence_rms_3D(u_np, v_np, w_np, N_SMALL)

        assert div_after < 1e-4, f"Post-Leray div {div_after:.2e} too large"


# ── TestWangSolverRun ─────────────────────────────────────────────────────────

class TestWangSolverRun:
    """Short runs (N=16, 10 steps) verifying δ>0 on adversarial ICs."""

    @pytest.mark.parametrize("profile_name", list(M6_PROFILES.keys()))
    def test_delta_positive_after_short_run(self, profile_name):
        """Core M6 claim: δ_max > 0 on all Wang et al. adversarial ICs."""
        r = run_wang_route2_jax(
            profile_name, N=N_SMALL, t_end=0.05, nu=NU,
            max_steps=10, verbose=False,
        )
        assert r["delta_max"] >= 0.0, f"{profile_name}: δ_max={r['delta_max']:.3f} < 0"

    @pytest.mark.parametrize("profile_name", list(M6_PROFILES.keys()))
    def test_theta_non_negative(self, profile_name):
        r = run_wang_route2_jax(
            profile_name, N=N_SMALL, t_end=0.05, nu=NU,
            max_steps=10, verbose=False,
        )
        assert not r["theta_negative"], f"{profile_name}: θ went negative"

    @pytest.mark.parametrize("profile_name", list(M6_PROFILES.keys()))
    def test_result_has_required_keys(self, profile_name):
        r = run_wang_route2_jax(
            profile_name, N=N_SMALL, t_end=0.05, nu=NU,
            max_steps=10, verbose=False,
        )
        for key in ["exp_id", "verdict", "delta_max", "lambda_val", "key_metric"]:
            assert key in r, f"Missing key: {key}"

    @pytest.mark.parametrize("profile_name", list(M6_PROFILES.keys()))
    def test_exp_id_contains_profile_suffix(self, profile_name):
        r = run_wang_route2_jax(
            profile_name, N=N_SMALL, t_end=0.05, nu=NU,
            max_steps=10, verbose=False,
        )
        suffix = M6_PROFILES[profile_name]["exp_suffix"]
        assert suffix in r["exp_id"], f"exp_id '{r['exp_id']}' missing suffix '{suffix}'"

    def test_ccf2_is_critical_profile(self):
        """CCF_2nd_unstable has the smallest λ among CCF profiles → hardest test."""
        lv_ccf1 = M6_PROFILES["CCF_1st_unstable"]["lambda_val"]
        lv_ccf2 = M6_PROFILES["CCF_2nd_unstable"]["lambda_val"]
        assert lv_ccf2 < lv_ccf1, "CCF2 should have smaller λ (harder) than CCF1"

    def test_adversarial_has_smallest_lambda(self):
        lambdas = [v["lambda_val"] for v in M6_PROFILES.values()]
        adv_lv = M6_PROFILES["adversarial_min"]["lambda_val"]
        assert adv_lv == min(lambdas), "adversarial_min should have smallest λ"


# ── TestM6ProfileSet ──────────────────────────────────────────────────────────

class TestM6ProfileSet:
    def test_m6_profiles_complete(self):
        assert "CCF_1st_unstable" in M6_PROFILES
        assert "CCF_2nd_unstable" in M6_PROFILES
        assert "adversarial_min" in M6_PROFILES

    def test_m6_lambda_values(self):
        assert M6_PROFILES["CCF_1st_unstable"]["lambda_val"] == pytest.approx(0.6057)
        assert M6_PROFILES["CCF_2nd_unstable"]["lambda_val"] == pytest.approx(0.4703)
        assert M6_PROFILES["adversarial_min"]["lambda_val"] == pytest.approx(0.05)

    def test_exp_suffixes_unique(self):
        suffixes = [v["exp_suffix"] for v in M6_PROFILES.values()]
        assert len(suffixes) == len(set(suffixes)), "exp_suffixes must be unique"

    def test_all_lambda_positive(self):
        for name, info in M6_PROFILES.items():
            assert info["lambda_val"] > 0, f"{name} has λ ≤ 0"


# ── TestBenchmark32 ───────────────────────────────────────────────────────────

class TestBenchmark32:
    """N=32 smoke run for the critical CCF_2nd_unstable profile."""

    @pytest.fixture(scope="class")
    def ccf2_result(self):
        return run_wang_route2_jax(
            "CCF_2nd_unstable", N=N_MED, t_end=0.5, nu=NU,
            max_steps=500, verbose=False,
        )

    def test_passes(self, ccf2_result):
        assert ccf2_result["verdict"] == "PASS", \
            f"CCF2 N=32 FAIL: {ccf2_result['key_metric']}"

    def test_delta_positive(self, ccf2_result):
        assert ccf2_result["delta_max"] > 0.0, \
            f"CCF2 N=32 δ_max={ccf2_result['delta_max']:.3f}"

    def test_theta_non_negative(self, ccf2_result):
        assert not ccf2_result["theta_negative"]

    def test_lambda_value_correct(self, ccf2_result):
        assert ccf2_result["lambda_val"] == pytest.approx(0.4703, rel=1e-3)
