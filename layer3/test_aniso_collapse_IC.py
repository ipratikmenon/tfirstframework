"""
test_aniso_collapse_IC.py — Test suite for the anisotropic collapsing-vortex
adversarial IC generator (layer3/aniso_collapse_IC.py).

Coverage:
  - aniso_collapse_ic_3D: shape, energy calibration, aspect-ratio scaling
  - Divergence after Leray projection
  - Short solver runs: delta_max > 0 on the adversarial IC (exact unforced
    equations, f = 0 throughout — route2_3D_jax.run_jax carries no forcing
    term at all)
"""

from __future__ import annotations

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from layer3.aniso_collapse_IC import (
    aniso_collapse_ic_3D,
    run_aniso_route2_jax,
    divergence_rms_3D,
    ANISO_SWEEP_VALUES,
)

N_SMALL = 16
NU = 1e-3


# ── TestAnisoCollapseIC3D ────────────────────────────────────────────────────

class TestAnisoCollapseIC3D:
    @pytest.mark.parametrize("aspect_ratio", ANISO_SWEEP_VALUES)
    def test_output_shape(self, aspect_ratio):
        ic = aniso_collapse_ic_3D(N=N_SMALL, aspect_ratio=aspect_ratio)
        assert ic["u"].shape == (N_SMALL, N_SMALL, N_SMALL)
        assert ic["v"].shape == (N_SMALL, N_SMALL, N_SMALL)
        assert ic["w"].shape == (N_SMALL, N_SMALL, N_SMALL)

    @pytest.mark.parametrize("aspect_ratio", ANISO_SWEEP_VALUES)
    def test_non_zero_energy(self, aspect_ratio):
        ic = aniso_collapse_ic_3D(N=N_SMALL, aspect_ratio=aspect_ratio)
        assert ic["E0_analytic"] > 1e-8

    def test_energy_matches_reference_across_aspect_ratios(self):
        """Gamma calibration should equalize E0 across aspect ratios (same
        adversarial severity by kinetic energy, only geometry differs)."""
        energies = [
            aniso_collapse_ic_3D(N=N_SMALL, aspect_ratio=ar)["E0_analytic"]
            for ar in ANISO_SWEEP_VALUES
        ]
        for E in energies:
            assert E == pytest.approx(energies[0], rel=1e-6)

    def test_ellz_scales_with_aspect_ratio(self):
        ic_1 = aniso_collapse_ic_3D(N=64, aspect_ratio=1.0)
        ic_2 = aniso_collapse_ic_3D(N=64, aspect_ratio=0.1)
        assert ic_2["ellz"] < ic_1["ellz"]
        assert ic_1["ellr"] == pytest.approx(ic_2["ellr"])

    def test_ellz_floor_at_grid_resolution(self):
        """ellz should never collapse below one grid cell."""
        ic = aniso_collapse_ic_3D(N=32, aspect_ratio=1e-6)
        assert ic["ellz"] >= (2.0 * np.pi / 32) - 1e-12

    def test_invalid_aspect_ratio_raises(self):
        with pytest.raises(ValueError, match="aspect_ratio"):
            aniso_collapse_ic_3D(N=N_SMALL, aspect_ratio=0.0)
        with pytest.raises(ValueError, match="aspect_ratio"):
            aniso_collapse_ic_3D(N=N_SMALL, aspect_ratio=-0.5)

    def test_metadata_present(self):
        ic = aniso_collapse_ic_3D(N=N_SMALL, aspect_ratio=0.5)
        assert ic["name"] == "aniso_collapse"
        assert ic["aspect_ratio"] == 0.5
        assert ic["N"] == N_SMALL
        assert "Gamma" in ic


# ── TestDivergence ────────────────────────────────────────────────────────────

class TestDivergence:
    @pytest.mark.parametrize("aspect_ratio", ANISO_SWEEP_VALUES)
    def test_leray_reduces_divergence(self, aspect_ratio):
        from layer3.route2_3D_jax import make_solver_jax, set_ic_jax

        ic = aniso_collapse_ic_3D(N=N_SMALL, aspect_ratio=aspect_ratio)
        solver = make_solver_jax(N=N_SMALL, nu=NU)
        solver = set_ic_jax(solver, ic)

        u_np = np.array(solver["u"])
        v_np = np.array(solver["v"])
        w_np = np.array(solver["w"])
        div_after = divergence_rms_3D(u_np, v_np, w_np, N_SMALL)

        assert div_after < 1e-3, f"Post-Leray div {div_after:.2e} too large"


# ── TestAnisoSolverRun ────────────────────────────────────────────────────────

class TestAnisoSolverRun:
    """Short runs (N=16, 10 steps) under the exact unforced Route 2
    equations (f = 0). Verifies delta >= 0 and no forcing artifacts."""

    @pytest.mark.parametrize("aspect_ratio", ANISO_SWEEP_VALUES)
    def test_delta_nonnegative_after_short_run(self, aspect_ratio):
        r = run_aniso_route2_jax(
            aspect_ratio, N=N_SMALL, t_end=0.05, nu=NU,
            max_steps=10, verbose=False,
        )
        assert r["delta_max"] >= 0.0, \
            f"aspect_ratio={aspect_ratio}: delta_max={r['delta_max']:.3f} < 0"

    @pytest.mark.parametrize("aspect_ratio", ANISO_SWEEP_VALUES)
    def test_theta_non_negative(self, aspect_ratio):
        r = run_aniso_route2_jax(
            aspect_ratio, N=N_SMALL, t_end=0.05, nu=NU,
            max_steps=10, verbose=False,
        )
        assert not r["theta_negative"], f"aspect_ratio={aspect_ratio}: theta went negative"

    @pytest.mark.parametrize("aspect_ratio", ANISO_SWEEP_VALUES)
    def test_result_has_required_keys(self, aspect_ratio):
        r = run_aniso_route2_jax(
            aspect_ratio, N=N_SMALL, t_end=0.05, nu=NU,
            max_steps=10, verbose=False,
        )
        for key in ("exp_id", "claim_id", "delta_max", "verdict", "key_metric"):
            assert key in r

    def test_exp_id_format(self):
        r = run_aniso_route2_jax(0.5, N=N_SMALL, t_end=0.02, nu=NU, max_steps=5)
        assert r["exp_id"].startswith("EXP-L3-R2-ANISO")
        assert r["exp_id"].endswith(f"-JAX-{N_SMALL:03d}")

    def test_claim_id_fixed(self):
        r = run_aniso_route2_jax(1.0, N=N_SMALL, t_end=0.02, nu=NU, max_steps=5)
        assert r["claim_id"] == "claim_l3_aniso_collapse"
