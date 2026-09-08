"""
test_mu_sweep_3D.py — Test suite for M8 Phase 4 μ(T)→ν Limit 3D.

Coverage:
  - make_T_initial_3D: shape, range, L∞ norm, mean, reproducibility
  - compute_A_min: positive for physical T, zero-handling
  - run_single_delta: A_min > 0, verdict PASS, convergence direction
  - analyse_linearity: slope near 1 for linear δ scaling
  - DELTA_SWEEP_VALUES completeness
  - Short sweep (2 δ values at N=16) for integration test
"""

from __future__ import annotations

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from layer3.mu_sweep_3D import (
    make_T_initial_3D,
    compute_A_min,
    run_single_delta,
    analyse_linearity,
    DELTA_SWEEP_VALUES,
    _get_A_ref,
    _T_MEAN,
)


# ── constants ──────────────────────────────────────────────────────────────────
N_SMALL = 16
NU      = 1e-3


# ── TestMakeTInitial3D ────────────────────────────────────────────────────────

class TestMakeTInitial3D:
    def test_shape(self):
        T = make_T_initial_3D(N_SMALL, _T_MEAN, delta=1.0)
        assert T.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_dtype_float32(self):
        T = make_T_initial_3D(N_SMALL, _T_MEAN, delta=1.0)
        assert T.dtype == np.float32

    def test_mean_close_to_T_mean(self):
        T = make_T_initial_3D(N_SMALL, _T_MEAN, delta=1.0)
        assert abs(float(np.mean(T)) - _T_MEAN) < 1.0, \
            f"mean(T) = {np.mean(T):.2f} not close to T_mean={_T_MEAN}"

    def test_range_bounded_by_delta(self):
        delta = 5.0
        T = make_T_initial_3D(N_SMALL, _T_MEAN, delta=delta)
        assert float(np.min(T)) >= _T_MEAN - delta - 1e-4
        assert float(np.max(T)) <= _T_MEAN + delta + 1e-4

    def test_linf_norm_is_delta(self):
        delta = 3.0
        T = make_T_initial_3D(N_SMALL, _T_MEAN, delta=delta)
        g = T - _T_MEAN
        g_linf = float(np.max(np.abs(g)))
        assert abs(g_linf - delta) < delta * 0.05, \
            f"‖g‖_L∞ = {g_linf:.4f}, expected ≈ {delta:.4f}"

    def test_small_delta_gives_near_uniform(self):
        delta = 0.001
        T = make_T_initial_3D(N_SMALL, _T_MEAN, delta=delta)
        std = float(np.std(T))
        assert std < 0.01, f"std(T) = {std:.4f} too large for δ={delta}"

    def test_zero_delta_is_uniform(self):
        T = make_T_initial_3D(N_SMALL, _T_MEAN, delta=0.0)
        np.testing.assert_allclose(T, _T_MEAN, atol=1e-5)

    def test_reproducibility(self):
        T1 = make_T_initial_3D(N_SMALL, _T_MEAN, delta=2.0, seed=7)
        T2 = make_T_initial_3D(N_SMALL, _T_MEAN, delta=2.0, seed=7)
        np.testing.assert_array_equal(T1, T2)

    def test_different_seeds_different(self):
        T1 = make_T_initial_3D(N_SMALL, _T_MEAN, delta=2.0, seed=0)
        T2 = make_T_initial_3D(N_SMALL, _T_MEAN, delta=2.0, seed=1)
        assert not np.allclose(T1, T2), "Different seeds should give different fields"

    def test_T_always_positive(self):
        # Even for large δ, T = T_mean ± δ and T_mean >> δ for physical values
        delta = 50.0
        T = make_T_initial_3D(N_SMALL, _T_MEAN, delta=delta)
        assert float(np.min(T)) > 0.0, "T must remain positive"


# ── TestComputeAMin ──────────────────────────────────────────────────────────

class TestComputeAMin:
    def test_uniform_T_gives_A_ref(self):
        """A_min for uniform T = T_mean should equal A_ref."""
        import jax.numpy as jnp
        from layer3.mu_sweep_3D import _d
        T_uniform = _d(jnp.full((N_SMALL, N_SMALL, N_SMALL), _T_MEAN, dtype=jnp.float32))
        A_min = compute_A_min(T_uniform)
        A_ref = _get_A_ref()
        assert abs(A_min - A_ref) / A_ref < 1e-3, \
            f"A_min={A_min:.4e} should be close to A_ref={A_ref:.4e} for uniform T"

    def test_A_min_positive_for_physical_T(self):
        import jax.numpy as jnp
        from layer3.mu_sweep_3D import _d
        T = _d(jnp.full((N_SMALL, N_SMALL, N_SMALL), 300.0, dtype=jnp.float32))
        assert compute_A_min(T) > 0.0

    def test_A_min_positive_for_varying_T(self):
        import jax.numpy as jnp
        from layer3.mu_sweep_3D import _d
        T_np = make_T_initial_3D(N_SMALL, _T_MEAN, delta=50.0)
        T = _d(jnp.array(T_np, dtype=jnp.float32))
        assert compute_A_min(T) > 0.0, "A_min must be > 0 (second law)"

    def test_A_min_increases_as_T_increases(self):
        """Higher T → higher k(T) → larger A(T) for Sutherland ideal gas."""
        import jax.numpy as jnp
        from layer3.mu_sweep_3D import _d
        T_low  = _d(jnp.full((N_SMALL,) * 3, 200.0, dtype=jnp.float32))
        T_high = _d(jnp.full((N_SMALL,) * 3, 400.0, dtype=jnp.float32))
        A_low  = compute_A_min(T_low)
        A_high = compute_A_min(T_high)
        assert A_high > A_low, f"A_high={A_high:.4e} should > A_low={A_low:.4e}"


# ── TestARef ─────────────────────────────────────────────────────────────────

class TestARef:
    def test_A_ref_positive(self):
        A_ref = _get_A_ref()
        assert A_ref > 0.0, f"A_ref = {A_ref:.4e} must be > 0"

    def test_A_ref_matches_M4_value(self):
        """A_ref should be ≈ 3.08e-5 (same as M4 2D result)."""
        A_ref = _get_A_ref()
        assert 2e-5 < A_ref < 5e-5, f"A_ref={A_ref:.4e} out of expected range [2e-5, 5e-5]"

    def test_A_ref_cached(self):
        """Calling _get_A_ref twice returns the same value."""
        A1 = _get_A_ref()
        A2 = _get_A_ref()
        assert A1 == A2


# ── TestRunSingleDelta ────────────────────────────────────────────────────────

class TestRunSingleDelta:
    @pytest.mark.parametrize("delta", [1.0, 0.1, 0.01])
    def test_A_min_positive(self, delta):
        r = run_single_delta(delta, exp_index=1, N=N_SMALL, t_end=0.05,
                             max_steps=10, verbose=False)
        assert r["A_min_global"] > 0.0, \
            f"δ={delta}: A_min_global={r['A_min_global']:.3e} must be > 0"

    @pytest.mark.parametrize("delta", [1.0, 0.1, 0.01])
    def test_verdict_pass(self, delta):
        r = run_single_delta(delta, exp_index=1, N=N_SMALL, t_end=0.05,
                             max_steps=10, verbose=False)
        assert r["verdict"] == "PASS", f"δ={delta}: verdict={r['verdict']}"

    def test_exp_id_format(self):
        r = run_single_delta(1.0, exp_index=3, N=N_SMALL, t_end=0.05,
                             max_steps=5, verbose=False)
        assert r["exp_id"] == "EXP-L3-R1-MULIMIT-003"

    def test_smaller_delta_smaller_eps(self):
        """ε(δ) = |A_min - A_ref| should decrease as δ decreases."""
        r_large = run_single_delta(1.0, exp_index=1, N=N_SMALL, t_end=0.05,
                                   max_steps=10, verbose=False)
        r_small = run_single_delta(0.01, exp_index=2, N=N_SMALL, t_end=0.05,
                                   max_steps=10, verbose=False)
        eps_large = abs(r_large["eps_delta"])
        eps_small = abs(r_small["eps_delta"])
        assert eps_small <= eps_large + 1e-8, \
            f"ε(0.01)={eps_small:.4e} should ≤ ε(1.0)={eps_large:.4e}"

    def test_required_keys(self):
        r = run_single_delta(0.1, exp_index=1, N=N_SMALL, t_end=0.05,
                             max_steps=5, verbose=False)
        for key in ["exp_id", "A_ref", "A_min_global", "eps_delta", "verdict", "key_metric"]:
            assert key in r, f"Missing key: {key}"

    def test_A_ref_consistent(self):
        """A_ref returned by experiment should match _get_A_ref()."""
        r = run_single_delta(0.1, exp_index=1, N=N_SMALL, t_end=0.05,
                             max_steps=5, verbose=False)
        assert abs(r["A_ref"] - _get_A_ref()) < 1e-12


# ── TestAnalyseLinearity ──────────────────────────────────────────────────────

class TestAnalyseLinearity:
    def _make_linear_results(self, deltas, A_ref, alpha=1.0):
        """Construct synthetic results with |ε| ∝ δ^alpha."""
        return {
            d: {"delta": d, "eps_delta": A_ref * d**alpha * 0.01}
            for d in deltas
        }

    def test_linear_scaling_detected(self):
        deltas = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
        A_ref = _get_A_ref()
        results = self._make_linear_results(deltas, A_ref, alpha=1.0)
        fit = analyse_linearity(results)
        assert fit["fit_ok"], f"Expected linear fit, got α={fit['alpha']:.3f}"
        assert 0.8 < fit["alpha"] < 1.2, f"α={fit['alpha']:.3f} not near 1"

    def test_quadratic_not_linear(self):
        deltas = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
        A_ref = _get_A_ref()
        results = self._make_linear_results(deltas, A_ref, alpha=2.0)
        fit = analyse_linearity(results)
        # α should be ≈ 2 for quadratic scaling
        assert fit["alpha"] > 1.5, f"Expected α≈2, got {fit['alpha']:.3f}"

    def test_too_few_points_returns_nan(self):
        results = {0.1: {"eps_delta": 1e-5}}
        fit = analyse_linearity(results)
        assert not fit["fit_ok"]

    def test_zero_eps_excluded(self):
        results = {
            1.0 : {"eps_delta": 0.01},
            0.5 : {"eps_delta": 0.005},
            0.001: {"eps_delta": 0.0},   # zero excluded
        }
        fit = analyse_linearity(results)
        assert fit["n_pts"] == 2


# ── TestDeltaSweepValues ──────────────────────────────────────────────────────

class TestDeltaSweepValues:
    def test_ten_values(self):
        assert len(DELTA_SWEEP_VALUES) == 10

    def test_decreasing_order(self):
        for i in range(len(DELTA_SWEEP_VALUES) - 1):
            assert DELTA_SWEEP_VALUES[i] > DELTA_SWEEP_VALUES[i + 1], \
                "DELTA_SWEEP_VALUES must be in decreasing order"

    def test_includes_large_and_small(self):
        assert max(DELTA_SWEEP_VALUES) >= 1.0
        assert min(DELTA_SWEEP_VALUES) <= 0.002

    def test_all_positive(self):
        assert all(d > 0 for d in DELTA_SWEEP_VALUES)


# ── TestIntegration ───────────────────────────────────────────────────────────

class TestIntegration:
    """2-point sweep at N=16 — verifies end-to-end pipeline."""

    @pytest.fixture(scope="class")
    def short_sweep(self):
        from layer3.mu_sweep_3D import run_mu_sweep_3D
        return run_mu_sweep_3D(
            N=N_SMALL, t_end=0.05, max_steps=10,
            verbose=False, delta_values=[1.0, 0.01],
        )

    def test_all_pass(self, short_sweep):
        for delta, r in short_sweep.items():
            assert r["verdict"] == "PASS", f"δ={delta}: {r['verdict']}"

    def test_A_min_positive_all(self, short_sweep):
        for delta, r in short_sweep.items():
            assert r["A_min_global"] > 0.0, f"δ={delta}: A_min={r['A_min_global']:.3e}"

    def test_eps_smaller_at_smaller_delta(self, short_sweep):
        eps_large = abs(short_sweep[1.0]["eps_delta"])
        eps_small = abs(short_sweep[0.01]["eps_delta"])
        assert eps_small <= eps_large + 1e-9, \
            f"ε(0.01)={eps_small:.4e} not ≤ ε(1.0)={eps_large:.4e}"

    def test_linearity_fit(self, short_sweep):
        fit = analyse_linearity(short_sweep)
        assert fit["n_pts"] >= 2
