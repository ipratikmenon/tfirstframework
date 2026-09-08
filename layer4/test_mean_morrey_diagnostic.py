"""Tests for mean_morrey_diagnostic.py

All tests use synthetic velocity fields — no NS solve required.
Run with: python -m pytest layer4/test_mean_morrey_diagnostic.py -v
"""
from __future__ import annotations

import sqlite3

import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from layer4.mean_morrey_diagnostic import (
    L_BOX,
    ball_mask,
    compute_mean_morrey_field,
    compute_mean_morrey_field_random,
    compute_spatial_mean,
    find_alpha_threshold,
    fit_morrey_scaling,
    init_db,
    log_result,
    run_mean_morrey_analysis,
)


# ─────────────────────────────────────────────────────────────────────────────
# Test 1: ball_mask correctness
# ─────────────────────────────────────────────────────────────────────────────

class TestBallMask:
    def test_symmetry_and_bounds(self):
        """Mask should be symmetric around centre, not all-True, not all-False."""
        N = 16
        cx, cy, cz = 8, 8, 8
        r_frac = 0.1   # small radius → well-contained ball

        mask = ball_mask(N, cx, cy, cz, r_frac)

        assert mask.shape == (N, N, N)
        assert mask.dtype == bool
        # Not degenerate
        assert mask.any(), "Expected at least one True voxel"
        assert not mask.all(), "Expected at least one False voxel"

    def test_centre_included(self):
        """Centre voxel must always be inside the ball."""
        N = 16
        mask = ball_mask(N, 8, 8, 8, r_frac=0.05)
        assert mask[8, 8, 8], "Centre voxel should be inside the ball"

    def test_all_true_in_ball_are_within_radius(self):
        """Every True voxel should be within r_frac * L_BOX of the centre."""
        N = 16
        cx, cy, cz = 8, 8, 8
        r_frac = 0.1
        r_phys = r_frac * L_BOX
        dx = L_BOX / N

        mask = ball_mask(N, cx, cy, cz, r_frac)
        idx = np.array(np.where(mask)).T   # shape (n_true, 3)

        for (ii, jj, kk) in idx:
            di = min(abs(ii - cx), N - abs(ii - cx)) * dx
            dj = min(abs(jj - cy), N - abs(jj - cy)) * dx
            dk = min(abs(kk - cz), N - abs(kk - cz)) * dx
            dist = np.sqrt(di**2 + dj**2 + dk**2)
            assert dist <= r_phys + 1e-10, (
                f"Voxel ({ii},{jj},{kk}) at dist {dist:.4f} exceeds r={r_phys:.4f}"
            )

    def test_periodic_distance(self):
        """A ball centred near the boundary should wrap around periodically."""
        N = 32
        # Centre at (0,0,0) — corner — with moderate radius
        mask_corner = ball_mask(N, 0, 0, 0, r_frac=0.15)
        # Centre at mid — for comparison
        mask_mid = ball_mask(N, N // 2, N // 2, N // 2, r_frac=0.15)
        # Both balls should have the same volume (periodic symmetry)
        assert mask_corner.sum() == mask_mid.sum(), (
            "Periodic ball should have same count as centred ball"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Test 2: Shear layer — mean is zero
# ─────────────────────────────────────────────────────────────────────────────

class TestShearLayerZeroMean:
    """
    Test with the zero velocity field u=0, which genuinely has zero spatial mean
    over any ball.  Note: a sinusoidal shear u=(sin(2πy/L),0,0) has zero mean
    over the full torus T³ but NOT over an arbitrary ball B_r(x₀) — the local
    mean depends on x₀ and r.  The NS shear-layer result (Thm 11.4) proves
    Claim A directly via parabolic regularity, not via mean Morrey.
    """
    def _make_zero_field(self, N: int) -> np.ndarray:
        """u = (0,0,0) — identically zero, mean is exactly zero everywhere."""
        return np.zeros((3, N, N, N))

    def test_spatial_mean_near_zero(self):
        """Spatial mean over any ball of the zero field is exactly 0."""
        N = 32
        u = self._make_zero_field(N)

        for cx, cy, cz in [(8, 8, 8), (16, 0, 16), (N // 2, N // 2, N // 2)]:
            _, mag = compute_spatial_mean(u, cx, cy, cz, r_frac=0.3)
            assert mag < 1e-14, (
                f"Zero field mean magnitude should be 0, got {mag:.2e}"
            )

    def test_run_analysis_morrey_near_zero(self):
        """run_mean_morrey_analysis with zero field: morrey_max should be 0."""
        N = 32
        u = self._make_zero_field(N)
        snapshots = [u, u, u]

        result = run_mean_morrey_analysis(
            exp_id="EXP-TEST-ZERO-001",
            u_snapshots=snapshots,
            t_values=[0.0, 0.1, 0.2],
            grid_N=N,
            ic_type="zero_field",
        )

        assert result["morrey_max"] < 1e-10, (
            f"Zero field morrey_max should be 0, got {result['morrey_max']:.2e}"
        )
        assert "verdict" in result


# ─────────────────────────────────────────────────────────────────────────────
# Test 3: Constant field — known scaling
# ─────────────────────────────────────────────────────────────────────────────

class TestConstantField:
    def _make_constant_field(self, N: int) -> np.ndarray:
        """u = (1, 0, 0) everywhere."""
        u = np.zeros((3, N, N, N))
        u[0] = 1.0
        return u

    def test_alpha_threshold_near_one(self):
        """For constant u=(1,0,0), alpha_threshold ≈ 1.0 (± 0.3 tolerance)."""
        N = 32
        u = self._make_constant_field(N)

        # Use a fine set of alphas straddling 1.0
        alpha_values = [0.1, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
        r_values = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2]

        mm, _ = compute_mean_morrey_field(
            u, r_values, alpha_values, n_centres=4
        )
        alpha_threshold, _ = find_alpha_threshold(r_values, mm, alpha_values)

        assert np.isfinite(alpha_threshold), "Should find a finite alpha_threshold"
        assert abs(alpha_threshold - 1.0) <= 0.3, (
            f"Expected alpha_threshold ≈ 1.0, got {alpha_threshold:.3f}"
        )

    def test_spatial_mean_magnitude_is_one(self):
        """Mean of constant u=(1,0,0) over any non-empty ball = 1."""
        N = 32
        u = self._make_constant_field(N)
        _, mag = compute_spatial_mean(u, N // 2, N // 2, N // 2, r_frac=0.2)
        assert abs(mag - 1.0) < 1e-10, f"Expected mean magnitude 1.0, got {mag:.6f}"


# ─────────────────────────────────────────────────────────────────────────────
# Test 4: TG-like divergence-free field
# ─────────────────────────────────────────────────────────────────────────────

class TestTaylorGreenField:
    def _make_tg_field(self, N: int, A: float = 1.0) -> np.ndarray:
        """2D Taylor-Green vortex extended to 3D: (A·sin(x)cos(y), -A·cos(x)sin(y), 0)."""
        x = np.linspace(0, L_BOX, N, endpoint=False)
        y = np.linspace(0, L_BOX, N, endpoint=False)
        xx, yy, _ = np.meshgrid(x, y, np.arange(N), indexing="ij")

        u = np.zeros((3, N, N, N))
        u[0] = A * np.sin(xx) * np.cos(yy)
        u[1] = -A * np.cos(xx) * np.sin(yy)
        # u[2] = 0 — divergence-free in 2D sense
        return u

    def test_large_ball_mean_small(self):
        """Spatial mean over a large ball of TG field should be small (< 0.1)."""
        N = 32
        u = self._make_tg_field(N, A=1.0)
        _, mag = compute_spatial_mean(u, N // 2, N // 2, N // 2, r_frac=0.4)
        assert mag < 0.1, (
            f"TG large-ball mean should be < 0.1 (zero-mean periodic), got {mag:.4f}"
        )

    def test_run_analysis_returns_verdict(self):
        """run_mean_morrey_analysis should return a dict with 'verdict' key."""
        N = 32
        u = self._make_tg_field(N)
        snapshots = [u]

        result = run_mean_morrey_analysis(
            exp_id="EXP-TEST-TG-001",
            u_snapshots=snapshots,
            t_values=[0.0],
            grid_N=N,
            ic_type="TG",
        )

        assert "verdict" in result
        assert result["verdict"] in ("PASS", "FAIL")
        assert "alpha_threshold" in result
        assert "morrey_max" in result
        assert "morrey_at_r_min" in result
        assert "scaling_exponent" in result


# ─────────────────────────────────────────────────────────────────────────────
# Test 5: Scaling test — constructed |⟨u⟩| ~ r^{1/2}
# ─────────────────────────────────────────────────────────────────────────────

class TestMorreyScalingFit:
    """
    Construct synthetic morrey_dict such that |⟨u⟩_{B_r}| = r^{1/2}.
    Then M_α(r) = r^{-1+α} * (r^{1/2})^3 = r^{-1+α+3/2} = r^{1/2+α}.
    So fit_morrey_scaling should return γ ≈ 1/2 + α.
    """

    def _make_synthetic_morrey_dict(
        self,
        r_values: list,
        alpha_values: list,
    ) -> dict:
        """Build morrey_dict where |mean| = r^{1/2} by construction."""
        morrey_dict = {}
        for r in r_values:
            mean_mag = r ** 0.5
            for alpha in alpha_values:
                exponent = -1.0 + alpha
                m_val = (r ** exponent) * (mean_mag ** 3)
                morrey_dict[(r, alpha)] = m_val
        return morrey_dict

    @pytest.mark.parametrize("alpha", [0.1, 0.5, 1.0, 1.5])
    def test_scaling_exponent(self, alpha):
        """fit_morrey_scaling should return γ ≈ 0.5 + alpha (tolerance 0.2)."""
        r_values = [0.05, 0.1, 0.2, 0.4, 0.8, 1.0, 1.5]
        alpha_values = [alpha]

        morrey_dict = self._make_synthetic_morrey_dict(r_values, alpha_values)
        C, gamma, R2 = fit_morrey_scaling(r_values, morrey_dict, alpha)

        expected_gamma = 0.5 + alpha
        assert abs(gamma - expected_gamma) < 0.2, (
            f"alpha={alpha}: expected γ ≈ {expected_gamma:.2f}, got {gamma:.4f}"
        )
        assert R2 > 0.95, f"Expected R² > 0.95 for exact power law, got {R2:.4f}"

    def test_positive_gamma_means_bounded(self):
        """γ ≥ 0 for alpha ≥ 0.5 when mean ~ r^{0.5} (since 0.5+0.5=1 > 0)."""
        r_values = [0.05, 0.1, 0.2, 0.4, 0.8]
        alpha_values = [0.5, 0.75, 1.0, 1.5, 2.0]
        morrey_dict = self._make_synthetic_morrey_dict(r_values, alpha_values)

        for alpha in alpha_values:
            _, gamma, _ = fit_morrey_scaling(r_values, morrey_dict, alpha)
            assert gamma >= 0.0, (
                f"alpha={alpha}: expected γ ≥ 0 (bounded), got {gamma:.4f}"
            )

    def test_negative_gamma_for_constant_mean(self):
        """
        For a CONSTANT mean (|⟨u⟩| = 1 independent of r):
          M_α(r) = r^{-1+α} * 1³ = r^{-1+α}  →  γ = -1+α.
        For α = 0.1: γ = -0.9 < 0 (diverges as r→0, seminorm unbounded).
        This uses a separate synthetic dict with constant mean magnitude.
        """
        r_values = [0.05, 0.1, 0.2, 0.4, 0.8]
        alpha = 0.1
        # Build morrey_dict with constant mean = 1
        morrey_dict = {(r, alpha): r ** (-1.0 + alpha) * 1.0 for r in r_values}
        _, gamma, _ = fit_morrey_scaling(r_values, morrey_dict, alpha)
        assert gamma < 0.0, (
            f"alpha=0.1 with mean~r^0.5: expected γ < 0, got {gamma:.4f}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Test 6: DB logging
# ─────────────────────────────────────────────────────────────────────────────

class TestDBLogging:
    def _make_result(self) -> dict:
        return {
            "exp_id":           "EXP-TEST-DB-001",
            "ic_type":          "TG",
            "grid_N":           32,
            "alpha_threshold":  1.0,
            "morrey_max":       3.14,
            "morrey_at_r_min":  0.5,
            "scaling_exponent": 1.5,
            "verdict":          "PASS",
            "notes":            "unit test row",
            "timestamp":        "2026-05-02T12:00:00+00:00",
        }

    def test_init_db_creates_table(self):
        """init_db should create mean_morrey_results table in-memory."""
        conn = sqlite3.connect(":memory:")
        init_db(conn)

        cursor = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='mean_morrey_results'"
        )
        assert cursor.fetchone() is not None, "Table mean_morrey_results should exist"
        conn.close()

    def test_log_result_inserts_row(self):
        """log_result should insert exactly one row."""
        conn = sqlite3.connect(":memory:")
        result = self._make_result()
        log_result(conn, result)

        cursor = conn.execute("SELECT COUNT(*) FROM mean_morrey_results")
        count = cursor.fetchone()[0]
        assert count == 1, f"Expected 1 row after log_result, got {count}"
        conn.close()

    def test_log_result_values_round_trip(self):
        """Logged values should match what was inserted."""
        conn = sqlite3.connect(":memory:")
        result = self._make_result()
        log_result(conn, result)

        row = conn.execute(
            "SELECT exp_id, verdict, grid_N, alpha_threshold FROM mean_morrey_results"
        ).fetchone()

        assert row[0] == "EXP-TEST-DB-001"
        assert row[1] == "PASS"
        assert row[2] == 32
        assert abs(row[3] - 1.0) < 1e-10
        conn.close()

    def test_log_result_idempotent_multiple_rows(self):
        """Calling log_result twice should produce two distinct rows."""
        conn = sqlite3.connect(":memory:")
        result = self._make_result()
        log_result(conn, result)
        log_result(conn, result)

        cursor = conn.execute("SELECT COUNT(*) FROM mean_morrey_results")
        assert cursor.fetchone()[0] == 2
        conn.close()


# ─────────────────────────────────────────────────────────────────────────────
# Test 7: Random centre sampling
# ─────────────────────────────────────────────────────────────────────────────

class TestRandomCentres:
    """Tests for compute_mean_morrey_field_random."""

    def _make_tg_field(self, N: int = 32, A: float = 1.0) -> np.ndarray:
        """2D Taylor-Green vortex extended to 3D."""
        x = np.linspace(0, L_BOX, N, endpoint=False)
        xx, yy, _ = np.meshgrid(x, x, np.arange(N), indexing="ij")
        u = np.zeros((3, N, N, N))
        u[0] = A * np.sin(xx) * np.cos(yy)
        u[1] = -A * np.cos(xx) * np.sin(yy)
        return u

    _R_VALUES = [0.1, 0.2, 0.4, 0.8]
    _A_VALUES = [0.5, 1.0, 1.5]

    def test_random_centres_reproducible(self):
        """Same seed → identical morrey_dict output."""
        N = 32
        u = self._make_tg_field(N)
        mm1, centres1 = compute_mean_morrey_field_random(
            u, self._R_VALUES, self._A_VALUES, n_centres=16, seed=42
        )
        mm2, centres2 = compute_mean_morrey_field_random(
            u, self._R_VALUES, self._A_VALUES, n_centres=16, seed=42
        )
        assert np.array_equal(centres1, centres2), (
            "Same seed should produce identical centre arrays"
        )
        for key in mm1:
            assert abs(mm1[key] - mm2[key]) < 1e-14, (
                f"Same seed: morrey values differ at {key}"
            )

    def test_random_centres_different_seeds(self):
        """Different seeds may produce different centre samples (probabilistic)."""
        N = 32
        u = self._make_tg_field(N)
        _, centres_a = compute_mean_morrey_field_random(
            u, self._R_VALUES, self._A_VALUES, n_centres=16, seed=1
        )
        _, centres_b = compute_mean_morrey_field_random(
            u, self._R_VALUES, self._A_VALUES, n_centres=16, seed=99
        )
        # It would be astronomically unlikely for 16 random draws in [0,31]^3
        # to be identical with different seeds — assert they differ
        assert not np.array_equal(centres_a, centres_b), (
            "Different seeds should (almost certainly) yield different centres"
        )

    def test_random_max_geq_fixed(self):
        """Max over N random centres should be >= the value at a single centre.

        This is a monotonicity property: taking max over a superset of centres
        can only equal or exceed the single-centre result.
        """
        N = 32
        u = self._make_tg_field(N)
        # Fix one centre at (N//2, N//2, N//2) — the box centre
        cx, cy, cz = N // 2, N // 2, N // 2
        r_frac = 0.2
        r_phys = r_frac * L_BOX
        _, single_mag = compute_spatial_mean(u, cx, cy, cz, r_frac)
        single_cube = single_mag ** 3

        # Random sampling with many centres
        mm_rand, _ = compute_mean_morrey_field_random(
            u, [r_phys], [1.0], n_centres=64, seed=7
        )
        rand_val = mm_rand[(r_phys, 1.0)]

        # M_1.0(r) = r^0 * |mean|^3 = |mean|^3 (alpha=1 → exponent=0)
        assert rand_val >= single_cube - 1e-12, (
            f"Random max ({rand_val:.4e}) should be >= single-centre value ({single_cube:.4e})"
        )

    def test_zero_field_random_zero(self):
        """Random sampling of a zero velocity field → all morrey values are 0."""
        N = 32
        u = np.zeros((3, N, N, N))
        mm, centres = compute_mean_morrey_field_random(
            u, self._R_VALUES, self._A_VALUES, n_centres=32, seed=0
        )
        assert centres.shape == (32, 3), (
            f"centres_xyz should have shape (32,3), got {centres.shape}"
        )
        for key, val in mm.items():
            assert abs(val) < 1e-14, (
                f"Zero field: morrey value at {key} should be 0, got {val:.2e}"
            )
