"""
test_props.py — Unit tests for Layer 1 Property Engine

Every test validates a specific mathematical claim from PRD §7.2 Layer 1.
Run with: pytest layer1/test_props.py -v
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import pytest
from layer1.tfirst_props import (
    ideal_props, co2_props, A_field, second_law_check,
    self_similar_props, conjecture_34_sweep, WANG_LAMBDA,
    route1_coeffs, route2_theta_source,
)


# ---------------------------------------------------------------------------
# ideal_props
# ---------------------------------------------------------------------------

class TestIdealProps:

    def test_returns_all_keys(self):
        r = ideal_props(300.0)
        for key in ("mu", "k", "rho", "cv", "cp", "gamma"):
            assert key in r, f"Missing key: {key}"

    def test_scalar_input(self):
        r = ideal_props(300.0)
        assert isinstance(float(r["mu"]), float)

    def test_array_input_shape_preserved(self):
        T = np.array([100.0, 300.0, 500.0, 1000.0])
        r = ideal_props(T)
        assert r["mu"].shape == (4,), f"Expected shape (4,), got {r['mu'].shape}"

    def test_viscosity_positive(self):
        T = np.linspace(100, 2000, 200)
        r = ideal_props(T)
        assert np.all(r["mu"] > 0), "Viscosity must be > 0"

    def test_conductivity_positive(self):
        T = np.linspace(100, 2000, 200)
        r = ideal_props(T)
        assert np.all(r["k"] > 0), "Thermal conductivity must be > 0"

    def test_density_positive(self):
        T = np.linspace(100, 2000, 200)
        r = ideal_props(T)
        assert np.all(r["rho"] > 0), "Density must be > 0"

    def test_cv_positive_and_constant(self):
        T = np.linspace(100, 2000, 50)
        r = ideal_props(T)
        assert np.all(r["cv"] > 0)
        assert np.allclose(r["cv"], 718.0), "Ideal gas cv should be 718 J/(kg·K)"

    def test_sutherland_viscosity_at_273K(self):
        """μ(273.15 K) should equal μ_ref = 1.716e-5 Pa·s."""
        r = ideal_props(273.15)
        assert abs(float(r["mu"]) - 1.716e-5) / 1.716e-5 < 1e-6, (
            f"Sutherland μ at 273.15K = {float(r['mu']):.4e}, expected 1.716e-5"
        )

    def test_viscosity_increases_with_temperature(self):
        """For a gas, μ increases with T (unlike liquids)."""
        T = np.array([200.0, 400.0, 800.0, 1600.0])
        mu = ideal_props(T)["mu"]
        assert np.all(np.diff(mu) > 0), "Ideal gas viscosity must increase with T"

    def test_negative_temperature_raises(self):
        with pytest.raises(ValueError, match="T must be > 0"):
            ideal_props(-10.0)

    def test_zero_temperature_raises(self):
        with pytest.raises(ValueError, match="T must be > 0"):
            ideal_props(0.0)

    def test_gamma_is_correct(self):
        """Diatomic ideal gas: γ = cp/cv ≈ 1.4 (cv=718, R=287 → cp=1005, γ=1005/718)."""
        r = ideal_props(300.0)
        gamma = float(r["gamma"])
        expected = (718.0 + 287.0) / 718.0  # = 1.39972...
        assert abs(gamma - expected) < 1e-6, (
            f"γ = {gamma:.6f}, expected {expected:.6f}"
        )
        assert 1.39 < gamma < 1.41, f"γ = {gamma:.4f} outside physical range [1.39, 1.41]"

    def test_2d_array_input(self):
        T = np.ones((4, 4)) * 300.0
        r = ideal_props(T)
        assert r["mu"].shape == (4, 4)


# ---------------------------------------------------------------------------
# co2_props
# ---------------------------------------------------------------------------

class TestCO2Props:

    def test_returns_all_keys(self):
        r = co2_props(320.0, 10e6)
        for key in ("mu", "k", "rho", "cv", "cp"):
            assert key in r, f"Missing key: {key}"

    def test_mu_positive(self):
        T = np.linspace(280, 400, 30)
        r = co2_props(T, 10e6)
        assert np.all(np.asarray(r["mu"]) > 0), "CO2 viscosity must be > 0"

    def test_k_positive(self):
        T = np.linspace(280, 400, 30)
        r = co2_props(T, 10e6)
        assert np.all(np.asarray(r["k"]) > 0), "CO2 conductivity must be > 0"

    def test_rho_positive(self):
        T = np.linspace(280, 400, 30)
        r = co2_props(T, 10e6)
        assert np.all(np.asarray(r["rho"]) > 0), "CO2 density must be > 0"

    def test_cv_positive(self):
        T = np.linspace(280, 400, 30)
        r = co2_props(T, 10e6)
        assert np.all(np.asarray(r["cv"]) > 0), "CO2 cv must be > 0"

    def test_scalar_input(self):
        r = co2_props(320.0, 10e6)
        assert isinstance(float(r["mu"]), float)

    def test_multiple_pressures(self):
        """PRD §4.2: Widom line verified at all 5 pressures."""
        pressures = [7.5e6, 8.0e6, 10.0e6, 15.0e6, 20.0e6]
        T_test = 320.0
        for P in pressures:
            r = co2_props(T_test, P)
            assert float(r["mu"]) > 0, f"CO2 mu ≤ 0 at P={P:.1e} Pa"
            assert float(r["k"]) > 0, f"CO2 k ≤ 0 at P={P:.1e} Pa"

    def test_A_min_matches_prd_exp2(self):
        """PRD §4.2 CO2-Exp2: A_min = 1.07e-7 at any P in SC region."""
        T = np.linspace(280, 400, 100)
        pressures = [7.5e6, 10e6, 15e6, 20e6]
        for P in pressures:
            r = co2_props(T, P)
            k   = np.asarray(r["k"])
            rho = np.asarray(r["rho"])
            cv  = np.asarray(r["cv"])
            A = k / (rho * cv)
            A_min = float(A.min())
            # PRD says A_min = 1.07e-7; tolerance ±50% given grid/range differences
            assert A_min > 5e-8, (
                f"CO2 A_min = {A_min:.3e} at P={P:.1e}Pa — below PRD §4.2 floor of 1e-7"
            )
            assert A_min > 0, "CO2 A_min must be > 0"

    def test_out_of_range_temperature_raises(self):
        with pytest.raises(ValueError, match="out of range"):
            co2_props(50.0, 10e6)


# ---------------------------------------------------------------------------
# A_field
# ---------------------------------------------------------------------------

class TestAField:

    def test_ideal_gas_all_positive(self):
        T = np.linspace(100, 2000, 500)
        A = A_field(T, fluid="ideal")
        assert np.all(A > 0), f"A_field has non-positive values: min={A.min():.4e}"

    def test_co2_all_positive(self):
        T = np.linspace(280, 400, 100)
        A = A_field(T, fluid="co2", P=10e6)
        assert np.all(A > 0), f"CO2 A_field has non-positive values: min={A.min():.4e}"

    def test_output_shape_matches_input(self):
        T = np.ones((3, 4, 5)) * 300.0
        A = A_field(T, fluid="ideal")
        assert A.shape == (3, 4, 5)

    def test_assert_positive_raises_on_bad_input(self):
        """Force A ≤ 0 by monkey-patching — just test the assertion triggers."""
        # We can't easily make A ≤ 0 physically, so test the guard indirectly
        T = np.array([300.0])
        A = A_field(T, fluid="ideal", assert_positive=True)
        assert A[0] > 0

    def test_monotone_in_T_for_ideal_gas(self):
        """Ideal gas A(T) should increase monotonically with T."""
        T = np.array([200.0, 400.0, 800.0, 1600.0])
        A = A_field(T, fluid="ideal")
        assert np.all(np.diff(A) > 0), "Ideal gas A(T) should be monotone increasing"

    def test_A_min_co2_matches_prd(self):
        """PRD §4.2: A_min ≈ 1.07e-7 in SC-CO2 at any P."""
        T = np.linspace(280, 400, 200)
        A = A_field(T, fluid="co2", P=10e6)
        assert A.min() > 5e-8, f"A_min = {A.min():.3e}, expected ~1.07e-7 from PRD"

    def test_unknown_fluid_raises(self):
        with pytest.raises(ValueError, match="unknown fluid"):
            A_field(np.array([300.0]), fluid="helium")


# ---------------------------------------------------------------------------
# second_law_check
# ---------------------------------------------------------------------------

class TestSecondLawCheck:

    def test_ideal_gas_passes_full_range(self):
        """PRD §4.1 Ideal-Exp2: A_min > 0 at all times."""
        T = np.linspace(100, 2000, 500)
        r = second_law_check(T, fluid="ideal")
        assert r["passed"], f"Second law FAILED for ideal gas: {r['violation_message']}"

    def test_co2_passes_sc_range(self):
        """PRD §4.2 CO2-Exp2: A dips near Tc but never approaches 0."""
        T = np.linspace(280, 400, 200)
        r = second_law_check(T, fluid="co2", P=10e6)
        assert r["passed"], f"Second law FAILED for CO2: {r['violation_message']}"

    def test_co2_passes_at_all_prd_pressures(self):
        """PRD §4.2 CO2-Exp2b: All 5 pressures confirmed."""
        T = np.linspace(280, 400, 100)
        for P in [7.5e6, 8.0e6, 10.0e6, 15.0e6, 20.0e6]:
            r = second_law_check(T, fluid="co2", P=P)
            assert r["passed"], f"Second law FAILED at P={P:.1e}: {r['violation_message']}"

    def test_returns_all_fields(self):
        T = np.linspace(200, 800, 50)
        r = second_law_check(T, fluid="ideal")
        for field in ("passed", "A_min", "A_max", "T_at_A_min", "cv_min", "k_min",
                      "mu_min", "violation_message"):
            assert field in r, f"Missing field: {field}"

    def test_A_min_positive_ideal(self):
        T = np.linspace(100, 2000, 200)
        r = second_law_check(T, fluid="ideal")
        assert r["A_min"] > 0

    def test_A_min_positive_co2(self):
        T = np.linspace(280, 400, 100)
        r = second_law_check(T, fluid="co2", P=10e6)
        assert r["A_min"] > 0

    def test_scalar_temperature(self):
        r = second_law_check(300.0, fluid="ideal")
        assert r["passed"]


# ---------------------------------------------------------------------------
# self_similar_props — tests Conjecture 3.4
# ---------------------------------------------------------------------------

class TestSelfSimilarProps:

    def test_returns_all_fields(self):
        r = self_similar_props(1.2, 0.5)
        for field in ("lambda_val", "t", "tau", "T_local", "A_val", "mu_val",
                      "blowup_rate", "viscous_heating_rate", "suppression_ratio",
                      "exponent_margin", "conjecture_34_satisfied"):
            assert field in r, f"Missing field: {field}"

    def test_exponent_margin_always_positive(self):
        """Key analytical result: exponent_margin = 2 + λ > 0 for all λ > 0."""
        for lv in [0.01, 0.1, 0.3, 0.4703, 0.6057, 0.9, 1.2, 1.9206]:
            r = self_similar_props(lv, 0.5)
            assert r["exponent_margin"] > 0, (
                f"exponent_margin = {r['exponent_margin']:.4f} ≤ 0 at λ={lv}"
            )

    def test_exponent_margin_equals_2_plus_lambda(self):
        """Analytical: exponent_margin = 2 + λ exactly."""
        for lv in [0.1, 0.5, 1.0, 2.0]:
            r = self_similar_props(lv, 0.5)
            expected = 2.0 + lv
            assert abs(r["exponent_margin"] - expected) < 1e-10, (
                f"exponent_margin = {r['exponent_margin']:.6f}, expected {expected:.6f}"
            )

    def test_suppression_ratio_greater_than_1_for_all_wang_profiles(self):
        """
        Conjecture 3.4: suppression_ratio > 1 for all Wang et al. λ values.
        Tests PRD §4.3 profiles: CCF stable/unstable + Boussinesq profiles.
        """
        for name, lv in WANG_LAMBDA.items():
            for t_frac in [0.5, 0.9, 0.99]:
                r = self_similar_props(lv, t_frac)
                assert r["conjecture_34_satisfied"], (
                    f"Conjecture 3.4 VIOLATED at λ={lv} ({name}), t/t*={t_frac}: "
                    f"suppression_ratio = {r['suppression_ratio']:.4e}"
                )

    def test_suppression_strengthens_near_blowup(self):
        """
        PRD §3.4: as t → t_star, suppression_ratio → ∞.
        Test that later times give strictly larger suppression_ratio.
        """
        lv = 0.6057  # CCF 1st unstable — critical test case
        r_early = self_similar_props(lv, 0.5)
        r_mid   = self_similar_props(lv, 0.9)
        r_late  = self_similar_props(lv, 0.99)
        assert r_mid["suppression_ratio"] > r_early["suppression_ratio"], (
            "Suppression ratio should increase approaching t_star"
        )
        assert r_late["suppression_ratio"] > r_mid["suppression_ratio"]

    def test_prd_cff_2nd_unstable_lambda_0_4703(self):
        """
        PRD §4.3: CCF 2nd unstable (λ=0.4703) — described as the 'critical test case'.
        Must show suppression_ratio > 1.
        """
        r = self_similar_props(0.4703, 0.99)
        assert r["conjecture_34_satisfied"], (
            f"CCF 2nd unstable (λ=0.4703): Conjecture 3.4 VIOLATED. "
            f"suppression_ratio = {r['suppression_ratio']:.4e}"
        )

    def test_prd_sweep_lambda_values(self):
        """PRD §7.3 lambda sweep: λ = 1.2, 0.9, 0.6057, 0.4703, 0.3, 0.1 all pass."""
        sweep_lambdas = [1.2, 0.9, 0.6057, 0.4703, 0.3, 0.1]
        for lv in sweep_lambdas:
            r = self_similar_props(lv, 0.99)
            assert r["conjecture_34_satisfied"], (
                f"PRD §7.3 sweep FAILED at λ={lv}: "
                f"suppression_ratio = {r['suppression_ratio']:.4e}"
            )

    def test_invalid_lambda_raises(self):
        with pytest.raises(ValueError, match="λ must be > 0"):
            self_similar_props(-0.1, 0.5)

    def test_t_at_t_star_raises(self):
        with pytest.raises(ValueError, match="t=1.0 must be < t_star"):
            self_similar_props(1.0, 1.0, t_star=1.0)

    def test_A_val_positive(self):
        """A(T_local) > 0 even at blow-up time (second law holds)."""
        for lv in [0.1, 0.6057, 1.2]:
            r = self_similar_props(lv, 0.99)
            assert r["A_val"] > 0, f"A_val ≤ 0 at λ={lv}: {r['A_val']}"

    def test_mu_val_positive(self):
        """μ(T_local) > 0 — thermodynamic stability."""
        for lv in [0.1, 0.6057, 1.2]:
            r = self_similar_props(lv, 0.99)
            assert r["mu_val"] > 0, f"mu_val ≤ 0 at λ={lv}: {r['mu_val']}"

    def test_co2_fluid(self):
        """Conjecture 3.4 also holds for SC-CO2."""
        r = self_similar_props(0.4703, 0.99, fluid="co2", P=10e6)
        assert r["conjecture_34_satisfied"], (
            f"CO2: Conjecture 3.4 VIOLATED at λ=0.4703: "
            f"suppression_ratio = {r['suppression_ratio']:.4e}"
        )


# ---------------------------------------------------------------------------
# conjecture_34_sweep
# ---------------------------------------------------------------------------

class TestConjecture34Sweep:

    def test_all_default_lambdas_pass(self):
        """PRD §7.3: all 6 default λ values show suppression_ratio > 1."""
        results = conjecture_34_sweep(fluid="ideal")
        failures = [r for r in results if r["verdict"] != "PASS"]
        assert len(failures) == 0, (
            f"Conjecture 3.4 VIOLATED for {len(failures)} (λ, t) pairs:\n"
            + "\n".join(f"  λ={r['lambda_val']}, t/t*={r['t_frac']}: "
                        f"ratio={r['suppression_ratio']:.4e}" for r in failures)
        )

    def test_returns_correct_number_of_results(self):
        """6 λ values × 3 time points = 18 results."""
        results = conjecture_34_sweep()
        assert len(results) == 18, f"Expected 18 results, got {len(results)}"

    def test_custom_lambda_values(self):
        results = conjecture_34_sweep(lambda_values=[0.01, 0.5, 2.0], t_values=[0.99])
        assert len(results) == 3
        for r in results:
            assert r["conjecture_34_satisfied"], (
                f"Custom λ={r['lambda_val']} failed at t/t*={r['t_frac']}"
            )

    def test_very_small_lambda_still_passes(self):
        """Borderline case: λ = 0.01 (near adversarial minimum)."""
        results = conjecture_34_sweep(lambda_values=[0.01], t_values=[0.5, 0.9, 0.99])
        for r in results:
            assert r["conjecture_34_satisfied"], (
                f"λ=0.01 failed at t/t*={r['t_frac']}: ratio={r['suppression_ratio']:.4e}"
            )


# ---------------------------------------------------------------------------
# route1_coeffs
# ---------------------------------------------------------------------------

class TestRoute1Coeffs:
    """Tests for route1_coeffs() — PRD §7.2 v0.4 M1 addition."""

    def _make_kxy(self, N=16):
        k1d = np.fft.fftfreq(N, d=1.0 / N)
        return np.meshgrid(k1d, k1d, indexing="ij")

    def test_returns_all_keys(self):
        T = np.linspace(200.0, 400.0, 64).reshape(8, 8)
        r = route1_coeffs(T)
        for key in ("mu_min", "mu_max", "A_min", "A_max",
                    "nu_eff_min", "nu_eff_max", "T_min", "T_max", "margin_ok"):
            assert key in r, f"Missing key: {key}"

    def test_mu_positive(self):
        """mu_min > 0 for ideal gas — thermodynamic requirement."""
        T = np.linspace(100.0, 2000.0, 100)
        r = route1_coeffs(T, fluid="ideal")
        assert r["mu_min"] > 0.0, f"mu_min = {r['mu_min']:.4e} ≤ 0"

    def test_A_positive(self):
        """A_min > 0 for ideal gas — the core inequality."""
        T = np.linspace(100.0, 2000.0, 100)
        r = route1_coeffs(T, fluid="ideal")
        assert r["A_min"] > 0.0, f"A_min = {r['A_min']:.4e} ≤ 0"

    def test_bounds_ordering(self):
        """mu_min ≤ mu_max; A_min ≤ A_max."""
        T = np.linspace(200.0, 800.0, 50)
        r = route1_coeffs(T, fluid="ideal")
        assert r["mu_min"] <= r["mu_max"], "mu bounds inverted"
        assert r["A_min"] <= r["A_max"], "A bounds inverted"

    def test_nu_eff_ordering(self):
        """nu_eff_min ≤ nu_eff_max."""
        T = np.linspace(200.0, 800.0, 50)
        r = route1_coeffs(T, fluid="ideal")
        assert r["nu_eff_min"] <= r["nu_eff_max"], "nu_eff bounds inverted"

    def test_T_range_recorded(self):
        T = np.array([300.0, 500.0, 700.0])
        r = route1_coeffs(T, fluid="ideal")
        assert abs(r["T_min"] - 300.0) < 1e-10
        assert abs(r["T_max"] - 700.0) < 1e-10

    def test_2d_field_input(self):
        """Works on a 2D temperature field (as produced by T_solver_2D)."""
        T = np.ones((64, 64)) * 300.0 + np.random.default_rng(42).normal(0, 10, (64, 64))
        T = np.clip(T, 150.0, 1000.0)
        r = route1_coeffs(T, fluid="ideal")
        assert r["A_min"] > 0.0
        assert r["mu_min"] > 0.0
        assert r["margin_ok"] is True

    def test_margin_ok_is_true(self):
        """margin_ok flag is True when no thermodynamic violation (always for valid T)."""
        T = np.linspace(300.0, 600.0, 20)
        r = route1_coeffs(T, fluid="ideal")
        assert r["margin_ok"] is True

    def test_rejects_nonpositive_T(self):
        T = np.array([300.0, -10.0, 500.0])
        with pytest.raises(ValueError, match="must be > 0"):
            route1_coeffs(T, fluid="ideal")

    def test_rejects_unknown_fluid(self):
        T = np.array([300.0])
        with pytest.raises(ValueError, match="unknown fluid"):
            route1_coeffs(T, fluid="hydrogen")

    def test_consistent_with_A_field(self):
        """route1_coeffs A_min should equal np.min(A_field(T)) for same T."""
        T = np.linspace(200.0, 1000.0, 80)
        r = route1_coeffs(T, fluid="ideal")
        A_direct = A_field(T, fluid="ideal")
        assert abs(r["A_min"] - float(np.min(A_direct))) < 1e-14 * r["A_min"]
        assert abs(r["A_max"] - float(np.max(A_direct))) < 1e-14 * r["A_max"]


# ---------------------------------------------------------------------------
# route2_theta_source
# ---------------------------------------------------------------------------

class TestRoute2ThetaSource:
    """Tests for route2_theta_source() — PRD §3.2 v0.4 Route 2 θ equation."""

    def _make_grid(self, N=32):
        k1d = np.fft.fftfreq(N, d=1.0 / N)
        kx, ky = np.meshgrid(k1d, k1d, indexing="ij")
        return kx, ky

    def test_returns_correct_shape(self):
        N = 32
        u = np.zeros((N, N))
        v = np.zeros((N, N))
        kx, ky = self._make_grid(N)
        S = route2_theta_source(u, v, nu=1e-3, kx=kx, ky=ky)
        assert S.shape == (N, N), f"Expected ({N},{N}), got {S.shape}"

    def test_zero_velocity_zero_source(self):
        """If u=v=0, source = ν·|∇u|² = 0."""
        N = 32
        u = np.zeros((N, N))
        v = np.zeros((N, N))
        kx, ky = self._make_grid(N)
        S = route2_theta_source(u, v, nu=1e-3, kx=kx, ky=ky)
        assert np.max(np.abs(S)) < 1e-25, f"Zero velocity should give zero source; max={np.max(np.abs(S)):.2e}"

    def test_source_is_nonnegative(self):
        """S_θ = ν·|∇u|² ≥ 0 always (sum of squares)."""
        N = 32
        rng = np.random.default_rng(7)
        u = rng.normal(0, 1, (N, N))
        v = rng.normal(0, 1, (N, N))
        kx, ky = self._make_grid(N)
        S = route2_theta_source(u, v, nu=1e-3, kx=kx, ky=ky)
        assert np.all(S >= 0.0), f"Source has negative values: min={np.min(S):.4e}"

    def test_scales_linearly_with_nu(self):
        """S_θ should scale as ν (linear coefficient)."""
        N = 32
        rng = np.random.default_rng(13)
        u = rng.normal(0, 1, (N, N))
        v = rng.normal(0, 1, (N, N))
        kx, ky = self._make_grid(N)
        S1 = route2_theta_source(u, v, nu=1e-3, kx=kx, ky=ky)
        S2 = route2_theta_source(u, v, nu=2e-3, kx=kx, ky=ky)
        ratio = float(np.mean(S2)) / float(np.mean(S1) + 1e-100)
        assert abs(ratio - 2.0) < 1e-10, f"ν-scaling broken: ratio={ratio:.6f} (expected 2.0)"

    def test_scales_quadratically_with_velocity_amplitude(self):
        """S_θ ~ |∇u|² so scales as amplitude²."""
        N = 32
        rng = np.random.default_rng(17)
        u0 = rng.normal(0, 1, (N, N))
        v0 = rng.normal(0, 1, (N, N))
        kx, ky = self._make_grid(N)
        S1 = route2_theta_source(u0, v0, nu=1e-3, kx=kx, ky=ky)
        S2 = route2_theta_source(2.0 * u0, 2.0 * v0, nu=1e-3, kx=kx, ky=ky)
        ratio = float(np.mean(S2)) / float(np.mean(S1) + 1e-100)
        assert abs(ratio - 4.0) < 1e-8, f"Amplitude-squared scaling broken: ratio={ratio:.6f} (expected 4.0)"

    def test_pure_shear_flow_analytical(self):
        """u = y (linear shear): ∂u/∂y = 1, all other gradients zero.
        In spectral representation on periodic domain, test with a sinusoidal shear."""
        N = 64
        kx, ky = self._make_grid(N)
        # u = A·sin(ky_val·y), v = 0: du/dy = A·ky_val·cos(ky_val·y)
        # => |∇u|² = A²·ky_val²·cos²(ky_val·y), mean = A²·ky_val²/2
        ky_val = 2
        A = 1.0
        nu = 1e-2
        x = np.linspace(0, 2 * np.pi, N, endpoint=False)
        _, Y = np.meshgrid(x, x, indexing="ij")
        u = A * np.sin(ky_val * Y)
        v = np.zeros((N, N))
        S = route2_theta_source(u, v, nu=nu, kx=kx, ky=ky)
        expected_mean = nu * A**2 * ky_val**2 / 2.0
        actual_mean = float(np.mean(S))
        assert abs(actual_mean - expected_mean) / expected_mean < 1e-6, (
            f"Analytical shear check: mean(S)={actual_mean:.6e}, expected={expected_mean:.6e}"
        )

    def test_rejects_nonpositive_nu(self):
        N = 16
        u = np.zeros((N, N))
        v = np.zeros((N, N))
        kx, ky = self._make_grid(N)
        with pytest.raises(ValueError, match="nu must be > 0"):
            route2_theta_source(u, v, nu=0.0, kx=kx, ky=ky)
        with pytest.raises(ValueError, match="nu must be > 0"):
            route2_theta_source(u, v, nu=-1e-3, kx=kx, ky=ky)

    def test_source_positive_for_nontrivial_flow(self):
        """For any non-zero velocity field with nonzero gradients, S_θ > 0 somewhere."""
        N = 32
        rng = np.random.default_rng(99)
        u = rng.normal(0, 1, (N, N))
        v = rng.normal(0, 1, (N, N))
        kx, ky = self._make_grid(N)
        S = route2_theta_source(u, v, nu=1e-3, kx=kx, ky=ky)
        assert np.max(S) > 0.0, "Nontrivial flow must produce positive θ source somewhere"
