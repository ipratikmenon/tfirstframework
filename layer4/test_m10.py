"""
test_m10.py — Test suite for M10 diagnostic tools.

Coverage:
  continuation_criterion: A/μ margin checks, full verdict, blowup result analysis
  gronwall_fitter: E1 extraction, fit on synthetic decay, history fit
  viscosity_scaling: blowup_rate, power law fit, suppression ratio, verdict
  mu_limit_tracker: eps computation, scaling fit, verdict, 2D/3D comparison
"""

from __future__ import annotations

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

# ── imports ───────────────────────────────────────────────────────────────────

from layer4.continuation_criterion import (
    check_A_min,
    check_mu_min,
    continuation_criterion_verdict,
    A_ref_at,
    mu_ref_at,
    A_of_T,
    mu_of_T,
    analyse_blowup_result,
)

from layer4.gronwall_fitter import (
    fit_gronwall,
    gronwall_verdict,
    fit_energy_timeseries,
    extract_E1_spectral,
)

from layer4.viscosity_scaling import (
    blowup_rate,
    fit_power_law_decay,
    suppression_ratio,
    viscosity_scaling_verdict,
    mu_timeseries_from_T,
)

from layer4.mu_limit_tracker import (
    compute_eps,
    fit_scaling_exponent,
    mu_limit_verdict,
    run_mu_limit_analysis,
    compare_2d_3d,
    load_from_sweep_results,
)


# ─────────────────────────────────────────────────────────────────────────────
# TestContinuationCriterion
# ─────────────────────────────────────────────────────────────────────────────

class TestContinuationCriterion:
    def setup_method(self):
        self.A_ref  = A_ref_at(300.0)
        self.mu_ref = mu_ref_at(300.0)

    # -- check_A_min -----------------------------------------------------------

    def test_not_triggered_when_A_min_well_above_threshold(self):
        A_series = [self.A_ref * 0.5, self.A_ref * 0.8, self.A_ref]
        r = check_A_min(A_series, self.A_ref, threshold_frac=0.01)
        assert not r["triggered"], f"Should not trigger: margin={r['margin']:.3f}"

    def test_triggered_when_A_min_near_zero(self):
        A_series = [self.A_ref, self.A_ref * 0.001]
        r = check_A_min(A_series, self.A_ref, threshold_frac=0.01)
        assert r["triggered"]

    def test_margin_is_ratio(self):
        A_series = [self.A_ref * 0.5]
        r = check_A_min(A_series, self.A_ref)
        assert abs(r["margin"] - 0.5) < 1e-6

    def test_empty_series_not_triggered(self):
        r = check_A_min([], self.A_ref)
        assert not r["triggered"]
        assert r["n_points"] == 0

    # -- check_mu_min ----------------------------------------------------------

    def test_mu_not_triggered_physical(self):
        mu_series = [self.mu_ref * 0.9, self.mu_ref * 0.95, self.mu_ref]
        r = check_mu_min(mu_series, self.mu_ref)
        assert not r["triggered"]

    def test_mu_triggered_near_zero(self):
        mu_series = [self.mu_ref * 0.001]
        r = check_mu_min(mu_series, self.mu_ref, threshold_frac=0.01)
        assert r["triggered"]

    # -- A_of_T / mu_of_T ------------------------------------------------------

    def test_A_of_T_positive(self):
        for T in [100.0, 200.0, 300.0, 500.0, 1000.0]:
            assert A_of_T(T) > 0.0

    def test_A_increases_with_T(self):
        assert A_of_T(400.0) > A_of_T(200.0)

    def test_mu_positive(self):
        for T in [150.0, 300.0, 600.0]:
            assert mu_of_T(T) > 0.0

    # -- full verdict ----------------------------------------------------------

    def test_verdict_pass_for_normal_history(self):
        history = [{"t": i * 0.1, "A_min": self.A_ref * (0.9 + 0.01 * i)}
                   for i in range(20)]
        r = continuation_criterion_verdict(history)
        assert r["verdict"] == "PASS"
        assert not r["criterion_triggered"]

    def test_verdict_fail_when_A_min_collapses(self):
        history = [{"t": 0.1, "A_min": self.A_ref},
                   {"t": 0.5, "A_min": self.A_ref * 0.001}]
        r = continuation_criterion_verdict(history, threshold_frac=0.01)
        assert r["verdict"] == "FAIL"
        assert r["criterion_triggered"]

    def test_empty_history_pass(self):
        r = continuation_criterion_verdict([])
        assert r["verdict"] == "PASS"

    def test_required_keys(self):
        history = [{"t": 0.1, "A_min": self.A_ref}]
        r = continuation_criterion_verdict(history)
        for k in ["verdict", "A_check", "mu_check", "criterion_triggered",
                  "A_ref", "mu_ref"]:
            assert k in r, f"Missing key: {k}"

    def test_analyse_blowup_result(self):
        result = {
            "history": [{"t": 0.1 * i, "A_min": self.A_ref * 0.9}
                        for i in range(10)]
        }
        r = analyse_blowup_result(result)
        assert r["verdict"] == "PASS"


# ─────────────────────────────────────────────────────────────────────────────
# TestGronwallFitter
# ─────────────────────────────────────────────────────────────────────────────

class TestGronwallFitter:
    def _synthetic_decay(self, c=1.5, E0=0.1, n=30, t_end=2.0, noise=0.0):
        t = np.linspace(0.0, t_end, n)
        E = E0 * np.exp(-c * t)
        if noise > 0.0:
            rng = np.random.default_rng(42)
            E = E * (1.0 + noise * rng.standard_normal(n))
            E = np.abs(E)
        return t, E

    def test_fit_perfect_decay(self):
        c_true = 1.5
        t, E = self._synthetic_decay(c=c_true, n=50)
        result = fit_gronwall(E, t)
        assert result["fit_ok"]
        assert abs(result["c"] - c_true) < 0.05, \
            f"c={result['c']:.4f} != c_true={c_true}"

    def test_fit_R2_near_one_for_clean_data(self):
        t, E = self._synthetic_decay(c=2.0, n=50)
        result = fit_gronwall(E, t)
        assert result["R2"] > 0.99

    def test_verdict_pass_for_positive_c(self):
        t, E = self._synthetic_decay(c=0.5, n=20)
        result = fit_gronwall(E, t)
        assert result["verdict"] == "PASS"

    def test_verdict_fail_for_negative_c(self):
        # Growing E: c < 0
        t = np.linspace(0.0, 1.0, 20)
        E = 0.1 * np.exp(2.0 * t)  # growing
        result = fit_gronwall(E, t)
        # c should be negative
        assert result["c"] < 0.0 or result["verdict"] == "FAIL"

    def test_insufficient_data(self):
        result = fit_gronwall([0.1, 0.09], [0.0, 0.1])
        assert not result["fit_ok"]
        assert result["verdict"] == "INSUFFICIENT_DATA"

    def test_gronwall_verdict_pass(self):
        assert gronwall_verdict(c=1.0, R2=0.95) == "PASS"

    def test_gronwall_verdict_fail_negative_c(self):
        assert gronwall_verdict(c=-0.1, R2=0.99) == "FAIL"

    def test_gronwall_verdict_fail_low_R2(self):
        assert gronwall_verdict(c=1.0, R2=0.3) == "FAIL"

    def test_fit_energy_timeseries(self):
        c_true = 1.0
        t, E = self._synthetic_decay(c=c_true, n=40)
        history = [{"t": float(ti), "E": float(Ei)} for ti, Ei in zip(t, E)]
        result = fit_energy_timeseries(history)
        assert result["fit_ok"]
        assert abs(result["c"] - c_true) < 0.1

    def test_E1_spectral_positive(self):
        N = 16
        x = np.linspace(0, 2 * np.pi, N, endpoint=False).astype(np.float32)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        u = np.sin(X).astype(np.float32)
        v = np.cos(Y).astype(np.float32)
        w = np.zeros_like(u)
        freqs = np.fft.fftfreq(N, d=1.0 / N).astype(np.float32)
        kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing="ij")
        E1 = extract_E1_spectral(u, v, w, kx, ky, kz)
        assert E1 > 0.0

    def test_E1_zero_for_zero_velocity(self):
        N = 8
        u = np.zeros((N, N, N), dtype=np.float32)
        v = np.zeros_like(u)
        w = np.zeros_like(u)
        freqs = np.fft.fftfreq(N, d=1.0 / N).astype(np.float32)
        kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing="ij")
        E1 = extract_E1_spectral(u, v, w, kx, ky, kz)
        assert abs(E1) < 1e-10


# ─────────────────────────────────────────────────────────────────────────────
# TestViscosityScaling
# ─────────────────────────────────────────────────────────────────────────────

class TestViscosityScaling:
    def test_blowup_rate_shape(self):
        t = np.linspace(0.0, 0.9, 10)
        rate = blowup_rate(t, T_star=1.0, lambda_val=0.6)
        assert rate.shape == (10,)

    def test_blowup_rate_diverges_near_T_star(self):
        t = np.array([0.0, 0.5, 0.9, 0.99])
        rate = blowup_rate(t, T_star=1.0, lambda_val=1.0)
        assert rate[-1] > rate[-2] > rate[0]

    def test_blowup_rate_nan_beyond_T_star(self):
        t = np.array([0.5, 1.0, 1.5])
        rate = blowup_rate(t, T_star=1.0, lambda_val=0.5)
        assert np.isnan(rate[2])
        assert np.isnan(rate[1])
        assert not np.isnan(rate[0])

    def test_fit_power_law_decay_synthetic(self):
        tau = np.array([0.5, 0.3, 0.2, 0.1, 0.05])
        alpha_true = -0.5
        f = 2.0 * tau ** alpha_true
        t = 1.0 - tau  # t = T_star - tau
        result = fit_power_law_decay(t, f, T_star=1.0)
        assert result["fit_ok"]
        assert abs(result["alpha"] - alpha_true) < 0.05

    def test_fit_power_law_insufficient_data(self):
        result = fit_power_law_decay([0.1, 0.2], [1.0, 2.0], T_star=1.0)
        assert not result["fit_ok"]

    def test_suppression_ratio_shape(self):
        t = np.linspace(0.0, 0.8, 10)
        mu = np.full(10, 1.716e-5)
        ratio = suppression_ratio(t, mu, T_star=1.0, lambda_val=0.5)
        assert ratio.shape == (10,)
        assert np.all(np.isfinite(ratio))

    def test_viscosity_scaling_verdict_constant_mu(self):
        t = np.linspace(0.0, 0.9, 20)
        mu = np.full(20, 1.716e-5)
        result = viscosity_scaling_verdict(t, mu, T_star=1.5, lambda_val=0.5)
        assert "verdict" in result
        assert "suppression_ok" in result

    def test_mu_timeseries_from_T_positive(self):
        T_series = [150.0, 200.0, 250.0, 300.0]
        mu = mu_timeseries_from_T(T_series)
        assert len(mu) == 4
        assert np.all(mu > 0)

    def test_mu_increases_with_T(self):
        mu_low  = mu_timeseries_from_T([200.0])[0]
        mu_high = mu_timeseries_from_T([400.0])[0]
        assert mu_high > mu_low


# ─────────────────────────────────────────────────────────────────────────────
# TestMuLimitTracker
# ─────────────────────────────────────────────────────────────────────────────

class TestMuLimitTracker:
    _A_REF = 3.0801e-5  # M4/M8 measured A_ref

    def _synthetic_linear(self, n=10):
        """ε = 0.01 * δ (exactly linear: α = 1)."""
        deltas = np.logspace(0, -3, n)
        A_min  = self._A_REF + 0.01 * deltas
        return deltas, A_min

    def test_compute_eps_shape(self):
        deltas, A_min = self._synthetic_linear()
        eps = compute_eps(deltas, A_min, self._A_REF)
        assert eps.shape == (len(deltas),)

    def test_compute_eps_positive(self):
        deltas, A_min = self._synthetic_linear()
        eps = compute_eps(deltas, A_min, self._A_REF)
        assert np.all(eps >= 0.0)

    def test_eps_decreases_with_delta(self):
        """Linear ε ~ δ: smaller δ → smaller ε."""
        deltas = np.array([1.0, 0.5, 0.1, 0.01])
        A_min  = self._A_REF + 0.02 * deltas
        eps    = compute_eps(deltas, A_min, self._A_REF)
        assert np.all(np.diff(eps) < 0), "ε should decrease as δ decreases"

    def test_fit_scaling_linear(self):
        deltas, A_min = self._synthetic_linear(n=8)
        eps = compute_eps(deltas, A_min, self._A_REF)
        fit = fit_scaling_exponent(deltas, eps)
        assert fit["fit_ok"]
        assert abs(fit["alpha"] - 1.0) < 0.1, \
            f"Expected α≈1.0 (linear), got {fit['alpha']:.3f}"

    def test_fit_scaling_quadratic(self):
        deltas = np.logspace(0, -2, 8)
        eps    = 0.01 * deltas ** 2
        fit = fit_scaling_exponent(deltas, eps)
        assert fit["fit_ok"]
        assert abs(fit["alpha"] - 2.0) < 0.1

    def test_fit_insufficient_data(self):
        fit = fit_scaling_exponent([0.1], [1e-5])
        assert not fit["fit_ok"]

    def test_mu_limit_verdict_linear_pass(self):
        v = mu_limit_verdict(alpha=0.98, eps_values=[0.01, 0.005, 0.001])
        assert v == "PASS"

    def test_mu_limit_verdict_super_linear_fail(self):
        v = mu_limit_verdict(alpha=2.5, eps_values=[0.01, 0.001, 0.0001],
                              alpha_tolerance=0.3)
        assert v == "FAIL"

    def test_mu_limit_verdict_all_zero_eps_pass(self):
        v = mu_limit_verdict(alpha=float("nan"), eps_values=[0.0, 0.0, 0.0])
        assert v == "PASS"

    def test_run_mu_limit_analysis_m8_data(self):
        """Use actual M8 measured values from SESSION-LOG S14."""
        deltas = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
        A_mins = [
            3.0636e-05, 3.0719e-05, 3.0768e-05, 3.0785e-05,
            3.0793e-05, 3.0798e-05, 3.0800e-05, 3.0800e-05,
            3.0801e-05, 3.0801e-05,
        ]
        result = run_mu_limit_analysis(deltas, A_mins, self._A_REF, label="3D M8")
        assert result["verdict"] == "PASS"
        assert abs(result["alpha"] - 1.0) < 0.15, \
            f"Expected α≈1.0, got {result['alpha']:.3f}"
        assert result["label"] == "3D M8"

    def test_run_mu_limit_analysis_required_keys(self):
        deltas, A_min = self._synthetic_linear(n=6)
        result = run_mu_limit_analysis(deltas, A_min, self._A_REF)
        for k in ["delta_values", "A_min_values", "A_ref", "eps_values",
                  "scaling_fit", "verdict", "A_min_plateau", "alpha"]:
            assert k in result, f"Missing key: {k}"

    def test_compare_2d_3d_consistent(self):
        d = [1.0, 0.5, 0.1, 0.01]
        A_m = [self._A_REF + 0.01 * dd for dd in d]
        r2d = run_mu_limit_analysis(d, A_m, self._A_REF, label="2D")
        r3d = run_mu_limit_analysis(d, A_m, self._A_REF, label="3D")
        cmp = compare_2d_3d(r2d, r3d)
        assert cmp["consistent"]
        assert cmp["verdict"] == "CONSISTENT"

    def test_load_from_sweep_results(self):
        sweep = {
            d: {"A_min_global": self._A_REF + 0.01 * d}
            for d in [1.0, 0.5, 0.1, 0.01]
        }
        result = load_from_sweep_results(sweep, self._A_REF, label="test")
        assert result["verdict"] == "PASS"
        assert len(result["delta_values"]) == 4
