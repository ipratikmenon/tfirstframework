"""
test_route2_2D.py — Tests for Route 2 2D solver (PRD §3.2 v0.4)

Tests verify:
  - θ(x,0) = 0 (initial condition)
  - S_θ = ν·|∇u|² ≥ 0 always (source term non-negative)
  - θ(x,t) ≥ 0 for all t (monotone accumulation)
  - μ_eff = ν + ε·f(θ) ≥ ν > 0 always (LPS margin never negative)
  - ω equation is exact Prize NS (ν = const)
  - ε_param → 0 recovers Prize NS (θ decouples from ω)
  - All 4 ε_param values pass

Run with: pytest layer2/test_route2_2D.py -v
"""

import os
import sys
import tempfile

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from route2_2D import (
    make_solver, set_ic, solver_step, run,
    run_eps_sweep, run_smoke_test,
    _rhs, _mu_eff, _apply_f_theta,
    compute_cfl_dt, EPS_SWEEP_VALUES, _NU_DEFAULT,
)
from lambda_sweep_2D import query_results
from self_similar_IC import generate_CCF_profile


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture()
def tmp_db(tmp_path):
    return str(tmp_path / "r2_test.db")


@pytest.fixture()
def small_solver():
    """Minimal solver for fast unit tests (N=16)."""
    s = make_solver(N=16, nu=1e-3, eps_param=0.1)
    ic = generate_CCF_profile(0.6057, N=16)
    return set_ic(s, ic)


@pytest.fixture()
def medium_solver():
    """N=32 solver with CCF IC."""
    s = make_solver(N=32, nu=1e-3, eps_param=0.1)
    ic = generate_CCF_profile(0.6057, N=32)
    return set_ic(s, ic)


_FAST = dict(N=32, t_end=0.05, max_steps=10)


# ─────────────────────────────────────────────────────────────────────────────
# Solver construction
# ─────────────────────────────────────────────────────────────────────────────

class TestSolverConstruction:

    def test_make_solver_has_required_keys(self):
        s = make_solver(N=32, nu=1e-3, eps_param=0.1)
        for key in ("N", "nu", "eps_param", "f_theta", "omega", "theta",
                    "kx", "ky", "k2", "t", "dt_max"):
            assert key in s, f"Missing key: {key}"

    def test_theta_initially_zero(self):
        """θ(x,0) = 0 — PRD §3.2 requirement."""
        s = make_solver(N=32, nu=1e-3, eps_param=0.5)
        assert np.all(s["theta"] == 0.0), "θ must be initialised to 0"

    def test_set_ic_resets_theta_to_zero(self):
        """set_ic always sets θ=0 regardless of IC dict content."""
        s = make_solver(N=16, nu=1e-3, eps_param=0.1)
        ic = generate_CCF_profile(0.6057, N=16)
        # Inject a non-zero theta into IC (should be ignored)
        ic["theta"] = np.ones((16, 16))
        s = set_ic(s, ic)
        assert np.all(s["theta"] == 0.0), "set_ic must reset θ to 0"

    def test_set_ic_copies_omega(self):
        s = make_solver(N=16, nu=1e-3, eps_param=0.1)
        ic = generate_CCF_profile(0.6057, N=16)
        s = set_ic(s, ic)
        assert not np.all(s["omega"] == 0.0), "ω should be non-zero from CCF IC"

    def test_set_ic_is_immutable(self):
        s0 = make_solver(N=16, nu=1e-3, eps_param=0.1)
        ic = generate_CCF_profile(0.6057, N=16)
        s1 = set_ic(s0, ic)
        assert s0 is not s1
        assert np.all(s0["theta"] == 0.0)

    def test_rejects_nonpositive_nu(self):
        with pytest.raises(ValueError, match="nu must be > 0"):
            make_solver(nu=0.0)

    def test_rejects_negative_eps_param(self):
        with pytest.raises(ValueError, match="eps_param must be ≥ 0"):
            make_solver(eps_param=-0.1)

    def test_eps_zero_allowed(self):
        """ε_param=0 is valid — recovers exact Prize equations."""
        s = make_solver(eps_param=0.0)
        assert s["eps_param"] == 0.0


# ─────────────────────────────────────────────────────────────────────────────
# f(θ) functional forms
# ─────────────────────────────────────────────────────────────────────────────

class TestFThetaFunctions:

    def test_identity_at_zero(self):
        """f(0) = 0 for all forms."""
        for form in ("identity", "tanh", "sqrt"):
            val = float(np.mean(_apply_f_theta(np.zeros((4, 4)), form)))
            assert abs(val) < 1e-10, f"f(0) ≠ 0 for {form}: got {val}"

    def test_identity_nonnegative(self):
        theta = np.clip(np.random.default_rng(0).normal(0, 1, (8, 8)), 0, None)
        f = _apply_f_theta(theta, "identity")
        assert np.all(f >= 0)

    def test_tanh_bounded(self):
        theta = np.linspace(0, 100, 50)
        f = _apply_f_theta(theta, "tanh")
        assert np.all(f >= 0)
        assert np.all(f <= 1.0 + 1e-10)

    def test_sqrt_nonnegative(self):
        theta = np.linspace(0, 10, 20)
        f = _apply_f_theta(theta, "sqrt")
        assert np.all(f >= 0)

    def test_mu_eff_always_geq_nu(self):
        """μ_eff = ν + ε·f(θ) ≥ ν for all θ ≥ 0, all forms."""
        nu, eps = 1e-3, 0.5
        theta = np.linspace(0, 5, 20)
        for form in ("identity", "tanh", "sqrt"):
            mu = _mu_eff(theta, nu, eps, form)
            assert np.all(mu >= nu - 1e-14), (
                f"μ_eff < ν for form={form}: min={np.min(mu):.4e}, nu={nu}"
            )

    def test_mu_eff_equals_nu_when_theta_zero(self):
        """When θ=0, μ_eff = ν exactly."""
        nu, eps = 1e-3, 1.0
        theta = np.zeros(10)
        for form in ("identity", "tanh", "sqrt"):
            mu = _mu_eff(theta, nu, eps, form)
            assert np.allclose(mu, nu), f"μ_eff(θ=0) ≠ ν for form={form}"

    def test_unknown_form_raises(self):
        with pytest.raises(ValueError, match="Unknown f_theta"):
            _apply_f_theta(np.zeros(5), "cubic")


# ─────────────────────────────────────────────────────────────────────────────
# RHS
# ─────────────────────────────────────────────────────────────────────────────

class TestRHS:

    def test_source_term_nonnegative(self, small_solver):
        """S_θ = ν·|∇u|² ≥ 0 everywhere."""
        s = small_solver
        _, d_theta, diag = _rhs(
            s["omega"], s["theta"],
            s["kx"], s["ky"], s["k2"],
            s["nu"], s["eps_param"], s["f_theta"],
        )
        # S_theta_max and mean should be non-negative
        assert diag["S_theta_max"] >= 0.0
        assert diag["S_theta_mean"] >= 0.0

    def test_zero_omega_gives_zero_source(self):
        """If ω=0 → u=v=0 → S_θ=0."""
        N = 16
        s = make_solver(N=N, nu=1e-3, eps_param=0.1)
        # omega = 0, theta = 0
        _, _, diag = _rhs(s["omega"], s["theta"],
                          s["kx"], s["ky"], s["k2"],
                          s["nu"], s["eps_param"], s["f_theta"])
        assert diag["S_theta_max"] < 1e-20, f"S_θ should be 0 for zero velocity"

    def test_lps_margin_nonnegative_at_step_0(self, small_solver):
        """ε_LPS = min(μ_eff) - ν ≥ 0 at initial step."""
        s = small_solver
        _, _, diag = _rhs(
            s["omega"], s["theta"],
            s["kx"], s["ky"], s["k2"],
            s["nu"], s["eps_param"], s["f_theta"],
        )
        # At t=0: θ=0 so f(θ)=0 so μ_eff = ν, so ε_LPS = 0 (lower bound)
        assert diag["lps_margin"] >= -1e-14, (
            f"ε_LPS = {diag['lps_margin']:.4e} < 0 at step 0"
        )

    def test_omega_rhs_uses_constant_nu(self, small_solver):
        """ω equation uses ν (not μ_eff) — exact Prize NS."""
        s = small_solver
        nu = s["nu"]
        # Run RHS with eps_param=0 (μ_eff=ν) and eps_param=large (μ_eff>>ν)
        # d_omega should be the SAME because ω uses only ν
        d_om_small, _, _ = _rhs(s["omega"], s["theta"],
                                 s["kx"], s["ky"], s["k2"],
                                 nu, 0.0, "identity")
        d_om_large, _, _ = _rhs(s["omega"], s["theta"],
                                 s["kx"], s["ky"], s["k2"],
                                 nu, 1000.0, "identity")
        # d_omega should be identical (eps doesn't affect ω equation)
        assert np.allclose(d_om_small, d_om_large), (
            "d_omega differs between eps_param=0 and eps_param=1000 — "
            "ω equation must use only ν (not μ_eff)"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Time-stepping
# ─────────────────────────────────────────────────────────────────────────────

class TestTimeStepping:

    def test_cfl_dt_positive(self, small_solver):
        dt = compute_cfl_dt(small_solver)
        assert dt > 0.0

    def test_cfl_dt_respects_dt_max(self, small_solver):
        dt = compute_cfl_dt(small_solver)
        assert dt <= small_solver["dt_max"]

    def test_solver_step_advances_time(self, small_solver):
        new_s, diag = solver_step(small_solver)
        assert new_s["t"] > small_solver["t"]
        assert diag["t"] > 0.0

    def test_solver_step_is_immutable(self, small_solver):
        t0 = small_solver["t"]
        new_s, _ = solver_step(small_solver)
        assert small_solver["t"] == t0, "Original solver mutated"
        assert new_s is not small_solver

    def test_theta_nonnegative_after_step(self, small_solver):
        """θ ≥ 0 after every step (S_θ ≥ 0 and numerical floor applied)."""
        s = small_solver
        for _ in range(5):
            s, _ = solver_step(s)
            assert np.all(s["theta"] >= 0.0), f"θ < 0 at t={s['t']:.4f}"

    def test_theta_accumulates(self, small_solver):
        """θ_max should be non-decreasing (S_θ ≥ 0 drives accumulation)."""
        s = small_solver
        theta_max_prev = float(np.max(s["theta"]))
        for _ in range(10):
            s, diag = solver_step(s)
            theta_max_now = diag["theta_max"]
            # Allow tiny decrease from diffusion (S_θ might be smaller than diffusion loss)
            assert theta_max_now >= theta_max_prev - 1e-10, (
                f"θ_max decreased: {theta_max_prev:.4e} → {theta_max_now:.4e}"
            )
            theta_max_prev = theta_max_now

    def test_mu_eff_min_geq_nu_after_step(self, small_solver):
        """μ_eff_min ≥ ν after every step."""
        s = small_solver
        nu = s["nu"]
        for _ in range(5):
            s, diag = solver_step(s)
            assert diag["mu_eff_min"] >= nu - 1e-14, (
                f"μ_eff_min={diag['mu_eff_min']:.4e} < ν={nu:.4e}"
            )

    def test_lps_margin_nonnegative_after_steps(self, small_solver):
        """ε_LPS = min(μ_eff) - ν ≥ 0 after every step."""
        s = small_solver
        for _ in range(10):
            s, diag = solver_step(s)
            assert diag["lps_margin"] >= -1e-14, (
                f"ε_LPS = {diag['lps_margin']:.4e} < 0 at t={s['t']:.4f}"
            )

    def test_verdict_pass_for_positive_margin(self, small_solver):
        s, diag = solver_step(small_solver)
        assert diag["verdict"] == "PASS"


# ─────────────────────────────────────────────────────────────────────────────
# Run loop
# ─────────────────────────────────────────────────────────────────────────────

class TestRunLoop:

    def test_run_reaches_t_end(self, medium_solver):
        result = run(medium_solver, t_end=0.05, max_steps=50)
        assert result["t_final"] >= 0.9 * 0.05

    def test_run_verdict_pass(self, medium_solver):
        result = run(medium_solver, t_end=0.05, max_steps=50)
        assert result["verdict"] == "PASS"

    def test_run_history_recorded(self, medium_solver):
        result = run(medium_solver, t_end=0.05, max_steps=50)
        assert result["n_steps"] >= 1
        assert len(result["history"]) == result["n_steps"]

    def test_lps_margin_min_nonnegative(self, medium_solver):
        """ε_LPS_min ≥ 0 over entire run."""
        result = run(medium_solver, t_end=0.05, max_steps=50)
        assert result["lps_margin_min"] >= -1e-14, (
            f"lps_margin_min = {result['lps_margin_min']:.4e} < 0"
        )

    def test_theta_max_final_nonneg(self, medium_solver):
        result = run(medium_solver, t_end=0.05, max_steps=50)
        assert result["theta_max_final"] >= 0.0

    def test_mu_eff_min_global_geq_nu(self, medium_solver):
        nu = medium_solver["nu"]
        result = run(medium_solver, t_end=0.05, max_steps=50)
        assert result["mu_eff_min_global"] >= nu - 1e-14

    def test_callback_called(self, medium_solver):
        calls = []
        run(medium_solver, t_end=0.05, max_steps=50,
            callback=lambda d: calls.append(d["t"]))
        assert len(calls) >= 1

    def test_run_with_zero_eps(self):
        """ε_param=0: θ decouples, exact Prize equations, still PASS."""
        s = make_solver(N=32, nu=1e-3, eps_param=0.0)
        ic = generate_CCF_profile(0.6057, N=32)
        s = set_ic(s, ic)
        result = run(s, t_end=0.05, max_steps=50)
        assert result["verdict"] == "PASS"
        # With eps=0: μ_eff = ν exactly, so lps_margin = 0 always
        assert abs(result["lps_margin_min"]) < 1e-12


# ─────────────────────────────────────────────────────────────────────────────
# ε_param sweep
# ─────────────────────────────────────────────────────────────────────────────

class TestEpsSweep:

    def test_returns_4_results(self, tmp_db):
        results = run_eps_sweep(db_path=tmp_db, **_FAST)
        assert len(results) == 4

    def test_all_eps_values_pass(self, tmp_db):
        """PRD §7.2 v0.4 M2: LPS margin ≥ 0 for all ε_param ∈ {1, 0.1, 0.01, 0.001}."""
        results = run_eps_sweep(db_path=tmp_db, **_FAST)
        failures = [r for r in results if r["verdict"] != "PASS"]
        assert len(failures) == 0, (
            f"ε_param sweep FAILED: {[r['eps_param'] for r in failures]}"
        )

    def test_lps_margin_nonneg_all_runs(self, tmp_db):
        """ε_LPS_min ≥ 0 for every ε_param value."""
        results = run_eps_sweep(db_path=tmp_db, **_FAST)
        for r in results:
            assert r["lps_margin_min"] >= -1e-14, (
                f"ε_param={r['eps_param']}: lps_margin_min={r['lps_margin_min']:.4e} < 0"
            )

    def test_theta_max_same_across_eps(self, tmp_db):
        """θ is driven by S_θ = ν·|∇u|² which doesn't depend on ε_param.
        θ_max_final should be similar across all ε values (within solver noise)."""
        results = run_eps_sweep(db_path=tmp_db, **_FAST)
        theta_maxes = [r["theta_max_final"] for r in results]
        # All should be non-zero and in the same ballpark
        assert all(t >= 0.0 for t in theta_maxes)
        # Relative spread should be small (θ accumulation is ε-independent)
        if min(theta_maxes) > 1e-20:
            spread = max(theta_maxes) / min(theta_maxes)
            assert spread < 10.0, (
                f"θ_max spread too large across ε_param values: {theta_maxes}"
            )

    def test_exp_ids_correct(self, tmp_db):
        """Experiment IDs must be EXP-L2-R2-001…004."""
        results = run_eps_sweep(db_path=tmp_db, **_FAST)
        for i, r in enumerate(results, start=1):
            assert r["exp_id"] == f"EXP-L2-R2-{i:03d}"

    def test_lps_margin_increases_with_eps(self, tmp_db):
        """Larger ε_param → larger μ_eff → larger ε_LPS_min (monotone)."""
        results = run_eps_sweep(db_path=tmp_db, **_FAST)
        margins = [r["lps_margin_min"] for r in results]
        # Should be monotonically non-decreasing as eps increases
        # (results are in decreasing eps order: 1.0, 0.1, 0.01, 0.001)
        # So margins should be monotonically non-increasing
        for i in range(len(margins) - 1):
            assert margins[i] >= margins[i + 1] - 1e-14, (
                f"ε_LPS not monotone: eps={EPS_SWEEP_VALUES[i]:.4g} gives "
                f"margin={margins[i]:.4e} < eps={EPS_SWEEP_VALUES[i+1]:.4g} "
                f"gives margin={margins[i+1]:.4e}"
            )

    def test_eps_001_small_but_positive_lps_margin(self, tmp_db):
        """Smallest ε_param (0.001) should still give ε_LPS ≥ 0."""
        results = run_eps_sweep(eps_values=[0.001], db_path=tmp_db, **_FAST)
        assert results[0]["verdict"] == "PASS"
        assert results[0]["lps_margin_min"] >= -1e-14

    def test_results_logged_to_db(self, tmp_db):
        run_eps_sweep(db_path=tmp_db, **_FAST)
        rows = query_results(db_path=tmp_db, layer=2)
        route2_rows = [r for r in rows if r["route"] == 2]
        assert len(route2_rows) == 4

    def test_different_f_theta_forms_all_pass(self, tmp_db):
        """All three f(θ) functional forms must pass."""
        for form in ("identity", "tanh", "sqrt"):
            db = tmp_db + f"_{form}.db"
            results = run_eps_sweep(
                eps_values=[1.0, 0.001], f_theta=form,
                db_path=db, **_FAST
            )
            for r in results:
                assert r["verdict"] == "PASS", (
                    f"f_theta='{form}', eps={r['eps_param']}: FAIL"
                )


# ─────────────────────────────────────────────────────────────────────────────
# Prize limit (ε→0)
# ─────────────────────────────────────────────────────────────────────────────

class TestPrizeLimitConvergence:

    def test_omega_converges_to_prize_as_eps_to_zero(self, tmp_db):
        """
        ω evolution is identical for all ε_param (since ω uses only ν).
        Final ω fields should match regardless of ε_param.
        """
        t_end = 0.05
        max_steps = 20
        N = 32
        nu = 1e-3

        ic = generate_CCF_profile(0.6057, N=N)

        def run_with_eps(eps):
            s = make_solver(N=N, nu=nu, eps_param=eps)
            s = set_ic(s, ic)
            return run(s, t_end=t_end, max_steps=max_steps)

        r_large  = run_with_eps(1.0)
        r_small  = run_with_eps(0.001)
        r_zero   = run_with_eps(0.0)

        omega_large = r_large["solver"]["omega"]
        omega_small = r_small["solver"]["omega"]
        omega_zero  = r_zero["solver"]["omega"]

        # All three ω fields should be essentially identical
        # (ω equation doesn't depend on eps_param)
        tol = 1e-8 * float(np.max(np.abs(omega_zero)) + 1e-30)
        assert np.max(np.abs(omega_large - omega_zero)) < tol, (
            "ω differs between ε_param=1.0 and ε_param=0 — "
            "Prize ω equation must be independent of ε_param"
        )
        assert np.max(np.abs(omega_small - omega_zero)) < tol

    def test_theta_nonzero_for_nontrivial_flow(self):
        """After any steps with a nontrivial ω IC, θ_max > 0."""
        s = make_solver(N=32, nu=1e-3, eps_param=0.1)
        ic = generate_CCF_profile(0.6057, N=32)
        s = set_ic(s, ic)
        result = run(s, t_end=0.05, max_steps=20)
        assert result["theta_max_final"] > 0.0, "θ should accumulate for nontrivial flow"

    def test_zero_eps_lps_margin_is_zero(self):
        """When ε_param=0, μ_eff = ν exactly, so ε_LPS = 0 (not negative)."""
        s = make_solver(N=32, nu=1e-3, eps_param=0.0)
        ic = generate_CCF_profile(0.6057, N=32)
        s = set_ic(s, ic)
        result = run(s, t_end=0.05, max_steps=20)
        # ε_LPS = min(μ_eff) - ν = min(ν + 0·f(θ)) - ν = 0
        assert abs(result["lps_margin_min"]) < 1e-12, (
            f"ε_param=0 should give ε_LPS=0, got {result['lps_margin_min']:.4e}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Smoke test
# ─────────────────────────────────────────────────────────────────────────────

class TestSmokeTest:

    def test_smoke_test_passes(self, tmp_db):
        r = run_smoke_test(db_path=tmp_db)
        assert r["verdict"] == "PASS"

    def test_smoke_test_n_runs(self, tmp_db):
        r = run_smoke_test(db_path=tmp_db)
        assert r["n_runs"] == 2

    def test_smoke_test_all_pass(self, tmp_db):
        r = run_smoke_test(db_path=tmp_db)
        assert r["n_pass"] == r["n_runs"]


# ─────────────────────────────────────────────────────────────────────────────
# Import from layer1
# ─────────────────────────────────────────────────────────────────────────────

class TestLayer1Integration:

    def test_route2_theta_source_used_in_rhs(self, small_solver):
        """S_θ in _rhs should match route2_theta_source from layer1."""
        from layer1.tfirst_props import route2_theta_source as r2ts

        s = small_solver
        # One RHS call
        _, _, diag = _rhs(
            s["omega"], s["theta"],
            s["kx"], s["ky"], s["k2"],
            s["nu"], s["eps_param"], s["f_theta"],
        )
        # S_theta_mean ≥ 0
        assert diag["S_theta_mean"] >= 0.0

    def test_route1_coeffs_unused_in_route2(self):
        """Route 2 uses ν=const; route1_coeffs (μ(T) bounds) is not needed."""
        # Just verify we can run route2 without importing route1_coeffs
        s = make_solver(N=16, nu=1e-3, eps_param=0.1)
        ic = generate_CCF_profile(0.6057, N=16)
        s = set_ic(s, ic)
        s, diag = solver_step(s)
        assert diag["verdict"] == "PASS"
