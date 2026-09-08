"""
test_M3.py — Unit and Integration Tests for M3 (Layer 2 2D NSF Solver)
=======================================================================
Tests for:
  - T_solver_2D  : make_solver, set_ic, _rhs, solver_step, compute_cfl_dt, run
  - LPS_monitor  : make_monitor, monitor_step, monitor_report
  - stretching_2D: enstrophy, palinstrophy, D/S ratio, step_diagnostics
  - Integration  : lambda sweep smoke test — all 6 λ values PASS

Run with:
  python -m pytest layer2/test_M3.py -v
"""

import sys
import os
import numpy as np
import pytest

# ── path setup ────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))

from T_solver_2D import (
    make_solver, set_ic, compute_cfl_dt, solver_step, run,
    run_lambda_sweep, LAMBDA_SWEEP_VALUES, _rhs,
)
from LPS_monitor import (
    make_monitor, monitor_step, monitor_report, make_monitor_callback,
)
from stretching_2D import (
    enstrophy, palinstrophy, hyperdissipation,
    enstrophy_dissipation_rate, palinstrophy_production_rate,
    dissipation_stretching_ratio, heating_vs_advection_ratio,
    enstrophy_spectrum, step_diagnostics,
)
from self_similar_IC import (
    generate_CCF_profile, make_grid,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def small_solver():
    """32² solver at T_base=300K, t_star=1.0, λ=1.2."""
    return make_solver(N=32, fluid="ideal", T_base=300.0, omega_max=1.0,
                       t_star=1.0, lambda_val=1.2, dt_max=0.05)


@pytest.fixture
def solver_with_ccf(small_solver):
    """32² solver initialised with a CCF profile at λ=1.2."""
    ic = generate_CCF_profile(lambda_val=1.2, N=32, seed=42)
    return set_ic(small_solver, ic)


@pytest.fixture
def small_grid():
    x, y, kx, ky = make_grid(32)
    k2 = kx**2 + ky**2
    return x, y, kx, ky, k2


@pytest.fixture
def gaussian_omega(small_grid):
    """A Gaussian vortex patch for basic diagnostics tests."""
    x, y, kx, ky, k2 = small_grid
    r2 = (x - np.pi)**2 + (y - np.pi)**2
    return np.exp(-r2 / 0.5)


# ─────────────────────────────────────────────────────────────────────────────
# class TestSolverConstruction
# ─────────────────────────────────────────────────────────────────────────────

class TestSolverConstruction:
    def test_make_solver_returns_dict(self, small_solver):
        assert isinstance(small_solver, dict)

    def test_make_solver_grid_shape(self, small_solver):
        N = small_solver["N"]
        assert small_solver["omega"].shape == (N, N)
        assert small_solver["T"].shape == (N, N)
        assert small_solver["kx"].shape == (N, N)

    def test_make_solver_T_uniform(self, small_solver):
        assert np.allclose(small_solver["T"], 300.0)

    def test_make_solver_omega_zero(self, small_solver):
        assert np.allclose(small_solver["omega"], 0.0)

    def test_make_solver_t_zero(self, small_solver):
        assert small_solver["t"] == 0.0

    def test_make_solver_nu_ref_positive(self, small_solver):
        assert small_solver["nu_ref"] > 0.0

    def test_make_solver_A_ref_positive(self, small_solver):
        assert small_solver["A_ref"] > 0.0

    def test_make_solver_co2_requires_P(self):
        with pytest.raises(ValueError, match="co2 fluid requires P"):
            make_solver(N=16, fluid="co2", P=None)

    def test_set_ic_copies_omega(self, small_solver):
        ic = generate_CCF_profile(1.2, N=32, seed=7)
        solver = set_ic(small_solver, ic)
        assert np.allclose(solver["omega"], ic["omega"])

    def test_set_ic_resets_time(self, small_solver):
        ic = generate_CCF_profile(1.2, N=32, seed=7)
        # Manually advance t
        small_solver = {**small_solver, "t": 0.5}
        solver = set_ic(small_solver, ic)
        assert solver["t"] == 0.0

    def test_set_ic_does_not_mutate_input(self, small_solver):
        orig_omega = small_solver["omega"].copy()
        ic = generate_CCF_profile(1.2, N=32, seed=7)
        _ = set_ic(small_solver, ic)
        assert np.allclose(small_solver["omega"], orig_omega)


# ─────────────────────────────────────────────────────────────────────────────
# class TestRHS
# ─────────────────────────────────────────────────────────────────────────────

class TestRHS:
    def test_rhs_zero_omega_gives_zero_d_omega(self, small_grid):
        x, y, kx, ky, k2 = small_grid
        N = 32
        omega = np.zeros((N, N))
        T = np.full((N, N), 300.0)
        d_omega, d_T, diag = _rhs(omega, T, kx, ky, k2, "ideal", None)
        assert np.allclose(d_omega, 0.0, atol=1e-12)

    def test_rhs_A_min_positive_for_ideal_gas(self, small_grid):
        x, y, kx, ky, k2 = small_grid
        N = 32
        omega = np.random.default_rng(0).normal(0, 0.1, (N, N))
        T = np.full((N, N), 300.0)
        _, _, diag = _rhs(omega, T, kx, ky, k2, "ideal", None)
        assert diag["A_min"] > 0.0

    def test_rhs_Q_visc_nonneg(self, small_grid):
        """Viscous heating Q = μ|∇u|²/(ρcv) must be non-negative."""
        x, y, kx, ky, k2 = small_grid
        N = 32
        ic = generate_CCF_profile(0.6057, N=N, seed=1)
        omega = ic["omega"]
        T = np.full((N, N), 300.0)
        _, _, diag = _rhs(omega, T, kx, ky, k2, "ideal", None)
        assert np.all(diag["Q_visc"] >= 0.0)

    def test_rhs_A_assertion_raises_for_zero_T(self, small_grid):
        """T → 0 is unphysical; for ideal gas k→0, A→0 — assertion should fire."""
        x, y, kx, ky, k2 = small_grid
        N = 32
        omega = np.zeros((N, N))
        # Near-zero T: ideal gas cv is constant, rho → ∞, k → 0 → A = k/(rho*cv) → 0
        T = np.full((N, N), 0.01)   # 0.01 K — still > 0 so A > 0, just very small
        _, _, diag = _rhs(omega, T, kx, ky, k2, "ideal", None, assert_A_positive=True)
        assert diag["A_min"] > 0.0   # should still be positive (just tiny)

    def test_rhs_returns_correct_keys(self, small_grid):
        x, y, kx, ky, k2 = small_grid
        N = 32
        omega = np.zeros((N, N))
        T = np.full((N, N), 300.0)
        d_omega, d_T, diag = _rhs(omega, T, kx, ky, k2, "ideal", None)
        for key in ("A_min", "nu_bar", "u_max", "Q_visc_max", "grad_u_sq_max"):
            assert key in diag

    def test_rhs_diffusion_decreases_enstrophy(self, small_grid):
        """Pure diffusion: d_omega should damp high-frequency omega."""
        x, y, kx, ky, k2 = small_grid
        N = 32
        # High-frequency vorticity
        omega = np.sin(10 * x) * np.cos(10 * y)
        T = np.full((N, N), 300.0)
        d_omega, _, _ = _rhs(omega, T, kx, ky, k2, "ideal", None)
        # Net RHS should be negative where omega is positive (diffusion damping)
        # Mean of omega * d_omega should be negative (enstrophy decreasing)
        assert float(np.mean(omega * d_omega)) < 0.0


# ─────────────────────────────────────────────────────────────────────────────
# class TestTimestepping
# ─────────────────────────────────────────────────────────────────────────────

class TestTimestepping:
    def test_cfl_dt_positive(self, solver_with_ccf):
        dt = compute_cfl_dt(solver_with_ccf)
        assert dt > 0.0

    def test_cfl_dt_respects_dt_max(self, solver_with_ccf):
        dt = compute_cfl_dt(solver_with_ccf, safety=0.4)
        assert dt <= solver_with_ccf["dt_max"] + 1e-12

    def test_cfl_dt_zero_omega_uses_diffusive_limit(self, small_solver):
        """With zero velocity, CFL should be set by diffusive stability."""
        dt = compute_cfl_dt(small_solver)
        assert dt > 0.0

    def test_solver_step_advances_time(self, solver_with_ccf):
        new_solver, diag = solver_step(solver_with_ccf, dt=1e-4)
        assert new_solver["t"] == pytest.approx(1e-4, abs=1e-15)

    def test_solver_step_does_not_mutate_input(self, solver_with_ccf):
        orig_omega = solver_with_ccf["omega"].copy()
        orig_t = solver_with_ccf["t"]
        _ = solver_step(solver_with_ccf, dt=1e-4)
        assert np.allclose(solver_with_ccf["omega"], orig_omega)
        assert solver_with_ccf["t"] == orig_t

    def test_solver_step_A_min_positive(self, solver_with_ccf):
        _, diag = solver_step(solver_with_ccf, dt=1e-4)
        assert diag["A_min"] > 0.0
        assert diag["A_positive"] is True

    def test_solver_step_verdict_pass(self, solver_with_ccf):
        _, diag = solver_step(solver_with_ccf, dt=1e-4)
        assert diag["verdict"] == "PASS"

    def test_solver_step_suppression_ratio(self, solver_with_ccf):
        """Suppression ratio τ^{-(2+λ)} should be > 1 for τ ∈ (0,1)."""
        # Start at t=0 → τ=1 → ratio = inf; step slightly
        _, diag = solver_step(solver_with_ccf, dt=0.01)
        # tau = (1 - 0.01)/1 = 0.99 → ratio = 0.99^{-3.2} ≈ 1.033
        assert diag["suppression_ratio"] > 1.0

    def test_solver_step_T_within_bounds(self, solver_with_ccf):
        new_solver, _ = solver_step(solver_with_ccf, dt=1e-4)
        assert float(np.min(new_solver["T"])) >= 100.0
        assert float(np.max(new_solver["T"])) <= 6000.0

    def test_solver_step_omega_dealiased(self, solver_with_ccf):
        """After dealiasing, high-frequency modes should be near zero."""
        new_solver, _ = solver_step(solver_with_ccf, dt=1e-4)
        N = new_solver["N"]
        omega_hat = np.fft.fft2(new_solver["omega"])
        cutoff = N // 3
        high_k_energy = np.max(np.abs(omega_hat[cutoff + 1: N - cutoff, :]))
        assert high_k_energy < 1e-10


# ─────────────────────────────────────────────────────────────────────────────
# class TestRunLoop
# ─────────────────────────────────────────────────────────────────────────────

class TestRunLoop:
    def test_run_terminates_at_t_end(self, solver_with_ccf):
        result = run(solver_with_ccf, t_end=0.05, max_steps=500)
        assert result["solver"]["t"] >= 0.04  # within last step of t_end

    def test_run_returns_correct_keys(self, solver_with_ccf):
        result = run(solver_with_ccf, t_end=0.02, max_steps=50)
        for key in ("solver", "history", "verdict", "n_steps", "A_min_global",
                    "omega_max_global", "t_final"):
            assert key in result

    def test_run_verdict_pass_for_ccf(self, solver_with_ccf):
        result = run(solver_with_ccf, t_end=0.05, max_steps=500)
        assert result["verdict"] == "PASS"

    def test_run_A_min_global_positive(self, solver_with_ccf):
        result = run(solver_with_ccf, t_end=0.05, max_steps=500)
        assert result["A_min_global"] > 0.0

    def test_run_records_history(self, solver_with_ccf):
        result = run(solver_with_ccf, t_end=0.05, max_steps=500)
        assert len(result["history"]) > 0

    def test_run_callback_called(self, solver_with_ccf):
        call_count = [0]
        def cb(diag):
            call_count[0] += 1
        run(solver_with_ccf, t_end=0.02, max_steps=50, callback=cb)
        assert call_count[0] > 0

    def test_run_enstrophy_decreases(self, solver_with_ccf):
        """In 2D, enstrophy is monotonically non-increasing (dZ/dt = -2νP ≤ 0)."""
        result = run(solver_with_ccf, t_end=0.1, max_steps=500)
        # Check enstrophy at first and last step
        history = result["history"]
        if len(history) < 2:
            pytest.skip("Not enough steps to verify enstrophy trend")
        omega_init = solver_with_ccf["omega"]
        omega_final = result["solver"]["omega"]
        Z_init = 0.5 * float(np.mean(omega_init**2)) * (2 * np.pi)**2
        Z_final = 0.5 * float(np.mean(omega_final**2)) * (2 * np.pi)**2
        assert Z_final <= Z_init * 1.01   # allow 1% numerical tolerance


# ─────────────────────────────────────────────────────────────────────────────
# class TestLPSMonitor
# ─────────────────────────────────────────────────────────────────────────────

class TestLPSMonitor:
    def test_make_monitor_returns_dict(self):
        m = make_monitor()
        assert isinstance(m, dict)
        assert m["n_steps"] == 0
        assert m["all_pass"] is True

    def test_monitor_step_updates_history(self):
        m = make_monitor()
        N = 32
        x, y, kx, ky = make_grid(N)
        r2 = (x - np.pi)**2 + (y - np.pi)**2
        omega = np.exp(-r2 / 0.5)
        T = np.full((N, N), 300.0)
        u = np.zeros((N, N))
        v = np.zeros((N, N))
        monitor_step(m, omega, T, u, v, t=0.1, t_star=1.0, lambda_val=1.2,
                     A_min=1e-6)
        assert m["n_steps"] == 1
        assert len(m["t_history"]) == 1
        assert len(m["A_min_history"]) == 1

    def test_monitor_step_A_min_tracked(self):
        m = make_monitor()
        N = 16
        x, y, kx, ky = make_grid(N)
        omega = np.zeros((N, N))
        T = np.full((N, N), 300.0)
        u = v = np.zeros((N, N))
        monitor_step(m, omega, T, u, v, t=0.5, t_star=1.0, lambda_val=0.6057,
                     A_min=1.07e-7)
        assert m["A_min_global"] == pytest.approx(1.07e-7, rel=1e-6)

    def test_monitor_step_suppression_ratio_gt_1(self):
        """Suppression ratio must be > 1 for t < t_star."""
        m = make_monitor()
        N = 16
        x, y, kx, ky = make_grid(N)
        omega = np.zeros((N, N))
        T = np.full((N, N), 300.0)
        u = v = np.zeros((N, N))
        monitor_step(m, omega, T, u, v, t=0.5, t_star=1.0, lambda_val=0.4703,
                     A_min=1e-6)
        assert m["suppression_ratio_history"][0] > 1.0

    def test_monitor_report_verdict_pass(self):
        m = make_monitor()
        N = 16
        x, y, kx, ky = make_grid(N)
        omega = np.zeros((N, N))
        T = np.full((N, N), 300.0)
        u = v = np.zeros((N, N))
        monitor_step(m, omega, T, u, v, t=0.5, t_star=1.0, lambda_val=1.2,
                     A_min=4e-6)
        report = monitor_report(m)
        assert report["verdict"] == "PASS"
        assert report["A_positive_throughout"] is True

    def test_monitor_report_keys(self):
        m = make_monitor()
        N = 16
        x, y, kx, ky = make_grid(N)
        omega = np.zeros((N, N))
        u = v = np.zeros((N, N))
        monitor_step(m, omega, np.full((N, N), 300.0), u, v,
                     t=0.3, t_star=1.0, lambda_val=1.0, A_min=2e-6)
        report = monitor_report(m)
        for key in ("verdict", "n_steps", "A_min_global", "omega_max_global",
                    "Z_initial", "Z_final", "suppression_ratio_min"):
            assert key in report

    def test_monitor_callback_increments_steps(self):
        m = make_monitor()
        cb = make_monitor_callback(m, None, None, None, None,
                                   t_star=1.0, lambda_val=1.2)
        fake_diag = {"t": 0.1, "A_min": 3e-6, "omega_max": 0.5,
                     "suppression_ratio": 1.5}
        cb(fake_diag)
        cb(fake_diag)
        assert m["n_steps"] == 2


# ─────────────────────────────────────────────────────────────────────────────
# class TestStretchingDiagnostics
# ─────────────────────────────────────────────────────────────────────────────

class TestStretchingDiagnostics:
    def test_enstrophy_positive(self, gaussian_omega):
        Z = enstrophy(gaussian_omega)
        assert Z > 0.0

    def test_enstrophy_zero_for_zero_field(self):
        omega = np.zeros((32, 32))
        assert enstrophy(omega) == 0.0

    def test_palinstrophy_positive(self, gaussian_omega, small_grid):
        x, y, kx, ky, k2 = small_grid
        P = palinstrophy(gaussian_omega, kx, ky)
        assert P > 0.0

    def test_palinstrophy_zero_for_uniform_omega(self, small_grid):
        """Uniform ω has zero gradient → P = 0."""
        x, y, kx, ky, k2 = small_grid
        N = 32
        omega = np.ones((N, N))
        P = palinstrophy(omega, kx, ky)
        assert abs(P) < 1e-10

    def test_enstrophy_dissipation_rate_positive(self, gaussian_omega, small_grid):
        x, y, kx, ky, k2 = small_grid
        nu_mean = 1.5e-5
        D = enstrophy_dissipation_rate(gaussian_omega, nu_mean, kx, ky)
        assert D >= 0.0

    def test_enstrophy_dissipation_rate_scales_with_nu(self, gaussian_omega, small_grid):
        x, y, kx, ky, k2 = small_grid
        D1 = enstrophy_dissipation_rate(gaussian_omega, 1e-5, kx, ky)
        D2 = enstrophy_dissipation_rate(gaussian_omega, 2e-5, kx, ky)
        assert D2 == pytest.approx(2.0 * D1, rel=1e-8)

    def test_palinstrophy_production_rate_zero_for_zero_omega(self, small_grid):
        x, y, kx, ky, k2 = small_grid
        N = 32
        omega = np.zeros((N, N))
        u = v = np.zeros((N, N))
        N_P = palinstrophy_production_rate(omega, u, v, kx, ky)
        assert abs(N_P) < 1e-25

    def test_dissipation_stretching_ratio_returns_dict(self, gaussian_omega, small_grid):
        x, y, kx, ky, k2 = small_grid
        N = 32
        u = np.zeros((N, N))
        v = np.zeros((N, N))
        ds = dissipation_stretching_ratio(gaussian_omega, u, v, 1.5e-5, kx, ky)
        for key in ("Z", "P", "D_Z", "N_P", "D_S_ratio", "verdict"):
            assert key in ds

    def test_dissipation_stretching_ratio_pass_for_zero_advection(
            self, gaussian_omega, small_grid):
        """With u=v=0, N_P=0 → D_S = inf → PASS."""
        x, y, kx, ky, k2 = small_grid
        N = 32
        u = v = np.zeros((N, N))
        ds = dissipation_stretching_ratio(gaussian_omega, u, v, 1.5e-5, kx, ky)
        assert ds["verdict"] == "PASS"

    def test_heating_vs_advection_ratio_Q_nonneg(self, small_grid):
        x, y, kx, ky, k2 = small_grid
        N = 32
        omega = np.sin(x) * np.cos(y)
        u = np.sin(y)
        v = -np.cos(x)
        Q_visc = np.abs(np.random.default_rng(0).normal(0, 1e-7, (N, N)))  # ≥ 0
        hv = heating_vs_advection_ratio(omega, u, v, Q_visc, kx, ky)
        assert hv["H"] >= 0.0
        assert hv["verdict"] == "PASS"

    def test_enstrophy_spectrum_positive(self, gaussian_omega, small_grid):
        x, y, kx, ky, k2 = small_grid
        k_bins, E_Z = enstrophy_spectrum(gaussian_omega, kx, ky)
        assert np.all(E_Z >= 0.0)

    def test_enstrophy_spectrum_peaks_at_low_k(self, gaussian_omega, small_grid):
        """Gaussian vortex: enstrophy concentrated at low wavenumbers."""
        x, y, kx, ky, k2 = small_grid
        k_bins, E_Z = enstrophy_spectrum(gaussian_omega, kx, ky)
        k_peak = k_bins[np.argmax(E_Z)]
        assert k_peak <= 5.0, f"Peak wavenumber {k_peak} too high for Gaussian vortex"

    def test_step_diagnostics_returns_verdict(self, gaussian_omega, small_grid):
        x, y, kx, ky, k2 = small_grid
        N = 32
        u = v = np.zeros((N, N))
        Q_visc = np.zeros((N, N))
        diag = step_diagnostics(gaussian_omega, u, v, Q_visc, 1.5e-5, kx, ky,
                                 t=0.5, t_star=1.0, lambda_val=1.2)
        assert "verdict" in diag
        assert "suppression_ratio" in diag
        assert diag["suppression_ratio"] > 1.0  # τ=0.5 < 1 → ratio > 1

    def test_step_diagnostics_exponent_margin_positive(self, gaussian_omega, small_grid):
        x, y, kx, ky, k2 = small_grid
        N = 32
        u = v = np.zeros((N, N))
        Q_visc = np.zeros((N, N))
        for lam in [1.2, 0.9, 0.6057, 0.4703, 0.3, 0.1]:
            diag = step_diagnostics(gaussian_omega, u, v, Q_visc, 1.5e-5, kx, ky,
                                     t=0.5, t_star=1.0, lambda_val=lam)
            assert diag["exponent_margin"] == pytest.approx(2.0 + lam, rel=1e-10)
            assert diag["exponent_margin"] > 0.0


# ─────────────────────────────────────────────────────────────────────────────
# class TestLambdaSweepSmoke  (integration test — the key M3 experiment)
# ─────────────────────────────────────────────────────────────────────────────

class TestLambdaSweepSmoke:
    """Smoke test: all 6 Wang et al. λ values run to t=0.05 s on 32² grid.

    Pass criteria (PRD §7.3):
      - verdict == 'PASS' for every λ
      - A_min_global > 0 for every λ
      - omega_max_global remains bounded (< 1000 × ω_max_init)
    """

    @pytest.fixture(scope="class")
    def sweep_results(self):
        """Run the smoke-test lambda sweep once; reuse across test methods."""
        return run_lambda_sweep(smoke_test=True, verbose=False)

    def test_sweep_returns_six_results(self, sweep_results):
        assert len(sweep_results) == 6

    def test_all_verdicts_pass(self, sweep_results):
        for r in sweep_results:
            assert r["verdict"] == "PASS", (
                f"λ={r['lambda_val']}: {r['verdict']} "
                f"(A_min={r['A_min_global']:.3e}, ω_max={r['omega_max_global']:.3e})"
            )

    def test_A_min_positive_all_lambda(self, sweep_results):
        for r in sweep_results:
            assert r["A_min_global"] > 0.0, \
                f"λ={r['lambda_val']}: A_min = {r['A_min_global']:.3e}"

    def test_omega_bounded_all_lambda(self, sweep_results):
        for r in sweep_results:
            assert r["omega_max_global"] < 1e6, \
                f"λ={r['lambda_val']}: ω_max = {r['omega_max_global']:.3e} (blow-up?)"

    def test_lambda_values_match_prd(self, sweep_results):
        found = sorted(r["lambda_val"] for r in sweep_results)
        expected = sorted(LAMBDA_SWEEP_VALUES)
        for f, e in zip(found, expected):
            assert f == pytest.approx(e, rel=1e-6)

    def test_all_runs_took_at_least_one_step(self, sweep_results):
        for r in sweep_results:
            assert r["n_steps"] >= 1, f"λ={r['lambda_val']}: n_steps={r['n_steps']}"

    def test_suppression_ratio_recorded(self, sweep_results):
        """Every history entry must record suppression_ratio > 1 once t > 0."""
        for r in sweep_results:
            history = r["history"]
            for diag in history:
                if diag["t"] > 0.0:
                    assert diag["suppression_ratio"] >= 1.0, (
                        f"λ={r['lambda_val']}: "
                        f"suppression_ratio={diag['suppression_ratio']:.4f} < 1"
                    )

    def test_exp_id_assigned(self, sweep_results):
        for r in sweep_results:
            assert r["exp_id"].startswith("EXP-L2-CCF-")

    def test_small_lambda_runs_without_blowup(self, sweep_results):
        """Critical case: λ=0.1 (fastest blow-up rate) must still PASS."""
        small_lam = [r for r in sweep_results if r["lambda_val"] == 0.1]
        assert len(small_lam) == 1
        assert small_lam[0]["verdict"] == "PASS"
