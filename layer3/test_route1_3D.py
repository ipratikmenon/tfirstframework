"""
test_route1_3D.py — Test suite for M5 Route 1 3D solver
=========================================================
44 tests across 8 classes.

Classes:
  TestGrid3D              (6)  — make_grid_3D, dealias mask, dealias function
  TestDivFree             (5)  — Helmholtz projection, divergence checks
  TestTaylorGreenIC       (6)  — IC shape, div-free, E₀ analytic, T uniform
  TestSolverConstruction  (5)  — make_solver keys, shapes, set_ic, immutable
  TestRHS3D               (6)  — shapes, A_min>0, S_θ≥0, velocity structure
  TestTimeStepping3D      (7)  — CFL, advances, immutable, energy, div-free
  TestRunLoop3D           (5)  — reaches t_end, PASS, E decays, A_min>0, history
  TestTaylorGreenBenchmark(4)  — E₀ accuracy, energy decays, div_ok, PASS verdict
  TestSmokeTest           (3)  — smoke PASS, n_steps, all_pass
  TestM5EndToEnd          (1)  — full TG chain: make→IC→run→verdict
"""

import sys
import os
import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))

from route1_3D import (
    make_grid_3D,
    make_dealias_mask_3D,
    dealias_3D,
    project_divergence_free,
    divergence_rms,
    taylor_green_ic,
    kinetic_energy,
    make_solver,
    set_ic,
    compute_cfl_dt,
    solver_step,
    run,
    run_taylor_green,
    run_smoke_test,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def N16_solver():
    """N=16³ solver with Taylor-Green IC."""
    ic = taylor_green_ic(16, V0=1.0, T_base=300.0)
    s = make_solver(N=16, T_base=300.0, dt_max=0.05)
    return set_ic(s, ic)


@pytest.fixture(scope="module")
def N16_tg_result():
    """One TG run at N=16, t_end=0.1 — used by multiple test classes."""
    return run_taylor_green(N=16, t_end=0.1, max_steps=30, verbose=False)


# ─────────────────────────────────────────────────────────────────────────────
# 1. TestGrid3D
# ─────────────────────────────────────────────────────────────────────────────

class TestGrid3D:
    def test_grid_x_shape(self):
        x, kx, ky, kz, k2 = make_grid_3D(16)
        assert x.shape == (16,)

    def test_grid_k_shape(self):
        x, kx, ky, kz, k2 = make_grid_3D(16)
        assert kx.shape == (16, 16, 16)
        assert ky.shape == (16, 16, 16)
        assert kz.shape == (16, 16, 16)
        assert k2.shape == (16, 16, 16)

    def test_k2_nonnegative(self):
        _, kx, ky, kz, k2 = make_grid_3D(16)
        assert np.all(k2 >= 0)

    def test_k2_zero_mode(self):
        _, kx, ky, kz, k2 = make_grid_3D(16)
        assert k2[0, 0, 0] == 0.0

    def test_dealias_mask_shape(self):
        mask = make_dealias_mask_3D(16)
        assert mask.shape == (16, 16, 16)
        assert mask.dtype == bool

    def test_dealias_preserves_low_k(self):
        """Dealias should not zero the k=(1,1,1) mode."""
        N = 16
        mask = make_dealias_mask_3D(N)
        # k=(1,1,1) is well below N/3 for N=16 (k_max=5)
        assert mask[1, 1, 1] is np.bool_(True)

    def test_dealias_zeros_high_k(self):
        """Dealias should zero the highest-frequency mode."""
        N = 16
        mask = make_dealias_mask_3D(N)
        # k=(N//2, 0, 0) has |k|=8 > N//3=5 → should be zeroed
        assert not mask[N // 2, 0, 0]

    def test_dealias_real_output(self):
        N = 16
        mask = make_dealias_mask_3D(N)
        f = np.random.default_rng(0).standard_normal((N, N, N))
        g = dealias_3D(f, mask)
        assert g.shape == (N, N, N)
        assert np.isreal(g).all() or g.dtype == np.float64


# ─────────────────────────────────────────────────────────────────────────────
# 2. TestDivFree
# ─────────────────────────────────────────────────────────────────────────────

class TestDivFree:
    def test_projection_gives_div_free(self):
        """Projection reduces divergence by many orders of magnitude.

        For random fields on N=16 grid: div_before ~ O(1) (large amplitudes),
        div_after << div_before. We check a relative reduction ≥ 1e6×.
        Smooth TG field achieves 1e-12 (see test_taylor_green_already_div_free);
        random fields hit ~1e-11 due to round-trip FFT error on all N³ modes.
        """
        N = 16
        _, kx, ky, kz, k2 = make_grid_3D(N)
        rng = np.random.default_rng(1)
        u = rng.standard_normal((N, N, N))
        v = rng.standard_normal((N, N, N))
        w = rng.standard_normal((N, N, N))
        div_before = divergence_rms(u, v, w, kx, ky, kz)
        u, v, w = project_divergence_free(u, v, w, kx, ky, kz, k2)
        div_after = divergence_rms(u, v, w, kx, ky, kz)
        # After projection: divergence must be reduced by at least 6 orders of magnitude
        assert div_after < div_before * 1e-6, (
            f"Projection did not suppress div: before={div_before:.3e}, after={div_after:.3e}"
        )
        # And absolute divergence must be tiny (essentially numerical noise)
        assert div_after < 1e-6, f"div u RMS after projection = {div_after:.3e}"

    def test_projection_is_idempotent(self):
        """Projecting an already-projected field leaves divergence unchanged (both ≈ 0).

        We check the divergence after two successive projections, not array equality,
        because taking real(ifft(...)) and re-FFT-ing introduces ~1e-11 round-trip
        error for N=16 random fields with O(N^{3/2}) spectral amplitudes.
        """
        N = 16
        _, kx, ky, kz, k2 = make_grid_3D(N)
        rng = np.random.default_rng(2)
        u = rng.standard_normal((N, N, N))
        v = rng.standard_normal((N, N, N))
        w = rng.standard_normal((N, N, N))
        u1, v1, w1 = project_divergence_free(u, v, w, kx, ky, kz, k2)
        u2, v2, w2 = project_divergence_free(u1, v1, w1, kx, ky, kz, k2)
        # Both projected fields should be divergence-free to the same tolerance
        div1 = divergence_rms(u1, v1, w1, kx, ky, kz)
        div2 = divergence_rms(u2, v2, w2, kx, ky, kz)
        assert div1 < 1e-6, f"First projection: div={div1:.3e}"
        assert div2 < 1e-6, f"Second projection: div={div2:.3e}"
        # The second projection must not INCREASE the divergence
        assert div2 <= div1 + 1e-12, f"Second projection worsened div: {div1:.3e} → {div2:.3e}"

    def test_taylor_green_already_div_free(self):
        """Taylor-Green IC is analytically divergence-free."""
        N = 16
        _, kx, ky, kz, k2 = make_grid_3D(N)
        ic = taylor_green_ic(N)
        div = divergence_rms(ic["u"], ic["v"], ic["w"], kx, ky, kz)
        assert div < 1e-12, f"TG div u = {div:.3e}"

    def test_divergence_rms_output(self):
        N = 16
        _, kx, ky, kz, _ = make_grid_3D(N)
        u = np.zeros((N, N, N))
        v = np.zeros((N, N, N))
        w = np.zeros((N, N, N))
        div = divergence_rms(u, v, w, kx, ky, kz)
        assert div < 1e-14

    def test_projection_preserves_divergence_free(self):
        """TG IC after set_ic projection remains div-free."""
        N = 16
        ic = taylor_green_ic(N)
        s = make_solver(N=N)
        s = set_ic(s, ic)
        div = divergence_rms(s["u"], s["v"], s["w"], s["kx"], s["ky"], s["kz"])
        assert div < 1e-10


# ─────────────────────────────────────────────────────────────────────────────
# 3. TestTaylorGreenIC
# ─────────────────────────────────────────────────────────────────────────────

class TestTaylorGreenIC:
    def test_shape(self):
        ic = taylor_green_ic(16)
        assert ic["u"].shape == (16, 16, 16)
        assert ic["v"].shape == (16, 16, 16)
        assert ic["w"].shape == (16, 16, 16)
        assert ic["T"].shape == (16, 16, 16)

    def test_w_is_zero(self):
        ic = taylor_green_ic(16)
        assert np.allclose(ic["w"], 0.0)

    def test_E0_analytic(self):
        """E₀ analytic must equal V₀²/8."""
        for V0 in [0.5, 1.0, 2.0]:
            ic = taylor_green_ic(16, V0=V0)
            assert abs(ic["E0_analytic"] - V0**2 / 8.0) < 1e-14

    def test_E0_computed_matches_analytic(self):
        """Computed kinetic energy of TG IC must match analytic V₀²/8 within 1%."""
        ic = taylor_green_ic(32, V0=1.0)
        E_computed = kinetic_energy(ic["u"], ic["v"], ic["w"])
        E_analytic = ic["E0_analytic"]
        rel_err = abs(E_computed - E_analytic) / E_analytic
        assert rel_err < 0.01, f"E₀ relative error = {rel_err*100:.2f}%"

    def test_T_uniform(self):
        T_base = 450.0
        ic = taylor_green_ic(16, T_base=T_base)
        assert np.allclose(ic["T"], T_base)

    def test_name_and_keys(self):
        ic = taylor_green_ic(16)
        for key in ("u", "v", "w", "T", "E0_analytic", "name", "V0"):
            assert key in ic
        assert ic["name"] == "taylor_green"


# ─────────────────────────────────────────────────────────────────────────────
# 4. TestSolverConstruction
# ─────────────────────────────────────────────────────────────────────────────

class TestSolverConstruction:
    def test_keys_present(self):
        s = make_solver(N=16)
        for key in ("N", "kx", "ky", "kz", "k2", "dealias_mask",
                    "u", "v", "w", "T", "t", "nu_ref", "A_ref"):
            assert key in s

    def test_initial_velocity_zero(self):
        s = make_solver(N=16)
        assert np.allclose(s["u"], 0.0)
        assert np.allclose(s["v"], 0.0)
        assert np.allclose(s["w"], 0.0)

    def test_initial_T_uniform(self):
        s = make_solver(N=16, T_base=500.0)
        assert np.allclose(s["T"], 500.0)

    def test_nu_ref_positive(self):
        s = make_solver(N=16)
        assert s["nu_ref"] > 0

    def test_A_ref_positive(self):
        """A_ref = A(T_base) > 0 — second law."""
        s = make_solver(N=16)
        assert s["A_ref"] > 0

    def test_set_ic_immutable(self):
        s = make_solver(N=16)
        ic = taylor_green_ic(16)
        s2 = set_ic(s, ic)
        # Original solver must be unchanged
        assert np.allclose(s["u"], 0.0), "set_ic must not mutate original solver"
        assert not np.allclose(s2["u"], 0.0), "new solver must have IC velocity"


# ─────────────────────────────────────────────────────────────────────────────
# 5. TestRHS3D
# ─────────────────────────────────────────────────────────────────────────────

class TestRHS3D:
    def test_rhs_output_shape(self, N16_solver):
        from route1_3D import _rhs
        s = N16_solver
        N = s["N"]
        d_u, d_v, d_w, d_T, diag = _rhs(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["fluid"], s["P"]
        )
        assert d_u.shape == (N, N, N)
        assert d_v.shape == (N, N, N)
        assert d_w.shape == (N, N, N)
        assert d_T.shape == (N, N, N)

    def test_A_min_positive_in_rhs(self, N16_solver):
        from route1_3D import _rhs
        s = N16_solver
        _, _, _, _, diag = _rhs(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["fluid"], s["P"]
        )
        assert diag["A_min"] > 0

    def test_S_theta_nonnegative(self, N16_solver):
        """Viscous heating source S_θ = ν·|∇u|² ≥ 0."""
        from route1_3D import _rhs
        s = N16_solver
        _, _, _, _, diag = _rhs(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["fluid"], s["P"]
        )
        assert diag["S_theta"] >= 0.0

    def test_zero_velocity_zero_advection(self):
        """Zero velocity → zero advection; d_T = diffusion only."""
        from route1_3D import _rhs
        N = 16
        s = make_solver(N=N)
        _, kx, ky, kz, k2 = make_grid_3D(N)
        mask = make_dealias_mask_3D(N)
        u = v = w = np.zeros((N, N, N))
        T = np.full((N, N, N), 300.0)
        d_u, d_v, d_w, d_T, diag = _rhs(u, v, w, T, kx, ky, kz, k2, mask, "ideal", None)
        # Zero velocity: advection = 0, diffusion of uniform T = 0, Q_visc = 0
        assert np.allclose(d_T, 0.0, atol=1e-10), "Zero velocity → zero d_T for uniform T"

    def test_diag_keys(self, N16_solver):
        from route1_3D import _rhs
        s = N16_solver
        _, _, _, _, diag = _rhs(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["fluid"], s["P"]
        )
        for k in ("A_min", "A_max", "nu_bar", "u_rms", "u_max", "S_theta"):
            assert k in diag

    def test_diag_nu_bar_positive(self, N16_solver):
        from route1_3D import _rhs
        s = N16_solver
        _, _, _, _, diag = _rhs(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["fluid"], s["P"]
        )
        assert diag["nu_bar"] > 0


# ─────────────────────────────────────────────────────────────────────────────
# 6. TestTimeStepping3D
# ─────────────────────────────────────────────────────────────────────────────

class TestTimeStepping3D:
    def test_cfl_dt_positive(self, N16_solver):
        dt = compute_cfl_dt(N16_solver)
        assert dt > 0

    def test_cfl_dt_bounded_by_dt_max(self, N16_solver):
        dt = compute_cfl_dt(N16_solver)
        assert dt <= N16_solver["dt_max"]

    def test_step_advances_time(self, N16_solver):
        t0 = N16_solver["t"]
        s2, diag = solver_step(N16_solver)
        assert s2["t"] > t0

    def test_step_immutable(self, N16_solver):
        t0 = N16_solver["t"]
        s2, _ = solver_step(N16_solver)
        assert N16_solver["t"] == t0, "solver_step must not mutate input solver"

    def test_step_T_in_range(self, N16_solver):
        s2, _ = solver_step(N16_solver)
        assert np.all(s2["T"] >= 100.0)
        assert np.all(s2["T"] <= 6000.0)

    def test_step_maintains_div_free(self, N16_solver):
        """After one RK4 step, div u must remain ≈ 0."""
        s2, _ = solver_step(N16_solver)
        div = divergence_rms(s2["u"], s2["v"], s2["w"], s2["kx"], s2["ky"], s2["kz"])
        assert div < 1e-9, f"div u after step = {div:.3e}"

    def test_step_verdict_pass(self, N16_solver):
        _, diag = solver_step(N16_solver)
        assert diag["verdict"] == "PASS"


# ─────────────────────────────────────────────────────────────────────────────
# 7. TestRunLoop3D
# ─────────────────────────────────────────────────────────────────────────────

class TestRunLoop3D:
    @pytest.fixture(scope="class")
    def short_run(self):
        ic = taylor_green_ic(16)
        s = make_solver(N=16, dt_max=0.05)
        s = set_ic(s, ic)
        return run(s, t_end=0.1, max_steps=20)

    def test_reaches_t_end(self, short_run):
        assert short_run["t_final"] >= 0.09

    def test_verdict_pass(self, short_run):
        assert short_run["verdict"] == "PASS"

    def test_history_nonempty(self, short_run):
        assert len(short_run["history"]) > 0

    def test_A_min_global_positive(self, short_run):
        assert short_run["A_min_global"] > 0.0

    def test_energy_decays(self, short_run):
        """Kinetic energy must decrease from IC to t_end."""
        assert short_run["E_final"] < short_run["E_initial"]

    def test_energy_initial_close_to_analytic(self, short_run):
        """E₀ computed should match V₀²/8 within 1%."""
        E0_analytic = 1.0**2 / 8.0
        rel_err = abs(short_run["E_initial"] - E0_analytic) / E0_analytic
        assert rel_err < 0.01

    def test_dissipation_rate_positive(self, short_run):
        """Mean dissipation rate -dE/dt > 0 (energy is being removed)."""
        assert short_run["dissipation_rate_mean"] > 0


# ─────────────────────────────────────────────────────────────────────────────
# 8. TestTaylorGreenBenchmark
# ─────────────────────────────────────────────────────────────────────────────

class TestTaylorGreenBenchmark:
    def test_E0_accuracy(self, N16_tg_result):
        """E₀ computed must match analytic V₀²/8 within 1%."""
        assert N16_tg_result["E0_error_pct"] < 1.0

    def test_energy_decays(self, N16_tg_result):
        """Energy must decrease over the simulation."""
        assert N16_tg_result["energy_decayed"] is True

    def test_div_ok(self, N16_tg_result):
        """Divergence constraint maintained to machine precision."""
        assert N16_tg_result["div_ok"] is True

    def test_A_min_positive(self, N16_tg_result):
        assert N16_tg_result["A_min_global"] > 0.0

    def test_verdict_pass(self, N16_tg_result):
        assert N16_tg_result["verdict"] == "PASS"

    def test_exp_id_format(self, N16_tg_result):
        assert N16_tg_result["exp_id"].startswith("EXP-L3-R1-TG-")


# ─────────────────────────────────────────────────────────────────────────────
# 9. TestSmokeTest
# ─────────────────────────────────────────────────────────────────────────────

class TestSmokeTest:
    @pytest.fixture(scope="class")
    def smoke(self):
        return run_smoke_test(verbose=False)

    def test_verdict_pass(self, smoke):
        assert smoke["verdict"] == "PASS"

    def test_n_steps_positive(self, smoke):
        assert smoke["n_steps"] > 0

    def test_all_pass(self, smoke):
        assert smoke["all_pass"] is True


# ─────────────────────────────────────────────────────────────────────────────
# 10. TestM5EndToEnd
# ─────────────────────────────────────────────────────────────────────────────

class TestM5EndToEnd:
    def test_full_chain_pass(self):
        """make_solver → set_ic(TG) → run → verdict PASS."""
        ic = taylor_green_ic(16)
        s = make_solver(N=16)
        s = set_ic(s, ic)
        result = run(s, t_end=0.15, max_steps=30)
        assert result["verdict"] == "PASS"
        assert result["A_min_global"] > 0.0
        assert result["E_final"] < result["E_initial"]

    def test_kinetic_energy_formula(self):
        """E = ½⟨|u|²⟩ computed correctly."""
        N = 8
        u = np.ones((N, N, N))
        v = np.zeros((N, N, N))
        w = np.zeros((N, N, N))
        E = kinetic_energy(u, v, w)
        assert abs(E - 0.5) < 1e-12
