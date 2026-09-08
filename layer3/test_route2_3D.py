"""
test_route2_3D.py — Test suite for Route 2 3D solver (Exact Prize NS + θ)
=========================================================================
Tests for layer3/route2_3D.py

Physics validated:
  - Leray projector (Nyquist-zeroed, B-005 fix)
  - Source S_θ = ν|∇u|² ≥ 0 always
  - θ(x,0) = 0 always (Route 2 initial condition)
  - θ ≥ 0 after steps (source non-negative; CN scheme)
  - Energy decay (viscous dissipation)
  - div u ≈ 0 maintained by Leray re-projection
  - δ_est ≥ 0
  - LPS margin ε_LPS = min(μ_eff) − ν ≥ 0 by construction
  - Taylor-Green E₀ matches analytic formula within 1%
  - Sobolev norms monotone in s for smooth fields
  - CZ integrability probe finite

Run:
  cd /Users/pratikmenon/Documents/Claude/tfirst
  python -m pytest layer3/test_route2_3D.py -v
"""

from __future__ import annotations

import numpy as np
import pytest
import sys
import os

# ── path setup ──────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

import route2_3D as r2


# =============================================================================
# Fixtures
# =============================================================================

N_SMALL = 16     # fast, for most unit tests
N_MED   = 32     # for integration / benchmark tests

NU_DEFAULT = 1.0e-3


@pytest.fixture(scope="module")
def grid_small():
    return r2.make_grid_3D(N_SMALL)


@pytest.fixture(scope="module")
def grid_med():
    return r2.make_grid_3D(N_MED)


@pytest.fixture(scope="module")
def solver_small():
    s = r2.make_solver(N=N_SMALL, nu=NU_DEFAULT, eps_param=0.1)
    ic = r2.taylor_green_ic(N_SMALL, V0=1.0)
    return r2.set_ic(s, ic)


@pytest.fixture(scope="module")
def solver_med():
    s = r2.make_solver(N=N_MED, nu=NU_DEFAULT, eps_param=0.1)
    ic = r2.taylor_green_ic(N_MED, V0=1.0)
    return r2.set_ic(s, ic)


# =============================================================================
# 1. Grid and mask tests
# =============================================================================

class TestGrid:
    def test_grid_shape(self, grid_small):
        kx, ky, kz, k2 = grid_small
        for arr in (kx, ky, kz, k2):
            assert arr.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_k2_nonnegative(self, grid_small):
        _, _, _, k2 = grid_small
        assert np.all(k2 >= 0.0)

    def test_k2_zero_at_origin(self, grid_small):
        _, _, _, k2 = grid_small
        assert k2[0, 0, 0] == 0.0

    def test_dealias_mask_shape(self):
        mask = r2.make_dealias_mask_3D(N_SMALL)
        assert mask.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_dealias_mask_binary(self):
        mask = r2.make_dealias_mask_3D(N_SMALL)
        unique_vals = np.unique(mask)
        assert set(unique_vals).issubset({0.0, 1.0})

    def test_dealias_dc_mode_kept(self):
        mask = r2.make_dealias_mask_3D(N_SMALL)
        assert mask[0, 0, 0] == 1.0

    def test_dealias_high_k_zeroed(self):
        mask = r2.make_dealias_mask_3D(N_SMALL)
        # Mode at Nyquist should be zeroed
        nq = N_SMALL // 2
        assert mask[nq, 0, 0] == 0.0

    def test_grid_symmetry(self, grid_small):
        """kx[i,j,k] and kx[-i,j,k] should be negatives (FFT frequency symmetry)."""
        kx, ky, kz, _ = grid_small
        # Skip DC and Nyquist (index 0 and N//2)
        assert kx[1, 0, 0] == -kx[-1, 0, 0]
        assert ky[0, 1, 0] == -ky[0, -1, 0]
        assert kz[0, 0, 1] == -kz[0, 0, -1]


# =============================================================================
# 2. Leray projector tests (B-005 fixed)
# =============================================================================

class TestLeray:
    def test_tg_is_div_free(self, solver_small):
        """TG IC after set_ic should be divergence-free."""
        u, v, w = solver_small["u"], solver_small["v"], solver_small["w"]
        kx, ky, kz = solver_small["kx"], solver_small["ky"], solver_small["kz"]
        div = r2.divergence_rms(u, v, w, kx, ky, kz)
        assert div < 1e-12, f"TG IC div_rms={div:.2e} (expected < 1e-12)"

    def test_random_field_projected_div_free(self, grid_small):
        """Leray projector should make any random field divergence-free."""
        kx, ky, kz, k2 = grid_small
        rng = np.random.default_rng(123)
        u = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        v = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        w = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))

        up, vp, wp = r2.project_divergence_free(u, v, w, kx, ky, kz, k2)
        div = r2.divergence_rms(up, vp, wp, kx, ky, kz)
        assert div < 1e-10, f"div_rms after Leray={div:.2e}"

    def test_leray_idempotent(self, grid_small):
        """Projecting an already div-free field should change it negligibly."""
        kx, ky, kz, k2 = grid_small
        rng = np.random.default_rng(456)
        u = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        v = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        w = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))

        up, vp, wp = r2.project_divergence_free(u, v, w, kx, ky, kz, k2)
        up2, vp2, wp2 = r2.project_divergence_free(up, vp, wp, kx, ky, kz, k2)

        diff = np.sqrt(np.mean((up - up2)**2 + (vp - vp2)**2 + (wp - wp2)**2))
        assert diff < 1e-10, f"Leray not idempotent: Δu_rms={diff:.2e}"

    def test_leray_preserves_div_free(self, grid_small):
        """A purely solenoidal field should be unchanged by Leray projection."""
        # curl of a vector potential → div-free by construction
        kx, ky, kz, k2 = grid_small
        rng = np.random.default_rng(789)
        # Make a smooth curl field
        A_hat = np.zeros((N_SMALL, N_SMALL, N_SMALL), dtype=complex)
        k_abs = np.sqrt(k2)
        low_k = (k_abs > 0) & (k_abs <= 3)
        A_hat[low_k] = rng.standard_normal(int(np.sum(low_k)))
        A = np.real(np.fft.ifftn(A_hat))
        # u = curl(A, A, A) — simplified: use known div-free field
        x = np.linspace(0, 2*np.pi, N_SMALL, endpoint=False)
        xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
        u = np.sin(xx) * np.cos(yy)
        v = -np.cos(xx) * np.sin(yy)
        w = np.zeros_like(u)
        up, vp, wp = r2.project_divergence_free(u, v, w, kx, ky, kz, k2)
        diff = np.sqrt(np.mean((u - up)**2 + (v - vp)**2 + (w - wp)**2))
        assert diff < 1e-10, f"Leray changed div-free field: Δ={diff:.2e}"


# =============================================================================
# 3. Source S_θ tests
# =============================================================================

class TestSourceTheta:
    def test_source_nonnegative_tg(self, solver_small):
        """S_θ = ν|∇u|² must be ≥ 0 everywhere for TG IC."""
        u, v, w = solver_small["u"], solver_small["v"], solver_small["w"]
        kx, ky, kz = solver_small["kx"], solver_small["ky"], solver_small["kz"]
        nu = solver_small["nu"]
        S = r2._source_theta(u, v, w, kx, ky, kz, nu)
        assert np.all(S >= 0.0), f"S_θ has negative values: min={np.min(S):.2e}"

    def test_source_zero_for_zero_velocity(self, grid_small):
        """S_θ = 0 when u = v = w = 0."""
        kx, ky, kz, k2 = grid_small
        u = np.zeros((N_SMALL, N_SMALL, N_SMALL))
        v = np.zeros_like(u)
        w = np.zeros_like(u)
        S = r2._source_theta(u, v, w, kx, ky, kz, nu=1e-3)
        assert np.max(np.abs(S)) < 1e-12, f"S_θ non-zero for zero u: max={np.max(np.abs(S)):.2e}"

    def test_source_nonnegative_random(self, grid_small):
        """S_θ ≥ 0 for random divergence-free field."""
        kx, ky, kz, k2 = grid_small
        rng = np.random.default_rng(111)
        u = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        v = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        w = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        S = r2._source_theta(u, v, w, kx, ky, kz, nu=1e-3)
        # Allow tiny numerical negatives due to floating point
        assert np.min(S) >= -1e-12, f"S_θ significantly negative: min={np.min(S):.2e}"

    def test_source_scales_with_nu(self, grid_small):
        """S_θ = ν|∇u|² should scale linearly with ν."""
        kx, ky, kz, k2 = grid_small
        rng = np.random.default_rng(222)
        u = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        v = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        w = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        nu1, nu2 = 1e-3, 2e-3
        S1 = r2._source_theta(u, v, w, kx, ky, kz, nu1)
        S2 = r2._source_theta(u, v, w, kx, ky, kz, nu2)
        ratio = np.mean(S2) / np.mean(S1)
        assert abs(ratio - 2.0) < 1e-10, f"S_θ does not scale with ν: ratio={ratio:.6f}"


# =============================================================================
# 4. Initial condition tests
# =============================================================================

class TestInitialConditions:
    def test_theta_zero_at_t0(self, solver_small):
        """θ must be exactly 0 at t=0 (PRD §3.3)."""
        theta = solver_small["theta"]
        assert np.max(np.abs(theta)) == 0.0, f"θ not zero at t=0: max|θ|={np.max(np.abs(theta)):.2e}"

    def test_t_zero_at_start(self, solver_small):
        """Solver time must be 0.0 after set_ic."""
        assert solver_small["t"] == 0.0

    def test_tg_e0_analytic_small(self):
        """TG E₀ = V₀²/8 for small N."""
        ic = r2.taylor_green_ic(N_SMALL, V0=1.0)
        E0_analytic = ic["E0_analytic"]
        assert abs(E0_analytic - 0.125) < 1e-12, f"TG E₀ analytic wrong: {E0_analytic}"

    def test_tg_e0_computed_matches_analytic(self, solver_small):
        """Computed E₀ should match TG analytic E₀ = V₀²/8 within 1%."""
        u, v, w = solver_small["u"], solver_small["v"], solver_small["w"]
        E0_computed = r2.kinetic_energy(u, v, w)
        E0_analytic = 0.125  # V0=1
        err = abs(E0_computed - E0_analytic) / E0_analytic
        assert err < 0.01, f"TG E₀ error {err*100:.2f}% > 1%"

    def test_tg_w_component_zero(self, solver_small):
        """TG vortex has w=0 (2D structure in 3D)."""
        w = solver_small["w"]
        assert np.max(np.abs(w)) < 1e-14, f"TG w not zero: max|w|={np.max(np.abs(w)):.2e}"

    def test_shear_layer_ic_structure(self):
        """Shear layer IC: u varies in z, v=w=0."""
        ic = r2.shear_layer_ic(N_SMALL, amp=1.0, width=0.5)
        assert np.max(np.abs(ic["v"])) < 1e-14
        assert np.max(np.abs(ic["w"])) < 1e-14
        # u should vary with z
        assert np.std(ic["u"]) > 0.1

    def test_random_div_free_ic(self):
        """Random IC should be approximately div-free after set_ic."""
        kx, ky, kz, k2 = r2.make_grid_3D(N_SMALL)
        ic = r2.random_div_free_ic(N_SMALL, amp=0.5, seed=99)
        div = r2.divergence_rms(ic["u"], ic["v"], ic["w"], kx, ky, kz)
        assert div < 1e-10, f"Random IC not div-free: div_rms={div:.2e}"

    def test_set_ic_returns_new_dict(self, solver_small):
        """set_ic must not mutate the input solver (immutable pattern)."""
        solver_orig = r2.make_solver(N=N_SMALL)
        ic = r2.taylor_green_ic(N_SMALL)
        solver_new = r2.set_ic(solver_orig, ic)
        # Original should still have zero velocity
        assert solver_orig is not solver_new
        assert np.max(np.abs(solver_orig["u"])) == 0.0


# =============================================================================
# 5. mu_eff and LPS margin tests
# =============================================================================

class TestMuEff:
    def test_mu_eff_ge_nu(self):
        """μ_eff = ν + ε·f(θ) must be ≥ ν everywhere (ε ≥ 0, f(θ) ≥ 0)."""
        nu = 1e-3
        eps = 0.5
        theta = np.random.default_rng(0).uniform(0, 1, (N_SMALL, N_SMALL, N_SMALL))
        mu = r2.mu_eff_field(theta, nu, eps, "identity")
        assert np.all(mu >= nu), f"μ_eff < ν at some point: min={np.min(mu):.4e}"

    def test_mu_eff_eps_zero_equals_nu(self):
        """At ε=0, μ_eff = ν everywhere (exact Prize equations)."""
        nu = 1e-3
        theta = np.ones((N_SMALL, N_SMALL, N_SMALL)) * 0.5
        mu = r2.mu_eff_field(theta, nu, eps_param=0.0, f_theta="identity")
        assert np.allclose(mu, nu), f"μ_eff ≠ ν at ε=0: max diff={np.max(np.abs(mu - nu)):.2e}"

    def test_mu_eff_identity(self):
        """f_theta='identity': μ_eff = ν + ε·θ at θ=1."""
        nu, eps = 1e-3, 0.1
        theta = np.ones((4, 4, 4))
        mu = r2.mu_eff_field(theta, nu, eps, "identity")
        expected = nu + eps * 1.0
        assert np.allclose(mu, expected, rtol=1e-10)

    def test_mu_eff_tanh(self):
        """f_theta='tanh': μ_eff = ν + ε·tanh(θ)."""
        nu, eps = 1e-3, 0.1
        theta = np.full((4, 4, 4), 0.5)
        mu = r2.mu_eff_field(theta, nu, eps, "tanh")
        expected = nu + eps * np.tanh(0.5)
        assert np.allclose(mu, expected, rtol=1e-10)

    def test_lps_margin_nonneg_by_construction(self):
        """LPS margin min(μ_eff) − ν ≥ 0 always (ε ≥ 0, f(θ) ≥ 0)."""
        nu, eps = 1e-3, 0.1
        theta = np.random.default_rng(7).uniform(0, 2, (N_SMALL, N_SMALL, N_SMALL))
        mu = r2.mu_eff_field(theta, nu, eps, "identity")
        eps_lps = float(np.min(mu)) - nu
        assert eps_lps >= 0.0, f"LPS margin negative: {eps_lps:.4e}"

    def test_unknown_f_theta_raises(self):
        """Unknown f_theta string should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown f_theta"):
            r2.mu_eff_field(np.zeros((4, 4, 4)), nu=1e-3, eps_param=0.1, f_theta="cubic")


# =============================================================================
# 6. Diagnostics tests (δ, CZ, Sobolev norms)
# =============================================================================

class TestDiagnostics:
    def test_delta_zero_theta_returns_zero(self, grid_small):
        """δ_est should be 0 when θ = 0 (early time)."""
        kx, ky, kz, _ = grid_small
        theta = np.zeros((N_SMALL, N_SMALL, N_SMALL))
        delta = r2.delta_from_theta(theta, kx, ky, kz)
        assert delta == 0.0

    def test_delta_nonneg_for_smooth_theta(self, grid_small):
        """δ_est ≥ 0 for any non-zero θ."""
        kx, ky, kz, _ = grid_small
        rng = np.random.default_rng(42)
        theta = np.clip(rng.standard_normal((N_SMALL, N_SMALL, N_SMALL)), 0, None)
        delta = r2.delta_from_theta(theta, kx, ky, kz)
        assert delta >= 0.0, f"δ_est < 0: {delta}"

    def test_cz_integrability_finite(self, grid_small):
        """CZ integrability should return a finite positive number."""
        rng = np.random.default_rng(99)
        source = np.abs(rng.standard_normal((N_SMALL, N_SMALL, N_SMALL)))
        cz = r2.cz_integrability(source, eps_probe=0.1)
        assert np.isfinite(cz) and cz > 0.0

    def test_cz_zero_source_is_zero(self):
        """Zero source → CZ norm = 0."""
        source = np.zeros((N_SMALL, N_SMALL, N_SMALL))
        cz = r2.cz_integrability(source, eps_probe=0.1)
        assert cz == 0.0

    def test_sobolev_norms_ordered(self, grid_small):
        """H^s norms should be monotonically increasing in s for non-smooth θ."""
        kx, ky, kz, _ = grid_small
        # A theta with energy at high k — H^s should grow with s
        rng = np.random.default_rng(55)
        theta = rng.standard_normal((N_SMALL, N_SMALL, N_SMALL))
        Hs = r2.sobolev_norms(theta, kx, ky, kz, s_list=(0.0, 0.5, 1.0, 1.5, 2.0))
        vals = [Hs[f"H{s:.1f}"] for s in [0.0, 0.5, 1.0, 1.5, 2.0]]
        for i in range(len(vals) - 1):
            assert vals[i] <= vals[i+1] + 1e-8, \
                f"H^{0.5*i} > H^{0.5*(i+1)}: {vals[i]:.4e} > {vals[i+1]:.4e}"

    def test_sobolev_norms_zero_theta(self, grid_small):
        """H^s norms should all be 0 for θ=0."""
        kx, ky, kz, _ = grid_small
        theta = np.zeros((N_SMALL, N_SMALL, N_SMALL))
        Hs = r2.sobolev_norms(theta, kx, ky, kz)
        for k, v in Hs.items():
            assert v < 1e-14, f"{k} norm non-zero for zero θ: {v:.2e}"

    def test_divergence_rms_formula(self, grid_small):
        """div_rms for a known field with known divergence."""
        kx, ky, kz, k2 = grid_small
        x = np.linspace(0, 2*np.pi, N_SMALL, endpoint=False)
        xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
        # sin(x) → ∂u/∂x = cos(x), so divf ≠ 0 in general
        u = np.sin(xx)
        v = np.zeros_like(u)
        w = np.zeros_like(u)
        div = r2.divergence_rms(u, v, w, kx, ky, kz)
        # ∂u/∂x = cos(x), so div_rms ≈ 1/√2 ≈ 0.707
        assert 0.3 < div < 0.9, f"div_rms = {div:.4f} (expected ~0.707)"

    def test_kinetic_energy_formula(self):
        """E = ½⟨|u|²⟩ should match manual computation."""
        u = np.ones((4, 4, 4)) * 2.0
        v = np.ones((4, 4, 4)) * 2.0
        w = np.zeros((4, 4, 4))
        E = r2.kinetic_energy(u, v, w)
        E_manual = 0.5 * (4.0 + 4.0)  # ½(2² + 2²)
        assert abs(E - E_manual) < 1e-12


# =============================================================================
# 7. Solver construction tests
# =============================================================================

class TestMakeSolver:
    def test_make_solver_default_fields_zero(self):
        """make_solver should initialise u,v,w,theta to zero."""
        s = r2.make_solver(N=N_SMALL)
        for field in ("u", "v", "w", "theta"):
            assert np.max(np.abs(s[field])) == 0.0

    def test_make_solver_nu_positive_required(self):
        """nu ≤ 0 must raise ValueError."""
        with pytest.raises(ValueError, match="nu must be > 0"):
            r2.make_solver(N=N_SMALL, nu=0.0)
        with pytest.raises(ValueError, match="nu must be > 0"):
            r2.make_solver(N=N_SMALL, nu=-1.0)

    def test_make_solver_eps_negative_raises(self):
        """eps_param < 0 must raise ValueError."""
        with pytest.raises(ValueError, match="eps_param must be ≥ 0"):
            r2.make_solver(N=N_SMALL, eps_param=-0.1)

    def test_make_solver_eps_zero_allowed(self):
        """eps_param = 0.0 must be allowed (exact Prize equations)."""
        s = r2.make_solver(N=N_SMALL, eps_param=0.0)
        assert s["eps_param"] == 0.0

    def test_make_solver_shapes(self):
        """All grid and field arrays should have shape (N, N, N)."""
        s = r2.make_solver(N=N_SMALL)
        for key in ("kx", "ky", "kz", "k2", "dealias", "u", "v", "w", "theta"):
            assert s[key].shape == (N_SMALL, N_SMALL, N_SMALL), \
                f"{key} shape {s[key].shape} ≠ ({N_SMALL},{N_SMALL},{N_SMALL})"

    def test_make_solver_t_zero(self):
        """Solver time should start at 0.0."""
        s = r2.make_solver(N=N_SMALL)
        assert s["t"] == 0.0


# =============================================================================
# 8. CFL timestep tests
# =============================================================================

class TestCFL:
    def test_cfl_positive(self, solver_small):
        """CFL timestep must be positive."""
        dt = r2.compute_cfl_dt(solver_small)
        assert dt > 0.0

    def test_cfl_bounded_by_dt_max(self, solver_small):
        """CFL dt must not exceed dt_max."""
        dt = r2.compute_cfl_dt(solver_small)
        assert dt <= solver_small["dt_max"]

    def test_cfl_zero_velocity_gives_dt_max(self):
        """Zero velocity field → no advection CFL → dt limited by dt_max and diffusion."""
        s = r2.make_solver(N=N_SMALL, nu=NU_DEFAULT, dt_max=0.05)
        # u=v=w=0 → U_max=1e-12 → dt_adv very large → capped at dt_max
        dt = r2.compute_cfl_dt(s)
        assert dt <= 0.05

    def test_cfl_scales_with_velocity(self, grid_small):
        """Higher velocity → smaller CFL dt."""
        s1 = r2.make_solver(N=N_SMALL, nu=NU_DEFAULT, dt_max=1.0)
        ic1 = r2.taylor_green_ic(N_SMALL, V0=1.0)
        s1 = r2.set_ic(s1, ic1)

        s2 = r2.make_solver(N=N_SMALL, nu=NU_DEFAULT, dt_max=1.0)
        ic2 = r2.taylor_green_ic(N_SMALL, V0=5.0)
        s2 = r2.set_ic(s2, ic2)

        dt1 = r2.compute_cfl_dt(s1)
        dt2 = r2.compute_cfl_dt(s2)
        assert dt2 < dt1, f"Higher V₀ should give smaller dt: dt1={dt1:.4e} dt2={dt2:.4e}"


# =============================================================================
# 9. Single step tests
# =============================================================================

class TestSolverStep:
    def test_step_returns_two_items(self, solver_small):
        """solver_step must return (new_solver, diagnostics)."""
        new_s, diag = r2.solver_step(solver_small, dt=1e-4)
        assert isinstance(new_s, dict)
        assert isinstance(diag, dict)

    def test_step_immutable(self, solver_small):
        """solver_step must not mutate the input solver."""
        u_orig = solver_small["u"].copy()
        t_orig = solver_small["t"]
        _ = r2.solver_step(solver_small, dt=1e-4)
        np.testing.assert_array_equal(solver_small["u"], u_orig)
        assert solver_small["t"] == t_orig

    def test_step_advances_time(self, solver_small):
        """Time must advance by dt after one step."""
        dt = 1e-4
        t0 = solver_small["t"]
        new_s, _ = r2.solver_step(solver_small, dt=dt)
        assert abs(new_s["t"] - (t0 + dt)) < 1e-12

    def test_step_theta_nonneg(self, solver_small):
        """θ must remain ≥ 0 after each step (source ≥ 0 + clipping)."""
        new_s, _ = r2.solver_step(solver_small, dt=1e-4)
        assert np.min(new_s["theta"]) >= 0.0, \
            f"θ negative after step: min={np.min(new_s['theta']):.2e}"

    def test_step_div_free_maintained(self, solver_small):
        """div u should remain < 1e-10 after one step."""
        new_s, _ = r2.solver_step(solver_small, dt=1e-4)
        kx, ky, kz = new_s["kx"], new_s["ky"], new_s["kz"]
        div = r2.divergence_rms(new_s["u"], new_s["v"], new_s["w"], kx, ky, kz)
        assert div < 1e-10, f"div_rms after step = {div:.2e} (expected < 1e-10)"

    def test_step_diag_has_required_keys(self, solver_small):
        """Step diagnostics must include all required keys."""
        _, diag = r2.solver_step(solver_small, dt=1e-4)
        required = {
            "t", "dt", "E", "div_rms",
            "theta_min", "theta_max", "theta_mean",
            "S_theta_max", "S_theta_mean",
            "mu_eff_min", "eps_lps", "delta_est", "cz_L1p",
        }
        for key in required:
            assert key in diag, f"Missing diagnostic key: '{key}'"

    def test_step_eps_lps_nonneg(self, solver_small):
        """LPS margin ε_LPS = min(μ_eff) − ν must be ≥ 0."""
        _, diag = r2.solver_step(solver_small, dt=1e-4)
        assert diag["eps_lps"] >= 0.0, f"ε_LPS = {diag['eps_lps']:.4e} < 0"

    def test_step_source_max_nonneg(self, solver_small):
        """S_θ_max must be ≥ 0."""
        _, diag = r2.solver_step(solver_small, dt=1e-4)
        assert diag["S_theta_max"] >= 0.0

    def test_step_delta_est_nonneg(self, solver_small):
        """δ_est must be ≥ 0."""
        _, diag = r2.solver_step(solver_small, dt=1e-4)
        assert diag["delta_est"] >= 0.0

    def test_step_auto_dt(self, solver_small):
        """solver_step with dt=None should use CFL dt automatically."""
        new_s, diag = r2.solver_step(solver_small, dt=None)
        assert diag["dt"] > 0.0
        assert new_s["t"] > solver_small["t"]


# =============================================================================
# 10. Multi-step integration tests
# =============================================================================

class TestIntegration:
    def test_energy_decays_over_time(self, solver_small):
        """Energy should decrease monotonically over ~20 steps (viscous dissipation)."""
        solver = solver_small
        E_prev = r2.kinetic_energy(solver["u"], solver["v"], solver["w"])
        for _ in range(20):
            solver, diag = r2.solver_step(solver, dt=1e-4)

        E_final = r2.kinetic_energy(solver["u"], solver["v"], solver["w"])
        assert E_final < E_prev, \
            f"Energy did not decay: E_initial={E_prev:.4e}, E_final={E_final:.4e}"

    def test_theta_grows_from_zero(self, solver_small):
        """θ should grow from 0 as S_θ = ν|∇u|² accumulates."""
        solver = solver_small
        theta_max_prev = 0.0
        for _ in range(10):
            solver, diag = r2.solver_step(solver, dt=1e-4)

        theta_max = float(np.max(solver["theta"]))
        assert theta_max > theta_max_prev, \
            f"θ did not grow from 0: theta_max={theta_max:.2e}"

    def test_theta_stays_nonneg_multi_step(self, solver_small):
        """θ ≥ 0 must hold throughout multi-step integration."""
        solver = solver_small
        for step in range(30):
            solver, diag = r2.solver_step(solver, dt=1e-4)
            assert np.min(solver["theta"]) >= 0.0, \
                f"θ < 0 at step {step}: min={np.min(solver['theta']):.2e}"

    def test_div_free_multi_step(self, solver_small):
        """div u < 1e-9 after 20 steps."""
        solver = solver_small
        for _ in range(20):
            solver, _ = r2.solver_step(solver, dt=1e-4)
        kx, ky, kz = solver["kx"], solver["ky"], solver["kz"]
        div = r2.divergence_rms(solver["u"], solver["v"], solver["w"], kx, ky, kz)
        assert div < 1e-9, f"div_rms after 20 steps = {div:.2e}"

    def test_eps_zero_gives_valid_run(self, grid_small):
        """ε=0 (exact Prize equations) should run without error."""
        s = r2.make_solver(N=N_SMALL, nu=NU_DEFAULT, eps_param=0.0)
        ic = r2.taylor_green_ic(N_SMALL, V0=0.5)
        s = r2.set_ic(s, ic)
        for _ in range(5):
            s, diag = r2.solver_step(s, dt=1e-4)
        assert diag["eps_lps"] == 0.0  # ε=0 → μ_eff = ν → margin = 0

    def test_run_function_verdict(self, solver_small):
        """run() should return a dict with 'verdict', 'history', and pass for TG."""
        result = r2.run(solver_small, t_end=0.002, max_steps=30, verbose=False)
        assert "verdict" in result
        assert "history" in result
        assert len(result["history"]) > 0
        assert result["verdict"] == "PASS"

    def test_run_energy_decayed_flag(self, solver_small):
        """energy_decayed flag should be True for TG (ν > 0)."""
        result = r2.run(solver_small, t_end=0.005, max_steps=100, verbose=False)
        assert result["energy_decayed"] is True

    def test_run_theta_negative_flag_false_for_tg(self, solver_small):
        """theta_negative flag should be False for TG IC."""
        result = r2.run(solver_small, t_end=0.002, max_steps=30, verbose=False)
        assert result["theta_negative"] is False

    def test_run_delta_max_nonneg(self, solver_small):
        """delta_max from run() must be ≥ 0."""
        result = r2.run(solver_small, t_end=0.002, max_steps=30, verbose=False)
        assert result["delta_max"] >= 0.0

    def test_run_eps_lps_min_nonneg(self, solver_small):
        """eps_lps_min from run() must be ≥ 0."""
        result = r2.run(solver_small, t_end=0.002, max_steps=30, verbose=False)
        assert result["eps_lps_min"] >= 0.0


# =============================================================================
# 11. Nyquist mask (B-005 fix) test
# =============================================================================

class TestNyquistZero:
    def test_zero_nyquist_kills_nyquist_planes(self):
        """_zero_nyquist should zero the three Nyquist planes."""
        N = 8
        f = np.ones((N, N, N), dtype=complex)
        f_zeroed = r2._zero_nyquist(f, N)
        nq = N // 2
        assert np.all(f_zeroed[nq, :, :] == 0.0)
        assert np.all(f_zeroed[:, nq, :] == 0.0)
        assert np.all(f_zeroed[:, :, nq] == 0.0)

    def test_zero_nyquist_preserves_other_modes(self):
        """_zero_nyquist should not touch non-Nyquist modes."""
        N = 8
        f = np.ones((N, N, N), dtype=complex) * 3.0
        f_zeroed = r2._zero_nyquist(f, N)
        # Mode [1, 1, 1] should be unchanged
        assert f_zeroed[1, 1, 1] == 3.0

    def test_zero_nyquist_does_not_mutate_input(self):
        """_zero_nyquist must return a copy, not mutate the input."""
        N = 8
        f = np.ones((N, N, N), dtype=complex)
        f_orig = f.copy()
        _ = r2._zero_nyquist(f, N)
        np.testing.assert_array_equal(f, f_orig)


# =============================================================================
# 12. Benchmark: Taylor-Green N=32 (EXP-L3-R2-TG-032)
# =============================================================================

class TestTaylorGreenBenchmark32:
    """Integration test: EXP-L3-R2-TG-032.

    Pass criteria (PRD §7.3):
      1. E₀ error < 1%
      2. θ ≥ 0 always
      3. Energy decays
      4. ε_LPS ≥ 0
    """

    @pytest.fixture(scope="class")
    def tg32_result(self):
        return r2.run_taylor_green(
            N=N_MED,
            t_end=0.5,
            nu=NU_DEFAULT,
            eps_param=0.1,
            max_steps=2000,
            verbose=False,
            db_path=None,   # no DB write in unit tests
        )

    def test_tg32_verdict_pass(self, tg32_result):
        assert tg32_result["verdict"] == "PASS", \
            f"EXP-L3-R2-TG-032 FAIL: {tg32_result.get('key_metric', '')}"

    def test_tg32_e0_error_lt_1pct(self, tg32_result):
        err = tg32_result["E0_error_pct"]
        assert err < 1.0, f"E₀ error = {err:.3f}% > 1%"

    def test_tg32_theta_nonneg(self, tg32_result):
        assert not tg32_result["theta_negative"], "θ went negative in TG32 run"

    def test_tg32_energy_decayed(self, tg32_result):
        assert tg32_result["energy_decayed"], "Energy did not decay in TG32 run"

    def test_tg32_eps_lps_min_nonneg(self, tg32_result):
        assert tg32_result["eps_lps_min"] >= 0.0, \
            f"ε_LPS_min = {tg32_result['eps_lps_min']:.4e} < 0"

    def test_tg32_n_steps_gt_zero(self, tg32_result):
        assert tg32_result["n_steps"] > 0

    def test_tg32_t_final_ge_target(self, tg32_result):
        assert tg32_result["t_final"] >= 0.45, \
            f"t_final={tg32_result['t_final']:.4f} < 0.45 (too early termination)"


# =============================================================================
# 13. f_theta variants
# =============================================================================

class TestFTheta:
    """Test that all f_theta variants satisfy f(θ) ≥ 0 and f(0) = 0."""

    VARIANTS = ["identity", "tanh", "sqrt"]

    def test_f_identity_at_zero(self):
        theta = np.zeros((4, 4, 4))
        mu = r2.mu_eff_field(theta, nu=1e-3, eps_param=0.5, f_theta="identity")
        assert np.allclose(mu, 1e-3)  # f(0)=0 → μ_eff=ν

    def test_f_tanh_at_zero(self):
        theta = np.zeros((4, 4, 4))
        mu = r2.mu_eff_field(theta, nu=1e-3, eps_param=0.5, f_theta="tanh")
        assert np.allclose(mu, 1e-3)  # tanh(0)=0

    def test_f_sqrt_at_zero(self):
        theta = np.zeros((4, 4, 4))
        mu = r2.mu_eff_field(theta, nu=1e-3, eps_param=0.5, f_theta="sqrt")
        assert np.allclose(mu, 1e-3, atol=1e-10)  # sqrt(ε_reg) - sqrt(ε_reg) ≈ 0

    @pytest.mark.parametrize("f_th", VARIANTS)
    def test_f_variant_mu_eff_ge_nu(self, f_th):
        nu, eps = 1e-3, 0.2
        theta = np.abs(np.random.default_rng(1).standard_normal((4, 4, 4)))
        mu = r2.mu_eff_field(theta, nu, eps, f_th)
        assert np.all(mu >= nu - 1e-12), f"μ_eff < ν for f_theta='{f_th}'"


# =============================================================================
# 14. Epsilon sweep correctness
# =============================================================================

class TestEpsSweep:
    def test_eps_zero_theta_max_matches_higher_eps(self):
        """θ_max should be identical for ε=0 and ε=0.1 (source is ε-independent)."""
        results = []
        for eps in [0.0, 0.1]:
            s = r2.make_solver(N=N_SMALL, nu=NU_DEFAULT, eps_param=eps)
            ic = r2.taylor_green_ic(N_SMALL, V0=0.5)
            s = r2.set_ic(s, ic)
            for _ in range(10):
                s, _ = r2.solver_step(s, dt=1e-4)
            results.append(float(np.max(s["theta"])))

        # θ evolution depends only on u and S_θ, not on ε
        assert abs(results[0] - results[1]) < 1e-8, \
            f"θ differs between ε=0 and ε=0.1: {results[0]:.6e} vs {results[1]:.6e}"


# =============================================================================
# 15. Smoke test: run_taylor_green with small N, db_path=None
# =============================================================================

class TestSmokeRun:
    def test_run_taylor_green_smoke(self):
        """Smoke test: run_taylor_green N=16 for t=0.05, no DB."""
        r = r2.run_taylor_green(
            N=N_SMALL,
            t_end=0.05,
            nu=NU_DEFAULT,
            eps_param=0.1,
            max_steps=200,
            verbose=False,
            db_path=None,
        )
        assert r["verdict"] == "PASS"
        assert r["n_steps"] > 0
        assert not r["theta_negative"]

    def test_run_shear_layer_smoke(self):
        """Smoke test: run_shear_layer N=16 for t=0.05, no DB."""
        r = r2.run_shear_layer(
            N=N_SMALL,
            t_end=0.05,
            nu=NU_DEFAULT,
            eps_param=0.1,
            max_steps=200,
            verbose=False,
            db_path=None,
        )
        assert r["verdict"] == "PASS"
        assert not r["theta_negative"]

    def test_exp_id_format_tg(self):
        """exp_id should follow EXP-L3-R2-TG-{N:03d} format."""
        r = r2.run_taylor_green(N=16, t_end=0.01, max_steps=20, verbose=False)
        assert r["exp_id"] == "EXP-L3-R2-TG-016"

    def test_exp_id_format_shear(self):
        """exp_id should follow EXP-L3-R2-SH-{N:03d} format."""
        r = r2.run_shear_layer(N=16, t_end=0.01, max_steps=20, verbose=False)
        assert r["exp_id"] == "EXP-L3-R2-SH-016"
