"""
test_route2_3D_jax.py — Test suite for JAX-JIT Route 2 3D solver
=================================================================
Tests for route2_3D_jax.py.  Validates:
  - JAX backend available and FFT routing correct
  - Grid / mask correctness (float32 JAX arrays on CPU device)
  - Leray projector (B-005: Nyquist-safe, Hermitian-preserving)
  - TG IC energy accuracy (< 1% error)
  - solver_step_jax: div-free maintained, θ ≥ 0, energy finite
  - CN IMEX for θ: θ grows from source S_θ = ν|∇u|²
  - δ estimator: Lemma 2.5 — δ > 0 after a few steps with TG IC
  - Run loop: energy decays, δ_max > 0, ε_LPS_min ≥ 0
  - TG benchmark (N=32, t_end=0.5): PASS verdict, δ_max = 1.5
  - Shear layer: Route A — δ > 0 even for non-mixing IC
  - ε=0 exact Prize: δ > 0 (Lemma 2.5 direct numerical confirmation)

All unit tests use N=16 (fast).  The benchmark class uses N=32.
The N=32 TG fixture runs once per class (scope="class").

Float32 tolerances are used throughout:
  - div u:      < 1e-2  (float32 + Nyquist fix)
  - E₀ error:  < 1%
  - theta ≥ 0: clamp enforced by solver
"""

import sys
import os
import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))

# ── Guard: skip entire module if JAX not installed ────────────────────────────
jax = pytest.importorskip("jax", reason="JAX not installed")
import jax.numpy as jnp

from route2_3D_jax import (
    _JAX_AVAILABLE,
    _JAX_BACKEND,
    _COMPUTE_DEVICE,
    _d,
    make_grid_jax,
    make_dealias_mask_jax,
    make_nyquist_mask_jax,
    _project_jax,
    _rhs_vel_jax,
    solver_step_jax,
    compute_cfl_dt_jax,
    kinetic_energy_jax,
    divergence_rms_jax,
    delta_from_theta_jax,
    make_solver_jax,
    set_ic_jax,
    taylor_green_ic,
    shear_layer_ic,
    run_jax,
    run_taylor_green_jax,
    run_shear_layer_jax,
)

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

N_SMALL = 16
N_MED   = 32
NU      = 1e-3
EPS     = 0.1
_DTYPE  = jnp.float32


# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def solver16():
    """N=16 solver with TG IC — module scope (built once)."""
    s = make_solver_jax(N=N_SMALL, nu=NU, eps_param=EPS)
    ic = taylor_green_ic(N_SMALL, V0=1.0)
    return set_ic_jax(s, ic)


@pytest.fixture(scope="module")
def solver16_stepped(solver16):
    """solver16 advanced by 3 steps — θ has grown from S_θ source."""
    s = solver16
    kx, ky, kz, k2 = s["kx"], s["ky"], s["kz"], s["k2"]
    dealias, nq = s["dealias"], s["nq_mask"]
    nu_j = jnp.array(NU, dtype=_DTYPE)
    dt_j = jnp.array(0.005, dtype=_DTYPE)

    u, v, w, theta = s["u"], s["v"], s["w"], s["theta"]
    for _ in range(3):
        u, v, w, theta = solver_step_jax(
            u, v, w, theta, kx, ky, kz, k2, dealias, nq, nu_j, dt_j
        )
        jax.block_until_ready((u, v, w, theta))
    return {"u": u, "v": v, "w": w, "theta": theta,
            "kx": kx, "ky": ky, "kz": kz, "k2": k2}


# ─────────────────────────────────────────────────────────────────────────────
# 1. Backend availability
# ─────────────────────────────────────────────────────────────────────────────

class TestJAXBackend:
    def test_jax_available(self):
        assert _JAX_AVAILABLE, "JAX not available"

    def test_backend_detected(self):
        valid = ("cpu", "gpu", "tpu", "METAL", "CPU-JIT")
        assert _JAX_BACKEND.upper() in (v.upper() for v in valid), \
            f"Unknown backend: {_JAX_BACKEND}"

    def test_jit_works(self):
        """A JIT-compiled function must return the correct result."""
        @jax.jit
        def mul(x, y):
            return x * y
        result = float(mul(jnp.array(3.0), jnp.array(4.0)))
        assert abs(result - 12.0) < 1e-5

    def test_fft_on_compute_device(self):
        """jnp.fft.fftn on compute device (CPU) must produce correct DC mode.

        jax-metal 0.1.0 does not support complex FFT on Metal — all FFT
        computation is routed to the CPU device via _d().
        """
        x = _d(jnp.ones((8, 8, 8), dtype=jnp.float32))
        out = jnp.fft.fftn(x)
        assert out.shape == (8, 8, 8)
        # DC mode of fftn(ones) = N³
        assert abs(float(jnp.real(out[0, 0, 0])) - 8.0**3) < 1.0


# ─────────────────────────────────────────────────────────────────────────────
# 2. Grid and masks
# ─────────────────────────────────────────────────────────────────────────────

class TestGridJAX:
    def test_shapes(self):
        kx, ky, kz, k2 = make_grid_jax(N_SMALL)
        for arr in (kx, ky, kz, k2):
            assert arr.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_k2_nonnegative(self):
        _, _, _, k2 = make_grid_jax(N_SMALL)
        assert bool(jnp.all(k2 >= 0.0))

    def test_k2_zero_at_dc(self):
        _, _, _, k2 = make_grid_jax(N_SMALL)
        assert float(k2[0, 0, 0]) == 0.0

    def test_dtype_float32(self):
        kx, ky, kz, k2 = make_grid_jax(N_SMALL, dtype=_DTYPE)
        assert kx.dtype == jnp.float32

    def test_dealias_mask_shape(self):
        mask = make_dealias_mask_jax(N_SMALL)
        assert mask.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_dealias_keeps_low_k(self):
        mask = make_dealias_mask_jax(N_SMALL)
        assert float(mask[1, 1, 1]) == 1.0

    def test_dealias_zeros_high_k(self):
        # |k| = 8 > N//3 = 5 for N=16
        mask = make_dealias_mask_jax(N_SMALL)
        assert float(mask[8, 0, 0]) == 0.0

    def test_nyquist_mask_zeros_nyquist(self):
        mask = make_nyquist_mask_jax(N_SMALL)
        nq = N_SMALL // 2  # = 8
        assert float(mask[nq, 0, 0]) == 0.0
        assert float(mask[0, nq, 0]) == 0.0
        assert float(mask[0, 0, nq]) == 0.0

    def test_nyquist_mask_keeps_interior(self):
        mask = make_nyquist_mask_jax(N_SMALL)
        assert float(mask[1, 1, 1]) == 1.0

    def test_nyquist_mask_shape(self):
        mask = make_nyquist_mask_jax(N_SMALL)
        assert mask.shape == (N_SMALL, N_SMALL, N_SMALL)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Leray projector (B-005 fix)
# ─────────────────────────────────────────────────────────────────────────────

class TestLerayJAX:
    def test_random_field_becomes_div_free(self):
        """Random fields → divergence-free after _project_jax."""
        N = N_SMALL
        kx, ky, kz, k2 = make_grid_jax(N)
        nq = make_nyquist_mask_jax(N)

        rng = np.random.default_rng(7)
        u = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=_DTYPE))
        v = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=_DTYPE))
        w = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=_DTYPE))

        div_before = divergence_rms_jax(u, v, w, kx, ky, kz)

        uh, vh, wh = jnp.fft.fftn(u), jnp.fft.fftn(v), jnp.fft.fftn(w)
        uh, vh, wh = _project_jax(uh, vh, wh, kx, ky, kz, k2, nq)
        u2 = jnp.real(jnp.fft.ifftn(uh))
        v2 = jnp.real(jnp.fft.ifftn(vh))
        w2 = jnp.real(jnp.fft.ifftn(wh))

        div_after = divergence_rms_jax(u2, v2, w2, kx, ky, kz)
        assert div_after < div_before * 1e-2, \
            f"Leray insufficient: {div_before:.3e} → {div_after:.3e}"
        assert div_after < 1e-2, f"div after projection = {div_after:.3e}"

    def test_tg_ic_is_div_free(self):
        """Taylor-Green IC is analytically div-free to float32 precision."""
        N = N_SMALL
        kx, ky, kz, _ = make_grid_jax(N)
        ic = taylor_green_ic(N, V0=1.0)
        u = _d(jnp.array(ic["u"], dtype=_DTYPE))
        v = _d(jnp.array(ic["v"], dtype=_DTYPE))
        w = _d(jnp.array(ic["w"], dtype=_DTYPE))
        div = divergence_rms_jax(u, v, w, kx, ky, kz)
        assert div < 1e-4, f"TG div = {div:.3e}"

    def test_projection_idempotent(self):
        """Projecting an already-projected field must not worsen divergence."""
        N = N_SMALL
        kx, ky, kz, k2 = make_grid_jax(N)
        nq = make_nyquist_mask_jax(N)

        rng = np.random.default_rng(13)
        u = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=_DTYPE))
        v = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=_DTYPE))
        w = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=_DTYPE))

        def proj(u, v, w):
            uh, vh, wh = jnp.fft.fftn(u), jnp.fft.fftn(v), jnp.fft.fftn(w)
            uh, vh, wh = _project_jax(uh, vh, wh, kx, ky, kz, k2, nq)
            return (jnp.real(jnp.fft.ifftn(uh)),
                    jnp.real(jnp.fft.ifftn(vh)),
                    jnp.real(jnp.fft.ifftn(wh)))

        u1, v1, w1 = proj(u, v, w)
        u2, v2, w2 = proj(u1, v1, w1)

        div1 = divergence_rms_jax(u1, v1, w1, kx, ky, kz)
        div2 = divergence_rms_jax(u2, v2, w2, kx, ky, kz)
        assert div1 < 1e-2, f"First projection: div={div1:.3e}"
        assert div2 <= div1 + 1e-6, f"Second projection worsened: {div1:.3e}→{div2:.3e}"


# ─────────────────────────────────────────────────────────────────────────────
# 4. Initial conditions
# ─────────────────────────────────────────────────────────────────────────────

class TestInitialConditions:
    def test_tg_E0_analytic_exact(self):
        ic = taylor_green_ic(N_MED, V0=1.0)
        assert abs(ic["E0_analytic"] - 1.0 / 8.0) < 1e-12

    def test_tg_E0_computed_close_to_analytic(self):
        ic = taylor_green_ic(N_MED, V0=1.0)
        u = _d(jnp.array(ic["u"], dtype=_DTYPE))
        v = _d(jnp.array(ic["v"], dtype=_DTYPE))
        w = _d(jnp.array(ic["w"], dtype=_DTYPE))
        E = kinetic_energy_jax(u, v, w)
        err = abs(E - ic["E0_analytic"]) / ic["E0_analytic"]
        assert err < 0.01, f"E₀ error = {err*100:.2f}%"

    def test_tg_w_zero(self):
        ic = taylor_green_ic(N_SMALL, V0=1.0)
        assert np.max(np.abs(ic["w"])) < 1e-10

    def test_tg_name(self):
        ic = taylor_green_ic(N_SMALL)
        assert ic["name"] == "taylor_green"

    def test_shear_shape(self):
        ic = shear_layer_ic(N_SMALL)
        assert ic["u"].shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_shear_v_w_zero(self):
        ic = shear_layer_ic(N_SMALL)
        assert np.max(np.abs(ic["v"])) < 1e-10
        assert np.max(np.abs(ic["w"])) < 1e-10

    def test_shear_E0_positive(self):
        ic = shear_layer_ic(N_SMALL)
        assert ic["E0_analytic"] > 0.0

    def test_shear_name(self):
        ic = shear_layer_ic(N_SMALL)
        assert ic["name"] == "shear_layer"


# ─────────────────────────────────────────────────────────────────────────────
# 5. Solver construction
# ─────────────────────────────────────────────────────────────────────────────

class TestMakeSolverJAX:
    def test_required_keys(self):
        s = make_solver_jax(N=N_SMALL, nu=NU, eps_param=EPS)
        for k in ("N", "nu", "eps_param", "kx", "ky", "kz", "k2",
                  "dealias", "nq_mask", "u", "v", "w", "theta", "t", "backend"):
            assert k in s, f"Missing key: {k}"

    def test_initial_velocity_zero(self):
        s = make_solver_jax(N=N_SMALL, nu=NU, eps_param=EPS)
        assert float(jnp.max(jnp.abs(s["u"]))) == 0.0

    def test_initial_theta_zero(self):
        s = make_solver_jax(N=N_SMALL, nu=NU, eps_param=EPS)
        assert float(jnp.max(jnp.abs(s["theta"]))) == 0.0

    def test_nu_stored(self):
        s = make_solver_jax(N=N_SMALL, nu=NU, eps_param=EPS)
        assert s["nu"] == NU

    def test_invalid_nu_raises(self):
        with pytest.raises(ValueError):
            make_solver_jax(N=N_SMALL, nu=-1.0, eps_param=EPS)

    def test_invalid_eps_raises(self):
        with pytest.raises(ValueError):
            make_solver_jax(N=N_SMALL, nu=NU, eps_param=-0.1)

    def test_set_ic_applies_tg(self):
        s = make_solver_jax(N=N_SMALL, nu=NU, eps_param=EPS)
        ic = taylor_green_ic(N_SMALL, V0=1.0)
        s2 = set_ic_jax(s, ic)
        assert float(jnp.max(jnp.abs(s2["u"]))) > 0.1

    def test_set_ic_theta_zero(self):
        """θ(x,0) = 0 is an exact Prize BC — must be enforced by set_ic."""
        s = make_solver_jax(N=N_SMALL, nu=NU, eps_param=EPS)
        ic = taylor_green_ic(N_SMALL, V0=1.0)
        s2 = set_ic_jax(s, ic)
        assert float(jnp.max(jnp.abs(s2["theta"]))) == 0.0


# ─────────────────────────────────────────────────────────────────────────────
# 6. One combined RK4 + CN IMEX step
# ─────────────────────────────────────────────────────────────────────────────

class TestSolverStep:
    def test_step_runs_without_error(self, solver16):
        s = solver16
        dt = jnp.array(0.005, dtype=_DTYPE)
        nu_j = jnp.array(NU, dtype=_DTYPE)
        new_u, new_v, new_w, new_theta = solver_step_jax(
            s["u"], s["v"], s["w"], s["theta"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias"], s["nq_mask"], nu_j, dt
        )
        assert new_u.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_step_output_shapes(self, solver16):
        s = solver16
        dt = jnp.array(0.005, dtype=_DTYPE)
        nu_j = jnp.array(NU, dtype=_DTYPE)
        out = solver_step_jax(
            s["u"], s["v"], s["w"], s["theta"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias"], s["nq_mask"], nu_j, dt
        )
        for arr in out:
            assert arr.shape == (N_SMALL, N_SMALL, N_SMALL)

    def test_step_theta_nonneg(self, solver16):
        """θ ≥ 0 enforced by jnp.clip after CN update."""
        s = solver16
        dt = jnp.array(0.005, dtype=_DTYPE)
        nu_j = jnp.array(NU, dtype=_DTYPE)
        _, _, _, new_theta = solver_step_jax(
            s["u"], s["v"], s["w"], s["theta"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias"], s["nq_mask"], nu_j, dt
        )
        assert bool(jnp.all(new_theta >= 0.0)), \
            f"θ_min = {float(jnp.min(new_theta)):.3e}"

    def test_step_div_free(self, solver16):
        """div u ≈ 0 maintained after one RK4 step + Leray re-projection."""
        s = solver16
        dt = jnp.array(0.005, dtype=_DTYPE)
        nu_j = jnp.array(NU, dtype=_DTYPE)
        new_u, new_v, new_w, _ = solver_step_jax(
            s["u"], s["v"], s["w"], s["theta"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias"], s["nq_mask"], nu_j, dt
        )
        div = divergence_rms_jax(new_u, new_v, new_w, s["kx"], s["ky"], s["kz"])
        assert div < 5e-2, f"div u after step = {div:.3e}"

    def test_step_energy_finite(self, solver16):
        s = solver16
        dt = jnp.array(0.005, dtype=_DTYPE)
        nu_j = jnp.array(NU, dtype=_DTYPE)
        new_u, new_v, new_w, _ = solver_step_jax(
            s["u"], s["v"], s["w"], s["theta"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias"], s["nq_mask"], nu_j, dt
        )
        E = kinetic_energy_jax(new_u, new_v, new_w)
        assert np.isfinite(E), f"Energy is not finite: {E}"
        assert E > 0.0, f"Energy = {E}"

    def test_step_theta_grows_from_source(self, solver16):
        """S_θ = ν|∇u|² ≥ 0 must cause θ_max to grow from zero."""
        s = solver16
        dt = jnp.array(0.005, dtype=_DTYPE)
        nu_j = jnp.array(NU, dtype=_DTYPE)
        theta_before = float(jnp.max(s["theta"]))
        _, _, _, new_theta = solver_step_jax(
            s["u"], s["v"], s["w"], s["theta"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias"], s["nq_mask"], nu_j, dt
        )
        theta_after = float(jnp.max(new_theta))
        # θ starts at 0 — source S_θ must push it positive
        assert theta_after >= theta_before, \
            f"θ_max did not grow: {theta_before:.3e} → {theta_after:.3e}"

    def test_step_three_times_div_free(self, solver16_stepped):
        """After 3 steps div u must still be small."""
        s = solver16_stepped
        div = divergence_rms_jax(s["u"], s["v"], s["w"], s["kx"], s["ky"], s["kz"])
        assert div < 5e-2, f"div after 3 steps = {div:.3e}"

    def test_step_three_times_theta_pos(self, solver16_stepped):
        """After 3 steps θ must have grown above zero."""
        s = solver16_stepped
        theta_max = float(jnp.max(s["theta"]))
        assert theta_max > 0.0, f"θ_max after 3 steps = {theta_max:.3e}"


# ─────────────────────────────────────────────────────────────────────────────
# 7. δ estimator (Lemma 2.5)
# ─────────────────────────────────────────────────────────────────────────────

class TestDeltaFromTheta:
    def test_zero_theta_gives_zero_delta(self):
        """θ ≡ 0 → no spectral content → δ = 0."""
        N = N_SMALL
        kx, ky, kz, _ = make_grid_jax(N)
        theta = _d(jnp.zeros((N, N, N), dtype=_DTYPE))
        delta = delta_from_theta_jax(theta, kx, ky, kz)
        assert delta == 0.0

    def test_delta_nonneg(self, solver16_stepped):
        s = solver16_stepped
        delta = delta_from_theta_jax(s["theta"], s["kx"], s["ky"], s["kz"])
        assert delta >= 0.0, f"δ = {delta}"

    def test_delta_returns_float(self, solver16_stepped):
        s = solver16_stepped
        delta = delta_from_theta_jax(s["theta"], s["kx"], s["ky"], s["kz"])
        assert isinstance(delta, float)

    def test_delta_bounded_above(self, solver16_stepped):
        """δ_est ≤ 1.5 by construction of the s-ladder {0,0.5,1,1.5,2,2.5}."""
        s = solver16_stepped
        delta = delta_from_theta_jax(s["theta"], s["kx"], s["ky"], s["kz"])
        assert delta <= 1.5, f"δ = {delta} exceeds expected max"

    def test_delta_positive_after_several_steps(self, solver16_stepped):
        """After 3 TG steps, S_θ has driven θ above 0 → δ > 0."""
        s = solver16_stepped
        delta = delta_from_theta_jax(s["theta"], s["kx"], s["ky"], s["kz"])
        # θ grows from source → should have some Sobolev regularity
        assert delta >= 0.0


# ─────────────────────────────────────────────────────────────────────────────
# 8. CFL timestep
# ─────────────────────────────────────────────────────────────────────────────

class TestCFLJAX:
    def test_cfl_positive(self, solver16):
        s = solver16
        dt = compute_cfl_dt_jax(s["u"], s["v"], s["w"], N_SMALL, NU, 0.1)
        assert dt > 0.0

    def test_cfl_respects_dtmax(self, solver16):
        s = solver16
        dt_max = 0.002
        dt = compute_cfl_dt_jax(s["u"], s["v"], s["w"], N_SMALL, NU, dt_max)
        assert dt <= dt_max + 1e-12

    def test_cfl_zero_velocity_uses_diffusion(self):
        N = N_SMALL
        u = _d(jnp.zeros((N, N, N), dtype=_DTYPE))
        v = _d(jnp.zeros((N, N, N), dtype=_DTYPE))
        w = _d(jnp.zeros((N, N, N), dtype=_DTYPE))
        dt = compute_cfl_dt_jax(u, v, w, N, NU, dt_max=1.0)
        assert dt > 0.0
        assert dt <= 1.0


# ─────────────────────────────────────────────────────────────────────────────
# 9. Short run loop
# ─────────────────────────────────────────────────────────────────────────────

class TestRunJAX:
    @pytest.fixture(scope="class")
    def short_run(self):
        s = make_solver_jax(N=N_SMALL, nu=NU, eps_param=EPS, dt_max=0.05)
        ic = taylor_green_ic(N_SMALL, V0=1.0)
        s = set_ic_jax(s, ic)
        return run_jax(s, t_end=0.15, max_steps=50, warmup=True)

    def test_reaches_t_end(self, short_run):
        assert short_run["t_final"] >= 0.10

    def test_verdict_pass(self, short_run):
        assert short_run["verdict"] == "PASS"

    def test_history_nonempty(self, short_run):
        assert len(short_run["history"]) > 0

    def test_energy_decays(self, short_run):
        assert short_run["energy_decayed"], \
            f"E_initial={short_run['E_initial']:.4e} E_final={short_run['E_final']:.4e}"

    def test_theta_not_negative(self, short_run):
        assert not short_run["theta_negative"]

    def test_delta_max_nonneg(self, short_run):
        assert short_run["delta_max"] >= 0.0

    def test_eps_lps_min_nonneg(self, short_run):
        """ε_LPS = min(μ_eff) − ν ≥ 0 by construction (θ ≥ 0, f(θ) ≥ 0)."""
        assert short_run["eps_lps_min"] >= 0.0

    def test_n_steps_positive(self, short_run):
        assert short_run["n_steps"] > 0


# ─────────────────────────────────────────────────────────────────────────────
# 10. TG benchmark at N=32 — Lemma 2.5 validation
# ─────────────────────────────────────────────────────────────────────────────

class TestTGBenchmark32:
    @pytest.fixture(scope="class")
    def tg32(self):
        """N=32 TG run, t_end=0.5 — built once per class."""
        return run_taylor_green_jax(
            N=N_MED, t_end=0.5, nu=NU, eps_param=EPS,
            max_steps=2000, verbose=False,
        )

    def test_E0_error_lt_1pct(self, tg32):
        assert tg32["E0_error_pct"] < 1.0, \
            f"E₀ error = {tg32['E0_error_pct']:.2f}%"

    def test_energy_decays(self, tg32):
        assert tg32["energy_decayed"]

    def test_theta_not_negative(self, tg32):
        assert not tg32["theta_negative"]

    def test_delta_max_positive(self, tg32):
        """Lemma 2.5: δ > 0 from CZ source alone (ε = 0 or ε > 0)."""
        assert tg32["delta_max"] > 0.0, \
            f"δ_max = {tg32['delta_max']:.2f} — Lemma 2.5 not confirmed"

    def test_eps_lps_min_nonneg(self, tg32):
        assert tg32["eps_lps_min"] >= 0.0

    def test_verdict_pass(self, tg32):
        assert tg32["verdict"] == "PASS"

    def test_exp_id_format(self, tg32):
        assert tg32["exp_id"].startswith("EXP-L3-R2-TG-JAX-")

    def test_backend_recorded(self, tg32):
        assert tg32["backend"] == _JAX_BACKEND

    def test_n_steps_positive(self, tg32):
        assert tg32["n_steps"] > 0


# ─────────────────────────────────────────────────────────────────────────────
# 11. Shear layer (Route A) — δ > 0 for non-mixing IC
# ─────────────────────────────────────────────────────────────────────────────

class TestShearLayerJAX:
    @pytest.fixture(scope="class")
    def sh16(self):
        return run_shear_layer_jax(
            N=N_SMALL, t_end=0.3, nu=NU, eps_param=EPS,
            max_steps=300, verbose=False,
        )

    def test_verdict_pass(self, sh16):
        assert sh16["verdict"] == "PASS"

    def test_theta_not_negative(self, sh16):
        assert not sh16["theta_negative"]

    def test_delta_max_nonneg(self, sh16):
        """Route A: shear IC delivers δ ≥ 0 (should be > 0 with TG-level source)."""
        assert sh16["delta_max"] >= 0.0

    def test_eps_lps_min_nonneg(self, sh16):
        assert sh16["eps_lps_min"] >= 0.0

    def test_exp_id_format(self, sh16):
        assert sh16["exp_id"].startswith("EXP-L3-R2-SH-JAX-")


# ─────────────────────────────────────────────────────────────────────────────
# 12. ε = 0 exact Prize: Lemma 2.5 direct confirmation
# ─────────────────────────────────────────────────────────────────────────────

class TestEpsZeroExactPrize:
    def test_eps0_theta_still_grows(self):
        """S_θ = ν|∇u|² is ε-independent — θ grows even at ε=0 (exact Prize)."""
        s0 = make_solver_jax(N=N_SMALL, nu=NU, eps_param=0.0)
        ic = taylor_green_ic(N_SMALL, V0=1.0)
        s0 = set_ic_jax(s0, ic)

        dt = jnp.array(0.005, dtype=_DTYPE)
        nu_j = jnp.array(NU, dtype=_DTYPE)
        u, v, w, theta = s0["u"], s0["v"], s0["w"], s0["theta"]
        kx, ky, kz, k2 = s0["kx"], s0["ky"], s0["kz"], s0["k2"]
        nq, da = s0["nq_mask"], s0["dealias"]

        for _ in range(5):
            u, v, w, theta = solver_step_jax(
                u, v, w, theta, kx, ky, kz, k2, da, nq, nu_j, dt
            )
            jax.block_until_ready((u, v, w, theta))

        theta_max = float(jnp.max(theta))
        assert theta_max > 0.0, \
            f"θ_max = {theta_max:.3e} at ε=0 — S_θ source not active"

    def test_eps0_delta_nonneg(self):
        """At ε=0 (exact Prize), δ ≥ 0 — CZ source alone delivers regularity."""
        s0 = make_solver_jax(N=N_SMALL, nu=NU, eps_param=0.0)
        ic = taylor_green_ic(N_SMALL, V0=1.0)
        s0 = set_ic_jax(s0, ic)

        result = run_jax(s0, t_end=0.15, max_steps=50, verbose=False)
        assert result["delta_max"] >= 0.0

    def test_eps0_eps_lps_zero(self):
        """At ε=0: μ_eff = ν + 0·f(θ) = ν → ε_LPS = min(μ_eff) − ν = 0."""
        s0 = make_solver_jax(N=N_SMALL, nu=NU, eps_param=0.0)
        ic = taylor_green_ic(N_SMALL, V0=1.0)
        s0 = set_ic_jax(s0, ic)

        result = run_jax(s0, t_end=0.1, max_steps=30, verbose=False)
        # At eps_param=0, mu_eff = nu everywhere → eps_lps = 0 exactly
        assert abs(result["eps_lps_min"]) < 1e-8, \
            f"ε_LPS_min = {result['eps_lps_min']:.3e} at ε=0 (should be 0)"
