"""
test_route1_3D_jax.py — Test suite for JAX Metal Route 1 3D solver
====================================================================
Tests for route1_3D_jax.py.  Validates:
  - JAX backend available
  - Grid / mask correctness (JAX arrays)
  - Leray projector (Nyquist-safe, Hermitian-preserving)
  - TG IC energy accuracy
  - One RK4 step: time advances, div-free maintained, T in range
  - Run loop: energy decays, A_min > 0
  - TG benchmark: PASS verdict

All tests use N=16 (fast) except the benchmark (N=16, t_end=0.1).
"""

import sys
import os
import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))

# ── Guard: skip entire module if JAX not installed ─────────────────────────────
jax = pytest.importorskip("jax", reason="JAX not installed")
# float32 throughout — must NOT enable x64 (Metal does not support float64).
import jax.numpy as jnp

from route1_3D_jax import (
    _JAX_AVAILABLE,
    _JAX_BACKEND,
    _COMPUTE_DEVICE,
    _d,
    make_grid_jax,
    make_dealias_mask_jax,
    make_nyquist_mask_jax,
    ideal_props_jax,
    _project_spectral_jax,
    divergence_rms_jax,
    kinetic_energy_jax,
    taylor_green_ic_jax,
    make_jax_solver,
    set_ic_jax,
    solver_step_jax,
    run_jax,
    run_taylor_green_jax,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _dtype():
    """float32 — jax-metal 0.1.0 does not support float64 on Metal."""
    return jnp.float32


def _atol():
    """Absolute tolerance for divergence-free tests (float32 CPU XLA)."""
    return 1e-5


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def solver16():
    """N=16 solver with TG IC — module scope (built once)."""
    s = make_jax_solver(N=16, T_base=300.0, dt_max=0.05)
    ic = taylor_green_ic_jax(16, V0=1.0, T_base=300.0, dtype=s["dtype"])
    return set_ic_jax(s, ic)


@pytest.fixture(scope="module")
def tg_result16():
    """Short TG run at N=16, t_end=0.1."""
    return run_taylor_green_jax(N=16, t_end=0.1, max_steps=30, verbose=False)


# ─────────────────────────────────────────────────────────────────────────────
# 1. Backend availability
# ─────────────────────────────────────────────────────────────────────────────

class TestJAXBackend:
    def test_jax_available(self):
        assert _JAX_AVAILABLE, "JAX not available"

    def test_backend_detected(self):
        # CPU-JIT = Metal detected but routing to CPU for FFT support
        valid = ("cpu", "gpu", "tpu", "METAL", "CPU-JIT")
        assert _JAX_BACKEND in valid, f"Unknown backend: {_JAX_BACKEND}"

    def test_jit_works(self):
        """A simple JIT-compiled function must return the correct result."""
        @jax.jit
        def add(x, y):
            return x + y
        result = float(add(jnp.array(1.0), jnp.array(2.0)))
        assert abs(result - 3.0) < 1e-6

    def test_fft_on_compute_device(self):
        """jnp.fft.fftn must run on the compute device (CPU XLA) without error.

        Note: Metal (jax-metal 0.1.0) does not support complex FFT.
        We always route FFT to the CPU device — this test verifies that works.
        """
        x = _d(jnp.ones((8, 8, 8), dtype=jnp.float32))
        out = jnp.fft.fftn(x)
        assert out.shape == (8, 8, 8)
        # DC mode of FFT of all-ones is N³ = 512
        assert abs(float(jnp.real(out[0, 0, 0])) - 8.0 ** 3) < 1.0


# ─────────────────────────────────────────────────────────────────────────────
# 2. Grid and masks
# ─────────────────────────────────────────────────────────────────────────────

class TestGridJAX:
    def test_shapes(self):
        kx, ky, kz, k2 = make_grid_jax(16)
        for arr in (kx, ky, kz, k2):
            assert arr.shape == (16, 16, 16)

    def test_k2_nonnegative(self):
        _, _, _, k2 = make_grid_jax(16)
        assert bool(jnp.all(k2 >= 0.0))

    def test_k2_zero_at_dc(self):
        _, _, _, k2 = make_grid_jax(16)
        assert float(k2[0, 0, 0]) == 0.0

    def test_dealias_mask_shape(self):
        mask = make_dealias_mask_jax(16)
        assert mask.shape == (16, 16, 16)

    def test_dealias_keeps_low_k(self):
        mask = make_dealias_mask_jax(16)
        assert float(mask[1, 1, 1]) == 1.0

    def test_dealias_zeros_high_k(self):
        mask = make_dealias_mask_jax(16)
        # k = (8, 0, 0): |k| = 8 > N//3 = 5 → should be zeroed
        assert float(mask[8, 0, 0]) == 0.0

    def test_nyquist_mask_zeros_nyquist(self):
        mask = make_nyquist_mask_jax(16)
        # All Nyquist-plane modes should be 0
        assert float(mask[8, 0, 0]) == 0.0
        assert float(mask[0, 8, 0]) == 0.0
        assert float(mask[0, 0, 8]) == 0.0

    def test_nyquist_mask_keeps_interior(self):
        mask = make_nyquist_mask_jax(16)
        assert float(mask[1, 1, 1]) == 1.0


# ─────────────────────────────────────────────────────────────────────────────
# 3. Fluid properties
# ─────────────────────────────────────────────────────────────────────────────

class TestPropsJAX:
    def test_mu_positive(self):
        T = jnp.array([200.0, 300.0, 500.0])
        p = ideal_props_jax(T)
        assert bool(jnp.all(p["mu"] > 0))

    def test_rho_positive(self):
        T = jnp.array([200.0, 300.0, 500.0])
        p = ideal_props_jax(T)
        assert bool(jnp.all(p["rho"] > 0))

    def test_A_positive(self):
        """A(T) = k/(ρ·cv) > 0 — second law."""
        T = jnp.linspace(200.0, 600.0, 10)
        p = ideal_props_jax(T)
        A = p["k"] / (p["rho"] * p["cv"])
        assert bool(jnp.all(A > 0))

    def test_field_shapes(self):
        T = jnp.ones((4, 4, 4)) * 300.0
        p = ideal_props_jax(T)
        assert p["mu"].shape == (4, 4, 4)
        assert p["rho"].shape == (4, 4, 4)


# ─────────────────────────────────────────────────────────────────────────────
# 4. Leray projector
# ─────────────────────────────────────────────────────────────────────────────

class TestLerayJAX:
    def test_projection_gives_div_free_random(self):
        """Random fields on compute device → divergence-free after Leray projection."""
        N = 16
        kx, ky, kz, k2 = make_grid_jax(N)
        nqmask = make_nyquist_mask_jax(N)

        # Create random fields on compute device (CPU when Metal is default)
        rng = np.random.default_rng(42)
        u = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=jnp.float32))
        v = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=jnp.float32))
        w = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=jnp.float32))

        div_before = divergence_rms_jax(u, v, w, kx, ky, kz)

        u_hat, v_hat, w_hat = (jnp.fft.fftn(u), jnp.fft.fftn(v), jnp.fft.fftn(w))
        u_hat, v_hat, w_hat = _project_spectral_jax(u_hat, v_hat, w_hat,
                                                     kx, ky, kz, k2, nqmask)
        u_p = jnp.real(jnp.fft.ifftn(u_hat))
        v_p = jnp.real(jnp.fft.ifftn(v_hat))
        w_p = jnp.real(jnp.fft.ifftn(w_hat))

        div_after = divergence_rms_jax(u_p, v_p, w_p, kx, ky, kz)

        # float32: project reduces divergence by many orders of magnitude
        assert div_after < div_before * 1e-3, (
            f"Leray projection insufficient: before={div_before:.3e}, after={div_after:.3e}"
        )
        assert div_after < 1e-3, (
            f"div u after projection = {div_after:.3e}"
        )

    def test_tg_stays_div_free(self):
        """TG IC is analytically divergence-free."""
        N = 16
        kx, ky, kz, k2 = make_grid_jax(N)
        ic = taylor_green_ic_jax(N, V0=1.0, T_base=300.0)
        div = divergence_rms_jax(ic["u"], ic["v"], ic["w"], kx, ky, kz)
        # float32 FFT: machine eps ~1e-7 per mode; N=16³=4096 modes → ~1e-5 total
        assert div < 1e-4, f"TG div = {div:.3e}"

    def test_projection_idempotent(self):
        """Projecting an already-projected field: div unchanged and tiny."""
        N = 16
        kx, ky, kz, k2 = make_grid_jax(N)
        nqmask = make_nyquist_mask_jax(N)

        rng = np.random.default_rng(99)
        u = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=jnp.float32))
        v = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=jnp.float32))
        w = _d(jnp.array(rng.standard_normal((N, N, N)), dtype=jnp.float32))

        def project_phys(u, v, w):
            uh, vh, wh = jnp.fft.fftn(u), jnp.fft.fftn(v), jnp.fft.fftn(w)
            uh, vh, wh = _project_spectral_jax(uh, vh, wh, kx, ky, kz, k2, nqmask)
            return jnp.real(jnp.fft.ifftn(uh)), jnp.real(jnp.fft.ifftn(vh)), jnp.real(jnp.fft.ifftn(wh))

        u1, v1, w1 = project_phys(u, v, w)
        u2, v2, w2 = project_phys(u1, v1, w1)

        div1 = divergence_rms_jax(u1, v1, w1, kx, ky, kz)
        div2 = divergence_rms_jax(u2, v2, w2, kx, ky, kz)

        assert div1 < 1e-3, f"First projection: div={div1:.3e}"
        assert div2 < 1e-3, f"Second projection: div={div2:.3e}"
        assert div2 <= div1 + 1e-12, f"Second projection worsened: {div1:.3e} → {div2:.3e}"


# ─────────────────────────────────────────────────────────────────────────────
# 5. Initial conditions and solver construction
# ─────────────────────────────────────────────────────────────────────────────

class TestICJAX:
    def test_tg_E0_analytic(self):
        ic = taylor_green_ic_jax(32, V0=1.0, T_base=300.0)
        # E0_analytic is a Python float (V0²/8); exact to floating-point
        assert abs(ic["E0_analytic"] - 1.0 / 8.0) < 1e-12

    def test_tg_E0_computed_matches_analytic(self):
        ic = taylor_green_ic_jax(32, V0=1.0, T_base=300.0)
        E = kinetic_energy_jax(ic["u"], ic["v"], ic["w"])
        rel_err = abs(E - ic["E0_analytic"]) / ic["E0_analytic"]
        assert rel_err < 0.01, f"E₀ relative error = {rel_err*100:.2f}%"

    def test_tg_w_zero(self):
        ic = taylor_green_ic_jax(16)
        # w = np.zeros(...) converted to float32 → exactly 0.0
        assert float(jnp.max(jnp.abs(ic["w"]))) < 1e-6

    def test_solver_keys(self):
        s = make_jax_solver(N=16)
        for k in ("N", "kx", "ky", "kz", "k2", "dealias_mask",
                  "nyquist_mask", "u", "v", "w", "T", "nu_ref", "A_ref"):
            assert k in s

    def test_A_ref_positive(self):
        s = make_jax_solver(N=16)
        assert s["A_ref"] > 0.0

    def test_set_ic_immutable(self):
        """set_ic_jax must not mutate the original solver."""
        s = make_jax_solver(N=16)
        ic = taylor_green_ic_jax(16)
        s2 = set_ic_jax(s, ic)
        assert float(jnp.max(jnp.abs(s["u"]))) == 0.0, "set_ic mutated original"
        assert float(jnp.max(jnp.abs(s2["u"]))) > 0.0, "set_ic did not apply IC"


# ─────────────────────────────────────────────────────────────────────────────
# 6. One RK4 step
# ─────────────────────────────────────────────────────────────────────────────

class TestStepJAX:
    def test_step_runs(self, solver16):
        """solver_step_jax must complete without error."""
        s = solver16
        dt = jnp.array(0.01, dtype=s["dtype"])
        new_u, new_v, new_w, new_T = solver_step_jax(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["nyquist_mask"], dt
        )
        assert new_u.shape == s["u"].shape

    def test_step_T_in_range(self, solver16):
        s = solver16
        dt = jnp.array(0.01, dtype=s["dtype"])
        _, _, _, new_T = solver_step_jax(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["nyquist_mask"], dt
        )
        assert bool(jnp.all(new_T >= 100.0))
        assert bool(jnp.all(new_T <= 6000.0))

    def test_step_div_free(self, solver16):
        """div u must stay ≈ 0 after one RK4 step."""
        s = solver16
        dt = jnp.array(0.01, dtype=s["dtype"])
        new_u, new_v, new_w, _ = solver_step_jax(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["nyquist_mask"], dt
        )
        div = divergence_rms_jax(new_u, new_v, new_w, s["kx"], s["ky"], s["kz"])
        tol = _atol() * 1000  # float32: ~1e-2 after projection; float64: ~1e-7
        assert div < tol, f"div u after step = {div:.3e}"

    def test_step_A_min_positive(self, solver16):
        s = solver16
        dt = jnp.array(0.01, dtype=s["dtype"])
        _, _, _, new_T = solver_step_jax(
            s["u"], s["v"], s["w"], s["T"],
            s["kx"], s["ky"], s["kz"], s["k2"],
            s["dealias_mask"], s["nyquist_mask"], dt
        )
        props = ideal_props_jax(new_T)
        A = props["k"] / (props["rho"] * props["cv"])
        assert float(jnp.min(A)) > 0.0


# ─────────────────────────────────────────────────────────────────────────────
# 7. Run loop
# ─────────────────────────────────────────────────────────────────────────────

class TestRunJAX:
    @pytest.fixture(scope="class")
    def short_run(self):
        s = make_jax_solver(N=16, dt_max=0.05)
        ic = taylor_green_ic_jax(16, V0=1.0, T_base=300.0, dtype=s["dtype"])
        s = set_ic_jax(s, ic)
        return run_jax(s, t_end=0.1, max_steps=20, warmup=False)

    def test_reaches_t_end(self, short_run):
        assert short_run["t_final"] >= 0.09

    def test_verdict_pass(self, short_run):
        assert short_run["verdict"] == "PASS"

    def test_history_nonempty(self, short_run):
        assert len(short_run["history"]) > 0

    def test_A_min_positive(self, short_run):
        assert short_run["A_min_global"] > 0.0

    def test_energy_decays(self, short_run):
        assert short_run["E_final"] < short_run["E_initial"]


# ─────────────────────────────────────────────────────────────────────────────
# 8. TG benchmark
# ─────────────────────────────────────────────────────────────────────────────

class TestTGBenchmarkJAX:
    def test_E0_accuracy(self, tg_result16):
        assert tg_result16["E0_error_pct"] < 1.0

    def test_energy_decays(self, tg_result16):
        assert tg_result16["energy_decayed"] is True

    def test_A_min_positive(self, tg_result16):
        assert tg_result16["A_min_global"] > 0.0

    def test_verdict_pass(self, tg_result16):
        assert tg_result16["verdict"] == "PASS"

    def test_exp_id_format(self, tg_result16):
        assert tg_result16["exp_id"].startswith("EXP-L3-R1-TG-JAX-")

    def test_backend_recorded(self, tg_result16):
        assert tg_result16["backend"] == _JAX_BACKEND
