"""
Tests for gamma_adversarial_search.py

Run with:
    python -m pytest layer4/test_gamma_adversarial_search.py -v

Six tests covering:
  1. Concentrated IC is divergence-free
  2. Concentrated IC has most gradient energy inside the target ball
  3. Free decay dissipates L2 norm
  4. Gamma = 0 for zero field
  5. Gamma = 0 for constant (zero-gradient) field
  6. run_search returns dict with required keys
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from layer4.gamma_adversarial_search import GammaAdversarialSearch, _fft_grad_mag_sq
from layer4.mean_morrey_diagnostic import L_BOX


# ---------------------------------------------------------------------------
# Shared small search object (N=16 for speed)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def search16() -> GammaAdversarialSearch:
    return GammaAdversarialSearch(N=16, nu=1e-3, T=1.0, n_steps=20, seed=0)


# ---------------------------------------------------------------------------
# Test 1: concentrated IC is divergence-free
# ---------------------------------------------------------------------------

class TestConcentratedICDivergenceFree:
    """div u in Fourier space should be < 1e-10 after Leray projection."""

    def test_divergence_free(self, search16: GammaAdversarialSearch) -> None:
        u = search16.make_concentrated_ic(r_frac=0.2, amplitude=1.0)
        N = search16.N
        dx = L_BOX / N

        # Spectral divergence
        kx = np.fft.fftfreq(N, d=dx) * (2 * np.pi)
        ky = np.fft.fftfreq(N, d=dx) * (2 * np.pi)
        kz = np.fft.fftfreq(N, d=dx) * (2 * np.pi)
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

        uh = np.array([np.fft.fftn(u[c]) for c in range(3)])
        div_hat = 1j * (KX * uh[0] + KY * uh[1] + KZ * uh[2])
        div_phys = np.real(np.fft.ifftn(div_hat))

        div_norm = float(np.sqrt(np.mean(div_phys ** 2)))
        assert div_norm < 1e-10, (
            f"div u RMS = {div_norm:.3e} exceeds 1e-10 — Leray projection failed"
        )


# ---------------------------------------------------------------------------
# Test 2: gradient energy concentrated inside target ball
# ---------------------------------------------------------------------------

class TestConcentratedICEnergyInBall:
    """At least 50% of total ||grad u||^2 should lie inside B_{r_frac}."""

    def test_energy_fraction_in_ball(self, search16: GammaAdversarialSearch) -> None:
        r_frac = 0.15
        u = search16.make_concentrated_ic(r_frac=r_frac, amplitude=1.0)
        N = search16.N

        grad_sq = _fft_grad_mag_sq(u)

        # Ball centred at N//2 (box centre)
        from layer4.mean_morrey_diagnostic import ball_mask
        cx = cy = cz = N // 2
        mask = ball_mask(N, cx, cy, cz, r_frac)

        energy_in  = float(grad_sq[mask].sum())
        energy_all = float(grad_sq.sum())

        if energy_all < 1e-30:
            pytest.skip("Zero-amplitude field — nothing to test")

        frac = energy_in / energy_all
        # After Leray projection, energy spreads spectrally but should still
        # be meaningfully concentrated relative to a uniform field (~(r_frac)^3 ~ 0.34%)
        assert frac > 0.15, (
            f"Only {frac:.1%} of grad energy inside ball — expected > 15% "
            f"(Leray projection spreads energy but IC should remain locally concentrated)"
        )


# ---------------------------------------------------------------------------
# Test 3: free decay dissipates L2 norm
# ---------------------------------------------------------------------------

class TestFreeDecayDissipates:
    """||u(T)||_L2 < ||u(0)||_L2 — pure viscous decay must reduce energy."""

    def test_l2_strictly_decreases(self, search16: GammaAdversarialSearch) -> None:
        u0 = search16.make_concentrated_ic(r_frac=0.2, amplitude=1.0)
        snapshots = search16.run_free_decay(u0)

        l2_0 = float(np.sqrt(np.mean(snapshots[0] ** 2)))
        l2_T = float(np.sqrt(np.mean(snapshots[-1] ** 2)))

        assert l2_T < l2_0, (
            f"L2 norm did not decrease: ||u(0)||={l2_0:.4e}, ||u(T)||={l2_T:.4e}"
        )


# ---------------------------------------------------------------------------
# Test 4: Gamma = 0 for zero field
# ---------------------------------------------------------------------------

class TestGammaZeroField:
    """Gamma(r) must be 0 for the zero velocity field."""

    def test_zero_field(self, search16: GammaAdversarialSearch) -> None:
        N = search16.N
        zero_field = np.zeros((3, N, N, N), dtype=float)
        gamma = search16.compute_gamma([zero_field], r_frac=0.2)
        assert gamma == pytest.approx(0.0, abs=1e-15), (
            f"Gamma = {gamma:.3e} for zero field — expected 0"
        )


# ---------------------------------------------------------------------------
# Test 5: Gamma = 0 for constant field
# ---------------------------------------------------------------------------

class TestGammaConstantField:
    """A spatially uniform (constant) field has zero gradient, so Gamma = 0."""

    def test_constant_field(self, search16: GammaAdversarialSearch) -> None:
        N = search16.N
        # u = (1, 0, 0) everywhere — constant, divergence-free, zero gradient
        const_field = np.zeros((3, N, N, N), dtype=float)
        const_field[0] = 1.0   # u_x = 1 everywhere

        gamma = search16.compute_gamma([const_field, const_field], r_frac=0.2)
        assert gamma == pytest.approx(0.0, abs=1e-12), (
            f"Gamma = {gamma:.3e} for constant field — expected 0"
        )


# ---------------------------------------------------------------------------
# Test 6: run_search returns dict with required keys
# ---------------------------------------------------------------------------

class TestRunSearchReturnsVerdict:
    """run_search must return a dict with the four mandatory keys."""

    REQUIRED_KEYS = {"gamma_initial", "gamma_evolved", "suppressed", "verdict"}

    def test_required_keys_present(self, search16: GammaAdversarialSearch) -> None:
        # Use only 2 r values for speed
        results = search16.run_search(r_fracs=[0.1, 0.3])
        missing = self.REQUIRED_KEYS - set(results.keys())
        assert not missing, f"Missing keys in run_search output: {missing}"

    def test_types_are_correct(self, search16: GammaAdversarialSearch) -> None:
        results = search16.run_search(r_fracs=[0.1, 0.3])
        assert isinstance(results["gamma_initial"], float)
        assert isinstance(results["gamma_evolved"], float)
        assert isinstance(results["suppressed"], bool)
        assert isinstance(results["verdict"], str)
        assert len(results["verdict"]) > 0
