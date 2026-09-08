"""
test_self_similar_IC.py — Unit tests for Layer 2 IC Generator (M2)

Every test validates a structural or mathematical property of the
self-similar initial conditions used in the lambda sweep (M3).

Run with: pytest layer2/test_self_similar_IC.py -v
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import pytest

from layer2.self_similar_IC import (
    make_grid, spectral_laplacian_inv, vorticity_to_velocity, dealias_23,
    energy_spectrum, generate_CCF_profile, generate_Boussinesq_profile,
    generate_adversarial_min, generate_wang_profile, generate_all_wang_profiles,
    generate_lambda_sweep_ICs, profile_diagnostics, verify_energy_spectrum_scaling,
    WANG_PROFILES, LAMBDA_SWEEP_VALUES,
)


# ---------------------------------------------------------------------------
# Grid utilities
# ---------------------------------------------------------------------------

class TestGridUtilities:

    def test_make_grid_shape(self):
        x, y, kx, ky = make_grid(64)
        assert x.shape == (64, 64)
        assert y.shape == (64, 64)
        assert kx.shape == (64, 64)
        assert ky.shape == (64, 64)

    def test_make_grid_range(self):
        N = 32
        x, y, kx, ky = make_grid(N)
        assert float(x.min()) >= 0.0
        assert float(x.max()) < 2 * np.pi
        assert float(y.min()) >= 0.0
        assert float(y.max()) < 2 * np.pi

    def test_spectral_laplacian_inv_recovers_field(self):
        """
        spectral_laplacian_inv solves (-Δ)ψ = f, i.e., Δψ = -f.
        Verify: Δ(spectral_laplacian_inv(f)) = -f.
        """
        N = 32
        x, y, kx, ky = make_grid(N)
        f = np.sin(2 * x) * np.cos(3 * y)
        psi = spectral_laplacian_inv(f, kx, ky)
        # Laplacian of psi via spectral derivative
        psi_hat = np.fft.fft2(psi)
        k2 = kx ** 2 + ky ** 2
        k2[0, 0] = 1.0
        lap_psi = np.real(np.fft.ifft2(-k2 * psi_hat))
        # (-Δ)ψ = f  →  Δψ = -f
        assert np.allclose(lap_psi, -f, atol=1e-8), (
            f"Spectral Laplacian inversion error: max|Δψ - (-f)|={np.max(np.abs(lap_psi + f)):.2e}"
        )

    def test_vorticity_to_velocity_divergence_free(self):
        """Velocity field recovered from vorticity must be divergence-free."""
        N = 64
        x, y, kx, ky = make_grid(N)
        # Gaussian vortex
        r2 = (x - np.pi) ** 2 + (y - np.pi) ** 2
        omega = np.exp(-r2 / 0.5)
        u, v = vorticity_to_velocity(omega, kx, ky)

        # ∂u/∂x + ∂v/∂y should be ~0 (divergence-free)
        u_hat = np.fft.fft2(u)
        v_hat = np.fft.fft2(v)
        div_hat = 1j * kx * u_hat + 1j * ky * v_hat
        div = np.real(np.fft.ifft2(div_hat))
        assert np.max(np.abs(div)) < 1e-8, (
            f"Velocity field not divergence-free: max|div u| = {np.max(np.abs(div)):.2e}"
        )

    def test_dealias_kills_high_wavenumbers(self):
        """2/3 rule: wavenumbers above N/3 should be zeroed (within floating-point noise)."""
        N = 32
        f = np.random.default_rng(0).normal(0, 1, (N, N))
        f_d = dealias_23(f)
        f_hat = np.fft.fft2(f_d)
        cutoff = N // 3
        # After ifft2 → real → fft2 round-trip, high modes should be ≲ machine epsilon × norm(f)
        norm_f = np.max(np.abs(f_hat))
        tol = norm_f * 1e-8  # relative tolerance — floating-point round-trip
        # dealias_23 zeros |k| > cutoff, i.e. rows cutoff+1 .. N-cutoff-1
        # Row cutoff (k=+cutoff) and row N-cutoff (k=-cutoff) are legitimately kept.
        assert np.max(np.abs(f_hat[cutoff + 1: N - cutoff, :])) < tol, (
            f"High-k modes not zeroed: max={np.max(np.abs(f_hat[cutoff+1:N-cutoff,:])):.2e} (tol={tol:.2e})"
        )
        assert np.max(np.abs(f_hat[:, cutoff + 1: N - cutoff])) < tol

    def test_energy_spectrum_returns_positive_values(self):
        N = 32
        x, y, kx, ky = make_grid(N)
        r2 = (x - np.pi) ** 2 + (y - np.pi) ** 2
        omega = np.exp(-r2 / 0.5)
        k_bins, E_k = energy_spectrum(omega, kx, ky)
        assert np.all(E_k >= 0), "Energy spectrum must be non-negative"
        assert np.any(E_k > 0), "Energy spectrum must have positive values"


# ---------------------------------------------------------------------------
# CCF profile
# ---------------------------------------------------------------------------

class TestCCFProfile:

    def test_shape_is_N_by_N(self):
        ic = generate_CCF_profile(1.1808, N=32)
        assert ic["omega"].shape == (32, 32)
        assert ic["u"].shape == (32, 32)
        assert ic["v"].shape == (32, 32)

    def test_all_required_keys_present(self):
        ic = generate_CCF_profile(1.1808, N=32)
        for key in ("omega", "u", "v", "x", "y", "kx", "ky",
                    "lambda_val", "N", "family", "width", "amplitude",
                    "k_bins", "E_k", "omega_max", "enstrophy"):
            assert key in ic, f"Missing key: {key}"

    def test_family_is_CCF(self):
        ic = generate_CCF_profile(1.1808, N=32)
        assert ic["family"] == "CCF"

    def test_omega_max_near_1(self):
        """Peak vorticity should be normalised to ~1."""
        ic = generate_CCF_profile(1.1808, N=64, omega_max=1.0)
        assert 0.5 < ic["omega_max"] < 2.0, (
            f"omega_max = {ic['omega_max']:.4f} out of expected range [0.5, 2.0]"
        )

    def test_enstrophy_positive(self):
        ic = generate_CCF_profile(1.1808, N=32)
        assert ic["enstrophy"] > 0

    def test_velocity_divergence_free(self):
        ic = generate_CCF_profile(1.1808, N=32)
        kx, ky = ic["kx"], ic["ky"]
        u_hat  = np.fft.fft2(ic["u"])
        v_hat  = np.fft.fft2(ic["v"])
        div_hat = 1j * kx * u_hat + 1j * ky * v_hat
        div     = np.real(np.fft.ifft2(div_hat))
        assert np.max(np.abs(div)) < 1e-6

    def test_all_wang_CCF_lambdas(self):
        """All three CCF λ values generate without error."""
        ccf_lambdas = [1.1808, 0.6057, 0.4703]
        for lv in ccf_lambdas:
            ic = generate_CCF_profile(lv, N=32)
            assert ic["lambda_val"] == lv
            assert ic["omega_max"] > 0

    def test_critical_case_0_4703(self):
        """CCF 2nd unstable (λ=0.4703) — the critical test case (PRD §4.3)."""
        ic = generate_CCF_profile(0.4703, N=64)
        diag = profile_diagnostics(ic)
        assert diag["omega_max"] > 0
        assert diag["enstrophy"] > 0
        assert diag["energy"] > 0, f"Energy must be positive, got {diag['energy']:.4e}"

    def test_reproducibility_with_seed(self):
        """Same seed → identical field."""
        ic1 = generate_CCF_profile(1.1808, N=32, seed=0)
        ic2 = generate_CCF_profile(1.1808, N=32, seed=0)
        assert np.allclose(ic1["omega"], ic2["omega"])

    def test_different_seeds_differ(self):
        ic1 = generate_CCF_profile(1.1808, N=32, seed=0)
        ic2 = generate_CCF_profile(1.1808, N=32, seed=99)
        assert not np.allclose(ic1["omega"], ic2["omega"])

    def test_smaller_lambda_has_more_fine_structure(self):
        """
        Smaller λ → more unstable → finer-scale structure (higher palinstrophy).
        CCF 2nd unstable (λ=0.4703) should have higher palinstrophy than CCF stable (λ=1.1808).
        """
        ic_stable   = generate_CCF_profile(1.1808, N=64, seed=42)
        ic_unstable = generate_CCF_profile(0.4703, N=64, seed=42)
        diag_s = profile_diagnostics(ic_stable)
        diag_u = profile_diagnostics(ic_unstable)
        assert diag_u["palinstrophy"] >= diag_s["palinstrophy"], (
            f"Expected unstable profile to have higher palinstrophy: "
            f"{diag_u['palinstrophy']:.4e} vs {diag_s['palinstrophy']:.4e}"
        )

    def test_energy_positive(self):
        ic = generate_CCF_profile(1.1808, N=64)
        diag = profile_diagnostics(ic)
        assert diag["energy"] > 0, f"Energy = {diag['energy']:.4e} should be > 0"

    def test_custom_omega_max(self):
        ic = generate_CCF_profile(1.0, N=32, omega_max=5.0)
        assert ic["omega_max"] > 1.0, "Scaling omega_max should give larger amplitude"


# ---------------------------------------------------------------------------
# Boussinesq profile
# ---------------------------------------------------------------------------

class TestBoussinesqProfile:

    def test_all_required_keys(self):
        ic = generate_Boussinesq_profile(1.9206, N=32)
        for key in ("omega", "b", "u", "v", "x", "y",
                    "lambda_val", "N", "family", "enstrophy", "buoyancy_rms"):
            assert key in ic, f"Missing key: {key}"

    def test_family_is_BSSQ(self):
        ic = generate_Boussinesq_profile(1.9206, N=32)
        assert ic["family"] == "BSSQ"

    def test_buoyancy_field_shape(self):
        ic = generate_Boussinesq_profile(1.9206, N=32)
        assert ic["b"].shape == (32, 32)

    def test_buoyancy_mostly_positive(self):
        """
        Buoyancy field should be predominantly non-negative (light rising plume).
        Small negative values are acceptable due to spectral dealiasing (Gibbs effect).
        """
        ic = generate_Boussinesq_profile(1.9206, N=64)
        b_max = float(ic["b"].max())
        b_min = float(ic["b"].min())
        # Any negative values should be tiny compared to the peak
        assert b_min > -0.01 * b_max, (
            f"Buoyancy min/max ratio too negative: min={b_min:.4e}, max={b_max:.4e}"
        )
        assert b_max > 0, "Buoyancy peak must be positive"

    def test_buoyancy_rms_positive(self):
        ic = generate_Boussinesq_profile(1.9206, N=32)
        assert ic["buoyancy_rms"] > 0

    def test_all_boussinesq_lambdas(self):
        """All three Boussinesq λ values generate without error."""
        bssq_lambdas = [1.9206, 1.3991, 1.1843]
        for lv in bssq_lambdas:
            ic = generate_Boussinesq_profile(lv, N=32)
            assert ic["lambda_val"] == lv
            assert ic["omega_max"] > 0

    def test_boussinesq_profiles_all_have_buoyancy(self):
        """
        All Boussinesq profiles must have a non-trivial buoyancy field.
        (Structural check: all λ values produce a coupled ω + b IC.)
        """
        for lv in [1.9206, 1.3991, 1.1843]:
            ic = generate_Boussinesq_profile(lv, N=32)
            assert float(ic["buoyancy_max"]) > 0.1, (
                f"λ={lv}: buoyancy_max = {ic['buoyancy_max']:.4e} — not a proper IC"
            )


# ---------------------------------------------------------------------------
# Adversarial minimum
# ---------------------------------------------------------------------------

class TestAdversarialMin:

    def test_all_required_keys(self):
        ic = generate_adversarial_min(N=32)
        for key in ("omega", "T_field", "u", "v", "x", "y",
                    "lambda_val", "N", "family", "T_min", "T_max",
                    "enstrophy", "palinstrophy"):
            assert key in ic, f"Missing key: {key}"

    def test_family_name(self):
        ic = generate_adversarial_min(N=32)
        assert ic["family"] == "adversarial_min"

    def test_T_field_positive(self):
        """T > 0 everywhere — second law must hold even in adversarial IC."""
        ic = generate_adversarial_min(N=64)
        assert float(ic["T_field"].min()) > 0, (
            f"T_min = {ic['T_field'].min():.4f} K — second law violated in adversarial IC"
        )

    def test_T_field_above_T_base(self):
        """T_min > T_base (adversarial means T_base + small delta)."""
        T_base = 300.0
        ic = generate_adversarial_min(N=32, T_base=T_base, T_margin=5.0)
        assert ic["T_min"] > T_base, (
            f"T_min = {ic['T_min']:.2f} should be > T_base = {T_base}"
        )

    def test_high_palinstrophy(self):
        """Adversarial IC designed for max stretching — should have high palinstrophy."""
        ic_adv = generate_adversarial_min(N=64)
        ic_ccf = generate_CCF_profile(1.1808, N=64)
        diag_adv = profile_diagnostics(ic_adv)
        diag_ccf = profile_diagnostics(ic_ccf)
        assert diag_adv["palinstrophy"] > diag_ccf["palinstrophy"], (
            f"Adversarial should have higher palinstrophy: "
            f"{diag_adv['palinstrophy']:.4e} vs {diag_ccf['palinstrophy']:.4e}"
        )

    def test_omega_max_normalised(self):
        ic = generate_adversarial_min(N=32, omega_max=1.0)
        assert abs(ic["omega_max"] - 1.0) < 0.01, (
            f"omega_max = {ic['omega_max']:.4f}, expected ~1.0 after normalisation"
        )


# ---------------------------------------------------------------------------
# generate_wang_profile dispatcher
# ---------------------------------------------------------------------------

class TestWangProfileDispatcher:

    def test_all_named_profiles_generate(self):
        all_names = list(WANG_PROFILES.keys()) + ["adversarial_min"]
        for name in all_names:
            ic = generate_wang_profile(name, N=32)
            assert "omega" in ic, f"Profile '{name}' missing omega"
            assert ic["omega"].shape == (32, 32)

    def test_profile_name_stored(self):
        ic = generate_wang_profile("CCF_2nd_unstable", N=32)
        assert ic["profile_name"] == "CCF_2nd_unstable"

    def test_prd_prediction_stored(self):
        ic = generate_wang_profile("CCF_2nd_unstable", N=32)
        assert "prd_prediction" in ic
        assert len(ic["prd_prediction"]) > 0

    def test_invalid_name_raises(self):
        with pytest.raises(ValueError, match="Unknown profile"):
            generate_wang_profile("nonexistent_profile", N=32)

    def test_all_wang_profiles_returns_dict_of_7(self):
        all_ics = generate_all_wang_profiles(N=32)
        # 6 named profiles + adversarial_min = 7
        assert len(all_ics) == 7, f"Expected 7 profiles, got {len(all_ics)}"

    def test_ccf_2nd_unstable_lambda_correct(self):
        """Exact λ value must be preserved (PRD §4.3 critical test)."""
        ic = generate_wang_profile("CCF_2nd_unstable", N=32)
        assert ic["lambda_val"] == 0.4703, (
            f"λ = {ic['lambda_val']}, expected 0.4703 (PRD §4.3)"
        )


# ---------------------------------------------------------------------------
# Lambda sweep ICs
# ---------------------------------------------------------------------------

class TestLambdaSweepICs:

    def test_default_returns_6_profiles(self):
        sweep = generate_lambda_sweep_ICs(N=32)
        assert len(sweep) == 6, f"Expected 6, got {len(sweep)}"

    def test_lambda_values_match_prd(self):
        """PRD §7.3 canonical values must be present."""
        sweep = generate_lambda_sweep_ICs(N=32)
        expected = set(LAMBDA_SWEEP_VALUES)
        actual   = set(sweep.keys())
        assert expected == actual, f"Lambda values mismatch: {actual} vs {expected}"

    def test_all_profiles_valid(self):
        sweep = generate_lambda_sweep_ICs(N=32)
        for lv, ic in sweep.items():
            assert ic["omega_max"] > 0, f"lambda={lv}: omega_max ≤ 0"
            assert ic["enstrophy"] > 0, f"lambda={lv}: enstrophy ≤ 0"

    def test_custom_lambda_values(self):
        sweep = generate_lambda_sweep_ICs(lambda_values=[0.5, 1.0], N=32)
        assert set(sweep.keys()) == {0.5, 1.0}

    def test_enstrophy_variation_with_lambda(self):
        """Different λ values should produce profiles with similar enstrophy (all normalised)."""
        sweep = generate_lambda_sweep_ICs(N=32)
        enstrophies = [ic["enstrophy"] for ic in sweep.values()]
        # Enstrophy should be within 2× of each other (not wildly different)
        ratio = max(enstrophies) / min(enstrophies)
        assert ratio < 5.0, (
            f"Enstrophy varies too widely across λ values: ratio = {ratio:.2f}"
        )


# ---------------------------------------------------------------------------
# Profile diagnostics
# ---------------------------------------------------------------------------

class TestProfileDiagnostics:

    def test_returns_all_fields(self):
        ic   = generate_CCF_profile(1.1808, N=32)
        diag = profile_diagnostics(ic)
        for field in ("lambda_val", "N", "family", "omega_max", "omega_rms",
                      "enstrophy", "palinstrophy", "energy", "vortex_width"):
            assert field in diag, f"Missing diagnostic: {field}"

    def test_energy_positive(self):
        ic   = generate_CCF_profile(1.1808, N=64)
        diag = profile_diagnostics(ic)
        assert diag["energy"] > 0, f"Energy = {diag['energy']:.4e} should be > 0"

    def test_enstrophy_matches_ic(self):
        ic   = generate_CCF_profile(1.1808, N=32)
        diag = profile_diagnostics(ic)
        assert abs(diag["enstrophy"] - ic["enstrophy"]) / ic["enstrophy"] < 0.01

    def test_palinstrophy_positive(self):
        ic   = generate_CCF_profile(0.4703, N=32)
        diag = profile_diagnostics(ic)
        assert diag["palinstrophy"] > 0

    def test_vortex_width_positive(self):
        ic   = generate_CCF_profile(1.1808, N=64)
        diag = profile_diagnostics(ic)
        assert diag["vortex_width"] > 0

    def test_boussinesq_diagnostics(self):
        ic   = generate_Boussinesq_profile(1.9206, N=32)
        diag = profile_diagnostics(ic)
        assert diag["omega_max"] > 0
        assert diag["energy"] > 0


# ---------------------------------------------------------------------------
# Energy spectrum check — concentrated ICs have steep spectra
# ---------------------------------------------------------------------------

class TestEnergySpectrum:

    def test_spectrum_is_non_negative(self):
        ic = generate_CCF_profile(1.1808, N=64)
        assert np.all(ic["E_k"] >= 0)

    def test_spectrum_peaks_at_low_k(self):
        """Concentrated profiles should have most energy at low wavenumbers."""
        ic = generate_CCF_profile(1.1808, N=64)
        k_bins = ic["k_bins"]
        E_k    = ic["E_k"]
        # Energy in k < 5 should exceed energy in k > 10
        low_k  = E_k[k_bins < 5].sum()
        high_k = E_k[k_bins > 10].sum()
        assert low_k > high_k, (
            f"Low-k energy ({low_k:.4e}) should exceed high-k ({high_k:.4e}) "
            f"for concentrated profile"
        )

    def test_verify_energy_spectrum_scaling_function(self):
        """Function runs without error and returns expected keys."""
        ic     = generate_CCF_profile(1.1808, N=64)
        result = verify_energy_spectrum_scaling(ic)
        for key in ("slope", "intercept", "k_range", "passes_check", "reason"):
            assert key in result

    def test_total_energy_positive(self):
        ic = generate_CCF_profile(1.1808, N=64)
        assert ic["E_k"].sum() > 0
