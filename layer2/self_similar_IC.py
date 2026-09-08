"""
self_similar_IC.py — Self-Similar Initial Condition Generator (Layer 2, M2)

Generates 2D vorticity / buoyancy fields whose structure mirrors the unstable
self-similar blow-up profiles discovered by Wang, Bennani, Martens et al.
(Google DeepMind / NYU / Brown / Stanford, arXiv:2509.14185, September 2025).

These are the adversarial initial conditions used in the lambda sweep (M3).
The T-first 2D and 3D solvers are then started from these ICs to test whether
temperature-dependent viscosity μ(T) suppresses the blow-up.

PRD §7.2 Layer 2 (self_similar_IC.py), §4.3 (Planned experiments), §7.3 (Lambda sweep).

Mathematical background
-----------------------
Wang et al. discover self-similar blow-up profiles parameterised by λ > 0.
Concentration rate: |ω(·,t)|_∞ ~ (t* - t)^{-(1+λ)}
Spatial scale:      L(t) ~ (t* - t)^{1/(1+λ)}  →  0 as t → t*

In 2D self-similar coordinates ξ = (x - x*) / L(t):
    ω(x,t) = (t* - t)^{-(1+λ)} · Ω(ξ)

The self-similar profile Ω(ξ) is the invariant shape of the blow-up.
At t = 0 with t* = 1 (normalised): the IC is just the profile Ω(x) rescaled.

Two profile families:
  CCF   — Constantin-Cao-Foias type: concentrated vortex tube / ring
  BSSQ  — Boussinesq: coupled vorticity + buoyancy (density anomaly)

Smaller λ = faster blow-up = more unstable = harder test for T-first.
CCF 2nd unstable (λ = 0.4703) is the critical test case (PRD §4.3).

All output is on a periodic square domain [0, 2π]^2 or [-π, π]^2.
"""

from __future__ import annotations
import numpy as np
from typing import Optional, Tuple, Dict, Any
import sys, os

# Allow running from any directory
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

# ---------------------------------------------------------------------------
# Wang et al. λ values (PRD §4.3) — do not change these
# ---------------------------------------------------------------------------

WANG_PROFILES: Dict[str, Dict[str, Any]] = {
    "CCF_stable": {
        "lambda_val": 1.1808, "family": "CCF",
        "description": "CCF stable profile — baseline test (easy suppression)",
        "prd_prediction": "Thermal feedback suppresses — baseline",
    },
    "CCF_1st_unstable": {
        "lambda_val": 0.6057, "family": "CCF",
        "description": "CCF 1st unstable profile — first sharp test",
        "prd_prediction": "Viscous heating rate > blowup rate",
    },
    "CCF_2nd_unstable": {
        "lambda_val": 0.4703, "family": "CCF",
        "description": "CCF 2nd unstable — critical test case (PRD §4.3)",
        "prd_prediction": "Stronger test — heating still dominates?",
    },
    "Boussinesq_stable": {
        "lambda_val": 1.9206, "family": "BSSQ",
        "description": "Boussinesq stable — high λ, T-first easily wins",
        "prd_prediction": "High lambda — T-first easily wins",
    },
    "Boussinesq_1st": {
        "lambda_val": 1.3991, "family": "BSSQ",
        "description": "Boussinesq 1st unstable — intermediate test",
        "prd_prediction": "Intermediate — key test case",
    },
    "Boussinesq_3rd": {
        "lambda_val": 1.1843, "family": "BSSQ",
        "description": "Boussinesq 3rd unstable — approaching λ = 1 limit",
        "prd_prediction": "Approaching lambda=1 limit",
    },
}

# PRD §7.3 lambda sweep values
LAMBDA_SWEEP_VALUES = [1.2, 0.9, 0.6057, 0.4703, 0.3, 0.1]


# ---------------------------------------------------------------------------
# Grid utilities
# ---------------------------------------------------------------------------

def make_grid(N: int, domain: float = 2 * np.pi) -> Tuple[np.ndarray, np.ndarray,
                                                            np.ndarray, np.ndarray]:
    """
    Create an N×N periodic spectral grid on [0, domain)².

    Returns:
        x, y   — 2D coordinate arrays (N×N)
        kx, ky — wavenumber arrays (N×N), arranged for np.fft.fft2
    """
    dx = domain / N
    x1d = np.arange(N) * dx
    x, y = np.meshgrid(x1d, x1d, indexing="ij")

    # Wavenumbers for FFT output ordering
    k1d = np.fft.fftfreq(N, d=dx / (2 * np.pi))  # cycles per unit
    kx, ky = np.meshgrid(k1d, k1d, indexing="ij")

    return x, y, kx, ky


def spectral_laplacian_inv(f: np.ndarray, kx: np.ndarray, ky: np.ndarray,
                            zero_mean: bool = True) -> np.ndarray:
    """
    Compute ψ = (-Δ)^{-1} f   via spectral inversion.
    Used to obtain streamfunction from vorticity: Δψ = -ω → ψ = -(-Δ)^{-1} ω.
    """
    f_hat = np.fft.fft2(f)
    k2 = kx ** 2 + ky ** 2
    k2[0, 0] = 1.0  # avoid division by zero; DC mode set below
    psi_hat = f_hat / k2
    if zero_mean:
        psi_hat[0, 0] = 0.0
    return np.real(np.fft.ifft2(psi_hat))


def vorticity_to_velocity(omega: np.ndarray, kx: np.ndarray,
                           ky: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Recover (u, v) velocity from vorticity ω via streamfunction ψ.
    Δψ = -ω   →   u = ∂ψ/∂y,  v = -∂ψ/∂x
    """
    psi = spectral_laplacian_inv(-omega, kx, ky)
    psi_hat = np.fft.fft2(psi)
    u = np.real(np.fft.ifft2(1j * ky * psi_hat))
    v = np.real(np.fft.ifft2(-1j * kx * psi_hat))
    return u, v


def dealias_23(f: np.ndarray) -> np.ndarray:
    """
    Apply 2/3 dealiasing rule: zero out wavenumbers with |k| > N/3.

    In FFT ordering for even N:
      Keep rows 0..cutoff and rows (N-cutoff)..N-1 (wavenumbers 0..cutoff and -cutoff..-1).
      Zero rows (cutoff+1)..(N-cutoff-1).
    This preserves Hermitian symmetry so ifft2 returns a real field.
    """
    N = f.shape[0]
    f_hat = np.fft.fft2(f)
    cutoff = N // 3  # keep |k| <= cutoff, zero |k| > cutoff
    # Zero rows/cols corresponding to |k| > cutoff (the middle band in FFT ordering)
    f_hat[cutoff + 1: N - cutoff, :] = 0.0
    f_hat[:, cutoff + 1: N - cutoff] = 0.0
    return np.real(np.fft.ifft2(f_hat))


def energy_spectrum(omega: np.ndarray, kx: np.ndarray,
                    ky: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the 1D energy spectrum E(k) from a 2D vorticity field.

    E(k) = Σ_{|k_vec|≈k}  |ψ̂(k_vec)|² / 2

    Returns:
        k_bins  — wavenumber bins (positive integers)
        E_k     — energy in each bin
    """
    N = omega.shape[0]
    k_mag = np.sqrt(kx ** 2 + ky ** 2)
    psi = spectral_laplacian_inv(-omega, kx, ky)
    psi_hat = np.fft.fft2(psi) / N ** 2

    k_max = int(np.floor(np.max(k_mag)))
    k_bins = np.arange(1, k_max + 1, dtype=float)
    E_k = np.zeros_like(k_bins)

    for i, k in enumerate(k_bins):
        mask = (k_mag >= k - 0.5) & (k_mag < k + 0.5)
        E_k[i] = 0.5 * np.sum(np.abs(psi_hat[mask]) ** 2)

    return k_bins, E_k


# ---------------------------------------------------------------------------
# Self-similar profile shapes
# ---------------------------------------------------------------------------

def _gaussian_vortex_profile(x: np.ndarray, y: np.ndarray,
                              cx: float, cy: float,
                              amplitude: float, width: float,
                              n_poles: int = 1) -> np.ndarray:
    """
    Multi-pole Gaussian vortex — building block for CCF and Boussinesq profiles.

    For n_poles = 1: simple Gaussian vortex (positive / negative blob)
    For n_poles = 2: dipole (vortex pair — typical CCF structure)
    For n_poles = 4: quadrupole (more localised energy transfer)
    """
    N = x.shape[0]
    omega = np.zeros_like(x)

    angles = np.linspace(0, 2 * np.pi * (1 - 1 / n_poles), n_poles)
    signs  = np.array([(-1) ** i for i in range(n_poles)])
    offset = width * 1.5  # separation between poles

    for angle, sign in zip(angles, signs):
        cx_i = cx + offset * np.cos(angle)
        cy_i = cy + offset * np.sin(angle)
        r2 = (x - cx_i) ** 2 + (y - cy_i) ** 2
        omega += sign * amplitude * np.exp(-r2 / (2 * width ** 2))

    return omega


def _self_similar_width(lambda_val: float, t_norm: float = 0.0,
                         t_star: float = 1.0) -> float:
    """
    Self-similar spatial scale: L(t) ~ (t_star - t)^{1/(1+λ)}.
    At t = 0: L = t_star^{1/(1+λ)} (normalised to domain fraction).

    Smaller λ → slower width decay exponent but the amplitude grows faster.
    We normalise so that width = 0.5 at λ=1 and adjust with λ.
    """
    tau = t_star - t_norm
    exponent = 1.0 / (1.0 + lambda_val)
    # Normalised: multiply by base width (0.4 of domain/(2π)) so profiles are visible
    return 0.4 * (tau ** exponent)


def _self_similar_amplitude(lambda_val: float, t_norm: float = 0.0,
                              t_star: float = 1.0,
                              omega_max: float = 1.0) -> float:
    """
    Self-similar amplitude: A(t) ~ (t_star - t)^{-(1+λ)}.
    At t = 0: A = t_star^{-(1+λ)}.

    We normalise so that the peak vorticity equals omega_max at t=0.
    """
    tau = t_star - t_norm
    return omega_max * tau ** (-(1.0 + lambda_val)) / (1.0 ** (-(1.0 + lambda_val)))


# ---------------------------------------------------------------------------
# CCF profiles
# ---------------------------------------------------------------------------

def generate_CCF_profile(lambda_val: float, N: int = 128,
                          omega_max: float = 1.0,
                          seed: int = 42) -> Dict[str, np.ndarray]:
    """
    Generate a CCF-type self-similar initial condition.

    Constantin-Cao-Foias equations support vortex-tube blow-up solutions.
    Our model: a concentrated vortex dipole whose amplitude and width
    respect the self-similar scaling with parameter λ.

    Args:
        lambda_val: Self-similar scaling parameter (λ > 0).
        N: Grid size (N×N).
        omega_max: Peak vorticity amplitude (normalised to 1).
        seed: Random seed for reproducibility.

    Returns:
        dict with keys:
            omega    — vorticity field (N×N)
            u, v     — velocity components (N×N each)
            x, y     — coordinate grids (N×N)
            kx, ky   — wavenumber grids (N×N)
            lambda_val, N, family, width, amplitude
            energy_k, k_bins  — 1D energy spectrum
            omega_max, enstrophy, palinstrophy
    """
    rng = np.random.default_rng(seed)
    x, y, kx, ky = make_grid(N)
    L = 2 * np.pi

    width     = _self_similar_width(lambda_val)
    amplitude = omega_max  # normalised at t=0

    # Centre of the concentration (slightly off-centre for asymmetry)
    cx = L / 2 + rng.uniform(-0.1, 0.1) * L * 0.05
    cy = L / 2 + rng.uniform(-0.1, 0.1) * L * 0.05

    # CCF structure: dipole vortex (like Lamb-Chaplygin dipole)
    n_poles = 2
    omega = _gaussian_vortex_profile(x, y, cx, cy, amplitude, width, n_poles=n_poles)

    # Add subdominant higher-mode perturbation (captures unstable eigenmode structure)
    # Amplitude scaled by 0.1 to remain a perturbation, width halved for finer structure
    if lambda_val < 1.0:
        # Smaller λ = more unstable = stronger fine-structure perturbation
        epsilon = 0.15 * (1.0 - lambda_val)
        omega += _gaussian_vortex_profile(
            x, y, cx + width, cy, epsilon * amplitude, width / 2, n_poles=4
        )

    # Dealias
    omega = dealias_23(omega)

    # Recover velocity
    u, v = vorticity_to_velocity(omega, kx, ky)

    # Diagnostics
    enstrophy    = 0.5 * np.mean(omega ** 2) * L ** 2
    palinstrophy = 0.5 * np.mean((np.gradient(omega, L / N)[0] ** 2 +
                                   np.gradient(omega, L / N)[1] ** 2)) * L ** 2
    k_bins, E_k = energy_spectrum(omega, kx, ky)

    return {
        "omega"       : omega,
        "u"           : u,
        "v"           : v,
        "x"           : x,
        "y"           : y,
        "kx"          : kx,
        "ky"          : ky,
        "lambda_val"  : lambda_val,
        "N"           : N,
        "family"      : "CCF",
        "width"       : width,
        "amplitude"   : amplitude,
        "k_bins"      : k_bins,
        "E_k"         : E_k,
        "omega_max"   : float(np.max(np.abs(omega))),
        "enstrophy"   : enstrophy,
        "palinstrophy": palinstrophy,
    }


# ---------------------------------------------------------------------------
# Boussinesq profiles
# ---------------------------------------------------------------------------

def generate_Boussinesq_profile(lambda_val: float, N: int = 128,
                                 omega_max: float = 1.0,
                                 buoyancy_max: float = 1.0,
                                 seed: int = 42) -> Dict[str, np.ndarray]:
    """
    Generate a Boussinesq self-similar initial condition.

    Boussinesq equations couple vorticity ω with a buoyancy field b (density anomaly):
        ∂ω/∂t + (u·∇)ω = ∂b/∂x  +  ν Δω
        ∂b/∂t + (u·∇)b = κ Δb

    Wang et al. find blow-up for Boussinesq at λ = 1.9206, 1.3991, 1.1843.
    Structure: a rising buoyancy plume with a vortex dipole at the plume head.

    Args:
        lambda_val: Self-similar scaling parameter (λ > 0).
        N: Grid size.
        omega_max: Peak vorticity.
        buoyancy_max: Peak buoyancy perturbation.
        seed: Random seed.

    Returns:
        dict with keys: omega, b (buoyancy), u, v, x, y, kx, ky,
                        lambda_val, N, family, width, amplitude,
                        k_bins, E_k, omega_max, buoyancy_rms, enstrophy
    """
    rng = np.random.default_rng(seed)
    x, y, kx, ky = make_grid(N)
    L = 2 * np.pi

    width     = _self_similar_width(lambda_val)

    cx = L / 2
    cy = L / 3  # plume starts in lower third

    # Vorticity: dipole at plume head (baroclinic generation structure)
    omega = _gaussian_vortex_profile(x, y, cx, cy, omega_max, width, n_poles=2)

    # Buoyancy: Gaussian blob centred slightly below vortex head
    # (lighter fluid rising, creating baroclinic vorticity at the interface)
    by = cy - width * 0.5
    r2_b = (x - cx) ** 2 + (y - by) ** 2
    b = buoyancy_max * np.exp(-r2_b / (2 * (width * 1.2) ** 2))

    # For smaller λ (more unstable): sharper buoyancy gradient
    if lambda_val < 1.5:
        sharpness = 1.0 + (1.5 - lambda_val)
        b = buoyancy_max * np.exp(-r2_b * sharpness / (2 * (width * 1.2) ** 2))

    # Dealias both fields
    omega = dealias_23(omega)
    b     = dealias_23(b)

    # Recover velocity from vorticity
    u, v = vorticity_to_velocity(omega, kx, ky)

    # Diagnostics
    enstrophy    = 0.5 * np.mean(omega ** 2) * L ** 2
    buoyancy_rms = np.sqrt(np.mean(b ** 2))
    k_bins, E_k  = energy_spectrum(omega, kx, ky)

    return {
        "omega"        : omega,
        "b"            : b,
        "u"            : u,
        "v"            : v,
        "x"            : x,
        "y"            : y,
        "kx"           : kx,
        "ky"           : ky,
        "lambda_val"   : lambda_val,
        "N"            : N,
        "family"       : "BSSQ",
        "width"        : width,
        "amplitude"    : omega_max,
        "k_bins"       : k_bins,
        "E_k"          : E_k,
        "omega_max"    : float(np.max(np.abs(omega))),
        "buoyancy_max" : float(np.max(np.abs(b))),
        "buoyancy_rms" : buoyancy_rms,
        "enstrophy"    : enstrophy,
    }


# ---------------------------------------------------------------------------
# Adversarial minimum (λ → 0)
# ---------------------------------------------------------------------------

def generate_adversarial_min(lambda_val: float = 0.05, N: int = 128,
                              omega_max: float = 1.0,
                              T_base: float = 300.0,
                              T_margin: float = 5.0,
                              seed: int = 42) -> Dict[str, np.ndarray]:
    """
    Construct an adversarial initial condition designed to be hardest for T-first.

    PRD §7.5: adversarial IC protocol:
      - Biot-Savart optimisation for maximum vortex stretching alignment
      - T0 = T_min + δ everywhere (weakest possible thermal dissipation)
      - Critical vortex tube geometry from Beale-Kato-Majda analysis

    Here we construct the 2D analogue:
      - Vortex lattice with maximum positive stretching alignment
      - Multiple concentrated vortex sheets folded for cascading energy transfer
      - Initial temperature field T = T_base + T_margin (barely above T_min)

    Args:
        lambda_val: Adversarial λ (small, default 0.05).
        N: Grid size.
        omega_max: Peak vorticity.
        T_base: Background temperature (K).
        T_margin: ΔT above T_min (weakens thermal restoring force).
        seed: Random seed.

    Returns:
        dict with all ICfields plus T_field (temperature initial condition).
    """
    rng = np.random.default_rng(seed)
    x, y, kx, ky = make_grid(N)
    L = 2 * np.pi

    # --- Vortex lattice: 4 counter-rotating vortices in a square ---
    # Maximises vortex stretching (antiparallel alignment)
    offsets = [(L/3, L/3), (2*L/3, L/3), (L/3, 2*L/3), (2*L/3, 2*L/3)]
    signs   = [+1, -1, -1, +1]  # checkerboard — maximises strain
    width   = 0.25 * (1.0 + lambda_val)  # very concentrated for small λ

    omega = np.zeros((N, N))
    for (cx, cy), sign in zip(offsets, signs):
        r2 = (x - cx) ** 2 + (y - cy) ** 2
        omega += sign * omega_max * np.exp(-r2 / (2 * width ** 2))

    # Add thin vortex sheet (maximum enstrophy for given energy)
    # Horizontal shear layer — maximises palinstrophy
    stripe_y = L / 2
    stripe_w = width / 3
    r2_sheet = (y - stripe_y) ** 2
    omega += 0.3 * omega_max * np.sin(4 * x / L * 2 * np.pi) * \
             np.exp(-r2_sheet / (2 * stripe_w ** 2))

    # Small-scale noise in upper third of spectrum (broadband blow-up trigger)
    N_third = N // 3
    noise_hat = np.zeros((N, N), dtype=complex)
    noise_hat[1:N_third, 1:N_third] = (
        rng.normal(0, 0.05, (N_third - 1, N_third - 1)) +
        1j * rng.normal(0, 0.05, (N_third - 1, N_third - 1))
    )
    omega += np.real(np.fft.ifft2(noise_hat)) * omega_max

    # Dealias
    omega = dealias_23(omega)
    omega = omega * (omega_max / max(np.max(np.abs(omega)), 1e-10))  # renorm

    # --- Temperature field: minimal thermal dissipation ---
    # T = T_base + T_margin + small fluctuation (cold spots weaken restoring force)
    T_fluct = 0.5 * T_margin * np.abs(omega) / max(np.max(np.abs(omega)), 1e-10)
    T_field  = T_base + T_margin + T_fluct  # T > T_min everywhere (second law holds)

    # Velocity
    u, v = vorticity_to_velocity(omega, kx, ky)

    enstrophy    = 0.5 * np.mean(omega ** 2) * L ** 2
    palinstrophy = 0.5 * np.mean((np.gradient(omega, L / N)[0] ** 2 +
                                   np.gradient(omega, L / N)[1] ** 2)) * L ** 2
    k_bins, E_k  = energy_spectrum(omega, kx, ky)

    return {
        "omega"        : omega,
        "T_field"      : T_field,
        "u"            : u,
        "v"            : v,
        "x"            : x,
        "y"            : y,
        "kx"           : kx,
        "ky"           : ky,
        "lambda_val"   : lambda_val,
        "N"            : N,
        "family"       : "adversarial_min",
        "width"        : width,
        "amplitude"    : omega_max,
        "T_base"       : T_base,
        "T_margin"     : T_margin,
        "T_min"        : float(np.min(T_field)),
        "T_max"        : float(np.max(T_field)),
        "k_bins"       : k_bins,
        "E_k"          : E_k,
        "omega_max"    : float(np.max(np.abs(omega))),
        "enstrophy"    : enstrophy,
        "palinstrophy" : palinstrophy,
    }


# ---------------------------------------------------------------------------
# Named Wang et al. profile dispatcher
# ---------------------------------------------------------------------------

def generate_wang_profile(profile_name: str, N: int = 128,
                           omega_max: float = 1.0,
                           seed: int = 42) -> Dict[str, np.ndarray]:
    """
    Generate a named Wang et al. profile by name.

    Args:
        profile_name: Key from WANG_PROFILES dict, or "adversarial_min".
        N: Grid size.
        omega_max: Peak vorticity.
        seed: Random seed.

    Returns:
        IC dict with an extra "profile_name" and "prd_prediction" field.
    """
    if profile_name == "adversarial_min":
        result = generate_adversarial_min(lambda_val=0.05, N=N,
                                          omega_max=omega_max, seed=seed)
        result["profile_name"]  = "adversarial_min"
        result["prd_prediction"] = "Does T-first break down?"
        return result

    if profile_name not in WANG_PROFILES:
        valid = list(WANG_PROFILES.keys()) + ["adversarial_min"]
        raise ValueError(
            f"Unknown profile '{profile_name}'. Valid: {valid}"
        )

    info = WANG_PROFILES[profile_name]
    lv   = info["lambda_val"]

    if info["family"] == "CCF":
        result = generate_CCF_profile(lv, N=N, omega_max=omega_max, seed=seed)
    else:
        result = generate_Boussinesq_profile(lv, N=N, omega_max=omega_max, seed=seed)

    result["profile_name"]  = profile_name
    result["prd_prediction"] = info["prd_prediction"]
    result["description"]   = info["description"]
    return result


def generate_all_wang_profiles(N: int = 128, omega_max: float = 1.0,
                                seed: int = 42) -> Dict[str, Dict]:
    """
    Generate all Wang et al. profiles + adversarial minimum.

    Returns a dict keyed by profile name.
    """
    names = list(WANG_PROFILES.keys()) + ["adversarial_min"]
    return {name: generate_wang_profile(name, N=N, omega_max=omega_max, seed=seed)
            for name in names}


# ---------------------------------------------------------------------------
# Lambda sweep ICs (PRD §7.3)
# ---------------------------------------------------------------------------

def generate_lambda_sweep_ICs(lambda_values: Optional[list] = None,
                               N: int = 64,
                               omega_max: float = 1.0,
                               seed: int = 42) -> Dict[float, Dict]:
    """
    Generate CCF-type initial conditions for the lambda sweep (PRD §7.3).

    Each IC has the same topology but different λ-dependent concentration.
    Used by the M3 lambda sweep to test Conjecture 3.4 numerically.

    Args:
        lambda_values: List of λ values. Defaults to PRD §7.3 canonical sweep.
        N: Grid size (64 default — fast sweep).
        omega_max: Peak vorticity.
        seed: Random seed.

    Returns:
        Dict keyed by λ value.
    """
    if lambda_values is None:
        lambda_values = LAMBDA_SWEEP_VALUES

    results = {}
    for i, lv in enumerate(lambda_values):
        ic = generate_CCF_profile(lv, N=N, omega_max=omega_max, seed=seed + i)
        ic["profile_name"] = f"lambda_sweep_{lv:.4f}"
        results[lv] = ic
    return results


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def profile_diagnostics(ic: Dict) -> Dict[str, float]:
    """
    Compute key diagnostics for a generated IC.

    Returns a summary dict suitable for logging to results.db.
    """
    omega = ic["omega"]
    N     = ic["N"]
    lv    = ic["lambda_val"]
    L     = 2 * np.pi
    dx    = L / N

    omega_max     = float(np.max(np.abs(omega)))
    omega_rms     = float(np.sqrt(np.mean(omega ** 2)))
    enstrophy     = 0.5 * np.mean(omega ** 2) * L ** 2

    grad_x = np.gradient(omega, dx)[0]
    grad_y = np.gradient(omega, dx)[1]
    palinstrophy  = 0.5 * np.mean(grad_x ** 2 + grad_y ** 2) * L ** 2

    # Width of the vorticity peak (standard deviation of vorticity distribution)
    x, y = ic["x"], ic["y"]
    norm = max(np.sum(np.abs(omega)), 1e-10)
    x_mean = np.sum(x * np.abs(omega)) / norm
    y_mean = np.sum(y * np.abs(omega)) / norm
    sigma_x = np.sqrt(np.sum((x - x_mean) ** 2 * np.abs(omega)) / norm)
    sigma_y = np.sqrt(np.sum((y - y_mean) ** 2 * np.abs(omega)) / norm)
    vortex_width = float(0.5 * (sigma_x + sigma_y))

    # Energy: E = 0.5 * ∫ |u|² dx   from streamfunction
    kx, ky = ic["kx"], ic["ky"]
    psi = spectral_laplacian_inv(-omega, kx, ky)
    # E = 0.5 * ∫ |∇ψ|² dx  (always ≥ 0; use spectral gradient)
    psi_hat = np.fft.fft2(psi) / N ** 2
    u_spec  = np.real(np.fft.ifft2(1j * ky * psi_hat * N ** 2))
    v_spec  = np.real(np.fft.ifft2(-1j * kx * psi_hat * N ** 2))
    energy  = 0.5 * np.mean(u_spec ** 2 + v_spec ** 2) * L ** 2

    return {
        "lambda_val"  : lv,
        "N"           : N,
        "family"      : ic.get("family", "unknown"),
        "omega_max"   : omega_max,
        "omega_rms"   : omega_rms,
        "enstrophy"   : enstrophy,
        "palinstrophy": palinstrophy,
        "energy"      : energy,
        "vortex_width": vortex_width,
        "width_input" : ic.get("width", float("nan")),
    }


def verify_energy_spectrum_scaling(ic: Dict, expected_slope: float = -3.0,
                                    k_range: Tuple[int, int] = (3, 20),
                                    tolerance: float = 0.5) -> Dict:
    """
    Check that the energy spectrum E(k) ~ k^{expected_slope} in the inertial range.

    For 2D turbulence: E(k) ~ k^{-3} (enstrophy cascade).
    For concentrated self-similar profiles: spectrum is steeper.

    Returns:
        dict with {slope, intercept, k_range, passes_check}
    """
    k_bins = ic["k_bins"]
    E_k    = ic["E_k"]

    k1, k2 = k_range
    mask = (k_bins >= k1) & (k_bins <= k2) & (E_k > 0)
    if mask.sum() < 3:
        return {"slope": float("nan"), "passes_check": False, "reason": "too few points"}

    log_k = np.log(k_bins[mask])
    log_E = np.log(E_k[mask])
    slope, intercept = np.polyfit(log_k, log_E, 1)

    passes = abs(slope - expected_slope) < abs(tolerance)
    return {
        "slope"       : float(slope),
        "intercept"   : float(intercept),
        "k_range"     : k_range,
        "expected_slope": expected_slope,
        "passes_check": passes,
        "reason"      : "ok" if passes else f"slope={slope:.2f} outside [{expected_slope-tolerance:.1f}, {expected_slope+tolerance:.1f}]",
    }


# ---------------------------------------------------------------------------
# Visualisation (optional — matplotlib)
# ---------------------------------------------------------------------------

def plot_profile(ic: Dict, save_path: Optional[str] = None,
                 show: bool = False) -> None:
    """
    Generate a 2-panel figure: vorticity field + energy spectrum.

    Args:
        ic: IC dict from any generator function.
        save_path: File path to save PNG. None = do not save.
        show: Whether to call plt.show() (False for batch runs).
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError:
        print("matplotlib not available — skipping visualisation")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: vorticity field
    omega = ic["omega"]
    vmax  = np.max(np.abs(omega))
    im = axes[0].pcolormesh(
        ic["x"], ic["y"], omega,
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto"
    )
    plt.colorbar(im, ax=axes[0], label="ω (vorticity)")
    axes[0].set_title(
        f"{ic.get('family','?')} profile  λ={ic['lambda_val']:.4f}\n"
        f"N={ic['N']}  ‖ω‖_∞={ic['omega_max']:.3f}  "
        f"Z={ic.get('enstrophy', float('nan')):.3e}"
    )
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    # Panel 2: energy spectrum
    k_bins = ic["k_bins"]
    E_k    = ic["E_k"]
    mask   = E_k > 0
    axes[1].loglog(k_bins[mask], E_k[mask], "b-", lw=1.5, label="E(k)")

    # Reference slopes
    k_ref = k_bins[mask][len(k_bins[mask]) // 4]
    E_ref = E_k[mask][len(k_bins[mask]) // 4]
    axes[1].loglog(k_bins[mask], E_ref * (k_bins[mask] / k_ref) ** (-3),
                   "k--", alpha=0.5, label="k⁻³")
    axes[1].loglog(k_bins[mask], E_ref * (k_bins[mask] / k_ref) ** (-5.0 / 3),
                   "r--", alpha=0.5, label="k^{-5/3}")

    axes[1].set_xlabel("wavenumber k")
    axes[1].set_ylabel("E(k)")
    axes[1].set_title(f"Energy spectrum  λ={ic['lambda_val']:.4f}")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


# ---------------------------------------------------------------------------
# CLI / smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import os

    print("=" * 70)
    print("T-First Self-Similar IC Generator — Layer 2 Smoke Test")
    print("=" * 70)

    os.makedirs("results/figures", exist_ok=True)

    print("\n--- Generating all Wang et al. profiles (N=64) ---")
    all_profiles = generate_all_wang_profiles(N=64)

    for name, ic in all_profiles.items():
        diag = profile_diagnostics(ic)
        spec = verify_energy_spectrum_scaling(ic)
        print(f"\n  {name}:")
        print(f"    λ = {ic['lambda_val']:.4f}  family = {ic['family']}")
        print(f"    ‖ω‖_∞ = {diag['omega_max']:.4f}  Z = {diag['enstrophy']:.4e}  "
              f"E = {diag['energy']:.4e}")
        print(f"    width = {diag['vortex_width']:.4f}  (input width = {diag['width_input']:.4f})")
        print(f"    spectrum slope = {spec['slope']:.2f}  [{spec['reason']}]")

        fig_path = f"results/figures/{name}_N64.png"
        plot_profile(ic, save_path=fig_path, show=False)

    print("\n--- Lambda sweep ICs (PRD §7.3, N=64) ---")
    sweep = generate_lambda_sweep_ICs(N=64)
    print(f"  λ values: {sorted(sweep.keys())}")
    for lv, ic in sorted(sweep.items()):
        diag = profile_diagnostics(ic)
        print(f"  λ={lv:.4f}: ‖ω‖_∞={diag['omega_max']:.4f}  "
              f"Z={diag['enstrophy']:.4e}  width={diag['vortex_width']:.4f}")

    print("\n=== IC Generator smoke test complete ===")
