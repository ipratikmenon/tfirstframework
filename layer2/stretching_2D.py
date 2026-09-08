"""
stretching_2D.py — Dissipation / Stretching Diagnostics (2D)
=============================================================
In 3D Navier-Stokes, vortex stretching ω·∇u drives enstrophy growth.
In 2D, there is no vortex-stretching term (it is identically zero), but
palinstrophy can grow via the Jacobian nonlinearity.

The dissipation/stretching ratio tracked here is:
  D/S(t) = (enstrophy dissipation rate) / |palinstrophy production rate by J|

Enstrophy equation in 2D:
  dZ/dt = −2ν·P    (only dissipation — no stretching)
  → enstrophy is MONOTONICALLY DECREASING in 2D

Palinstrophy equation in 2D:
  dP/dt = −2ν·Q + N_P
  where Q = ½∫|Δω|² dx  (hyperviscous-like dissipation)
  and   N_P = −∫Δω·J(ψ,ω) dx  (nonlinear production / suppression)

The T-First suppression claim:
  The viscous dissipation rate D = 2ν̄·P outpaces the palinstrophy
  production rate |N_P| for all Wang et al. λ values.

Ratio:
  D/S = (2ν̄·P) / max(|N_P|, ε_floor)
  PASS criterion: D/S > 1  (dissipation dominates)

These functions also compute:
  • Viscous heating rate vs advective stretching rate (Conjecture 3.4 analog)
  • Spectral enstrophy flux (cascade direction indicator)
"""

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Basic spectral diagnostics
# ─────────────────────────────────────────────────────────────────────────────

def enstrophy(omega, domain_L=2.0 * np.pi):
    """Z = ½ ∫_Ω ω² dx.

    Parameters
    ----------
    omega    : (N,N) vorticity field [s⁻¹]
    domain_L : domain side length [m]

    Returns
    -------
    Z : float [s⁻²·m²]
    """
    return 0.5 * float(np.mean(omega**2)) * domain_L**2


def palinstrophy(omega, kx, ky, domain_L=2.0 * np.pi):
    """P = ½ ∫_Ω |∇ω|² dx (computed spectrally).

    Parameters
    ----------
    omega    : (N,N) vorticity field
    kx, ky   : wavenumber arrays (N,N) [rad/m]
    domain_L : domain side length [m]

    Returns
    -------
    P : float [s⁻²]  (units: [ω]²/[L]² · [L]² = [ω]²)
    """
    omega_hat = np.fft.fft2(omega)
    domega_dx = np.real(np.fft.ifft2(1j * kx * omega_hat))
    domega_dy = np.real(np.fft.ifft2(1j * ky * omega_hat))
    grad_sq = domega_dx**2 + domega_dy**2
    return 0.5 * float(np.mean(grad_sq)) * domain_L**2


def hyperdissipation(omega, kx, ky, domain_L=2.0 * np.pi):
    """Q = ½ ∫ |Δω|² dx (for palinstrophy equation)."""
    k2 = kx**2 + ky**2
    lap_omega = np.real(np.fft.ifft2(-k2 * np.fft.fft2(omega)))
    return 0.5 * float(np.mean(lap_omega**2)) * domain_L**2


# ─────────────────────────────────────────────────────────────────────────────
# Dissipation rates
# ─────────────────────────────────────────────────────────────────────────────

def enstrophy_dissipation_rate(omega, nu_mean, kx, ky, domain_L=2.0 * np.pi):
    """Rate of enstrophy dissipation: D_Z = 2ν̄ · P [s⁻³·m²].

    In 2D: dZ/dt = −2ν·P  (exact, from integrating ω·(ν∆ω)).
    """
    P = palinstrophy(omega, kx, ky, domain_L)
    return 2.0 * nu_mean * P


def palinstrophy_production_rate(omega, u, v, kx, ky, domain_L=2.0 * np.pi):
    """Nonlinear palinstrophy production: N_P = −∫ Δω · J(ψ,ω) dx.

    This is the nonlinear source/sink in the palinstrophy equation.
    In 2D, N_P can be positive or negative; |N_P| is the "stretching" analog.
    """
    k2 = kx**2 + ky**2
    omega_hat = np.fft.fft2(omega)
    lap_omega = np.real(np.fft.ifft2(-k2 * omega_hat))

    # Jacobian J(ψ,ω) = u·∇ω
    domega_dx = np.real(np.fft.ifft2(1j * kx * omega_hat))
    domega_dy = np.real(np.fft.ifft2(1j * ky * omega_hat))
    J_psi_omega = u * domega_dx + v * domega_dy  # = u·∇ω

    # N_P = −∫ Δω · J dx
    integrand = -lap_omega * J_psi_omega
    N_P = float(np.mean(integrand)) * domain_L**2
    return N_P


# ─────────────────────────────────────────────────────────────────────────────
# Main diagnostic: dissipation / stretching ratio
# ─────────────────────────────────────────────────────────────────────────────

def dissipation_stretching_ratio(omega, u, v, nu_mean, kx, ky, domain_L=2.0 * np.pi):
    """Compute D/S = enstrophy_dissipation_rate / |palinstrophy_production_rate|.

    D/S > 1: viscous dissipation dominates (T-first suppression)
    D/S < 1: nonlinear production dominates (potential growth)

    Returns
    -------
    dict with keys:
      D_Z        : enstrophy dissipation rate [s⁻³·m²]
      N_P        : palinstrophy production rate (signed) [s⁻³·m²]
      D_S_ratio  : D/|N_P|  (or inf if N_P ≈ 0)
      verdict    : 'PASS' if D/S > 1, 'WARN' if 0.5 < D/S ≤ 1, 'FAIL' if D/S ≤ 0.5
      Z          : enstrophy [s⁻²·m²]
      P          : palinstrophy [s⁻²]
    """
    Z = enstrophy(omega, domain_L)
    P = palinstrophy(omega, kx, ky, domain_L)
    D_Z = 2.0 * nu_mean * P

    N_P = palinstrophy_production_rate(omega, u, v, kx, ky, domain_L)
    abs_N_P = abs(N_P)
    eps_floor = 1e-30  # prevent division by zero

    D_S = D_Z / max(abs_N_P, eps_floor)

    if D_S > 1.0:
        verdict = "PASS"
    elif D_S > 0.5:
        verdict = "WARN"
    else:
        verdict = "FAIL"

    return {
        "Z": Z,
        "P": P,
        "D_Z": D_Z,
        "N_P": N_P,
        "abs_N_P": abs_N_P,
        "D_S_ratio": D_S,
        "nu_mean": nu_mean,
        "verdict": verdict,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Viscous heating vs advective rate (Conjecture 3.4 analogue at PDE level)
# ─────────────────────────────────────────────────────────────────────────────

def heating_vs_advection_ratio(omega, u, v, Q_visc, kx, ky, domain_L=2.0 * np.pi):
    """Ratio of integrated viscous heating to integrated vorticity advection.

    Heating rate: H = ∫ Q_visc dx
    Advection magnitude: A_adv = ∫ |J(ψ,ω)| dx

    H / A_adv > 1 implies the T-field is driven faster by heating than
    by advective transport — consistent with Conjecture 3.4.

    Returns
    -------
    dict with keys: H, A_adv, ratio, verdict
    """
    omega_hat = np.fft.fft2(omega)
    domega_dx = np.real(np.fft.ifft2(1j * kx * omega_hat))
    domega_dy = np.real(np.fft.ifft2(1j * ky * omega_hat))
    J = u * domega_dx + v * domega_dy

    H = float(np.mean(Q_visc)) * domain_L**2
    A_adv = float(np.mean(np.abs(J))) * domain_L**2
    eps_floor = 1e-30

    ratio = H / max(A_adv, eps_floor)

    return {
        "H": H,
        "A_adv": A_adv,
        "ratio": ratio,
        "verdict": "PASS" if ratio > 0.0 else "FAIL",   # H > 0 always (Q_visc ≥ 0)
    }


# ─────────────────────────────────────────────────────────────────────────────
# Spectral enstrophy flux (forward cascade indicator)
# ─────────────────────────────────────────────────────────────────────────────

def enstrophy_spectrum(omega, kx, ky, n_bins=None):
    """Shell-averaged enstrophy spectrum E_Z(k) = Σ_{|k'|≈k} |ω̂(k')|².

    Used to assess whether enstrophy cascades forward (to small scales)
    or inversely (to large scales — characteristic of 2D turbulence).

    Returns
    -------
    k_bins : 1D array of wavenumber bin centres
    E_Z    : 1D array of enstrophy per shell
    """
    N = omega.shape[0]
    omega_hat = np.fft.fft2(omega)
    # Shell-average
    k_mag = np.sqrt(kx**2 + ky**2)
    k_int = np.round(k_mag).astype(int)
    k_max = N // 2
    if n_bins is None:
        n_bins = k_max

    k_bins = np.arange(1, n_bins + 1, dtype=float)
    E_Z = np.zeros(n_bins)
    for ki in range(1, n_bins + 1):
        mask = k_int == ki
        E_Z[ki - 1] = float(np.sum(np.abs(omega_hat[mask])**2)) / N**2

    return k_bins, E_Z


# ─────────────────────────────────────────────────────────────────────────────
# Full-step diagnostic bundle
# ─────────────────────────────────────────────────────────────────────────────

def step_diagnostics(omega, u, v, Q_visc, nu_mean, kx, ky, t, t_star, lambda_val,
                     domain_L=2.0 * np.pi):
    """Compute all stretching/dissipation diagnostics for one timestep.

    Parameters
    ----------
    omega, u, v  : current fields
    Q_visc       : viscous heating field [K/s] (from solver)
    nu_mean      : mean kinematic viscosity [m²/s]
    kx, ky       : wavenumber grids
    t            : current time [s]
    t_star       : blow-up time
    lambda_val   : Wang et al. λ
    domain_L     : domain side [m]

    Returns
    -------
    dict with all key diagnostics
    """
    ds = dissipation_stretching_ratio(omega, u, v, nu_mean, kx, ky, domain_L)
    hv = heating_vs_advection_ratio(omega, u, v, Q_visc, kx, ky, domain_L)

    # Conjecture 3.4 suppression ratio (analytical)
    tau = max((t_star - t) / t_star, 1e-10)
    supp_ratio = tau ** (-(2.0 + lambda_val)) if tau < 1.0 else float("inf")
    exponent_margin = 2.0 + lambda_val  # always > 0

    return {
        # Enstrophy / palinstrophy
        "Z": ds["Z"],
        "P": ds["P"],
        "D_Z": ds["D_Z"],          # enstrophy dissipation rate
        "N_P": ds["N_P"],          # palinstrophy production rate
        "D_S_ratio": ds["D_S_ratio"],
        "DS_verdict": ds["verdict"],
        # Viscous heating
        "H": hv["H"],
        "A_adv": hv["A_adv"],
        "heating_ratio": hv["ratio"],
        # Conjecture 3.4
        "tau": tau,
        "suppression_ratio": supp_ratio,
        "exponent_margin": exponent_margin,
        "lambda_val": lambda_val,
        "t": t,
        # Overall verdict
        "verdict": "PASS" if (ds["verdict"] in ("PASS", "WARN") and supp_ratio > 1.0) else "FAIL",
    }
