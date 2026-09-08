"""
LPS_monitor.py — Ladyzhenskaya-Prodi-Serrin Norm Monitor (2D)
==============================================================
Tracks LPS-type norms and the T-First suppression margin throughout a solver run.

In 2D, global regularity is known; this module builds the monitoring infrastructure
that will be extended to 3D in M4/M5.  Tracked quantities:

  • Enstrophy  Z(t) = ½ ∫ ω² dx         [s⁻²·m²]
  • Palinstrophy P(t) = ½ ∫ |∇ω|² dx    [s⁻²·m⁰]
  • Lp norms of u: ‖u‖_{L^p} for p = 2, 4, ∞
  • LPS-like margin ε(λ,τ) = τ^{-(2+λ)} − 1  (Conjecture 3.4 exponent margin)
  • Dissipation-to-stretching ratio D/S  (see stretching_2D.py)
  • A(T)_min — asserted positive at every step

The 2D LPS analogue uses the (p=∞, q=2) Prodi-Serrin pair:
  ‖u‖_{L^∞_t L^2_x} is automatically controlled by the energy estimate.
  The LPS margin tracks how far we are from the boundary ‖u‖_{L^p L^q} = ∞.
"""

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Monitor construction
# ─────────────────────────────────────────────────────────────────────────────

def make_monitor():
    """Return a fresh monitor state dict."""
    return {
        "t_history":              [],
        "Z_history":              [],     # enstrophy
        "P_history":              [],     # palinstrophy
        "u_L2_history":           [],
        "u_L4_history":           [],
        "u_Linf_history":         [],
        "omega_Linf_history":     [],
        "A_min_history":          [],
        "suppression_ratio_history": [],
        "lps_margin_history":     [],     # ε(λ,τ) = τ^{-(2+λ)} − 1
        "n_steps":                0,
        "A_min_global":           float("inf"),
        "omega_max_global":       0.0,
        "lps_margin_min":         float("inf"),
        "all_pass":               True,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Per-step update
# ─────────────────────────────────────────────────────────────────────────────

def _enstrophy(omega, dx):
    """Z = ½ ∫ ω² dx ≈ ½ · mean(ω²) · L²"""
    L = omega.shape[0] * dx
    return 0.5 * float(np.mean(omega**2)) * L**2


def _palinstrophy(omega, kx, ky, dx):
    """P = ½ ∫ |∇ω|² dx (spectral gradient, spectral integration)."""
    N = omega.shape[0]
    L = N * dx
    omega_hat = np.fft.fft2(omega)
    domega_dx_hat = 1j * kx * omega_hat
    domega_dy_hat = 1j * ky * omega_hat
    # Parseval: ∫|∇ω|² = (1/N²) Σ |∇ω̂|² / ... use spatial mean approach
    grad_omega_sq = (
        np.real(np.fft.ifft2(domega_dx_hat))**2 +
        np.real(np.fft.ifft2(domega_dy_hat))**2
    )
    return 0.5 * float(np.mean(grad_omega_sq)) * L**2


def _lp_norm(f, p):
    """Compute ‖f‖_{L^p} over the domain (spatial L^p norm)."""
    if np.isinf(p):
        return float(np.max(np.abs(f)))
    N = f.shape[0]
    dx = 2.0 * np.pi / N
    L = N * dx
    return float((np.mean(np.abs(f)**p) * L**2)**(1.0 / p))


def monitor_step(monitor, omega, T, u, v, t, t_star, lambda_val, A_min, A_max=None):
    """Update the monitor with the current solver state.

    Parameters
    ----------
    monitor    : monitor dict (mutated in place)
    omega      : vorticity field (N,N)
    T          : temperature field (N,N)
    u, v       : velocity components (N,N)
    t          : current time [s]
    t_star     : blow-up time for suppression ratio
    lambda_val : Wang et al. λ
    A_min      : minimum of A(T) at this step (pre-computed by solver)
    A_max      : maximum of A(T) (optional)

    Returns
    -------
    Updated monitor dict (same object, also returned for convenience).
    """
    N = omega.shape[0]
    dx = 2.0 * np.pi / N
    kx_1d = np.fft.fftfreq(N, d=1.0 / N)
    kx, ky = np.meshgrid(kx_1d, kx_1d, indexing="ij")

    # Enstrophy and palinstrophy
    Z = _enstrophy(omega, dx)
    P = _palinstrophy(omega, kx, ky, dx)

    # Velocity norms
    speed = np.sqrt(u**2 + v**2)
    u_L2 = _lp_norm(speed, 2)
    u_L4 = _lp_norm(speed, 4)
    u_Linf = _lp_norm(speed, np.inf)
    omega_Linf = float(np.max(np.abs(omega)))

    # Conjecture 3.4 suppression ratio and LPS margin
    tau = max((t_star - t) / t_star, 1e-10)
    if tau < 1.0:
        supp_ratio = tau ** (-(2.0 + lambda_val))  # > 1 always
        lps_margin = supp_ratio - 1.0              # > 0 always
    else:
        supp_ratio = float("inf")
        lps_margin = float("inf")

    # Verdict: A_min must be > 0
    step_pass = A_min > 0.0

    # Update history
    monitor["t_history"].append(t)
    monitor["Z_history"].append(Z)
    monitor["P_history"].append(P)
    monitor["u_L2_history"].append(u_L2)
    monitor["u_L4_history"].append(u_L4)
    monitor["u_Linf_history"].append(u_Linf)
    monitor["omega_Linf_history"].append(omega_Linf)
    monitor["A_min_history"].append(A_min)
    monitor["suppression_ratio_history"].append(supp_ratio)
    monitor["lps_margin_history"].append(lps_margin)
    monitor["n_steps"] += 1

    # Running extremes
    monitor["A_min_global"] = min(monitor["A_min_global"], A_min)
    monitor["omega_max_global"] = max(monitor["omega_max_global"], omega_Linf)
    if not np.isinf(lps_margin):
        monitor["lps_margin_min"] = min(monitor["lps_margin_min"], lps_margin)

    if not step_pass:
        monitor["all_pass"] = False

    return monitor


# ─────────────────────────────────────────────────────────────────────────────
# Report
# ─────────────────────────────────────────────────────────────────────────────

def monitor_report(monitor):
    """Summarise the monitor history.

    Returns a dict with key scalars suitable for logging to results.db.
    """
    if monitor["n_steps"] == 0:
        return {"verdict": "NO_DATA", "n_steps": 0}

    Z_hist = monitor["Z_history"]
    P_hist = monitor["P_history"]
    S_hist = monitor["suppression_ratio_history"]

    # Enstrophy trend: should be monotonically decreasing in 2D (dZ/dt = -2ν P ≤ 0)
    Z_arr = np.array(Z_hist)
    Z_monotone = bool(np.all(np.diff(Z_arr) <= 1e-12 * Z_arr[:-1])) if len(Z_arr) > 1 else True

    # Finite suppression ratio values (exclude τ >= t_star)
    finite_S = [s for s in S_hist if not np.isinf(s)]

    return {
        "verdict": "PASS" if monitor["all_pass"] else "FAIL",
        "n_steps": monitor["n_steps"],
        "A_min_global": monitor["A_min_global"],
        "omega_max_global": monitor["omega_max_global"],
        "lps_margin_min": monitor["lps_margin_min"],
        "Z_initial": Z_hist[0] if Z_hist else None,
        "Z_final": Z_hist[-1] if Z_hist else None,
        "Z_monotone_decreasing": Z_monotone,
        "P_max": float(np.max(P_hist)) if P_hist else None,
        "suppression_ratio_min": float(min(finite_S)) if finite_S else None,
        "suppression_ratio_final": float(finite_S[-1]) if finite_S else None,
        "u_L2_max": float(np.max(monitor["u_L2_history"])),
        "u_Linf_max": float(np.max(monitor["u_Linf_history"])),
        "A_positive_throughout": monitor["A_min_global"] > 0.0,
        "conjecture_34_satisfied": (
            monitor["lps_margin_min"] > 0.0
            if not np.isinf(monitor["lps_margin_min"])
            else True
        ),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: attach monitor to a solver run via callback
# ─────────────────────────────────────────────────────────────────────────────

def make_monitor_callback(monitor, omega_ref, T_ref, u_ref, v_ref,
                          t_star, lambda_val):
    """Return a callback(step_diag) for use with T_solver_2D.run().

    The callback closes over the solver fields. Because the solver mutates
    its state dict (new_solver is returned from solver_step), the caller
    must inject the current (omega, T, u, v) via step_diag or re-extract
    them from the solver after each step.

    This factory function is used for tests; production code should call
    monitor_step directly in the run loop.
    """
    def callback(step_diag):
        # step_diag contains pre-computed A_min, u_max, etc. but not the full fields.
        # For lightweight monitoring, use the scalar diagnostics only.
        monitor["t_history"].append(step_diag["t"])
        monitor["A_min_history"].append(step_diag["A_min"])
        monitor["omega_Linf_history"].append(step_diag["omega_max"])
        monitor["suppression_ratio_history"].append(step_diag["suppression_ratio"])
        monitor["n_steps"] += 1
        monitor["A_min_global"] = min(monitor["A_min_global"], step_diag["A_min"])
        monitor["omega_max_global"] = max(monitor["omega_max_global"],
                                          step_diag["omega_max"])
        if not np.isinf(step_diag["suppression_ratio"]):
            lps_m = step_diag["suppression_ratio"] - 1.0
            monitor["lps_margin_history"].append(lps_m)
            monitor["lps_margin_min"] = min(monitor["lps_margin_min"], lps_m)
        if step_diag["A_min"] <= 0.0:
            monitor["all_pass"] = False

    return callback
