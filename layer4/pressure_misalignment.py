"""
layer4/pressure_misalignment.py
§23 pressure misalignment analysis.

Measures the axial (tube-axis) component of ∇p for nearly aligned vortex fields
to test Open Problem 23.1: does |∂_z p_glob / M^{3/2}|_{B_{r*}} ≥ c₀ > 0?

Key decomposition:
  ∇p = (∂_r p)e_r  +  (∂_z p)e_z    (in cylindrical coords about tube axis ê)
  ∂_r p ~ M^{3/2}   (centripetal — does NOT create misalignment)
  ∂_z p             (axial — creates S_{rz} strain → misalignment)
"""

import numpy as np
from layer4.geometric_disorder import (
    spectral_wavenumbers,
    compute_M,
    compute_e_hat,
    compute_A_loc,
    _spectral_grad,
)


# ─────────────────────────────────────────────────────────────────────────────
# Pressure solver: -Δp = ∂_i ∂_j (u_i u_j) = tr((∇u)²)
# ─────────────────────────────────────────────────────────────────────────────

def solve_pressure_spectral(u, kx, ky, kz):
    """
    Solve -Δp = ∂_i ∂_j (u_i u_j) spectrally on T³.

    RHS = -∂_i ∂_j (u_i u_j) = tr((∇u)^T ∇u) - div(div) correction
    Uses: -Δp = ∂_x² (u_x²) + ∂_y² (u_y²) + ∂_z² (u_z²)
                + 2∂_x∂_y(u_x u_y) + 2∂_x∂_z(u_x u_z) + 2∂_y∂_z(u_y u_z)

    Returns p on the grid (mean-zero).
    """
    N = u.shape[1]
    dx = 2 * np.pi / N

    # Build RHS = ∂_i ∂_j (u_i u_j) in spectral space
    rhs_hat = np.zeros((N, N, N), dtype=complex)

    # Diagonal terms
    for i, (ki, comp) in enumerate([(kx, 0), (ky, 1), (kz, 2)]):
        f_hat = np.fft.fftn(u[comp] ** 2) * dx**3
        rhs_hat += -(ki**2) * f_hat  # -∂²/∂x_i² → -k_i² in Fourier

    # Off-diagonal cross terms (factor 2 each)
    pairs = [(kx, ky, 0, 1), (kx, kz, 0, 2), (ky, kz, 1, 2)]
    for ki, kj, ci, cj in pairs:
        f_hat = np.fft.fftn(u[ci] * u[cj]) * dx**3
        rhs_hat += -2 * ki * kj * f_hat

    # Solve -Δp = rhs → p_hat = rhs_hat / (k² )  with p_hat[0] = 0
    k2 = kx**2 + ky**2 + kz**2
    with np.errstate(divide="ignore", invalid="ignore"):
        p_hat = np.where(k2 > 0, rhs_hat / k2, 0.0)

    p = np.real(np.fft.ifftn(p_hat)) / dx**3
    return p


def compute_pressure_gradient(p, kx, ky, kz):
    """Return ∇p as (3, N, N, N) array via spectral differentiation."""
    N = p.shape[0]
    dx = 2 * np.pi / N
    p_hat = np.fft.fftn(p) * dx**3
    gx = np.real(np.fft.ifftn(1j * kx * p_hat)) / dx**3
    gy = np.real(np.fft.ifftn(1j * ky * p_hat)) / dx**3
    gz = np.real(np.fft.ifftn(1j * kz * p_hat)) / dx**3
    return np.array([gx, gy, gz])


# ─────────────────────────────────────────────────────────────────────────────
# Axial / radial decomposition at x*
# ─────────────────────────────────────────────────────────────────────────────

def decompose_pressure_gradient(grad_p, e_hat, idx):
    """
    At location idx, decompose ∇p into axial (along ê) and radial (⊥ ê) parts.

    Returns:
        dp_axial   : scalar, component of ∇p along ê (tube axis) at x*
        dp_radial  : scalar, magnitude of ∇p ⊥ ê at x*
        dp_total   : scalar, |∇p| at x*
    """
    ix, iy, iz = idx
    grad_at_x = np.array([grad_p[c, ix, iy, iz] for c in range(3)])
    e_at_x    = np.array([e_hat[c, ix, iy, iz]  for c in range(3)])

    dp_axial_vec   = np.dot(grad_at_x, e_at_x) * e_at_x
    dp_radial_vec  = grad_at_x - dp_axial_vec

    dp_axial  = abs(np.dot(grad_at_x, e_at_x))
    dp_radial = np.linalg.norm(dp_radial_vec)
    dp_total  = np.linalg.norm(grad_at_x)
    return dp_axial, dp_radial, dp_total


# ─────────────────────────────────────────────────────────────────────────────
# Ball-averaged decomposition
# ─────────────────────────────────────────────────────────────────────────────

def ball_averaged_pressure_decomposition(grad_p, e_hat, mag, idx, r_star, N):
    """
    Average the axial/radial pressure gradient over B_{r*}(x*).

    Returns mean_axial, mean_radial, mean_total (all normalised per unit volume).
    """
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")

    ix, iy, iz = idx
    x0 = np.array([x1d[ix], x1d[iy], x1d[iz]])

    # Periodic distance
    dX = np.minimum(np.abs(X - x0[0]), 2 * np.pi - np.abs(X - x0[0]))
    dY = np.minimum(np.abs(Y - x0[1]), 2 * np.pi - np.abs(Y - x0[1]))
    dZ = np.minimum(np.abs(Z - x0[2]), 2 * np.pi - np.abs(Z - x0[2]))
    dist = np.sqrt(dX**2 + dY**2 + dZ**2)
    ball = dist <= r_star

    axials, radials, totals = [], [], []
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if ball[i, j, k] and mag[i, j, k] > 1e-8:
                    gv = np.array([grad_p[c, i, j, k] for c in range(3)])
                    ev = np.array([e_hat[c, i, j, k]  for c in range(3)])
                    norm_ev = np.linalg.norm(ev)
                    if norm_ev < 1e-12:
                        continue
                    ev /= norm_ev
                    ax = abs(np.dot(gv, ev))
                    rad = np.linalg.norm(gv - np.dot(gv, ev) * ev)
                    tot = np.linalg.norm(gv)
                    axials.append(ax)
                    radials.append(rad)
                    totals.append(tot)

    if not axials:
        return 0.0, 0.0, 0.0

    return np.mean(axials), np.mean(radials), np.mean(totals)


# ─────────────────────────────────────────────────────────────────────────────
# Main diagnostic: measure axial pressure ratio for various fields
# ─────────────────────────────────────────────────────────────────────────────

def pressure_misalignment_diagnostic(u, omega, kx, ky, kz):
    """
    Full §23 pressure diagnostic.

    Returns dict with:
      M, sigma, r_star,
      dp_axial_x*   : |∂_z p| at x* (tube axis direction)
      dp_radial_x*  : |∇p|_⊥ at x*
      dp_total_x*   : |∇p| at x*
      ratio_axial   : dp_axial / M^{3/2}  (should be ≥ c₀ > 0 if Open Problem holds)
      ratio_radial  : dp_radial / M^{3/2}
      ball_axial    : mean |∂_z p| over B_{r*}
      ball_ratio    : ball_axial / M^{3/2}
    """
    N = omega.shape[1]
    dx = 2 * np.pi / N

    M, mag, idx = compute_M(omega)
    if M < 1e-10:
        return {"M": 0.0, "error": "zero vorticity"}

    r_star = M ** (-0.5)
    e_hat, _ = compute_e_hat(omega, mag=mag)
    A_loc, sigma, _, _ = compute_A_loc(omega, kx, ky, kz, M=M, idx=idx, mag=mag)

    p = solve_pressure_spectral(u, kx, ky, kz)
    grad_p = compute_pressure_gradient(p, kx, ky, kz)

    dp_axial, dp_radial, dp_total = decompose_pressure_gradient(grad_p, e_hat, idx)
    ball_axial, ball_radial, ball_total = ball_averaged_pressure_decomposition(
        grad_p, e_hat, mag, idx, r_star, N
    )

    M32 = M ** 1.5

    return {
        "M":               M,
        "sigma":           sigma,
        "r_star":          r_star,
        "dp_axial_xstar":  dp_axial,
        "dp_radial_xstar": dp_radial,
        "dp_total_xstar":  dp_total,
        "ratio_axial":     dp_axial  / M32 if M32 > 0 else 0.0,
        "ratio_radial":    dp_radial / M32 if M32 > 0 else 0.0,
        "ratio_total":     dp_total  / M32 if M32 > 0 else 0.0,
        "ball_axial":      ball_axial,
        "ball_ratio":      ball_axial / M32 if M32 > 0 else 0.0,
        "ball_radial":     ball_radial,
    }
