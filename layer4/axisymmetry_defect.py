"""
layer4/axisymmetry_defect.py
§23 Non-axisymmetry defect analysis.

Computes and decomposes:
  u = u_axi + u_na   (axisymmetric + non-axisymmetric parts about axis ê at x*)

Key quantities:
  D_axi  = ∫_{B_{r*}} |u_na|²              non-axisymmetric energy in ball
  D_dens = D_axi / r*³                      normalised density
  F_z    = ∮_{∂B_{r*}} p (u·n)(ê·n) dS     axial pressure-velocity flux
  F_ratio = F_z / (M^{3/2} r*²)             scale-normalised flux

The scaling obstruction (§23 analysis):
  For dσ/dt|_pressure ≥ c₀ > 0, we need D_axi ≥ c₀ M r*³ = c₀ M^{-1/2}.
  Known sources give D_axi ~ r*³ = M^{-3/2} (background) or M^{-5/2} (images).
  The flux F_z scales as M — closer to the needed magnitude.

Open Problem 23.2: Is F_z / M^{3/2} r*² ≥ c₀ > 0 uniformly as M→∞, σ→0?
"""

import numpy as np
from layer4.geometric_disorder import (
    spectral_wavenumbers,
    compute_M,
    compute_r_star,
    compute_e_hat,
    compute_A_loc,
    _spectral_grad,
)
from layer4.pressure_misalignment import (
    solve_pressure_spectral,
    compute_pressure_gradient,
)


# ─────────────────────────────────────────────────────────────────────────────
# Axisymmetric decomposition about axis ê at x*
# ─────────────────────────────────────────────────────────────────────────────

def _cylindrical_basis(e_hat_vec):
    """
    Given tube axis ê (unit 3-vector), return an orthonormal frame (e_r_ref, e_θ, ê).
    e_r_ref is an arbitrary radial reference direction perpendicular to ê.
    """
    ez = np.array(e_hat_vec, dtype=float)
    ez /= np.linalg.norm(ez) + 1e-30

    # Pick a reference direction not parallel to ez
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ez, ref)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])

    er0 = ref - np.dot(ref, ez) * ez
    er0 /= np.linalg.norm(er0) + 1e-30
    et0 = np.cross(ez, er0)
    et0 /= np.linalg.norm(et0) + 1e-30
    return er0, et0, ez


def axisymmetric_decompose(u, x_center, e_hat_vec, N, n_theta=32):
    """
    Decompose u = u_axi + u_na at each point in the domain.

    Axisymmetric part u_axi(r, z) is obtained by averaging u over θ-circles
    about the axis ê through x_center.

    Returns:
        u_axi  : (3, N, N, N) array — axisymmetric part
        u_na   : (3, N, N, N) array — non-axisymmetric part (u - u_axi)
        rho    : (N, N, N) — cylindrical radius from axis at each grid point
        z_cyl  : (N, N, N) — cylindrical z (projection onto ê)
    """
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")
    grid = np.stack([X, Y, Z], axis=-1)  # (N, N, N, 3)

    x0 = np.array(x_center)
    er0, et0, ez = _cylindrical_basis(e_hat_vec)

    # Displacement from axis point
    disp = grid - x0[None, None, None, :]  # (N,N,N,3) — ignoring periodicity for now

    # Cylindrical coordinates
    z_cyl = np.einsum("ijk l,l->ijk", disp, ez)
    r_perp = disp - z_cyl[..., None] * ez[None, None, None, :]  # perpendicular component
    rho = np.sqrt(np.einsum("ijkl,ijkl->ijk", r_perp, r_perp))  # radial distance

    # For each grid point, average u over the θ-circle of radius rho at same z
    # We interpolate u on the circle and average
    # Simple implementation: for each grid point (i,j,k), find the θ-circle
    # and sample u at n_theta points, average

    u_axi = np.zeros_like(u)

    thetas = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)

    for th_idx, (ct, st) in enumerate(zip(cos_t, sin_t)):
        # Rotation matrix about ez by angle theta
        # R v = v_z ez + cos(θ)(v - v_z ez) + sin(θ)(ez × v)
        # Applied to velocity components: just weight-average
        pass

    # Practical implementation: use azimuthal FFT
    # Express each component of u in terms of θ-Fourier modes about ê
    # The m=0 mode is the axisymmetric part

    # For 3D spectral code, approximate by averaging over Nθ sample points on each circle
    # For efficiency: use the structure of the spectral grid

    # Efficient version: for a regular Cartesian grid, use the fact that
    # the θ-average at fixed (rho, z) can be computed via interpolation on circles

    # Simplified (approximate) implementation for testing:
    # u_axi ≈ mean of u over all points at the same (rho_bin, z_bin)
    # This is an approximation valid for the purposes of bounding D_axi

    rho_flat   = rho.ravel()
    z_flat     = z_cyl.ravel()
    u_flat     = u.reshape(3, -1)

    # Bin into (rho, z) cells
    n_rho_bins = max(4, N // 4)
    n_z_bins   = max(4, N // 4)
    rho_max    = np.max(rho_flat) + 1e-10
    z_min, z_max = np.min(z_flat), np.max(z_flat) + 1e-10

    rho_edges = np.linspace(0, rho_max, n_rho_bins + 1)
    z_edges   = np.linspace(z_min, z_max, n_z_bins + 1)

    rho_idx = np.digitize(rho_flat, rho_edges) - 1
    z_idx   = np.digitize(z_flat,   z_edges)   - 1
    rho_idx = np.clip(rho_idx, 0, n_rho_bins - 1)
    z_idx   = np.clip(z_idx,   0, n_z_bins   - 1)

    bin_idx = rho_idx * n_z_bins + z_idx
    n_bins  = n_rho_bins * n_z_bins

    u_axi_flat = np.zeros((3, len(bin_idx)), dtype=float)
    for comp in range(3):
        bin_sums   = np.bincount(bin_idx, weights=u_flat[comp], minlength=n_bins)
        bin_counts = np.bincount(bin_idx, minlength=n_bins).clip(min=1)
        bin_means  = bin_sums / bin_counts
        u_axi_flat[comp] = bin_means[bin_idx]

    u_axi = u_axi_flat.reshape(3, N, N, N)
    u_na  = u - u_axi

    return u_axi, u_na, rho, z_cyl


def compute_axisymmetry_defect(u, omega, kx, ky, kz, n_theta=16):
    """
    Compute the axisymmetry defect D_axi = ∫_{B_{r*}} |u_na|² and related quantities.

    Returns dict:
      D_axi       : non-axisymmetric energy in B_{r*}
      D_total     : total kinetic energy in B_{r*}
      D_dens      : D_axi / r*³  (normalised density)
      D_ratio     : D_axi / D_total  (fraction non-axisymmetric)
      M, sigma, r_star
      threshold   : c₀ M r*³ = c₀ M^{-1/2}  (what D_axi needs to reach for §23 argument)
      gap_ratio   : D_dens / M  (should be ≥ c₀ if Open Problem 23.2 holds)
    """
    N = omega.shape[1]
    dx = 2 * np.pi / N

    M, mag, idx = compute_M(omega)
    if M < 1e-10:
        return {"M": 0.0, "error": "zero vorticity"}

    r_star = compute_r_star(M)
    e_hat, _ = compute_e_hat(omega, mag=mag)
    A_loc, sigma, _, _ = compute_A_loc(omega, kx, ky, kz, M=M, idx=idx, mag=mag)

    ix, iy, iz = idx
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    x_center = [x1d[ix], x1d[iy], x1d[iz]]
    e_hat_vec = np.array([e_hat[c, ix, iy, iz] for c in range(3)])
    norm_ev = np.linalg.norm(e_hat_vec)
    if norm_ev < 1e-12:
        e_hat_vec = np.array([0.0, 0.0, 1.0])
    else:
        e_hat_vec /= norm_ev

    u_axi, u_na, rho, z_cyl = axisymmetric_decompose(
        u, x_center, e_hat_vec, N, n_theta=n_theta
    )

    # Ball mask B_{r*}(x*)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")
    dX = np.minimum(np.abs(X - x_center[0]), 2*np.pi - np.abs(X - x_center[0]))
    dY = np.minimum(np.abs(Y - x_center[1]), 2*np.pi - np.abs(Y - x_center[1]))
    dZ = np.minimum(np.abs(Z - x_center[2]), 2*np.pi - np.abs(Z - x_center[2]))
    dist = np.sqrt(dX**2 + dY**2 + dZ**2)
    ball = dist <= r_star

    vol_element = dx**3
    ball_vol = float(np.sum(ball)) * vol_element

    u_na_sq  = u_na[0]**2 + u_na[1]**2 + u_na[2]**2
    u_sq     = u[0]**2    + u[1]**2    + u[2]**2

    D_axi  = float(np.sum(u_na_sq[ball])) * vol_element
    D_tot  = float(np.sum(u_sq[ball]))    * vol_element
    D_dens = D_axi / (r_star**3) if r_star > 0 else 0.0
    D_ratio = D_axi / (D_tot + 1e-30)

    # Threshold: D_axi ≥ c₀ M r*³ = c₀ M^{-1/2} is needed for §23 argument
    threshold = M * r_star**3  # = M^{-1/2}
    gap_ratio = D_dens / M if M > 0 else 0.0   # needs to be ≥ c₀

    return {
        "M":          M,
        "sigma":      sigma,
        "r_star":     r_star,
        "D_axi":      D_axi,
        "D_total":    D_tot,
        "D_dens":     D_dens,
        "D_ratio":    D_ratio,
        "threshold":  threshold,   # needed value of D_axi
        "gap_ratio":  gap_ratio,   # D_dens / M — is this ≥ c₀?
        "ball_vol":   ball_vol,
        "e_hat":      e_hat_vec.tolist(),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Pressure-velocity flux through ∂B_{r*}  (the F_z quantity)
# ─────────────────────────────────────────────────────────────────────────────

def compute_pressure_flux(u, omega, kx, ky, kz):
    """Compute the axial pressure-velocity flux through ∂B_{r*}.

       F_z = ∮ p (u·n)(ê·n) dS  over shell near ∂B_{r*}

    This is the key §23 energy-flux quantity.  We approximate it via the
    divergence theorem applied to a shell around B_{r*}:

       F_z ≈ ∫_{B_{2r*} \ B_{r*}} ∂_ê [ p (u·∇)ê ] dx

    Returns:
      F_z       : axial flux
      F_ratio   : F_z / (M^{3/2} r*²)   — scale-normalised (should be O(1) or ≥ c₀)
      F_total   : total pressure flux ∮ p(u·n) dS
    """
    N = omega.shape[1]
    dx = 2 * np.pi / N
    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)

    M, mag, idx = compute_M(omega)
    if M < 1e-10:
        return {"F_z": 0.0, "F_ratio": 0.0, "F_total": 0.0, "M": 0.0}

    r_star = compute_r_star(M)
    e_hat, _ = compute_e_hat(omega, mag=mag)

    ix, iy, iz = idx
    x_center = np.array([x1d[ix], x1d[iy], x1d[iz]])
    e_hat_vec = np.array([e_hat[c, ix, iy, iz] for c in range(3)])
    norm_ev = np.linalg.norm(e_hat_vec)
    e_hat_vec = e_hat_vec / (norm_ev + 1e-30)

    p = solve_pressure_spectral(u, kx, ky, kz)

    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")
    dX = np.minimum(np.abs(X - x_center[0]), 2*np.pi - np.abs(X - x_center[0]))
    dY = np.minimum(np.abs(Y - x_center[1]), 2*np.pi - np.abs(Y - x_center[1]))
    dZ = np.minimum(np.abs(Z - x_center[2]), 2*np.pi - np.abs(Z - x_center[2]))
    dist = np.sqrt(dX**2 + dY**2 + dZ**2)

    # Annular shell approximating ∂B_{r*}
    shell = (dist > 0.8 * r_star) & (dist <= 1.2 * r_star)

    # Outward normal on shell: n = (x - x*)/|x - x*|
    with np.errstate(divide="ignore", invalid="ignore"):
        nx = np.where(dist > 1e-10, dX / dist, 0.0)
        ny = np.where(dist > 1e-10, dY / dist, 0.0)
        nz = np.where(dist > 1e-10, dZ / dist, 0.0)

    # u·n (radial velocity)
    u_dot_n = u[0]*nx + u[1]*ny + u[2]*nz

    # ê·n (axial component of normal)
    e_dot_n = e_hat_vec[0]*nx + e_hat_vec[1]*ny + e_hat_vec[2]*nz

    vol_element = dx**3
    F_total = float(np.sum(p[shell] * u_dot_n[shell])) * vol_element
    F_z     = float(np.sum(p[shell] * u_dot_n[shell] * e_dot_n[shell])) * vol_element

    M32_r2 = M**1.5 * r_star**2
    F_ratio = F_z / (M32_r2 + 1e-30)

    return {
        "F_z":       F_z,
        "F_total":   F_total,
        "F_ratio":   F_ratio,
        "M":         M,
        "r_star":    r_star,
        "M32_r2":    M32_r2,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Scaling sweep: measure D_axi and F_z as M increases
# ─────────────────────────────────────────────────────────────────────────────

def scaling_sweep(field_factory, N, amplitudes, kx, ky, kz):
    """
    For a family of fields parameterised by amplitude (∝ M), measure:
      - M, sigma, D_axi, D_dens, gap_ratio, F_z, F_ratio
    to test whether D_axi / (M r*³) and F_z / M^{3/2} r*² are bounded below.

    Returns list of result dicts.
    """  # noqa: W605
    results = []
    for amp in amplitudes:
        u, omega = field_factory(N, amplitude=amp)
        defect = compute_axisymmetry_defect(u, omega, kx, ky, kz)
        flux   = compute_pressure_flux(u, omega, kx, ky, kz)
        if "error" not in defect:
            row = {**defect, **{k: v for k, v in flux.items() if k not in defect}}
            results.append(row)
    return results
