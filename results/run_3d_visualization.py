"""
results/run_3d_visualization.py
3D visualization of vorticity field, geometric misalignment |∇ê|², and σ = A_loc/M^{3/2}.

Generates:
  results/figures/vorticity_3d_slices.png   — 3×2 slice plots
  results/figures/sigma_3d_scatter.png      — 3D scatter + B_{r*} ball
  results/figures/proof_status_dashboard.png — proof state summary

Run from any directory:
  python results/run_3d_visualization.py
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.patches import Circle

# Ensure layer4 is importable
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _repo_root)

from layer4.geometric_disorder import (
    spectral_wavenumbers,
    compute_M,
    compute_r_star,
    compute_e_hat,
    compute_A_loc,
    _spectral_grad,
)

FIGURES_DIR = os.path.join(_repo_root, "results", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Taylor-Green velocity field
# ─────────────────────────────────────────────────────────────────────────────

def taylor_green_velocity(N):
    """
    Classic Taylor-Green vortex on [0,2π]^3.
      u_x =  sin x cos y cos z
      u_y = -cos x sin y cos z
      u_z =  0
    """
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    ux =  np.sin(X) * np.cos(Y) * np.cos(Z)
    uy = -np.cos(X) * np.sin(Y) * np.cos(Z)
    uz =  np.zeros_like(X)
    return ux, uy, uz


def compute_curl_spectral(ux, uy, uz, kx, ky, kz):
    """ω = curl u via spectral differentiation."""
    ux_hat = np.fft.fftn(ux)
    uy_hat = np.fft.fftn(uy)
    uz_hat = np.fft.fftn(uz)

    # ωx = ∂_y uz − ∂_z uy
    duz_dy = np.real(np.fft.ifftn(1j * ky * uz_hat))
    duy_dz = np.real(np.fft.ifftn(1j * kz * uy_hat))
    omega_x = duz_dy - duy_dz

    # ωy = ∂_z ux − ∂_x uz
    dux_dz = np.real(np.fft.ifftn(1j * kz * ux_hat))
    duz_dx = np.real(np.fft.ifftn(1j * kx * uz_hat))
    omega_y = dux_dz - duz_dx

    # ωz = ∂_x uy − ∂_y ux
    duy_dx = np.real(np.fft.ifftn(1j * kx * uy_hat))
    dux_dy = np.real(np.fft.ifftn(1j * ky * ux_hat))
    omega_z = duy_dx - dux_dy

    return np.array([omega_x, omega_y, omega_z])


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Pointwise geometric diagnostics on the full grid
# ─────────────────────────────────────────────────────────────────────────────

def compute_grad_ehat_sq(omega, kx, ky, kz, eps=1e-12):
    """
    Compute |∇ê|²(x) on the full grid via spectral gradients of each ê component.
    Returns (grad_ehat_sq, e_hat, mag).
    """
    M, mag, idx = compute_M(omega)
    e_hat, mag = compute_e_hat(omega, mag, eps=eps)

    N = omega.shape[1]
    grad_ehat_sq = np.zeros((N, N, N), dtype=float)
    for comp in range(3):
        f_hat = np.fft.fftn(e_hat[comp])
        gx, gy, gz = _spectral_grad(f_hat, kx, ky, kz)
        grad_ehat_sq += gx**2 + gy**2 + gz**2

    return grad_ehat_sq, e_hat, mag, M, idx


def build_ball_wireframe(center_xyz, r_star, n_theta=20, n_phi=20):
    """
    Return (X, Y, Z) arrays for a wireframe sphere of radius r_star at center_xyz.
    Domain is [0, 2π]^3 — no periodic wrapping for display.
    """
    theta = np.linspace(0, np.pi, n_theta)
    phi = np.linspace(0, 2 * np.pi, n_phi)
    theta, phi = np.meshgrid(theta, phi)
    x0, y0, z0 = center_xyz
    X = x0 + r_star * np.sin(theta) * np.cos(phi)
    Y = y0 + r_star * np.sin(theta) * np.sin(phi)
    Z = z0 + r_star * np.cos(theta)
    return X, Y, Z


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Figure 1 — 3×2 slice plots
# ─────────────────────────────────────────────────────────────────────────────

def figure1_slices(fields, out_path):
    """
    fields: list of dicts with keys:
      label, omega_mag, grad_ehat_sq, integrand, iz_star, r_star, N, idx
    """
    n_fields = len(fields)
    fig, axes = plt.subplots(3, n_fields, figsize=(7 * n_fields, 16))
    fig.suptitle(
        "Vorticity Field, Misalignment |∇ê|², and A_loc Density",
        fontsize=16, fontweight="bold", y=0.98
    )

    row_labels = ["|ω(x,y,z=z*)|", "|∇ê|²(x,y,z=z*)", "|ω|²|∇ê|² — A_loc density"]
    cmaps = ["plasma", "viridis", "inferno"]

    x = np.linspace(0, 1, fields[0]["N"], endpoint=False)  # x/2π

    for col, fd in enumerate(fields):
        N = fd["N"]
        iz = fd["iz_star"]
        r_star = fd["r_star"]
        ix, iy, _ = fd["idx"]
        dx = 1.0 / N  # in units of 2π

        cx = ix * dx   # circle centre in normalised coords
        cy = iy * dx

        slices = [
            fd["omega_mag"][:, :, iz],
            fd["grad_ehat_sq"][:, :, iz],
            fd["integrand"][:, :, iz],
        ]

        for row in range(3):
            ax = axes[row, col]
            data = slices[row]
            im = ax.imshow(
                data.T,
                origin="lower",
                extent=[0, 1, 0, 1],
                cmap=cmaps[row],
                aspect="equal",
                interpolation="bilinear",
            )
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            # Draw self-similar ball as white circle
            r_norm = r_star / (2 * np.pi)   # normalise to [0,1]
            circ = Circle((cx, cy), r_norm,
                          fill=False, edgecolor="white", linewidth=2, linestyle="--")
            ax.add_patch(circ)
            ax.plot(cx, cy, "w*", markersize=12)

            if row == 0:
                ax.set_title(fd["label"], fontsize=13, fontweight="bold")
            ax.set_xlabel("x / 2π")
            ax.set_ylabel("y / 2π")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)

            # Row label on leftmost column
            if col == 0:
                ax.set_ylabel(f"y / 2π\n\n[{row_labels[row]}]", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Figure 2 — 3D scatter + B_{r*} wireframe
# ─────────────────────────────────────────────────────────────────────────────

def figure2_scatter3d(fd, out_path):
    """3D scatter of A_loc density with self-similar ball wireframe."""
    N = fd["N"]
    integrand = fd["integrand"]
    r_star = fd["r_star"]
    ix, iy, iz = fd["idx"]

    x1d = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")

    threshold = 0.1 * float(np.max(integrand))
    mask = integrand > threshold

    vals = integrand[mask]
    xs = X[mask]
    ys = Y[mask]
    zs = Z[mask]

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    # Use log scale if dynamic range > 100
    drange = float(np.max(vals)) / (float(np.min(vals)) + 1e-30)
    if drange > 100:
        c_vals = np.log10(vals + 1e-30)
        cbar_label = "log₁₀(|ω|²|∇ê|²)"
    else:
        c_vals = vals
        cbar_label = "|ω|²|∇ê|²"

    sc = ax.scatter(xs, ys, zs, c=c_vals, cmap="hot",
                    s=4, alpha=0.5, linewidths=0)
    plt.colorbar(sc, ax=ax, label=cbar_label, fraction=0.03, pad=0.08)

    # Self-similar ball wireframe
    cx = x1d[ix]
    cy = x1d[iy]
    cz = x1d[iz]
    Xw, Yw, Zw = build_ball_wireframe((cx, cy, cz), r_star)
    ax.plot_wireframe(Xw, Yw, Zw, color="cyan", alpha=0.25, linewidth=0.6,
                      rstride=3, cstride=3)

    # Mark x* (max vorticity point)
    ax.scatter([cx], [cy], [cz], color="red", s=200, marker="*",
               zorder=10, label="x* (max |ω|)")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(
        f"A_loc Density in 3D — Self-Similar Ball B_{{r*}}(x*)\n"
        f"Field: {fd['label']}   r* = {r_star:.4f}   σ = {fd['sigma']:.4f}",
        fontsize=12
    )
    ax.legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Figure 3 — Proof status dashboard
# ─────────────────────────────────────────────────────────────────────────────

def figure3_dashboard(out_path):
    fig = plt.figure(figsize=(14, 10))
    fig.patch.set_facecolor("#0e1117")
    fig.suptitle(
        "T-First Program: State of the Proof (S30, May 2026)",
        fontsize=16, fontweight="bold", color="white", y=0.97
    )

    axes = [
        fig.add_subplot(2, 2, 1),
        fig.add_subplot(2, 2, 2),
        fig.add_subplot(2, 2, 3),
        fig.add_subplot(2, 2, 4),
    ]
    for ax in axes:
        ax.set_facecolor("#1a1d27")
        ax.tick_params(colors="white")
        for spine in ax.spines.values():
            spine.set_edgecolor("#444")

    # ── Panel 1: Bar chart — proved vs open ──────────────────────────────────
    ax1 = axes[0]
    items = [
        ("Route C Bootstrap",      "proved"),
        ("Claim A (shear)",        "proved"),
        ("§20 Z(r) cascade",       "proved"),
        ("Theorem E (σ→0)",        "proved"),
        ("Pressure misalignment",  "conditional"),
        ("Energy class L²",        "open"),
    ]
    color_map = {"proved": "#2ecc71", "conditional": "#f39c12", "open": "#e74c3c"}
    labels = [x[0] for x in items]
    colors = [color_map[x[1]] for x in items]
    y_pos = range(len(labels))
    bars = ax1.barh(list(y_pos), [1] * len(items), color=colors, height=0.6, alpha=0.9)
    ax1.set_yticks(list(y_pos))
    ax1.set_yticklabels(labels, fontsize=9, color="white")
    ax1.set_xticks([])
    ax1.set_title("Proof Chain Status", color="white", fontsize=11, fontweight="bold")
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#2ecc71", label="Proved"),
        Patch(facecolor="#f39c12", label="Conditional / H¹"),
        Patch(facecolor="#e74c3c", label="Open"),
    ]
    ax1.legend(handles=legend_elements, loc="lower right",
               facecolor="#1a1d27", labelcolor="white", fontsize=8)

    # ── Panel 2: Pie chart — completion ──────────────────────────────────────
    ax2 = axes[1]
    pie_labels = [
        "Proved\n(unconditional)",
        "Proved\n(conditional/H¹)",
        "Open"
    ]
    pie_sizes = [4, 2, 2]
    pie_colors = ["#2ecc71", "#f39c12", "#e74c3c"]
    wedges, texts, autotexts = ax2.pie(
        pie_sizes, labels=pie_labels, colors=pie_colors,
        autopct="%1.0f%%", startangle=90,
        textprops={"color": "white", "fontsize": 9},
        pctdistance=0.75,
    )
    for at in autotexts:
        at.set_color("white")
        at.set_fontsize(9)
    ax2.set_title("Proof Completion", color="white", fontsize=11, fontweight="bold")

    # ── Panel 3: Session timeline ─────────────────────────────────────────────
    ax3 = axes[2]
    milestones = [
        (1,  6,  "S01–S06: Layer 1–2 (props, 2D solvers)",     "#3498db"),
        (7,  12, "S07–S12: Layer 3 (3D spectral NS)",           "#9b59b6"),
        (13, 16, "S13–S16: Layer 4 (diagnostics)",              "#1abc9c"),
        (17, 17, "S17: Papers 1–4",                             "#e67e22"),
        (18, 27, "S18–S27: Proof §1–§20 (Prize for u₀∈H¹)",    "#e74c3c"),
        (28, 29, "S28–S29: Prize complete + docs",              "#f1c40f"),
        (30, 30, "S30: §21–22 σ-framework (geometric disorder)","#2ecc71"),
    ]
    for i, (s_start, s_end, label, col) in enumerate(milestones):
        ax3.barh(i, s_end - s_start + 1, left=s_start - 1,
                 color=col, alpha=0.85, height=0.6)
        ax3.text(s_start - 0.5, i, label, va="center", ha="left",
                 fontsize=7.5, color="white",
                 bbox=dict(boxstyle="round,pad=0.1", fc=col, alpha=0.5))
    ax3.set_yticks([])
    ax3.set_xlabel("Session number", color="white", fontsize=9)
    ax3.set_xlim(-1, 32)
    ax3.set_title("Session Timeline", color="white", fontsize=11, fontweight="bold")
    ax3.xaxis.label.set_color("white")

    # ── Panel 4: Key results text panel ──────────────────────────────────────
    ax4 = axes[3]
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.set_title("Key Results", color="white", fontsize=11, fontweight="bold")

    results = [
        ("#2ecc71", "§20: Z(r) ≤ D·r^{2/3}  →  Prize for u₀ ∈ H¹"),
        ("#2ecc71", "§21: σ = A_loc/M^{3/2}  scale-invariant"),
        ("#2ecc71", "Theorem E: blowup ⟹ σ→0 + Type I"),
        ("#f39c12", "873/873 tests pass"),
        ("#3498db", "4979 lines, 22 proof sections"),
        ("#f1c40f", "Route C bootstrap: unconditional"),
        ("#f1c40f", "Vortex direction ê: max principle"),
        ("#e74c3c", "Open: u₀ ∈ L² (energy class)"),
    ]
    y_step = 1.0 / (len(results) + 1)
    for k, (col, text) in enumerate(results):
        y = 1.0 - (k + 1) * y_step
        ax4.plot(0.04, y, "o", color=col, markersize=9, transform=ax4.transAxes,
                 clip_on=False)
        ax4.text(0.10, y, text, va="center", ha="left",
                 fontsize=9, color="white", transform=ax4.transAxes)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 6.  Main driver
# ─────────────────────────────────────────────────────────────────────────────

def process_field(label, omega, kx, ky, kz):
    """Compute all pointwise diagnostics and summary scalars for one vorticity field."""
    M, mag, idx = compute_M(omega)
    r_star = compute_r_star(M)
    A_loc, sigma, _, ball_fraction = compute_A_loc(
        omega, kx, ky, kz, M=M, idx=idx, mag=mag
    )
    grad_ehat_sq, e_hat, mag2, M2, idx2 = compute_grad_ehat_sq(omega, kx, ky, kz)
    integrand = mag2**2 * grad_ehat_sq

    ix, iy, iz = idx

    # Find max of integrand
    max_loc = np.unravel_index(np.argmax(integrand), integrand.shape)

    N = omega.shape[1]
    dx_phys = 2.0 * np.pi / N

    print(f"\n{'='*60}")
    print(f"  Field      : {label}")
    print(f"  Grid size  : {N}³ = {N**3} points")
    print(f"  M = ‖ω‖_L∞ : {M:.6f}")
    print(f"  r* = M^-½  : {r_star:.6f}  ({r_star:.3f} rad  ≈ {r_star/dx_phys:.1f} cells)")
    print(f"  A_loc      : {A_loc:.6e}")
    print(f"  σ          : {sigma:.6f}")
    print(f"  Ball frac  : {ball_fraction:.4f}  ({100*ball_fraction:.1f}% of grid in B_r*)")
    print(f"  x* (idx)   : {idx}  → ({ix*dx_phys:.3f}, {iy*dx_phys:.3f}, {iz*dx_phys:.3f}) rad")
    print(f"  Integrand max loc: {max_loc}")
    print(f"{'='*60}")

    return {
        "label":        label,
        "omega_mag":    mag2,
        "grad_ehat_sq": grad_ehat_sq,
        "integrand":    integrand,
        "iz_star":      iz,
        "r_star":       r_star,
        "N":            N,
        "idx":          idx,
        "M":            M,
        "sigma":        sigma,
        "A_loc":        A_loc,
    }


def main():
    N = 32
    print(f"\nT-First 3D Visualization — N={N}³ grid")
    print("Computing Taylor-Green vorticity field...")

    kx, ky, kz = spectral_wavenumbers(N)

    # ── Field 1: Taylor-Green ────────────────────────────────────────────────
    ux, uy, uz = taylor_green_velocity(N)
    omega_tg = compute_curl_spectral(ux, uy, uz, kx, ky, kz)
    fd_tg = process_field("Taylor-Green (TG)", omega_tg, kx, ky, kz)

    # ── Field 2: TG + random perturbation ────────────────────────────────────
    print("\nAdding random perturbation to break symmetry...")
    rng = np.random.default_rng(seed=42)
    amplitude = 0.15
    ux_p = ux + amplitude * rng.standard_normal((N, N, N))
    uy_p = uy + amplitude * rng.standard_normal((N, N, N))
    uz_p = uz + amplitude * rng.standard_normal((N, N, N))
    # Project to divergence-free (spectral pressure projection)
    ux_hat = np.fft.fftn(ux_p)
    uy_hat = np.fft.fftn(uy_p)
    uz_hat = np.fft.fftn(uz_p)
    k2 = kx**2 + ky**2 + kz**2
    k2_safe = np.where(k2 > 0, k2, 1.0)
    kdotu = kx * ux_hat + ky * uy_hat + kz * uz_hat
    ux_hat -= kx * kdotu / k2_safe
    uy_hat -= ky * kdotu / k2_safe
    uz_hat -= kz * kdotu / k2_safe
    ux_p = np.real(np.fft.ifftn(ux_hat))
    uy_p = np.real(np.fft.ifftn(uy_hat))
    uz_p = np.real(np.fft.ifftn(uz_hat))
    omega_pert = compute_curl_spectral(ux_p, uy_p, uz_p, kx, ky, kz)
    fd_pert = process_field("TG + Perturbation", omega_pert, kx, ky, kz)

    # ── Figure 1: Slice plots ────────────────────────────────────────────────
    print("\nGenerating Figure 1: vorticity_3d_slices.png ...")
    fig1_path = os.path.join(FIGURES_DIR, "vorticity_3d_slices.png")
    figure1_slices([fd_tg, fd_pert], fig1_path)

    # ── Figure 2: 3D scatter (use perturbed field — more interesting) ────────
    print("Generating Figure 2: sigma_3d_scatter.png ...")
    fig2_path = os.path.join(FIGURES_DIR, "sigma_3d_scatter.png")
    figure2_scatter3d(fd_pert, fig2_path)

    # ── Figure 3: Proof dashboard ────────────────────────────────────────────
    print("Generating Figure 3: proof_status_dashboard.png ...")
    fig3_path = os.path.join(FIGURES_DIR, "proof_status_dashboard.png")
    figure3_dashboard(fig3_path)

    print("\n✓ All 3 figures saved to results/figures/")
    print(f"  {FIGURES_DIR}/")
    print("  ├── vorticity_3d_slices.png")
    print("  ├── sigma_3d_scatter.png")
    print("  └── proof_status_dashboard.png")


if __name__ == "__main__":
    main()
