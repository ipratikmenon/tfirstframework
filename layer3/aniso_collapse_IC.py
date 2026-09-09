"""
aniso_collapse_IC.py — Anisotropic Collapsing-Vortex Adversarial IC (Layer 3)
==============================================================================
Inspired by the geometry of OpenAI's forced Navier-Stokes blowup construction
("Finite Time Blowup for Navier-Stokes", Sep 2026, cdn.openai.com) — but this
module introduces NO forcing term. f = 0 everywhere; the equations run are
the exact, unforced Route 2 Prize equations (route2_3D_jax.run_jax), same as
every other experiment in this layer. Nothing here diverts from the Prize
equations: this is purely a new *initial condition* shape.

What we borrowed and what we didn't
------------------------------------
OpenAI's construction (Theorem 1.1) needs an external smooth force f to
sustain an axisymmetric vortex core collapsing at DIFFERENT rates in the
radial and axial directions (Re_theta -> infinity while Re_r = O(1)). Their
force is structurally load-bearing — it seeds the amplifying pulses, cancels
the base flow's singular residual, and supplies the compact-support cutoff.
None of that is reproducible or desired here.

What we DO borrow is the qualitative shape: an axisymmetric swirl-plus-
axial-outflow vortex column whose radial core scale ℓr and axial core scale
ℓz can be set independently (aspect_ratio = ℓz/ℓr). We use this as a new
*adversarial initial condition family*, complementary to the Wang et al.
self-similar profiles (layer2/self_similar_IC.py), and ask the same question
M6 already asks of those profiles: does δ > 0 (Lemma 2.5) — or, run under
Route 1, does A(T) > 0 — survive contact with this geometry under the real,
unforced, dissipative equations, as the anisotropy is made more severe?

This is exploratory adversarial IC construction, not a literal embedding of
OpenAI's analytic profile (which requires solving their Appendix A-C
ODE/analytic-continuation system — out of scope here).

Claim: claim_l3_aniso_collapse
  The CZ/delta-gain mechanism (Route 2) survives an anisotropic
  swirl-plus-outflow adversarial IC at increasing anisotropy severity
  (aspect_ratio -> 0), under the exact unforced equations (f = 0).

PASS: delta_max > 0 and theta_negative == False for every aspect_ratio
      in the sweep, at the tested resolution N.
FAIL: either diagnostic fails at some aspect_ratio before the resolution
      limit is reached — a genuine new adversarial threshold, not a
      resolution artifact.

Resolution caveat (same epistemic status as M9 blowup_search_3D.py):
finite N cannot resolve arbitrarily thin cores. A PASS at N grid points
down to a given aspect_ratio does not rule out failure at ratios thinner
than the grid can resolve (ell_z below ~2 grid cells).
"""

from __future__ import annotations

import os
import sys
import time
from datetime import datetime

import numpy as np

# ── path setup ─────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from layer3.route2_3D_jax import (
    make_solver_jax,
    set_ic_jax,
    run_jax,
    kinetic_energy_jax,
    _log_result,
    _JAX_BACKEND,
    _NU_DEFAULT,
)

# Reference energy: match the Wang et al. adversarial_min IC's kinetic energy
# at the same N, so this family is not a weaker test by construction.
from layer2.self_similar_IC import generate_adversarial_min


# ── aspect-ratio sweep (PRD-style adversarial sweep, analogous to LAMBDA_SWEEP_VALUES) ──

ANISO_SWEEP_VALUES = [1.0, 0.5, 0.2, 0.1]          # N=64 production sweep
ANISO_SWEEP_VALUES_HIRES = [0.05]                   # requires N=128 to resolve


def _reference_energy(N: int, seed: int = 42) -> float:
    """Kinetic energy of the existing adversarial_min Wang profile at this N."""
    ic_2d = generate_adversarial_min(lambda_val=0.05, N=N, seed=seed)
    u2, v2 = ic_2d["u"], ic_2d["v"]
    return float(np.mean(u2**2 + v2**2) / 2.0)


def aniso_collapse_ic_3D(
    N: int,
    aspect_ratio: float,
    seed: int = 42,
    target_E0: float | None = None,
) -> dict:
    """Axisymmetric swirl-plus-axial-outflow vortex column, no forcing.

    Geometry (qualitative surrogate for the OpenAI collapse geometry):
      u_theta(r,z) = Gamma * (r/ellr) * exp(-(r/ellr)^2 - (z/ellz)^2)
      u_z(r,z)     = Gamma * bias_slope * (z/ellz + z0_bias)
                            * exp(-(r/ellr)^2 - (z/ellz)^2)
      u_r          recovered by enforcing div(u) = 0 (Leray projection,
                   applied again explicitly by set_ic_jax downstream)

    ellz = aspect_ratio * ellr. aspect_ratio = 1.0 is an isotropic Gaussian
    vortex ring (baseline); aspect_ratio -> 0 is the severe anisotropic
    "slender collapsing column" geometry that OpenAI's forced construction
    needed. Gamma is calibrated so the total kinetic energy matches the
    existing adversarial_min Wang profile at the same N (target_E0), unless
    overridden.

    Returns dict with keys u, v, w (each shape (N,N,N)), plus metadata,
    compatible with route2_3D_jax.set_ic_jax() — same schema as wang_ic_3D().
    """
    if aspect_ratio <= 0.0:
        raise ValueError(f"aspect_ratio must be > 0, got {aspect_ratio}")

    domain = 2.0 * np.pi
    x1d = np.linspace(0.0, domain, N, endpoint=False)
    xx, yy, zz = np.meshgrid(x1d, x1d, x1d, indexing="ij")

    x0 = y0 = z0 = domain / 2.0
    r = np.sqrt((xx - x0) ** 2 + (yy - y0) ** 2)
    z = zz - z0
    theta_ang = np.arctan2(yy - y0, xx - x0)

    ellr = domain / 8.0
    ellz = aspect_ratio * ellr
    ellz = max(ellz, domain / N)  # never thinner than one grid cell

    envelope = np.exp(-((r / ellr) ** 2) - ((z / ellz) ** 2))

    u_theta = (r / ellr) * envelope
    # small asymmetric axial bias, mirroring the paper's deliberate breaking
    # of exact z -> -z reflection symmetry (Sec 2.1: "small upward bias and
    # nonzero velocity at z = 0")
    bias = 0.15
    u_z = ((z / ellz) + bias) * envelope

    # Cylindrical -> Cartesian for the swirl component
    u_swirl = -np.sin(theta_ang) * u_theta
    v_swirl = np.cos(theta_ang) * u_theta
    w_raw = u_z

    Gamma_raw = 1.0
    u = Gamma_raw * u_swirl
    v = Gamma_raw * v_swirl
    w = Gamma_raw * w_raw

    E_raw = float(np.mean(u**2 + v**2 + w**2) / 2.0)
    if target_E0 is None:
        target_E0 = _reference_energy(N, seed=seed)
    Gamma = np.sqrt(target_E0 / max(E_raw, 1e-30))

    u *= Gamma
    v *= Gamma
    w *= Gamma

    E0 = float(np.mean(u**2 + v**2 + w**2) / 2.0)

    return {
        "u": u,
        "v": v,
        "w": w,
        "E0_analytic": E0,
        "name": "aniso_collapse",
        "aspect_ratio": aspect_ratio,
        "ellr": ellr,
        "ellz": ellz,
        "Gamma": float(Gamma),
        "N": N,
    }


def _exp_suffix(aspect_ratio: float) -> str:
    return f"ANISO{aspect_ratio:.3f}".replace(".", "")


def run_aniso_route2_jax(
    aspect_ratio: float,
    N: int = 64,
    t_end: float = 1.0,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    max_steps: int = 3000,
    seed: int = 42,
    verbose: bool = False,
    db_path: str | None = None,
) -> dict:
    """Run Route 2 JAX (exact unforced Prize equations, f=0) from an
    anisotropic collapsing-vortex 3D IC.

    Pass criterion:
      delta_max > 0 and theta_negative == False  (Lemma 2.5 holds)

    Returns experiment result dict (same schema as run_wang_route2_jax).
    """
    exp_id = f"EXP-L3-R2-{_exp_suffix(aspect_ratio)}-JAX-{N:03d}"

    ic = aniso_collapse_ic_3D(N=N, aspect_ratio=aspect_ratio, seed=seed)
    solver = make_solver_jax(N=N, nu=nu, eps_param=eps_param)
    solver = set_ic_jax(solver, ic)

    E0_computed = float(kinetic_energy_jax(solver["u"], solver["v"], solver["w"]))

    if verbose:
        print(f"\n  Route 2 JAX aniso-collapse  aspect_ratio={aspect_ratio:.4f}")
        print(f"  backend={_JAX_BACKEND}  N={N}^3  nu={nu}  eps={eps_param}")
        print(f"  E0={E0_computed:.6f}  ellr={ic['ellr']:.4f}  ellz={ic['ellz']:.4f}")

    t_wall0 = time.perf_counter()
    result = run_jax(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)
    wall = time.perf_counter() - t_wall0

    passed = not result["theta_negative"] and result["delta_max"] > 0.0

    exp = {
        "exp_id": exp_id,
        "claim_id": "claim_l3_aniso_collapse",
        "timestamp": datetime.utcnow().isoformat(),
        "backend": _JAX_BACKEND,
        "ic": f"aniso_collapse_ar{aspect_ratio:.3f}",
        "N": N,
        "nu": nu,
        "eps_param": eps_param,
        "t_end": t_end,
        "n_steps": result["n_steps"],
        "t_final": result["t_final"],
        "E0_computed": E0_computed,
        "E_final": result["E_final"],
        "energy_decayed": result["energy_decayed"],
        "theta_negative": result["theta_negative"],
        "delta_max": result["delta_max"],
        "eps_lps_min": result["eps_lps_min"],
        "verdict": "PASS" if passed else "FAIL",
        "wall_time_s": wall,
        "compile_time_s": result.get("compile_time_s", 0.0),
        "step_ms_mean": result["step_ms_mean"],
        "key_metric": (
            f"aspect_ratio={aspect_ratio:.3f}  "
            f"delta_max={result['delta_max']:.2f}  "
            f"eps_LPS_min={result['eps_lps_min']:.4f}  "
            f"step={result['step_ms_mean']:.1f}ms"
        ),
    }

    if verbose:
        print(f"\n  Results: {exp_id}")
        print(f"    aspect_ratio: {aspect_ratio:.3f}")
        print(f"    theta >= 0:   {'OK' if not result['theta_negative'] else 'FAIL'}")
        print(f"    delta_max:    {result['delta_max']:.2f}  "
              f"{'OK (>0)' if result['delta_max'] > 0 else 'FAIL (=0)'}")
        print(f"    eps_LPS_min:  {result['eps_lps_min']:.4f}")
        print(f"    Wall time:    {wall:.1f} s")
        print(f"    Verdict:      {exp['verdict']}")

    if db_path is not None:
        _log_result(db_path, exp)

    return exp


def run_aniso_sweep_route2_jax(
    aspect_ratios: list[float] | None = None,
    N: int = 64,
    t_end: float = 1.0,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    max_steps: int = 3000,
    seed: int = 42,
    verbose: bool = True,
    db_path: str | None = None,
) -> dict[float, dict]:
    """Run the full anisotropy sweep at fixed N. Returns dict keyed by aspect_ratio."""
    if aspect_ratios is None:
        aspect_ratios = ANISO_SWEEP_VALUES

    results: dict[float, dict] = {}
    for ar in aspect_ratios:
        r = run_aniso_route2_jax(
            ar, N=N, t_end=t_end, nu=nu, eps_param=eps_param,
            max_steps=max_steps, seed=seed, verbose=verbose, db_path=db_path,
        )
        results[ar] = r

    if verbose:
        print("\n" + "=" * 60)
        print(f"Anisotropic Collapse Sweep Summary (N={N})")
        print("=" * 60)
        print(f"{'aspect_ratio':>12} {'delta_max':>10} {'Verdict'}")
        print("-" * 40)
        for ar, r in results.items():
            print(f"{ar:>12.3f} {r['delta_max']:>10.2f}  {r['verdict']}")
        all_pass = all(r["verdict"] == "PASS" for r in results.values())
        print("-" * 40)
        print(f"Overall: {'ALL PASS' if all_pass else 'SOME FAIL'}")

    return results


# ── diagnostics ────────────────────────────────────────────────────────────────

def divergence_rms_3D(u: np.ndarray, v: np.ndarray, w: np.ndarray, N: int) -> float:
    """Compute RMS divergence of (u,v,w) using spectral differentiation."""
    kx1d = np.fft.fftfreq(N, d=1.0 / N)
    kx, ky, kz = np.meshgrid(kx1d, kx1d, kx1d, indexing="ij")

    uh = np.fft.fftn(u)
    vh = np.fft.fftn(v)
    wh = np.fft.fftn(w)

    div_hat = 1j * kx * uh + 1j * ky * vh + 1j * kz * wh
    div_field = np.real(np.fft.ifftn(div_hat))
    return float(np.sqrt(np.mean(div_field**2)))


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=64)
    p.add_argument("--t-end", type=float, default=1.0)
    p.add_argument("--nu", type=float, default=_NU_DEFAULT)
    p.add_argument("--eps-param", type=float, default=0.1)
    p.add_argument("--max-steps", type=int, default=3000)
    p.add_argument("--db", default="results/results.db")
    args = p.parse_args()

    run_aniso_sweep_route2_jax(
        N=args.N, t_end=args.t_end, nu=args.nu, eps_param=args.eps_param,
        max_steps=args.max_steps, verbose=True, db_path=args.db,
    )
