"""
wang_profiles_3D.py — Wang et al. Self-Similar Profiles as 3D Initial Conditions
==================================================================================
M6 deliverable (PRD v1.0 FINAL §7.3).

Embeds 2D Wang et al. blow-up profiles (from layer2/self_similar_IC.py) into 3D
as vortex-tube initial conditions for the Route 2 3D JAX solver.

Embedding strategy (vortex-tube + z-perturbation):
  - u_3D(x,y,z) = u_2D(x,y)          (z-independent vortex tube)
  - v_3D(x,y,z) = v_2D(x,y)
  - w_3D(x,y,z) = amp * ũ(x,y) * sin(n * z)   (weak 3D perturbation)
  - Leray projection (inside set_ic_jax) restores exact div u = 0

Physics test (M6 central question):
  Does δ > 0 (Lemma 2.5) survive on Wang et al. adversarial ICs in 3D?
  The 2D profiles are designed to maximise blow-up rate for a given λ.
  If δ_max > 0 here, the CZ mechanism defeats the adversarial structure.

Profiles tested:
  CCF_1st_unstable   λ=0.6057  EXP-L3-R2-CCF1-JAX-064
  CCF_2nd_unstable   λ=0.4703  EXP-L3-R2-CCF2-JAX-064  (critical — PRD §4.3)
  adversarial_min    λ=0.05    EXP-L3-R2-ADV3D-JAX-064  (hardest test)
"""

from __future__ import annotations

import os
import sqlite3
import sys
import time
from datetime import datetime

import numpy as np

# ── path setup ─────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from layer2.self_similar_IC import (
    generate_CCF_profile,
    generate_adversarial_min,
    WANG_PROFILES,
)
from layer3.route2_3D_jax import (
    make_solver_jax,
    set_ic_jax,
    run_jax,
    delta_from_theta_jax,
    kinetic_energy_jax,
    _log_result,
    _JAX_BACKEND,
    _NU_DEFAULT,
)

try:
    import jax.numpy as jnp
    _JAX_AVAILABLE = True
except ImportError:
    _JAX_AVAILABLE = False


# ── profile constants ──────────────────────────────────────────────────────────

# Profiles used in M6 experiments (PRD §7.3)
M6_PROFILES = {
    "CCF_1st_unstable": {"lambda_val": 0.6057, "exp_suffix": "CCF1"},
    "CCF_2nd_unstable": {"lambda_val": 0.4703, "exp_suffix": "CCF2"},
    "adversarial_min":  {"lambda_val": 0.05,   "exp_suffix": "ADV3D"},
}


# ── embedding ──────────────────────────────────────────────────────────────────

def embed_2D_to_3D(
    u_2d: np.ndarray,
    v_2d: np.ndarray,
    N: int,
    perturb_amp: float = 0.02,
    perturb_modes: int = 2,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Embed a 2D velocity field as a 3D vortex tube plus z-perturbation.

    The 2D (u,v) field is broadcast along z (z-independent vortex tube).
    A small w-component proportional to u_2D * sin(m*z) is added to break
    the 2D symmetry while keeping the perturbation physically meaningful.

    Leray projection (applied later by set_ic_jax) restores exact div u = 0.

    Returns: u_3d, v_3d, w_3d  each of shape (N, N, N)
    """
    assert u_2d.shape == (N, N), f"Expected ({N},{N}), got {u_2d.shape}"
    assert v_2d.shape == (N, N), f"Expected ({N},{N}), got {v_2d.shape}"

    z = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)

    # Broadcast 2D → 3D  (shape: N,N,1 tiled to N,N,N)
    u_3d = np.broadcast_to(u_2d[:, :, np.newaxis], (N, N, N)).copy()
    v_3d = np.broadcast_to(v_2d[:, :, np.newaxis], (N, N, N)).copy()

    # z-perturbation in w: amp * |u_2D| * sin(m * z)
    # Scaled to perturb_amp of kinetic energy so it doesn't overwhelm the 2D structure
    rng = np.random.default_rng(seed)
    E_u = float(np.mean(u_2d**2 + v_2d**2))
    scale = perturb_amp * np.sqrt(max(E_u, 1e-10))

    envelope = np.sqrt(u_2d**2 + v_2d**2)   # shape (N, N)
    sign = rng.choice([-1.0, 1.0])

    # Stack z-modes
    w_3d = np.zeros((N, N, N))
    for m in range(1, perturb_modes + 1):
        coeff = scale / perturb_modes
        # shape (N, N, N):  envelope[:,:,None] * sin(m * z[None,None,:])
        w_3d += sign * coeff * envelope[:, :, np.newaxis] * np.sin(m * z[np.newaxis, np.newaxis, :])

    return u_3d, v_3d, w_3d


def wang_ic_3D(
    profile_name: str,
    N: int = 64,
    perturb_amp: float = 0.02,
    seed: int = 42,
) -> dict:
    """Build a 3D IC dict from a named Wang et al. profile.

    Returns dict with keys: u, v, w  (each shape (N,N,N)), plus metadata.
    Compatible with route2_3D_jax.set_ic_jax().
    """
    # Generate 2D profile
    if profile_name == "adversarial_min":
        ic_2d = generate_adversarial_min(lambda_val=0.05, N=N, seed=seed)
        lambda_val = 0.05
        family = "adversarial_min"
    elif profile_name in WANG_PROFILES:
        info = WANG_PROFILES[profile_name]
        lambda_val = info["lambda_val"]
        family = info["family"]
        ic_2d = generate_CCF_profile(lambda_val, N=N, seed=seed)
    else:
        valid = list(WANG_PROFILES.keys()) + ["adversarial_min"]
        raise ValueError(f"Unknown profile '{profile_name}'. Valid: {valid}")

    u_2d = ic_2d["u"]
    v_2d = ic_2d["v"]

    u_3d, v_3d, w_3d = embed_2D_to_3D(
        u_2d, v_2d, N, perturb_amp=perturb_amp, seed=seed
    )

    E0 = float(np.mean(u_3d**2 + v_3d**2 + w_3d**2) / 2.0)

    return {
        "u"            : u_3d,
        "v"            : v_3d,
        "w"            : w_3d,
        "E0_analytic"  : E0,
        "name"         : profile_name,
        "lambda_val"   : lambda_val,
        "family"       : family,
        "omega_max_2d" : float(np.max(np.abs(ic_2d.get("omega", np.zeros((N, N)))))),
        "N"            : N,
    }


# ── experiments ────────────────────────────────────────────────────────────────

def run_wang_route2_jax(
    profile_name: str,
    N: int = 64,
    t_end: float = 1.0,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    max_steps: int = 3000,
    perturb_amp: float = 0.02,
    seed: int = 42,
    verbose: bool = False,
    db_path: str | None = None,
) -> dict:
    """Run Route 2 JAX from a Wang et al. 3D IC.

    Pass criterion:
      δ_max > 0  (Lemma 2.5 holds on adversarial IC)

    Returns experiment result dict (same schema as run_taylor_green_jax).
    """
    info = M6_PROFILES.get(profile_name, {"lambda_val": 0.0, "exp_suffix": profile_name})
    exp_id = f"EXP-L3-R2-{info['exp_suffix']}-JAX-{N:03d}"

    ic = wang_ic_3D(profile_name, N=N, perturb_amp=perturb_amp, seed=seed)
    solver = make_solver_jax(N=N, nu=nu, eps_param=eps_param)
    solver = set_ic_jax(solver, ic)

    E0_computed = float(kinetic_energy_jax(solver["u"], solver["v"], solver["w"]))

    if verbose:
        lv = ic["lambda_val"]
        print(f"\n  Route 2 JAX Wang-3D  profile={profile_name}  λ={lv:.4f}")
        print(f"  backend={_JAX_BACKEND}  N={N}³  ν={nu}  ε={eps_param}")
        print(f"  E₀={E0_computed:.6f}  ω_max_2D={ic['omega_max_2d']:.4f}")

    t_wall0 = time.perf_counter()
    result = run_jax(solver, t_end=t_end, max_steps=max_steps, verbose=verbose)
    wall = time.perf_counter() - t_wall0

    passed = not result["theta_negative"] and result["delta_max"] > 0.0

    exp = {
        "exp_id"         : exp_id,
        "claim_id"       : f"claim_m6_wang_3D_{profile_name}",
        "timestamp"      : datetime.utcnow().isoformat(),
        "backend"        : _JAX_BACKEND,
        "ic"             : profile_name,
        "lambda_val"     : ic["lambda_val"],
        "family"         : ic["family"],
        "N"              : N,
        "nu"             : nu,
        "eps_param"      : eps_param,
        "t_end"          : t_end,
        "n_steps"        : result["n_steps"],
        "t_final"        : result["t_final"],
        "E0_computed"    : E0_computed,
        "E_final"        : result["E_final"],
        "energy_decayed" : result["energy_decayed"],
        "theta_negative" : result["theta_negative"],
        "delta_max"      : result["delta_max"],
        "eps_lps_min"    : result["eps_lps_min"],
        "verdict"        : "PASS" if passed else "FAIL",
        "wall_time_s"    : wall,
        "compile_time_s" : result.get("compile_time_s", 0.0),
        "step_ms_mean"   : result["step_ms_mean"],
        "key_metric"     : (
            f"λ={ic['lambda_val']:.4f}  "
            f"δ_max={result['delta_max']:.2f}  "
            f"ε_LPS_min={result['eps_lps_min']:.4f}  "
            f"step={result['step_ms_mean']:.1f}ms"
        ),
    }

    if verbose:
        print(f"\n  Results: {exp_id}")
        print(f"    λ:           {ic['lambda_val']:.4f}")
        print(f"    θ ≥ 0:       {'✓' if not result['theta_negative'] else '✗'}")
        print(f"    δ_max:       {result['delta_max']:.2f}  {'✓ (>0)' if result['delta_max'] > 0 else '✗ (=0)'}")
        print(f"    ε_LPS_min:   {result['eps_lps_min']:.4f}")
        print(f"    Step time:   {result['step_ms_mean']:.1f} ms")
        print(f"    Wall time:   {wall:.1f} s")
        print(f"    Verdict:     {exp['verdict']}")

    if db_path is not None:
        _log_result(db_path, exp)

    return exp


def run_all_m6_experiments(
    N: int = 64,
    t_end: float = 1.0,
    nu: float = _NU_DEFAULT,
    eps_param: float = 0.1,
    max_steps: int = 3000,
    verbose: bool = True,
    db_path: str | None = None,
) -> dict[str, dict]:
    """Run all M6 Wang et al. 3D experiments.

    Returns dict keyed by profile_name.
    """
    results = {}
    for profile_name in M6_PROFILES:
        r = run_wang_route2_jax(
            profile_name, N=N, t_end=t_end, nu=nu, eps_param=eps_param,
            max_steps=max_steps, verbose=verbose, db_path=db_path,
        )
        results[profile_name] = r

    if verbose:
        print("\n" + "=" * 60)
        print("M6 Wang et al. 3D Summary")
        print("=" * 60)
        print(f"{'Profile':<22} {'λ':>7} {'δ_max':>6} {'Verdict'}")
        print("-" * 50)
        for name, r in results.items():
            print(f"{name:<22} {r['lambda_val']:>7.4f} {r['delta_max']:>6.2f}  {r['verdict']}")
        all_pass = all(r["verdict"] == "PASS" for r in results.values())
        print("-" * 50)
        print(f"M6 overall: {'✓ ALL PASS' if all_pass else '✗ SOME FAIL'}")

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


def profile_energy_3D(u: np.ndarray, v: np.ndarray, w: np.ndarray) -> float:
    return float(np.mean(u**2 + v**2 + w**2) / 2.0)


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--profile", default="CCF_2nd_unstable",
                   choices=list(M6_PROFILES.keys()) + ["all"])
    p.add_argument("--N", type=int, default=64)
    p.add_argument("--t-end", type=float, default=1.0)
    p.add_argument("--nu", type=float, default=_NU_DEFAULT)
    p.add_argument("--eps-param", type=float, default=0.1)
    p.add_argument("--max-steps", type=int, default=3000)
    p.add_argument("--db", default="results/results.db")
    args = p.parse_args()

    if args.profile == "all":
        run_all_m6_experiments(
            N=args.N, t_end=args.t_end, nu=args.nu,
            eps_param=args.eps_param, max_steps=args.max_steps,
            verbose=True, db_path=args.db,
        )
    else:
        run_wang_route2_jax(
            args.profile, N=args.N, t_end=args.t_end, nu=args.nu,
            eps_param=args.eps_param, max_steps=args.max_steps,
            verbose=True, db_path=args.db,
        )
