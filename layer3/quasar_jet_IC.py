"""
layer3/quasar_jet_IC.py — Bipolar Swirling Jet IC and Strain Sweep
===================================================================
A bipolar ("quasar-jet") adversarial initial condition: axisymmetric swirl
plus opposing axial outflow either side of a mid-plane, with radial inflow
supplied by incompressibility. NO forcing term is introduced anywhere -- the
IC is fed to the existing, unmodified Route 2 solver (exact unforced Prize
equations, f = 0), exactly as `layer3/aniso_collapse_IC.py` does.

Motivation
----------
This geometry is not arbitrary. In its steady, unbounded form it is the
classical **Burgers(-Rott) vortex**: swirl + radial inflow + axial outflow
held by an axial strain gamma, an EXACT steady solution of incompressible
Navier-Stokes which is provably regular, with core radius delta = sqrt(4 nu/gamma)
-- viscosity balances stretching exactly. It is also, structurally, the inner
core of OpenAI's forced blowup construction (their Sec. 2.1: fluid spirals
inward and "flows upward and downward on opposite sides of a dividing layer
close to z = 0"), which required forcing to sustain.

The interesting regime for THIS program is that a collimated jet is a
maximally *aligned* configuration. In the perfectly aligned limit ehat is
constant, so w = |grad ehat|^2 = 0 and hence sigma* = 0 -- precisely the
corner the program's own results say a singularity must approach (S30
Thm E; S38 proved sigma*-decay NECESSARY). Yet the Burgers vortex sits at
sigma* = 0 and is perfectly regular, because M stays bounded. So this family
probes the aligned corner directly, and the sweep asks whether driving the
strain balance away from the Burgers value degrades the diagnostics there.

Design decisions worth stating
------------------------------
1. **A tilt perturbation is mandatory, not cosmetic.** A perfectly
   axisymmetric aligned jet has ehat constant, so w == 0 identically and
   every correlation diagnostic in this program degenerates to its
   zero-by-convention branch. `tilt` breaks the alignment so the diagnostics
   are non-degenerate. Setting tilt=0 gives the degenerate reference case
   and is retained for exactly that purpose.
2. **Energy is held FIXED across the sweep; the swept parameter is the
   strain-to-swirl RATIO.** Sweeping the raw strain amplitude would confound
   "more strain" with "more energy" and make every downstream diagnostic
   move for the trivial reason. Fixing E0 isolates the geometry.
3. Radial inflow is not imposed by hand: axial outflow away from the
   mid-plane forces it through incompressibility, and the Leray projection
   inside `set_ic` supplies it exactly.

SCOPE CAVEAT: per S37's standing caveat, numerics on a globally-regular
decaying flow cannot confirm or refute the analytical conjectures; this can
only surface adversarial structure. And the real physical collimator of a
quasar jet is magnetic, which incompressible NS does not have -- the
analogy licenses the GEOMETRY, not a sustaining mechanism.
"""

from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

# Strain-to-swirl ratios. 0.0 = pure swirl (no jet); large = strain-dominated.
JET_SWEEP_RATIOS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]

_DEFAULT_TILT = 0.05


def quasar_jet_ic(
    N: int,
    strain_ratio: float,
    tilt: float = _DEFAULT_TILT,
    core_frac: float = 8.0,
    target_E0: float | None = None,
    seed: int = 42,
) -> dict:
    """Bipolar swirling-jet IC on [0, 2pi)^3, energy-normalised.

    Structure (cylindrical about the box centre, mid-plane at z = pi):
      u_theta = (r/delta) exp(-(r/delta)^2)              swirl
      u_z     = strain_ratio * sin(z - pi) * exp(-(r/delta)^2)
                                                          bipolar outflow
      u_r     supplied by the Leray projection in set_ic (incompressibility)
      plus a small transverse tilt breaking exact alignment.

    `strain_ratio` is the swept shape parameter: the amplitude of the axial
    jet relative to the swirl. Total kinetic energy is normalised to
    `target_E0` (default: match `layer3.aniso_collapse_IC`'s reference so the
    two adversarial families are directly comparable), so the sweep varies
    geometry at fixed energy.

    Returns the same dict schema as `wang_ic_3D` / `aniso_collapse_ic_3D`.
    """
    if strain_ratio < 0.0:
        raise ValueError(f"strain_ratio must be >= 0, got {strain_ratio}")
    if core_frac <= 0.0:
        raise ValueError(f"core_frac must be > 0, got {core_frac}")

    domain = 2.0 * np.pi
    x1d = np.linspace(0.0, domain, N, endpoint=False)
    xx, yy, zz = np.meshgrid(x1d, x1d, x1d, indexing="ij")

    x0 = y0 = z0 = domain / 2.0
    rx, ry = xx - x0, yy - y0
    r = np.sqrt(rx**2 + ry**2)
    zeta = zz - z0                      # mid-plane at zeta = 0
    theta = np.arctan2(ry, rx)

    delta = domain / core_frac
    env = np.exp(-((r / delta) ** 2))

    # Swirl (azimuthal), vanishing on the axis
    u_theta = (r / delta) * env
    u = -np.sin(theta) * u_theta
    v = np.cos(theta) * u_theta

    # Bipolar axial outflow: positive above the mid-plane, negative below
    w = strain_ratio * np.sin(zeta) * env

    # Tilt: breaks exact axisymmetry/alignment so that ehat varies and the
    # correlation diagnostics are non-degenerate (see module docstring).
    if tilt != 0.0:
        u = u + tilt * np.cos(zeta) * env

    if target_E0 is None:
        from layer3.aniso_collapse_IC import _reference_energy
        target_E0 = _reference_energy(N, seed=seed)

    E_raw = float(np.mean(u**2 + v**2 + w**2) / 2.0)
    scale = np.sqrt(target_E0 / max(E_raw, 1e-30))
    u, v, w = u * scale, v * scale, w * scale

    return {
        "u": u, "v": v, "w": w,
        "E0_analytic": float(np.mean(u**2 + v**2 + w**2) / 2.0),
        "name": "quasar_jet",
        "strain_ratio": strain_ratio,
        "tilt": tilt,
        "delta": delta,
        "N": N,
    }


def _tag(strain_ratio: float) -> str:
    return f"JET{strain_ratio:.2f}".replace(".", "")


def run_jet_strain_sweep(
    strain_ratios: list[float] | None = None,
    N: int = 64,
    n_steps: int = 200,
    nu: float | None = None,
    tilt: float = _DEFAULT_TILT,
    record_every: int = 10,
    seed: int = 42,
    verbose: bool = True,
    db_path: str | None = None,
) -> dict[float, dict]:
    """Sweep the strain-to-swirl ratio, tracking scale-r* diagnostics.

    Reuses `layer4.scale_diagnostics.run_scale_experiment` unchanged via its
    `ic_dict_override` hook -- the solver, the diagnostics and the logging
    are all existing, tested code; only the IC is new.
    """
    from layer4.scale_diagnostics import run_scale_experiment, l_dependence_summary

    if strain_ratios is None:
        strain_ratios = JET_SWEEP_RATIOS

    results: dict[float, dict] = {}
    for sr in strain_ratios:
        ic = quasar_jet_ic(N=N, strain_ratio=sr, tilt=tilt, seed=seed)
        if verbose:
            print(f"\n  jet strain_ratio={sr:.2f} (tilt={tilt}, N={N})")
        r = run_scale_experiment(
            ic="jet", N=N, n_steps=n_steps, nu=nu, record_every=record_every,
            seed=seed, verbose=verbose, db_path=db_path,
            exp_id_override=f"EXP-L4-{_tag(sr)}-{N:03d}",
            claim_id_override=f"jet-strain-sweep-{sr:.2f}",
            ic_dict_override=ic,
        )
        r["strain_ratio"] = sr
        results[sr] = r

    if verbose:
        print("\n" + "=" * 78)
        print(f"Jet strain sweep (N={N}, tilt={tilt})")
        print("=" * 78)
        print(f"{'ratio':>7} {'|E|/|B|':>9} {'|Bd|/|B|':>9} {'c_K':>8} "
              f"{'c_K^band':>9} {'sigma*':>11} {'L':>7}  verdict")
        print("-" * 78)
        for sr, r in results.items():
            print(f"{sr:>7.2f} {r['frac_E_rstar_median']:>9.4f} "
                  f"{r['frac_band_rstar_median']:>9.4f} "
                  f"{r['c_K_rstar_median']:>8.3f} {r['c_K_band_median']:>9.3f} "
                  f"{r['sigma_star_median']:>11.3e} {r['L_median']:>7.3f}  "
                  f"{r['verdict']}")
        summ = l_dependence_summary(list(results.values()))
        print("-" * 78)
        print(f"pooled L-regression: c_K slope={summ['c_K_vs_L']['slope']:.3f} "
              f"(r2={summ['c_K_vs_L']['r2']:.3f}, n={summ['c_K_vs_L']['n_used']}), "
              f"c_K^band slope={summ['c_K_band_vs_L']['slope']:.3f} "
              f"(r2={summ['c_K_band_vs_L']['r2']:.3f}, n={summ['c_K_band_vs_L']['n_used']})")
        print("  slope ~ -1 => decaying form c_K <= C/L; slope ~ 0 => boundedness only")

    return results
