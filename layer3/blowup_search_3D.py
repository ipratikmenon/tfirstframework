"""
blowup_search_3D.py — M9 Adversarial Blowup Search (Route 1, 3D)
=================================================================
Strongest single test of the T-First programme: anti-parallel Gaussian
vortex tubes (Kerr/Hou-Li type) + T₀=T_min to maximise vortex stretching
and minimise thermal resistance.

Physics (Route 1):
  ∂_t u + (u·∇)u = −∇p + ∇·(μ(T)∇u)    div u = 0
  ∂_t T + u·∇T   = A(T)·ΔT + μ(T)|∇u|²/ρcv

Adversarial IC:
  Two anti-parallel Gaussian vortex tubes displaced in y.
  Velocity recovered from vorticity via Biot-Savart: û = i·k×ω̂/|k|².
  T₀ = T_base (set low, e.g. 150 K, to minimise A(T) = k/(ρcv) ∝ T^1.82).

Monitor:
  A_min(t) — must stay > 0 (second law; central claim)
  Z(t) = ½‖ω‖²  — enstrophy
  P(t) = ½‖∇ω‖² — palinstrophy
  ‖u‖_{L∞}(t)    — max velocity

PASS criterion:
  A_min_global > 0  AND  blowup_detected = False

Experiments:
  EXP-L3-R1-ADV-064  — N=64³   smoke test  (~2 s wall)
  EXP-L3-R1-ADV-128  — N=128³  production  (~30 s wall)
  EXP-L3-R1-ADV-256  — N=256³  high-res    (~10 min wall)
"""

from __future__ import annotations

import os
import sys
import time
import json
import sqlite3
import numpy as np
from datetime import datetime

# ── JAX import (identical pattern to route1_3D_jax.py) ───────────────────────
try:
    import jax
    import jax.numpy as jnp
    from jax import jit

    _JAX_AVAILABLE = True
    _JAX_BACKEND_RAW = jax.default_backend()

    if _JAX_BACKEND_RAW == "METAL":
        _COMPUTE_DEVICE = jax.devices("cpu")[0]
        _JAX_BACKEND = "CPU-JIT"
    else:
        _COMPUTE_DEVICE = jax.devices()[0]
        _JAX_BACKEND = _JAX_BACKEND_RAW

except ImportError:
    _JAX_AVAILABLE = False
    _JAX_BACKEND_RAW = "unavailable"
    _JAX_BACKEND = "unavailable"
    _COMPUTE_DEVICE = None


def _d(arr):
    """Move a JAX array to the compute device."""
    if _COMPUTE_DEVICE is None:
        return arr
    return jax.device_put(arr, _COMPUTE_DEVICE)


# ── path setup ────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, ".."))

from layer3.route1_3D_jax import (
    ideal_props_jax,
    make_grid_jax,
    make_dealias_mask_jax,
    make_nyquist_mask_jax,
    make_jax_solver,
    set_ic_jax,
    run_jax,
    kinetic_energy_jax,
    divergence_rms_jax,
    _project_spectral_jax,
)

# ── Default adversarial parameters ───────────────────────────────────────────
_T_ADV_DEFAULT = 150.0   # K — low T → low A(T) → minimum thermal resistance
_AMP_DEFAULT   = 2.0     # vortex tube peak vorticity
_TUBE_SIGMA    = 0.6     # Gaussian half-width (radians)
_TUBE_SEP      = 1.5     # separation between tube centres (radians)
_BLOWUP_FACTOR = 20.0    # Z(t)/Z(0) threshold for blowup alarm

EXPERIMENT_IDS = {
    64:  "EXP-L3-R1-ADV-064",
    128: "EXP-L3-R1-ADV-128",
    256: "EXP-L3-R1-ADV-256",
}


# ─────────────────────────────────────────────────────────────────────────────
# Vorticity field: anti-parallel Gaussian tubes
# ─────────────────────────────────────────────────────────────────────────────

def make_vortex_omega(
    N: int,
    amp: float = _AMP_DEFAULT,
    sigma: float = _TUBE_SIGMA,
    separation: float = _TUBE_SEP,
    seed: int = 42,
) -> tuple:
    """Anti-parallel Gaussian vortex tubes.

    Two counter-rotating tubes aligned with the z-axis, displaced ±separation/2
    in the y-direction.  A small random perturbation in x breaks symmetry to
    allow the tubes to interact and stretch.

    Returns: (omega_x, omega_y, omega_z) each shape (N, N, N), NumPy float32.
    """
    rng = np.random.default_rng(seed)
    x = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False).astype(np.float32)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    cy = np.pi           # y-centre of domain
    d  = separation / 2  # half-separation

    # Tube 1: centre at (π, cy − d, π), ω in +z direction
    r1_sq = (X - np.pi) ** 2 + (Y - (cy - d)) ** 2
    profile1 = amp * np.exp(-r1_sq / (2.0 * sigma ** 2)).astype(np.float32)

    # Tube 2: centre at (π, cy + d, π), ω in −z direction (anti-parallel)
    r2_sq = (X - np.pi) ** 2 + (Y - (cy + d)) ** 2
    profile2 = amp * np.exp(-r2_sq / (2.0 * sigma ** 2)).astype(np.float32)

    # Small x-axis perturbation breaks the planar symmetry so tubes interact
    noise_amp = 0.05 * amp
    noise = noise_amp * rng.standard_normal((N, N, N)).astype(np.float32)

    omega_x = noise
    omega_y = np.zeros((N, N, N), dtype=np.float32)
    omega_z = profile1 - profile2     # anti-parallel in z

    return omega_x, omega_y, omega_z


# ─────────────────────────────────────────────────────────────────────────────
# Biot-Savart: recover velocity from vorticity
# ─────────────────────────────────────────────────────────────────────────────

def biot_savart_spectral(
    omega_x: np.ndarray,
    omega_y: np.ndarray,
    omega_z: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    kz: np.ndarray,
    k2: np.ndarray,
) -> tuple:
    """Recover velocity from vorticity via Biot-Savart in spectral space.

    û = i·k × ω̂ / |k|²

    All inputs may be NumPy arrays (called before JAX JIT context).
    Returns (u, v, w) as float32 NumPy arrays.
    """
    ox_hat = np.fft.fftn(omega_x)
    oy_hat = np.fft.fftn(omega_y)
    oz_hat = np.fft.fftn(omega_z)

    k2_safe = np.where(k2 > 0.0, k2, 1.0)

    # k × ω̂ = (ky*oz - kz*oy, kz*ox - kx*oz, kx*oy - ky*ox)
    cross_x = ky * oz_hat - kz * oy_hat
    cross_y = kz * ox_hat - kx * oz_hat
    cross_z = kx * oy_hat - ky * ox_hat

    valid = (k2 > 0.0).astype(np.float32)

    u_hat = 1j * cross_x / k2_safe * valid
    v_hat = 1j * cross_y / k2_safe * valid
    w_hat = 1j * cross_z / k2_safe * valid

    u = np.real(np.fft.ifftn(u_hat)).astype(np.float32)
    v = np.real(np.fft.ifftn(v_hat)).astype(np.float32)
    w = np.real(np.fft.ifftn(w_hat)).astype(np.float32)

    return u, v, w


# ─────────────────────────────────────────────────────────────────────────────
# Full adversarial IC
# ─────────────────────────────────────────────────────────────────────────────

def make_adversarial_ic_jax(
    N: int,
    T_base: float = _T_ADV_DEFAULT,
    amp: float = _AMP_DEFAULT,
    sigma: float = _TUBE_SIGMA,
    separation: float = _TUBE_SEP,
    seed: int = 42,
    dtype=None,
) -> dict:
    """Anti-parallel vortex tubes IC as JAX arrays on the compute device.

    Returns IC dict compatible with set_ic_jax:
      u, v, w, T — JAX arrays shape (N, N, N)
      name, N, T_base, amp, sigma, separation, Z0
    """
    if dtype is None:
        dtype = jnp.float32

    # Build vorticity in NumPy
    kx_np, ky_np, kz_np = np.meshgrid(
        np.fft.fftfreq(N, d=1.0 / N).astype(np.float32),
        np.fft.fftfreq(N, d=1.0 / N).astype(np.float32),
        np.fft.fftfreq(N, d=1.0 / N).astype(np.float32),
        indexing="ij",
    )
    k2_np = kx_np ** 2 + ky_np ** 2 + kz_np ** 2

    omega_x, omega_y, omega_z = make_vortex_omega(
        N, amp=amp, sigma=sigma, separation=separation, seed=seed
    )

    # Biot-Savart: ω → u
    u_np, v_np, w_np = biot_savart_spectral(
        omega_x, omega_y, omega_z, kx_np, ky_np, kz_np, k2_np
    )

    # Temperature field: uniform T_base (adversarial: low T → low A(T))
    T_np = np.full((N, N, N), T_base, dtype=np.float32)

    # Initial enstrophy Z₀ = ½ ‖ω‖²_{L²} (used for blowup alarm threshold)
    Z0 = 0.5 * float(np.mean(omega_x ** 2 + omega_y ** 2 + omega_z ** 2))

    return {
        "u": _d(jnp.array(u_np, dtype=dtype)),
        "v": _d(jnp.array(v_np, dtype=dtype)),
        "w": _d(jnp.array(w_np, dtype=dtype)),
        "T": _d(jnp.array(T_np, dtype=dtype)),
        "name": "adversarial_antitube",
        "N": N,
        "T_base": T_base,
        "amp": amp,
        "sigma": sigma,
        "separation": separation,
        "Z0": Z0,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Diagnostic quantities
# ─────────────────────────────────────────────────────────────────────────────

def compute_enstrophy_jax(u, v, w, kx, ky, kz) -> float:
    """Z = ½ ‖ω‖²_{L²} where ω = ∇ × u."""
    u_hat = jnp.fft.fftn(u)
    v_hat = jnp.fft.fftn(v)
    w_hat = jnp.fft.fftn(w)

    def deriv_hat(f_hat, k): return 1j * k * f_hat

    ox = jnp.real(jnp.fft.ifftn(deriv_hat(w_hat, ky) - deriv_hat(v_hat, kz)))
    oy = jnp.real(jnp.fft.ifftn(deriv_hat(u_hat, kz) - deriv_hat(w_hat, kx)))
    oz = jnp.real(jnp.fft.ifftn(deriv_hat(v_hat, kx) - deriv_hat(u_hat, ky)))

    return float(0.5 * jnp.mean(ox ** 2 + oy ** 2 + oz ** 2))


def compute_palinstrophy_jax(u, v, w, kx, ky, kz, k2) -> float:
    """P = ½ ‖∇ω‖²_{L²}."""
    u_hat = jnp.fft.fftn(u)
    v_hat = jnp.fft.fftn(v)
    w_hat = jnp.fft.fftn(w)

    def deriv_hat(f_hat, k): return 1j * k * f_hat

    ox_hat = deriv_hat(w_hat, ky) - deriv_hat(v_hat, kz)
    oy_hat = deriv_hat(u_hat, kz) - deriv_hat(w_hat, kx)
    oz_hat = deriv_hat(v_hat, kx) - deriv_hat(u_hat, ky)

    def grad_sq_hat(f_hat):
        return jnp.real(jnp.fft.ifftn(-k2 * f_hat)) ** 2

    return float(0.5 * jnp.mean(
        grad_sq_hat(ox_hat) + grad_sq_hat(oy_hat) + grad_sq_hat(oz_hat)
    ))


def compute_u_linf_jax(u, v, w) -> float:
    """‖u‖_{L∞} = max |u|."""
    return float(jnp.max(jnp.sqrt(u ** 2 + v ** 2 + w ** 2)))


def check_blowup(Z: float, Z0: float, factor: float = _BLOWUP_FACTOR) -> bool:
    """Return True if enstrophy has grown by more than `factor` from Z0."""
    if Z0 <= 0.0:
        return False
    return Z > factor * Z0


# ─────────────────────────────────────────────────────────────────────────────
# Main blowup search runner
# ─────────────────────────────────────────────────────────────────────────────

def run_blowup_search(
    N: int = 128,
    T_base: float = _T_ADV_DEFAULT,
    amp: float = _AMP_DEFAULT,
    sigma: float = _TUBE_SIGMA,
    separation: float = _TUBE_SEP,
    t_end: float = 1.0,
    max_steps: int = 5000,
    verbose: bool = False,
    db_path: str = None,
    seed: int = 42,
    blowup_factor: float = _BLOWUP_FACTOR,
) -> dict:
    """Run adversarial blowup search on N³ grid.

    Constructs anti-parallel vortex tube IC + low-T thermal resistance.
    Runs Route 1 JAX solver, monitoring A_min, Z, P, ‖u‖_∞ at each step.

    PASS criterion: A_min_global > 0 AND blowup_detected = False.

    Returns result dict with keys:
      exp_id, verdict, A_min_global, Z0, Z_max, P_max, u_linf_max,
      blowup_detected, n_steps, t_final, wall_time_s, key_metric, history
    """
    exp_id = EXPERIMENT_IDS.get(N, f"EXP-L3-R1-ADV-{N:03d}")

    if verbose:
        print(f"\n{'='*60}")
        print(f"M9 Blowup Search — {exp_id}")
        print(f"N={N}³  T_base={T_base}K  amp={amp}  t_end={t_end}")
        print(f"{'='*60}")

    # ── Build IC ──────────────────────────────────────────────────────────────
    solver = make_jax_solver(N=N, T_base=T_base, dt_max=0.05)
    ic = make_adversarial_ic_jax(
        N, T_base=T_base, amp=amp, sigma=sigma,
        separation=separation, seed=seed, dtype=solver["dtype"]
    )
    solver = set_ic_jax(solver, ic)
    Z0 = ic["Z0"]

    E_initial = kinetic_energy_jax(solver["u"], solver["v"], solver["w"])
    Z_initial = compute_enstrophy_jax(
        solver["u"], solver["v"], solver["w"],
        solver["kx"], solver["ky"], solver["kz"]
    )

    if verbose:
        print(f"  IC: E₀={E_initial:.4e}  Z₀={Z_initial:.4e}  T_base={T_base}K")

    # ── Run with extended diagnostics ─────────────────────────────────────────
    from layer3.route1_3D_jax import (
        compute_cfl_dt_jax, solver_step_jax, ideal_props_jax,
    )

    kx, ky, kz, k2 = solver["kx"], solver["ky"], solver["kz"], solver["k2"]
    dealias_mask = solver["dealias_mask"]
    nyquist_mask = solver["nyquist_mask"]
    u, v, w, T = solver["u"], solver["v"], solver["w"], solver["T"]

    A_min_global = float("inf")
    Z_max        = Z_initial
    P_max        = 0.0
    u_linf_max   = 0.0
    blowup_detected = False
    history = []
    verdict = "PASS"

    t_wall_start = time.perf_counter()
    t = 0.0

    for step_idx in range(max_steps):
        if t >= t_end:
            break

        dt = compute_cfl_dt_jax(u, v, w, T, N, dt_max=solver["dt_max"])
        dt = min(dt, t_end - t)
        dt_jnp = jnp.array(dt, dtype=solver["dtype"])

        new_u, new_v, new_w, new_T = solver_step_jax(
            u, v, w, T, kx, ky, kz, k2, dealias_mask, nyquist_mask, dt_jnp
        )
        jax.block_until_ready((new_u, new_v, new_w, new_T))
        u, v, w, T = new_u, new_v, new_w, new_T
        t += dt

        # A_min diagnostic
        props = ideal_props_jax(T)
        A_arr = props["k"] / (props["rho"] * props["cv"])
        A_min = float(jnp.min(A_arr))
        A_min_global = min(A_min_global, A_min)

        if A_min <= 0.0:
            verdict = "FAIL"
            if verbose:
                print(f"  [FAIL] Step {step_idx}: A_min={A_min:.3e} ≤ 0")
            break

        # Enstrophy, palinstrophy, L∞
        Z = compute_enstrophy_jax(u, v, w, kx, ky, kz)
        P = compute_palinstrophy_jax(u, v, w, kx, ky, kz, k2)
        u_linf = compute_u_linf_jax(u, v, w)

        Z_max = max(Z_max, Z)
        P_max = max(P_max, P)
        u_linf_max = max(u_linf_max, u_linf)

        # Blowup alarm
        if check_blowup(Z, Z0, factor=blowup_factor):
            blowup_detected = True
            verdict = "FAIL"
            if verbose:
                print(f"  [BLOWUP] Step {step_idx}: Z={Z:.3e} > {blowup_factor}×Z₀={blowup_factor*Z0:.3e}")
            break

        diag = {
            "t": t, "dt": dt,
            "A_min": A_min, "Z": Z, "P": P, "u_linf": u_linf,
            "E": float(0.5 * jnp.mean(u**2 + v**2 + w**2)),
        }
        history.append(diag)

        if verbose and step_idx % max(1, max_steps // 10) == 0:
            print(
                f"  step {step_idx:4d}  t={t:.4f}  A_min={A_min:.3e}  "
                f"Z={Z:.3e}  Z/Z₀={Z/max(Z0,1e-20):.2f}  u∞={u_linf:.3e}"
            )

    t_wall = time.perf_counter() - t_wall_start
    n_steps = len(history)

    if verbose:
        print(f"\n  Results:")
        print(f"    A_min_global = {A_min_global:.4e}")
        print(f"    Z_max / Z0   = {Z_max / max(Z0, 1e-20):.3f}")
        print(f"    P_max        = {P_max:.4e}")
        print(f"    ‖u‖_∞ max   = {u_linf_max:.4e}")
        print(f"    blowup_detected = {blowup_detected}")
        print(f"    verdict      = {verdict}")
        print(f"    wall time    = {t_wall:.1f} s  ({n_steps} steps)")

    key_metric = (
        f"A_min={A_min_global:.4e}  Z_max/Z0={Z_max/max(Z0,1e-20):.3f}  "
        f"P_max={P_max:.4e}  u∞={u_linf_max:.4e}  blowup={blowup_detected}"
    )

    result = {
        "exp_id": exp_id,
        "claim_id": "claim_m9_blowup_search_3D",
        "timestamp": datetime.utcnow().isoformat(),
        "backend": _JAX_BACKEND,
        "N": N,
        "T_base": T_base,
        "amp": amp,
        "sigma": sigma,
        "separation": separation,
        "t_end": t_end,
        "n_steps": n_steps,
        "t_final": t,
        "Z0": Z0,
        "E_initial": E_initial,
        "A_min_global": A_min_global,
        "Z_max": Z_max,
        "P_max": P_max,
        "u_linf_max": u_linf_max,
        "blowup_detected": blowup_detected,
        "verdict": verdict,
        "wall_time_s": t_wall,
        "key_metric": key_metric,
        "history": history,
        "params_json": json.dumps({
            "N": N, "T_base": T_base, "amp": amp, "sigma": sigma,
            "separation": separation, "t_end": t_end, "blowup_factor": blowup_factor,
        }),
    }

    if db_path is not None:
        _log_result_m9(db_path, result)

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Multi-resolution sweep
# ─────────────────────────────────────────────────────────────────────────────

def run_m9_sweep(
    resolutions: list = None,
    t_end: float = 1.0,
    verbose: bool = True,
    db_path: str = None,
    **kwargs,
) -> dict:
    """Run blowup search at multiple N values.

    Default: [64, 128] (256 requires explicit request due to ~10 min wall time).
    Returns dict mapping N → result.
    """
    if resolutions is None:
        resolutions = [64, 128]

    results = {}
    for N in resolutions:
        if verbose:
            print(f"\n{'─'*50}")
            print(f"Running N={N}³ blowup search...")
        results[N] = run_blowup_search(
            N=N, t_end=t_end, verbose=verbose, db_path=db_path, **kwargs
        )

    if verbose:
        print(f"\n{'='*60}")
        print("M9 SWEEP SUMMARY")
        print(f"{'='*60}")
        for N, r in results.items():
            print(
                f"  N={N:3d}³  {r['verdict']:4s}  "
                f"A_min={r['A_min_global']:.4e}  "
                f"Z_max/Z0={r['Z_max']/max(r['Z0'],1e-20):.3f}  "
                f"blowup={r['blowup_detected']}"
            )

    return results


# ─────────────────────────────────────────────────────────────────────────────
# SQLite results logging
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_m9_db(db_path: str):
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS m9_blowup_search (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id           TEXT UNIQUE NOT NULL,
            claim_id         TEXT NOT NULL,
            timestamp        TEXT NOT NULL,
            backend          TEXT,
            N                INTEGER NOT NULL,
            T_base           REAL,
            amp              REAL,
            sigma            REAL,
            separation       REAL,
            t_end            REAL,
            n_steps          INTEGER,
            t_final          REAL,
            Z0               REAL,
            E_initial        REAL,
            A_min_global     REAL,
            Z_max            REAL,
            P_max            REAL,
            u_linf_max       REAL,
            blowup_detected  INTEGER,
            verdict          TEXT NOT NULL,
            wall_time_s      REAL,
            key_metric       TEXT,
            params_json      TEXT
        )
    """)
    conn.commit()
    return conn


def _log_result_m9(db_path: str, result: dict) -> None:
    conn = _ensure_m9_db(db_path)
    cols = [
        "exp_id", "claim_id", "timestamp", "backend", "N",
        "T_base", "amp", "sigma", "separation", "t_end",
        "n_steps", "t_final", "Z0", "E_initial",
        "A_min_global", "Z_max", "P_max", "u_linf_max",
        "blowup_detected", "verdict", "wall_time_s",
        "key_metric", "params_json",
    ]
    values = tuple(
        int(result[c]) if isinstance(result.get(c), bool) else result.get(c)
        for c in cols
    )
    ph = ", ".join("?" for _ in cols)
    conn.execute(
        f"INSERT OR REPLACE INTO m9_blowup_search ({', '.join(cols)}) VALUES ({ph})",
        values,
    )
    conn.commit()
    conn.close()


def query_m9_results(db_path: str) -> list:
    if not os.path.exists(db_path):
        return []
    conn = sqlite3.connect(db_path)
    cur = conn.execute("SELECT * FROM m9_blowup_search ORDER BY timestamp DESC")
    cols = [d[0] for d in cur.description]
    rows = [dict(zip(cols, r)) for r in cur.fetchall()]
    conn.close()
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def _default_db() -> str:
    d = os.path.join(_HERE, "..", "results")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, "results.db")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="M9 Blowup Search 3D")
    parser.add_argument("--N", type=int, default=128)
    parser.add_argument("--t-end", type=float, default=1.0)
    parser.add_argument("--T-base", type=float, default=_T_ADV_DEFAULT)
    parser.add_argument("--amp", type=float, default=_AMP_DEFAULT)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--sweep", action="store_true",
                        help="Run N=64 and N=128 sweep")
    parser.add_argument("--with-256", action="store_true",
                        help="Include N=256 in sweep (~10 min)")
    parser.add_argument("--query", action="store_true")
    parser.add_argument("--db", type=str, default=None)
    args = parser.parse_args()

    db_path = args.db or _default_db()

    if args.query:
        rows = query_m9_results(db_path)
        for r in rows:
            print(
                f"{r['exp_id']}  N={r['N']}  {r['verdict']}  "
                f"A_min={r['A_min_global']:.4e}  blowup={r['blowup_detected']}"
            )
    elif args.sweep:
        res_list = [64, 128]
        if args.with_256:
            res_list.append(256)
        run_m9_sweep(resolutions=res_list, t_end=args.t_end, verbose=True,
                     T_base=args.T_base, amp=args.amp, db_path=db_path)
    else:
        result = run_blowup_search(
            N=args.N, T_base=args.T_base, amp=args.amp,
            t_end=args.t_end, verbose=args.verbose, db_path=db_path,
        )
        print(f"\nVerdict: {result['verdict']}")
        print(f"  {result['key_metric']}")
