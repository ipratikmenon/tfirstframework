"""
Adversarial Gamma(r) Search — layer4/gamma_adversarial_search.py

Constructs synthetic velocity fields designed to maximise
    Gamma(r) = integral_0^T (fint_{B_r} |grad u|^2)^{1/2} dt
then checks whether NS evolution (modelled as free viscous decay) suppresses them.

Goal:
  (a) Find a flow where Gamma(r) ~ r^{-beta} for some beta > 0
      (violates Claim A conjecture -> gap may be real), OR
  (b) Confirm empirically that NS dynamics suppress Gamma even for
      adversarially designed ICs.

Experiment ID: EXP-L4-GAM-ADV-001
"""
from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Import shared utilities from mean_morrey_diagnostic
# ---------------------------------------------------------------------------
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from layer4.mean_morrey_diagnostic import ball_mask, compute_spatial_mean, L_BOX

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DB_PATH = Path(__file__).parent.parent / "results" / "results.db"

_CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS gamma_adversarial_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    exp_id          TEXT    NOT NULL,
    timestamp       TEXT,
    grid_N          INTEGER,
    nu              REAL,
    T_final         REAL,
    n_steps         INTEGER,
    seed            INTEGER,
    r_fracs_json    TEXT,
    gamma_initial   REAL,
    gamma_evolved   REAL,
    gamma_exponent  REAL,
    suppressed      INTEGER,
    verdict         TEXT,
    notes           TEXT
)
"""

_INSERT_RESULT = """
INSERT INTO gamma_adversarial_results
    (exp_id, timestamp, grid_N, nu, T_final, n_steps, seed,
     r_fracs_json, gamma_initial, gamma_evolved, gamma_exponent,
     suppressed, verdict, notes)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


# ---------------------------------------------------------------------------
# Helper: FFT-based gradient magnitude squared, averaged over a ball
# ---------------------------------------------------------------------------

def _fft_grad_mag_sq(u: np.ndarray) -> np.ndarray:
    """Compute |grad u|^2 = sum_{i,j} (partial_i u_j)^2 on an N^3 grid.

    Uses spectral differentiation (FFT). Returns shape (N, N, N).
    """
    N = u.shape[1]
    # Wavenumber array for a grid of length L_BOX
    kx = np.fft.fftfreq(N, d=L_BOX / (2 * np.pi * N))
    ky = np.fft.fftfreq(N, d=L_BOX / (2 * np.pi * N))
    kz = np.fft.fftfreq(N, d=L_BOX / (2 * np.pi * N))
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

    grad_sq = np.zeros((N, N, N), dtype=float)
    for c in range(3):
        uh = np.fft.fftn(u[c])
        for K in (KX, KY, KZ):
            # spectral derivative: multiply by i*k then IFFT
            du = np.real(np.fft.ifftn(1j * K * uh))
            grad_sq += du ** 2
    return grad_sq


def _ball_avg_grad_sq(u: np.ndarray, r_frac: float,
                      x0_frac: Tuple[float, float, float] = (0.5, 0.5, 0.5)) -> float:
    """Average |grad u|^2 over ball B_{r_frac}(x0_frac) using FFT gradients.

    Parameters
    ----------
    u       : shape (3, N, N, N)
    r_frac  : radius as fraction of L_BOX
    x0_frac : centre as fraction of L_BOX in each dimension

    Returns
    -------
    ball_avg : float — (1/|B_r|) integral_{B_r} |grad u|^2 dx
    """
    N = u.shape[1]
    cx = int(round(x0_frac[0] * N)) % N
    cy = int(round(x0_frac[1] * N)) % N
    cz = int(round(x0_frac[2] * N)) % N

    grad_sq = _fft_grad_mag_sq(u)
    mask = ball_mask(N, cx, cy, cz, r_frac)
    n_pts = int(mask.sum())
    if n_pts == 0:
        return 0.0
    return float(grad_sq[mask].mean())


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class GammaAdversarialSearch:
    """Build adversarial ICs that concentrate gradient energy and test NS suppression."""

    def __init__(
        self,
        N: int = 32,
        nu: float = 1e-3,
        T: float = 1.0,
        n_steps: int = 100,
        seed: int = 42,
    ) -> None:
        self.N = N
        self.nu = nu
        self.T = T
        self.n_steps = n_steps
        self.seed = seed
        self._rng = np.random.default_rng(seed)

    # ------------------------------------------------------------------
    # Method 1: make_concentrated_ic
    # ------------------------------------------------------------------

    def make_concentrated_ic(
        self, r_frac: float, amplitude: float = 1.0
    ) -> np.ndarray:
        """Divergence-free velocity field with gradient energy concentrated in B_{r_frac}.

        Uses a Gaussian bump psi(x) = exp(-|x - x0|^2 / (2 sigma^2)):
            u[0] =  amplitude * d_y psi
            u[1] = -amplitude * d_x psi
            u[2] = 0
        Then projects to divergence-free subspace (Leray projector in Fourier space).

        Parameters
        ----------
        r_frac    : concentration radius as fraction of L_BOX
        amplitude : overall amplitude scaling

        Returns
        -------
        u : shape (3, N, N, N)
        """
        N = self.N
        sigma = r_frac * L_BOX          # Gaussian width = concentration radius
        dx = L_BOX / N

        # Grid coordinates (periodic box [0, L_BOX]^3)
        x = np.linspace(0.0, L_BOX, N, endpoint=False)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

        # Centre at box midpoint
        x0 = L_BOX / 2.0

        # Periodic minimum-image distance squared
        def _pbc_dist_sq(V: np.ndarray) -> np.ndarray:
            dV = V - x0
            # Wrap into (-L/2, L/2]
            dV = dV - L_BOX * np.round(dV / L_BOX)
            return dV ** 2

        r_sq = _pbc_dist_sq(X) + _pbc_dist_sq(Y) + _pbc_dist_sq(Z)
        psi = np.exp(-r_sq / (2.0 * sigma ** 2))

        # Spectral derivatives of psi
        kx_arr = np.fft.fftfreq(N, d=dx) * (2 * np.pi)
        ky_arr = np.fft.fftfreq(N, d=dx) * (2 * np.pi)
        kz_arr = np.fft.fftfreq(N, d=dx) * (2 * np.pi)
        KX, KY, KZ = np.meshgrid(kx_arr, ky_arr, kz_arr, indexing="ij")

        psi_hat = np.fft.fftn(psi)
        dpsi_dx = np.real(np.fft.ifftn(1j * KX * psi_hat))
        dpsi_dy = np.real(np.fft.ifftn(1j * KY * psi_hat))

        u = np.zeros((3, N, N, N), dtype=float)
        u[0] =  amplitude * dpsi_dy    # curl to ensure div u ~ 0
        u[1] = -amplitude * dpsi_dx
        u[2] = np.zeros((N, N, N))

        # Leray projection: remove curl-free (gradient) part in Fourier space
        u = self._leray_project(u)
        return u

    def _leray_project(self, u: np.ndarray) -> np.ndarray:
        """Project u onto divergence-free subspace via spectral Leray projector.

        P u = u - grad(Delta^{-1} div u)

        Returns new array, does not mutate input.
        """
        N = self.N
        # Wavenumber convention: k_j = 2*pi*j/L_BOX  (radians / metre)
        # fftfreq(N, d=L_BOX/(2*pi*N)) achieves this correctly.
        kx_arr = np.fft.fftfreq(N, d=L_BOX / (2 * np.pi * N))
        ky_arr = np.fft.fftfreq(N, d=L_BOX / (2 * np.pi * N))
        kz_arr = np.fft.fftfreq(N, d=L_BOX / (2 * np.pi * N))
        KX, KY, KZ = np.meshgrid(kx_arr, ky_arr, kz_arr, indexing="ij")

        K2 = KX ** 2 + KY ** 2 + KZ ** 2
        K2[0, 0, 0] = 1.0   # avoid division by zero; zero mode forced to 0 below

        # Compute divergence in Fourier space
        uhat = np.array([np.fft.fftn(u[c]) for c in range(3)])
        div_hat = 1j * (KX * uhat[0] + KY * uhat[1] + KZ * uhat[2])

        # Subtract gradient of pressure: phi_hat = div_hat / (-K2)
        phi_hat = div_hat / (-K2)

        proj_uhat = np.zeros_like(uhat)
        proj_uhat[0] = uhat[0] - 1j * KX * (-phi_hat)
        proj_uhat[1] = uhat[1] - 1j * KY * (-phi_hat)
        proj_uhat[2] = uhat[2] - 1j * KZ * (-phi_hat)

        # Zero the mean mode (enforce zero net flow)
        proj_uhat[:, 0, 0, 0] = 0.0

        # Zero Nyquist planes to enforce Hermitian symmetry
        proj_uhat[:, N // 2, :, :] = 0.0
        proj_uhat[:, :, N // 2, :] = 0.0
        proj_uhat[:, :, :, N // 2] = 0.0

        return np.array([np.real(np.fft.ifftn(proj_uhat[c])) for c in range(3)])

    # ------------------------------------------------------------------
    # Method 2: compute_gamma
    # ------------------------------------------------------------------

    def compute_gamma(
        self,
        u_snapshots: List[np.ndarray],
        r_frac: float,
        x0_frac: Tuple[float, float, float] = (0.5, 0.5, 0.5),
    ) -> float:
        """Compute Gamma(r) = integral_0^T (fint_{B_r} |grad u|^2)^{1/2} dt.

        Parameters
        ----------
        u_snapshots : list of (3, N, N, N) arrays at uniform time intervals
        r_frac      : radius as fraction of L_BOX
        x0_frac     : centre of ball as fraction of L_BOX

        Returns
        -------
        Gamma : float (time-integrated local gradient seminorm)
        """
        if len(u_snapshots) == 0:
            return 0.0

        n = len(u_snapshots)
        dt = self.T / max(n - 1, 1)

        integrand = np.zeros(n, dtype=float)
        for i, u_t in enumerate(u_snapshots):
            ball_avg = _ball_avg_grad_sq(u_t, r_frac, x0_frac)
            integrand[i] = float(np.sqrt(max(ball_avg, 0.0)))

        # Trapezoidal rule
        gamma = float(np.trapezoid(integrand, dx=dt))
        return gamma

    # ------------------------------------------------------------------
    # Method 3: run_free_decay
    # ------------------------------------------------------------------

    def run_free_decay(self, u0: np.ndarray) -> List[np.ndarray]:
        """Simulate free viscous decay: d_t u = nu * Delta u (heat equation).

        Exact in Fourier space: uhat(k, t) = uhat(k, 0) * exp(-nu |k|^2 t).
        The nonlinear term is omitted (conservative estimate — nonlinearity
        typically adds more dissipation).

        Parameters
        ----------
        u0 : shape (3, N, N, N) initial condition

        Returns
        -------
        snapshots : list of (3, N, N, N) arrays at n_steps+1 times in [0, T]
        """
        N = self.N
        kx_arr = np.fft.fftfreq(N, d=L_BOX / N) * (2 * np.pi * N / L_BOX)
        ky_arr = np.fft.fftfreq(N, d=L_BOX / N) * (2 * np.pi * N / L_BOX)
        kz_arr = np.fft.fftfreq(N, d=L_BOX / N) * (2 * np.pi * N / L_BOX)
        KX, KY, KZ = np.meshgrid(kx_arr, ky_arr, kz_arr, indexing="ij")
        K2 = KX ** 2 + KY ** 2 + KZ ** 2   # |k|^2

        # Initial Fourier coefficients
        u0hat = np.array([np.fft.fftn(u0[c]) for c in range(3)])

        times = np.linspace(0.0, self.T, self.n_steps + 1)
        snapshots: List[np.ndarray] = []

        for t in times:
            decay = np.exp(-self.nu * K2 * t)
            u_t = np.array([
                np.real(np.fft.ifftn(u0hat[c] * decay))
                for c in range(3)
            ])
            snapshots.append(u_t)

        return snapshots

    # ------------------------------------------------------------------
    # Method 4: run_search
    # ------------------------------------------------------------------

    def run_search(
        self,
        r_fracs: Optional[List[float]] = None,
    ) -> dict:
        """Adversarial Gamma(r) search across multiple concentration radii.

        For each r in r_fracs:
          - Construct adversarial IC concentrated at scale r
          - Compute Gamma_initial(r) from IC alone
          - Run free decay -> get snapshots
          - Compute Gamma_evolved(r) from evolved snapshots
        Then fit Gamma vs r to power law Gamma ~ C * r^gamma and record verdict.

        Parameters
        ----------
        r_fracs : list of concentration radii as fractions of L_BOX
                  Default: [0.05, 0.1, 0.2, 0.3, 0.4]

        Returns
        -------
        dict with keys:
          gamma_initial  : power-law exponent of Gamma_IC(r)
          gamma_evolved  : power-law exponent of Gamma_NS(r)
          suppressed     : bool (True if gamma_evolved > gamma_initial)
          verdict        : str
          r_fracs        : list of r values tested
          Gamma_IC       : list of Gamma values for IC
          Gamma_NS       : list of Gamma values after decay
          details        : per-r dict
        """
        if r_fracs is None:
            r_fracs = [0.05, 0.1, 0.2, 0.3, 0.4]

        Gamma_IC: List[float] = []
        Gamma_NS: List[float] = []
        details: dict = {}

        for r in r_fracs:
            # Adversarial IC concentrated at scale r; amplitude ~ r^{-1.5}
            # to maximally concentrate gradient energy (energy-class limit)
            amplitude = 1.0 / max(r, 1e-6) ** 1.5
            # Cap amplitude to avoid numerical overflow
            amplitude = min(amplitude, 1e6)

            u0 = self.make_concentrated_ic(r_frac=r, amplitude=amplitude)

            # Gamma from IC alone (single snapshot)
            g_ic = self.compute_gamma([u0], r_frac=r)
            Gamma_IC.append(g_ic)

            # Run free decay and compute Gamma from all snapshots
            snapshots = self.run_free_decay(u0)
            g_ns = self.compute_gamma(snapshots, r_frac=r)
            Gamma_NS.append(g_ns)

            details[r] = {"Gamma_IC": g_ic, "Gamma_NS": g_ns}

        # Fit power laws Gamma ~ C * r^gamma
        gamma_initial = self._fit_power_law(r_fracs, Gamma_IC)
        gamma_evolved = self._fit_power_law(r_fracs, Gamma_NS)

        # Suppression: gamma_evolved > gamma_initial means faster decay with r
        # i.e., NS dynamics push Gamma closer to zero as r -> 0
        suppressed = bool(gamma_evolved > gamma_initial)

        if suppressed:
            verdict = (
                "SUPPRESSED: NS evolution increases Gamma exponent "
                f"({gamma_initial:.3f} -> {gamma_evolved:.3f}). "
                "Claim A conjecture supported — no adverse scaling found."
            )
        else:
            verdict = (
                "NOT SUPPRESSED: NS evolution does not increase Gamma exponent "
                f"({gamma_initial:.3f} -> {gamma_evolved:.3f}). "
                "Possible gap in Claim A conjecture for adversarial ICs."
            )

        return {
            "gamma_initial": float(gamma_initial),
            "gamma_evolved": float(gamma_evolved),
            "suppressed":    suppressed,
            "verdict":       verdict,
            "r_fracs":       list(r_fracs),
            "Gamma_IC":      Gamma_IC,
            "Gamma_NS":      Gamma_NS,
            "details":       details,
        }

    # ------------------------------------------------------------------
    # Internal: power-law fit
    # ------------------------------------------------------------------

    def _fit_power_law(
        self, r_vals: List[float], gamma_vals: List[float]
    ) -> float:
        """Fit Gamma(r) ~ C * r^exponent in log-log space.

        Returns exponent (positive = decays to 0 as r->0 = bounded; negative = diverges).
        """
        r_arr = np.array(r_vals, dtype=float)
        g_arr = np.array(gamma_vals, dtype=float)

        valid = (r_arr > 0.0) & (g_arr > 0.0) & np.isfinite(g_arr)
        if valid.sum() < 2:
            return 0.0

        log_r = np.log(r_arr[valid])
        log_g = np.log(g_arr[valid])

        A = np.column_stack([np.ones(valid.sum()), log_r])
        coeffs, _, _, _ = np.linalg.lstsq(A, log_g, rcond=None)
        return float(coeffs[1])   # exponent


# ---------------------------------------------------------------------------
# Database logging
# ---------------------------------------------------------------------------

def _init_db(conn: sqlite3.Connection) -> None:
    conn.execute(_CREATE_TABLE)
    conn.commit()


def _log_result(conn: sqlite3.Connection, exp_id: str, search: GammaAdversarialSearch,
                results: dict) -> None:
    _init_db(conn)
    import json
    conn.execute(
        _INSERT_RESULT,
        (
            exp_id,
            datetime.now(timezone.utc).isoformat(),
            search.N,
            search.nu,
            search.T,
            search.n_steps,
            search.seed,
            json.dumps(results["r_fracs"]),
            results["gamma_initial"],
            results["gamma_evolved"],
            results.get("gamma_exponent", results["gamma_evolved"]),
            int(results["suppressed"]),
            results["verdict"],
            f"N={search.N}, nu={search.nu}, T={search.T}",
        ),
    )
    conn.commit()


# ---------------------------------------------------------------------------
# Experiment runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    EXP_ID = "EXP-L4-GAM-ADV-001"

    search = GammaAdversarialSearch(N=32, nu=1e-3, T=1.0, n_steps=50, seed=42)
    results = search.run_search(r_fracs=[0.05, 0.1, 0.15, 0.2, 0.3, 0.4])

    print("=" * 60)
    print("ADVERSARIAL Gamma(r) SEARCH — RESULTS")
    print("=" * 60)
    print(f"  Experiment ID:                {EXP_ID}")
    print(f"  Grid: N={search.N}^3, nu={search.nu}, T={search.T}, steps={search.n_steps}")
    print()
    print(f"  gamma_initial (adversarial IC):  {results['gamma_initial']:.3f}")
    print(f"  gamma_evolved (after NS decay):  {results['gamma_evolved']:.3f}")
    print(f"  Suppression observed:             {results['suppressed']}")
    print()
    print("  Per-radius breakdown:")
    for r in results["r_fracs"]:
        d = results["details"][r]
        print(f"    r={r:.2f}: Gamma_IC={d['Gamma_IC']:.4e}, Gamma_NS={d['Gamma_NS']:.4e}")
    print()
    print(f"  Verdict: {results['verdict']}")
    print("=" * 60)

    # Log to DB
    db_path = DB_PATH
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    _log_result(conn, EXP_ID, search, results)
    conn.close()
    print(f"\n  Result logged to {db_path} (table: gamma_adversarial_results)")
    print(f"  claim_id={EXP_ID} | result: {'PASS' if results['suppressed'] else 'PARTIAL'}")
    print(f"  key_metric: gamma_initial={results['gamma_initial']:.3f}, "
          f"gamma_evolved={results['gamma_evolved']:.3f}")
