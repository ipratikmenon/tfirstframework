"""
test_layer4.py — Tests for Layer 4 diagnostic tools
=====================================================
Covers: delta_extractor, cz_integrability, lps_monitor_3d

Run:
  cd /Users/pratikmenon/Documents/Claude/tfirst
  python -m pytest layer4/test_layer4.py -v
"""

from __future__ import annotations

import numpy as np
import pytest
import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.join(_HERE, "..")
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_ROOT, "layer3"))
sys.path.insert(0, _ROOT)

import delta_extractor as de
import cz_integrability as cz
import lps_monitor_3d as lps

N = 16          # small grid for all tests — fast
NU = 1.0e-3
EPS = 0.1
T_END = 0.02    # short runs — enough to seed source; not slow


# =============================================================================
# delta_extractor tests
# =============================================================================

class TestDeltaExtractor:
    @pytest.fixture(scope="class")
    def tg_result(self):
        return de.extract_delta_timeseries(
            N=N, nu=NU, eps_param=EPS, ic="tg",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )

    def test_returns_dict(self, tg_result):
        assert isinstance(tg_result, dict)
        assert "summary" in tg_result
        assert "timeseries" in tg_result
        assert "verdict" in tg_result

    def test_timeseries_nonempty(self, tg_result):
        assert len(tg_result["timeseries"]) > 0

    def test_verdict_pass_for_tg(self, tg_result):
        assert tg_result["verdict"] == "PASS"

    def test_delta_max_nonneg(self, tg_result):
        assert tg_result["summary"]["delta_max"] >= 0.0

    def test_delta_positive_frac_in_range(self, tg_result):
        frac = tg_result["summary"]["delta_positive_frac"]
        assert 0.0 <= frac <= 1.0

    def test_timeseries_has_required_keys(self, tg_result):
        row = tg_result["timeseries"][0]
        for key in ("t", "dt", "delta_est", "E", "theta_max", "H0.0", "H1.0"):
            assert key in row, f"Missing key: '{key}'"

    def test_time_monotone_increasing(self, tg_result):
        ts = [r["t"] for r in tg_result["timeseries"]]
        for i in range(len(ts) - 1):
            assert ts[i+1] > ts[i], f"Time not increasing at step {i}"

    def test_sobolev_norms_nonneg(self, tg_result):
        for row in tg_result["timeseries"]:
            for s in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5):
                key = f"H{s:.1f}"
                assert row[key] >= 0.0, f"{key} < 0 at t={row['t']:.4f}"

    def test_shear_ic_returns_valid(self):
        r = de.extract_delta_timeseries(
            N=N, nu=NU, eps_param=EPS, ic="shear",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )
        assert r["verdict"] == "PASS"
        assert r["summary"]["delta_max"] >= 0.0

    def test_random_ic_returns_valid(self):
        r = de.extract_delta_timeseries(
            N=N, nu=NU, eps_param=EPS, ic="random",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )
        assert r["verdict"] == "PASS"

    def test_eps_zero_still_passes(self):
        """ε=0 (exact Prize) must also return PASS — central test."""
        r = de.extract_delta_timeseries(
            N=N, nu=NU, eps_param=0.0, ic="tg",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )
        assert r["verdict"] == "PASS"

    def test_unknown_ic_raises(self):
        with pytest.raises(ValueError, match="Unknown ic"):
            de.extract_delta_timeseries(N=N, ic="blah", t_end=T_END)

    def test_delta_eps_sweep_returns_list(self):
        results = de.delta_eps_sweep(
            N=N, nu=NU, ic="tg", t_end=T_END,
            eps_vals=(0.1, 0.0), max_steps=30, verbose=False, db_path=None,
        )
        assert len(results) == 2

    def test_delta_eps_sweep_all_pass(self):
        results = de.delta_eps_sweep(
            N=N, nu=NU, ic="tg", t_end=T_END,
            eps_vals=(0.1, 0.01, 0.0), max_steps=30, verbose=False, db_path=None,
        )
        for r in results:
            assert r["verdict"] == "PASS", \
                f"ε={r['eps_param']} FAIL: delta_max={r['delta_max']}"

    def test_exp_id_format(self):
        r = de.extract_delta_timeseries(
            N=N, nu=NU, ic="tg", t_end=T_END, max_steps=20, db_path=None,
        )
        assert f"EXP-L4-DELTA-TG-{N:03d}" == r["summary"]["exp_id"]

    def test_exp_id_override(self):
        r = de.extract_delta_timeseries(
            N=N, nu=NU, ic="tg", t_end=T_END, max_steps=20, db_path=None,
            exp_id_override="MY-CUSTOM-ID",
        )
        assert r["summary"]["exp_id"] == "MY-CUSTOM-ID"


# =============================================================================
# cz_integrability tests
# =============================================================================

class TestCZNorms:
    def test_lp_norms_nonneg(self):
        """L^p norms must be ≥ 0."""
        rng = np.random.default_rng(1)
        src = np.abs(rng.standard_normal((N, N, N)))
        norms = cz.cz_lp_norms(src)
        for p, v in norms.items():
            assert v >= 0.0, f"L^{p} norm negative: {v}"

    def test_lp_norms_monotone_in_p(self):
        """‖f‖_{L^p} should be monotonically increasing in p for non-constant f."""
        rng = np.random.default_rng(2)
        src = np.abs(rng.standard_normal((N, N, N))) + 0.1
        norms = cz.cz_lp_norms(src, (1.0, 1.1, 1.5, 2.0))
        vals = [norms[p] for p in sorted(norms)]
        for i in range(len(vals) - 1):
            assert vals[i] <= vals[i+1] + 1e-8

    def test_lp_norms_zero_source(self):
        """Zero source → all L^p norms = 0."""
        src = np.zeros((N, N, N))
        norms = cz.cz_lp_norms(src)
        for v in norms.values():
            assert v == 0.0

    def test_eps_cz_zero_for_zero_source(self):
        """ε_CZ = 0 when source is zero (undefined)."""
        norms = {p: 0.0 for p in cz._P_VALS}
        assert cz.estimate_eps_cz(norms) == 0.0

    def test_eps_cz_nonneg(self):
        """ε_CZ must be ≥ 0."""
        rng = np.random.default_rng(3)
        src = np.abs(rng.standard_normal((N, N, N)))
        norms = cz.cz_lp_norms(src)
        eps = cz.estimate_eps_cz(norms)
        assert eps >= 0.0

    def test_eps_cz_uniform_source(self):
        """Uniform source (constant) has L^p / L^1 = 1 → controlled → ε_CZ = max(p-1)."""
        src = np.ones((N, N, N))
        norms = cz.cz_lp_norms(src, (1.0, 1.1, 1.5, 2.0))
        # All ratios = 1.0 < threshold → ε_CZ = 2.0 - 1.0 = 1.0
        eps = cz.estimate_eps_cz(norms, threshold=10.0)
        assert eps == pytest.approx(1.0, abs=1e-10)


class TestCZTracker:
    @pytest.fixture(scope="class")
    def tg_result(self):
        return cz.track_cz_integrability(
            N=N, nu=NU, eps_param=EPS, ic="tg",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )

    def test_returns_dict(self, tg_result):
        assert isinstance(tg_result, dict)
        for key in ("summary", "timeseries", "verdict"):
            assert key in tg_result

    def test_verdict_pass(self, tg_result):
        assert tg_result["verdict"] == "PASS"

    def test_timeseries_nonempty(self, tg_result):
        assert len(tg_result["timeseries"]) > 0

    def test_eps_cz_nonneg_all_steps(self, tg_result):
        for row in tg_result["timeseries"]:
            assert row["eps_cz"] >= 0.0, \
                f"ε_CZ < 0 at t={row['t']:.4f}: {row['eps_cz']}"

    def test_s_max_nonneg(self, tg_result):
        for row in tg_result["timeseries"]:
            assert row["S_max"] >= -1e-12, f"S_max negative: {row['S_max']:.3e}"

    def test_s_mean_nonneg(self, tg_result):
        for row in tg_result["timeseries"]:
            assert row["S_mean"] >= -1e-12

    def test_required_keys_in_summary(self, tg_result):
        s = tg_result["summary"]
        for key in ("eps_cz_max", "eps_cz_min", "eps_cz_mean",
                    "claim_a_frac", "verdict", "key_metric"):
            assert key in s, f"Missing summary key: '{key}'"

    def test_claim_a_frac_in_range(self, tg_result):
        frac = tg_result["summary"]["claim_a_frac"]
        assert 0.0 <= frac <= 1.0

    def test_shear_ic_cz_valid(self):
        r = cz.track_cz_integrability(
            N=N, nu=NU, eps_param=EPS, ic="shear",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )
        assert r["verdict"] == "PASS"

    def test_unknown_ic_raises(self):
        with pytest.raises(ValueError):
            cz.track_cz_integrability(N=N, ic="bad", t_end=T_END)

    def test_exp_id_format(self, tg_result):
        exp_id = tg_result["summary"]["exp_id"]
        assert f"EXP-L4-CZ-TG-{N:03d}" == exp_id


# =============================================================================
# lps_monitor_3d tests
# =============================================================================

class TestLPSNorms:
    def setup_method(self):
        import route2_3D as r2
        x = np.linspace(0, 2*np.pi, N, endpoint=False)
        xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
        self.u = np.sin(xx) * np.cos(yy)
        self.v = -np.cos(xx) * np.sin(yy)
        self.w = np.zeros((N, N, N))
        kx, ky, kz, _ = r2.make_grid_3D(N)
        self.kx, self.ky, self.kz = kx, ky, kz

    def test_spatial_lp_norm_nonneg(self):
        for p in (3, 4, 6):
            val = lps.spatial_lp_norm(self.u, self.v, self.w, p)
            assert val >= 0.0

    def test_spatial_lp_monotone(self):
        """‖f‖_{L^p} should increase with p (for non-constant f)."""
        vals = [lps.spatial_lp_norm(self.u, self.v, self.w, p) for p in (3, 4, 6)]
        for i in range(len(vals) - 1):
            assert vals[i] <= vals[i+1] + 1e-8

    def test_linf_ge_lp(self):
        """‖f‖_{L^∞} ≥ ‖f‖_{L^6} always."""
        linf = lps.linf_norm(self.u, self.v, self.w)
        lp6 = lps.spatial_lp_norm(self.u, self.v, self.w, 6)
        assert linf >= lp6 - 1e-10

    def test_enstrophy_nonneg(self):
        Z = lps.enstrophy(self.u, self.v, self.w, self.kx, self.ky, self.kz)
        assert Z >= 0.0

    def test_enstrophy_zero_velocity(self):
        u = np.zeros((N, N, N))
        Z = lps.enstrophy(u, u, u, self.kx, self.ky, self.kz)
        assert Z < 1e-20

    def test_palinstrophy_nonneg(self):
        P = lps.palinstrophy(self.u, self.v, self.w, self.kx, self.ky, self.kz)
        assert P >= 0.0

    def test_lp_norm_zero_velocity(self):
        u = np.zeros((N, N, N))
        for p in (3, 4, 6):
            val = lps.spatial_lp_norm(u, u, u, p)
            assert val < 1e-20


class TestLPSMonitor:
    @pytest.fixture(scope="class")
    def tg_result(self):
        return lps.monitor_lps(
            N=N, nu=NU, eps_param=EPS, ic="tg",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )

    def test_returns_dict(self, tg_result):
        for key in ("summary", "timeseries", "verdict"):
            assert key in tg_result

    def test_verdict_pass_tg(self, tg_result):
        assert tg_result["verdict"] == "PASS"

    def test_timeseries_nonempty(self, tg_result):
        assert len(tg_result["timeseries"]) > 0

    def test_enstrophy_nonneg_all_steps(self, tg_result):
        for row in tg_result["timeseries"]:
            assert row["Z"] >= 0.0, f"Z < 0 at t={row['t']:.4f}"

    def test_eps_lps_nonneg_all_steps(self, tg_result):
        for row in tg_result["timeseries"]:
            assert row["eps_lps"] >= -1e-12, f"ε_LPS < 0 at t={row['t']:.4f}"

    def test_lp6_nonneg_all_steps(self, tg_result):
        for row in tg_result["timeseries"]:
            assert row["Lp6"] >= 0.0

    def test_div_rms_stays_small(self, tg_result):
        """div u must remain < 1e-8 throughout."""
        for row in tg_result["timeseries"]:
            assert row["div_rms"] < 1e-8, \
                f"div_rms = {row['div_rms']:.2e} at t={row['t']:.4f}"

    def test_summary_required_keys(self, tg_result):
        s = tg_result["summary"]
        for key in ("Z_max", "lp6_max", "eps_lps_min",
                    "blowup_flag", "verdict", "key_metric"):
            assert key in s

    def test_no_blowup_flag_tg(self, tg_result):
        assert not tg_result["summary"]["blowup_flag"]

    def test_eps_lps_min_nonneg_summary(self, tg_result):
        assert tg_result["summary"]["eps_lps_min"] >= -1e-12

    def test_shear_ic_valid(self):
        r = lps.monitor_lps(
            N=N, nu=NU, eps_param=EPS, ic="shear",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )
        assert r["verdict"] == "PASS"

    def test_eps_zero_exact_prize(self):
        """ε=0 (exact Prize equations) should still PASS."""
        r = lps.monitor_lps(
            N=N, nu=NU, eps_param=0.0, ic="tg",
            t_end=T_END, max_steps=50, verbose=False, db_path=None,
        )
        assert r["verdict"] == "PASS"
        assert r["summary"]["eps_lps_min"] == pytest.approx(0.0, abs=1e-12)

    def test_unknown_ic_raises(self):
        with pytest.raises(ValueError):
            lps.monitor_lps(N=N, ic="bad_ic", t_end=T_END)

    def test_exp_id_format(self, tg_result):
        assert f"EXP-L4-LPS-TG-{N:03d}" == tg_result["summary"]["exp_id"]

    def test_time_monotone(self, tg_result):
        ts = [r["t"] for r in tg_result["timeseries"]]
        for i in range(len(ts) - 1):
            assert ts[i+1] > ts[i]


# =============================================================================
# Cross-tool consistency
# =============================================================================

class TestCrossToolConsistency:
    """delta_extractor and lps_monitor should agree on δ_est values."""

    def test_delta_and_lps_agree_on_energy(self):
        """Both tools should report non-growing energy (viscous dissipation)."""
        r_de = de.extract_delta_timeseries(
            N=N, nu=NU, eps_param=EPS, ic="tg",
            t_end=0.1, max_steps=200, verbose=False, db_path=None,
        )
        r_lps = lps.monitor_lps(
            N=N, nu=NU, eps_param=EPS, ic="tg",
            t_end=0.1, max_steps=200, verbose=False, db_path=None,
        )
        # Energy must not grow (viscosity always dissipates)
        E_de  = [r["E"] for r in r_de["timeseries"]]
        E_lps = [r["E"] for r in r_lps["timeseries"]]
        assert E_de[-1] <= E_de[0] + 1e-10, "delta_extractor: energy grew"
        assert E_lps[-1] <= E_lps[0] + 1e-10, "lps_monitor: energy grew"

    def test_cz_and_lps_both_pass_shear(self):
        """CZ tracker and LPS monitor should both PASS for shear IC."""
        r_cz = cz.track_cz_integrability(
            N=N, nu=NU, ic="shear", t_end=T_END, max_steps=30, db_path=None,
        )
        r_lps = lps.monitor_lps(
            N=N, nu=NU, ic="shear", t_end=T_END, max_steps=30, db_path=None,
        )
        assert r_cz["verdict"] == "PASS"
        assert r_lps["verdict"] == "PASS"
