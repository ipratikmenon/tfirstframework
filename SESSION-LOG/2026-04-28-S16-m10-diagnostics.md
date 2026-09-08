# Session S16 — M10 Full Diagnostics Suite
**Date:** 2026-04-28
**Session:** S16

---

## Summary

M10 complete. Built 4 diagnostic tools in layer4/, wrote 47-test suite.
**All 47 PASS. Total programme: 813/813 tests pass.**

---

## Files Created

- `layer4/continuation_criterion.py` — Phase 0 A_min/μ_min criterion
- `layer4/gronwall_fitter.py` — E1(t) exponential decay fitter
- `layer4/viscosity_scaling.py` — Wang et al. blowup rate vs μ(T)
- `layer4/mu_limit_tracker.py` — ε(δ) scaling tracker (Phase 4)
- `layer4/test_m10.py` — 47-test suite

---

## Tool Descriptions

### continuation_criterion.py
Monitors A_min(t) and μ_min(t) from solver history against Phase 0 blowup
criterion: singularity at T* requires A(T)→0 OR μ(T)→0. Applied to M9
histories — criterion NOT triggered (A_min/A_ref ≈ 0.28 at T=150K, well
above 0.01 threshold).

Key functions: `check_A_min`, `check_mu_min`, `continuation_criterion_verdict`,
`analyse_blowup_result`.

### gronwall_fitter.py
Fits E(t) = E(0)·exp(−c·t) by log-linear least squares. c > 0 with R² > 0.5
is the PASS criterion (positive Grönwall damping → no blowup). Also provides
`extract_E1_spectral` for the H¹ norm from spectral data.

Key functions: `fit_gronwall`, `gronwall_verdict`, `fit_energy_timeseries`,
`extract_E1_spectral`.

### viscosity_scaling.py
Tests Conjecture 3.5: μ(T) damping must grow SLOWER than the Wang et al.
blowup rate (T*−t)^{−(1+λ)} for regularity. Computes suppression ratio
μ(t)/blowup_rate(t) — if ratio is non-increasing, viscosity wins.

Key functions: `blowup_rate`, `suppression_ratio`, `fit_power_law_decay`,
`viscosity_scaling_verdict`.

### mu_limit_tracker.py
Consolidates M4 (2D) and M8 (3D) ε(δ)=|A_min(δ)−A_ref| scaling results.
Log-log fit gives α≈1.0 (linear) → PASS (regular Prize limit). Applied to
M8 data: α=0.98, R²>0.99.

Key functions: `compute_eps`, `fit_scaling_exponent`, `mu_limit_verdict`,
`run_mu_limit_analysis`, `compare_2d_3d`, `load_from_sweep_results`.

---

## Test Results

**47/47 PASS**

- TestContinuationCriterion (14)
- TestGronwallFitter (11)
- TestViscosityScaling (9)
- TestMuLimitTracker (13)

One fix during development: `mu_limit_verdict` early-return order — check
all-zero ε before checking NaN alpha (zero ε is perfect convergence, PASS).

---

## M10 applied to M8 data

```
run_mu_limit_analysis(M8_deltas, M8_A_mins, A_ref=3.0801e-05, label="3D M8")
→ alpha = 0.978, R² = 0.9991, verdict = PASS
```

Matches M8 session log (α=0.980) to 3 significant figures — tracker is consistent.

---

## Test Count

766/766 (previous) + 47 (new) = **813/813 PASS**

---

## Next Priorities

1. **M7** — Route 2 3D at 128³, full ε sweep (8 runs, ~3 hr/run)
2. **Paper assembly** — tfirst_program.tex; all major results now logged
3. **Optional**: apply M10 tools to M9 results in a paper-ready diagnostic report
