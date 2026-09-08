# TASKS.md — T-First Computational Program (v1.0 FINAL)
# Canonical task breakdown for the numerical verification program (PRD §7).
# PROGRESS.md tracks completion status. TASKS.md defines what needs to be done.
# Updated: 2026-04-25 (Session S07) — reflects PRD v1.0 FINAL: Route 2 PRIMARY, δ-gain central.

# KEY CHANGE FROM v0.4: Route 2 (exact Prize NS + auxiliary θ) is now PRIMARY.
# Central proof target: θ ∈ L²(H^{1+δ}) for any δ > 0 (Lemma 2.5).
# Primary proof mechanism: CZ source structure of traceless strain (Route A, no mixing needed).
# Numerical finding: δ ≥ 0.5 observed at λ_max ~ 0 — CZ mechanism confirmed.
# Route 1 (incomp. NS + μ(T)) remains as INDEPENDENT BACKUP.

---

## Computational Architecture Overview (PRD v1.0 FINAL §7)

Route 2 is PRIMARY. Route 1 is independent backup. Four layers.

```
Route 2 (PRIMARY):  Exact Prize NS + auxiliary scalar θ → measure δ > 0
Route 1 (BACKUP):   Incompressible NS + μ(T) + scalar T → A(T) > 0 mechanism
Route 3 (COMPLETE): Compressible NSF + Ma→0 limit (original, archived)

Layer 1 — Property Engine        layer1/tfirst_props.py
          ↓
Layer 2 — 2D Spectral Solver     T_solver_2D.py, route2_2D.py, LPS_monitor.py,
                                  self_similar_IC.py, stretching_2D.py,
                                  lambda_sweep_2D.py, mu_limit_2D.py
          ↓
Layer 3 — 3D Spectral NS         route1_3D.py (TG baseline), route2_3d.py (PRIMARY),
                                  wang_profiles_3D.py, mu_sweep_3D.py,
                                  blowup_search_3D.py
          ↓
Layer 4 — Diagnostic Tools       delta_extractor.py, cz_integrability.py,
                                  lyapunov_3d.py, lps_monitor_3d.py,
                                  strain_anisotropy.py, nonmixing_test.py,
                                  enstrophy_tracker.py, blowup_detector.py
```

**Central proof target (PRD v1.0 §3.3 Lemma 2.5):**
θ ∈ L²(0,T; H^{1+δ}(R³)) for any δ > 0 ⟹ Prodi-Serrin ⟹ global smooth

**Primary proof route (Route A — CZ source structure):**
div-free, traceless strain S → ν|∇u|² ∈ L^{1+ε} → θ ∈ L²(H^{1+δ}) → Lemma 2.5 → Prize

**Key numerical evidence:** δ ≥ 0.5 observed at λ_max ~ 0 in 2D Route 2 experiments.
Gain does NOT require mixing — it is algebraic (CZ structure of the NS nonlinearity).

**Every experiment must:**
1. Have an experiment ID in format `EXP-L{N}-{FLUID}-{NNN}` or `EXP-L{N}-R{route}-{type}-{NNN}`
2. Log claim_id, result (PASS/FAIL/PARTIAL), key_metric, figure_path to results.db
3. Be reproducible (seed + parameters saved)

---

## Session Workflow

Every milestone follows this pattern:

```
Step 1 — SPEC   Read the PRD section targeting this milestone
Step 2 — BUILD  Implement the module; every function has docstring with claim reference
Step 3 — TEST   pytest passes; smoke test with canonical parameters
Step 4 — LOG    First real experiment run; result logged to results.db
Step 5 — UPDATE Update PROGRESS.md; mark milestone complete
```

---

## M0 — Route 1 2D Solver ✅ COMPLETE (Session S03, Apr 19 2026)

**Files:** `layer2/T_solver_2D.py`, `layer2/LPS_monitor.py`, `layer2/stretching_2D.py`
**Tests:** `layer2/test_M3.py` — 64/64 PASS

### What was built

- [x] `T_solver_2D.py` — Route 1 2D: vorticity-stream function (incompressible, div u=0 exact, ρ=const)
  - Spectral RK4, 2/3 dealiasing, adaptive CFL
  - T equation: ∂_t T + u·∇T = A(T)·ΔT + Q_visc
  - All properties via layer1/tfirst_props.py — no inline formulas
  - `run_lambda_sweep()` function covers all 6 Wang et al. λ values
- [x] `LPS_monitor.py` — tracks Z, P, ‖u‖_{Lp}, A_min, suppression ratio
- [x] `stretching_2D.py` — D_Z, N_P, D/S ratio, enstrophy spectrum
- [x] `test_M3.py` — 64/64 PASS (solver construction, RHS, timestepping, run loop,
      LPS monitor, stretching diagnostics, lambda sweep smoke test)

---

## M1 — Property Engine ✅ COMPLETE (Session S01, Apr 18 2026) — 2 additions pending

**File:** `layer1/tfirst_props.py`  ·  **Tests:** `layer1/test_props.py`  ·  **52/52 PASS**

### Completed

- [x] `ideal_props(T)` — Sutherland viscosity + power-law k
- [x] `co2_props(T, P)` — CoolProp SC-CO2
- [x] `A_field(T_arr)` — thermal diffusivity k/(ρ·cv), assert > 0
- [x] `second_law_check()` — ideal + CO2, 5 pressures
- [x] `self_similar_props(λ, t)` — Conjecture 3.5 suppression ratio
- [x] `conjecture_34_sweep()` — all 6 Wang et al. λ values, 18/18 PASS
- [x] Results logged: EXP-L1-IDEAL-001, EXP-L1-CO2-001

### [v0.4 NEW — pending]

- [ ] `route1_coeffs(T_field, fluid='ideal', P=None) → dict`
  - Returns `{mu_min, mu_max, A_min, A_max, nu_eff_min, nu_eff_max}`
  - Asserts mu_min > 0, A_min > 0 (both are thermodynamic theorems)
  - Used by LPS_monitor to compute margin bounds at each timestep
  - PRD §7.2 v0.4: coefficient bounds needed for Route 1 continuation criterion
- [ ] `route2_theta_source(u, v, nu, kx, ky) → np.ndarray`
  - Computes θ viscous source term: S_θ = ν·|∇u|² (spectral accuracy)
  - θ equation: ∂_t θ + u·∇θ = ν·Δθ + S_θ
  - This is the key coupling — exact Prize NS generates θ automatically
  - Used by route2_2D.py; in layer1 because it depends only on ν (constant) + velocity
- [ ] Add tests for both new functions to `layer1/test_props.py`

**M1 v0.4 complete when:** 2 new functions added + tests pass.

---

## IC Infrastructure ✅ COMPLETE (Session S02, Apr 19 2026)

**File:** `layer2/self_similar_IC.py`  ·  **Tests:** `layer2/test_self_similar_IC.py`  ·  **53/53 PASS**

All 7 Wang et al. profiles + adversarial_min. Used by M0, M3, M5, M6.

---

## M2 — Route 2 2D: Exact Prize NS + Auxiliary Scalar θ (Month 2 — 128²)

**Target:** `layer2/route2_2D.py`
**PRD reference:** §7.2 Route 2, §3.2
**Pass criterion:** LPS margin ε(ε_param) > 0 for all ε_param ∈ {1, 0.1, 0.01, 0.001}.

Route 2 works with exact Prize equations (ν=const, div u=0, ρ=const) and introduces
auxiliary scalar θ that evolves alongside u. No limit is needed — this directly
answers Tao's supercriticality objection.

### Tasks

- [ ] `layer2/route2_2D.py`
  - Exact Prize NS: ∂_t u + (u·∇)u = −∇p + ν·Δu; div u = 0; ρ=const
  - θ equation: ∂_t θ + u·∇θ = ν·Δθ + ν·|∇u|²
  - μ_eff = ν + ε_param · f(θ); f smooth, non-negative
  - Spectral RK4; 2/3 dealiasing; same grid setup as T_solver_2D.py
  - θ source term via `layer1.tfirst_props.route2_theta_source()`
  - Run 4 values: ε_param ∈ {1.0, 0.1, 0.01, 0.001}
  - Measure LPS margin ε at end of each run
- [ ] Tests: `layer2/test_route2_2D.py`
  - θ source term is non-negative (ν·|∇u|² ≥ 0)
  - θ monotone increasing (viscous source only adds)
  - LPS margin positive for all ε_param values
  - ε_param→0 limit: margin doesn't collapse
- [ ] Results logged: EXP-L2-R2-001 through EXP-L2-R2-004
- [ ] Key result: does ε(ε_param) plateau > 0 as ε_param → 0?

**M2 complete when:** All 4 ε_param runs done; LPS margin verified > 0; results in results.db.

---

## M3 — Lambda Sweep 2D (Month 3 — 64²×12 runs)

**Target:** `layer2/lambda_sweep_2D.py`
**PRD reference:** §7.3 — key Conjecture 3.5 numerical test
**Pass criterion:** All λ values show suppression (no λ_c found). PASS = Conjecture 3.5 confirmed.

Route 1 solver already built (M0). Need the sweep driver with results.db logging.

### Tasks

- [ ] `layer2/lambda_sweep_2D.py`
  - Import `run_lambda_sweep()` from T_solver_2D.py
  - Add results.db logging for each run (claim_id, key_metric, verdict)
  - Experiment IDs: EXP-L2-R1-CCF-001 … EXP-L2-R1-CCF-006
  - Generate comparison figure: A_min(t) and D/S ratio per λ
- [ ] Run all 6 λ values: {1.2, 0.9, 0.6057, 0.4703, 0.3, 0.1}
  - Grid: 64², t_end=0.5, fluid='ideal'
  - Smoke test first: N=32, t_end=0.05
- [ ] Boussinesq profiles: same λ sweep (EXP-L2-R1-BOUS-001…003)
- [ ] Adversarial min IC: λ→0 (EXP-L2-R1-ADV-001)
- [ ] Conjecture 3.5 verdict: PASS = all λ suppressed, no λ_c found
- [ ] Add smoke-test run to CI (N=32, 3 λ values, fast)

**M3 complete when:** All 12+ runs complete; conjecture verdict logged; comparison figure saved.

---

## M4 — μ(T)→ν Limit 2D (Month 4 — 128²×10 runs)

**Target:** `layer2/mu_limit_2D.py`
**PRD reference:** §7.4 Route 1 Phase 4 test
**Pass criterion:** LPS margin ε(δ) remains > 0 as δ → 0.

This is the 2D version of Phase 4: as T₀ → uniform T̄ (i.e., δ → 0),
μ(T) → const ≡ ν. Does the LPS margin survive?

### Tasks

- [ ] `layer2/mu_limit_2D.py`
  - IC: T₀(x) = T̄ + δ · f(x), with f a smooth test function
  - 10 δ values: {1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001}
  - Measure LPS margin ε at t=t_final for each δ
  - Fit ε(δ) vs δ: plateau (constant limit exists) or collapse?
- [ ] Results: EXP-L2-R1-MULIMIT-001 … EXP-L2-R1-MULIMIT-010
- [ ] Phase 4 Prize deliverable: ε(μ) > 0 as μ → ν
- [ ] Figure: ε(δ) vs δ on log-log scale

**M4 complete when:** All 10 δ runs done; ε(δ) curve shows plateau > 0.

---

## M5 — 3D Baseline + Route 2 3D (PRD v1.0 §7.3)

**Target:** `layer3/route1_3D.py` (TG baseline) + `layer3/route2_3d.py` (PRIMARY deliverable)
**PRD reference:** §7.3 Steps 1–7
**Pass criterion:** 
- Step 1: TG benchmark energy decay within 5% of reference
- Steps 2–7: δ(t) > 0 in 3D, ε_LPS > 0 at ε=0 in 3D, δ > 0 in non-mixing IC

### Step 1 — 3D Baseline (Route 1 TG) [IN PROGRESS — bug B-005 open]

- [x] `layer3/route1_3D.py` — full 3D incompressible NS + scalar T, RK4, dealiased
- [ ] **Fix `project_divergence_free` bug** — Leray projector not zeroing div for random fields
  - Current: div_after=1.465 (should be ~0)
  - 54/56 tests pass; TestDivFree (2) failing
- [ ] All 56/56 tests pass in `layer3/test_route1_3D.py`
- [ ] TG benchmark run at N=32 → EXP-L3-R1-TG-001 (smoke)
- [ ] TG benchmark run at N=64 → EXP-L3-R1-TG-002 (full)

### Step 2 — Theta Equation in 3D

- [ ] `layer3/route2_3d.py` — Route 2 in 3D
  - Exact Prize NS: ∂_t u + (u·∇)u = −∇p + ν·Δu + ∇·(ε·f(θ)·∇u); div u = 0
  - θ equation with IMEX: implicit ν·Δθ, explicit advection and source
  - Source: ν·Σᵢⱼ(∂ᵢuⱼ)² = 2ν·S:S (full strain tensor squared, traceless + div-free)
  - Verify: θ_min ≥ −1e-10 always; θ grows tracking viscous dissipation
- [ ] Tests: `layer3/test_route2_3d.py` — θ ≥ 0, θ source ≥ 0, energy decay, div-free

### Step 3 — δ Extraction (KEY DIAGNOSTIC — PRD §7.3 Step 3)

- [ ] `layer4/delta_extractor.py`
  - Compute ‖θ(t)‖_{H^s} for s = 1.0, 1.05, 1.1, 1.2, 1.3, 1.5 at each saved step
  - Fit δ: largest s−1 where ∫₀ᵀ ‖θ‖²_{H^{1+δ}} dt stays bounded
  - Track δ(t) as timeseries — stable positive δ confirms Route A
  - **Pass criterion: δ(t) > 0 throughout run**

### Step 4 — CZ Structure Verification (CLAIM A TEST)

- [ ] `layer4/cz_integrability.py`
  - Compute ‖ν|∇u|²‖_{L^{1+ε}} directly from u field at each step
  - Test whether ε > 0 — this is the key Claim A to verify numerically
  - Compare: does ε correlate with traceless strain anisotropy?
  - **Pass criterion: ε > 0 even in non-mixing (shear) flows**

### Step 5 — Lyapunov Exponent in 3D

- [ ] `layer4/lyapunov_3d.py`
  - Place 200 tracer particles, track pair separations, fit λ_max(t)
  - Critical: is δ > 0 in runs where λ_max ~ 0? (Tests Route A vs Route B)
  - **Pass criterion: δ > 0 with λ_max = 0 → Route A confirmed**

### Step 6 — Non-Mixing Control Experiment (DECISIVE)

- [ ] Run with shear flow IC: u₀ = (f(y), 0, 0) — explicitly non-mixing, λ_max = 0
  - Measure δ and ε_LPS in this run
  - **Pass criterion: δ > 0 even without mixing → Route A is primary**

### Step 7 — LPS Margin at ε=0 in 3D (MOST PRIZE-RELEVANT)

- [ ] Run exact Prize NS (ε=0) in 3D; measure ε_LPS(t)
  - Extends 2D result: ε_LPS > 0 (min=0.155, mean=0.162)
  - **Pass criterion: ε_LPS > 0 throughout 3D run at ε=0**
  - EXP-L3-R2-LPS-001

**M5 complete when:** Steps 1–7 all PASS; δ(t) > 0 and ε_LPS > 0 logged to results.db.

---

## M6 — Wang et al. 3D (Month 7 — 128³×6 runs)

**Target:** `layer3/wang_profiles_3D.py`
**PRD reference:** §7.3 3D lambda sweep
**Pass criterion:** CCF λ=0.4703 (2nd unstable) — T-first suppresses blow-up in 3D.

### Tasks

- [ ] `layer3/wang_profiles_3D.py`
  - Extends 2D Wang et al. profiles to 3D initial conditions
  - All 6 λ values + adversarial min
- [ ] Lambda sweep 3D: all 6 λ × 128³, results logged
- [ ] `layer4/viscosity_scaling.py` — measures μ(T) along profile vs (1−t)^{-(1+λ)}
- [ ] Key result: CCF λ=0.4703 in 3D — PASS or FAIL?
- [ ] Logs: EXP-L3-R1-CCF-001 … EXP-L3-R1-CCF-006

**M6 complete when:** All 3D lambda sweep runs complete; 3D Conjecture 3.5 verdict.

---

## M7 — Route 2 3D (Month 8 — 128³×8 runs)

**Target:** `layer3/route2_3D.py`
**PRD reference:** §7.2 Layer 3 Route 2
**Pass criterion:** LPS margin positive for all ε_param values in 3D.

### Tasks

- [ ] `layer3/route2_3D.py` — exact Prize NS + θ in 3D; ε_param sweep (same as M2 in 3D)
- [ ] 8 runs (4 ε_param × 2 ICs); logs EXP-L3-R2-001 … EXP-L3-R2-008

**M7 complete when:** 3D Route 2 LPS margin verified > 0.

---

## M8 — μ Limit 3D (Month 10 — 64³×15 runs)

**Target:** `layer3/mu_sweep_3D.py`
**PRD reference:** §7.4 Phase 4 3D
**Pass criterion:** ε(δ) plateau > 0 in 3D (Phase 4 Prize deliverable).

### Tasks

- [ ] `layer3/mu_sweep_3D.py` — Phase 4 δ sweep in 3D (15 δ values)
- [ ] `layer4/mu_limit_tracker.py` — fits ε(δ) vs δ, identifies plateau

**M8 complete when:** 3D ε(δ) curve shows plateau above zero.

---

## M9 — Blowup Search (Month 12 — 256³)

**Target:** `layer3/blowup_search_3D.py`
**PRD reference:** §7.5
**Pass criterion:** No singularity found (strongest evidence for T-first).

### Tasks

- [ ] `layer3/blowup_search_3D.py`
  - Adversarial IC: max vortex stretching, T₀ = T_min + δ
  - Wang et al. CCF and Boussinesq unstable profiles
  - 256³ grid; ~12hr; ~$200 cloud
- [ ] `layer4/blowup_detector.py` — automated norm monitoring + alarm
- [ ] Results logged; paper-quality figures

**M9 complete when:** Blowup search complete; verdict published in results.db.

---

## M10 — Full Diagnostics (Month 14)

**Target:** All Layer 4 tools complete and tested.

### Tasks

- [ ] `layer4/lps_classifier.py` — velocity field → (p,q) Prodi-Serrin pair
- [ ] `layer4/gronwall_fitter.py` — E1(t) → damping constant c
- [ ] `layer4/continuation_criterion.py` — tracks distance to continuation failure
- [ ] Automated result pipeline:
  - `run_experiment(route, fluid, grid, IC_type, lambda_val, duration)` wrapper
  - Calls all diagnostics; generates LaTeX table row; appends to results.db
- [ ] Full pipeline smoke test end-to-end

**M10 complete when:** Any experiment auto-generates a paper-ready table row.

---

## Publication Targets (PRD v1.0 FINAL §6)

| Paper | Content | Target journal | Month |
|---|---|---|---|
| Note 1 | Route 2 2D verification: exact Prize NS + θ, δ≥0.5 | arXiv | 1 |
| Paper 1 | Scalar θ in exact NS: max principle, Lemma 2.5 reduction | Arch. Rat. Mech. Anal. | 8 |
| Paper 2 | δ-gain closes LPS: any H^{1+δ} improvement suffices | Comm. Math. Phys. | 10 |
| Paper 3 | δ > 0 via CZ source structure — unconditional (Claim A) | Invent. Math. | 20 |
| Paper 4 | Mixing-enhanced δ > 0 — conditional result (Route B) | Comm. Pure Appl. Math. | 18 |
| Paper 5 | Global regularity for Route 2 NS — assembly | Ann. Math. | 30 |
| Paper 6 | Removing ε: global regularity of exact Prize NS | Ann. Math. | 42 |
| Paper 7 | Route 1: global regularity for incompressible NS with μ(T) | Invent. Math. | 28 |
| Paper 8 | SC-CO2 turbulence: T-first engineering predictions | J. Fluid Mech. | 14 |

---

## Computational Resource Planning (PRD v1.0 §7.2 — JAX Metal on M5 Pro)

| Milestone | Grid | NumPy time | JAX Metal time | Cloud cost |
|---|---|---|---|---|
| M0–M1 | up to 64² | Minutes | Seconds | $0 |
| M2–M4 | 64²–128² | Minutes | Minutes | $0 |
| M5 Step 1 (TG) | 64³ dev | ~15 min | ~2 min | $0 |
| M5 Steps 2–7 | 64³ | ~15 min/run | ~2 min/run | $0 |
| M6 Wang et al. 3D | 128³ × 8 | ~4 hr/run | ~15-30 min/run | $0 |
| M7 LPS 3D | 128³ × 8 | ~4 hr/run | ~15-30 min/run | $0 |
| M8 μ limit 3D | 64³ × 15 | ~1 hr/run | ~5 min/run | $0 |
| M9 Blowup search | 256³ | Overnight | ~4-6 hr | ~$200 AWS spot |
| M10 Diagnostics | N/A | Diagnostic only | N/A | $0 |
