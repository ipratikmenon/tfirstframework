# T-First — Root CLAUDE.md
# Read this file at the start of every session, every worktree, every agent.
# This is the single source of truth for program-wide rules and conventions.
# AUTHORITATIVE PRD: tfirst_prd_v4.docx  (v0.4 — April 2026)

---

## FIRST ACTION EVERY SESSION

Before doing anything else, read PROGRESS.md in the root of this repo.
It tells you exactly where the program is, what was last verified, what is in progress,
and what to do next. Do not proceed with any work until you have read it.

```
cat PROGRESS.md
```

After every completed milestone — no matter how small — update PROGRESS.md before stopping.

At the END of every session, write a session log file:
```
SESSION-LOG/YYYY-MM-DD-SNN-description.md
```
Both PROGRESS.md and the session log must be written before committing.
A session without a log file is an incomplete session.

---

## Program Identity

**T-First** is a mathematical research program attacking the Clay Millennium Prize:
the existence and global smoothness of solutions to the 3D incompressible Navier-Stokes equations.

**Central thesis (v0.4):** Temperature T(x,t) is a **scalar** field — the master variable.
Every other fluid quantity (ρ, p, μ, k, velocity u) is downstream of T. The LPS regularity
gap is a thermodynamic artifact that closes through **A(T) = k(T)/[ρ(T)·cv(T)] > 0**,
a direct consequence of the second law. Because T satisfies a quasilinear parabolic
**scalar** PDE, the **strong maximum principle** applies — a tool unavailable to
vector-based approaches working with u or ω directly.

**v0.4 pivot:** The primary system is now **Route 1 — incompressible NS with μ(T)**,
staying in Prize geometry (div u = 0, ρ = const) from the outset. Compressible NSF
(Route 3) is now an independent completeness argument, not the primary object.

**Primary document:** `tfirst_prd_v4.docx` — read before any architectural decisions.

**Program authors:** Pratik Menon + Claude (Anthropic)

**Current phase:** Numerical Verification — building the computational program (PRD §7).

---

## Three Routes to the Prize

| Route | System | Prize geometry? | Final limit | Status |
|---|---|---|---|---|
| **Route 1** (primary) | Incomp. NS + μ(T) + scalar T | ✅ from start | μ(T) → ν (simple thermal relaxation) | Active — building now |
| **Route 2** (backup) | Exact Prize NS + auxiliary scalar θ | ✅ exact | ε → 0 (elementary) | Planned M2 |
| **Route 3** (completeness) | Compressible NSF | ❌ ρ variable | Ma → 0 (hard) | Previous approach; independent |

### Route 1 — Primary System (PRD §3.2)

Incompressible NS with temperature-dependent viscosity; Prize geometry throughout:

```
ρ·(∂_t u + (u·∇)u) = −∇p + ∇·(μ(T)∇u)     div u = 0
ρ·cv·∂_t T = ∇·(k(T)∇T) − ρ·cv·(u·∇)T + μ(T)|∇u|²
μ(T) = μ₀·(T/T₀)^(3/2)     k(T) = k₀·(T/T₀)^(3/2)     A(T) = k(T)/(ρ·cv) > 0
```

Prize equations recovered by μ(T) → ν (uniform temperature → thermal relaxation). Simple limit.

### Route 2 — Backup System (PRD §3.3)

Exact Prize equations + auxiliary scalar θ slaved to velocity:

```
∂_t θ + u·∇θ = ν·Δθ + ν·|∇u|²     θ(x,0) = 0
μ_eff = ν + ε·f(θ)   →   Prize equations recovered as ε → 0
```

Most direct answer to Tao's objection: works inside exact Prize equations, no limit required.

---

## Six-Phase Proof Architecture (v0.4)

| Phase | Goal | Status |
|---|---|---|
| **Phase 0** (Months 1–3) | Local existence for Route 1; continuation criterion: singularity ⟹ A(T)→0 or μ(T)→0 | Planned |
| **Phase 1** (Months 1–8) | T scalar master PDE globally well-posed; T ∈ [T_min, T_max]; A(T), μ(T) uniformly bounded | Planned |
| **Phase 2** (Months 4–18) | Four regularisation mechanisms; Lemma 2.2 strengthened for Wang et al. (Conjecture 3.5) | Planned |
| **Phase 3** (Months 12–28) | Global strong solutions for Route 1; continuation criterion never triggered | Planned |
| **Phase 4** (Months 24–36) | μ(T)→ν limit; LPS margin ε survives; Prodi-Serrin ⟹ Prize | Planned |
| **Phase 5** (Months 32–48) | Assemble, 3D adversarial experiments, peer review, submit | Planned |

The **computational program** (PRD §7) runs alongside all phases.
Every numerical result targets a specific mathematical claim with a measurable PASS/FAIL criterion.

---

## Directory Structure (v0.4)

```
tfirst/
├── CLAUDE.md                   ← You are here
├── PROGRESS.md                 ← READ THIS FIRST EVERY SESSION
├── TASKS.md                    ← Milestone task breakdown
├── tfirst_prd_v4.docx          ← Primary specification (authoritative — v0.4)
├── tfirst_prd_v3.docx          ← Previous version (archived)
├── tfirst_program.tex          ← LaTeX research document
│
├── verified/                   ← Locked baseline experiments (PRD §4)
│   ├── tfirst_verify.py        ← Ideal gas, 128² spectral, 5 experiments — ALL PASS
│   └── sco2_verify.py          ← SC-CO2 + CoolProp, 7 claims — ALL PASS
│
├── layer1/                     ← Property Engine (M1)
│   └── tfirst_props.py         ← scalar_T_props, route1_coeffs, route2_theta,
│                                   A_field, second_law_check, self_similar_props
│
├── layer2/                     ← 2D Spectral Solvers (M0–M4)
│   ├── self_similar_IC.py      ← Wang et al. profiles as ICs (M0 infrastructure) ✅
│   ├── T_solver_2D.py          ← Route 1 2D solver: incomp. NS + scalar T (M0) ✅
│   ├── LPS_monitor.py          ← Track LPS norms and suppression margin ✅
│   ├── stretching_2D.py        ← Dissipation/stretching ratio diagnostics ✅
│   ├── route2_2D.py            ← Route 2: exact Prize NS + auxiliary θ (M2)
│   ├── mu_limit_2D.py          ← μ(T)→ν sweep; ε(μ) vs μ (M4)
│   └── lambda_sweep_2D.py      ← Lambda sweep driver + results logger (M3)
│
├── layer3/                     ← 3D Spectral NS (M5–M9)
│   ├── route1_3D.py            ← Route 1 3D: incomp. NS + scalar T (M5, 128³)
│   ├── route2_3D.py            ← Route 2 3D: Prize NS + θ (M7, 128³×8)
│   ├── wang_profiles_3D.py     ← CCF/Boussinesq unstable profiles as 3D IC (M6)
│   ├── mu_sweep_3D.py          ← μ(T)→ν in 3D; ε(μ) (M8, 64³×15)
│   ├── blowup_search_3D.py     ← Adversarial IC at 256³ (M9)
│   └── continuation_monitor.py ← A(T)_min and μ(T)_min continuously tracked
│
├── layer4/                     ← Diagnostic Tools (M10)
│   ├── scalar_T_check.py       ← Confirm T has no directional dependence
│   ├── lps_classifier.py       ← u field → (p,q) Prodi-Serrin pair
│   ├── continuation_criterion.py ← A_min and μ_min tracking; Phase 0 criterion
│   ├── gronwall_fitter.py      ← E1(t) → damping constant c
│   ├── viscosity_scaling.py    ← μ(T) along Wang et al. profile vs (1-t)^{-(1+λ)}
│   ├── lambda_threshold.py     ← Lambda sweep → λ_c or confirm none
│   ├── mu_limit_tracker.py     ← ε(μ) as μ(T) → ν
│   └── blowup_detector.py      ← Alarm if norm exceeds threshold
│
├── results/                    ← All experiment outputs
│   ├── results.db              ← SQLite results database (auto-generated)
│   └── figures/                ← Generated plots
│
└── SESSION-LOG/                ← Session logs (one per session)
```

---

## Experiment ID Convention

Every experiment must have a unique ID: `EXP-LAYER-ROUTE-NNN`

Examples:
- `EXP-L1-R1-001` — Layer 1, Route 1, first experiment
- `EXP-L2-R1-CCF-003` — Layer 2, Route 1, CCF profile, third experiment
- `EXP-L2-R2-001` — Layer 2, Route 2 (theta system), first experiment
- `EXP-L3-R1-ADV-001` — Layer 3, Route 1, adversarial IC

Every result is logged as: `claim_id | result: PASS/FAIL/PARTIAL | key_metric | figure_path`

---

## Key Mathematical Objects — Naming Convention

| Symbol | Python variable | Description |
|---|---|---|
| T(x,t) | `T_field` | **SCALAR** temperature field (master variable) — max principle applies |
| A(T) | `A_field` | Thermal diffusivity k/(ρ·cv) — MUST be > 0 — master coefficient |
| μ(T) | `mu` | Dynamic viscosity — always > 0 |
| k(T) | `k_th` | Thermal conductivity |
| ρ | `rho` | Density — **constant** in Route 1 (ρ = const, div u = 0 exact) |
| cv | `cv` | Specific heat at constant volume |
| θ | `theta` | Auxiliary scalar (Route 2 only): ∂_t θ + u·∇θ = ν·Δθ + ν·|∇u|² |
| μ_eff | `mu_eff` | Route 2 effective viscosity: ν + ε·f(θ) |
| ε(μ) | `lps_margin` | LPS margin — must stay > 0 as μ(T) → ν |
| λ | `lambda_val` | Wang et al. self-similar scaling parameter |
| λ_c | `lambda_c` | Hypothetical T-first failure threshold (Conjecture 3.5: λ_c = 0) |
| E0, E1, E2 | `E0, E1, E2` | Energy hierarchy (L², H¹, H²) |
| S_ij | `S_ij` | Strain rate tensor — **downstream** of u (NOT a primary variable) |
| τ_ij | `tau_ij` | Viscous stress 2μ(T)S_ij — **downstream** |

**CRITICAL:** T is a scalar. S_ij, τ_ij, q are downstream derived tensors/vectors.
They do NOT appear as primary variables in the T master equation.

---

## Architecture Rules — Absolute

### Property Engine (Layer 1)
- ALL thermophysical property lookups go through `layer1/tfirst_props.py`
- NEVER inline property formulas in solver files — call the property engine
- `A_field()` must return strictly positive values — assert this on every call in debug mode
- `second_law_check()` must pass for every fluid before any solver run
- `route1_coeffs(T_field)` returns μ_min, μ_max, A_min, A_max — uniform bounds for LPS margin
- `route2_theta(u_field, nu)` computes the auxiliary scalar θ for Route 2

### Solvers (Layers 2–3)
- **Route 1**: ρ = const and div u = 0 must be enforced **exactly** (not approximately)
- **Route 2**: ν = const, θ evolves, μ_eff = ν + ε·f(θ)
- Spectral methods only (FFT-based) — finite difference is not acceptable
- Dealiased via 2/3 rule for all runs
- RK4 time integration throughout; adaptive CFL
- NEVER hardcode fluid properties — always call Layer 1

### Results and Reproducibility
- Every run saves: parameters JSON + raw output + figures + PASS/FAIL verdict
- Random seeds always set and logged
- Grid resolution always logged in results
- Results appended to `results/results.db`
- A result without a logged claim_id is not a result

### Wang et al. Interface (PRD §3.7, §7.3)
- Self-similar profiles parameterised by λ
- CCF stable: λ = 1.1808 · CCF 1st unstable: λ = 0.6057 · CCF 2nd unstable: λ = 0.4703
- Boussinesq stable: λ = 1.9206 · Boussinesq 1st: λ = 1.3991 · Boussinesq 3rd: λ = 1.1843
- Lambda sweep values (PRD §7.3): {1.2, 0.9, 0.6057, 0.4703, 0.3, 0.1}
- `self_similar_IC.py` is the single source for all Wang et al. profiles as ICs
- Conjecture 3.5: μ(T) scaling defeats self-similar blow-up for **all** λ > 0 — no λ_c threshold

---

## What NOT To Do

- Do NOT hardcode thermophysical properties — always call `layer1/tfirst_props.py`
- Do NOT use ρ as variable in Layer 2 Route 1 solver — ρ is constant
- Do NOT use finite-difference methods in the spectral solvers
- Do NOT run Layer 3 experiments without Layer 1 and Layer 2 validated first
- Do NOT edit existing verified experiments in `verified/` — those results are locked
- Do NOT claim PASS without a logged key_metric in results.db
- Do NOT log a result without its experiment ID
- Do NOT guess Wang et al. λ values — use exact values from PRD §4.3
- Do NOT refer to "Conjecture 3.4" — it is now Conjecture 3.5 in PRD v0.4

---

## Test Commands

```bash
# Layer 1: Property engine unit tests
python -m pytest layer1/ -v

# Layer 2: Route 1 2D solver smoke test
python -m pytest layer2/test_M3.py -v

# Layer 2: Lambda sweep (smoke test, 32² grid, 10 steps)
python layer2/T_solver_2D.py --smoke-test --verbose

# Layer 2: Full test suite
python -m pytest layer2/ -v

# Full suite (all layers)
python -m pytest layer1/ layer2/ -q

# Verified experiments (should always pass — locked)
python verified/tfirst_verify.py
python verified/sco2_verify.py
```

---

## Priority Open Questions (PRD §10 — v0.4)

1. **[PRIORITY]** Phase 3 large-data regularity: Can global strong solutions be proved for Route 1 (incomp. NS + μ(T)) for large data? This is open in the literature. Phase 1-2 provide the tools.
2. **[PRIORITY]** Phase 4 μ limit: Does LPS margin ε(μ) remain > 0 as μ(T) → ν? If yes, Prize follows immediately.
3. **[PRIORITY]** Route 2 ε limit: Does the θ-system converge to exact Prize equations with LPS margin inherited as ε → 0?
4. Lambda threshold: Does Conjecture 3.5 hold for all λ > 0, or does λ_c exist? Lambda sweep answers this.
5. Boundary-free λ minimum: Do Wang et al. boundary-free profiles have λ ≥ λ_min ≥ 1?
6. Is A(T) > 0 provable from thermodynamic stability alone, independent of equation of state?
7. Does scalar T degeneration at T_λ = 2.17K constitute a rigorous prediction of the He-4 λ-transition?

---

## Commit Rules

- One commit per completed milestone
- Commit message format: `[M{N}] verb: description`
  - e.g. `[M1] feat: layer 1 property engine with route1_coeffs and route2_theta`
  - e.g. `[M0] feat: route1_2D solver — incompressible NS + scalar T, 64/64 tests pass`
  - e.g. `[M3] verify: lambda sweep 2D — Conjecture 3.5 numerically confirmed`
- ALWAYS update PROGRESS.md before committing
- NEVER commit with failing tests
- NEVER commit without a results.db entry for the completed experiment

---

## Session Protocol

1. Read `PROGRESS.md` — understand exactly where the program stands
2. Read the relevant PRD section from `tfirst_prd_v4.docx` for the current milestone
3. Do the work
4. Run the appropriate tests — they must pass
5. Log the result to results.db with claim_id + key_metric + verdict
6. Update `PROGRESS.md`
7. Write SESSION-LOG entry
8. Commit

If unsure about any mathematical claim, go to `tfirst_prd_v4.docx`. Do not invent claims.
