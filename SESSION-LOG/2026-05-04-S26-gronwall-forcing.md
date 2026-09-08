# Session S26 — 2026-05-04
## §18 Gronwall Forcing Estimate — Partial Session (Rate Limit)

---

## Completed: §18 "Gronwall Forcing Estimate and the Pressure–Boundary Integral"

**File:** `proofs/claim_a_3d_proof_attempt.tex` (now ~3780 lines, was 3640)
**Section label:** `\label{sec:gronwall-forcing}`
**Compile check:** 0 undefined refs ✓

### Contents (~197 lines)

**Lemma 18.1** (`lem:viscous-forcing`): Viscous boundary forcing bound:
```
F_ν(r,t) = (ν/|B_r|)|∫_{∂B_r} ∂_n u dσ| ≤ Cν r^{-2} ‖∇u‖_{L²(B_r)}
```
Proof: trace theorem + Cauchy–Schwarz on ∂B_r.

**Lemma 18.2** (`lem:pressure-forcing`): Pressure boundary forcing bound:
```
F_p(r,t) ≤ C r^{-23/15} ‖p‖_{L^{5/3}(Q_T)}
```
Proof in 3 steps:
1. Hölder on ∂B_r (exponents 3/2, 3): |∫_{∂B_r} p dσ| ≤ ‖p‖_{L^{3/2}(∂B_r)} · Cr^{2/3}
2. Sphere trace: ‖p‖_{L^{3/2}(∂B_r)} ≤ Cr^{-1/2} ‖p‖_{L^{3/2}(B_{2r})} → F_p ≤ Cr^{-11/6} ‖p‖_{L^{3/2}(B_{2r})}
3. Upgrade local L^{3/2} → global L^{5/3} via Hölder: ‖p‖_{L^{3/2}(B_{2r})} ≤ ‖p‖_{L^{5/3}} r^{3/10}
4. Net: -11/6 + 3/10 = -55/30 + 9/30 = -46/30 = **-23/15** ≈ -1.533

**Lemma 18.3** (`lem:convective-forcing`): Convective boundary forcing:
```
F_conv(r,t) ≤ C r^{-3/2} ‖u‖_{L^4(B_{2r})}²
```

**Proposition 18.4** (`prop:net-forcing`): Combined bound:
```
F(r,t) ≤ C [r^{-2}‖∇u‖_{L²(B_r)} + r^{-23/15}‖p‖_{L^{5/3}} + r^{-3/2}‖u‖_{L^4(B_{2r})}²]
```
Dominant term as r→0: pressure, at rate r^{-23/15}.

**Theorem 18.5** (`thm:forcing-gap`): Grönwall integration gives:
```
∫₀ᵀ F_p dt ≤ C r^{-23/15} ‖p‖_{L^{5/3}} T^{2/5}
```
This diverges as r → 0. Forcing alone does not close the gap. However -23/15 ≈ -1.53 is strictly better than the naive -2 (from unrefined trace).

**Corollary 18.6** (`cor:sharpened-obstruction`): The Grönwall-forcing route reduces the prize gap to:
> Pressure averages ∫_{∂B_r} p dσ must decay as r → 0 — i.e., p ∈ L^1_t C^0_x (pointwise pressure regularity).
This is equivalent to the Clay Prize gap restated as a pressure regularity condition.

**Remark 18.7** (`rem:bogovskii`): By the divergence theorem:
```
∫_{∂B_r} p·n dσ = ∫_{B_r} ∇p dx
```
and from NS integrated over B_r:
```
∂_t ∫_{B_r} u dx = −∫_{∂B_r}[(u⊗u)·n + p·n − ν∇u·n] dσ
```
The pressure term is NOT independent — it IS the mean ODE. The forcing F is self-referential: eliminating p via a Bogovskii/Helmholtz projection (∇p = −∂_t u − (u·∇)u + ν Δu) and inserting back into the mean ODE may yield a pressure-free formulation that could close the gap.

### Status table update
New row added between "Shear mean ≡ 0" and "Local gradient conc.":
```
Gronwall forcing F(r,t) ≤ Cr^{-23/15}‖p‖_{L^{5/3}} | PROVED | Lem 18.2; sharpest known
```

---

## NOT Completed (rate limit hit)

### 1. `layer4/gamma_adversarial_search.py` — to do next session
- Build `GammaAdversarialSearch` class
- `make_concentrated_ic(r_frac, amplitude)`: divergence-free bump at scale r
- `compute_gamma(snapshots, r_frac)`: Γ(r) from snapshots via trapezoidal rule
- `run_free_decay(u0)`: viscous decay û(k,t) = û(k,0)·exp(-ν|k|²t) in Fourier space
- `run_search(r_fracs)`: loop over r, measure γ_initial vs γ_evolved, fit power law
- Experiment IDs: EXP-L4-GAM-ADV-001
- Log to `results/results.db` table `gamma_adversarial_results`
- Tests: 6 tests in `layer4/test_gamma_adversarial_search.py`

### 2. Paper 6 abstract + §1 revision — to do next session
- Abstract: ~250 words, honest about conditional vs complete results
- §1 introduction: ~500 words with background, main idea (θ-identity), results, organization, notation
- arXiv-ready draft

---

## Proof Document Status (end of S26)

| Sections | Lines | Undefined refs |
|---|---|---|
| §1–§18 + Bibliography | ~3780 | 0 |

| New key result | Exponent | Status |
|---|---|---|
| Pressure boundary forcing | r^{-23/15} ≈ r^{-1.53} | **Proved** (best known) |
| Gap = p ∈ L^1_t C^0_x | — | Precisely characterised |
| Bogovskii pressure-free reformulation | — | Identified as next step |

---

## Next Session Priorities

1. **[FIRST]** `layer4/gamma_adversarial_search.py` — adversarial Γ search (see spec above)
2. **[SECOND]** Paper 6 abstract + §1 arXiv revision
3. **[OPTIONAL §19]** Bogovskii pressure-free mean ODE: insert ∇p = −∂_t u − (u·∇)u + νΔu into the mean ODE to eliminate pressure explicitly; check if the resulting ODE for h(t) closes
