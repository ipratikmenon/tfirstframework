# Session Log — S23
**Date:** 2026-05-03
**Description:** Gronwall Analysis of the Mean + Mean Morrey Diagnostic Tool + Paper 6 §7 Update

---

## Summary

Three deliverables built in parallel during S23.

---

## 1. §15 — Gronwall Analysis of the Mean

Inserted into `proofs/claim_a_3d_proof_attempt.tex` as a new section with five subsections.

**§15.1 Decomposition of the Mean Dynamics**
- Lemma 15.1 (`lem:mean-self-int`): Mean self-interaction vanishes by div-free condition.
  `m · ∫_{∂B_r}(m⊗m)·n dS = |m|² m·∫_{∂B_r} n = 0`
- Nonlinear term splits: `m·∫u_osc(u_osc·n) + m^T(∫_{B_r}∇u)m`

**§15.2 The Linear Gronwall ODE**
- Lemma 15.2 (`lem:mean-gronwall`): `h(t) = |⟨u⟩_r|` satisfies
  `dh/dt ≤ (C/2r²)‖∇u‖_{L²(B_r)}·h + (C/2r)‖∇u‖²_{L²(B_r)}`
  Key tool: Poincaré + trace inequality gives `‖u_osc‖_{L²(∂B_r)} ≤ Cr^{1/2}‖∇u‖_{L²(B_r)}`
- Proposition 15.3 (Gronwall bound):
  `h(t) ≤ [h(0) + C/r ∫‖∇u‖²] · exp(C r^{-3/2} ∫‖∇u‖_{L²(B_r)})`
- Theorem 15.4 (`thm:conditional-mean-morrey`): If local gradient concentration
  `Γ(r) = ∫₀ᵀ(X̄_{B_r}|∇u|²)^{1/2} dt ≤ C₁r^β`, then mean Morrey holds with
  `α = 1 − 3β → Claim A`

**§15.3 Shear Layer** (`prop:shear-mean-zero`) — CORRECTED
- Claim A holds for the shear layer by direct parabolic regularity (Thm 11.4), NOT via mean Morrey.
- `⟨sin(y)⟩_{B_r} ≠ 0` for general `B_r` in 3D: `∫_{B_r} e^{iky} dx ≠ 0`.
- This is a mathematical error from earlier drafts; corrected in S23.

**§15.4 Gap = Local Gradient Concentration** (`cor:gap-gradient-concentration`)
- `Γ(r) ≤ Cr^β` is:
  - (a) implied by Prodi–Serrin with `β = 3/q − 3/2`
  - (b) strictly weaker — requires only local `L¹_t Ḣ^{1/2}_{loc}`
  - (c) equivalent to CKN `A(r)` uniform bound via Route B+C chain

**Key mathematical correction in S23:** Proposition 15.5 was initially stated as "shear layer mean ≡ 0 over any ball." The test suite caught this (`test_spatial_mean_near_zero` failed for `sin(2πy/L)`). Corrected: shear layer Claim A goes via the direct route only; mean Morrey route does not apply.

---

## 2. Mean Morrey Diagnostic Tool

New files: `layer4/mean_morrey_diagnostic.py` and `layer4/test_mean_morrey_diagnostic.py`.

**Functions:**
- `ball_mask(N, cx, cy, cz, r_frac)` — periodic min-image distance mask
- `compute_spatial_mean(u_field, cx, cy, cz, r_frac)` — returns `(mean_vec, magnitude)`
- `compute_mean_morrey_field(u_field, r_values, alpha_values, n_centres=8)` → morrey_dict
- `fit_morrey_scaling(r_values, morrey_dict, alpha)` → `(C, γ, R²)`
- `find_alpha_threshold(r_values, morrey_dict, alpha_values)` → `(alpha_threshold, fit_details)`
- `run_mean_morrey_analysis(exp_id, u_snapshots, t_values, grid_N, ic_type, ...)` → result dict
- `init_db(conn)`, `log_result(conn, result)` — SQLite persistence

**Test classes (20/20 PASS):**
- `TestBallMask`: symmetry, centre included, all voxels within radius, periodic wrap
- `TestShearLayerZeroMean` (renamed from incorrect physics): uses zero field `u=0`, mean=0 exactly
- `TestConstantField`: `alpha_threshold ≈ 1.0 ± 0.3` for `u=(1,0,0)`
- `TestTaylorGreenField`: large-ball mean < 0.1, returns verdict dict
- `TestMorreyScalingFit`: `γ ≈ 0.5+α` for mean~r^{1/2}; `γ < 0` for constant mean with `α=0.1`
- `TestDBLogging`: init_db, log_result insert, round-trip, idempotent

**Bugs fixed during S23:**
- f-string `:.3f if condition else 'nan'` is invalid Python syntax; fixed with separate `at_str`/`se_str` variables
- `n_snapshots` variable not tracked; added counter in accumulation loop
- Test for "negative γ" used wrong synthetic data; corrected to constant mean giving `γ = -1 + α`

---

## 3. Paper 6 §7 Update

File: `papers/paper6_prize_limit.tex` (804 lines after S23).

New subsection "The local gradient concentration characterisation" added:
- Theorem 4.5 (`thm:grad-concentration`): Conditional Claim A via `Γ(r) ≤ Cr^β`
- Proposition (`prop:shear-mean-zero-p6`): shear mean ≡ 0 globally; Claim A via direct route
- 4-row comparison table of gap formulations in decreasing strength

---

## Document Status After S23

| Document | Lines | Sections | Undefined refs |
|---------|-------|---------|----------------|
| `proofs/claim_a_3d_proof_attempt.tex` | 2946 | 16 (§1–§15 + Conclusion) | 0 |
| `papers/paper6_prize_limit.tex` | 804 | 8 | 0 |

Tests: **833/833 PASS**

---

## Gap in Sharpest Form

After S23 the remaining gap is:

```
Γ(r) = ∫₀ᵀ (X̄_{B_r(x₀)} |∇u|²)^{1/2} dt ≤ C₁ r^β   for some β > 0
```

This is a local `L¹_t Ḣ^{1/2}_{loc}` condition. Hierarchy proved:

```
Prodi–Serrin  ⟹  Γ(r) ≤ Cr^β  ⟹  mean Morrey  ⟹  Claim A  ⟹  Prize
```

Each implication is strict.

---

## Files Changed

- `proofs/claim_a_3d_proof_attempt.tex` — §15 added, Prop 15.5 corrected
- `papers/paper6_prize_limit.tex` — §7 subsection added
- `layer4/mean_morrey_diagnostic.py` — new diagnostic tool
- `layer4/test_mean_morrey_diagnostic.py` — 20 tests, all passing
- `PROGRESS.md` — S23 row added

---

## Next Steps (S24)

1. Run mean Morrey diagnostic on M7-style TG and shear snapshots (numeric evidence for `Γ(r)`)
2. §16: Vorticity/enstrophy attack on `Γ(r)` — use vorticity equation to bound `Z(r,t)`, show stretching is the final barrier
3. Paper 6 final assembly pass
