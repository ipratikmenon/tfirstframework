# Session S35 — 2026-07-19 — Second Rung, K-Self-Absorption Identity, Reynolds Sweep

## Purpose

Advance H1 on both fronts per the S34 plan: (analytical) handle the critical
quadratic and attack K-absorption structurally; (numerical) test whether c_K trends
with Reynolds number. Also: answer the user's direct question "can we prove
existence with this?" honestly — recorded both in chat and in the proof doc itself.

## Honest answer to "does this prove existence?" — NO (recorded in doc)

Weak-solution global existence is classical (Leray 1934); the open question is
smoothness. Current dependency chain: Type-I + H1 + H2 ⟹ Type-I singularities need
persistent disorder; + σ*-decay (open) ⟹ no Type-I blowup; + Type-II exclusion
(wide open, outside current methods — every marginal estimate uses Type-I scaling)
⟹ Prize. A new "Honest distance to the Prize" remark and post-S35 ledger are now in
the proof doc (§secondrung-s35). Near-term realistic target: the conditional
disorder-barrier theorem — a real contribution, not the Prize.

## Analytical results — §secondrung-s35 (proof doc lines 5923–6083)

**Proposition (second rung) — PROVED given the first rung + ε₀-smallness.**
The critical quadratic ν∬w²η² (the reclassified item (iii)) is handled by the
pointwise split w² = f^{10/3}·w^{1/3} (f = |∇ê|η), the parabolic embedding
L∞L² ∩ L²H¹ ↪ L^{10/3}(Q) (LSU), and Young:
  ν∬w²η² ≤ δν(r*)³‖w‖_∞ + Cδ^{-1/2}νE_V^{5/2}(r*)^{-3/2},
dimensionally verified (every term ℓ¹). With σ* ≤ ε₀ the second term is
O((𝓛ε₀)^{5/2}r*) — superlinear in ε₀, harmless; the first absorbs into the Moser
sup-bound. The entire nonlinearity of the direction equation is subordinate to
σ*-smallness.

**Proposition (K-self-absorption identity) — DERIVED, PROVISIONAL (C1–C2 pending).**
Testing the magnitude equation with the weight mwη² (instead of mη²) puts BOTH
ν∬|∇m|²wη² (= the K-density) AND ν∬m²w²η² (the critical quadratic, favorably
weighted) on the LHS, controlled by the angular dissipation ν∬m²|∇²ê|²η² with an
absolute constant — no logarithm, no correlation hypothesis. The viscous cross-term
closes by Cauchy–Schwarz via |∇w| ≤ 2|∇²ê|w^{1/2}. Two couplings flagged for S36
verification: C1 (the D_tw substitution, including the ∇log m·∇w re-entry whose
constant must come out < 1) and C2 (the α-weighted marginal terms, 𝓛-degraded
constant chase). If C1–C2 verify: the c_K hypothesis becomes UNNECESSARY, and H1
reduces to σ* ≤ ε₀ + routine items (i),(ii),(iv).

## Numerical results — Reynolds sweep (6 runs, N=64, 400 steps, seed=42)

New layer4/kcorr_reynolds_sweep.py + 18 tests; layer4 278/278; program total
987/987, zero regressions. Session-limit interruption mid-sweep; agent resumed with
per-run incremental logging; all 6 runs completed synchronously.

| exp_id | ν | ic | c_K med | f_ang | tail_frac | res_warn | verdict |
|---|---|---|---|---|---|---|---|
| NU-TG-001 | 1e-3 | tg | 1.063 | 1.00 | 7.2e-2 | YES | PASS |
| NU-TG-002 | 5e-4 | tg | 1.021 | 1.00 | 3.3e-1 | YES | PASS |
| NU-TG-003 | 2.5e-4 | tg | 1.020 | 1.00 | 5.9e-1 | YES | PASS |
| NU-ADV-001 | 1e-3 | adv | 1.156 | 0.67 | 3.6e-6 | no | PASS |
| NU-ADV-002 | 5e-4 | adv | 1.120 | 0.84 | 1.7e-5 | no | PASS |
| NU-ADV-003 | 2.5e-4 | adv | 1.113 | 0.85 | 6.2e-5 | no | PASS |

Regression slopes of median c_K vs log(1/ν): TG −0.031, ADV −0.031 — flat with a
slight negative tilt.

**Headline:** c_K does NOT grow with Reynolds number across a 4× range; it sits
pinned at ~1.0–1.16 with a mild decreasing trend. This removes the one numerical
failure mode (c_K ≫ 1 positive correlation) that could defeat K-absorption, cleanly
supports O(1) absorbability, and is directionally consistent with (but too narrow a
range to confirm) the stronger 1/𝓛 form.

**Resolution caveat (prominent):** all three TG runs are under-resolved at N=64 for
these viscosities (spectral tail fraction 0.07–0.59) — the TG trend is indicative
only. The adversarial runs are well-resolved (tail ≤ 6e-5) and independently show
the same flat c_K ≈ 1.1; they are the trustworthy half. A resolved TG sweep needs
N=128+ (S36 optional item).

## Actions

- proofs/claim_a_3d_proof_attempt.tex: §secondrung-s35 inserted (6154 lines,
  490/490 balanced, all refs verified).
- layer4/kcorr_reynolds_sweep.py + test file new; 6 rows in layer4_kcorr_sweep.
- PROGRESS.md updated for S35 (Sonnet).
- Test total: 987/987.
- No commit — not a git repository.

## Next session (S36)

1. Verify couplings C1–C2 of the K-self-absorption identity (the decisive step:
   if they close, H1 = ε₀-smallness + routine).
2. Discharge routine items (i),(ii),(iv) of the first rung.
3. Optional: resolved TG sweep at N=128; extend ν range downward for the 1/𝓛 trend.
4. Then: H2 (localized CF at shrinking r*) becomes the program's front line.
