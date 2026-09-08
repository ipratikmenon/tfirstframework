# Session S34 — 2026-07-19 — Magnitude Caccioppoli: A Priori σ* Bound + K-Correlation Numerics

## Purpose

Attack the K-absorption target (the sharpened remainder of Bridge Lemma Obstacle (i))
via the magnitude equation's good-sign sink, per the S33 plan. Fable derivation;
Sonnet execution.

## Corrigendum to §S33 (self-audit)

A dimensional re-audit of my own §S33 section — before building on it — found and
fixed four defects:
- Lemma BS-split / Prop first-rung inhomogeneous terms corrected:
  ν⁻¹𝓛M²(r*)³ → ν⁻¹𝓛M²(r*)⁵ and C_I²(r*)⁻¹ → C_I²·r* (every term of the first-rung
  inequality scales as ℓ¹ in parabolic units). The corrected terms VANISH as r*→0 —
  the fix is favorable, not damaging.
- K-absorption conjecture allowance likewise corrected to C_K(1+σ*)·r*.
- Routine item (iii) RECLASSIFIED: the νwê term yields the critical quadratic ν∬w²
  — the standard critical nonlinearity requiring σ* ≤ ε₀ smallness at the second
  rung (GN interpolation), not bookkeeping. This is where the Bridge Lemma's
  smallness hypothesis enters the scheme, as it should.

## New analytical results — §magnitude-s34 (proof doc, lines 5766–5922)

**Theorem (a priori coherence bound under Type I) — PROVED.**
Testing the magnitude equation D_t|ω| = α|ω| + ν(Δ|ω| − |ω|w) with mη² puts the
good-sign sink ν∬m²wη² on the LHS. Every RHS term is O(M^{1/2}𝓛): the key device is
the pointwise orthogonal-split inequality m²|∇ê|² = |∇ω|² − |∇m|² ≤ |∇ω|², which
sends every transition-annulus term to the local enstrophy-gradient bound with NO
absorption argument needed. Result:

    ⟨σ*⟩ = (r*)⁻³ ∬_{Q(r*/2)} w 𝟙_{|ω|≥M/2} ≤ C(C_I)·𝓛/ν.

The coherence deficit near a Type-I singularity is a priori bounded (order log).
The gap between "automatic" and "needed" for the Bridge Lemma is now the explicit
factor C𝓛/(νε₀), turning the program's smallness requirement into a decay question
for the σ*-dynamics. First unconditional-within-Type-I quantitative output of the
pincer program.

**Proposition (crude K-absorption fails by exactly one log) — PROVED (sharp).**
Hölder + enstrophy bound give νK_w ≤ C𝓛‖w‖_∞(r*)³ — top-order re-entry times 𝓛.

**Conjecture (refined K-absorption via equipartition) — OPEN, MEASURABLE.**
c_K = ⟨|∇m|²w⟩_E/(⟨|∇m|²⟩_E⟨w⟩_E) ≤ C/𝓛 near Type-I blowup. Mechanism: the
orthogonal split makes |∇m|² and m²w competitors for the same |∇ω|² budget; only
strong POSITIVE correlation (c_K ≫ 1) could defeat absorption.

## Numerical results — c_K measured for the first time

New module layer4/k_correlation.py (360 lines) + 19 tests; alignment_bridge.py wired
to record c_K/f_ang per snapshot; escalated-amplitude adversarial IC support added.
Tests: layer4 260/260; whole program 969/969, zero regressions.

| exp_id | c_K med/max/final | f_ang med | verdict |
|---|---|---|---|
| EXP-L4-KCORR-TG-001 (ν=1e-3) | 1.063 / 5.30 / 1.17 | 1.00 | PASS |
| EXP-L4-KCORR-SH-001 (ν=1e-3) | 1.935 / 2.90 / 1.70 | 0.823 | PASS |
| EXP-L4-KCORR-ADV-001 (ν=1e-3) | 1.156 / 1.35 / 1.15 | 0.670 | PASS |
| EXP-L4-BRIDGE-ADV-002 (amp 3×, ν=5e-4) | 1.092 / 1.31 / 1.09 | 0.851 | PASS |

Escalated run bridge stats: R max/med/final = 1.521/0.995/1.025, σ* ∈ [1.28, 2.77]
— saturation band again, never violated.

**Headline:** c_K ≈ 1 (decorrelated) in every flow, drifting TOWARD 1 late-time,
never diverging; escalation slightly LOWERED c_K. No positive-correlation mechanism
that could defeat K-absorption is visible in any tested flow. f_ang (0.67–1.0) shows
the angular term dominates the enstrophy-gradient budget on the high set — the
magnitude sink −ν|ω||∇ê|² is strong there, as the absorption argument wants.

**Honest reading vs the conjecture:** measured c_K ≈ 1 confirms decorrelation
(absorbability with O(1) constant) but cannot yet distinguish c_K = O(1) from the
stronger conjectured c_K ≤ C/𝓛, since 𝓛 = O(1–10) at these Reynolds numbers.
Distinguishing requires a Reynolds sweep (c_K vs 𝓛 trend). Shear's c_K (~1.9) should
be read qualitatively (worst Gibbs split_residual, up to 6.3); TG/adv signals robust.

Session hiccup: the numerics agent twice stopped after detaching runs in the
background (runs died with it); resumed via message with synchronous-execution
instructions both times. Note for future prompts: forbid detached runs explicitly.

## Actions

- proofs/claim_a_3d_proof_attempt.tex: 4 corrigendum edits + §magnitude-s34 inserted
  (481/481 balanced; all refs verified; no label collisions).
- layer4/k_correlation.py + test_k_correlation.py new; alignment_bridge.py extended.
- results.db: 4 kcorr experiment rows + 164 timeseries rows + 1 bridge row.
- PROGRESS.md updated for S34 (Sonnet).
- Test total now 969/969 (recorded-count reconciliation done: 950 → 969 with S34).
- No commit — not a git repository.

## Next session (S35) candidates

1. The second rung: GN interpolation + σ* ≤ ε₀ handling of the critical quadratic
   ν∬w² (the reclassified item (iii)) — assemble first + second rungs into a
   conditional sup-bound modulo K.
2. Reynolds sweep for c_K (ν ∈ {1e-3, 5e-4, 2.5e-4} at N=64/128): does c_K trend
   with 1/𝓛? Distinguishes the conjecture's O(1/log) form from mere O(1).
3. Attempt the equipartition-based proof of c_K ≤ C: can the magnitude Caccioppoli
   itself bound the correlation (its sink penalizes exactly the c_K-large
   configurations)?
