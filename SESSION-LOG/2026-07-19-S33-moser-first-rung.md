# Session S33 — 2026-07-18/19 — Level-Set Moser First Rung + Adversarial Bridge Diagnostics

## Purpose

Attack Bridge Lemma Obstacle (i) — the singular drift 2ν∇log|ω| in the direction
equation — via a level-set Moser formulation (Fable derivation), and escalate the
numerical bridge diagnostic to N=64³ including the adversarial anti-parallel-tube IC
(Sonnet execution). Session spanned a session-limit interruption (both agents killed
mid-task on Jul 18 evening; proof-doc insertion had already completed; numerics agent
resumed from transcript on Jul 19 and finished).

## Analytical result — new proof-doc section §moser-s33

"Level-Set Moser Scheme for the Bridge Lemma: The First Rung" inserted after
§pincer-s32 (environments balanced 470/470). Content (all Fable-derived, honestly
labeled):

**Proved:**
- Non-circularity remark: the pincer is a contradiction argument at a putative first
  singularity, so all estimates run on the smooth interval [0,T) — S31's Defect 4
  (Leray–Hopf circularity) cannot recur on this route by design.
- Composite level-set cutoff η = ψ(x,t)·χ(|ω|/M) with derivative bounds
  (|∇χ| ≤ 8|∇|ω||/M; D_tχ bounded via the magnitude equation, with the Δ|ω| term
  never used pointwise).
- The three integrations by parts (Lemma): every second-derivative-of-magnitude term
  — (a) ∇²log|ω| from the drift, (b) Δ|ω| inside D_tχ, (c) χ′-cutoff terms — lands
  on first-order factors ≤ 4|∇|ω||/M and is absorbed by Young against ν|∇²ê|².
  Level-set boundary contributes nothing (χ vanishes below M/4).
- Localized stretching control under Type I: local enstrophy-gradient bound
  ν∬|∇ω|² ≲ M²ρ³·L (L = log factor), and Biot–Savart near/far splitting (CZ-L²
  local + far-field ≤ C·C_I(r*)⁻³).

**Proved modulo enumerated routine verifications (Prop first-rung):**
First-rung Caccioppoli for the angular energy — sup-slice + ν∬|∇²ê|² controlled by
C·L/(ρ′−ρ)²·∬w + Cν·K(ρ′) + Type-I critical terms. The four remaining routine items
are enumerated explicitly (transport term, a commutator constant, the νwê term, a
time-slice boundary term) — none involves ∇²log|ω| or nonlocality. Bookkeeping for
S34, not new estimates.

**Open (the sharpened target):**
K-absorption conjecture: absorb K(r*) = M⁻²∬|∇|ω||²(w + …) into the LHS with factor
θ_K < 1. This is STRICTLY WEAKER than Obstacle (i) (integrated coupling vs pointwise
drift bound). Two attack routes named; the natural S34 start is route (2): the
magnitude equation's own Caccioppoli carries −ν|ω||∇ê|² as a good-sign sink — the
radial equation dissipates exactly the coupling the angular equation needs absorbed.

## Numerical result — bridge diagnostic at N=64³ + adversarial IC

layer4/alignment_bridge.py extended (585 lines): adv_vortex_tubes_ic reusing M9's
make_vortex_omega + biot_savart_spectral unchanged; solver remains route2_3D (exact
Prize NS); IC passes through the solver's own Leray projector. 10 new tests → 46/46;
full layer4 241/241, zero regressions.

Production (N=64³, 400 steps, seed=42, ν=1e-3, ~2 min wall each):

| exp_id | R max/med/final | σ* range | M growth | verdict |
|---|---|---|---|---|
| EXP-L4-BRIDGE-TG-002 | 1.264 / 0.764 / 0.855 | [0.48, 177.9] | 10.5× | PASS |
| EXP-L4-BRIDGE-SH-002 | 0.763 / 0.587 / 0.546 | [4.0, 128.1] | 1.0× | PASS |
| EXP-L4-BRIDGE-ADV-001 | 1.364 / 1.059 / 0.981 | [3.18, 9.38] | 1.002× | PASS |

**Key scientific finding:** the adversarial anti-parallel-tube run NEARLY SATURATES
the bridge inequality (R ≈ 1 in a tight band [0.92, 1.36] for the whole run) without
ever violating it — pointwise and averaged coherence in near-exact lockstep at a
direction-reversal interface. This is precisely the configuration where a Bridge
Lemma counterexample would appear at this Re, and none does. R peaked (1.36) during
the closest tube-interaction window and then relaxed. Complementarity: ADV stresses
coherence structure (M flat), TG stresses stretching (M 10.5×, R ≤ 1.26) — both keep
R = O(1). No resolution-driven divergence N=32 → N=64.

Caveat: at ν=1e-3, amp=2 the tube collision does not produce strong stretching in
exact Prize NS; a higher-amp/higher-Re combined stress test is the natural escalation.

Logged: 3 rows layer4_bridge_experiments + 123 rows layer4_bridge_timeseries.

## Actions

- proofs/claim_a_3d_proof_attempt.tex: §moser-s33 inserted (Sonnet; completed before
  the session-limit interruption, verified after).
- layer4/alignment_bridge.py + test_alignment_bridge.py extended (Sonnet, resumed).
- results.db: 3 experiments + timeseries.
- PROGRESS.md updated for S33 (Sonnet).
- No commit — not a git repository.

## Next session (S34)

1. K-absorption via route (2): couple the magnitude-equation Caccioppoli (good-sign
   sink −ν|ω||∇ê|²) to the angular first rung; attempt the absorption with explicit
   θ_K.
2. Discharge the four enumerated routine verifications of Prop first-rung.
3. Optional numerics: higher-amp/higher-Re adversarial run (stretching + alignment
   combined stress test); N=128³ if warranted.
