# Session S31 — 2026-07-18 — Independent Audit of §20–§22: Prize Claim Withdrawn

## Purpose

The user opened the session asking to explore a new possibility for proving global
existence for 3D NS. Before selecting a new route, an independent line-by-line audit
of the §20 vorticity scale recursion (the basis of the S29 "PRIZE PROVED for u₀ ∈ H¹"
claim) and the §21–22 σ-framework (S30) was performed. The audit found the Prize
claim invalid. This session records the findings, propagates the downgrade across all
program documents, and halts arXiv preparation.

## Audit findings — §20 (vorticity scale recursion)

**Defect 1 — Theorem 20.4 induction is arithmetically impossible for large data.**
The inductive absorption requires C_v·D^{5/3}·2^{10/9}·r^{10/9} ≤ (D/2)·r^{2/3},
i.e. D^{2/3} ≤ 2^{-19/9}/C_v (D small). The base case requires D ≥ ‖u₀‖²_{H¹}·T
(D large). The stated choice D = max(C₀, (2C_v·2^{10/9})³) satisfies neither.
This is the standard smallness obstruction for superlinear recursions; the argument
reduces to the CKN small-energy regime — the same limitation as §19. Notably, the
"Precise gap" remark in §21 (added S30) already acknowledged this internally, but the
acknowledgment was never propagated to §20's statements, PROGRESS.md, or paper 6.

**Defect 2 — Lemma 20.2 Step 3: invalid Hölder triple.**
|stretching| ≤ ‖∇u‖_{L³}‖ω‖²_{L²} uses exponents 1/3 + 1/2 + 1/2 = 4/3 ≠ 1.
The correct pairing ‖∇u‖_{L³}‖ω‖²_{L³} changes all downstream exponents
(Young then produces ‖ω‖⁶_{L²}, not ‖ω‖^{10/3}_{L²}); the 5/3 superlinearity
does not survive.

**Defect 3 — Lemma 20.2 Step 5: Jensen reversed.**
∫g^{5/3}dt ≤ r^{-2/3}(∫g dt)^{5/3} with g = ‖ω‖²_{L²(B_2r)} is the false direction
of Jensen's inequality (Jensen gives ∫g^{5/3} ≥ r^{-4/3}·(∫g)^{5/3}·r^{2/3}).
The conversion of the stretching term into C·Z(2r)^{5/3} is invalid.

**Defect 4 — Leray–Hopf circularity and norm mixing.**
The local enstrophy inequality (Lemma 20.1) tests the vorticity equation against
φ²ω — justified for smooth solutions only; ∇ω ∈ L²_loc is not known for Leray–Hopf
solutions (suitable weak solutions carry a local *energy* inequality, not enstrophy).
Steps 1–2 also use global T³ norms where local B_2r quantities are required.

## Audit findings — §21–22 (σ-framework, S30)

**Defect 5 — Theorem sigma-gronwall is vacuous as a regularity route.** Its hypothesis
∫₀ᵀ M dt < ∞ is precisely the Beale–Kato–Majda criterion; regularity already follows
classically under it.

**Defect 6 — Two-regime dichotomy needs large viscosity.** β > 0 requires ν > C₁C/c,
a condition on the fixed physical viscosity vs universal constants; not arrangeable
by scaling.

**Defect 7 — Regime II ODE conclusion is false.** dM/dt ≤ CM^{5/4}(log M)^{1/2} is
superlinear; the comparison solution blows up at finite t* = 4M(0)^{-1/4}/C. The
claim "stays finite on any fixed time interval" is an ODE error.

**Defect 8 — Lemma A-evol coercivity unproved.** The −cν∫|∇²ω|²φ² diffusion
contribution is asserted, not derived (ê singular on {ω=0}; sign-indefinite cross
terms). Also ‖∇u‖_{L∞(B_r*)} ≤ CM is false without a logarithmic factor (CZ operators
unbounded on L∞).

## What survives

- Route C bootstrap (Thm routeC) — unconditional. Intact.
- Claim A + regularity for shear ICs (Thm 11.4) and near-shear (Thm 12.3). Intact.
- λ_max = 0 for shear (Prop 11.6). Intact.
- Algebraic identity |∇ω|² = |∇|ω||² + |ω|²|∇ê|²; exact scale invariance of
  σ = A_loc/M^{3/2}. Intact — sound foundation for a quantitative alignment program.
- All numerical results: 873/873 tests, all experiment PASSes (evidence, not proof).

## New-route feasibility brainstorm (Fable, this session)

1. **Alignment pincer (σ + Constantin–Fefferman)** — most promising new direction.
   Gap to bridge: Theorem-E-style alignment is L²-averaged; CF needs pointwise/local
   direction coherence. A quantitative averaged→pointwise bridge is a genuine open
   research question, not a restatement of Claim A. Prerequisite: repair Defects 5–8.
2. **De Giorgi–Nash–Moser on θ** — obstruction identified up front: source
   ν|∇u|² ∈ L¹ only; scalar theory then gives θ ∈ L^q, q < 5/3, and no Hölder gain.
   Any improvement must exploit source = dissipation of the drift's own energy —
   which is Claim A restated. Worth one bounded session, documented as such.
3. **Probabilistic route** — a.s. global well-posedness for randomized H¹ data
   (Nahmod–Pavlović–Staffilani lineage). Publishable, but a different theorem from
   the Prize.
4. **§20 repair** — must re-derive the stretching exponent with valid Hölder/Jensen
   steps in a suitable-weak-solution framework; expected honest outcome is a
   conditional (small-data/CKN-type) statement.

## Actions taken (S31)

- Audit section "Independent Audit of §20–§22 (Session S31)" appended to
  proofs/claim_a_3d_proof_attempt.tex; status table updated (Sonnet agent).
- PROGRESS.md: Current State, Proof Pathway, priorities, session log downgraded
  (Sonnet agent).
- papers/paper6_prize_limit.tex: abstract downgraded; Theorem D → conjecture with
  withdrawal remark; vorticity + σ sections annotated (Sonnet agent).
- tfirst_program.tex: §20/Prize claims downgraded (Sonnet agent).
- arXiv submission HALTED.
- Test verification: layer1/layer2/layer4 → 513 passed; layer3 → 391 passed (exit 0).
  Total observed: 904, vs 873 recorded at S30 — the delta appears to be test additions
  never counted into PROGRESS.md; reconcile next session. No failures anywhere.
- No code changes; numerical layers untouched.

## Next session

- Decide repair-vs-reformulate for §20 (Defects 2–4).
- Begin the alignment-pincer program: audit-hardened restatement of the σ evolution
  (fixing Defects 5–8), then the averaged→pointwise coherence bridge.
- Optional bounded DGNM-θ session.

## Orchestration note

Thinking/planning/audit: Claude Fable 5 (this session's main context). Document
propagation: three parallel Sonnet subagents on disjoint files. Division per user
instruction.
