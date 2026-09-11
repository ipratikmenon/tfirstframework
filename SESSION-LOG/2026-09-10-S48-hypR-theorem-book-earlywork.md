# Session S48 — 2026-09-10/11 — WP2 settles Hypothesis (R); the book begins

Two parallel Opus subagents under `HANDOVER-S46-OPUS.md`. Both merged only
after independent re-derivation of their load-bearing steps.

## Headline

**Hypothesis (R) is a theorem.** WP2 returns verdict (i), and it is the first
unambiguously positive structural result since S38. **Wall 2 is withdrawn.**

## WP2 — Hypothesis (R)

`sec:hypR-s48`, proof document 10,462 → 11,309 lines, 152 pages, compiles
clean. `thm:s48-R` proves (R) in a *stronger* form than the document stated:

    ‖∇^k ∂_s^l u_λ(·,s)‖_∞ ≤ C_{k,l}(c_I) · |s|^{-(k+1)/2 - l}

for all k, l ≥ 0, uniformly in λ ∈ (0,1], with **no** (1+R) factor and **no**
S₂ dependence. Taking k=1, l=0 is exactly (R).

**Why the named trap never bit.** The handover warned that ε-regularity needs
*smallness* of scaled local energy while Type-I supplies only *boundedness*.
That trap is real, and this route never touches it. The bound actually
available on the rescaled family is an L^∞_t L^∞_x **velocity** bound, and the
pair (s′,s) = (∞,∞) gives 2/s′ + 3/s = 0 < 1 — the **strict**, subcritical
range of Serrin's criterion, which Serrin proved in 1962 and where boundedness
alone suffices. No smallness, no pressure hypothesis, no ε-regularity.

**Mechanism.** A restart argument on KNSS 2009 Prop. 4.1, transferred to
Λ_λ = λ^{-1}𝕋³ by a periodised Oseen kernel whose only estimate
(∫|K^{(L)}| ≤ Ct^{-1/2}) is uniform in the period; parasitic solutions are
excluded outright because a harmonic function on a closed manifold is
constant.

**Downstream.** `prop:transfer-audit(d)` first clause: Proved-mod-(R) →
**Proved**. Wall 2: **open → withdrawn**. The profile route does not relocate
the CZ logarithm; on that branch it never passes through the
Calderón–Zygmund operator at all.

**Correction to the document's own reasoning.** (d)'s stated mechanism — "U
has, by construction, no structure below the scale |s|^{1/2}" — is an unproved
heuristic and is *wrong*. The conclusion stands for the different reason above
(`rem:s48-mechanism`). The correction strengthens the result.

### Residue — four items, none of them (R)

1. **Saturation, now the sharpest gap on this branch.** Part (d) has *two*
   clauses. The second (‖∇U‖_∞ ≤ C·M_U) carries a **pre-existing** proviso,
   "whenever the Type-I upper bound is saturated". Part (e) uses that form
   (r* = M_U^{-1/2}), so (e) moves Proved-mod-(R) → **Proved-mod-saturation**,
   not Proved. Verified as genuinely pre-existing, not introduced by the agent.
2. **The ω→u bridge is permanently closed.** ‖v‖_∞ ≲ ‖ω‖_∞^{3/5}‖v‖_{L²}^{2/5}
   is log-free (the 2^{-j} sum converges, which is exactly where the CZ
   logarithm is not paid), but 3/5 is *forced by scaling* and 3/5 > 1/2. From
   vorticity Type-I plus bounded energy the velocity bound is unreachable — a
   scaling identity, not a missing technique. No future session should retry.
3. **(R) as *used* exceeded (R) as *stated*.** `prop:profile-extraction`
   invoked "the local pressure bound supplied by (R)"; (R) supplies none.
   Genuine internal gap, repaired without any pressure bound.
4. **Nothing else moves.** (N), (N′), (A-up), `target:disorder-depletion`,
   Wall 1, Wall 3, Type-II all untouched. `rem:s39-distance` stands verbatim:
   the ceiling is still "no Type-I singularity", strictly below the Prize.

### Recorded, deliberately not acted on

`obs:s48-wall1`: the same estimate applies to u itself on 𝕋³, giving
‖∇u(·,t)‖_∞ ≤ C(c_I)(T−t)^{-1}, logarithm-free, under the document's standing
hypotheses. This does **not** remove Wall 1 — §33–§37 need ∇u measured against
M(t), not (T−t)^{-1}, and the two agree only on the saturated set. The agent
recorded it, used it nowhere, and did not revise Wall 1. **Strategic
implication: saturation, not the logarithm, may be what both branches
actually need.** Test that before attacking either wall separately.

### Document errors surfaced

- `\cite{Seregin2012}` **miscited twice** (lines 6810, 6848). The paper is an
  L³-endpoint blow-up result; it contains neither the "standard caveats" nor
  any Type-I-improved ε-regularity theorem. (N) is unaffected — still assumed.
- `\cite{LadyzhenskayaSeregin1999}` loosely attributed at 6798 as a
  "sharpening by local-pressure theory". It is an ε-regularity paper for
  *suitable* weak solutions requiring p ∈ L^{3/2}_loc and a smallness
  hypothesis. Not the tool; no local-pressure theory is needed at all.
- (R)'s own caveat (ii) is **void**: it worries that u_λ has no uniform global
  energy bound. True premise, invalid inference — the L^∞-based local theory
  uses no energy anywhere.

### My verification before merging

- Re-derived the restart arithmetic from scratch: with A = c_I S^{-1/2} and
  τ ≍ S, KNSS k=1 gives C·A·τ^{-1/2} ≍ C(c_I)S^{-1}. Matches `thm:s48-R`.
- Re-derived the forced 3/5 exponent independently by dimensional analysis
  (a − 3b/2 = 0 from length scaling, a + b = 1 from amplitude ⟹ a = 3/5).
- Confirmed the line-5565 standing assumption verbatim — it assumes **both**
  M(t) ≤ C_I/(T−t) **and** sup|u| ≤ C_I(T−t)^{-1/2}.
- Confirmed the saturation proviso and all three citation defects at their
  exact lines; confirmed the diff is purely additive outside the Wall 2
  subsection and the summary-table row.
- Ran `check_proof.py` (exit 0) and two `pdflatex` passes (exit 0, 152 pp).

**One agent imprecision corrected.** It reported "0 LaTeX warnings". True for
`LaTeX Warning` proper, but the build carried 113 hyperref bookmark warnings
against a baseline of 109 — four new, from one unguarded math subsection
heading. Fixed with `\texorpdfstring`; now back to baseline.

**A steer of mine was refuted, and rightly.** My WP2 pre-read asserted that the
programme's Type-I hypothesis is on the *vorticity*, making an ω→u bridge the
real work. The document assumes the velocity bound outright (line 5565), so the
premise was false; and the bridge is impossible anyway (residue 2). I had read
the summaries rather than the document — the exact failure mode the house rules
exist to catch. Recorded because the erroneous steer is in the scratchpad and
in my own messages.

## The book — Parts I–III and the Landscape chapter

Owner directive: the Claim-A proof attempt is **not** to be removed or quietly
superseded. It is retained and presented as **Early Work**, read forward.
Master file restructured accordingly: `\part{The Early Work}` (chs 08–11) and
`\part{The Landscape}` (ch 17).

Thirteen chapters written, 4,920 lines, 108 pages, 0 errors, 0 LaTeX warnings,
no undefined references. Chapters 08/09/10 each carry the fixed four-part
template — what we assumed / what we considered / what it produced / what
survived. Ch 11 "Inheritance" is the hinge, with `tab:cz-life` tracing the
L^∞ CZ bound from a clause in Audit Defect D8 to the programme's central
barrier, rediscovered three times as an apparently new difficulty. Ch 17
covers OpenAI's forced construction against the four Clay alternatives.

### Repairs during verification

- **The master's `record` environment was malformed.** Splitting
  `\colorbox{c}{` across the begin- and end-code of `\newenvironment` leaves
  the begin-code brace-unbalanced; TeX mis-parses it at definition time and
  fails only on first *use* — which is exactly why it survived the
  empty-chapter build in the previous session. Rewritten via `lrbox`. Worth
  carrying forward: an environment that is never exercised is never tested.
- `tab:cz-life` amended with an S48 row. The chapter was drafted before WP2
  landed and stated Wall 2 as "relocated, not removed" — correct at drafting,
  stale within the hour. Both cautions kept in the new row.
- That row overflowed the float; table converted to `longtable`.

### Not done

Chapters 12–16, 18–22 and appendices A1–A4 are unwritten — the agent was cut
off by a session limit. The book master `\include`s them; LaTeX notes the
missing files and builds anyway.

## Status against the handover

WP0 ✅ · WP1 ✅ · **WP2 ✅ (verdict (i))** · WP3 (axisymmetric-with-swirl
analytical pivot) and WP4 (numerics redirect) remain. **But S48 has changed
the priority ordering**: `obs:s48-wall1` suggests saturation is the common
requirement of both walls, and part (e)'s saturation proviso is a statement
about one scalar function of one variable — a cheaper and better-defined
target than either wall, and one that did not exist as a named question
before this session.

## Files changed

- `proofs/claim_a_3d_proof_attempt.tex` (+847 lines; `sec:hypR-s48`; Wall 2
  reclassified; S39 ledger amended; `prop:profile-extraction` repaired)
- `proofs/claim_a_3d_proof_attempt.{pdf,log}` regenerated from current source
  (retained, not removed, per the owner's directive on the early work)
- `book/tfirst-book.tex` (restructured; `record` environment repaired)
- `book/chapters/` — 13 new chapters
- `PROGRESS.md` (KEY FINDING S48; two decision-log entries), this session log
