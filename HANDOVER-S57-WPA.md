# Handover — WP-A: price Grujić's hypothesis against this programme's own no-concentration family

Written 2026-09-13, immediately after S56 (the strategic audit). Owner's
direction: launch the audit's own recommended gated session.

## What this is, and what it is not

**This is a pricing exercise, not an adoption.** The question is whether
Grujić's hypothesis (arXiv:2607.08866v3) is a genuinely new, cheaper route
around Wall 1, or a restatement of this programme's own standing σ*-decay
gap in different vocabulary. **The most likely outcome, stated in advance
and not to be treated as a disappointing one, is the latter** — a fifth
independent confirmation (after S48, S50, S52, S53) that the obstruction is
real and does not move under relabeling. That outcome is a complete,
valuable result and ends the session cleanly. Do not manufacture a route
around Wall 1 that is not there.

## Background, established at S56 and to be re-verified, not re-trusted

- **The paper**: Zoran Grujić, "Logarithmic Depletion of Vortex Stretching
  and Singularity Evasion in the 3D Navier-Stokes Equations,"
  arXiv:2607.08866 (v1: 9 Jul 2026; v3: 9 Sep 2026). Claims: if the
  vorticity direction $\hat e=\omega/|\omega|$ lies locally in
  $\mathrm{bmo}_{1/|\log r|}$ (a logarithmically-weighted BMO space failing
  the Dini condition) near a critical-point singularity with
  $\omega\in L^\infty_tL^{3/2,\infty}_x$ concentration, then vortex
  stretching is depleted enough to avert blow-up (Theorem 7.4). Mechanism:
  recasts the stretching term as a singular-integral commutator, bounded via
  a localized Coifman–Rochberg–Weiss estimate. Companion paper on
  propagating the hypothesis: arXiv:2609.05720.
- **Verified independently at S56**: the paper's direction-equation terms
  match this document's `lem:direction-eq` (line 5349) term for term, and
  its critical regime ($\omega\sim|x|^{-2}$, i.e. $L^{3/2,\infty}$) matches
  this document's own criticality $M\cdot(r^*)^2=1$. **Re-verify this
  yourself against both the live document and the primary arXiv source —
  do not take S56's match on faith.**
- **The counterweight, also established at S56 and load-bearing**: unrefereed,
  single author. The same author's earlier claim in adjacent territory
  (vorticity-direction sparseness reduces the "scaling gap," arXiv:1704.05546,
  1911.00974) was scrutinised and found not to hold up by Albritton–Bradshaw
  (arXiv:2110.02187): *"the scaling gap is not reduced."* Read
  Albritton–Bradshaw's argument, not just its conclusion, before starting.
- **The cheap dictionary sketched at S56, not proved**: Poincaré's inequality
  gives $\fint_{B_{r^*}}|\hat e-\hat e_B|\lesssim r^*\bigl(\fint|\nabla\hat
  e|^2\bigr)^{1/2}\approx(\sigma^*)^{1/2}$ at the critical scale, so
  Grujić's hypothesis, translated to this document's own quantities at
  $r^*$, plausibly reads as $\sigma^*\lesssim1/(\log M)^2$ — a quantified
  form of the standing σ*-decay gap (open since S38).

## Deliverable 1 — read both papers in full, from the primary sources, and record proved vs. assumed

Fetch both arXiv papers (abstracts already confirmed at S56; you need the
full text, especially the proof of Theorem 7.4 and the propagation argument
in the companion paper). For each load-bearing step, classify: proved in
the paper, assumed as a hypothesis, or asserted without proof. **Do this
before forming any opinion about the paper's correctness** — the goal here
is an accurate map of what is actually established, not a verdict yet.

## Deliverable 2 — prove or refute the dictionary, precisely

Prove or refute the Poincaré translation above as a genuine lemma, handling
what S56 flagged as unresolved:

- The restriction to $\{|\omega|\geq M/4\}$ that $\sigma^*$ carries by
  definition (`def:sigma-star`) — Grujić's hypothesis, as stated, is not
  obviously restricted to this level set. Does the translation still hold
  restricted, or does it need the hypothesis unrestricted, and if so is that
  a strictly stronger requirement?
- Grujić demands the bmo condition **at every scale**, not just at $r^*$.
  Does this correspond to a scale-family of $\sigma^*$ conditions, or does
  it collapse to the single critical-scale statement this programme already
  works with? Get this right — it is exactly the kind of place a marginal
  argument goes wrong (the S31 audit's lesson, cited at every prior session
  this thread for good reason).
- State the resulting implication(s) explicitly in both directions between
  $\hat e\in\mathrm{bmo}_{1/|\log r|}$ and $\sigma^*$-decay, and between the
  hypothesis and this document's other no-concentration objects: (AX)
  (S53's structural hypothesis on the axis), and log-free saturation (S50).
  If the hypothesis turns out to be strictly stronger or strictly weaker
  than any of these — not merely "the same" — say precisely which and why;
  a genuine strict weakening would be the significant finding S56 did not
  rule out, and a genuine strict strengthening would be worth recording too.

## Deliverable 3 — the one isolable, checkable technical core

Independently check the single step S56 identified as carrying everything:
**the recasting of $\alpha=\hat e\cdot S\hat e$ as a singular-integral
commutator with a BMO symbol, and whether the Coifman–Rochberg–Weiss
estimate applies as claimed** in this specific geometric setting (not the
classical CRW setting — check the hypotheses of CRW itself and whether
Grujić's application satisfies them). **This step is valuable to verify in
isolation regardless of the rest of the paper**: if it is correct, it is a
genuine new technique available to this programme even if Theorem 7.4's
overall conclusion does not hold or does not help; if it is wrong, that is
important to know before anyone builds on this paper.

## Deliverable 4 — a clear verdict

State one of:

- **The hypothesis is (provably) equivalent to, or strictly weaker/stronger
  than, an existing object in this document** (σ*-decay, (AX), saturation).
  If equivalent: this is F2 from the audit, a legitimate fifth confirmation.
  If strictly weaker: this is a genuine, significant finding and should be
  written up as carefully as S50's equivalence theorem was.
- **The translation cannot be completed** (Deliverable 2 fails for a
  specific, stated reason) — report exactly where and why.
- **The commutator step (Deliverable 3) does not hold as claimed** —
  report exactly which hypothesis of CRW fails and why, independent of
  Deliverables 1-2's outcome.

## Forbidden, stated in advance (from S56, do not relitigate)

- **Do not attempt to prove Grujić's bmo hypothesis holds for actual NS
  solutions.** That is the open problem itself in different notation, not a
  sub-task of this pricing exercise.
- **Do not write anything into `proofs/claim_a_3d_proof_attempt.tex` that
  depends on Theorem 7.4 or its companion being correct.** Both are
  unrefereed. You may write the dictionary/comparison (Deliverables 2-3) as
  this programme's own independent result, with the source's claim reported
  as a claim, not imported as a fact.
- **Do not spend the session re-deriving Grujić's own supporting lemmas**
  from scratch before reading his proofs (S52's failure mode 1, which fired
  for real in that session). Read first, verify second, attribute correctly.

## Standing traps, carried forward

- **Session-tag collision has happened four times.** `git log --oneline -1`
  and `ls SESSION-LOG/` before writing anything; re-check at merge time.
- **Fetch primary sources; do not recall S56's summary of them.** S56 itself
  is one layer removed from the primary sources on several points.
- **If any proof-document edit is made**: `check_proof.py` (exit 0) and a
  three-pass `pdflatex` compile (0 errors, 0 new LaTeX warnings, bookmark
  warnings at or below the current baseline — check it yourself, do not
  assume the number from S56).
- Do not touch `book/`, `layer3/`, `layer4/`, `results/`. `proofs/` and
  `PROGRESS.md`/`SESSION-LOG/` only, and only if Deliverable 2 or 3 produces
  a genuine, independently-owned result worth recording — a pure "F2,
  confirmed equivalent" outcome may be recorded as a KEY FINDING in
  `PROGRESS.md` without necessarily needing new `proofs/` content, at your
  judgement.
- Do not commit or push. Report your worktree path (if any proof-document
  content was written) for independent verification before merge.
