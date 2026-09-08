# Session S37 — 2026-07-19 — Localized CF Depletion, H2/σ*-Decay Unification

## Purpose

Derive H2 (localized Constantin–Fefferman depletion) as a consequence of H1 (Bridge
Lemma), using the good-sign structures found in S36. Continued immediately from S36
in the same session per user instruction ("yes for next step already").

## Analytical result — §cf-s37 (proof doc, now 6507 lines, 511/511 balanced)

**Derivation (Constantin's exact geometric kernel, classical near/far optimization):**
using the pointwise cancellation D(e,e,σ)=0 in Constantin's stretching-rate
representation, Lipschitz continuity of ê (from H1's Bridge Lemma output) controls
the near field, and Cauchy–Schwarz against the global enstrophy Ω=‖ω‖²_{L²}
controls the far field. Optimizing the cutoff radius gives

    α(x*) ≲ (KM)^{3/5} Ω^{1/5},   K = C_B(σ*)^{1/2}/r* = C_B(σ*)^{1/2}M^{1/2}.

**Self-audit catch during derivation:** my first-pass sketch assumed Ω = O(1)
(bounded enstrophy) for free — WRONG, since bounded enstrophy near a putative
singularity is not free (it's essentially assuming no enstrophy blowup, which for a
genuine Type-I profile with self-similar scaling ω ~ (T-t)^{-1}Ω₀(y/√(T-t)) actually
gives Ω(t) ~ M(t)^{1/2}, NOT O(1)). Redone correctly with this scaling substituted:

**Regime A (Ω = O(1)):** α(x*) ≲ C(σ*)^{3/10}M^{9/10} — θ=2/5<1/2, UNCONDITIONAL
depletion given H1, using only the a priori bound σ* ≤ C𝓛/ν (S34), no decay needed.

**Regime B (Ω ~ M^{1/2}, self-similar/generic):** α(x*) ≲ C(σ*)^{3/10}M — exactly
borderline (θ=1/2). The near/far optimization gains NOTHING beyond a constant
prefactor; strict depletion needs σ*(t) → 0, not merely bounded.

**Central finding — unification:** in the generic (Regime B or worse) case, H2's
remaining content and the standalone σ*-decay gap are THE SAME open problem. The
program's open analytical content reduces from three nominally independent items
(H1's annulus step, H2, σ*-decay) to TWO: the annulus level-iteration (mechanical)
and a single deep gap — a mechanism forcing σ*(t)→0 along any putative Type-I
blowup sequence. Regime A remains a genuine alternative escape route (a bounded-
enstrophy argument, independent of σ* entirely) if it can be established.

## Numerical result — enstrophy growth exponent

New layer4/enstrophy_exponent.py (399 lines, 18 tests; program 319/319, zero
regressions). Measured p in Ω(t) ~ M(t)^p on the TG stretching run (M grows 15×,
same seed/params as EXP-L4-DEPLETION-TG-001 for direct comparability):

- p_fit_all (M>5, n=33): 1.226, r²=0.624
- **p_fit_growth_phase (pre-peak stretching only, n=19): 0.891, r²=0.914**
- ADV run: uninformative (M never grows enough; correctly reported as NaN/2-point
  artifact rather than forced)

**Honest interpretation (important, must not be overstated):** p≈0.89 lands ABOVE
both Regime A (p≈0) and Regime B (p≈0.5) — nominally "worse than borderline." BUT
this measures ordinary transient vortex stretching in a globally-regular, DECAYING
flow (M peaks then falls), not the approach to an actual singularity — which by
definition cannot be produced by a solver known to stay regular. The measurement is
real and well-fit (r²=0.91) but its relevance to the Regime A/B question — which is
specifically about behavior NEAR a hypothetical blowup — is limited. Correct
reading: this is mild evidence AGAINST expecting Regime A (bounded enstrophy) to be
a generic free lunch even in ordinary dynamics, which sharpens the case that
σ*-decay (not enstrophy-boundedness) is the program's real, unavoidable core
difficulty — a clarifying, if sobering, result. Does not directly test or falsify
the theory.

## Actions

- proofs/claim_a_3d_proof_attempt.tex: §cf-s37 inserted (511/511 balanced, all refs
  verified, Constantin1994 bibitem added).
- layer4/enstrophy_exponent.py + tests new; 2 experiments logged.
- N=128 TG resolution sweep (from S36) completed during this session: TG-101
  (ν=1e-3) c_K=2.235, TG-102 (ν=5e-4) c_K=2.238 — both notably higher than the
  under-resolved N=64 values (~1.02-1.06), confirming N=64 TG numbers WERE
  resolution-biased as flagged; TG-103 pending at session-log time.
- PROGRESS.md update dispatched.
- No commit — not a git repository.

## Where the program stands after S37

The genuinely open analytical core is now precisely ONE deep problem: forcing
σ*(t)→0 along a putative Type-I blowup sequence (equivalently, defeating Regime B).
Everything else in H1's and H2's estimate layers is proved or mechanical. This is
real progress in problem STRUCTURE even though the deep problem itself remains
unsolved — the program went from "three separately-hard open items" to "one hard
item plus bookkeeping."

## Next session (S38) — recommended direction

Attack σ*-decay directly. Two candidate mechanisms worth deriving:
1. Use the Gram-negativity structure found in S36 (good-signed drift-curvature
   coupling) directly in the dσ*/dt evolution itself, rather than in the K-coupling
   — the S31-audited dσ/dt inequality was wrong; a redone version using the S34-S36
   toolkit (magnitude Caccioppoli, Gram negativity) may fare better.
2. Investigate whether Regime A (bounded local, not global, enstrophy near x* at
   scale r*) is achievable independently of σ*, using the same localized
   Caccioppoli machinery already built (S33's local enstrophy-gradient bound) rather
   than the global Ω used in this session's derivation — a local-enstrophy version
   of Regime A might be attainable without touching σ* at all.
