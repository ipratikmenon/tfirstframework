# Session S38 — 2026-07-19 — Regime A Impossibility, σ* Trend Replication, Closure Target

## Purpose

Attack σ*-decay directly, per the two S37-recommended directions. Direction 2
("local Regime A") was checked FIRST as a fact-check before deriving anything —
and found to be a dead end, cleanly. Direction 1 (assembling the proved Caccioppoli
chain into a genuine recursion) was deliberately NOT force-closed, given the
session's demonstrated error pattern on exactly this kind of multi-term scaling
argument (S31, S35→S36 correction, S37 self-correction).

## Analytical result — §closure-s38 (proof doc, now 6667 lines, balanced)

**Proposition (Regime A is impossible) — PROVED, standard fact.** Via the
standard H¹ continuation criterion (Fujita–Kato local well-posedness: solutions
with bounded H¹ norm extend for a time depending only on that norm), any genuine
blowup at time T forces limsup_{t→T⁻}‖u(t)‖_{H¹} = ∞; since ‖u‖_{L²} is
non-increasing (energy inequality), this forces ‖∇u‖_{L²} → ∞, hence (Biot–Savart
isomorphism on T³) Ω(t) → ∞. **Consequence:** Regime A (bounded enstrophy) is
IMPOSSIBLE for any genuine singularity — not merely unlikely. This closes off
"local Regime A" as an independent escape route entirely. σ*-decay is not one of
several paths to H2 — given H1, it is NECESSARY.

**Honest numerical remark (non-monotonic σ*, not simple decay or growth):**
recorded the S37-adjacent finding that σ* rises then partially recedes as M
approaches its true peak in the TG run, with the caveat that late-time (post-peak,
M-decaying) σ* growth is likely a definitional artifact of a less-selective
threshold set, not physical.

**Closure target — deliberately left OPEN, not forced.** Assembling the proved
first-rung + second-rung + amended-K-absorption chain (S33–S36) into a genuine
scale recursion for σ*(ρ) as ρ→0: whether the resulting effective coefficient is
subcritical (<1, self-improving regardless of starting σ*) or critical/
supercritical (≥1, needs a priori smallness — same fate as §19/§20) is NOT
determined. A preliminary scaling check (M-power consistency between the K-term
and the first-rung dissipation) is necessary-but-not-sufficient and does not
settle the question. Recorded honestly as unresolved — precisely to avoid a fourth
mid-session correction on this exact failure mode.

**Literature pointer:** Type-I blowup-profile compactness (Nečas–Růžička–Šverák
1996, Escauriaza–Seregin–Šverák 2003) is named as the natural longer-range tool
for attacking σ*-decay at an ACTUAL singularity — distinct from what generic
numerics (which by construction samples only non-singular dynamics) can test.

## Numerical result — σ* trend replication across existing runs

Pure data-mining pass (no new solver runs) across all 6 logged bridge experiments.
Only 2 have real M dynamic range: TG-001 (N=32, M ratio 12.4×) and TG-002 (N=64,
M ratio 15.0×) — shear and adversarial runs lack sufficient growth phases to test.

**Both TG runs replicate the same specific signature:** σ* hits a local maximum
well before M peaks, then drops sharply (1.88× for TG-001, 2.42× for TG-002) as M
reaches its true maximum, and both show strong negative σ*-M correlation in the
decay phase (−0.94, −0.74). Growth-phase correlations are also consistent between
runs (weak positive, 0.37 vs 0.39).

**Honest scope:** confirmed within a family of n=2 TG runs at different
resolutions — this rules out "N=64-specific artifact" but does NOT establish
generality across ICs. Shear/adversarial are uninformative (no counter-evidence,
simply no dynamic range), not confirmatory. More IC diversity with real M growth
needed before calling this general.

## Actions

- proofs/claim_a_3d_proof_attempt.tex: §closure-s38 inserted (balanced, all refs
  verified, 2 new bibliography entries added).
- No new solver code — this session was pure derivation + data mining.
- PROGRESS.md update dispatched.
- No commit — not a git repository.

## Where the program stands after S38

The picture is now maximally sharp: **exactly one deep, unavoidable open
problem** — forcing σ*(t)→0 near a putative Type-I singularity — plus one
mechanical bookkeeping step (the annulus iteration). No remaining escape routes.
The closure attempt (assembling proved pieces into a decay mechanism) remains
genuinely open and is honestly labeled as such rather than forced. Numerics
consistently shows generic (non-singular) dynamics keeps σ* away from small
values, which is consistent with — not contradictory to — the disorder-barrier
picture: if blowup requires σ* small and ordinary turbulence doesn't produce
small σ*, a genuine singularity (if any) would have to be a special, rare
configuration.

## Next session (S39) — recommended direction

1. Complete the closure-target bookkeeping (the S38-flagged open item) with full
   care: nondimensionalize every term in the assembled first-rung/second-rung/
   K-amended chain and determine sub vs supercritical, learning explicitly from
   the S31/S35 failure modes (check every Hölder exponent, every Jensen
   direction, verify no circular regularity assumption).
2. If the closure turns out supercritical (as §19/§20 did), pivot fully to the
   Type-I compactness literature route rather than another scale-recursion
   attempt — a different tool for a problem that may need one.
3. Optional: broaden the σ* numerical dataset with a diverse-IC growth-phase run
   (something between shear and TG in character) to test whether the hump
   pattern is TG-specific or general.
