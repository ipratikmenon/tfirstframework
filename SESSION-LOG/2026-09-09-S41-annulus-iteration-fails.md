# Session S41 — 2026-09-09 — Annulus Level Iteration: Executed, Found Not to Close

On the S39/H1 critical path (unlike S40, which was an independent numerical side-track
earlier the same day). Delegated to an Opus subagent per user direction to pick a
specific, well-scoped item and attack it: the "annulus mechanical step," the last
technical item standing in H1's estimate layer per the S36 consolidated status table.

## Orchestration

Given this program's history of over-claimed results later withdrawn on audit (S31:
arithmetic-impossible induction, reversed Jensen's inequality, circularity; S35→S36:
a claimed "bonus" term retracted after self-audit), the subagent was launched in an
**isolated git worktree** rather than editing the shared proof document directly, and
was explicitly instructed to (a) trace every symbol back to its real definition rather
than guess notation, (b) actually carry out the convergent-series computation behind
"geometric beats geometric" rather than assert it, (c) report a clean failure if the
argument didn't close rather than paper over a gap, and (d) stay scoped to only this
one item.

After the subagent returned its self-report, **I independently reviewed the actual
diff before merging anything** — not just the self-report:

- Programmatically checked brace balance (final depth 0, never negative) and that the
  multiset of `\begin{...}`/`\end{...}` environments match exactly.
- Confirmed every environment used is either a standard LaTeX/amsmath environment or
  declared via `\newtheorem` in the preamble (lines 23-34), except two pre-existing
  undeclared environments (`observation`, `openproblem` at lines 5039/5117/5144) that
  the subagent correctly flagged as **not its own** — verified these predate this
  session's diff (S39 found this general error class; worth a separate fix someday,
  not done here, out of scope).
- Confirmed `git diff --stat` shows exactly one file changed, 864 insertions / 6
  deletions, all inside three locations: the `rem:level-iteration` addendum, the S36
  consolidated status table, and the new section — nothing else in the 8600-line
  document touched.
- **Hand-re-derived the two load-bearing computations myself**, rather than trusting
  them: the Young's-inequality step `4XZ ≤ δX² + 4δ⁻¹Z²` (verified via
  `2ab≤a²+b²` with `a=√δ X`, `b=2Z/√δ`: gives `δX²+4Z²/δ ≥ 2·√δX·2Z/√δ = 4XZ` ✓), and
  the `Λ_k = 2^8·4^k + 2^8·δ_k^{-1}·4^k` formula, confirming the "no contraction for any
  δ_k" conclusion follows immediately and correctly once that formula is granted.
- Read the full new section once as a skeptical referee, watching specifically for the
  error classes previous audits caught (invalid Hölder exponents, reversed Jensen,
  circular reasoning, dimensional inconsistency, unjustified steps). Found none.

Only after this independent pass did I copy the file from the worktree into the main
tree and remove the worktree.

## What was proved

New `\section{Execution of the Level Iteration for the Annulus Residual: an Exact
Inventory and a Structural Obstruction (Session S41)}`, `\label{sec:annulus-s41}`,
lines ~7666-8480 of `proofs/claim_a_3d_proof_attempt.tex` (7737 → 8600 lines).

**Setting:** Proposition `prop:K-amended` (S36) bounds the K-coupling term with an
absolute constant, no logarithm, no correlation hypothesis — *modulo* an annulus
residual `\mathcal{A}_{\mathrm{ann}}` (Definition `def:annulus-residual`), cutoff-
derivative terms supported on the transition band `{M/8 ≤ |ω| ≤ M/4}`. Remark
`rem:level-iteration` (S36) sketched, but did not execute, a De Giorgi nested-level
scheme meant to dispose of it, classifying it "open-mechanical... it involves no new
estimate."

**Positive result:** an exact inventory of the 8 cutoff-derivative families arising
from the weighted magnitude identity shows 4 (labeled G1-G4) are good-signed and can
be discarded outright, and 2 more (I1-I2) fall in the inhomogeneous class already
allowed by `prop:K-amended`, via the level-free enstrophy budget. This leaves an
explicit irreducible core of exactly 2 families (S1, S2).

**Negative result (the main finding):** executing the sketched iteration on this core
gives a recursion `Y_k ≤ G_k + Λ_k·Y_{k+1}` with
`Λ_k = 2^8·4^k + 2^8·δ_k^{-1}·4^k`. This is proved (`prop:s41-nocontraction`) to satisfy
`Λ_k ≥ 2^8·4^k > 1` for **every** choice of the Young parameter `δ_k` — the (S1) family
is not a product of a dissipative factor with anything else, so Young's inequality has
no purchase on it regardless of parameter choice. The originally proposed
`δ_k = δ_0·8^{-k}` makes it strictly worse (`Λ_k ~ 32^k`), and demanding
`Σδ_k ≤ 1/2` (needed elsewhere in the argument) forces `δ_k^{-1} ≥ 2`, so both
geometric factors reinforce rather than cancel — "geometric beats geometric" has no
realization here. Unrolling the recursion (`prop:s41-vacuous`) gives an inequality
that is weaker than the trivial `Y_0 ≤ Y_0` for every N, because the remainder
coefficient `∏Λ_i` diverges while `Y_N ≥ Y_0` by monotonicity.

**Root cause, identified precisely** (`rem:s41-whichlemma`): neither classical
iteration lemma applies. The linear (Ladyzhenskaya-Uraltseva / Giaquinta-Giusti
hole-filling) lemma needs a contraction factor `< 1`, which is proved impossible here.
The nonlinear (Stampacchia-De Giorgi) lemma needs the coarse level bounded by the fine
one with a superlinear gain and small initial data — three simultaneous failures:
the recursion runs coarse-to-fine (opposite of standard De Giorgi truncation, where
superlevel sets shrink to measure zero), it's homogeneous of degree one (no superlinear
gain), and smallness of `Y_0` is the sought conclusion, not an available hypothesis.
The structural reason: De Giorgi truncation gains because the truncated sets shrink;
here the truncation is in an **auxiliary variable** (`|ω|`), confined to a fixed band —
no set shrinks, no measure decays, nothing for a nonlinear gain to feed on.

**Confirmed by an independent argument** (`prop:s41-coarea`, co-area identity): as
`k → ∞`, the nested cutoffs converge weakly-* to a Dirac mass at the boundary level,
so the residual converges to a literal **flux through the level set** `{|ω| = M/8}` —
exactly the pointwise information H1 is trying to produce in the first place. More
generally, by the scale-free normalization `∫|χ'| ds = 1`, *any* admissible level
cutoff produces a residual that is a level density of the underlying measure, bounded
only by a fixed multiple of the band mass — a clean pigeonhole check confirms no
redistribution of levels (uniform, geometric, or optimally chosen) can beat this: with
N equally-spaced levels, `|χ'| ≤ 16N` against a best-annulus mass `≤ 1/N` cancels
exactly.

**What survives unconditionally:** the core residual admits a crude bound
(`prop:s41-crude`) that reproduces `prop:crude-K`'s bound term-for-term — i.e. it is
short of closure by **exactly the same single logarithm** already identified (S34) as
obstructing K-absorption. This is not a new kind of obstruction; it's the known one
resurfacing. The end target is stated precisely (`rem:s41-endtarget`,
eq:s41-endtarget) and restated as `Conjecture target:band-absorption` — a
decorrelation/equipartition hypothesis on the transition band, of the same analytic
type as the existing `Conjecture conj:K-refined` (the `c_K` conjecture, numerically
supported at `c_K ≈ 1` in S34-S35).

## Net effect on the program

`\mathcal{A}_{\mathrm{ann}}` is reclassified from "open-mechanical" (S36's belief that
it was free bookkeeping) to genuinely open (Conjecture target:band-absorption). H1
still rests on exactly two things — (a) `σ* ≤ ε₀` and (b) this conjecture — unchanged
in *count*, but (b) is now known to be real analytical content rather than a
formality. This is not progress toward closing H1, and is reported as such; it is,
however, a correction of a false expectation, and the honest thing to log. Nothing
else in the document is touched: H2, σ*-decay, Type-II exclusion, the withdrawn
§20-22 material, and the overall Prize/Claim-A conclusion are all unaffected.

## Process note

Given the program's history (S31 audit, S35 retraction), this session establishes the
standing bar for future analytical work delegated to a subagent: isolate in a
worktree, then independently re-derive the load-bearing computations by hand before
merging anything into the shared proof document — a self-report, however carefully
written, is not itself verification.

## Files changed

- `proofs/claim_a_3d_proof_attempt.tex` (+864/-6 lines, one new section + two amended
  locations)
- `PROGRESS.md` (new KEY FINDING S41 block, Last/Second-to-last session pointers,
  Session Log table row)
- `SESSION-LOG/2026-09-09-S41-annulus-iteration-fails.md` (this file)

No Python/layer1-4 code touched; no test suite affected (pure LaTeX/derivation
session).
