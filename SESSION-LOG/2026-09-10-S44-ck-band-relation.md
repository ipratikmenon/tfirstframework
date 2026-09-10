# Session S44 — 2026-09-10 — Are `c_K` and Band-Absorption the Same Conjecture?

Continues the S41 core-term thread. Analytical; touches only
`proofs/claim_a_3d_proof_attempt.tex`. No code, no numerics run.

## Motivation

S41 proved the annulus level iteration cannot close and left
`target:band-absorption` as the remaining core term. S42-S43 stress-tested it
numerically (supported, not proved). The proof document repeatedly describes it as
"of the same analytic type" as S34's `conj:K-refined` (the `c_K` decorrelation
conjecture), and notes that their crude bounds fail by *exactly the same single
logarithm* (`prop:s41-crude` "reproduces, term for term, the bound of
`prop:crude-K`"). Three coincidences, and nobody had ever checked whether the
resemblance was substantive. If one implied the other, or both followed from a
common hypothesis, the program's open-item count would drop — the same kind of win
as S37's H2/σ*-decay unification.

## Orchestration

Opus subagent in an isolated git worktree (the S41 pattern), briefed with the exact
statements rather than paraphrases, and told explicitly that a negative result was
first-class and not to manufacture a link. Two rounds: an initial submission, then a
correction round after my review found a load-bearing error (below). I verified both
rounds before merging anything.

## Results

### (a), (b) — Neither implies the other. Settled.

`prop:s44-a-fails` exhibits admissible configurations with `c_K = 0` — the strongest
possible form of the c_K hypothesis — in which *every* admissible level pair of width
≥ 1/16 violates `eq:s41-bandY` by a factor `(r*/ℓ)^5 → ∞`. The adversary is `|ω|`
oscillating rapidly **inside** the band, so `|∇m|²` and `w` co-concentrate exactly
where the c_K hypothesis has no purchase.

`prop:s44-b-fails` goes the other way: configurations where both sides of
`eq:s41-bandY`/`eq:s41-bandD` vanish identically (band-absorption holds with `C=0`)
while `c_K → ∞`. So band-absorption implies **no upper bound on `c_K` whatsoever**.

The structural reason (`rem:s44-disjoint`) is simple once stated and is the real
content of the finding: **the two hypotheses are supported on disjoint regions of the
|ω|-level decomposition.** `E ⊆ {m ≥ M/4}`; `Band ⊆ {m ≤ M/4}`. They meet only on the
level set `{m = M/4}`, which carries no measure whenever `M(t)/4` is a regular value.
Neither hypothesis constrains any integral supported in the other's region. "Same
analytic type" is a statement about the **form** of the two hypotheses, not their
content.

(Independent corroboration from this program's own code: `layer4/test_band_concentration.py`
contains `test_band_and_E_are_essentially_disjoint`, written in S42 without any of this
in view, asserting exactly this disjointness numerically.)

### (c) — A common target exists, and the better one was found second.

Round 1 proposed `conj:jle-s44` (**JLE**), built on a *volume-normalised* coupling
excess `Ξ[ϖ]` that divides by the ambient `|Q(r*)|` rather than by the measure of the
region carrying the weight. `prop:s44-jle-implies` shows (JLE-1) yields K-absorption
outright — no `D_ang`, no `‖w‖_∞`, no annulus residual — and also `eq:s41-bandY`.

Round 2 produced the better one. **`conj:lc-s44` (LC)** is *verbatim* `conj:K-refined`
— same normalisation, same `C/𝓛` decay — merely asserted on `{m ≥ M/16}` instead of
`{m ≥ M/4}`. That region contains both `E` and `Band`, and the integrand is
non-negative, so `prop:s44-lc-implies` bounds both terms **with no volume hypothesis**.
JLE, by contrast, degrades like `|Q(r*)|/|S|` under concentration on a set `S`. LC is
therefore the target the program should aim at; JLE is retained only because it gives
the stronger conclusion via the budget route.

Honest non-claims (`rem:s44-lc-reading`): LC does **not** imply `conj:K-refined` —
correlation ratios are not monotone under restriction, so LC supplies the *conclusion*
both conjectures were designed to give, not either conjecture as a statement. The
converse fails too. This is the legitimate content of "same analytic type": both
targets are instances of one hypothesis on adjacent regions.

### (d) — A first-draft claim was wrong, caught in review, and inverts

**This is the part worth reading.**

Round 1 reported that `target:band-absorption` is strictly *weaker* than
`conj:K-refined` — needing only an O(1) equipartition constant where `c_K` needs the
decaying `C/𝓛` — and concluded that the program's measured `c_K ≈ 1` already supported
it directly. That would have been genuinely good news.

Verification found `lem:s44-xi-vs-cK` stated

    Ξ[1_E] = c_K·|E|/|Q(r*)| ≤ c_K        ← WRONG, reciprocal

whereas both quantities differ only in a leading volume factor, giving

    Ξ[1_E] = c_K·|Q(r*)|/|E| ≥ c_K        ← correct

Checks I ran before challenging it: (i) the algebra directly from `conj:K-refined`'s
own definition and `def:s44-excess`; (ii) a numerical check on four random
configurations — the stated form matches nothing, the corrected form matches to machine
precision; (iii) the sanity case — `|∇m|²` and `w` constant on `E` and zero off it give
`c_K = 1` (uncorrelated, as they should) but `Ξ = |Q|/|E| ≫ 1`. The error also
contradicted **two of the agent's own results** (`rem:s44-jle-status`, which correctly
computes `Ξ = |Q|/|S|` under concentration, and `prop:s44-transport`), and that internal
inconsistency is what located it.

Corrected consequences:

- A `Ξ`-bound is **strictly stronger** than the corresponding correlation bound, by
  exactly the volume ratio.
- `conj:K-refined` yields only `Ξ[1_E] ≤ (C/𝓛)·|Q(r*)|/|E|`, which is O(1) *only*
  under a volume lower bound on `E`. Round 1 had correctly flagged exactly this kind of
  hypothesis **(NB)** on the band side and missed its symmetric analogue on the bulk
  side; that is now stated as **(NB-E)**, `eq:s44-NBE`.
- **Corrected (d), via the new `prop:s44-linfty-route`:** routed through `‖w‖_∞` — the
  route requiring no measure information — `target:band-absorption` needs
  `c_K^band ≤ C/𝓛`, the *exact analogue* of `conj:K-refined` on the band. So
  band-absorption is **not weaker**, the logarithm is **not an artefact** but the price
  of that route, and in this narrow sense the document's "same analytic type" language
  is **vindicated** — precisely the opposite of what round 1 reported. The O(1)
  alternative exists only in the volume-normalised statistic, and is strictly stronger.

I re-derived `prop:s44-linfty-route` by hand: the `|R|` cancels exactly against
`∬_R w ψ² ≤ ‖w‖_∞|R|`; the enstrophy budget `eq:s41-budget` contributes `𝓛M^{1/2}`;
`M^{-3/2} = (r*)³`; and the hypothesised `C₁/𝓛` cancels the logarithm, landing exactly
in the class `prop:second-rung` absorbs. Sound.

`rem:s44-corrigendum` records the inversion and explicitly withdraws all three
first-draft conclusions, in the style of the S34 corrigendum and the S35 retraction.

### (e) — Independent as stated

Over admissible test configurations (A1)-(A3), with `rem:s44-scope` explicitly bounding
the claim: this is independence over the budgets and pointwise structure the two
statements are phrased in, **not** over actual Navier-Stokes solutions.

### (c′) — A correction to the surrounding text

`rem:s44-bandD`: `eq:s41-bandD` (the D-half of the closure target) contains no product
of two densities and no correlation at all — it is a bare level non-concentration
statement for `ν|∇²ê|²ψ²`. **No decorrelation hypothesis in this family bears on it.**
So "same analytic type as `conj:K-refined`" is accurate for the Y-half and *inaccurate*
for the D-half. Future sessions should not assume the D-half comes along for free.

## Net effect on the program

The open-item count is **not** reduced — both conjectures remain open, and are now
proved independent rather than suspected equivalent. What improved is the targeting:
one sharply-stated hypothesis (LC) covers both Y-halves instead of two separate ones,
the D-half is explicitly flagged as covered by nothing, and the `O(1)` vs `C/𝓛`
confusion is resolved (it was a route distinction, not a strength distinction).

## What was NOT done, and what to measure next

`rem:s44-numerics` names three quantities that bear on this and were **not** measured:
`|E|/|Q(r*)|` and `|Band|/|Q(r*)|` **at scale r\***; `c_K^band`; and the `𝓛`-dependence
of both correlation constants — since `prop:s44-linfty-route` needs *decay*, not
boundedness, and S34-S35 explicitly could not distinguish `c_K = O(1)` from
`c_K ≤ C/𝓛`. It also warns that S42/S43's recorded volume fractions are taken over the
whole periodic box rather than over `Q(r*)`, so they are motivation for that
measurement and **not** evidence about it.

## Verification I performed (not merely relayed)

- Re-derived the corrected `lem:s44-xi-vs-cK` and `prop:s44-linfty-route` by hand.
- Independent numerical check of the `Ξ`/`c_K` identity on four random configurations.
- Brace balance (depth 0, never negative); `\begin`/`\end` multiset equality; zero
  duplicate labels; zero **new** unresolved refs (the same four pre-existing ones
  remain: `eq:NS-vorticity`, `lem:GN-ehat`, `prop:blowup-alignment`,
  `thm:vorticity-uniform`); all S44 environments declared in the preamble; no
  recurrence of the malformed `\end{env>` typo class the agent self-reported fixing.
- Confirmed the merge is strictly additive against the live file: 1192 insertions,
  **0** deletions or alterations, with S41 intact.

**Not** independently verified: the internal constants of the two counterexample
constructions (`prop:s44-a-fails`, `prop:s44-b-fails`) were re-checked line by line by
the agent but not re-derived from scratch by me. Their *structural* basis — the
disjointness of `E` and `Band` — I did verify, and it is independently corroborated by
this program's own S42 test suite.

## Process note

The S41 session established the bar: isolate in a worktree, then independently
re-derive the load-bearing computations before merging. That bar is what caught this
one — a single inverted volume factor that would otherwise have entered the proof
document as a favourable-sounding headline result and propagated into future sessions'
reasoning about which conjecture to target.

## Files changed

- `proofs/claim_a_3d_proof_attempt.tex` (+1192 lines, one new section
  `sec:relation-s44` plus minimal S41 cross-references)
- `PROGRESS.md` (KEY FINDING S44, session pointers, Session Log row, date)
- `SESSION-LOG/2026-09-10-S44-ck-band-relation.md` (this file)

No Python, no `results.db`, no test suites affected.
