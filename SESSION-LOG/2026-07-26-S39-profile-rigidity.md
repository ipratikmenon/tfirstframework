# Session S39 — 2026-07-26 — Profile Rigidity (Type-I Compactness Reframe)

Orchestration note: at the user's direction the analytical derivation was delegated to
an **Opus** subagent (owning the proof document exclusively), with Sonnet tracks in
parallel on new layer4 files. Fable specified the four-step program, named the two
expected failure modes in advance, verified the returned mathematics, and wrote this
log. See the companion log 2026-07-26-S39-aposteriori-verification.md for the
numerical track.

## Motivation

Every direct attempt to force σ*(t) → 0 as a *dynamical decay* statement has failed the
same way (§19, §20 per the S31 audit; the S38 closure target left open). The reframe:
use Type-I blowup-profile compactness to convert dynamical decay into **structural
rigidity** of a limit object, where σ* is an invariant rather than a decaying quantity.
S38 had already proved σ*-decay is *necessary* (Regime A impossible), so this is the
only remaining line.

## New section: §sec:rigidity-s39 (proof doc 6667 → 7737 lines, 582/582 balanced)

### Proved

- **lem:typeI-scaleinv-bounds** — the Type-I ancient class 𝒜(c_I) is closed under the NS
  scaling group with the *same* constant; σ*(U_μ)(s) = σ*(U)(μ²s). This is what makes a
  second compactness extraction legal later.
- **prop:transfer-audit(a)–(c)** — direction equation, magnitude equation, orthogonal
  split, Gram negativity all hold on any profile. Important framing: these are
  **re-derived from scratch on U**, not transferred as limits of inequalities (U is a
  smooth solution with scale-invariant Type-I bounds; every proof is local/pointwise).
- **lem:sigma-semicontinuity** — lower semicontinuity σ*(U) ≤ liminf σ*(u_λj), with the
  explicit statement that no reverse inequality holds.
- **lem:alignment-local** — σ*(U)(s) = 0 ⟹ ê constant on the relevant component of the
  high-vorticity set.
- **lem:harmonic-reduction** — *the session's genuinely new mathematics.* If ω_U ∥ e₃
  globally then ∂₃U ≡ 0, U₃ = U₃(s), and U_h solves 2D NS exactly. Proof: div ω = 0 gives
  ∂₃ζ = 0, so f := ∂₃U is curl-free and divergence-free, hence componentwise harmonic;
  Type-I gives boundedness, so Liouville forces f = f(s); then U = f(s)y₃ + Ũ with U
  bounded in y₃ forces f ≡ 0. **Verified independently by Fable, all five steps.** This is
  the step that makes any Liouville theorem applicable at all — general 3D
  bounded-ancient Liouville is open, so no version of the aligned horn can work without a
  real dimensional reduction.
- **prop:no-summable-budget** — the naive infinite-past budget argument is **provably
  impossible**: dyadic past shells and self-similar cylinders are commensurate, and every
  cylinder carries the same dimensionless bound with no k-decay (identical, for an exactly
  self-similar profile, by scale invariance), so Σ_k ⟨σ*⟩ = ∞ is forced by the class.
  This converts the trap Fable named in advance into a theorem, and dictates the only
  viable framework: a bounded quantity **monotone in logarithmic time**.
- **prop:logtime-necessary** — with τ = log|s| and Φ(τ) = |s|M_U(s) ∈ (0, c_I], the
  magnitude equation at the vorticity maximum (Hamilton's trick, Δ|ω| ≤ 0 there) gives
  d(logΦ)/dτ ≥ 1 + |s|w_max − |s|α_max; boundedness of Φ over infinite τ-length yields
  limsup of the log-time average of [|s|α_max − |s|w_max] ≥ 1. **Computation re-derived
  and confirmed independently by Fable.** Good sign structure: disorder drives Φ up, and
  Φ is bounded.

### Proved modulo stated hypotheses

prop:profile-extraction (mod R, N); prop:transfer-audit(d),(e) (mod R) — including the
notable finding that **profiles are logarithm-free**: the 𝓛 = log(e+‖u‖_{H³}/M) factor
degrading every estimate in §S33–§S36 becomes an absolute constant C(c_I) on the profile,
since scale invariance leaves the log-divergence nowhere to come from; lem:sigma-equality
(mod R, N); lem:harmonic-reduction (mod R); thm:aligned-horn (mod R, A-up, plus a flagged
mild-solution class caveat on the cited KNSS theorem); lem:logtime /
prop:logtime-necessary (mod R and attainment of the spatial supremum, with two honest
repairs given, one adopted); prop:two-step-reduction (mod R, N′).

Hypotheses: **(R)** scale-uniform interior regularity of the rescaled family (⟹ C^∞_loc
convergence); **(N)/(N′)** non-degeneracy/nontriviality of the limit; **(A-up)** global
alignment upgrade.

### Step-2 verdict (the question Fable flagged in advance as the crux)

The σ* transfer **works in the direction the argument needs — but the cost was relocated
into an explicitly labelled hypothesis, not solved.** Specifics:
1. The derivative-level obstruction does not appear as a defect inside Step 2; it appears
   as hypothesis (R). With (R), σ* transfers with **equality** (the indicator
   1_{|ω|≥M/4} bounds the denominator of |∇ê|² = (|∇ω|² − |∇|ω||²)/|ω|² away from zero).
   Without (R): lower semicontinuity with a defect measure.
2. Lower semicontinuity is the **useful** direction — smallness descends to the limit,
   which is what prop:two-step-reduction uses to absorb the oscillatory middle case
   (liminf σ* = 0 but never zero) into the aligned horn via a second compactness.
3. **Neither horn is broken by the direction**, because the dichotomy is intrinsic to U.
   What semicontinuity kills outright is the hybrid strategy "σ*(U) = 0 ⟹ σ*(u_λj) ≤ ε₀
   ⟹ apply H1/H2 to the original sequence" — no reverse inequality exists. Recorded as
   an error-marker at rem:s39-nontransfer(ii).
4. The horn that genuinely depends on (R) is **Step 4 (disordered)**, not Step 3: with
   only weak convergence the disorder could live entirely in the defect measure, so the
   disordered horn would be attacking the wrong object.

### New obstruction found (constrains prior work)

**rem:s39-nontransfer(i): prop:localized-cf (§S37) does NOT survive passage to the
profile.** Its far-field step needs ‖ω‖_{L²(ℝ³)} < ∞, which an ancient Type-I solution
need not satisfy; the crude substitute |ω_U| ≤ M_U makes the far-field kernel integral
logarithmically divergent at large |y|. So the entire §S37 H2/Regime-B analysis is
unavailable on profiles as written. This is the honest reason Step 4 does not close, and
any future session invoking prop:localized-cf on a profile is in error.

### Open targets

- **(A-up), the global alignment upgrade — the largest new gap.** Local, one-time,
  one-component alignment must be upgraded to global ω ∥ e₀ on ℝ³×(−∞,0), because the
  harmonic Liouville step has no ball version. Three mechanisms discussed
  (rem:A-up-discussion): unique continuation for the w-inequality (reduces to the
  program's own K-family Kato-class problem — mildly positive); ESS backward uniqueness
  (lacks the spatial decay); scale exhaustion (needs σ* ≡ 0 for *all* s, which the
  dichotomy does not deliver).
- **target:disorder-depletion** — the disordered horn's remaining inequality, sharply
  posed: persistent σ* ≥ δ must force the log-time average ≤ 1 − c(δ). Flagged as
  **marginal** (must beat the constant 1 exactly — the regime where §19/§20 errors hid);
  mitigating difference is that it compares two averages by a fixed constant rather than
  being a recursion needing bootstrapped smallness. Decomposed into (T4a) average-to-point
  for w at the maximum (natural mechanism: weak Harnack, plugs into the existing K-family
  machinery) and (T4b) profile depletion of α (the harder half, and now known to need a
  genuinely different far-field control per the obstruction above).

## Toolchain finding (process, important)

There is **no TeX toolchain on this machine** — the proof document has never been
compiled. The Opus agent found that `conjecture` and `fact` were used as environments
across several prior sessions without ever being declared via `\newtheorem` (3 of the 5
sites predate S39), so the file would have failed on first compile. Declarations added at
preamble lines 33–34. Every prior "verified" check in this program was grep-based brace
counting, which cannot catch this error class. **Install LaTeX and do a real compile
before any arXiv step.**

## Distance to the goal (unchanged in kind, sharper in detail)

Even with both horns closed and (R),(N),(N′),(A-up) discharged, the conclusion is "no
Type-I singularity at a non-degenerate point" — strictly less than the Prize. Type-II
exclusion is now **more** decisively outside these methods: without Type-I rates the
rescaled family has no uniform bound, so there is no limit object at all and the reframe
is vacuous against Type-II. (N) is not a technicality either: CKN ε-regularity gives the
right scale-invariant lower bound on Q₁(0), but that cylinder touches the terminal slice
s = 0 where compactness fails, so pushing mass off it needs a Type-I-improved
ε-regularity theorem (Seregin) that is not proved here.

## Next session (S40) candidates

1. **(A-up)** via unique continuation for the w-inequality — the most promising of the
   three mechanisms, and it reuses the program's existing K-family machinery.
2. **(T4a)** weak Harnack for average-to-point control of w at the vorticity maximum.
3. **Install a TeX toolchain and compile the document** — non-negotiable before any
   external submission; likely to surface more latent errors of the `\newtheorem` class.
4. Do NOT attempt (T4b) using prop:localized-cf — proved unavailable on profiles.
