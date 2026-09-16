# S62 — CD-2: fractional forcing-floor extension (negative, H2) + second-instance search

Date: 2026-09-16. Follow-up to S60/S61 (`HANDOVER-S61-CD2.md`), owner
direction: "Yes please." Executed by a subagent under full orchestrator
supervision; every claim below was independently re-verified by the
orchestrator against primary sources before being recorded here (see
Verification section).

## Scope reminder

Still Thread 1 of the (C)/(D) line only. This session studies forcing
requirements in *other people's* published constructions; it does not
introduce forcing into this programme's own arguments anywhere, and it is
not advocacy for (C)/(D) as a goal. Thread 2 (vanishing viscosity) stayed
closed per S60, not reopened.

## Deliverable 1 — fractional extension of the swirl maximum principle: FAILS (H2), located precisely

For a pure-swirl axisymmetric field `u = u_θ(r,z) e_θ(x)`, the azimuthal
component of the fractional Laplacian `(-Δ)^{α/2}u` (defined component-wise
in the fixed Cartesian frame, `0<α<2`) works out to

```
g_θ(x) = c_{3,α} P.V. ∫ [u_θ(x) − u_θ(y) cos(θ_x−θ_y)] |x−y|^{-3-α} dy
```

(independently re-derived by hand from the definition and `e_θ(x)·e_θ(y) =
cos(θ_x−θ_y)`, and confirmed to match the subagent's formula exactly). Unlike
the classical Córdoba–Córdoba scalar setting — the actual PNAS 2003 paper
(via PMC307564) was fetched and its Theorem 1, `2θΛ^αθ ≥ Λ^α(θ²)`, confirmed
to apply only to **scalars** — the bracket here is not `u_θ(x)−u_θ(y)`; it
carries the angle-dependent, sign-changing weight `cos(θ_x−θ_y)`.

The classical (`α=2`) proof works because the maximum principle is applied
to `Γ = r u_θ`, not `u_θ` directly, and the vector-Laplacian identity
`(Δu)_θ = Δu_θ − u_θ/r²` converts the equation for `u_θ` into a genuine
closed scalar equation for `Γ`. The counterexample found (and independently
re-checked by hand here): fix `x` at an interior maximizer of `Γ`, take `y`
at the **same angle** (`cos(θ_x−θ_y)=1`) with `r_y` slightly less than
`r_x`. Since `Γ'(r_x)=0` at the max, `Γ(y) ≈ Γ(x)` to leading order, but
`u_θ(y) = Γ(y)/r_y > Γ(x)/r_x = u_θ(x)` purely because `r_y<r_x` — so the
bracket `u_θ(x)−u_θ(y)` is **strictly negative** even at zero angular
separation. The mechanism that makes the classical case work (the local,
`y→x` degeneration of the angular weight into one algebraic curvature term
that cancels exactly at a critical point of `Γ`) has no analogue once `α<2`
makes the operator genuinely long-range: `g_θ` cannot be re-expressed as a
nonlocal operator on `Γ` alone with a sign property at `Γ`'s maximum.

**Verdict: H2, precisely located** — the swirl maximum principle's classical
proof mechanism fails for `0<α<2` because the natural Cartesian-component
comparison quantity is `u_θ`, not `Γ`, and these need not share a
maximizer. This is a negative result about *this specific proof route*, not
a claim that no fractional swirl maximum principle exists by any means (an
added radial-monotonicity hypothesis on `u_θ` would block the
counterexample). No prior result addressing this specific question was
found in the literature search run. Because D1 fails, no fractional
`cor:cd1-rate` exists to apply in D2.

## Deliverable 2 — second-instance search

Two candidates named in `HANDOVER-S61-CD2.md`, both checked from primary
sources and both **ruled out on axisymmetry grounds** (not existence):

- **arXiv:2309.08495** (Córdoba–Martínez-Zoroa, forced 3D Euler,
  `C^{1,1/2-ε}∩L²` force): confirmed to exist and be the paper meant. Its own
  abstract states "a specific focus on **non-axisymmetric** solutions,"
  contrasting against axisymmetric-without-swirl work by others. Does not
  satisfy `thm:cd1-general-floor`'s hypotheses.
- **arXiv:2407.06776** (Córdoba–Martínez-Zoroa–Zheng, hypodissipative forced
  NS, already used at S60 for its `α₀≈0.0927` threshold but never checked
  for axisymmetry against Thread 1): fetched the full PDF directly; §1.2.5
  "Symmetry of the construction" states verbatim that the vorticity is
  **axial** (even in its own coordinate, odd in the other two — reflection
  symmetry across coordinate planes) and the velocity **polar** — a Cartesian
  vortex-layer symmetry, not azimuthal rotation. Confirmed non-axisymmetric,
  independent of the fractional-dissipation question. This completes, rather
  than corrects, S60 — S60 never claimed this paper was axisymmetric, it
  simply never checked applicability to the Thread-1 theorem at all.

**New instance found, not named in the handover: Qi S. Zhang, arXiv:2311.12306,
"A blow up solution of the Navier-Stokes equations with a super critical
forcing term."** Confirmed genuinely axisymmetric, pure-swirl, forced: its
system (1.4) is literally the ASNS system `lem:wpt3-forced-swirl`/
`lem:wpt3-forced-max` are built for, `ν=1`, finite cylinder domain (making
the geometric compactness hypothesis automatic). Theorem 1.1 gives two
explicit closed-form self-similar solutions.

Independently re-derived by the orchestrator (not just re-checked): with
`v_θ(r,t) = u_θ(r,t) + αr`, `u_θ(r,t) = φ_0(ρ)/√(2(T−t))`, `ρ = r/√(2(T−t))`
(paper's eq. 1.12–1.15),

```
Γ(r,t) = r v_θ = r u_θ + αr² = ρ φ_0(ρ) + 2α(T−t)ρ²
```

using `r = ρ√(2(T−t))`. The first term is an exact function of `ρ` alone
(no residual `t`-dependence); the second is an explicit `O(T−t)` correction
that vanishes at any fixed `ρ` as `t→T`. Since `φ_0(r) = O(r/(r²+1))` decays
(paper's own bound), `ρφ_0(ρ)` attains a finite maximum at some `ρ*=O(1)`,
whose physical location `r_m(t) = ρ*√(2(T−t)) → 0` — i.e. **exponent
`b=1/2`** in this programme's `(G1)/(G2)` notation, matching the subagent's
claim. At that maximizer the correction term is `O(T−t)→0`, so **`M(t) =
max Γ → ρ*φ_0(ρ*)`, a nonzero constant, to leading order as `t→T` — exponent
`a=0`.**

**Caveat the orchestrator adds to the subagent's framing**: the subagent's
report describes this as "reduces identically... with no residual
t-dependence," which overstates it slightly — the `2α(T−t)ρ²` term is a
genuine residual, just one that is asymptotically subleading, not absent.
The substance (`a=0` in the asymptotic-rate sense `cor:cd1-rate` actually
uses) is correct and independently confirmed by the derivation above.

Because `a=0`, `cor:cd1-rate`'s floor is only the trivial `K≥0` here — this
construction satisfies the theorem's hypotheses but sits in a degenerate
corner the rate corollary isn't built to say anything sharp about. Precise
statement of where this leaves Deliverable 2: **the class of constructions
satisfying `thm:cd1-general-floor`'s hypotheses now has (at least) two
confirmed members** (OpenAI's leading axisymmetric field; Zhang's ASNS
construction), but **the class for which `cor:cd1-rate` yields a nontrivial
numeric floor (`a>0`) still has exactly one known member** (OpenAI). This is
a sharper, more accurate version of H1 than a flat "one member," not a full
S55-style three-way comparison (there is nothing nontrivial to compare for
Zhang's construction).

Part 2 of Zhang's theorem (log-transformed variant, bounded `Γ` while raw
`v_θ` blows up) was reported by the subagent with lower confidence
("moderate confidence... please re-derive by hand") and was **not**
independently re-derived by the orchestrator in this session — flagged here
as unverified detail, not load-bearing for the session's conclusions (both
outcomes give `a=0`).

## Verification performed by the orchestrator (independent of the subagent)

- Fetched arXiv:2309.08495 abstract directly: confirmed non-axisymmetric
  focus, quoted verbatim, matches subagent's claim exactly.
- Fetched arXiv:2407.06776 abstract and author list directly (Córdoba,
  Martínez-Zoroa, Zheng — confirmed); downloaded the full PDF and extracted
  text via `pdftotext`; located and quoted §1.2.5 verbatim — matches the
  subagent's claim exactly, independent of the subagent's own reading.
- Fetched arXiv:2311.12306 abstract and full PDF text directly; located
  Theorem 1.1 and equations (1.1)–(1.24) verbatim; independently re-derived
  `Γ(r,t) = ρφ_0(ρ) + 2α(T−t)ρ²` from the paper's own formulas (not copied
  from the subagent's derivation) and confirmed the `a=0`, `b=1/2` claims,
  with the stated caveat about "no residual t-dependence" being an
  overstatement of an otherwise-correct asymptotic claim.
- Confirmed via `PROGRESS.md`'s S60 KEY FINDING that S60 used
  arXiv:2407.06776 only for its Thread-2 hypodissipative-ceiling comparison,
  never for Thread-1 axisymmetry — so the subagent's "correction" is
  correctly characterized as completing an unchecked gap, not fixing an
  error.
- Did not independently re-verify the Córdoba–Córdoba PNAS citation's exact
  bibliographic details beyond confirming the inequality's standard form is
  consistent with the programme's existing knowledge of that literature; did
  not independently re-derive the fractional-operator counterexample from
  scratch before reading the subagent's version, but did re-derive it
  step-by-step by hand afterward and found no gap.
- Did not re-verify Zhang's Part 2 (log-transformed) computation — flagged
  above as unverified, non-load-bearing detail.

No arithmetic or attribution error was found in S60's own prior derivations
on the parts re-examined this session.

## Recommendation

**Close the CD-2 line here.** Deliverable 1 is a complete, precisely-located
negative result. Deliverable 2's search was genuine: both named candidates
were independently verified and correctly ruled out on axisymmetry grounds
(not existence), and one additional, unnamed instance was found, checked,
and shown to sit in the theorem's degenerate (`a=0`) corner rather than
extending its reach. Nothing here produces a new number worth its own
`proofs/` theorem. No proof-document edit was made this session — none was
warranted at the bar this programme holds such edits to; a single-sentence
future remark noting the second (degenerate) class member in
`sec:cd1-general-floor` is a bounded, low-effort option if ever wanted, not
a reason to reopen this line now. No further (C)/(D) forcing-floor session
is recommended.

## Standing discipline confirmed

- Session tag: `git log --oneline -1` at S62 start showed tip `08ae4b9`
  ([S61] plan commit); `ls SESSION-LOG/` showed S60 as the highest existing
  log file; S61 has a handover but no log (by design — it was a
  planning-only commit); **S62 is correctly the next tag.**
- No forcing introduced into this programme's own arguments (H5) —
  confirmed on inspection of everything written above.
- `book/`, `layer3/`, `layer4/`, `results/` untouched.
- No worktree was created (no proof-document edit made); nothing to merge.
