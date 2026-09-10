# HANDOVER — S46 → Opus Execution Programme

**From:** Fable (strategic audit, session S46, 2026-09-10)
**To:** Opus (execution, sessions S47+)
**Status of this document:** the authoritative work plan. `PROGRESS.md` points here.
Read `PROGRESS.md` first (CLAUDE.md protocol), then this file end to end before
touching anything.

---

## 0. The decision this plan serves

After the S46 audit the programme owner decided:

1. **Every route tried and every wall hit gets documented** — precisely, with
   line references, distinguishing *proved impossible* from *open*.
2. **The goal does not change.** We stay with the exact Prize equations
   (unforced incompressible NS, `f ≡ 0`, ℝ³ and 𝕋³) and their criteria, and we
   keep trying to prove global smoothness. No forcing terms, no modified
   equations, no "weak internal forces" — those were considered and rejected
   (S40, S46).
3. **The main line pivots** from the general alignment pincer to the
   **axisymmetric-with-swirl** subproblem, with helicity/twist as the new
   structural input, because that is where (a) Type-I is already excluded by
   the literature, (b) the Calderón–Zygmund logarithm that has blocked every
   route since S34 is expected to weaken or vanish, and (c) the programme's
   numerical intuitions and infrastructure already live.

The audit's findings that force this are in §1. Do not re-litigate them; do
extend them if you find them wrong, in the S31/S44 corrigendum tradition.

---

## 1. What the audit established (verified against the document)

Each item below was checked against `proofs/claim_a_3d_proof_attempt.tex`
(9787 lines) in S46. Labels and line numbers are given so you can re-verify.

**W1 — The logarithm is structural, not bookkeeping.** `lem:local-enstrophy-typeI`
(line ~5664) introduces `𝓛 = log(e + ‖u‖_{H³}/M)` through
`‖∇u‖_∞ ≲ M·𝓛`. That is the L^∞ endpoint failure of the Biot–Savart operator
(a Calderón–Zygmund operator, unbounded on L^∞). Every enstrophy budget in the
document inherits it. It is the *same* log that blocks crude K-absorption
(`prop:crude-K`), the annulus core (`rem:s41-shortfall`) and the corrected band
route (`prop:s44-linfty-route`). Any estimate routed through a sup-norm of a CZ
image pays it. This architecture routes everything through one.

**W2 — The profile route does not evade W1; it relocates it.**
`prop:transfer-audit(d)` ("profiles are logarithm-free") is
*Proved-modulo-(R)*, and its proof says the log disappears because "the CZ step
is replaced by the interior estimate (R)". Hypothesis **(R)** — scale-uniform
interior regularity of the rescaled Type-I family — *is* an assumed log-free
gradient bound `‖∇U‖_∞ ≤ C|s|⁻¹`. Status table: "Conditional (Serrin-type)".
An earlier S46 draft claimed the profile branch structurally escaped the
barrier; that claim was withdrawn on reading the proof. **Whether (R) follows
from existing Type-I ε-regularity is the single cheapest decisive question in
the programme** → WP2.

**W3 — The ceiling of the entire S32–S45 line is below the Prize.**
`rem:s39-distance`: if *every* open item closed — (A-up),
`target:disorder-depletion`, (R), (N), (N′) — the conclusion is "no Type-I
singularity", which the document itself calls "strictly less than the Clay
problem". Every estimate uses Type-I rates *essentially*; a Type-II singularity
"would rescale to nothing". The general pincer is blind to Type-II.

**W4 — The level iteration is proved impossible** (S41): `prop:s41-nocontraction`
(Λ_k ≥ 2⁸ for every Young parameter), `prop:s41-vacuous`, `prop:s41-coarea`
(no level selection helps, by a co-area identity).

**W5 — `conj:K-refined` and `target:band-absorption` are independent** (S44):
counterexamples both ways (`prop:s44-a-fails`, `prop:s44-b-fails`); structural
cause is disjoint supports (`rem:s44-disjoint`). A common weaker target
`conj:lc-s44` exists; corrected, it needs the same decaying constant `≤ C/𝓛`.
A first-draft strength claim was inverted and withdrawn (`rem:s44-corrigendum`).

**W6 — The c_K numerics were never discriminating.** S38 proved σ\*-decay is
*necessary for blowup*. A regular solver never blows up, hence never shows the
decay, hence measures `c_K ≈ 1` — which both hypotheses predict. Four sessions
(S34, S35, S42–43, S45) confirmed the non-information. S45 measured the
𝓛-slope ≈ 0 (not −1) across 𝓛 ∈ [1.35, 4.68]. **Do not measure c_K again to
test the conjecture.** The suite remains valid for *mechanism* questions (WP4).

**W7 — Verification substrate.** Six sessions (S31, S34, S35, S39, S41, S44)
had consequential *mathematical* errors caught only by review, on a document
that has never been compiled. The discipline that caught S41 and S44 —
isolated worktree, then independent re-derivation of load-bearing steps before
merge — is mandatory below.

---

## 2. Global rules for every work package

- **Isolation:** every edit to `proofs/claim_a_3d_proof_attempt.tex` happens in
  a git worktree. Nothing merges until the orchestrator has independently
  re-derived every load-bearing computation. A self-report is not verification.
- **Labelling:** every mathematical statement carries one of *Proved*,
  *Proved mod (X)*, *Open*, *Proved impossible*, *Conditional on …*. No bare
  "we have shown".
- **Literature first, always.** Before any new analytical claim, situate it
  against the named prior work in that WP. The programme has a history of
  rediscovering known results and of over-claiming marginal targets
  (`rem:s39-distance` item 3). Rediscovery is fine if labelled; unlabelled it
  is an error.
- **Numerics:** tests first; calibrate every new instrument on planted
  synthetic cases before real data (S45 practice); time-split every result and
  discard the early transient (S43 lesson: every spike lived at t<3); state the
  cell count of every region measured (S45: TG's r\*-ball had 81 cells);
  never let a core span < 4 grid cells (S40 floor lesson). House rules from
  CLAUDE.md hold: spectral, 2/3-dealiased, RK4, adaptive CFL, all properties
  via layer1, every result has `claim_id + key_metric + verdict` in
  `results/results.db`.
- **Logging:** per CLAUDE.md — update `PROGRESS.md` and write
  `SESSION-LOG/YYYY-MM-DD-SNN-*.md` before every commit; one commit per
  completed WP milestone; never commit with failing tests.
- **Cost discipline:** the owner hit a spend limit on 2026-09-10. Prefer one
  validated run to three speculative ones. Prefer N=64 with a stated
  resolution caveat over N=128 by default; go to N=128 only to settle a
  specific resolution question.

---

## 3. Work packages, in dependency order

### WP0 — Verification substrate  *(do first; cheap; gates everything)*

**Goal.** A compiling proof document and a mechanical check script, so the
S39/S44 class of syntactic error cannot recur silently.

**Deliverables.**
1. Install a TeX toolchain (`apt-get install texlive-latex-base
   texlive-latex-extra texlive-science` or minimal equivalent — `apt-get`
   works in this container; poppler-utils was installed this way in S40).
   Compile `proofs/claim_a_3d_proof_attempt.tex`. If full compilation is
   impossible, enumerate every error with line numbers — that list is itself
   the deliverable.
2. Fix the **known pre-existing defects**: undeclared environments
   `observation` and `openproblem` (lines ~5039, ~5117, ~5144 — declare them
   via `\newtheorem` in the preamble, lines ~23–34); four unresolved refs
   `eq:NS-vorticity`, `lem:GN-ehat`, `prop:blowup-alignment`,
   `thm:vorticity-uniform` (lines ~4740–5040 — locate the intended targets or
   mark them explicitly as dangling with a red note; do not invent labels).
3. `scripts/check_proof.py`: brace balance (comment-stripped, `\{ \} \\`
   handled), `\begin`/`\end` multiset equality, undeclared environments,
   duplicate labels, unresolved `\ref`/`\eqref`/`\cite`, and the malformed
   `end{[a-z]*>` typo class found in S44. Exit non-zero on any failure. Wire
   it as a pre-commit hook or document its invocation in CLAUDE.md's Test
   Commands.

**Acceptance.** `check_proof.py` passes; compilation succeeds or its failures
are enumerated; no mathematical content changed.

**Named failure mode.** "Fixing" an unresolved ref by guessing a plausible
label. If the target cannot be identified with certainty, leave a visible
dangling-reference note.

---

### WP1 — The Walls Document

**Goal.** One authoritative, dry, referenced record of every route and its
obstruction, so no future session repeats a proved-impossible step and so the
material exists for a standalone survey.

**Deliverable.** New section `\section{Routes and Walls: A Referenced Record
(Session S47)}`, `\label{sec:walls-s47}`, appended before the bibliography.
One subsection per wall, **fixed template**:
*(a) the route and what it hoped to prove; (b) the exact obstruction with
label/line refs; (c) classification — PROVED IMPOSSIBLE / OPEN / SUPERSEDED /
WITHDRAWN; (d) what it forbids or permits for future work.*

Walls to cover, minimum: S31 audit (`§sec:audit-s31`, all seven grounds); W1
the CZ log with its three appearances; W2 the (R) relocation; W3 the Type-II
ceiling; W4 level iteration; W5 c_K/band independence + the S44 corrigendum;
W6 non-discriminating numerics (give the *logical* argument via S38, not just
the data); the S35 retracted "bonus term"; the S39 `prop:localized-cf`
non-transfer; the S40 finding that OpenAI's forced construction cannot be
de-forced (so forcing-based intuitions are out of scope). Close with a status
table.

**Acceptance.** Every claim in the section resolves to a label or line in the
document; classification column complete; `check_proof.py` passes.

**Named failure mode.** Narrative drift — the section becoming a defence of
the programme rather than a record. Every sentence should survive being read
by the S31 auditor.

---

### WP2 — Settle (R) against the literature  *(decisive gate; one session)*

**Goal.** Determine whether hypothesis (R) — scale-uniform interior regularity
of the rescaled Type-I family, `‖∇U‖_∞ ≤ C(c_I)|s|⁻¹` — is already a theorem.
This decides whether the profile route is genuinely log-free (W2).

**Deliverable.** New subsection in `§sec:rigidity-s39` (or a short new
section) resolving (R) as exactly one of:
- **(i) (R) follows from [named theorem].** Then `prop:transfer-audit(d)`
  becomes *Proved*, the profile route is log-free, and the aligned horn
  reduces to (A-up) alone. State the exact hypotheses under which the cited
  theorem applies and confirm the rescaled family satisfies them.
- **(ii) (R) is equivalent to / a strict part of a named open problem.** Say
  which, and how close.
- **(iii) (R) is false in general.** Exhibit why.

**Literature to check, with exact statements (use WebSearch if needed):**
Escauriaza–Seregin–Šverák 2003 (L³,∞ regularity, Type-I compactness);
Seregin 2012 (Type-I-improved ε-regularity; already cited for (N));
Koch–Nadirashvili–Seregin–Šverák 2009 (Liouville, already cited);
Nečas–Růžička–Šverák 1996 and Tsai 1998 (self-similar exclusion);
Caffarelli–Kohn–Nirenberg 1982 (ε-regularity).

**Named failure modes — these are the traps:**
1. **Smallness vs boundedness.** CKN-type ε-regularity needs the *scaled local
   energy to be small*. Type-I gives it *bounded by c_I*, not small. Deriving
   (R) from ε-regularity without a smallness argument is the exact error to
   avoid. If smallness is available only for `c_I` small, say so — that is
   outcome (ii), not (i).
2. **Circularity.** (R) is close to "Type-I blowup has smooth limiting
   profiles", which is a large part of what the aligned horn tries to prove.
   Check that no cited theorem secretly assumes the conclusion.
3. **Interior vs global.** (R) is *interior* regularity at fixed rescaled
   distance from the singular point — weaker than regularity at the point.
   Do not conflate the two in either direction.

**Acceptance.** A labelled verdict with the applied theorem's hypotheses
checked line by line, or a precise statement of the gap.

**Branch after WP2.** Outcome (i): WP3 also reopens the aligned horn on
profiles (target: (A-up)). Outcome (ii)/(iii): WP3 proceeds purely
axisymmetric; the general profile route is parked, and W2 is recorded as
closed.

---

### WP3 — Axisymmetric-with-swirl: the analytical pivot

**Goal.** Rebuild the programme's direction-field machinery in the
axisymmetric class, situate it against the known axisymmetric literature, and
aim at a *swirl-only* regularity criterion. This is the **Prize equations
exactly**, restricted to a symmetry class that is **still open**.

**Setting (verify each line before using).** Cylindrical `(r, θ, z)`,
`u = u_r e_r + u_θ e_θ + u_z e_z`, `∂_θ ≡ 0`, swirl `u_θ ≢ 0`. Vorticity:
`ω_r = −∂_z u_θ`, `ω_θ = ∂_z u_r − ∂_r u_z`, `ω_z = (1/r)∂_r(r u_θ)`.
Swirl (angular momentum) `Γ := r u_θ` satisfies
`∂_t Γ + u_r ∂_r Γ + u_z ∂_z Γ = ν(Δ − (2/r)∂_r)Γ`
and hence a **maximum principle**: `‖Γ(t)‖_∞ ≤ ‖Γ(0)‖_∞`. The quantity
`ω_θ/r` satisfies
`∂_t(ω_θ/r) + u·∇(ω_θ/r) = ν(Δ + (2/r)∂_r)(ω_θ/r) + ∂_z(Γ²)/r⁴`.
**The entire difficulty of the swirl case is the single source term
`∂_z(Γ²)/r⁴`.** Without swirl it vanishes, `ω_θ/r` is materially conserved up
to diffusion, and regularity follows (Ladyzhenskaya 1968; Ukhovskii–Yudovich
1968).

**Known results the work must be situated against (check exact statements):**
- No swirl ⟹ global regularity: Ladyzhenskaya 1968; Ukhovskii–Yudovich 1968.
- **Type-I excluded with swirl:** Chen–Strain–Tsai–Yau 2008; Koch–Nadirashvili–
  Seregin–Šverák 2009 (Liouville, already `KNSS2009` in the bibliography);
  Seregin–Šverák. *This is not our result; cite it and use it.*
- Swirl-based regularity criteria under smallness/decay of `Γ` near the axis:
  Lei–Zhang 2017 ("Criticality of the axially symmetric Navier–Stokes
  equations"), Chen–Fang–Zhang 2017, Wei 2016; earlier Neustupa–Pokorný,
  Chae–Lee. **Read their exact hypotheses before claiming anything new.**
- One-component criteria in general: Kukavica–Ziane; Cao–Titi.

**Deliverables.**
1. `§sec:axisym-s48` (label to match session): the direction field `ê`, `w`,
   `σ\*`, the magnitude/direction equations and the Gram-negativity lemma
   re-derived in axisymmetric variables. **Identify exactly what happens to
   the CZ log (W1)**: in this class ∇u is recovered from `(ω_θ, Γ)` by
   2D-type elliptic problems in `(r, z)`; determine whether the L^∞ Biot–Savart
   bound still costs `𝓛`, costs a weaker factor, or is log-free for the
   poloidal part. This is the central question of the pivot — answer it
   first, honestly, before building anything on it.
2. State what H1 (Bridge) and H2 (localised CF) *become* in this class, and
   which parts are already discharged by the Type-I literature.
3. **Target:** a regularity criterion in terms of `Γ` alone (or `Γ` plus one
   scale-invariant quantity), stated with hypotheses strictly compared to
   Lei–Zhang / Chen–Fang–Zhang. A criterion that is *new* only in its
   constant is not a result; a criterion that replaces a smallness condition
   by a boundedness or structural condition is.
4. Ledger: what remains for *Type-II axisymmetric* — the true residual open
   problem of this class.

**Named failure modes.**
1. Rediscovering Lei–Zhang and calling it new. Situate first.
2. Reporting Type-I exclusion as own work.
3. Accidentally reducing to the no-swirl case (e.g. by an assumption that
   forces `Γ ≡ 0`).
4. Claiming the log is gone without a written estimate showing where the
   L^∞ bound on ∇u comes from in `(r, z)`.
5. Marginal targets. If a criterion needs a constant to beat exactly 1, flag
   it as marginal (`rem:s39-distance` item 3) and do not build on it until
   the constant is nailed.

---

### WP4 — Numerics redirect  *(independent of WP3; run in parallel)*

**Goal.** Point the existing, tested suite at *mechanism* questions in the
axisymmetric/helicity class — questions whose answers carry information in a
regular flow (unlike c_K, W6).

**4a. First pass — no new solver.** Use the existing 3D periodic Route-2 JAX
solver with axisymmetric ICs and add **axisymmetric diagnostics**:
`Γ = r u_θ` (check its max principle numerically — a free correctness test),
`ω_θ/r`, and the source term `∂_z(Γ²)/r⁴` localised in `(r, z)`. Compare
against the exact Burgers vortex as a calibration (it is an exact steady
solution; the solver should hold it).

**4b. `layer3/twisted_ring_IC.py` — the helicity structure.** A twisted vortex
ring parameterised by rotational transform `ι` (the stellarator parameter:
number of poloidal twists per toroidal circuit). Diagnostics in a new
`layer4/helicity_diagnostics.py`: helicity `H = ∫ u·ω` and local density
`h = u·ω`; twist density along the core; stretching rate `α = ê·S·ê`; and the
**pointwise correlation between twist density and α on `{|ω| ≥ M/4}`**.
*The question: does twist locally deplete stretching?* That is a lemma
candidate, not a blowup proxy — which is why this measurement is
information-bearing where c_K was not. Sweep `ι`; time-split; report the
correlation's sign and its stability across `ι`.

**4c. Reconnection detection.** Flag topology change (helicity jump,
core-distance minimum) so 4b can distinguish pre- and post-reconnection
regimes.

**4d. Later, gated on 4a results:** a true `(r, z)` axisymmetric solver
(Chebyshev in `r` is still spectral and so house-rule-compliant; state the
amendment explicitly in CLAUDE.md if adopted), enabling Hou-style
computer-assisted profile search near the axis. **Do not attempt in the first
pass.**

**Named failure modes.**
1. A twisted ring with *zero net helicity by symmetry* (e.g. a mirror-symmetric
   twist). Check `H ≠ 0` for `ι ≠ 0` in the tests before any run.
2. Reading the twist–stretching correlation in the early transient (S43).
3. Ring core under 4 cells (S40).
4. A bespoke finite-difference solver (house-rule violation).
5. Reporting a correlation as a mechanism without the sign being stable across
   `ι` and across at least two resolutions.

---

## 4. Sequencing and parallelism

```
WP0 ──► WP1 ──► WP2 (decisive gate)
                  │
        ┌─────────┴──────────┐
        ▼                    ▼
   WP3 (analytical,     WP4 (numerical,
   proof doc, worktree)  layer3/layer4)
        │                    │
        └──────┬─────────────┘
               ▼
   WP4d + Hou-style search (gated)
```

WP0→WP1→WP2 are sequential and each is roughly one session. WP3 and WP4 are
independent and may run as parallel worktrees/agents. WP2's outcome branches
WP3 (see WP2).

---

## 5. What is explicitly out of scope

- Forcing terms, modified equations, "weak internal forces", magnetism
  analogies as *mechanisms* (geometries borrowed from them are fine — S40, S46).
- Measuring `c_K`, `c_K^band`, or the band ratios again to test the
  conjectures (W6). They may be computed as *context* in WP4 but not as
  evidence.
- Claiming anything about Type-II from Type-I machinery (W3).
- Reopening the general (non-axisymmetric) level iteration (W4).
- The string-theory / fluid–gravity reframing as a proof method: it is a
  long-wavelength effective correspondence, UV-blind at the singular scale;
  recorded in S46 as a reformulation, not a route.

---

## 6. Definition of done for this handover cycle

- WP0–WP2 complete and logged; WP2 verdict stated with its label.
- WP3 deliverable 1 answered (what happens to the log in axisymmetric flow)
  with a written estimate — this single answer determines whether the pivot
  is real.
- WP4a–4c built, calibrated on Burgers + planted synthetic twist, one `ι`
  sweep run and time-split.
- `PROGRESS.md` and session logs current; `check_proof.py` green; all tests
  green; everything pushed to `claude/shared-conversation-jgoavo`.
- A one-paragraph honest statement in `PROGRESS.md` of whether the
  axisymmetric pivot has a live path to a new criterion, or has hit its own
  wall — in which case that wall goes into `§sec:walls-s47` like the others.
