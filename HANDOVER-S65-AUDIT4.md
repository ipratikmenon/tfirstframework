# Handover — fourth strategic audit: two more untried pathways toward (A)/(B)

Written 2026-10-07. Owner's direction, verbatim: "Lets have an another round
of audit on some more possible pathways?" — a fourth strategic audit,
requested after three consecutive negatives (S56: six angles, S59: three
angles, S64: three angles). **Pre-assigned tag: S65.** Verify with
`git log --oneline -1` (expect tip `b163b4d`, "[S63] docs: book chapter...")
and `ls SESSION-LOG/` (highest existing tag should be S64) before writing
anything.

## Be honest about where this stands before starting

This is the fourth audit in a row. Twelve angles have been tried and all
twelve have returned negative, for reasons ranging from "reduces to Wall 1
in different vocabulary" to "addresses a different axis than (A)/(B) asks
about" to "is itself the unresolved problem restated." **A fourth negative
here is a live, fully legitimate, expected outcome — do not manufacture a
positive finding to avoid a fourth negative in a row, and do not pad a thin
report to make the session look more productive than it was.** If both
angles below reduce to Wall 1 or to "the unresolved problem restated" like
everything before them, say so plainly and say it fast — this session does
not need to be long to be complete.

## Read first, so this audit does not repeat any of the twelve angles already tried

Read, in full: `HANDOVER-S56-AUDIT-FINDINGS.md`, `HANDOVER-S58-AUDIT2.md`,
`HANDOVER-S59-AUDIT2-FINDINGS.md`, `HANDOVER-S64-AUDIT3.md`, and the S56/
S57/S59/S64 KEY FINDING blocks in `PROGRESS.md`. The twelve angles already
closed: direct BKM/CZ-logarithm weakening; recent ancient-Euler Liouville
theorems beyond Seregin's four papers; computer-assisted/interval-arithmetic
profile exclusion; Tao/Barker quantitative regularity criteria;
energy-flux/cascade arguments; helicity-based mechanisms; Grujić's
log-weighted BMO mechanism (priced separately at S57); the programme's own
founding thermodynamic thesis (Route 1/Route 2); stochastic Lagrangian
representations (Constantin–Iyer) and rough path theory/regularity
structures; kinetic/hydrodynamic (Boltzmann-to-NS) limits; random-data/
probabilistic well-posedness (Nahmod–Pavlović–Staffilani, Deng–Nahmod–Yue
lineage, Földes–Sy); anisotropic Littlewood-Paley regularity criteria
beyond BKM (Cao–Titi and continuations). **Do not re-run any of these.**

**Only revisit a closed angle if you find something specifically dated
after Oct 7 2026 minus... no — specifically dated after S64 (16 Sept 2026)
that materially changes the picture.** A three-week gap has now passed
since the last recency check, wider than any prior pass — treat the
recency angle below as genuinely needing more than a token search this
time.

## The two new angles

### 1. Escauriaza–Seregin–Šverák's L³ criterion and the backward-uniqueness/unique-continuation machinery behind it

The ESS theorem (Escauriaza, Seregin, Šverák, 2003: *"L_{3,\infty}-solutions
of Navier-Stokes equations and backward uniqueness"*) proves that a
Leray-Hopf weak solution blowing up at time `T` must have
`limsup_{t\to T} \|u(t)\|_{L^3}=\infty` — strengthening the critical space
beyond the Ladyzhenskaya-Prodi-Serrin scale. Crucially, **the proof
technique is not a direct Calderón-Zygmund/Biot-Savart estimate on `∇u` at
all** — it proceeds via blow-up compactness, reduction to a bounded ancient
solution, and a **backward uniqueness theorem for the heat/Stokes operator**
(built on Carleman estimates), a genuinely different piece of machinery
from anything this programme's internal toolkit or any of the twelve closed
angles has used directly.

Investigate, from primary sources (the original ESS paper and any
2020-2026 continuations or simplifications of the backward-uniqueness
technique):

- Does the backward-uniqueness/unique-continuation route supply, anywhere
  in its proof, an estimate that could serve as a genuinely different
  substitute for the specific `L^\infty`-type bound on `∇u` recovered from
  `ω` via Biot-Savart that Wall 1 (`lem:local-enstrophy-typeI`) obstructs —
  or does the blow-up compactness step itself secretly require exactly that
  kind of control to extract a nontrivial ancient-solution limit (the same
  way S53's Type-II reconnaissance found rescaling toward a supercritical
  rate forces an Euler limit where this programme's tools are void)?
- Is there a *quantitative* form of ESS (explicit rate of `L^3`-norm
  blow-up, not just qualitative divergence) in the literature, and if so,
  does its proof technique reveal anything about whether backward-uniqueness
  methods are fundamentally compatible with removing the logarithm, or
  whether they simply relocate the same difficulty into the compactness
  step?
- State plainly whether this is a different toolkit hitting the same wall
  in different vocabulary (the S57/S64-angle-2 pattern), a toolkit that
  cannot reach (A)/(B) at all for a structural reason (the S59/S64-angle-1
  pattern), or something that has not been tried in a way this programme
  can evaluate.

### 2. Concentration-compactness / the "critical element" (minimal blow-up solution) method

Originating in critical dispersive PDE (Kenig–Merle for NLS/wave) and
adapted to Navier–Stokes by **Kenig–Koch** (*"An alternative approach to
regularity for the Navier-Stokes equations in critical spaces,"* 2011) and
**Gallagher–Koch–Planchon** (profile decomposition for NS in critical
Besov spaces, several papers 2011-2016). The method: assume for
contradiction that (A)/(B) is false, extract a "minimal" blow-up solution
via profile decomposition and concentration-compactness in a critical
space, derive strong compactness/rigidity properties of this minimal
element, and seek a contradiction. This is a genuinely different technical
toolkit (functional-analytic compactness methods, not direct PDE estimates
on `∇u`) that this programme has never engaged.

Investigate, from primary sources:

- What do Kenig–Koch and Gallagher–Koch–Planchon actually prove — do they
  reach unconditional global regularity, or only a reduction (e.g. "if
  blow-up occurs, a minimal blow-up solution with property X exists," with
  X itself an open question)? Read for the actual theorem statements, not
  the method's reputation.
- At the point where the minimal-element argument needs to rule out the
  extracted critical element, does it require exactly the same kind of
  `L^\infty(\nabla u)` or Biot-Savart-recovery control this programme's
  Wall 1 obstructs — i.e., does concentration-compactness relocate the
  difficulty into "ruling out the minimal element" the same way profile
  rigidity did for this programme's own blow-up-profile work at S39 (a
  structural parallel worth checking explicitly, since this programme has
  direct experience with exactly this kind of reduction)?
- Any 2023-2026 continuations of this specific method for 3D NS
  specifically (not just well-posedness in critical spaces generally,
  which is a separate and larger literature already adjacent to what S59
  covered via Koch-Tataru-type small-data results).
- State plainly whether this reduces to Wall 1, reduces to a different
  open rigidity question analogous to this programme's own, or offers
  something genuinely new.

### 3. Recency rescan — wider window this time, treat it as a real pass, not a token check

Three weeks have passed since S64's recency check (16 Sept 2026 to 7 Oct
2026) — wider than any prior gap in this programme's history between audits.
Check specifically, with dated search:

- Any claimed weakening or removal of the BKM/CZ logarithm itself.
- Any new Grujić-lineage, Tao/Barker-lineage, or Cao-Titi-lineage
  (one-component/anisotropic criteria) publication.
- Any update to the OpenAI/Clay situation (the `navier-stokes.org` tracker
  or equivalent) — it has now been roughly a month since the Sept 8 2026
  announcement and Sept 11 2026 Clay response; check whether anything has
  materially progressed (official verification, a retraction, a new
  construction) that bears on this programme's record, bearing in mind
  this is about (C)/(D) and does not itself change Wall 1's status.
- Any new preprint specifically on ESS-type L³ criteria or
  concentration-compactness for 3D NS, dated in this window, that might be
  directly relevant to angles 1-2 above.

## Deliverable

Same shape as S56/S59/S64: for each of the two new angles, a clear finding
(something real and checkable found, or a clean, well-justified "no"), with
primary sources fetched and quoted, not recalled. The recency rescan gets a
short, dated report. End with one of:

- A concrete, scoped next work package, stated with named failure modes in
  advance — **only if a genuine candidate survives scrutiny**, not
  manufactured to avoid a fourth negative in a row.
- An honest statement that both angles reduce to Wall 1 (or to a structural
  mismatch with what (A)/(B) actually asks, or to a differently-named
  version of the same open problem) in different vocabulary, consistent
  with the pattern of the prior three audits, and that the recency rescan
  found nothing new.

**Do not manufacture a positive recommendation.** This programme's most
valuable sessions have been precise negative results, and this audit is
explicitly requested with that understanding.

## Standing discipline

- **Session-tag collision has happened six times historically.** Pre-assign
  S65 here; verify with `git log --oneline -1` and `ls SESSION-LOG/` before
  writing anything, and re-check at merge time since this is a solo
  (non-parallel) session this round.
- Fetch every source directly. Do not recall S56/S59/S64's summaries of
  what they found — re-read those handovers and session logs directly, and
  fetch any new external source fresh.
- This is a **literature-and-strategy session**. Do not write proof-document
  content, numerics, or touch `book/`, `layer3/`, `layer4/`, `results/`.
- **This programme still never introduces forcing into its own regularity
  arguments**, and this audit is about (A)/(B), not (C)/(D) — the (C)/(D)
  forcing-floor line (S60-S62) is closed and not reopened by anything here,
  including the OpenAI/Clay recency check (report status only, do not
  re-engage the pricing exercise).
- Do not commit or push. Report findings and recommendation directly.
