# Handover — pivoting the active line to (C)/(D): how much forcing, and is it removable?

Written 2026-09-13, immediately after S59. Owner's direction: pivot the
programme's active research question from (A)/(B) to (C)/(D).

## What changes, and what does not

**What changes.** The active question is no longer "prove global regularity
for the unforced Prize equations." It is now: **for the constructions that
claim or approach forced finite-time singularity formation, how much
forcing is actually required, and can that requirement be reduced or
removed?** This is explicitly a step toward understanding (C)/(D), not
(A)/(B) — S56 §3.3 and S59's recommendation both named this as the honest
alternative once the cheapest moves against (A)/(B) were exhausted (internal
toolkit: S52, S53; best external candidate: S57; the programme's own
founding thesis: S59).

**What does not change.** This programme still never introduces a forcing
term into its *own* attempts at a regularity proof — that rule (Wall 9,
CLAUDE.md's standing prohibition) was about not smuggling forcing into an
(A)/(B) argument. It says nothing about *studying* forced constructions that
already exist in the literature, which is what this new line does. Do not
confuse "we study forcing requirements" with "we force our own equations" —
the former is now the programme's job; the latter remains forbidden, always.

## Background already in hand — read before anything else

- **S40 / Wall 9**: OpenAI's forced NS construction cannot be de-forced by
  deletion; the forcing is defined as the construction's own residual.
- **S54 / S55**: this programme already built the only quantitative
  de-forcing tool that exists anywhere in this line — `lem:wpt3-forced-max`,
  a rigorous forced generalisation of the axisymmetric swirl maximum
  principle, giving a genuine lower bound `|f_θ| ≳ τ^{-(A+1)}` on the forcing
  the OpenAI profile's swirl growth requires, confirmed against the primary
  source's own stress-exponent scaling by an independent route.
  **This tool currently applies to one profile. Generalising it is
  Deliverable 1 below.**
- **S53 / S56**: OpenAI has a *separate*, **unforced** Euler paper (Theorem
  1.1, `u_0\in C^\infty_{c,\sigma}(\mathbb R^3)`, verified from the primary
  PDF at S56) — the first object in this programme's history whose geometry
  is not behind a de-forcing prohibition, because it was never forced in the
  first place.
- **S56 §3**: Tao's own account of the Alpöge–Buckmaster programme —
  removing the forcing from their porous-medium/Boussinesq/Euler
  constructions "looks very feasible... in the near future," and the method
  "has a high likelihood of also extending to Navier-Stokes as well." This
  is the field's own stated direction of travel.

## The two concrete threads — treat this as a reconnaissance session, in the tradition of S53 (WP-T1)

**Do not jump straight to heavy technical work.** This is genuinely new
territory for this programme (general forcing-floor theorems; the
vanishing-viscosity literature). Spend real effort mapping what is
tractable before committing to a specific proof target, the same discipline
that made S53 productive before S54/S55 committed to a specific calculation.

### Thread 1 — generalise the forcing floor from one profile to a class

`lem:wpt3-forced-max` currently answers: *given* a specific swirl growth
rate and core geometry (OpenAI's), how much forcing is needed. The natural
generalisation: **for the whole class of axisymmetric constructions with a
prescribed (Type-II) swirl growth rate on a prescribed core-radius
geometry, is there a general lower bound on the required forcing, stated
once and for all rather than per-profile?**

- Re-derive `lem:wpt3-forced-max`'s proof (Danskin's envelope theorem
  applied to `M(t) = \max\Gamma`) and identify exactly which hypotheses were
  specific to the OpenAI profile and which are generic to any axisymmetric
  field with a given growth rate and radius. Write the general statement
  explicitly, with its own hypotheses, before checking it against any
  specific example.
- If a genuine general theorem results, this is real, new, general
  mathematics — a "cost of forcing" lower bound — of a kind this programme
  has not produced before (previous results were checks against one
  external object). State plainly if the generalisation is not clean (e.g.
  if the OpenAI-specific self-similar exactness that made S55's floor exact
  rather than merely asymptotic does not generalise) — that is itself
  useful information about how special that one profile is.

### Thread 2 — the OpenAI Euler paper and the vanishing-viscosity question

Their construction is Euler, unforced, from compact smooth data. The
natural, well-precedented mathematical question this raises, entirely
independent of anything this programme built: **does the singular behaviour
survive as viscosity is added and taken to zero, or does it fail to survive
any positive viscosity?**

This connects to the mature vanishing-viscosity and anomalous-dissipation
literature (Onsager's conjecture and its resolution, Constantin–E–Titi-type
work, and more recent work on inviscid limits for singular Euler solutions).
**This is a literature survey first, not a computation.** Establish:

- What, if anything, is known in general about whether Euler singularities
  of this type (self-similar-adjacent, from smooth compact data via
  amplification of localized oscillation — read the actual mechanism, not
  just the theorem statement, from the primary source used at S56) survive
  or are destroyed by vanishing viscosity. Fetch primary sources; do not
  recall.
- Whether the specific mechanism in this construction (iterated
  amplification across scales via wave-geometry transfer, per S56's
  reading) has any structural feature that makes it obviously
  viscosity-fragile or viscosity-robust — this is a question about the
  construction's own geometry, checkable without new heavy machinery, in
  the same spirit as S54's swirl-radius-cap argument.
- State plainly whether this is currently a tractable question for this
  programme or a genuinely open, hard question in its own right (much of
  the vanishing-viscosity literature is itself unresolved for general data).

## Deliverable

End with:

- A clear report on Thread 1: either a genuine general theorem (stated and
  proved, with hypotheses), or a precise statement of why the
  generalisation fails or is not clean, with the specific obstruction named.
- A clear report on Thread 2: what is known, whether this construction's
  specific geometry offers a tractable angle, and an honest assessment of
  whether pursuing it further is worth a dedicated follow-up session.
- A recommendation: a scoped next work package for whichever thread (or
  both) looks tractable, with named failure modes stated in advance — or an
  honest "neither is tractable at this programme's current scale" if that
  is what the reconnaissance finds. **Do not manufacture a target.**

## Standing discipline, unchanged

- **Session-tag collision has happened six times.** `git log --oneline -1`
  and `ls SESSION-LOG/` before writing anything — note that `S58` is used by
  a commit message but has no session log (S59 already flagged and
  corrected this); check the true state, not the assumption that S60 is
  free.
- Fetch every source; do not recall.
- **This programme still never introduces forcing into its own regularity
  arguments.** This handover's entire content is about *studying* forcing
  requirements in constructions that are not this programme's own — keep
  that line sharp in anything written up.
- If any proof-document edit is made: `check_proof.py` (exit 0) and a
  three-pass `pdflatex` compile (0 errors, 0 new LaTeX warnings, bookmark
  warnings at or below the current baseline — measure it yourself).
- Do not touch `book/`, `layer3/`, `layer4/`, `results/` unless a genuine
  general theorem from Thread 1 is produced and belongs in `proofs/` — in
  which case follow the standard isolated-worktree discipline.
- Do not commit or push. Report your worktree path, if any, for independent
  verification before merge.
