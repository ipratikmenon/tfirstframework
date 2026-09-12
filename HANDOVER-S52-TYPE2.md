# Handover — Type-II reconnaissance (post-S52)

Written 2026-09-12, after all four work packages of `HANDOVER-S46-OPUS.md`
completed (WP2/S48, WP2b/S50, WP4/S51, WP3/S52). Owner's direction, given
directly: **attack Type-II next.** This is the first session of a new line,
not a continuation of the S46 handover, which is exhausted.

## Why this, and why now

Every route this programme has built since S32 — the alignment pincer, the
profile-rigidity reframe, the saturation question, the axisymmetric pivot —
proves things about **Type-I** singularities only, and each closing document
says so in terms. `rem:s39-distance` (line 7703) states the reason
structurally, not as a gap to be filled later: the profile-extraction
rescaling `\eqref{eq:rescaling}` needs `M(t)\lesssim(T-t)^{-1}` to get a
**uniform bound** on the rescaled family at all; without it there is no limit
object, and "a Type-II singularity would rescale to nothing." Wall 3
(`subsec:wall3-typeII`, line 12668) records the same thing as a named,
`\textsc{open}` wall: the whole S32–S45 line's ceiling, even if every open
item in it closed, is "no Type-I singularity" — strictly less than the Clay
problem. S52's axisymmetric ledger (`rem:s51-ledger`, item 2) reconfirms the
vacuity transfers unchanged into that symmetry class.

**This is therefore not a gap this programme has been slowly closing. It is
untouched ground.** No session has ever attempted a Type-II argument, a
Type-II-adapted rescaling, or even a literature survey of what is known.
That is the entire justification for spending a session on reconnaissance
before committing to a technical approach.

## What "Type-II" means here, stated precisely so it is not re-derived wrong

Under the document's convention (`rem:s39-distance`, Wall 3): a putative
first singularity at time $T$ is **Type-I** if $M(t) := \|\omega(\cdot,t)\|_\infty
\leq C_I (T-t)^{-1}$ for some constant $C_I$ — the rate matching the
equation's natural scaling. It is **Type-II** if no such bound holds, i.e.
$\limsup_{t\uparrow T} (T-t)M(t) = \infty$ — blow-up strictly faster than the
scaling-critical rate. (Verify this against the literature rather than
trusting this paragraph; get the exact definitions used by the sources you
fetch, and flag immediately if any source's convention differs from the
document's.)

## The work package: WP-T1, Type-II reconnaissance (literature and framing only)

**This is a one-session, literature-and-writing package. Do not attempt a
new theorem, do not write or run numerics, and do not commit to a technical
programme before the survey below is done and reported.** The failure mode
to avoid is the opposite of S48–S52's: those sessions had a concrete
technical target handed to them. This one does not, and manufacturing one
prematurely (before knowing what exists) is worse than reporting "here is
what is known, and here is why it is hard."

### Deliverable 1 — what is actually known about Type-II for 3D Navier–Stokes

Fetch and quote exactly, do not recall:

- Is **any** partial result — construction, exclusion, or conditional
  criterion — known for Type-II blow-up of the true 3D incompressible NS
  equations specifically (not a toy/shell/modified model)? Survey the recent
  regularity-criteria literature (Escauriaza–Seregin–Šverák and its
  descendants, Seregin's later papers, Barker, Albritton, Palasek — several
  of whom already appear in this document's bibliography for other reasons)
  for anything that bears on Type-II, even partially or under extra
  hypotheses.
- Is Type-II blow-up known or believed to be **possible** for NS at all —
  is there a numerical or heuristic scenario in the literature (Kerr,
  Hou–Li-type, or others) that is explicitly Type-II rather than Type-I?
  The OpenAI/Alpöge–Buckmaster forced constructions this programme already
  catalogued (ch. 17 of `book/`; `sec:s49-two-strands`) produced a Type-II
  profile in one case (`u_\theta \sim \tau^{-(1/2+h)}`, $h>0$) — but under
  **forcing**. Is anything known for the **unforced** case?
- Distinguish sharply between "Type-II is open" (nobody has excluded it) and
  "Type-II is expected" (heuristic or numerical evidence it occurs) — these
  are different states of knowledge and the document must not conflate them.

### Deliverable 2 — is the vacuity fundamental, or specific to this construction?

`rem:s39-distance` states that *this programme's* rescaling is vacuous
against Type-II because it requires the Type-I rate for a uniform bound.
Determine, and prove or cite rather than assert:

- Does **any** compactness-based profile-extraction approach to blow-up
  exclusion face the same vacuity against Type-II — i.e., is "no uniform
  bound ⟹ no limit object" a fundamental obstruction to this entire *style*
  of argument (parabolic rescaling + compactness + rigidity of the limit),
  regardless of which quantity (vorticity alignment, coherence deficit, or
  anything else) the rescaled family carries? If so, Type-II is not a
  further target for *this programme's tools in any adapted form* — it needs
  a different style of argument entirely, and that should be stated as
  plainly as `rem:s39-distance` states the vacuity itself.
- Or, is there a **different, rate-adapted rescaling** — e.g. one that
  rescales by the *actual* $M(t)$ at each reference time rather than assuming
  it decays like $(T-t)^{-1}$, in the manner of Giga–Kohn-type or
  Merle–Zaag-type modulation analysis for other critical/supercritical
  PDEs (semilinear heat equations, harmonic map heat flow, mean curvature
  flow) — that produces a well-posed limit object even when the blow-up rate
  is unknown or faster than critical? If such a technique exists and
  transfers, that is a genuine second work package; if it does not transfer,
  say exactly why (most such techniques rely on an explicit stationary or
  self-similar "bubble" profile to linearize around, which the harmonic-map
  and semilinear-heat literature has and NS may not).

### Deliverable 3 — an honest verdict and, if warranted, a scoped WP-T2

End with one of:

- **(i)** A concrete, well-posed technical target for a follow-up session —
  stated as precisely as `target:disorder-depletion` or hypothesis (R) were
  stated, with named failure modes in advance, in the tradition of
  `HANDOVER-S46-OPUS.md`.
- **(ii)** A reasoned conclusion that Type-II is not tractable by any
  adaptation of this programme's existing machinery, with the specific
  obstruction stated precisely enough that a future session does not have to
  re-derive it. This is a legitimate and useful outcome — the programme's
  most valuable sessions (S41, S44, S48's residue, S50, S52) have been
  exactly this kind of precise negative result.
- **(iii)** A recommendation to consult a different existing branch of
  mathematics not yet touched by this programme (name it specifically, with
  the papers that would need to be read), if the honest assessment is that
  Type-II for NS requires tools this programme has no foothold in yet.

**Do not force a verdict of (i) if the evidence supports (ii) or (iii).**
This programme's credibility rests on exactly that discipline, four times
running now.

## Standing traps, carried forward

- **The rescale-to-nothing claim is a document assertion, not (yet, to this
  session's knowledge) an independently re-derived fact.** Re-derive it
  yourself from `\eqref{eq:rescaling}` before relying on it, the same way
  S48 re-derived the restart arithmetic and S52 re-derived the swirl-block
  algebra rather than trusting the brief. If you find the vacuity is less
  total than stated — e.g. it applies to the *specific* normalisation used
  in `\eqref{eq:rescaling}` but not to every possible rescaling — that is
  itself a significant finding and should be reported prominently, not
  buried.
- **Do not conflate this document's Type-I/Type-II with a different
  convention.** Some axisymmetric literature (flagged at S52,
  `rem:s51-ledger` and its citation corrections) uses "Type I" for a
  *spatial* decay bound rather than the *temporal* rate bound this document
  uses. Check every source's convention explicitly before quoting a result
  as bearing on this document's Type-II question.
- **Session-tag collision has now happened three times running** (S50/S51
  independently claimed by parallel agents, S51/S52 independently claimed by
  WP3 and WP4). Before writing any session log or `PROGRESS.md` entry:
  `git log --oneline -1` on your worktree's base AND `ls SESSION-LOG/` to
  find the true highest existing tag, not the one this handover assumes.
  Verify the current tip and highest tag at merge time too, since other
  agents may run concurrently.

## Working discipline (unchanged from every prior work package)

- **Isolated git worktree.** Copy the live file in; verify
  `wc -l proofs/claim_a_3d_proof_attempt.tex` against the current tip before
  editing anything.
- **Verify against the document, never against a session log or this
  handover, where they might disagree.** This handover's own line-number
  citations should be re-verified, not trusted — every prior handover in
  this programme has had at least one stale line number by the time it was
  read.
- **Fetch literature; do not recall it.** This is a pure literature session;
  the entire deliverable is worthless if built on misremembered statements.
- **Do not overstate.** "Open" and "no known result" are acceptable,
  valuable conclusions. Do not manufacture a technical target to have
  something to hand off.
- Before finishing, if any proof-document edit was made: `check_proof.py`
  (exit 0 required) and a real `pdflatex` compile (three passes, 0 errors,
  0 new LaTeX warnings, bookmark warnings at or below the current baseline —
  check the current baseline yourself, do not assume 109).
- Do not touch `book/`, `layer3/`, `layer4/`, `results/`. This is a
  `proofs/` and `PROGRESS.md`/`SESSION-LOG/`-only session.
- Do not commit or push. Report the worktree path for independent
  verification before merge.
