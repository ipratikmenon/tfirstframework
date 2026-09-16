# Handover — CD-2: does the forcing-floor theorem have a second instance?

Written 2026-09-16, following directly from S60's own recommendation.
Owner's direction: launch it.

## What this is

**S60 built one instance of a class theorem** (`thm:cd1-general-floor`,
`cor:cd1-rate`) and verified it against exactly one external object (the
OpenAI axisymmetric profile). A theorem checked against one instance is not
yet known to be a useful *class* theorem — it might simply restate what was
already known about that one profile in more general language. **CD-2's job
is to find out whether the class has a second member**, and, on the way,
to close the one gap S60 flagged as the step most likely to fail: whether
the underlying swirl maximum principle extends from the classical Laplacian
to fractional dissipation at all.

**This is still Thread 1 only.** Thread 2 (vanishing-viscosity fate of the
unforced Euler construction) is closed per S60's own recommendation — do
not reopen it. This is still a documentation/study exercise about existing
forced constructions in the literature, not a step toward (A)/(B), and this
programme still never introduces forcing into its own arguments anywhere.

## Do the fractional extension first — it is the one step that can fail

S60 named this explicitly: attempt the generalisation of the swirl maximum
principle to `|\nabla|^\alpha`-type (fractional) dissipation **before**
looking for a second example, because if it fails, there is no point
searching for examples that would need it.

- `lem:wpt3-forced-swirl` and `lem:wpt3-forced-max` are built on the
  classical Laplacian's sign behaviour at an interior maximum
  (`(\Delta-2/r\partial_r)\Gamma\leq0`). Investigate whether an analogous
  sign property holds for the relevant fractional operator at an interior
  maximum — the Córdoba–Córdoba inequality (`\int|x-y|^{-3-\alpha}(\theta(x)
  -\theta(y))\,dy \geq c|\nabla\theta(x)|^{...}` at the point where `\theta`
  attains its max, or the appropriate maximum-principle statement for the
  fractional Laplacian on a scalar — read the actual Córdoba–Córdoba paper
  and any refinements, do not reconstruct it from memory) survives being
  adapted to the specific singular structure of the axisymmetric swirl
  equation (the `-2\nu/r` term, and whatever the fractional analogue of the
  diffusion operator does to it).
- **If this fails, that is a complete, valuable, reportable outcome** (S60's
  named failure mode H2: "the fractional swirl equation has no clean scalar
  form"). Say precisely where it fails — which structural feature of the
  classical case does not survive — rather than a bare "no."
- **If it succeeds**, you have the fractional generalisation of
  `lem:wpt3-forced-max` itself, which is worth recording as its own result
  independent of whether Deliverable 2 below finds an example to apply it
  to.

## Then look for a second instance

**S60's report names two candidate papers, `arXiv:2309.08495` and
`arXiv:2407.06776`. Only `2407.06776` (Córdoba–Martínez-Zoroa–Zheng) was
independently verified in this programme's prior sessions. `2309.08495` is
unverified — confirm it exists and is the paper S60 meant before using it
for anything; if it does not resolve to a real, relevant paper, do not
guess a substitute — say so and search properly instead.**

For whichever genuine forced axisymmetric (or, if the fractional extension
in Deliverable 1 succeeds, fractionally-dissipative) construction you find
in the literature with an explicit growth rate and core geometry:

- Extract its stated growth exponent and radius scaling the way S54/S55 did
  for the OpenAI profile, from the primary source, not from a secondary
  description.
- Apply `cor:cd1-rate` (or its fractional generalisation, if built) and
  check whether the construction's own stated forcing matches, exceeds, or
  falls short of the derived floor — the same three-way comparison S55 ran.
- **The most likely outcome, named in advance and a completely legitimate
  place to stop (S60's H1): the relevant external class has exactly one
  member (the OpenAI profile), and no second instance exists in the current
  literature to test the theorem against.** If that is what a genuine search
  finds, report it as such. Do not stretch an unrelated construction to fit.

## Named failure modes, carried forward from S60 and binding

- **H1** — the class has one member. Legitimate stopping point, not a
  failure of this session.
- **H2** — the fractional swirl equation has no clean scalar form. Also
  legitimate; report precisely where it breaks.
- **H3** — rediscovery: re-deriving results already in the fractional
  Calderón–Zygmund / Córdoba–Córdoba literature without attribution. Read
  first, derive second, attribute correctly.
- **H4** — drifting into advocacy for (C)/(D) as a research goal in its own
  right. This programme studies forcing requirements; it does not campaign
  for anyone to prove (C)/(D), and does not need to.
- **H5** — blurring the standing prohibition on introducing forcing into
  this programme's own regularity arguments. Everything here concerns
  forcing in *other people's* constructions, examined from outside.

## Standing discipline, unchanged

- **Session-tag collision has happened six times, plus one commit tagged
  without a session log (S58).** `git log --oneline -1` and
  `ls SESSION-LOG/` before writing anything; re-check at merge time.
- Fetch every source; do not recall, including S60's own claims about
  `2309.08495` — verify it directly.
- If any proof-document edit is made: isolated worktree, copy the live file
  in, verify its line count against the current tip first; `check_proof.py`
  (exit 0) and a three-pass `pdflatex` compile (0 errors, 0 new LaTeX
  warnings, bookmark warnings at or below the current baseline — measure it
  yourself).
- Do not touch `book/`, `layer3/`, `layer4/`, `results/`.
- Do not commit or push. Report your worktree path, if any, for independent
  verification before merge.

## Deliverable

- Deliverable 1: the fractional extension, proved or precisely refuted.
- Deliverable 2: a second instance found and checked, or a precise, honest
  report that none exists in the literature currently available.
- A final recommendation: whether this line has anything further worth a
  session, or should be closed here — do not manufacture a reason to
  continue if the honest answer is that it is done.
