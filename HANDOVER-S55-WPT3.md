# Handover — WP-T3: a quantitative floor on the forcing in Wall 9

Written 2026-09-12, immediately after S54 (WP-T2). Owner's direction: turn
Wall 9's qualitative finding ("the forcing is load-bearing, not decorative")
into a derived number — a lower bound on how much forcing the OpenAI
construction must supply, checked against what the construction actually
uses.

## What this is, and what it is not

**This is documentation-sharpening, not a step toward (A)/(B) or a challenge
to (C)/(D).** It makes Wall 9 more precise and, if the numbers check out, more
citable. It does not change Wall 9's classification (`\textsc{proved
impossible} as an import route`), does not touch Wall 1, Wall 3, or the S54
axisymmetric-swirl exclusion, and says nothing new about whether OpenAI's
Theorem 1.1 is correct. State this plainly in the deliverable, the way every
work package since S52 has stated its own scope.

## The one clean lever: the forced swirl equation has no pressure term

This is why the sharpening is tractable at all. `lem:s51-gamma` (already in
the document, S52) gives, for axisymmetric flow, the *pressure-free* equation
$$\partial_t\Gamma+u_r\partial_r\Gamma+u_z\partial_z\Gamma
  =\nu\Bigl(\Delta-\frac2r\partial_r\Bigr)\Gamma,$$
whose unforced maximum principle S54 used. Adding a forcing term
$f=(f_r,f_\theta,f_z)$ to the momentum equations adds exactly
$r f_\theta$ to the right-hand side (since $\Gamma=ru_\theta$ and the
$\theta$-momentum equation gains $f_\theta$). **Re-derive this from the
momentum equations yourself; do not take it on the strength of this
paragraph.** No other component of the swirl equation is affected by $f_r$
or $f_z$ at all — this is a real structural fact about axisymmetric flow, not
an approximation, and it is what makes the swirl direction the only one
where a clean forced floor is derivable this cheaply.

## Deliverable 1 — the general forced floor, as a genuine lemma

Prove, rigorously (not by scaling heuristic), a forced generalisation of
`lem:s51-gamma`: if $M_\Gamma(t):=\|\Gamma(\cdot,t)\|_{L^\infty}$ is achieved
at an interior point $x_m(t)$ and is differentiable there (or work with the
Dini derivative if it is not), then at the maximiser the diffusion term is
$\leq0$ (this is the entire content of why the unforced maximum principle
holds), giving
$$\frac{d}{dt}M_\Gamma(t)\ \leq\ \bigl(rf_\theta\bigr)\bigl(x_m(t),t\bigr).$$
State the precise hypotheses this needs (smoothness of $M_\Gamma$, or the
correct one-sided/viscosity-sense statement if it is not everywhere
differentiable — do this carefully, it is exactly the kind of place a
marginal argument goes wrong, per the S31 audit's own lesson about where
errors hide) and prove it as a standalone lemma with its own label, before
applying it to anything.

## Deliverable 2 — apply it to the specific profile, get a number

Using $M_\Gamma(t)\gtrsim\tau^{-h}$ (S54's finding, already independently
verified — re-confirm the exponent yourself rather than re-trusting it) and
the core radius $\ell_r\asymp\tau^{1/2}$, derive the resulting floor on
$f_\theta$ as an explicit power of $\tau$ near the core. Show your work
fully; do not just state the exponent. Be careful about exactly which norm
and exactly which region this floor applies to (pointwise at the running
maximiser is what Deliverable 1 actually proves — an $L^\infty$ statement
over a whole region is a different, stronger claim and needs its own
argument if you want it).

## Deliverable 3 — check the floor against the actual construction

Fetch the primary source again (`cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`,
already used and verified at S54) and find what it states about the
magnitude, decay, or scaling of its own forcing term $f$ — Theorem 1.1 states
$f\in C^\infty_c$, but the technical sections proving the construction almost
certainly give a quantitative bound or an explicit scaling for the force
used to seed pulses and absorb the residual (S40 and Wall 9 already
identify three places forcing enters: pulse-seeding, residual absorption,
compact-support cutoff — the residual-absorption piece is the one this
lemma's floor should be checked against, since that is the piece defined as
*equal to* the Navier–Stokes residual of the swirl component). Compare:

- If the source's stated forcing scaling **meets or exceeds** the derived
  floor: this is a real, satisfying consistency check — report it as
  confirming, from an independent angle, that S40's "not decorative" finding
  was correctly diagnosed, now with a number attached.
- If the source's stated forcing scaling is **below** the derived floor: this
  is either an error in Deliverable 1/2's derivation (check first, most
  likely explanation) or a genuine discrepancy worth flagging carefully and
  precisely rather than either silently resolved or triumphantly announced —
  the standing discipline in this programme is to state exactly what was
  found and let the mathematics speak.
- If the source does not state an explicit enough bound to compare against:
  say so, and report the derived floor as a standalone quantitative addition
  to Wall 9 regardless — it is still worth having, even unconfronted by a
  comparison number.

## Standing traps, carried forward

- **Session-tag collision has happened four times now.** Before writing
  anything: `git log --oneline -1` and `ls SESSION-LOG/` for the true current
  tip and highest tag. This handover's own guess (S55) may already be stale.
- **Re-derive, do not re-trust.** Every exponent and formula carried over
  from S52/S53/S54 above should be independently re-confirmed against the
  document and the primary source before being used in a new derivation —
  this is not a formality, it is how this programme's last four sessions
  each caught a real error in the one before it.
- **The pointwise-vs-uniform distinction in Deliverable 1 is the most likely
  place for a marginal error.** Be explicit about which one you have proved.
- **Do not let this become a critique of OpenAI's paper.** If a discrepancy
  is found, it is reported as a fact to be checked, not as a refutation —
  their construction is not fully re-derivable or re-verifiable from outside
  in one session, and this programme's own S52/S54 exclusions apply to
  idealised reductions, not to their actual (forced, non-axisymmetric)
  object.

## Working discipline (unchanged)

- Isolated git worktree; copy the live `proofs/claim_a_3d_proof_attempt.tex`
  in; verify its line count against the current tip before editing anything.
- Fetch literature and the primary source; do not recall either.
- If any proof-document edit is made: `check_proof.py` (exit 0) and a real
  three-pass `pdflatex` compile (0 errors, 0 new LaTeX warnings, bookmark
  warnings at or below the current baseline — check it yourself).
- Do not touch `book/`, `layer3/`, `layer4/`, `results/`. `proofs/` and
  `PROGRESS.md`/`SESSION-LOG/` only.
- Do not commit or push. Report the worktree path for independent
  verification before merge.
