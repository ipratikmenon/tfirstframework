# Handover — WP-T2: does the OpenAI profile survive as an ancient Euler solution?

Written 2026-09-12, immediately after S53 (`HANDOVER-S52-TYPE2.md`, WP-T1).
Owner's direction, given directly after discussing the idea: scope and launch
this specific follow-up, not the Burgers-vortex-barrier alternative also
discussed (that idea is recorded below, in "Explicitly out of scope", in case
a future session returns to it).

## What this is, and what it is not

**This is a narrow, well-bounded exclusion attempt for one specific candidate
singularity profile.** If fully successful, its result is the same *shape* as
Nečas–Růžička–Šverák (1996) excluding self-similar Type-I blow-up: one more
profile family ruled out. **It is not, and must not be reported as, a step
toward (A)/(B).** The full logical chain to global regularity needs (i) Wall 1
closed on the Type-I side (untouched, still open) and (ii) a *general*
Type-II profile-extraction theorem plus a *general*, unconditional Liouville
theorem for ancient Euler solutions (recognized as hard open problems in their
own right — Seregin's own papers, read in full at S53, are conditional on an
assumed growth scenario for exactly this reason). WP-T2 touches none of that
general machinery. It asks a single, checkable question about a single,
explicit object. Say this plainly in the deliverable's own verdict section,
the way S52 and S53 said plainly what their negative results did and did not
achieve.

## Background this session needs, already established and verified

- **S40 / Wall 9**: the OpenAI forced construction cannot be de-forced by
  simply deleting $f$ — the forcing is structurally load-bearing (it seeds
  cancelling pulses and cancels a base residual), and since the construction
  starts from rest ($u^0=0$), $f\equiv0$ trivially gives $u\equiv0$ forever.
  **Direct import is closed. This work package does not reopen it.** What it
  does instead is treat the profile as an object to test against exclusion
  theorems, independent of whether the forced construction that produced it
  can be stripped of its forcing.
- **S49 / `sec:s49-two-strands`**: OpenAI's is one of two strands; only theirs
  (not Alpöge–Buckmaster) claims the full 3D NS equations. Their overview
  states the construction does not address the unforced case.
- **S53 / `lem:wpt1-nueff`, `cor:wpt1-trichotomy`**: for the *isotropic*
  two-parameter rescaling $v(y,s)=A\,u(x^*+Ly,T+ALs)$, the limit solves NS
  with $\nu_{\mathrm{eff}}=\nu A/L$. Normalising a rate $\sup|u|\asymp(T-t)^{-\beta}$
  forces $A=L^{\beta/(1-\beta)}$; for $\tfrac12<\beta<1$, $\nu_{\mathrm{eff}}\to0$
  as $L\downarrow0$ and the limit solves **Euler**. This is proved and was
  independently re-derived before merge — trust it, but note it is derived for
  a *single* isotropic length scale and a *single* amplitude, which the
  profile below does not have.
- **The profile itself, as transcribed from a screenshot in conversation, NOT
  yet independently verified against the primary source:**
  $$u^{(0)} = \frac{V_0(X,\eta)}{r}e_r + q^{-A}E(X,\eta)e_\theta + q^{-A}U(X,\eta)e_z,
  \qquad p^{(0)}=q^{-2A}\Pi(X,\eta),$$
  with $\tau=1-t$, $A=\tfrac12+h$, $D=\tfrac12-h$, $\tau=q(1-\eta^2)$,
  $z=q^D\eta$, $X=r^2/(2q)$, along the ray
  $x_\tau=(\sqrt{2X_{\mathrm{in}}\tau},0,0)$, giving $u_\theta(x_\tau,1-\tau)
  =\tau^{-A}(e_0+O(\tau^{2h}))\to\infty$, with $h<1/100$ reported at S53.
  **Do not trust these formulas as ground truth. Re-derive them from a
  primary source before using them for anything (Deliverable 0 below).**

## Deliverable 0 — get the real profile, from a primary source

Before any mathematics: fetch the actual OpenAI paper (or the most reliable
technical primary source available — a preprint, an official technical
writeup, or the Lean-adjacent technical note if one exists) and extract the
precise leading-order profile, its exact scaling exponents, the precise
domain of the similarity variables $(X,\eta)$, and — critically — **whether
the profile equations for $(V_0,E,U,\Pi)$ are stated explicitly anywhere**,
or only their asymptotic consequence ($u_\theta\sim\tau^{-A}$). If the source
only states the consequence and not the governing profile ODEs/PDEs, say so:
that changes Deliverable 2 from "test an explicit solution" into "reconstruct
the profile equations from the stated asymptotics before testing them,"
which is a harder and different task, and should be reported as such rather
than glossed over.

## Deliverable 1 — an anisotropic rescaling, derived, not assumed

The profile above needs at least these degrees of freedom beyond
`lem:wpt1-nueff`'s two: separate radial and axial length scales
($\ell_r\sim\tau^{1/2}$, $\ell_z\sim\tau^{1/2-h}$ — genuinely different
exponents), and separate amplitude behaviour by velocity component ($u_r$
scaling like $1/r\sim\tau^{-1/2}$, i.e. Type-I; $u_\theta,u_z$ scaling like
$\tau^{-(1/2+h)}$, i.e. Type-II). **This is real work, not a formality.**
Concretely:

- Write the rescaling $v_r(y,s)=A_r\,u_r(\ldots)$,
  $v_\theta(y,s)=A_{\theta z}\,u_\theta(\ldots)$,
  $v_z(y,s)=A_{\theta z}\,u_z(\ldots)$, with independent radial scale $L_r$ and
  axial scale $L_z$ in the spatial arguments, and derive — the same way
  `lem:wpt1-nueff` was derived, by direct substitution into the NS equations
  in cylindrical coordinates, not by analogy — what equation (if any) the
  rescaled family solves in the limit.
- **Check first whether this produces a well-posed limiting equation at
  all.** The Laplacian and the nonlinear advection term will not scale
  uniformly when $L_r\neq L_z$ and the components carry different
  amplitudes; it is entirely possible the naive limit is degenerate,
  ill-posed, or simply not a recognisable PDE. If so, that is the answer to
  Deliverable 1, and Deliverable 2 becomes moot — report this plainly rather
  than forcing a limit to exist.
- If a well-posed limit does exist, identify it exactly (Euler? An anisotropic
  or boundary-layer-type reduced system, e.g. something in the family of
  Prandtl-type or thin-domain limits, given the two very different length
  scales? A hybrid?). Do not assume it is Euler merely because S53's isotropic
  case gave Euler — that conclusion may not transfer.

## Deliverable 2 — test against ancient-solution exclusion, on whatever equation Deliverable 1 actually produces

If Deliverable 1 produces a genuine ancient solution (of Euler or of
whatever the correct limiting system is): check its properties (energy class,
symmetry, non-degeneracy, swirl) against the hypotheses of the four Seregin
papers already read in full at S53 (arXiv:2304.04045, 2402.13229, 2507.08733,
2606.29468) and any other Liouville-type theorem for ancient solutions that
appears relevant once the correct limiting equation is known. Determine
whether the profile:

- **(i) is excluded** by an existing theorem — state which one, and check
  every hypothesis line by line, the way S48 checked Serrin's criterion and
  S52 checked Seregin–Švérák's;
- **(ii) is not excluded**, and state precisely which hypothesis of which
  theorem it fails to satisfy, so a future session does not have to
  re-discover this;
- **(iii) never reaches this stage**, because Deliverable 0 or Deliverable 1
  stopped the chain first. This is a legitimate, complete outcome.

## Explicitly out of scope (do not attempt this session)

The Burgers-vortex-as-barrier/comparison-principle idea discussed alongside
this one is a **different** technical question (whether a nonlocal vorticity
equation admits any useful comparison/maximum-principle argument against an
explicit viscous-core profile) and is not part of WP-T2. Do not merge it in
or attempt it as a shortcut if Deliverable 1 or 2 stalls.

## Standing traps, carried forward

- **Session-tag collision has now happened three times.** Before writing
  anything: `git log --oneline -1` on your worktree's base, and
  `ls SESSION-LOG/` for the true highest tag. This handover's own guess (S54)
  may be wrong by the time you start — verify, and re-verify at merge time.
- **The profile transcription above is unverified.** Do not build Deliverables
  1–2 on it without first completing Deliverable 0.
- **Do not force a verdict.** "The rescaling doesn't produce a well-posed
  limit" and "the profile evades every existing exclusion theorem" are both
  complete, valuable, reportable outcomes — this programme's best sessions
  have been exactly this kind of precise negative result.
- **Do not conflate this with a route to (A)/(B).** State the narrow scope of
  any positive result in the same terms this handover uses above.

## Working discipline (unchanged)

- Isolated git worktree; copy the live `proofs/claim_a_3d_proof_attempt.tex`
  in; verify its line count against the current tip before editing anything.
- Fetch literature; do not recall it. This applies with extra force to
  Deliverable 0.
- Verify against the document, never against a session log or this handover
  where they might disagree — re-derive rather than trust any cited line
  number, exponent, or claim carried over from S53 or from this handover.
- If any proof-document edit is made: `check_proof.py` (exit 0) and a real
  three-pass `pdflatex` compile (0 errors, 0 new LaTeX warnings, bookmark
  warnings at or below the current baseline — check it yourself).
- Do not touch `book/`, `layer3/`, `layer4/`, `results/`. `proofs/` and
  `PROGRESS.md`/`SESSION-LOG/` only.
- Do not commit or push. Report the worktree path for independent
  verification before merge.
