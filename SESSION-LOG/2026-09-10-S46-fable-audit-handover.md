# Session S46 — 2026-09-10 — Fable Audit, Strategic Pivot, Opus Handover

Strategic session. No proof-document edits, no code. Deliverable is
`HANDOVER-S46-OPUS.md` plus this record.

## Why an audit

S45 returned the sharpest negative result of the recent programme: the 𝓛-decay
that `prop:s44-linfty-route` requires is absent (slopes ≈ 0, not −1, across a
3.5× range of 𝓛), and the alternative volume route's hypothesis (NB) failed
outright in a measured case. Both of S44's closure routes were then
unsupported by measurement. The owner asked for a Fable audit of how to proceed
and what methods remain for the Prize equations.

A first attempt to run the audit as a Fable subagent failed twice on a monthly
spend limit (HTTP 429). An interim audit was written by the Opus session model
from memory. The owner then switched the session model to Fable and asked for
the audit again; this session is that audit, performed against the document
rather than from memory. One conclusion of the interim audit was overturned by
that check (see W2 below) — which is the reason to audit from sources.

## Findings, verified against `proofs/claim_a_3d_proof_attempt.tex`

**W1 — The logarithm is structural.** `lem:local-enstrophy-typeI` (line ~5664)
introduces `𝓛 = log(e+‖u‖_{H³}/M)` through `‖∇u‖_∞ ≲ M·𝓛`. That is the L^∞
endpoint of the Biot–Savart operator, a Calderón–Zygmund operator unbounded on
L^∞. Every enstrophy budget inherits it, which is why crude K-absorption
(`prop:crude-K`), the annulus core (`rem:s41-shortfall`) and the corrected band
route (`prop:s44-linfty-route`) all miss by the *same* single log. Any estimate
routed through a sup-norm of a CZ image pays it; this architecture routes
everything through one. The BKM endpoint problem, in the programme's own
notation.

**W2 — The profile route relocates the barrier; it does not evade it.**
*This overturns the interim audit.* The interim claimed S39's "profiles are
logarithm-free" was the one branch that structurally escaped W1. Reading the
proof of `prop:transfer-audit(d)`: it is *Proved-modulo-(R)*, and the text says
the log disappears because "the CZ step is replaced by the interior estimate
(R), which is logarithm-free." Hypothesis (R) — scale-uniform interior
regularity of the rescaled Type-I family, `‖∇U‖_∞ ≤ C(c_I)|s|⁻¹` — *is* an
assumed log-free gradient bound. The status table lists it as "Conditional
(Serrin-type)". The heuristic beneath ("no structure below scale |s|^{1/2}") is
stated and unproved. So the barrier was moved from an estimate into a
hypothesis. Whether (R) is already a theorem (via ESS 2003, Seregin 2012, KNSS
2009, NRS 1996, Tsai 1998, CKN 1982) is the cheapest decisive question in the
programme and is WP2 of the handover. Named trap: ε-regularity needs
*smallness* of scaled local energy; Type-I gives *boundedness* by `c_I`.

**W3 — The ceiling is below the Prize, by the document's own statement.**
`rem:s39-distance`: if (A-up), `target:disorder-depletion`, (R), (N), (N′) all
closed, the conclusion is "no Type-I singularity", "strictly less than the Clay
problem". Every estimate uses Type-I rates *essentially*; "a Type-II
singularity would rescale to nothing." The S32–S45 machinery is blind to
Type-II. The S39 status table (~25 rows) has almost every load-bearing row
marked "mod (R)", "mod (R),(N)", "mod (R),(A-up)" or "Conditional on 5 inputs".

**W6 — The c_K numerics were never discriminating.** S38 proved σ\*-decay is
*necessary for blowup*. A regular solver never blows up, so it never exhibits
the decay, so it measures `c_K ≈ 1` — exactly what both hypotheses predict on
any flow the solver can produce. Four sessions (S34, S35, S42–43, S45)
confirmed the non-information. This indicts work directed under this
session's own earlier model instances and is recorded as such.

W4 (level iteration proved impossible, S41) and W5 (c_K/band independence,
S44) were re-confirmed from their labels and are unchanged.

## On the string-theory suggestion

Assessed on its merits. **Fluid/gravity duality** (Bhattacharyya–Hubeny–
Minwalla–Rangamani 2008; Bredberg–Keeler–Lysov–Strominger 2011, *From
Navier–Stokes to Einstein*) is a genuine correspondence: incompressible NS in
*d* dimensions is dual to a vacuum Einstein solution in *d+2*, so an NS
singularity would correspond to horizon singularity formation and regularity
reframes as a cosmic-censorship question. It is not a method: the duality is a
long-wavelength gradient expansion, valid only where a singularity is *not*
forming; it is UV-blind by construction. Cosmic censorship is itself unproven,
and Bizoń–Rostworowski 2011 shows AdS is nonlinearly unstable to black-hole
formation. The scale-criticality that blocks the programme reappears there as
the breakdown of the expansion. Recorded as a reformulation.

The one physics import judged more than metaphor: in renormalisation-group
language, 𝓛 is a *marginal-dimension logarithm* (3D NS at its critical
dimension, as 4D gauge theory). The Prize question "does the marginal coupling
resum harmlessly or pile up" is precisely the `c_K ≤ C/𝓛` question. Not a
technique — the correct mental model of what a genuine resummation argument
would have to do.

## The stellarator suggestion

The owner proposed the stellarator's twisted-torus geometry as a "specific
structure" in the sense of Tao's barrier. Decoded: a twisted torus read as a
vortex tube is the geometric embodiment of **helicity** (linking = twist +
writhe, Călugăreanu–White; the rotational transform is the twist term). The
three dynamics described — counter-phase oscillation, decaying frequency,
unbounded axial stretch — are the Crow instability, the Kelvin-wave cascade,
and the vortex-stretching term: together, the leading *blowup* candidate of
the last thirty years (Kerr 1993, Hou–Li 2006, Yao–Hussain). Studying the most
dangerous configuration is the correct regularity strategy, and twist along a
tube is swirl, which supplies centrifugal/pressure support against core
collapse (the Burgers balance). Warning recorded: Kerr's anti-parallel tubes
are helicity-neutral by symmetry, so global helicity conservation alone is not
a defence; the claim worth chasing is *local* — twist density depleting the
stretching rate `α = ê·S·ê` pointwise. That is a mechanism question, hence
information-bearing in regular flows, unlike c_K. It became WP4b.

## Decision (owner's)

1. Document every route and every wall, precisely and with references.
2. Keep the exact Prize equations and the regularity goal. No forcing, no
   modified equations, no "weak internal forces" (considered and rejected).
3. Pivot the main line to **axisymmetric-with-swirl**: the Prize equations
   restricted to a symmetry class that is still open; Type-I already excluded
   there by the literature (CSTY 2008; KNSS 2009, already cited); the whole
   difficulty is the single source term `∂_z(Γ²)/r⁴`, with `Γ = r·u_θ` obeying
   a maximum principle; the CZ log expected to weaken because ∇u is recovered
   by 2D-type elliptic problems in (r,z) — an expectation WP3 must prove or
   refute with a written estimate before anything is built on it.

## Handover

`HANDOVER-S46-OPUS.md` — five work packages in dependency order (WP0
verification substrate; WP1 walls document; WP2 settle (R); WP3 analytical
pivot; WP4 numerics redirect with `twisted_ring_IC` and helicity diagnostics),
each with deliverables, acceptance criteria and *named failure modes in
advance*; the mandatory verification bar (isolated worktree, then independent
re-derivation of load-bearing steps before merge — the discipline that caught
S41 and S44); cost discipline; an explicit out-of-scope list; and a definition
of done that requires an honest statement of whether the pivot has a live path
or has hit its own wall.

## Files changed

- `HANDOVER-S46-OPUS.md` (new — authoritative work plan)
- `PROGRESS.md` (active milestone redirected to the handover; KEY FINDING S46;
  session pointers; Session Log row)
- `SESSION-LOG/2026-09-10-S46-fable-audit-handover.md` (this file)

No proof-document or code changes.
