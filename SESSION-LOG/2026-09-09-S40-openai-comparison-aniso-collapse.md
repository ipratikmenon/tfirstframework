# Session S40 — 2026-09-09 — OpenAI Forced-NS Comparison + Anisotropic-Collapse Adversarial IC

Independent side-track. No bearing on the S31-S39 analytical proof-rigidity program
(profile rigidity reframe, hypothesis (R), gap (A-up), target:disorder-depletion all
untouched). This session did not read or modify proofs/claim_a_3d_proof_attempt.tex.

## Motivation

User surfaced a public claim that OpenAI's internal model produced a proof of
finite-time blowup for Navier-Stokes (widely reported alongside a credit dispute with
Tristan Buckmaster (NYU) and Levent Alpöge (Anthropic), who separately posted forced
finite-time blowup results for the incompressible porous medium equation, 2D Boussinesq,
and 3D Euler, extending the Córdoba–Martínez-Zoroa program). Asked to (1) determine
whether this bears on T-First's actual target, and (2) if there's anything usable, spec
and build it — explicitly without introducing a forcing term, since that would divert
from the Prize equations T-First targets.

## What was read

- `openai/NavierStokesAndEuler` GitHub repo README (Lean 4 formalization): confirms the
  NS result establishes Clay alternatives (C) and (D) — forced blowup on R^3 and T^3 for
  any fixed viscosity — not (A)/(B), the unforced conditions the $1M Prize covers.
- The actual writeup, `navier-stokes.pdf` (166 pages, cdn.openai.com, URL supplied by
  the user), fetched and extracted with `pdftotext -layout` (poppler-utils and `pypdf`
  were not pre-installed in this environment; `apt-get install poppler-utils` was run).
  Read Sections 1-3 (Introduction, Physical description, Proof outline) in full.

## Key facts extracted

- **Theorem 1.1 (exact statement):** for every ν>0, there exist smooth compactly
  supported (space-time) forcing f, a compact K, and smooth (u,p) on R^3×[0,1) solving
  the FORCED NS equations, u(·,0)=0, with kinetic energy uniformly bounded but
  limsup_{t↑1} ‖u(t)‖_L∞ = ∞. The paper states outright that this establishes
  alternative (C) (and, via compact support, (D)) in Fefferman's Clay problem
  statement — not (A)/(B).
- **Mechanism:** an axisymmetric, anisotropic self-similar collapsing vortex core.
  Similarity scales ℓr≍τ^{1/2}, ℓz≍τ^{1/2−h} (0<h<1/100, τ=1−t). Angular Reynolds
  number Re_θ≍τ^{−h}→∞ while radial Reynolds number Re_r=O(1) — viscosity stays in
  full control radially; only the azimuthal/spin-up direction escapes it. The base
  collapsing profile alone leaves an unbounded momentum residual as t→1 (it is not by
  itself a solution of anything); spatially oscillatory ring-shaped "pulses" whose
  nonlinear self-advection generates a mean momentum flux cancel the leading singular
  part of that residual (extending the Córdoba–Martínez-Zoroa vortex-layer-amplification
  program, cited directly as refs [6]-[8]); successive corrections handle the remainder;
  a compact-support cutoff (Section 10) handles localization.
- **Forcing is structurally load-bearing at three points, not decorative:** (a) "an
  exponentially small external force seeds each pulse"; (b) after pulse cancellation,
  what remains of the residual becomes part of f; (c) the compact-support cutoff
  introduces its own force term. f is *defined* as the NS residual R(u,p) — setting
  f≡0 in this construction is not an engineering simplification, it is asking whether
  this exact profile solves the actual unforced Prize equations, which the paper does
  not claim and states it does not establish.
- **No mention of temperature or thermodynamics anywhere** in the 166-page document
  (grepped in full: zero hits for temperature/thermal/thermodynam). ν is a fixed
  constant throughout. No structural overlap with Route 1's μ(T)/k(T)/A(T) framework.
- No citation to Alpöge; only Buckmaster reference is the unrelated 2019
  Buckmaster–Vicol nonuniqueness paper — consistent with secondary reporting that the
  2026 Buckmaster/Alpöge forced-blowup preprints are not acknowledged in this writeup.

## Conclusion on T-First relevance

No effect on any open item in the S31-S39 program (H1 annulus step, σ*-decay/H2,
hypothesis (R), gap (A-up)). Different equation class entirely: forced vs. unforced,
and the opposite question (constructing blowup vs. ruling it out). The one point of
mild interest: even to blow up a *fixed-viscosity* velocity field, OpenAI's
construction needed an externally injected, precisely tuned force to sustain the
anisotropic collapse against dissipation — i.e., the collapse does not sustain itself
under real viscosity alone. Not evidence for or against T-First's claims, just a
data point consistent with the general intuition that dissipation opposes this specific
collapse mode absent external forcing.

## What was built (explicitly not reusing their forcing)

New file `layer3/aniso_collapse_IC.py`: a new adversarial 3D initial-condition family,
inspired only by the *geometry* of OpenAI's construction (an axisymmetric,
independently-tunable-anisotropy swirl-plus-outflow vortex column), never by its
forcing mechanism. No forcing term is introduced anywhere — the IC is fed into the
existing, completely unmodified `route2_3D_jax.run_jax`, which carries f=0 throughout
(same solver used by every other Layer-3 Route 2 experiment, e.g. `wang_profiles_3D.py`).

- `aniso_collapse_ic_3D(N, aspect_ratio, seed, target_E0)` — builds
  `u_theta(r,z) = (r/ellr)*exp(-(r/ellr)^2-(z/ellz)^2)`,
  `u_z(r,z) = ((z/ellz)+bias)*exp(...)` (small asymmetric axial bias, mirroring the
  paper's deliberate breaking of exact z→−z reflection symmetry), `ellz = aspect_ratio
  * ellr`; radial component and exact div-free-ness recovered by the existing Leray
  projection in `set_ic_jax`. Amplitude Γ calibrated so total kinetic energy matches
  the existing `adversarial_min` Wang profile at the same N (not a weaker test by
  construction). Same return schema as `wang_ic_3D()`.
- `run_aniso_route2_jax(aspect_ratio, N, ...)` — models `run_wang_route2_jax()`
  exactly: builds the solver, sets the IC, runs `run_jax`, logs
  `claim_id="claim_l3_aniso_collapse"` to `results/results.db` via the existing
  `_log_result` / `route2_jax_experiments` table.
- `run_aniso_sweep_route2_jax(...)` — runs the full aspect-ratio sweep.
- `ANISO_SWEEP_VALUES = [1.0, 0.5, 0.2, 0.1]` (N=64 production);
  `ANISO_SWEEP_VALUES_HIRES = [0.05]` (flagged as needing N=128 to resolve; not run
  this session).

## Environment setup

- JAX was not pre-installed in this container; installed `jax[cpu]` (0.10.2) — CPU
  backend confirmed working.
- `poppler-utils` was not pre-installed; installed for `pdftotext`.
- `pytest` was not pre-installed; installed.

## Tests

`layer3/test_aniso_collapse_IC.py` — 31/31 PASS. Coverage: output shapes, energy
calibration equality across aspect ratios, ellz scaling with aspect_ratio, ellz floor
at grid resolution, invalid-aspect_ratio rejection, metadata, post-Leray divergence
(<1e-3 at N=16), short-run delta_max>=0 and theta non-negativity at N=16 across the
full sweep, exp_id / claim_id format checks.

No regressions: `layer3/test_wang_profiles_3D.py` (48/48) re-run and still passes
unmodified; no existing file was edited.

## Production experiments (N=64, nu=1e-3, eps_param=0.1, t_end=1.0, max_steps=3000)

Logged to `results/results.db`, table `route2_jax_experiments`, `claim_id =
"claim_l3_aniso_collapse"`:

| exp_id | aspect_ratio | delta_max | theta_negative | verdict |
|---|---|---|---|---|
| EXP-L3-R2-ANISO1000-JAX-064 | 1.0 | 1.50 | False | PASS |
| EXP-L3-R2-ANISO0500-JAX-064 | 0.5 | 1.50 | False | PASS |
| EXP-L3-R2-ANISO0200-JAX-064 | 0.2 | 0.50 | False | PASS |
| EXP-L3-R2-ANISO0100-JAX-064 | 0.1 | 0.50 | False | PASS |

**Result: 4/4 PASS.** delta_max stays strictly positive across the full sweep, though
it weakens monotonically as anisotropy increases (1.50 → 1.50 → 0.50 → 0.50) — a real
trend worth tracking in any follow-up, not just noise (each run's `key_metric` and full
diagnostic history is in `results.db`).

## Honest caveats

- Same epistemic status as `blowup_search_3D.py` (M9): finite N cannot resolve
  arbitrarily thin cores. A PASS down to aspect_ratio=0.1 at N=64 does not rule out
  failure at ratios thinner than the grid resolves. The natural follow-up,
  `EXP-L3-R2-ANISO005-JAX-128` (aspect_ratio=0.05 at N=128), is spec'd
  (`ANISO_SWEEP_VALUES_HIRES`) but was not run this session.
- This module does not reproduce OpenAI's actual analytic profile (their E, U, Π
  functions), which requires solving their Appendix A-C ODE/analytic-continuation
  system — out of scope here. It is a qualitative surrogate sharing the one structural
  feature (independently-tunable radial/axial anisotropy) that their construction
  exploits, used purely as an adversarial IC shape.
- This is a standalone robustness check against the existing Route 2 CZ/delta-gain
  mechanism. It does not feed into, and is not gated by, the S31-S39 analytical
  program's open items.

## Files changed

- `layer3/aniso_collapse_IC.py` (new)
- `layer3/test_aniso_collapse_IC.py` (new)
- `results/results.db` (4 new rows, `route2_jax_experiments` table)
- `PROGRESS.md` (new KEY FINDING S40 block, Last/Second-to-last session pointers
  updated, Session Log table row added, Last updated date)
- `SESSION-LOG/2026-09-09-S40-openai-comparison-aniso-collapse.md` (this file)
