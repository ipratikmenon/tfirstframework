# S56 — Strategic audit: findings, assessment, recommendation

**Date:** 2026-09-12
**Tip at start and at finish:** `fc32bf6` ([S55] plan: S56 strategic audit). Highest
session log = S55. **S56 is correct**, verified twice, working tree clean, no other
agent branch present.
**Scope:** literature and strategy only. Nothing written to `proofs/`, `book/`,
`layer3/`, `layer4/`, `results/`. Nothing committed or pushed.
**Read before writing:** `CLAUDE.md`; `HANDOVER-S56-AUDIT.md` in full;
`sec:walls-s47` in its current state (lines 14315–15029 of
`proofs/claim_a_3d_proof_attempt.tex`, Walls 0–10 plus status table and the two
closing remarks); session logs S48, S51, S52, S53, S54, S55; the S50 KEY FINDING
in `PROGRESS.md` (S49 and S50 have no log files — confirmed, as S55 recorded).

---

## 0. Executive summary

1. **Nobody has removed or weakened the logarithm unconditionally.** Confirmed by an
   exhaustive metadata sweep of the BKM literature and by targeted searches. That was
   the expected answer and it holds.
2. **But one line, invisible to this programme until now, targets exactly the
   logarithm, by a technique this programme has never used, in exactly this
   programme's own variables.** Zoran Grujić, `arXiv:2607.08866v3` (9 Jul 2026,
   revised **9 Sep 2026 — three days ago**) plus companion `arXiv:2609.05720v2`
   (4 Sep 2026, revised 9 Sep 2026). The mechanism recasts the stretching eigenvalue
   `α = ê·Sê` as a **singular-integral commutator** and pays for it with a
   Coifman–Rochberg–Weiss estimate, so that the vorticity **direction** needs only
   log-weighted `bmo` regularity — a Dini-failing condition, far below
   Constantin–Fefferman. Its direction equation is **term-for-term identical** to this
   document's `eq:dir-eq`. Status: unrefereed preprints, one author, sweeping
   conclusion, and this author's previous headline claim in the same area was analysed
   in the literature and found not to deliver. **Priced, not endorsed.**
3. **Nothing new at all on Liouville theorems for ancient Euler solutions.** The Type-II
   literature is still exactly the four Seregin papers S53 read. That branch has not
   moved, and it has quietly got *worse* for an unrelated reason: an unforced 3D Euler
   finite-time singularity from `C_c^∞` data is now claimed (and Lean-certified), which
   means the unconditional Liouville theorem S53 named as the owner of the Type-II
   question cannot exist in class-free form.
4. **The field moved while this programme was looking at walls.** Tao, 7 Sep 2026, on
   Alpöge–Buckmaster: forced blowup for IPM, 2D Boussinesq and **3D Euler**, and
   removing the forcing "looks very feasible … in the near future", with the method
   having "a high likelihood of also extending to Navier-Stokes as well". Every
   AI-assisted effort found — OpenAI, Alpöge–Buckmaster, DeepMind, Ganeshram et al. —
   targets **blowup**. None targets regularity. The handover's guess about that
   asymmetry is correct and there is a structural reason for it.
5. **The internal toolkit is exhausted** in the sense S52 and S53 used the word: three
   named residues remain, all of them capped below the Prize by Wall 3, and none of
   them cheap.
6. **Recommendation: one scoped session (WP-A), gated, with a real chance of a clean
   negative** — price Grujić's hypothesis against this programme's own no-concentration
   family, with the S50 discipline (a new target must be *proved* cheaper, not assumed
   cheaper). Plus an honest statement of the different research question the owner may
   or may not want, stated concretely in §4.3.

---

## 1. Deliverable 1 — the literature sweep

Everything below was fetched. Nothing is recalled. Where I could not read the primary
text I say so.

### 1.1 The highest-value question: has anyone weakened or removed 𝓛?

**Unconditionally: no.** Evidence, in the order it was gathered:

- An arXiv metadata sweep of the entire Beale–Kato–Majda corpus (`all:"Beale-Kato-Majda"`,
  **75 items total**, sorted newest first) returns, for 2024–2026, only BKM-type
  criteria for *other* systems: Hall/electron-MHD (`2407.04314`), free-boundary Euler
  with surface tension (`2507.10032`), Oldroyd-B, barotropic compressible NS, rigid
  body in perfect fluid. **No item addresses the endpoint logarithm in
  `‖∇u‖_∞ ≲ ‖ω‖_∞ log(e+‖u‖_{H^s}/‖ω‖_∞)`.**
- Targeted searches on Kozono–Taniuchi refinements return only *logarithmically
  improved regularity criteria* — a different object (weakening a Serrin exponent by a
  log), not the removal of this log.
- The one relevant optimality statement remains Kanamaru, *Optimality of logarithmic
  interpolation inequalities and extension criteria to the Navier–Stokes and Euler
  equations in Vishik spaces*, **J. Evol. Equ. 20 (2020), 1381–1397**, DOI
  10.1007/s00028-020-00559-0: the Vishik space "may be the largest normed space that
  satisfies the logarithmic interpolation inequality". That is 2020 and it points the
  wrong way for us.

This is exactly what Wall 1 predicts, and this programme already owns the sharpest
statement of *why*: `prop:s51-poloidal-sharp` (S52) exhibits a smooth compactly
supported axisymmetric divergence-free field, supported inside `B_{r*}` about the axis,
with `‖∇u‖_∞ ≥ cN‖ω‖_∞` and `‖u‖_{H³} ≤ C 2^{2N}`. **The inequality is saturated for
general fields, so "removing the logarithm" is not a thing anyone can do.** The only
routes are the two Wall 1(d) already permits: use structure special to NS solutions, or
stop taking a sup-norm of a Calderón–Zygmund image.

**And one line does exactly the second thing. This is the find of the audit.**

#### Grujić, `arXiv:2607.08866v3` — "Logarithmic Depletion of Vortex Stretching and Singularity Evasion in the 3D Navier-Stokes Equations"

Submitted 9 Jul 2026; **v3 dated 9 Sep 2026**; math.AP; no journal reference; comment
field: "cross-referencing a companion paper arXiv:2609.05720 plus a couple of
clarifications". Abstract, verbatim:

> "We present a geometric-analytic mechanism for the suppression of finite-time
> singularities in the 3D incompressible (unforced) Navier-Stokes equations for
> critical point singularities exhibiting $L^{3/2,\infty}$ spatial concentration of
> vorticity. We demonstrate that if the vorticity direction resides locally in a
> logarithmically weighted space of bounded mean oscillations, $\mathrm{bmo}_{1/|\log
> r|}$ — a space failing the Dini condition and thus permitting wild oscillatory
> defects — the non-linear vortex stretching is fundamentally depleted. By isolating a
> unidirectional geometric cancellation, we recast the stretching eigenvalue as a
> singular integral commutator. Utilizing a localized Coifman-Rochberg-Weiss estimate
> coupled with dyadic BMO tail bounds, we prove the stretching potential vanishes as a
> logarithmic envelope on shrinking super-level sets. This depletion forces the
> vorticity magnitude into a sub-critical Lorentz-Zygmund space via interpolated De
> Giorgi energy method. The logarithmic gain is subsequently transferred to the velocity
> field, forcing the geometric scale of local 1D sparseness below the uniform radius of
> spatial analyticity, ultimately averting the finite-time blow-up via the harmonic
> measure maximum principle."

Main theorem as rendered (Theorem 7.4):

> "Let $u_0\in L^\infty(\mathbb R^3)$ and consider a unique, spatially analytic solution
> $u$ on the interval $(0,T_*)$ where $T_*$ is the first possible singular time. Suppose
> that the vorticity conforms to the critical concentration profile uniformly on
> $(T_*-\epsilon,T_*)$, in particular $\omega\in L^\infty((T_*-\epsilon,T_*);
> L^{3/2,\infty}(\mathbb R^3))$. If the vorticity direction satisfies
> $\xi\in L^\infty((T_*-\epsilon,T_*);\mathrm{bmo}_{1/|\log r|}(\mathbb R^3))$, then
> $T_*$ is not a singular time and the blow-up is averted."

with

> "$\|f\|_{\mathrm{bmo}_\phi}=\|f\|_{L^\infty}+\sup_{x,r}\frac1{\phi(r)}\frac1{|B_r(x)|}
> \int_{B_r(x)}|f(y)-c_{B_r(x)}|\,dy$, where $\phi(r)=1/|\log r|$",

and the commutator step rendered as: the Coifman–Rochberg–Weiss theorem gives that
"singular integral commutators with BMO symbols are bounded on $L^p$ for $1<p<\infty$",
extended to `L^{3/2,∞}` by interpolation, "yielding the localized logarithmic smoothing
bound `‖α‖_{L^{3/2,∞}(B_R)} ≤ C_0/|log R|`".

#### Companion, `arXiv:2609.05720v2` — "On Decay of the Local Mean Oscillations of the Vorticity Direction in Critical Navier-Stokes Flows"

Submitted 4 Sep 2026; revised 9 Sep 2026. Abstract, in part, verbatim:

> "We isolate and analyze the geometric PDE governing the evolution of the vorticity
> direction in the 3D incompressible (unforced) Navier-Stokes equations (NSE),
> restricted to the case of a critical spatial point singularity where the vorticity
> magnitude concentrates as $O(|x|^{-2})$, inhabiting the critical Lorentz space
> $L^{3/2,\infty}$. The PDE consists of the Harmonic Map Heat Flow (HMHF) into the
> sphere supplemented with the fluid transport, cross-diffusion and tangential strain.
> The question is whether the NSE mechanics can propagate logarithmic decay of the local
> mean oscillations of the direction — the condition $\xi\in\mathrm{bmo}_{1/|\log r|}$
> which (in this setting) was shown in the companion paper to prevent finite time
> blow-up. … This yields a transfer theorem … Crucially, the strain enters the direction
> equation only through its tangential component $P_{\xi^\perp}S\xi$, which vanishes
> precisely when the direction is an eigenvector of the strain tensor. Since this
> includes the eigenvector carrying the maximal stretching, the result is in contrast to
> the classical geometric regularity criteria which are built on depleting the full
> vortex-stretching term and thus confined to configurations of weak stretching."

**Why this is the find, stated precisely and checked against our own document:**

1. **It is our own equation.** `lem:direction-eq` (line 5349) states
   `D_t ê = νΔê + 2ν(∇log|ω|·∇)ê + ν|∇ê|²ê + P_{ê⊥}((ê·∇)u)`. Grujić's decomposition
   — HMHF (`νΔê + ν|∇ê|²ê`) + transport (`D_t`) + cross-diffusion
   (`2ν(∇log|ω|·∇)ê`) + tangential strain (`P_{ê⊥}Sê`) — is that equation, term for
   term. (`P_{ê⊥}((ê·∇)u) = P_{ê⊥}(Sê)` exactly, since the antisymmetric part
   contributes `½ω×ê = 0`.) The statement that the tangential strain vanishes when `ê`
   is an eigenvector of `S` is the general-geometry form of this document's own
   `prop:s51-aligned` (S52) and `cor:s51-e0`.
2. **It is our scaling regime.** `ω ~ |x|^{-2}` spatial concentration is
   `M·(r*)² = 1` — the criticality this programme recorded at S32 and uses to define
   `r* = M^{-1/2}`.
3. **It attacks the wall where Wall 1(d) says attack is permitted.** Wall 1's
   Consequence clause reads: *"Permitted: any reformulation in which ∇u is not recovered
   by an L^∞ Calderón–Zygmund bound."* A commutator/`H¹`–BMO pairing is precisely such a
   reformulation. Grujić's related note `arXiv:2511.00725v3` ("On taming Moffatt-Kimura
   vortices of doom in the viscous case", rev. 10 Jun 2026) names the tool in the open:
   the mechanism "originates in certain analytic cancellation properties of the
   vortex-stretching term in the sense of **compensated compactness in Hardy spaces**".
4. **It is the only such line in existence.** arXiv metadata search
   `all:"commutator" AND all:"vortex stretching"` returns **three** items, of which
   `2607.08866` is the only serious one; `all:"Hardy space" AND all:"vorticity" AND
   all:"Navier-Stokes"` returns one unrelated item. This technique has one practitioner.

**Now the discipline, which matters more than the excitement.**

- These are **unrefereed preprints by a single author**, with a conclusion ("averting
  the finite-time blow-up") of the class this programme has learned to distrust. No
  independent verification was found.
- **The same author's previous headline claim in the same area was analysed and did not
  survive.** Albritton–Bradshaw, `arXiv:2110.02187v2`, exists precisely to "analyze the
  claims in [`arXiv:1704.05546`, `arXiv:1911.00974`] that *a priori* estimates on the
  sparseness of the vorticity and higher velocity derivatives reduce the 'scaling gap'".
  Their verdict, as rendered from the text: *"While the above interpretation is
  suggestive, still, the relationship (1.21) is informal and should not be taken
  seriously. Rather, to better demonstrate that the scaling gap is not reduced, we
  analyze the sparseness of ∇^k u in concrete blow-up scenarios"*, concluding from their
  Proposition 4.1: *"Within the class of examples we present here, membership in
  Z^{(k)}_{ᾱ_k} is not stronger than the condition of finite kinetic energy
  u ∈ L^∞_t L²_x."* The `2607.08866` argument routes its final step through exactly that
  sparseness/analyticity machinery ("forcing the geometric scale of local 1D sparseness
  below the uniform radius of spatial analyticity … via the harmonic measure maximum
  principle"). **This is a load-bearing prior and it is negative.**
- **The conclusion is conditional on a hypothesis in the family this programme has
  already been burned by twice.** `bmo_{1/|log r|}`-regularity of the direction is a
  no-concentration hypothesis. S50 proved that saturation, advertised as cheaper than
  Wall 1, *is* Wall 1; S52 proved that (AX) is "a no-concentration hypothesis of the
  same family as saturation" and made no cheapness claim. The correct first move is not
  to adopt Grujić's hypothesis but to **price it** against `σ*`-decay, (LC) and (AX).

**A dictionary that makes the pricing a one-session job** (sketch, to be verified, not
claimed): with `σ*(t) = (r*)^{-1}∫_{B_{r*}(x*)}|∇ê|²1_{|ω|≥M/4}` (`def:sigma-star`),
Poincaré on `B_{r*}` gives
`⨍_{B_{r*}}|ê − ê_{B_{r*}}| ≤ C r*(⨍_{B_{r*}}|∇ê|²)^{1/2} ≈ C(σ*)^{1/2}`,
since `∫_{B_{r*}}|∇ê|² = r*σ*` and `|B_{r*}| ≈ r*³`. So **the local mean oscillation of
the direction at the critical scale is controlled by `(σ*)^{1/2}`**, and Grujić's
hypothesis *at that one scale* reads `σ* ≲ 1/(log M)²` — a quantified form of this
programme's `σ*`-decay, which S37/S38 established is the single deep gap and is
*necessary*. The two hypotheses live on one axis; which is weaker, and at what scales,
is a finite question. The caveats are real and must be handled, not waved: `σ*` carries
the indicator `1_{|ω|≥M/4}` while `bmo` does not, and Grujić demands the bound at *every*
scale `r`, not only at `r*`.

### 1.2 The second highest-value question: ancient-Euler Liouville theory

**Nothing new. This is an honest negative and it is the expected one.**

- `abs:"Type II" AND abs:"Navier-Stokes"` returns 14 items; the only relevant ones are
  `2304.04045` (2023), `2402.13229` (2024), `2507.08733` (2025), `2606.29468` (Jun 2026)
  — **exactly the four Seregin papers S53 read in full.** Nothing since.
- `abs:"ancient solutions"` (204 items, newest first) is, for 2024–2026, essentially
  entirely geometric flows: mean curvature flow, Ricci flow, Yamabe flow, fast
  diffusion, parabolic Monge–Ampère. **No ancient-Euler work at all**, under that
  vocabulary or the "eternal solution", "steady Euler rigidity" or "vortex dynamics"
  vocabularies I also tried. Searches on rigidity of stationary Euler return 2D
  sector-domain results only (`2512.18700`).
- Adjacent and genuinely new, but not this: **Pineau–Vicol, `arXiv:2607.09619v2`**
  (10 Jul 2026, rev. 6 Aug 2026), "On rotated backwards self-similar solutions of the
  incompressible 3D Navier-Stokes equations" — Liouville-type triviality for *rotated*
  backwards self-similar solutions under a **Type-I upper bound**, extending
  Nečas–Růžička–Šverák '96 and Tsai '98 off `α=0`, plus a genuinely useful local
  criterion: *"if the solution satisfies a Type I upper bound in a unit parabolic
  cylinder, and there is a single time-slice at which the solution is locally
  approximately self-similar, then the top-center of the parabolic cylinder is a regular
  point"*, all via *"a robust weighted-L² framework … quantitative and … not sensitive
  to whether the Bernoulli head pressure satisfies a maximum principle, which was a key
  obstruction in previous works."* This is Type-I, so Wall 3 caps it, but it is newer and
  stronger technology than the NRS/ESS pointer `sec:closure-s38` records, and it is the
  right citation if the §39 profile branch is ever revived.

**But the branch got worse, for a reason S53 could not have known.** From the primary
PDF (`https://cdn.openai.com/pdf/315b36cd-ec98-4023-8342-93345194ece1/euler.pdf`,
fetched and text-extracted this session), **"Finite Time Blowup for the Euler Equation",
OpenAI, Theorem 1.1**:

> "There exists $u_0\in C^\infty_{c,\sigma}(\mathbb R^3)$ such that
> $0<T_*(u_0)<\infty$. Its smooth Euler solution satisfies
> $\limsup_{t\uparrow T_*}\|\nabla u(t)\|_{L^\infty}=\infty$,
> $\int_0^{T_*}\|\operatorname{curl}u(t)\|_{L^\infty}\,dt=\infty$."

with `C^∞_{c,σ}(ℝ³) = {u_0 ∈ C_c^∞(ℝ³;ℝ³) : ∇·u_0 = 0}` and the equations stated
**unforced**: `∂_t u + (u·∇)u + ∇p = 0`, `∇·u = 0`. Mechanism, from its table of
contents and §2: "a localized oscillation over a smooth Euler flow", "amplification and
transfer of the wave geometry", "choosing the scales and iterating the construction" —
i.e. the Córdoba–Martínez-Zoroa iterated-high-frequency scheme, **not** a self-similar
profile.

Two consequences for our ledger, stated with care:

1. **Wall 9 does not govern this paper.** Wall 9 is "forced blowup constructions cannot
   be de-forced", classified *proved impossible as an import route* because "the residual
   is the forcing by construction". The Euler paper is unforced. If a future session
   touches it, it must not cite Wall 9 as covering it. (Wall 9 continues to govern
   `\cite{OpenAI2026}`, the NS paper, exactly as written.)
2. **The tool S53 named cannot exist in class-free form.** S53's verdict was that closing
   Type-II in general requires "an unconditional Liouville theorem for ancient Euler
   solutions". If Euler singularities exist from `C_c^∞` data, then Euler dynamics
   supports the very degeneration such a theorem would have to exclude, so any such
   theorem must be sharply class-restricted — and the classes that arise from real Euler
   singularities are, on the evidence of this construction, iterated-oscillation objects
   across infinitely many scales, which is the worst possible input for a compactness or
   rigidity argument. **I do not claim that a nontrivial ancient Euler solution has been
   produced** — extracting one needs a rate and uniform bounds this construction does not
   advertise. What I claim is the weaker and safe statement: the hoped-for unconditional
   theorem is now known to be unavailable in the form S53's verdict implicitly wanted, and
   informed pessimism about the class-restricted versions is justified.

**Status caveat, mandatory.** OpenAI's results are Lean-formalized
(`github.com/openai/NavierStokesAndEuler`, two commits, `Euler.lean`,
`NavierStokes.lean`, `formalization.yaml`) but **not independently verified**, the Clay
Institute has not accepted them, and there is a live priority dispute with Alpöge and
Buckmaster, whose forced-Euler blowup (Lean-verified, posted 22 Aug 2026) preceded them.
Lean certifies a formal statement, not that the statement is the intended one — the same
caution this programme applied to its own compile at S47.

### 1.3 What actually changed: the field's centre of mass

Primary source, read in full: Terence Tao, **7 September 2026**, *"Finite time blowup
with smooth forcing term for the incompressible porous medium, Boussinesq, and
incompressible Euler equations"*
(`https://terrytao.wordpress.com/2026/09/07/finite-time-blowup-with-smooth-forcing-term-for-the-incompressible-porous-medium-boussinesq-and-incompressible-euler-equations/`).
Verbatim:

> "It is now widely expected that it should be possible to construct smooth initial data
> and smooth forcing term that would make these equations develop singularities in
> finite time; and it should even be possible to do without the forcing term. While
> these authors do not quite achieve these goals yet, they have made enough of a
> breakthrough that it looks very feasible to complete these goals in the near future."

> "(The first of these equations was already handled by Córdoba and Martínez-Zoroa, but
> Alpöge and Buckmaster found a variant of their method that also extended to the other
> two equations, and **has a high likelihood of also extending to Navier-Stokes as
> well**.)"

> "the actual solving of these problems is only a proxy goal for the primary goal of
> developing mathematical understanding and insight."

and, in an edit to the same post, the third AI-assisted line:

> "there is now also an independent preprint of Ganeshram, Duruisseaux, and Anandkumar
> that has made a significant advance on the other major approach to finite time blowup
> … For the Euler equations (with no forcing term or boundary), they have used a
> physics-informed neural network (PINN) to locate a numerically stable candidate
> solution; though actually establishing its stability … remains a major challenging
> task".

That is `arXiv:2609.10867` and `arXiv:2609.10860`, both **9 Sep 2026**, which reduce the
nonlinear-stability proof "to a large but finite collection of explicit estimates and
computable constants."

**This is the single most decision-relevant fact in the audit.** It is not a proof that
(A)/(B) is false. It is a statement about where the informed expectation and essentially
all of the current effort now sit. A programme whose stated goal is (A)/(B) should decide
deliberately whether to continue, rather than by default.

### 1.4 The four methodologies the handover named

**(a) Computer-assisted / interval-arithmetic proofs.** Mature and accelerating on the
*blowup* side: Chen–Hou for 3D Euler in an axially periodic cylinder with an impermeable
wall, smooth data (cited from the OpenAI Euler paper's own historical section, and
published *Multiscale Model. Simul.*, 2025); DeepMind's `arXiv:2509.14185` "Discovery of
Unstable Singularities" (17 Sep 2025; Wang, …, Buckmaster, Georgiev, Gómez-Serrano, Lai
— i.e. the *Wang et al.* of this repo's `CLAUDE.md` §7.3), which reaches "near
double-float machine precision" and states that "This level of precision meets the
requirements for rigorous mathematical validation via computer-assisted proofs", and
reports "a simple empirical asymptotic formula relating the blow-up rate to the order of
instability" — i.e. the successor of this programme's own λ-family (`λ = 1.1808`,
`0.6057`, `0.4703`, …); then `arXiv:2511.22819` and the two Ganeshram–Duruisseaux–
Anandkumar papers above.
  On the *regularity* side it exists and is honest about its ceiling:
  **Brunk–Giesselmann–Tscherpel, `arXiv:2509.25105`** (29 Sep 2025), "A posteriori
  existence of strong solutions to the Navier-Stokes equations in 3D" — a fully
  computable a posteriori criterion, built on the ESS `L^∞L³` criterion, that verifies
  existence of strong solutions by ruling out blow-up **on a time interval**, "limited to
  short time intervals". **Assessment: computer-assisted methods can exclude blowup for
  fixed data on a fixed interval; they cannot prove global regularity.** Handover item
  (a) is therefore not a route to (A)/(B), and the programme should not build there.

**(b) Tao-style quantitative regularity.** Alive and moving, and it is the *only* other
methodology found that does not route through a sup-norm of a CZ image (it routes
through Carleman inequalities):
  - Barker, `arXiv:2510.20757v3` (23 Oct 2025, rev. 11 Aug 2026) — first quantitative
    classification of potentially singular solutions "at *any* given time in the region
    of potential blow-up times" for approximately axisymmetric data, "in principle
    amenable to numerical testing", motivated explicitly by Hou's numerical candidate;
    also "the first such result establishing sufficient conditions for blow-up from the
    right in terms of blow-up from the left."
  - Barker–Popkin, `arXiv:2602.09951` (10 Feb 2026) — localisation of a slightly
    supercritical Orlicz criterion; new quantitative estimates for the *forced*
    equations.
  - Hu–Nguyen–Nguyen–Zhang, `arXiv:2411.06483v4` (rev. 20 Aug 2026) — quantitative
    bounds in endpoint critical Besov spaces, building on Tao's triple-log.
  **Assessment:** genuinely different in shape, genuinely active, and after seven years
  it has moved the critical-norm lower bound from triple-log to triple-log-with-better-
  constants. It is not a fast route, and entering it means competing with Barker,
  Palasek and Prange on their own ground — the same verdict S52 reached about entering
  the axisymmetric Γ literature.

**(c) Energy-flux / cascade arguments. Retire this item; it is closed by citation.**
Palasek, `arXiv:2605.13827` (13 May 2026), "Finite-time blow-up in an elementary model of
the 3D Navier-Stokes equations": finite-time blowup in an Obukhov-type shell model
`X_k' = −νN_k²X_k + N_{k-1}^α X_{k-1}X_k − N_k^α X_{k+1}² + f_k` with smooth data and
forcing (Theorem 1.3, `α>2`), **and inviscid unforced blowup (Theorem 1.8, `ν=0`,
`α≥1`, `f=0`)**. The point that matters for us is his §4 comparison with Tao's averaged
NS: this model "contains only nonlinear interactions that are organically represented in
the true Euler/Navier-Stokes nonlinearity", both the high-high-low and the Obukhov
interaction corresponding directly to standard PDE terms, whereas Tao's model has "a
complicated array of nonlinear interactions with different weights that do not have any
clear counterpart in the PDE nonlinearity". **So: an argument that sees only the energy
budget and the cascade cannot distinguish 3D NS from an object that blows up, and the
gap between model and equation is now much narrower than it was for Katz–Pavlović or
Tao.** Any global-inequality/flux route must import spatial-geometric information — which
is mildly supportive of the *kind* of mechanism this programme chose at S32, and fatal to
the alternative the handover asked about.

**(d) Helicity.** No live literature. `all:"helicity" AND all:"regularity" AND
all:"Navier-Stokes"` returns 7 items, newest relevant `2410.00813` (helicity conservation
by weak Euler/MHD solutions, not a regularity mechanism). **WP4's negative twist-depletion
measurement (S51) stands, and there is nothing external to build on.** Close the item.

### 1.5 An additional constraint nobody in this programme has recorded

**Cheskidov–Dai–Palasek, `arXiv:2511.09556v2`** (12 Nov 2025, rev. 12 Jan 2026),
"Instantaneous Type I blow-up and non-uniqueness of smooth solutions of the Navier-Stokes
equations". For *any* smooth divergence-free data they construct a solution exhibiting
Type-I blow-up of `‖u‖_{L^∞}` at `T_*>0` while remaining "smooth in space and time on
`𝕋^d×([0,T]∖{T_*})`", with rates as rendered: `‖u(t)‖_{L^∞} ≤ C/√(t−T_*)`,
`‖∇u(t)‖_{L^∞} ≤ C/(t−T_*)` — **saturating the critical scaling**, in all dimensions
`d ≥ 2`. Two rendered statements matter:

- Remark 1.4: the solution "lies in neither `L^∞_t L²_x` nor `L²_t H¹_x` near the
  blow-up time", and (Theorem 2.4) energy blow-up is *necessary* for instantaneous
  Type-I blowup, not an artefact.
- These solutions **satisfy logarithmically-weakened BKM and LPS conditions**:
  `∫_0^T ‖∇u‖_∞ /(1+(loglog(e+‖∇u‖_∞))^c)\,dt < ∞` and the corresponding `‖u‖_∞²` form.

**Lesson, directly applicable to Wall 1 and to the pincer's target.** Any chain of
reasoning toward "no Type-I singularity" that is built only out of BKM/LPS-flavoured
integrability of `‖∇u‖_∞` cannot succeed, because objects satisfying those criteria and
blowing up at exactly the Type-I rate exist. What excludes them is membership in
`L^∞_t L²_x ∩ L²_t H¹_x`. This programme's standing hypotheses are on `𝕋³` with the
energy inequality in force, so nothing here breaks anything — but it is a sharp statement
of *which* ingredient is load-bearing, and it belongs in the walls document if anyone ever
writes there again. I did **not** verify against the document that every estimate in
§33–§44 uses the energy class essentially; that is a check, not a finding.

### 1.6 The AI asymmetry — confirmed, with the reason

The handover asked whether anyone, including via AI-assisted proof search, has targeted
*regularity*. **Answer: no. Every AI-assisted effort found targets blowup.**

| Effort | Target | Certificate |
|---|---|---|
| OpenAI, NS (Sept 2026) | blowup, **forced** | Lean 4 |
| OpenAI, Euler (Sept 2026) | blowup, **unforced**, `C_c^∞` data | Lean 4 |
| Alpöge–Buckmaster (Aug 2026) | blowup, forced: IPM, 2D Boussinesq, 3D Euler | Lean 4 |
| DeepMind `2509.14185` (Sept 2025) | unstable self-similar blowup profiles | machine-precision numerics for CAP |
| Ganeshram–Duruisseaux–Anandkumar `2609.10867/60` (Sept 2026) | stable Euler blowup profile on `ℝ³` | PINN + stability framework |

There is a structural reason and it should be stated rather than lamented: **blowup is an
existence statement and admits a certificate; global regularity is a universally
quantified a priori estimate and admits none.** A formalization pipeline is worth
enormous amounts to the first and almost nothing to the second, until the human-level
argument already exists. So the asymmetry is not a fashion and will not reverse.

### 1.7 A hygiene note

Several 2026 arXiv preprints surfaced by keyword search in this area — titles of the
form "Critical Ledgers and Scale-Defect Cascades", "Finite-Chain CKN-Bad Scale Counting",
"Invisible Defect Cascades", "Endpoint Energy Atoms Force Local Pressure Concentration",
"Stable Finite-Time Singularity Formation for 3D Navier–Stokes via 5D-Lifted Axisymmetric
Reductions" — do not come from identifiable established groups and were **not used**.
I make no claim about them beyond that I did not rely on them. The practical consequence
for this programme is that keyword-only sweeps of math.AP in 2026 have a much lower
signal-to-noise ratio than the ones S46–S53 ran, and author/venue provenance now has to
be checked before a citation is entered.

---

## 2. Deliverable 2 — is anything left inside the toolkit?

**Substantially no, and I decline to manufacture a target.** The full inventory, against
the current `sec:walls-s47` and the S48–S55 residues:

| Item | Status | Why it is not a live target |
|---|---|---|
| `σ*`-decay ≡ H2 (the one deep gap) | **open** | S37 unified them; S38 proved decay *necessary*. Its estimate layer runs through Wall 1. |
| Annulus residual `A_ann` / (LC) | **open** | Wall 4 killed every level-selection route by co-area; Wall 5 killed the implications; `prop:s44-linfty-route` needs *decay*, and S45 measured no onset of it. |
| Localized CF substitute on `ℝ³` | **open and explicitly untried** — Wall 8(c): "Whether some localized substitute exists on `ℝ³` is open and untried" | This is the one genuinely untried internal item in the document. Its payoff is H2/Regime-B on the profile, i.e. inside the S32–S45 line, whose ceiling is Wall 3 ("no Type-I"), which in the one class where the programme has traction is Seregin–Šverák 2009. **Capped below the Prize by construction.** |
| `prop:s50-typeI-scale` | Proved-modulo-**routine** | Housekeeping; `rem:s50-not-e` already says it does not deliver (e). |
| (N), (N′), (A-up), `target:disorder-depletion` | open | Untouched since S39; `target:disorder-depletion` must beat the constant 1 exactly — the marginal-target pattern the S31 audit named as where errors hide. |
| Axisymmetric Γ apparatus | usable, not a target | `prop:s51-boundedness`: any "Γ bounded ⟹ regular" criterion *is* global regularity for the class. |
| `c_K` / band numerics | **proved impossible as evidence** (Wall 6) | — |

**The one internal move that is new, and it exists only because of §1.1:** the
`σ*` ↔ local-mean-oscillation dictionary. It is internal (it uses `def:sigma-star`,
Poincaré, and nothing else), it is cheap, and it is the correct first response to an
external claim of the S50 type. It is a *comparison*, not a theorem, and I am labelling
it as such.

Everything else in the toolkit has either been proved impossible, withdrawn, or reduced
to the same deep gap. That is the same answer S52 and S53 gave, and one more session of
looking has not changed it.

---

## 3. Deliverable 3 — recommendation

### 3.1 The recommendation, in one sentence

**Spend exactly one session (WP-A) pricing Grujić's hypothesis against this programme's
own no-concentration family, and gate everything else on the outcome.** If the answer is
"fifth avatar of the same requirement" — which is the most likely single outcome and
would be a good session — record it and stop. If the answer is "strictly weaker, and the
commutator step is sound", the programme has, for the first time since S31, an external
technique aimed at the exact wall.

### 3.2 WP-A, scoped

**Deliverables.**
1. Read `arXiv:2607.08866v3` and `arXiv:2609.05720v2` **in full from the primary PDFs**
   (not the HTML renderings this audit used for the abstracts and Theorem 7.4), and
   record what is actually proved versus assumed.
2. Settle the dictionary: prove or refute
   `⨍_{B_{r*}}|ê − ê_{B_{r*}}| ≲ (σ*)^{1/2}`, handling the `1_{|ω|≥M/4}` restriction
   honestly, and state the implication in both directions between
   `ξ ∈ bmo_{1/|log r|}` and (i) `σ*`-decay, (ii) `conj:lc-s44` (LC), (iii) (AX) of
   `prop:s51-ax`, (iv) saturation. **The S50 standard applies verbatim: a new target
   must be *proved* cheaper, not asserted cheaper.**
3. Independently check the one step that carries everything: the recasting of
   `α = ê·Sê` as a singular-integral commutator with a BMO symbol, and whether the
   Coifman–Rochberg–Weiss bound applies to it as claimed. This is the step that converts
   `𝓛` into a gain, it is the only genuinely new technique the sweep found anywhere, and
   it is checkable in isolation from the rest of the argument. **Even if the rest of the
   paper is wrong, this step being right is valuable to this programme.**
4. Record the Albritton–Bradshaw prior (`arXiv:2110.02187`) alongside, so that no future
   session adopts the sparseness/analyticity endgame without knowing that the earlier
   scaling-gap claims in the same lineage were analysed and found not to deliver.

**Named failure modes, in advance, in the S46 tradition.**
- **F1 — unadjudicable.** The preprints may be wrong in a way one session cannot settle.
  *Mitigation:* the deliverable is the comparison and the isolated commutator check, both
  of which are valid whichever way the paper falls. Do not attempt a full referee report.
- **F2 — fifth avatar.** The hypothesis turns out to be the same no-concentration
  requirement as `σ*`-decay/(AX)/saturation. *This is the expected outcome, it is a
  legitimate result, and the session should end there*, exactly as S50 did.
- **F3 — rediscovery.** The programme re-derives Grujić's own lemmas and reports them as
  its own, which is S52's failure mode 1 and it fired for real there. *Mitigation:* read
  first, derive second, and attribute.
- **F4 — scope creep.** Attempting to *prove* the `bmo_{1/|log r|}` hypothesis. That is
  the whole problem in new clothing (and is what the companion paper attempts).
  **Forbidden in advance.**
- **F5 — importing the conclusion.** Writing anything into `proofs/` that depends on an
  unrefereed preprint. Forbidden; the session's output is a comparison section that cites
  the preprint as a preprint, or nothing at all.

**Out of scope:** numerics of any kind; `book/`; any new conjecture; any statement about
whether Grujić's theorem is true.

**Definition of done:** a written logical relation with a proof or a counterexample, plus
a one-paragraph verdict on the commutator step, plus a one-line ledger entry. Nothing
else.

### 3.3 The honest (b): the different research question, stated concretely

If WP-A comes back F2 — and it probably will — then the audit's answer to the handover's
question is the one it was prepared to give: **the cheapest remaining moves against
(A)/(B) have all been tried, and what is left needs an external breakthrough this
programme cannot manufacture.** In that case, here is the different research question,
stated concretely enough to accept or decline:

> **How much forcing does a blowup construction actually require, and is the requirement
> removable?**

The reasons this is the right different question *for this programme specifically*:

- It is where the field now is. Córdoba–Martínez-Zoroa and Alpöge–Buckmaster build
  blowup by iterating high-frequency corrections "in a manner that makes the solution
  more singular towards the blowup time **while keeping the forcing term well behaved**"
  (Tao, 7 Sep 2026). Removing the forcing is the stated next step, and Tao calls it "very
  feasible … in the near future".
- This programme has already built, and verified against a primary source, the only
  quantitative de-forcing obstruction I am aware of: `lem:wpt3-forced-swirl`,
  `lem:wpt3-forced-max`, `prop:wpt3-floor` (`|f_θ| ≳ τ^{-(A+1)}`, full leading order in
  its own momentum equation), `cor:wpt3-annulus`, and `cor:wpt2-radius-cap` (an unforced
  axisymmetric core of radius `ℓ_r ≍ τ^{1/2}` can carry only Type-I swirl). S54 and S55
  demonstrated the capability on the hardest available target and got exponent-exact
  agreement with the source by an independent route.
- It is honest about direction of travel: **this is work toward understanding (C)/(D),
  not toward (A)/(B)**, and the owner should accept it as such or decline it. Wall 9's
  prohibition ("forbidden: reading a forced blowup result as evidence about the unforced
  equations in either direction") stays in force and is exactly what makes the *measured
  floor* the interesting object rather than the construction itself.
- Concrete first target if accepted: the OpenAI **Euler** paper, which is **unforced** and
  therefore outside Wall 9 entirely — the first object in this programme's history where
  the geometry can be read directly rather than through a de-forcing prohibition.

**I am not recommending §3.3. I am stating it so the owner can decide.** Starting it
without an explicit decision would be exactly the manufactured-target failure S53 and S55
refused.

### 3.4 What I am not recommending, and why

- **Not** the axisymmetric line (S52: Wall 1 sharp there, Type-I already Seregin–Šverák
  2009, cost of entry is moving off `𝕋³`).
- **Not** Type-II / ancient Euler (S53's verdict stands unchanged; §1.2 makes it worse,
  not better).
- **Not** computer-assisted proof (§1.4(a): structurally cannot deliver global
  regularity).
- **Not** energy-flux/cascade (§1.4(c): closed by Palasek's shell model).
- **Not** helicity (§1.4(d): no literature, and WP4 already returned negative).
- **Not** the Tao/Barker quantitative programme (§1.4(b): real, but entering it means
  competing with its authors using their tools, which is the S52 verdict re-run).

---

## 4. Verification record and non-claims

**Read as primary text this session:** `sec:walls-s47` in full at its current state;
`lem:direction-eq`/`eq:dir-eq` and `def:sigma-star` at lines 5349–5399; the bibliography
key list (59 entries — confirming that Grujić, Bradshaw, Farhat, Albritton, Palasek,
Cheskidov, Pineau, Vicol, Alpöge, Buckmaster, Córdoba and Martínez-Zoroa are **all
absent**); session logs S48, S52, S53, S54, S55 in full and the S50 KEY FINDING; Tao's
blog post of 7 Sep 2026 in full raw text; the front matter and Theorem 1.1 of the OpenAI
Euler PDF, text-extracted from the primary PDF; abstracts and metadata (dates, comments,
journal refs) fetched from the arXiv API for every arXiv item cited here; rendered main
theorems for `2607.08866` (Thm 7.4) and `2110.02187` (Defn 1.1, Thm 1.2, §1.2 conclusion),
and for `2511.09556` (Rem 1.4, Thms 2.3/2.4), `2605.13827` (Thms 1.3/1.8 and §4).

**Derived here, not taken from any source:** the term-for-term identification of Grujić's
direction PDE with `eq:dir-eq`, including `P_{ê⊥}((ê·∇)u) = P_{ê⊥}(Sê)`; the observation
that his `|x|^{-2}`/`L^{3/2,∞}` concentration is this document's `M(r*)²=1`; the
`σ*` ↔ mean-oscillation Poincaré dictionary (**sketch, explicitly not proved**); the
observation that Wall 9 does not govern the OpenAI Euler paper; the structural reason for
the AI blowup/regularity asymmetry.

**Explicit non-claims.**
1. I have **not** verified Grujić's proofs, and this audit takes no position on whether
   Theorem 7.4 is true. It is an unrefereed preprint.
2. I have **not** verified OpenAI's Euler theorem beyond reading Theorem 1.1 and the
   equation display from the primary PDF. It is unverified, disputed on priority, and not
   accepted by Clay.
3. I did **not** read Alpöge–Buckmaster; their preprints are not on arXiv and I relied on
   Tao's post, quoted verbatim, for their content.
4. The `σ*`↔`bmo` dictionary is a one-line heuristic with two named caveats (the
   `1_{|ω|≥M/4}` restriction; single scale versus all scales). It is the *object* of
   WP-A, not an input to it.
5. The claim in §1.2 that "the unconditional Liouville theorem cannot exist in class-free
   form" is an inference from a claimed and unverified theorem, and is flagged as
   informed pessimism, not as a result.
6. Absence of evidence: my searches for ancient-Euler Liouville work covered the arXiv
   metadata index and general web search under several vocabularies. A paper indexed
   under none of them could exist.
