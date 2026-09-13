# SESSION LOG — S59 (2026-09-13)

## Second strategic audit: the founding thesis, stochastic representations, kinetic limits

**Type:** literature and mathematical analysis only. No numerics, no proof-document
content, no `book/`, `layer3/`, `layer4/`, `results/`. Nothing committed, nothing pushed.

**Brief:** `HANDOVER-S58-AUDIT2.md`, executed as written.

**Full findings:** `HANDOVER-S59-AUDIT2-FINDINGS.md`.

---

## Session-tag verification (done twice, at start and before writing)

- `git log --oneline -1` → `26ba1e1  [S58] doc: book chapter -- pressure-testing the
  landscape (S53-S57)`.
- `ls SESSION-LOG/` → highest **logged** session is **S57**; `grep -c S58 PROGRESS.md`
  → **0**.
- So `S58` is **used** (commit tag, book-chapter session) but unlogged. **This session is
  S59.** Verified at start and again before any file was written. Tip unchanged
  throughout; working tree clean apart from this session's two new files.

---

## What was done

Three angles, per the handover. Every source fetched; nothing recalled.

### Angle 1 — the programme's own founding thesis (the main work)

Wrote out Route 1 and Route 2 exactly as `CLAUDE.md` and PRD v0.4 §§3.2–3.3 state them
(PRD text-extracted from the `.docx` and quoted), then determined where a regularity
argument for either would need an $L^\infty$ estimate of the same shape as
`lem:local-enstrophy-typeI`.

**Result: both inherit Wall 1, and by construction, not by accident.** Six independent
strands, all worked out in the findings document:

1. **Kinematic transfer.** Wall 1 is a property of $\omega\mapsto\nabla u$, i.e. of
   $\operatorname{div}u=0$ alone. `CLAUDE.md` makes exact enforcement of
   $\operatorname{div}u=0$ an *absolute* rule for both routes. So the recovery operator
   is identical, and `prop:s51-poloidal-sharp` — a **static, smooth, compactly
   supported divergence-free field**, no dynamics, no viscosity, no temperature — is a
   counterexample for Route 1 and Route 2 **verbatim, unmodified**. The v0.4 pivot's
   selling point ("Prize geometry from the first line") is exactly the property that
   guarantees the wall transfers.
2. **The maximum principle is on the demand side.** Free half: $T\ge T_{\min}\Rightarrow
   \mu\ge\mu_{\min}>0$, i.e. exactly what NS already grants. Useful half: $T\le T_{\max}$
   needs $\mu|\nabla u|^2\in L^{5/2+}$ (parabolic $L^\infty$ threshold $q>(n+2)/2$,
   re-derived by scaling), i.e. $\nabla u\in L^{5+}$ — and $\nabla u\in L^5$ already gives
   $u\in L^5_tL^{15/2}_x$, a Prodi–Serrin class ($2/5+2/5=4/5\le1$), hence regularity.
3. **The scalar's a priori content is the energy.** $\tfrac12\|u\|_{L^2}^2+\int T$ is
   conserved, so $\|T\|_{L^1}$ is bounded — and that is all. Boccardo–Gallouët for the
   $L^1$ source returns $T\in L^{q}(Q)$, $q<5/3$, which is *exactly* $u\in L^{10/3}(Q)$
   translated through $T\sim|u|^2$. **No gain.**
4. **Route 2 exactly.** At $\varepsilon=0$,
   $\int\theta(t)\,dx=\tfrac12(\|u_0\|_{L^2}^2-\|u(t)\|_{L^2}^2)$ — the auxiliary scalar's
   entire a priori content *is* the energy identity, which is precisely the structure
   Tao's averaged equation keeps.
5. **Extra terms, not fewer.** The Route 1 vorticity equation gains
   $(\nabla\mu\cdot\nabla)\omega$, $\nabla\mu\times\Delta u$ and
   $\varepsilon_{ijk}(\partial_j\partial_l\mu)(\partial_lu_k)$, requiring $L^\infty$
   bounds on $\nabla T$ and $\nabla^2T$. Count of $L^\infty$-of-recovered-quantity
   requirements goes 1 → 3.
6. **Both limits ARE the Prize problem.** At $\mu(T)\equiv\nu$ and at $\varepsilon=0$ the
   scalar decouples and the momentum equation is the Prize equation. So Phase 4 and Route
   2's $\varepsilon\to0$ — `CLAUDE.md` Priority Open Questions 2 and 3 — are restatements
   of the Millennium problem, not sub-problems.

**Two errors found in the authoritative PRD.** §3.5(c) has the sign of $d\!\int\!T/dt$
backwards ($\int T$ *increases*, at exactly the dissipation rate). §3.6's claim that LSU
theory for bounded smooth coefficients gives global regularity **proves too much**: its
$\mu\equiv\nu$ instance is the Millennium problem.

**Steelman done, not skipped.** Conjecture 3.5's scaling was worked out properly
(including the $A(T)\sim T^{3/2}$ feedback the PRD omits): self-consistent balance gives
$T\sim M^{4/(4+3\delta)}$, so at $\delta=1$ the core Reynolds number $\to0$ and the Type-I
ansatz is inconsistent — the conjecture is plausible *for Route 1*. It still does not
help: it needs a pointwise **lower** bound on $T$ (the direction no max principle gives),
its conclusion is Type-I exclusion (Wall 3, capped below the Prize), and the mechanism's
strength $M^{6\delta/(4+3\delta)}\to1$ as $\delta\to0$.

**Literature anchor.** Route 1 is a studied system (incompressible Navier–Stokes–Fourier
with temperature-dependent coefficients). Bulíček–Feireisl–Málek, *Nonlinear Anal. RWA*
10 (2009) 992–1015: large-data **weak** solutions in 3D. Bulíček–Málek–Shilkin,
*Nonlinear Anal. RWA* 19 (2014) 89–104: classical solutions only in **2D**, and only with
a **Ladyzhenskaya-type extra stress**.

### Angle 2 — stochastic / probabilistic representations

Constantin–Iyer `math/0511067v4` read **in full from the primary PDF**. Three findings:

- **The Leray–Hodge projection $\mathbb{P}$ is inside the representation formula**
  ($u=\mathbb{E}\,\mathbb{P}[(\nabla^tA)(u_0\circ A)]$, their Theorem 2.2 / Prop. 2.1).
  Same Calderón–Zygmund operator, same $L^\infty$ endpoint failure. The CZ step is
  relocated, not removed — so this is **not** in the class Wall 1(d) permits.
- **The stochastic Kelvin circulation theorem is pathwise**, and the authors write that
  its proof "is exactly the same as a proof showing circulation is conserved in inviscid
  flows" — no dissipative content. Euler-level identities cannot exclude blowup.
- **The authors state the limitation themselves**: known NS criteria "can be *translated*
  in criteria for the average of the stochastic flow map." Translated, not improved.
  Local existence is viscosity-independent (Euler-strength); the only global result in
  this line is Iyer's **small Reynolds number** theorem (`math/0702506`, Ann. IHP 2009),
  and the Le Jan–Sznitman/Bhattacharya-et-al cascade line is likewise **small data**.

The 2026 continuation `arXiv:2608.16915` (unrefereed, single unknown author — used only
for its self-declared negatives) states "**No claim is made regarding global
regularity**", proves with explicit counterexamples that
$\|\mathbb{E}\nabla A\|\not\Rightarrow\mathbb{E}\|\nabla A\|$ (Props 6.1, 6.3) — the
probabilistic avatar of Wall 1's own averaged-vs-sup failure — and its extracted criterion
(Thm 7.1) is **equivalent to the Serrin norm**. Fourth relabelling in this programme's
history (after S50, S52, S57).

Rough paths / regularity structures: Hairer's own text read from the primary PDF. Founding
hypothesis is **local subcriticality** — "at small scales, all nonlinear terms formally
disappear" — the exact negation of NS's difficulty, and **Lemma 8.10 proves it necessary**
(an iff). Output is local: "we do not claim that the solutions constructed here are
global" (Remark 1.18). Every NS-*regularity* application found modifies the equation:
Flandoli–Luo (PTRF 2021) needs added transport noise, "with **high probability**",
"sufficiently large noise intensity"; Flandoli–Hofmanová–Luo–Nilssen (AAP 2022) adds a
deterministic transport term and gets data "**outside arbitrary small sets**". Neither is
a universally quantified statement about the unperturbed equations.

### Angle 3 — kinetic and hydrodynamic limits

Clean negative, and the information flows the wrong way.

- Golse–Saint-Raymond, *J. Math. Pures Appl.* 91 (2009) 508–552, read from the primary
  PDF: the limit is governed by **Leray solutions**, and the authors write that "**the
  regularity of Leray solutions … is not known**", so the smooth branch is "limited to
  either local (in time) solutions, or to … small" data. Their Appendix C even calls the
  Leray projection's non-locality "one annoying difficulty" — the same CZ obstruction,
  inside the kinetic proof.
- Carrapatoso–Gallagher–Tristani `arXiv:2503.12046` (2025): "the Navier–Stokes equation
  can be solved in a **lower regularity setting** … than kinetic equations"; their
  Theorem 1 **assumes** an NS solution ($H^{1/2}$ data) and constructs kinetic solutions
  converging to it. NS → kinetic, not the reverse.
- Imbert–Silvestre, *EMS Surv. Math. Sci.* 7 (2020): Boltzmann's own global smoothness is
  "an outstanding open problem", and their conditional result holds "**provided that**"
  the mass, energy and entropy densities remain bounded — the macroscopic bounds are the
  **hypothesis**.
- The H-theorem's limit is the Leray energy inequality, which S56 §1.4(c) already closed
  by citation (Palasek's shell model satisfies the energy budget and blows up).

### Recency check

Nothing dated after S56 changes the picture. No BKM-logarithm progress, no new
ancient-Euler Liouville work, no Grujić version after S56's v3. Tao's blog has four
post-7-Sep posts in this area; the only fluid one (10 Sep, guest post on
Ganeshram–Duruisseaux–Anandkumar) is the item S56 already recorded, is about **Euler**,
and targets **blowup** — the AI asymmetry holds.

**Honest methodological caveat:** the arXiv API was rate-limited throughout this session
("Rate exceeded" on every attempt over ~40 minutes and two backoff strategies), so the
systematic metadata sweeps S56 ran were **not possible**. The recency check rests on
targeted web search plus direct primary fetches, and is correspondingly weaker evidence.
Recorded, not papered over.

---

## Verdict

**No new work package. All three angles reduce to Wall 1, or to "a different toolkit that
meets the same wall in its own vocabulary."** This is the **sixth** independent
confirmation after S48, S50, S52, S53, S57 — and the first that reaches the programme's
own founding thesis rather than an external candidate.

Options put to the owner (§5.4 of the findings): (1) stop adding to (A)/(B) and
consolidate; (2) accept S56 §3.3's different question, still not recommended; (3) stand
by for an external technique in Wall 1(d)'s permitted class, which the programme can now
price in one session, having done it twice.

Five failure modes named in advance (G1–G5), including two prohibitions: **G2 forbids a
Route 1/Route 2 numerical campaign as a response to this audit** (Wall 6 already
classifies this programme's numerics as unable to evidence an a priori estimate), and
**G3 forbids upgrading §1.8's formal scaling into `proofs/`**.

---

## Files written

- `HANDOVER-S59-AUDIT2-FINDINGS.md` (new)
- `SESSION-LOG/2026-09-13-S59-audit2-founding-thesis-stochastic-kinetic.md` (this file)
- `PROGRESS.md` (KEY FINDING block + Active-milestone line + session-log row)

Nothing committed, nothing pushed, per the handover.

## Tests

None. Literature and derivation session; no code touched.
