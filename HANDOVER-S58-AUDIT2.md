# Handover — second strategic audit: beyond what S56 already covered

Written 2026-09-13, immediately after S57 closed the Grujić thread. Owner's
direction: close that thread, then audit again for new methods toward
(A)/(B). **This audit must not simply repeat S56.** Read
`HANDOVER-S56-AUDIT-FINDINGS.md` and `PROGRESS.md`'s S56 and S57 KEY
FINDINGs in full first, so you know exactly what has already been searched
and found, and do not spend this session re-covering it.

## What S56 already closed, and is out of scope here unless something new and dated has appeared since

- Direct progress on weakening/removing the BKM/Calderón–Zygmund
  logarithm: searched exhaustively, nothing found.
- Recent (2024-2026) progress on ancient-Euler Liouville theorems beyond
  Seregin's four papers: searched, nothing found; the branch was shown to
  plausibly get *worse* given OpenAI's separate unforced Euler
  construction.
- Computer-assisted / interval-arithmetic profile exclusion: assessed,
  cannot deliver a universally-quantified regularity statement.
- Tao/Barker-style quantitative regularity criteria: assessed, alive but
  is the S52-verdict re-run (competing with specialists on their own tools).
- Energy-flux/cascade arguments: closed by citation (Palasek's shell model).
- Helicity-based mechanisms: no live literature; WP4's own negative stands.
- Grujić's log-weighted BMO mechanism: found at S56, priced and refuted at
  S57 (five/six/seventh confirmation Wall 1 is real).

**Only revisit any of the above if you find something specifically dated
after S56 (7-8 Sep 2026) that materially changes the picture** — the field
in this timeline has moved fast enough (three papers this thread engaged
with were revised within days of being found) that a short, targeted
recency check on each of the above is worthwhile, but should not consume
more than a fraction of this session.

## What this audit should actually spend its time on — three genuinely unexplored angles

### 1. The programme's own founding thesis — has it been properly re-examined?

`CLAUDE.md` states the **central thesis** of this entire programme is that
temperature $T$ is a scalar master variable, with $A(T)=k(T)/(\rho c_v)>0$
following from the second law, and that a **scalar** parabolic PDE admits a
**strong maximum principle unavailable to vector-based approaches** working
with $u$ or $\omega$ directly. Route 1 (incompressible NS with $\mu(T)$) and
Route 2 (exact Prize NS with an auxiliary scalar $\theta$) are named as the
primary and backup systems in the PRD.

**Every session from S18 onward has worked entirely inside the vector/
vorticity-alignment picture (Claim A, the pincer, $\sigma^*$, $\hat e$).**
Route 1 and Route 2 have not been actively pursued as a route to (A)/(B)
since the early (pre-S18) sessions, which were exploratory and numerical,
not aimed at a rigorous global-regularity proof.

**The question this audit must answer honestly**: does the scalar-$T$
picture's strong maximum principle actually supply anything Wall 1 lacks?
Wall 1 is fundamentally about recovering $\nabla u$ from $\omega$ via
Biot–Savart, an operation on a **vector** field with no maximum principle of
its own. Route 1's $\mu(T)$ enters the momentum equation as a coefficient,
not as a replacement for that recovery — so it is not obvious the scalar
maximum principle for $T$ ever touches the specific place Wall 1 bites. Work
this out rigorously rather than assuming either answer:

- Write out the Route 1 system precisely (as CLAUDE.md and the PRD state it)
  and identify exactly where, if anywhere, a regularity argument for it
  would need an estimate of the same shape as `lem:local-enstrophy-typeI`
  (an $L^\infty$ bound on a Biot–Savart-recovered quantity). If the answer
  is "nowhere, because $\mu(T)$ boundedness plus the standard NS estimates
  already close the argument," that would be a striking, high-value finding
  and should be checked with real care, not asserted. If the answer is "the
  same logarithm reappears when recovering $\nabla u$ from $\omega$
  regardless of $\mu(T)$," that is the expected, still-valuable negative
  finding — Route 1 inherits Wall 1 rather than avoiding it.
- Do the same for Route 2's auxiliary $\theta$ satisfying
  $\partial_t\theta+u\cdot\nabla\theta=\nu\Delta\theta+\nu|\nabla u|^2$ with
  $\mu_{\mathrm{eff}}=\nu+\varepsilon f(\theta)$. This equation's own
  right-hand side contains $|\nabla u|^2$ — check explicitly whether closing
  an estimate on $\theta$ requires exactly the same $L^\infty(\nabla u)$
  control Wall 1 obstructs.
- **Do not spend more than one focused pass on this.** If it clearly
  reduces to "Wall 1 all over again," say so precisely (which term, which
  estimate) and move on — this is meant to be a decisive, not open-ended,
  re-examination of the programme's own origin.

### 2. Probabilistic and stochastic representations of Navier–Stokes — genuinely untried by this programme

This programme has never engaged the stochastic-analysis literature on NS
regularity at all. Investigate, with real primary-source reading, not just
keyword collection:

- **Constantin–Iyer's stochastic Lagrangian (circulation) representation**
  of the Navier–Stokes equations and its use in regularity criteria (the
  original Constantin–Iyer papers and any 2023-2026 continuations or
  applications). Does this representation supply any estimate on $\nabla u$
  or $\omega$ that avoids the specific L^\infty Biot–Savart recovery Wall 1
  is built on? Read for the actual estimates, not just the existence of the
  representation.
- **Rough path theory / regularity structures** applied to fluid equations
  — search specifically for recent work applying these frameworks to 3D
  NS regularity (not just well-posedness of stochastically-forced NS, which
  is a different, larger literature and not directly relevant unless it
  bears on the deterministic regularity question).
- State plainly whether either genuinely offers a route around Wall 1's
  specific obstruction, or whether — as with the quantitative-regularity and
  energy-flux methodologies S56 already assessed — it is simply a different
  toolkit that runs into the same wall in its own vocabulary.

### 3. Kinetic and hydrodynamic-limit approaches — also untried

Boltzmann-to-Navier–Stokes hydrodynamic limits are a large, mature
literature (Golse and others). Investigate whether any regularity criterion
for the limiting NS equation has ever been derived *through* the kinetic
description rather than directly — i.e., whether kinetic-level structure
(entropy, the H-theorem, moment hierarchies) supplies any control on
$\nabla u$ that the direct PDE approach does not. This is a genuine
long shot; treat it as such, and report a clean negative if that is what a
real search finds rather than padding the report.

## Deliverable

Same shape as S56: for each of the three angles above, a clear finding
(something real and checkable found, or a clean, well-justified "no"), with
primary sources fetched and quoted, not recalled. End with one of:

- A concrete, scoped next work package, stated with named failure modes in
  advance, in the S46/S56 tradition.
- An honest statement that all three angles reduce to Wall 1 in different
  vocabulary, or offer nothing checkable, and that the programme's
  productive angles are genuinely exhausted for now.

**Do not manufacture a positive recommendation.** This programme's most
valuable sessions have been precise negative results.

## Standing discipline

- **Session-tag collision has happened five times running.** `git log
  --oneline -1` and `ls SESSION-LOG/` before writing anything; re-check
  before finalizing.
- Fetch every source. Do not recall S56's summary of what it found, and do
  not recall the contents of CLAUDE.md's Route 1/Route 2 description —
  re-read it directly.
- This is a **literature-and-strategy session**. Do not write proof-document
  content, numerics, or touch `book/`, `layer3/`, `layer4/`, `results/`. A
  second, concurrent session is writing a new book chapter — do not touch
  `book/` for any reason, including to check something; read the session
  logs and PROGRESS.md instead.
- Do not commit or push. Report findings and recommendation directly.
