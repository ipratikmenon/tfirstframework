# S59 — Second strategic audit: findings, worked calculation, recommendation

**Date:** 2026-09-13
**Tip at start and at finish:** `26ba1e1` (`[S58] doc: book chapter — pressure-testing the
landscape (S53–S57)`). **Session tag:** the book-chapter session took `S58` in its commit
message but wrote no `SESSION-LOG/` file and no `PROGRESS.md` entry, so the highest
*logged* session is S57 while `S58` is *used*. **This session is S59**, verified twice.
Working tree clean at start.
**Scope:** literature and mathematical analysis only. Nothing written to `proofs/`,
`book/`, `layer3/`, `layer4/`, `results/`. Nothing committed or pushed.
**Brief:** `HANDOVER-S58-AUDIT2.md`, read in full, executed as written.
**Read before writing:** `CLAUDE.md` in full; `HANDOVER-S58-AUDIT2.md` in full;
`HANDOVER-S56-AUDIT-FINDINGS.md` in full; the S56 and S57 KEY FINDING blocks in
`PROGRESS.md` (lines 986–1142); `lem:local-enstrophy-typeI` and `lem:BS-split`
(`proofs/claim_a_3d_proof_attempt.tex`, lines 5701–5735); `subsec:wall1-cz` (lines
15141–15213); `prop:s51-poloidal-sharp` (lines 11676–11705); PRD v0.4 §§3.1–3.7,
text-extracted from `tfirst_prd_v4.docx`.

---

## 0. Executive summary

1. **Angle 1 is a clean, decisive negative, and it is stronger than the handover
   anticipated.** Route 1 and Route 2 do not merely *inherit* Wall 1 — they inherit it
   **by construction**, because Wall 1 is a statement about the kinematic operator
   $\omega\mapsto\nabla u$, which is fixed by $\operatorname{div}u=0$ alone, and
   `CLAUDE.md` makes exact enforcement of $\operatorname{div}u=0$ an *absolute
   architecture rule* for both routes. The programme's own sharpness witness
   `prop:s51-poloidal-sharp` is a static, smooth, compactly supported divergence-free
   field: it contains no viscosity, no temperature and no dynamics, so it refutes
   $\|\nabla u\|_\infty\le K\|\omega\|_\infty$ for Route 1 and Route 2 *verbatim, with
   no adaptation whatsoever*.
2. **The scalar maximum principle sits on the wrong side of the ledger.** Its free half
   ($T\ge T_{\min}\Rightarrow\mu\ge\mu_{\min}>0$) delivers exactly the hypothesis the
   Prize equations already grant. Its useful half ($T\le T_{\max}$) requires, by the
   parabolic $L^\infty$ threshold, $\nabla u\in L^{5+}_{t,x}$ — and $\nabla u\in
   L^5_{t,x}$ already implies $u\in L^5_tL^{15/2}_x$, a Prodi–Serrin class, hence full
   regularity. **The maximum principle consumes exactly what Wall 1 fails to supply.**
3. **Route 2's auxiliary scalar carries no information beyond the energy identity.** At
   $\varepsilon=0$, $\int\theta(t)\,dx=\tfrac12(\|u_0\|_{L^2}^2-\|u(t)\|_{L^2}^2)$
   *exactly*. And Boccardo–Gallouët parabolic theory for the $L^1$ source returns
   $\theta\in L^{q}(Q)$ for $q<5/3$, which is *precisely* the Ladyzhenskaya–Prodi
   $u\in L^{10/3}(Q)$ bound translated through $\theta\sim|u|^2$. The scalar equation
   gives back the energy and nothing else.
4. **Both routes' final limits are the Prize problem, not perturbations of it.** At
   $\mu(T)\equiv\nu$ (Phase 4's endpoint) and at $\varepsilon=0$ (Route 2's endpoint),
   the scalar **decouples** and the momentum equation *is* the Prize equation.
   Uniformity in the limit parameter is therefore not "elementary" (PRD §3.2) and not a
   question about an LPS margin — **it is logically equivalent to the Prize theorem
   itself.** CLAUDE.md Priority Open Questions 2 and 3 are restatements of the whole
   problem.
5. **PRD §3.6 is false as stated**, and it is the load-bearing step of the whole
   architecture. It claims global regularity of $u$ follows from bounded smooth
   coefficients by LSU theory; the constant-coefficient instance $\mu\equiv\nu$ of that
   claim is the Millennium problem. The claim proves too much.
6. **Angle 2 is a clean negative with an unusually sharp form.** The Constantin–Iyer
   representation still contains the **Leray–Hodge projection $\mathbb{P}$** inside the
   formula $u=\mathbb{E}\,\mathbb{P}[(\nabla^tA)(u_0\circ A)]$ — the same
   Calderón–Zygmund operator, with the same $L^\infty$ endpoint failure. The stochastic
   Kelvin circulation theorem is pathwise and its proof "is exactly the same as a proof
   showing circulation is conserved in inviscid flows" (authors' words) — it therefore
   carries no dissipative information. Constantin and Iyer say themselves that known NS
   criteria "**can be translated in** criteria for the average of the stochastic flow
   map" — translation, not improvement. And the one 2026 continuation found states
   flatly "**No claim is made regarding global regularity**", proves the key exchange
   $\|\mathbb{E}\nabla A\|\not\Rightarrow\mathbb{E}\|\nabla A\|$ is **false** with
   explicit counterexamples, and its central criterion (Thm 7.1) is *equivalent to the
   Serrin norm*.
7. **Rough paths / regularity structures: structurally inapplicable, in Hairer's own
   terms.** The theory's founding hypothesis is **local subcriticality** — "at small
   scales, all nonlinear terms formally disappear" — which is the exact negation of 3D
   NS's difficulty; Hairer proves (Lemma 8.10) that this hypothesis is *necessary*, not
   just convenient; and he states "we do not claim that the solutions constructed here
   are global." Everything found applying these tools to NS regularity works on a
   **modified** system (added transport noise, added deterministic drift) and yields
   either high-probability or off-a-small-set statements — never a universally
   quantified statement about the unperturbed equations.
8. **Angle 3 is a clean negative and the information flows the wrong way.** The
   Boltzmann→NS limit produces **Leray solutions** (Golse–Saint-Raymond, verbatim), the
   authors state that "the regularity of Leray solutions … is not known", the 2025 state
   of the art *assumes* an NS solution as input and constructs kinetic solutions
   converging to it, and the Boltzmann equation's *own* global regularity is open, with
   its best conditional results (Imbert–Silvestre) **hypothesising the macroscopic
   bounds**. The kinetic level is strictly downstream.
9. **Recency check: nothing new.** No BKM-logarithm progress after S56. No new
   ancient-Euler Liouville work. Tao's blog has four posts after 7 Sep 2026 relevant to
   this area; the only fluid one (10 Sep, guest post on Ganeshram–Duruisseaux–
   Anandkumar) is the item S56 already recorded via the 7 Sep post's edit, is about
   **Euler**, and targets **blowup**.
10. **Recommendation: no new work package. All three angles reduce to Wall 1 or to
    "different toolkit, same wall."** This is the **sixth** independent confirmation
    after S48, S50, S52, S53, S57 — and the first that reaches back to the programme's
    own founding thesis. §5 states what that implies and offers the owner two honest
    options, neither of them a manufactured target.

---

## 1. Angle 1 — the founding thesis, worked out

The handover asked for the calculation, not an opinion. Here it is. Throughout, set
$\rho=1$ and $c_v=1$ (both constants in Route 1 by fiat), and work on $\mathbb{T}^3$ so
every divergence integrates to zero. I write $Q$ for the space-time domain.

### 1.1 The two systems, exactly as `CLAUDE.md` and PRD v0.4 §3.2–§3.3 state them

**Route 1** (PRD §3.2, verbatim from the extracted `.docx`: "Momentum equation
(incompressible, variable viscosity): `rho·(del_t u + (u·del)u) = −del p + del·(mu(T)
del u)     div u = 0`; Temperature equation (scalar master PDE): `rho·cv·del_t T =
del·(k(T) del T) − rho·cv·(u·del)T + mu(T)|del u|^2`"):

$$
\partial_t u+(u\cdot\nabla)u=-\nabla p+\nabla\!\cdot\!\big(\mu(T)\nabla u\big),
\qquad \operatorname{div}u=0, \tag{R1a}
$$
$$
\partial_t T+u\cdot\nabla T=\nabla\!\cdot\!\big(k(T)\nabla T\big)+\mu(T)|\nabla u|^2,
\tag{R1b}
$$
$$
\mu(T)=\mu_0(T/T_0)^{3/2},\quad k(T)=k_0(T/T_0)^{3/2},\quad A(T)=k(T)/(\rho c_v)>0.
$$

**Route 2** (PRD §3.3, verbatim: "`del_t theta + u·del theta = nu·Delta theta +
nu·|del u|^2     theta(x,0) = 0`", "`mu_eff = nu + epsilon·f(theta)`"):

$$
\partial_t u+(u\cdot\nabla)u=-\nabla p+\nabla\!\cdot\!\big(\mu_{\rm eff}\nabla u\big),
\qquad \operatorname{div}u=0,\qquad \mu_{\rm eff}=\nu+\varepsilon f(\theta), \tag{R2a}
$$
$$
\partial_t\theta+u\cdot\nabla\theta=\nu\Delta\theta+\nu|\nabla u|^2,\qquad
\theta(\cdot,0)=0. \tag{R2b}
$$

`CLAUDE.md`'s architecture rules make one property **absolute** for both:

> "**Route 1**: ρ = const and div u = 0 must be enforced **exactly** (not
> approximately)"

That rule is the whole finding of §1.2.

### 1.2 Where Wall 1 bites — and it bites in exactly the same place, for a kinematic reason

Wall 1 (`subsec:wall1-cz`, item (a)) is:

> "Every enstrophy-budget estimate … controls the vortex stretching by a sup-norm of
> $\nabla u$, recovered from $\omega$ by Biot–Savart. Closing any of them requires that
> sup-norm to cost $O(M)$."

and item (b) locates the failure at `lem:local-enstrophy-typeI`, whose statement is

> "the vorticity Caccioppoli inequality … with
> $\|\nabla u\|_{L^\infty}\lesssim M\log(e+\|u\|_{H^3}/M)$ gives
> $\nu\iint_{Q(2\rho)}|\nabla\omega|^2\lesssim M^2\rho^3(1+M\rho^2\mathcal L)\lesssim
> M^2\rho^3\mathcal L$."

**Step 1. The recovery operator.** For any divergence-free $u$ on $\mathbb{T}^3$ with
$\omega=\operatorname{curl}u$,
$$
u=\operatorname{curl}(-\Delta)^{-1}\omega,\qquad
\nabla u=\underbrace{\nabla\operatorname{curl}(-\Delta)^{-1}}_{=:K}\omega ,
$$
$K$ being a fixed matrix of second-order Riesz-transform compositions
$R_iR_j$. **$K$ depends on nothing but $\operatorname{div}u=0$.** It does not see the
momentum equation, the viscosity, the pressure, the temperature, or the time variable.

**Step 2. Route 1 and Route 2 use the identical $K$.** Both enforce
$\operatorname{div}u=0$ *exactly* — by CLAUDE.md's own absolute rule, and by PRD §3.2's
"div u = 0 holds exactly throughout". Therefore the map $\omega\mapsto\nabla u$ in
(R1a)/(R2a) is *the same operator, symbol for symbol*, as in the Prize equations. There
is no $\mu$, no $T$, no $\theta$ and no $\varepsilon$ anywhere in it.

**Step 3. The programme's own counterexample transfers with zero modification.**
`prop:s51-poloidal-sharp` produces, for each $N$, "a smooth, compactly supported,
divergence-free *axisymmetric* field $u^{(N)}$, supported in the ball $B_{r^*}$" with
$$
\|\nabla u^{(N)}\|_{L^\infty}\ge cN\|\omega^{(N)}\|_{L^\infty},\qquad
\|u^{(N)}\|_{H^3}\le C(r^*)2^{2N}\|\omega^{(N)}\|_{L^\infty}.
$$
This is a **static kinematic object**. It is not a solution of anything; it involves no
viscosity, no temperature and no time. It is therefore a counterexample to
$\|\nabla u\|_\infty\le K\|\omega\|_\infty$ **for the Route 1 and Route 2 velocity
fields verbatim**, because those fields are divergence-free.

> **Conclusion (angle 1, part (i)).** Route 1 and Route 2 do not inherit Wall 1 as an
> accident of technique. They inherit it because Wall 1 is a property of
> $\operatorname{div}u=0$, and both routes were *designed* to enforce
> $\operatorname{div}u=0$ exactly. **The v0.4 pivot's central selling point — "Prize
> geometry from the first line" — is precisely the property that guarantees Wall 1
> transfers with no adaptation.** (Route 3, compressible NSF, would at least have a
> *different* kinematic operator, since $\operatorname{div}u\neq0$ there — a worse one,
> but a different one. The pivot chose the side that inherits the wall exactly.)

**Step 4. The vorticity budget, term by term, confirms the location is identical.**
Taking $\operatorname{curl}$ of (R1a), using
$\big(\nabla\!\cdot(\mu\nabla u)\big)_i=\mu\Delta u_i+(\nabla\mu\cdot\nabla)u_i$ and
$\nabla\times(fV)=f\nabla\times V+\nabla f\times V$:
$$
\partial_t\omega+(u\cdot\nabla)\omega
=\underbrace{(\omega\cdot\nabla)u}_{\text{(S)}}
+\mu\Delta\omega
+\underbrace{(\nabla\mu\cdot\nabla)\omega}_{\text{(E1)}}
+\underbrace{\nabla\mu\times\Delta u}_{\text{(E2)}}
+\underbrace{\big[\varepsilon_{ijk}(\partial_j\partial_l\mu)(\partial_lu_k)\big]_i}_{\text{(E3)}} .
\tag{R1-vort}
$$
(Route 2 is the same with $\mu\to\mu_{\rm eff}$.) Term **(S)** is *literally the term
from the Prize equations*, unmodified — $\mu$ multiplies the diffusion, never the
stretching. So the local enstrophy Caccioppoli for (R1-vort) is
$$
\tfrac12\partial_t\!\!\int\!|\omega|^2\eta^2+\mu_{\min}\!\!\iint\!|\nabla\omega|^2\eta^2
\;\le\;\|\nabla u\|_{L^\infty(\operatorname{supp}\eta)}\!\iint\!|\omega|^2\eta^2
\;+\;(\text{E1})+(\text{E2})+(\text{E3})+\text{cutoff},
$$
and the first right-hand term is *exactly* the term
`lem:local-enstrophy-typeI` pays $\mathcal{L}=\log(e+\|u\|_{H^3}/M)$ for. **Same term,
same operator, same logarithm, same line of the argument.** Nothing about $\mu(T)$
touches it.

**Step 5. The new terms make it strictly worse, not better.** Using
$\Delta u=-\nabla\times\omega$ (valid because $\operatorname{div}u=0$):

- **(E1)** and **(E2)** are absorbable into $\mu_{\min}\iint|\nabla\omega|^2\eta^2$ by
  Young, at cost $C\mu_{\min}^{-1}\|\nabla\mu\|_{L^\infty}^2\iint|\omega|^2\eta^2$.
  This needs $\|\nabla\mu\|_{L^\infty}=\|\mu'(T)\nabla T\|_{L^\infty}$, i.e.
  **an $L^\infty$ bound on $\nabla T$.**
- **(E3)** needs $\nabla^2\mu=\mu''|\nabla T|^2+\mu'\nabla^2T$ paired against
  $\nabla u$ and $\omega$: **an $L^\infty$-type bound on $\nabla^2 T$**, against a first
  derivative of $u$.

Neither requirement exists in the Prize problem. Route 1 therefore requires the
inherited $L^\infty(\nabla u)$ estimate **plus two new $L^\infty$-type estimates on
derivatives of the master scalar**, each of which (§1.3) itself needs $L^\infty(\nabla
u)$-strength control. The count of $L^\infty$-of-a-recovered-quantity requirements goes
**up**, from one to three.

### 1.3 What the strong maximum principle actually delivers

`CLAUDE.md`'s central claim is that "Because T satisfies a quasilinear parabolic
**scalar** PDE, the **strong maximum principle** applies — a tool unavailable to
vector-based approaches." Both halves of that principle are available. Here is exactly
what each buys.

**Free half (minimum principle).** In (R1b) the source $\mu(T)|\nabla u|^2\ge0$, so
$T\ge\min T_0=T_{\min}>0$ for all time, hence
$$
\mu(T)\;\ge\;\mu_0(T_{\min}/T_0)^{3/2}=:\mu_{\min}>0 .
$$
This is unconditional and correct. **What it delivers is a strictly positive lower bound
on the viscosity — which is exactly the hypothesis the Prize equations already grant
($\nu>0$).** Zero new information. Identically for Route 2: $\theta\ge0$, hence
$\mu_{\rm eff}\ge\nu$, which is true by inspection.

**Useful half (maximum principle).** The maximum principle for (R1b) gives
$$
\frac{d}{dt}\max_x T\;\le\;\max_x\big(\mu(T)|\nabla u|^2\big)
\;\le\;\mu_{\max}\|\nabla u(t)\|_{L^\infty}^2,
\qquad\text{so}\qquad
T_{\max}(t)\le T_{\max}(0)+\mu_{\max}\!\!\int_0^t\!\|\nabla u\|_{L^\infty}^2\,ds .
$$
The input is $\nabla u\in L^2_tL^\infty_x$ — a BKM/Prodi–Serrin-endpoint quantity. **The
maximum principle's input is the very object Wall 1 obstructs.**

**Can the diffusion do better than the bare maximum principle?** Yes, a little, and this
is the sharpest way to price it. For a scalar parabolic equation
$\partial_tT+u\cdot\nabla T-\nabla\!\cdot\!(A\nabla T)=S$ in $n$ space dimensions,
boundedness of $T$ requires $S\in L^q$ with $q>(n+2)/2$; for $n=3$, $q>5/2$. (Standard
Moser/LSU–Aronson–Serrin theory; and it is forced by scaling: under
$T_\lambda(x,t)=T(\lambda x,\lambda^2t)$, $S_\lambda=\lambda^2S(\lambda x,\lambda^2 t)$,
one has $\|S_\lambda\|_{L^q(Q_1)}=\lambda^{2-(n+2)/q}\|S\|_{L^q(Q_\lambda)}$ while
$\|T_\lambda\|_\infty=\|T\|_\infty$, so the estimate is scale-consistent exactly at
$q=(n+2)/2$ and subcritical above it.)

So $T\in L^\infty$ needs $\mu|\nabla u|^2\in L^{5/2+}(Q)$, i.e.
$$
\boxed{\ \nabla u\in L^{5+\epsilon}_{t,x}. \ }
$$

**And that hypothesis already gives regularity outright.** On $\mathbb{T}^3$ (or
$\mathbb{R}^3$), $W^{1,5}\hookrightarrow L^{15/2}$, so $\nabla u\in L^5_{t,x}$ gives
$u\in L^5_tL^{15/2}_x$, and
$$
\frac{2}{5}+\frac{3}{15/2}=\frac25+\frac25=\frac45\;\le\;1 ,
$$
which is a **Prodi–Serrin–Ladyzhenskaya class with room to spare** ($q=15/2>3$). Hence:

> **Conclusion (angle 1, part (ii)).** Phase 1's deliverable ($T\in[T_{\min},T_{\max}]$,
> hence $\mu$ uniformly bounded above) is **not a step toward regularity**. It is a
> *consequence* of regularity, obtainable only from a hypothesis strictly stronger than
> one that already closes the problem. The strong maximum principle is on the **demand**
> side of the Wall 1 ledger, not the supply side.

### 1.4 The scaling ledger: $T$ and $\theta$ scale like $|u|^2$, and their only a priori bound *is* the energy

With coefficients frozen (constant $\mu,k$ — which is exactly Route 2 at $\varepsilon=0$,
and Route 1 up to the symmetry-breaking constant $T_0$), the pair (R1a)–(R1b) is
**exactly invariant** under the Navier–Stokes scaling
$$
u\mapsto\lambda u(\lambda x,\lambda^2t),\qquad
p\mapsto\lambda^2p(\lambda x,\lambda^2t),\qquad
T\mapsto\lambda^2T(\lambda x,\lambda^2t).
$$
(Check: every term of (R1b) scales as $\lambda^4$, including $\mu|\nabla u|^2$.) So
**$T$ has the scaling dimension of $|u|^2$.** Consequences, one-to-one with the velocity
ledger:

| velocity | scalar ($T$ or $\theta$) | status |
|---|---|---|
| critical $u\in L^\infty_tL^3_x$ | critical $T\in L^\infty_tL^{3/2}_x$ | what regularity needs |
| energy $u\in L^\infty_tL^2_x$ | **$T\in L^\infty_tL^1_x$** | what is available |
| $u\in L^{10/3}(Q)$ | $T\in L^{5/3}(Q)$ | Ladyzhenskaya–Prodi |
| $\|u\|_{L^\infty}$ supercritical ($\lambda^1$) | $\|T\|_{L^\infty}$ supercritical ($\lambda^2$) | what the max principle would bound |

And the available bound is *exactly* the energy, not an analogue of it. Integrating
(R1b) over $\mathbb{T}^3$ (both the diffusion and the transport are divergences):
$$
\frac{d}{dt}\int T\,dx=\int\mu(T)|\nabla u|^2dx,
\qquad
\frac{d}{dt}\,\tfrac12\!\int|u|^2dx=-\int\mu(T)|\nabla u|^2dx ,
$$
so $\tfrac12\|u(t)\|_{L^2}^2+\int T(t)\,dx$ is **conserved**, giving
$$
\|T(t)\|_{L^1}\;\le\;\tfrac12\|u_0\|_{L^2}^2+\|T_0\|_{L^1}.
$$

Two remarks that matter.

- **PRD §3.5(c) is backwards.** It asserts "Total thermal energy $E_\theta=\rho c_v\int
  T\,dx$ is dissipated (not created) by the system." The computation above shows
  $\int T$ is **monotonically increasing**, at exactly the rate of viscous dissipation.
  The *conclusion* PRD wants ($\int T$ bounded) is nonetheless true, by the conservation
  law. What is a non sequitur is the next sentence, "Therefore $T_{\max}(t)$ is
  controlled by initial energy $E_0$": an $L^1$ bound does not control $L^\infty$, and
  the gap between them is the whole difficulty.
- **The parabolic $L^1$ theory returns exactly the energy and nothing more.** With
  source only in $L^1(Q)$, Boccardo–Gallouët theory gives $T\in L^q(Q)$ for
  $q<(n+2)/n=5/3$ and $\nabla T\in L^r(Q)$ for $r<(n+2)/(n+1)=5/4$. And $5/3$ is
  *precisely* the exponent one gets by translating the Ladyzhenskaya–Prodi bound
  $u\in L^{10/3}(Q)$ through $T\sim|u|^2$. **The scalar master equation's own regularity
  theory reproduces the velocity energy estimate exactly, with no gain.**

This is why the maximum principle cannot help: $\|T\|_\infty$ is a *supercritical*
quantity (dimension $\lambda^2$) and the data are *subcritical* ($L^1$, dimension
$\lambda^{-1}$). A maximum principle is a tool for converting a controlled source into
an $L^\infty$ bound; it does not change scaling dimensions, and here the source is not
controlled in any norm above $L^1$.

### 1.5 Route 2's $\theta$: its entire a priori content is the energy identity

Route 2 makes this exact rather than analogous. Take $\varepsilon=0$ in (R2a)–(R2b)
(the case the whole route is a perturbation of). Then $\int\nu\Delta\theta=0$ and
$\int u\cdot\nabla\theta=0$, so
$$
\frac{d}{dt}\int\theta\,dx=\nu\!\int|\nabla u|^2dx=-\frac{d}{dt}\,\tfrac12\!\int|u|^2dx,
$$
and with $\theta(\cdot,0)=0$,
$$
\boxed{\ \int_{\mathbb{T}^3}\theta(x,t)\,dx\;=\;\tfrac12\big(\|u_0\|_{L^2}^2-\|u(t)\|_{L^2}^2\big)\ }
$$
— **identically the dissipated energy.** For $\varepsilon>0$ the same computation gives
$\int\theta(t)=\nu\iint|\nabla u|^2\le\iint\mu_{\rm eff}|\nabla u|^2=
\tfrac12(\|u_0\|^2-\|u(t)\|^2)\le\tfrac12\|u_0\|_{L^2}^2$.

So the auxiliary scalar is not new information encoded in a new variable; it is the
energy identity written as a PDE. PRD §3.3's claim that "$\theta$ encodes thermodynamic
heating from strain — it is the finer structure Tao's averaging removes" is exactly
wrong in the operative sense: **$\theta$'s a priori content is the energy budget, which
is precisely the structure Tao's averaged equation *keeps***. Tao's construction respects
the energy identity; that is why it is a counterexample. A variable whose only a priori
bound is the energy cannot distinguish NS from Tao's averaged NS. (S56 §1.4(c) recorded
the same verdict for energy-flux/cascade arguments, by citation to Palasek's shell
model; this is the same verdict reached internally, by identity rather than by
citation.)

And to *use* $\theta$ one needs $f(\theta)$ **large where the strain concentrates** —
i.e. a **pointwise lower bound** on $\theta$ at the concentration. A maximum principle
gives upper bounds from bounded sources and lower bounds from bounded-below data; the
lower bound it gives here is $\theta\ge0$, which is the trivial one. Getting a
*nontrivial* pointwise lower bound at a moving concentration requires a parabolic
Harnack inequality for $\partial_t+u\cdot\nabla-\nu\Delta$ with a critical-size drift —
available only in the divergence-free $L^\infty_tBMO^{-1}$ class, at a Type-I-borderline
constant, and its conclusion would be a *Type-I exclusion*, i.e. Wall 3 territory, capped
below the Prize.

### 1.6 PRD §3.6 is false as stated, and it is the load-bearing step

PRD §3.6, verbatim:

> "With T in [T_min, T_max] for all time, mu(T) and k(T) are uniformly bounded: 0 <
> mu_min <= mu(T) <= mu_max < infinity. … This is NS with bounded, smooth, T-dependent
> coefficients. By Ladyzhenskaya-Uraltseva-Solonnikov theory for parabolic systems with
> bounded smooth coefficients, global regularity of u follows from T regularity. In
> particular u satisfies the LPS condition u in Lp(Lq) at the critical exponent, closing
> the gap."

**This proves too much.** The class "NS with bounded, smooth coefficients" contains the
*constant* coefficient $\mu\equiv\nu$. For that member the system is the Prize equations
verbatim, and the claimed implication would settle the Millennium problem. LSU theory
gives regularity for *linear* parabolic problems (and for quasilinear scalar problems
under structure conditions); it does not, and cannot, absorb the quadratic convective
nonlinearity $(u\cdot\nabla)u$ in three dimensions. Granting Phase 1 *in full and for
free* returns the enstrophy budget
$$
\tfrac12\tfrac{d}{dt}\|\omega\|_{L^2}^2+\mu_{\min}\|\nabla\omega\|_{L^2}^2
\;\le\;\int\omega\cdot\nabla u\cdot\omega\;+\;C(\|\nabla\mu\|_\infty)\|\omega\|\|\nabla\omega\| ,
$$
whose first term is the classical supercritical term and whose closure requires exactly
$\|\nabla u\|_\infty\lesssim M$ — Wall 1. **Phase 1 → Phase 3 is not a step; it is the
whole problem, plus two extra terms.**

### 1.7 The two limits (Phase 4; Route 2's $\varepsilon\to0$) *are* the Prize problem

`CLAUDE.md` Priority Open Questions 2 and 3 ask whether the LPS margin
$\varepsilon(\mu)$ survives $\mu(T)\to\nu$, and whether the $\theta$-system converges
with the margin inherited as $\varepsilon\to0$. PRD §3.2 calls the first limit "simple"
and "elementary."

Interpolate Route 1 by $\delta\in[0,1]$:
$\mu_\delta(T)=\nu(T/T_0)^{3\delta/2}$, $k_\delta(T)=k_0(T/T_0)^{3\delta/2}$.

- At $\delta=1$ this is Route 1.
- At $\delta=0$, $\mu\equiv\nu$ identically, so (R1a) **is the Prize equation**, and
  (R1b) becomes a *passive scalar* advection–diffusion equation driven by $u$, with **no
  back-coupling whatsoever**.

Identically for Route 2 at $\varepsilon=0$: $\mu_{\rm eff}\equiv\nu$, (R2a) is the Prize
equation, and $\theta$ is a slaved passive scalar.

> **Conclusion (angle 1, part (iii)).** The limit system is not a degenerate
> approximation of the Prize equations — **it is the Prize equations**, with an inert
> scalar attached. Therefore "does the estimate survive the limit?" is not a question
> about a margin: an estimate for the family that is *uniform down to the endpoint* is,
> at the endpoint, an estimate for the Prize equations. **Priority Open Questions 2 and 3
> are restatements of the Millennium problem**, not sub-problems of it, and Phase 4 is
> not elementary.

This is the same structural defect as the well-studied Ladyzhenskaya/Smagorinsky
regularisation family (extra stress $\propto|\nabla u|^{p-2}\nabla u$), which yields
global smooth solutions for each fixed coefficient and loses them uniformly in the limit.
It is also, in miniature, the reason Wall 9 forbids reading forced constructions as
evidence about the unforced equations: the parameter that carries the mechanism is the
parameter the Prize problem sets to zero.

### 1.8 Steelman: Conjecture 3.5's scaling argument — and why it still does not help

The one place Route 1 has a genuine idea is PRD §3.7 / Conjecture 3.5: near a
concentration, $|\nabla u|^2$ is enormous, so $T$ heats, so $\mu(T)$ grows, so
dissipation wins. I ran this properly rather than dismissing it, because if it worked at
the level the handover flagged it would be a real finding.

Under a Type-I ansatz ($\tau=T_*-t$, $M\sim\tau^{-1}$, core radius $r^*=M^{-1/2}\sim
\tau^{1/2}$), suppose $T\sim M^a$ at the core. Then $\mu,A\sim M^{3\delta a/2}$; the heat
generated in the core per unit time is $\sim\mu M^2(r^*)^3=\mu M^{1/2}$; the diffusion
length over the remaining time is $\ell_D\sim(A\tau)^{1/2}\sim M^{3\delta a/4-1/2}$.
Because $\ell_D/r^*\sim M^{3\delta a/4}\to\infty$ for $\delta a>0$, **the heat leaves the
core faster than it accumulates**, so the correct balance spreads it over $\ell_D^3$:
$$
M^{a}\sim\frac{\mu M^{-1/2}}{\ell_D^{3}}=M^{1-3\delta a/4}
\quad\Longrightarrow\quad
a=\frac{4}{4+3\delta},\qquad
\mu_{\rm core}\sim M^{6\delta/(4+3\delta)} .
$$
At $\delta=1$: $T\sim M^{4/7}$, $\mu_{\rm core}\sim M^{6/7}$, local Reynolds number
$\mathrm{Re}_{\rm loc}=|u|r^*/\mu\sim\mu^{-1}\to0$ — the core is asymptotically Stokes and
the ansatz is inconsistent. **So Conjecture 3.5's *conclusion* is plausible at the level
of formal scaling, and the PRD's own exponent bookkeeping (which ignores that $A(T)$
blows up too, and so overestimates $T$ as $\sim M$ rather than $M^{4/7}$) is too
generous but not fatally so.** Three reasons it still does not help:

1. It requires a **pointwise lower bound** on $T$ at the core (§1.5) — the direction no
   maximum principle supplies, and one whose only known route is a parabolic Harnack
   inequality for a *degenerate quasilinear* operator ($A(T)\sim T^{3/2}$ blowing up)
   with critical drift. Nothing in the literature covers that combination.
2. Its conclusion, if obtained, is "**no Type-I singularity for Route 1**" — Wall 3
   territory, which `subsec:wall1-cz`'s neighbours already classify as capped below the
   Prize, and which in the one class where this programme has traction is
   Seregin–Šverák 2009.
3. The mechanism's strength is $M^{6\delta/(4+3\delta)}$, which is $\to1$ (no
   enhancement) as $\delta\to0$ at every fixed $M$. **The mechanism *is* the
   $\mu$-variation**, and Phase 4 is defined as removing it. This is §1.7 again, now
   with an exponent attached.

*(§1.8's scaling is formal — a self-consistent balance, not a theorem. Its qualitative
conclusion in item 3 does not depend on the details: at $\delta=0$ the coupling is
identically absent, which is exact.)*

### 1.9 The Route 1 system is a studied object, and its 3D state of the art is weak solutions

Route 1 is not a new system. It is the incompressible Navier–Stokes–Fourier system with
temperature-dependent material coefficients and viscous heating. Primary references
found and checked:

- **M. Bulíček, E. Feireisl, J. Málek**, *A Navier–Stokes–Fourier system for
  incompressible fluids with temperature dependent material coefficients*, **Nonlinear
  Anal. Real World Appl. 10 (2009), no. 2, 992–1015**, DOI
  `10.1016/j.nonrwa.2007.11.018`. Established: "long-time and large-data existence of
  **weak** solutions for the incompressible three-dimensional Navier–Stokes–Fourier
  system with the viscosity and the heat conductivity depending on the temperature
  (internal energy)." Weak solutions, in 3D, large data — and the temperature equation is
  handled through an entropy/global-balance formulation precisely because its right-hand
  side is only $L^1$.
- **M. Bulíček, J. Málek, T. N. Shilkin**, *On the Regularity of Two-Dimensional Unsteady
  Flows of Heat-Conducting Generalized Newtonian Fluids*, **Nonlinear Anal. Real World
  Appl. 19 (2014), 89–104**, DOI `10.1016/j.nonrwa.2014.03.003`: classical solvability
  for heat-conducting incompressible flows is obtained **in two dimensions**, and even
  there only "under certain structural assumptions on the Cauchy stress that include
  generalizations of the **Ladyzhenskaya** or power-law like models."

So: for the exact system Route 1 names, the specialists get *weak* solutions in 3D, and
*classical* solutions only in 2D and only with a Ladyzhenskaya-type extra stress. This
is consistent with everything in §§1.2–1.7 and inconsistent with PRD §3.6.

### 1.10 Angle 1 — verdict

**Expected negative, obtained, in a sharper form than expected.**

- Wall 1 transfers **verbatim**, for a kinematic reason that is a *consequence of the
  v0.4 pivot's own design goal*, and the programme's own sharpness witness is a
  counterexample for both routes without modification.
- The scalar maximum principle's free half returns what NS already has; its useful half
  requires a hypothesis strictly stronger than one that already implies regularity.
- The scalar's only a priori bound is the energy — exactly, as an identity, for Route 2.
- Both routes add $L^\infty$-type requirements on derivatives of the scalar that the
  Prize problem does not have.
- Both routes' limits *are* the Prize problem, so their Phase 4 / $\varepsilon\to0$
  questions are restatements of it.

**`CLAUDE.md`'s central-thesis sentence** — "Because T satisfies a quasilinear parabolic
**scalar** PDE, the **strong maximum principle** applies — a tool unavailable to
vector-based approaches working with u or ω directly" — **is true and irrelevant.** The
tool is genuinely unavailable to vector methods and genuinely available here; it simply
does not act at the place Wall 1 bites, and the quantity it would bound ($\|T\|_\infty$)
is supercritical, so it cannot be reached from the subcritical data available.

---

## 2. Angle 2 — stochastic and probabilistic representations

Everything below was fetched and read. Where I read only an abstract I say so.

### 2.1 Constantin–Iyer: the Leray projection is still in the formula

**P. Constantin, G. Iyer**, *A stochastic Lagrangian representation of the
3-dimensional incompressible Navier–Stokes equations*, `arXiv:math/0511067v4`
(3 Nov 2005, final 31 Aug 2006), **Comm. Pure Appl. Math. 61 (2008), 330–345**. Read
in full from the primary PDF.

Theorem 2.2, verbatim, is the system
$$
dX=u\,dt+\sqrt{2\nu}\,dW,\qquad A=X^{-1},\qquad
u=\mathbb{E}\,\mathbb{P}\big[(\nabla^{t}A)(u_0\circ A)\big],\qquad X(a,0)=a,
$$
with, from Proposition 2.1, "**Here $\mathbb{P}$ is the Leray–Hodge projection … on
divergence free vector fields.**"

**This is the finding.** $\mathbb{P}=I+\nabla(-\Delta)^{-1}\operatorname{div}$ is a matrix
of Riesz transforms — a Calderón–Zygmund operator of order zero, unbounded on $L^\infty$,
with exactly the endpoint failure Wall 1 is built on. The representation does not remove
the CZ recovery step; it **relocates** it, from "recover $\nabla u$ from $\omega$ by
Biot–Savart" to "recover $u$ from the stochastic Weber field by Leray projection", and
then one must still differentiate. Wall 1's Consequence clause permits "any reformulation
in which $\nabla u$ is **not** recovered by an $L^\infty$ Calderón–Zygmund bound"; this
reformulation is not one.

**The authors state the limitation themselves.** §2, verbatim:

> "Although it is evident from equations (2.9) and (2.13), we explicitly point out that
> the source of growth in the velocity and vorticity fields arises from the gradient of
> the noisy flow map $X$. The Beale-Kato-Majda [1] criterion guarantees if the vorticity
> $\omega$ stays bounded, then no blow up can occur in the Euler equations. In the case
> of the Navier-Stokes equations, well known criteria for regularity exist and they **can
> be translated in criteria for the average of the stochastic flow map**."

Translated, not improved. That is the authors' own assessment, and it is exactly the
"different toolkit, same wall" verdict S56 reached for quantitative regularity.

**The stochastic Kelvin circulation theorem is pathwise, hence inviscid in content.**
Proposition 2.9: with $\tilde u=\mathbb{P}[(\nabla^tA)(u_0\circ A)]$ (the *un-averaged*
stochastic velocity, $u=\mathbb{E}\tilde u$),
$$
\oint_{X(\Gamma)}\tilde u\cdot dr=\oint_\Gamma u_0\cdot dr ,
$$
and the authors add:

> "We remark that the above proof is **exactly the same as a proof showing circulation is
> conserved in inviscid flows**."

So the circulation identity holds realisation by realisation and carries **no dissipative
information**: viscosity enters only through $\mathbb{E}$. Any estimate derived from the
circulation theorem alone is an Euler-level estimate — and Euler now has a claimed
unforced finite-time singularity from $C_c^\infty$ data (OpenAI, recorded at S56), so
Euler-level identities cannot exclude blowup. To extract dissipation one must exchange
$\mathbb{E}$ with the nonlinearity, which is precisely where the difficulty sits (§2.2).

**What the method actually delivers.** Theorem 2.6 is *local* existence with an existence
time "**independent of viscosity**" — and the authors note "The theorem and proof also
work when the viscosity $\nu=0$", i.e. it is an Euler-strength local theory. The one
global result in this line is **G. Iyer**, *A stochastic Lagrangian proof of global
existence of the Navier-Stokes equations for flows with small Reynolds number*,
`arXiv:math/0702506`, **Ann. IHP Anal. Non Linéaire 26 (2009), 181–189**: "**If the
Reynolds number is small enough** we provide an elementary short proof of the existence of
global in time Hölder continuous solutions." Small data. The parallel probabilistic line
(Le Jan–Sznitman cascades; Bhattacharya–Chen–Dobson–Guenther–Orum–Ossiander–Thomann–
Waymire, *Majorizing kernels and stochastic cascades…*, **Trans. Amer. Math. Soc. 355
(2003), 5003–5040**) yields "unique global solutions … under **small initial data
conditions**". Every global result in the probabilistic literature is small-data.

### 2.2 The 2026 continuation says so explicitly, and quantifies the obstruction

`arXiv:2608.16915v1`, *Exact mean-covariance dynamics of the Weber field in the
stochastic Lagrangian representation of the 3D Navier–Stokes equations*, single author
(T. Mahithitarmmatorn), submitted 24 Jul 2026, 12 pp, math.PR/math.AP, no journal
reference. **Provenance caveat, per S56 §1.7: single unknown author, unrefereed. I use it
only for its explicitly-stated negatives, which are the useful part regardless of the
paper's fate.** Abstract, verbatim in the operative parts:

> "**No claim is made regarding global regularity; the missing steps are stated
> explicitly as open problems.**"

> "We further isolate several obstructions, each backed by an explicit counterexample —
> **the norm of an expected deformation gradient does not control the expected norm**,
> even for genuine stochastic flows of smooth divergence-free drifts; a spatial Hölder
> bound on the mean does not imply parabolic Campanato decay — and we prove a
> gauge-invariant quotient criterion **equivalent to the Serrin norm** of the velocity …"

Rendered specifics: Proposition 6.1 gives volume-preserving diffeomorphisms with
$\|\mathbb{E}\nabla A\|_{C^\alpha}=1$ but $\mathbb{E}\|\nabla A\|_{C^\alpha}\ge c_\alpha N$;
Proposition 6.3 gives a genuine stochastic flow with smooth divergence-free drift and
$\|\mathbb{E}\nabla A_{t_0}\|_{C^\alpha}\le C_\alpha$ while
$\mathbb{E}\|\nabla A_{t_0}\|_{L^\infty}\ge c\sqrt{\nu t_0}K^{1-\alpha}\to\infty$.
Theorem 7.1's criterion is a two-sided equivalence with $\|u\|_{L^p(I;L^q)}$, and the
quotient is gauge-invariant because it "**discards precisely the gradient part that the
Leray projection cannot see**."

Two consequences for this programme, stated conservatively:

1. **The exchange step is provably false in general.** Controlling
   $\|\mathbb{E}\nabla A\|$ (which is what the representation naturally bounds) does not
   control $\mathbb{E}\|\nabla A\|$ (which is what the nonlinearity needs). That is the
   probabilistic avatar of exactly the failure Wall 1 records: a bound in an averaged /
   integral sense does not upgrade to the sup-norm the budget requires.
2. **The one criterion extracted is equivalent to Serrin.** Fourth time this programme
   has met that pattern (S50 saturation ≡ Wall 1; S52 (AX); S57 $\mathrm{bmo}_\phi$; now
   the Serrin–Weber quotient). It is a relabelling, and the paper says so.

### 2.3 Rough paths and regularity structures: structurally aimed at a different difficulty

**M. Hairer**, *A theory of regularity structures* (the June 2015 revision of the
Inventiones paper), read from the primary PDF. Three verbatim quotes settle this.

> "Our main assumption will be that the equation described by (1.1) is **locally
> subcritical** (see Assumption 8.3 below). Roughly speaking, this means that if one
> rescales (1.1) in a way that keeps both $Lu$ and $\xi$ invariant then, at small scales,
> **all nonlinear terms formally disappear**."

> **Lemma 8.10.** "… our assumption of local subcriticality, Assumption 8.3, is really
> the **correct** assumption for the theory developed in this article to apply" (the
> lemma is an *if and only if*).

> Remark 1.18: "**Again, we do not claim that the solutions constructed here are
> global.** Indeed, the convergence holds in the space $C([0,T],\mathcal{C}^\alpha)$, but
> only up to some possibly finite explosion time."

and, from §1.1: "We furthermore consider parabolic problems, where we need to deal with
the problem of initial conditions and **local (rather than global) solutions**."

The three-line verdict:

1. **Different difficulty.** Regularity structures exists to make sense of a nonlinearity
   applied to a distribution too rough to multiply — a *local well-posedness* problem. 3D
   NS with smooth data has no such difficulty; its nonlinearity is classically defined.
   The Prize difficulty is a **global a priori bound**.
2. **Opposite hypothesis.** Local subcriticality ("at small scales all nonlinear terms
   formally disappear") is the precise negation of the supercriticality that makes 3D NS
   hard (Tao's objection, PRD §2). And Lemma 8.10 shows the hypothesis is *necessary* for
   the machinery, not merely convenient.
3. **Wrong output.** Local, not global.

Everything found actually applying these tools *to NS regularity* modifies the equation:

- **F. Flandoli, D. Luo**, *High mode transport noise improves vorticity blow-up control
  in 3D Navier–Stokes equations*, **Probab. Theory Relat. Fields 180 (2021), 309–363**:
  "a suitable multiplicative noise of transport type has a regularizing effect …
  provides a bound on vorticity which gives well posedness, **with high probability**.
  The result holds for **sufficiently large noise intensity and sufficiently high
  spectrum** of the noise."
- **F. Flandoli, M. Hofmanová, D. Luo, T. Nilssen**, *Global well-posedness of the 3D
  Navier–Stokes equations perturbed by a deterministic vector field*, `arXiv:2004.07528`,
  **Ann. Appl. Probab. 32 (2022), 2568–2586**: an **additional transport-type term** is
  added to the vorticity formulation, and global well-posedness is obtained "for large
  initial data **outside arbitrary small sets**", by "probabilistic methods, rough path
  theory, and a new Wong–Zakai approximation result."

Both are real theorems and neither is about the Prize equations. (A)/(B) is a
universally-quantified statement about the unperturbed system; "with high probability"
and "outside arbitrarily small sets" are the two ways this literature is structurally
unable to deliver it. The remaining rough-path/regularity-structures literature for 3D NS
is about *space-time white noise forcing* (Zhu–Zhu; Hofmanová–Zhu–Zhu, global existence
and **non-uniqueness**), which the handover correctly excluded as a different question.

### 2.4 Angle 2 — verdict

**Clean negative, with an unusually crisp reason.** The stochastic Lagrangian
representation keeps the Calderón–Zygmund operator (as $\mathbb{P}$), its circulation
theorem is a pathwise Euler identity with no dissipative content, its authors say known
criteria are *translated* rather than improved, its only global results are small-data,
and its 2026 continuation proves the crucial exchange step false and its extracted
criterion Serrin-equivalent. Rough paths and regularity structures are aimed at
*sub*critical ill-posedness and deliver *local* solutions, and every regularity-flavoured
application modifies the equation. **Neither offers a route around Wall 1's specific
obstruction; both are "different toolkit, same wall" in the S56 sense.**

---

## 3. Angle 3 — kinetic theory and hydrodynamic limits

Treated as the long shot the handover said it was. It is one, and the search returns a
clean negative with an extra structural reason.

### 3.1 The limit produces the *weakest* solution class, by design

**F. Golse, L. Saint-Raymond**, *The incompressible Navier–Stokes limit of the Boltzmann
equation for hard cutoff potentials*, **J. Math. Pures Appl. 91 (2009), 508–552**. Read
from the primary PDF. Abstract, verbatim:

> "The present paper proves that all limit points of sequences of renormalized solutions
> of the Boltzmann equation in the limit of small, asymptotically equivalent Mach and
> Knudsen numbers are governed by **Leray solutions** of the Navier–Stokes equations."

And §1, verbatim — the authors' own statement of why:

> "While the results above holds for global solutions of the Boltzmann equation without
> restriction on the size (or symmetries) of its initial data, earlier results had been
> obtained in the regime of smooth solutions [7,5]. **Since the regularity of Leray
> solutions of the Navier–Stokes equations in 3 space dimensions is not known at the time
> of this writing, such results are limited to either local (in time) solutions, or to
> solutions with initial data that are small in some appropriate norm.**"

So: the DiPerna–Lions → Leray programme is *constructed around* the NS regularity gap. It
inherits the gap; it does not touch it. And the smooth-solution branch of the same
literature is confined to local-in-time or small data — the same ceiling as §2.

**A small confirming detail.** Appendix C of the same paper is titled "Some regularity
results for the Leray projection", opening:

> "One annoying difficulty in handling incompressible or weakly compressible models is
> the **nonlocal** nature of the Leray projection $P$ … However, $P$ is **not continuous
> on local spaces of the type $L^p_{\rm loc}(dx)$**."

The same Calderón–Zygmund nonlocality that is Wall 1 shows up as an obstacle *inside* the
kinetic proof, in the authors' own words. The kinetic level does not dissolve it.

### 3.2 The modern state of the art runs the implication backwards

**K. Carrapatoso, I. Gallagher, I. Tristani**, *The Navier–Stokes limit of kinetic
equations for low regularity data*, `arXiv:2503.12046` (15 Mar 2025). Abstract, verbatim:

> "The main purpose of this work is to be as accurate as possible in terms of functional
> spaces. More precisely, **it is well-known that the Navier–Stokes equation can be
> solved in a lower regularity setting (in the space variable) than kinetic equations.**
> Our main result allows to get a rigorous link between solutions to the Navier–Stokes
> equation with such low regularity data and kinetic equations."

Their Theorem 1 **takes as input** a Navier–Stokes solution with data in
$H^{1/2}(\mathbb{T}^3)$ on $[0,T)$, and constructs kinetic solutions converging to it.
The information flows **NS → kinetic**. That is the opposite of what angle 3 would need.

### 3.3 Boltzmann's own regularity problem is open, and its conditional results assume the fluid bounds

**C. Imbert, L. Silvestre**, *Regularity for the Boltzmann equation conditional to
macroscopic bounds*, **EMS Surv. Math. Sci. 7 (2020), 1–56**. Abstract, verbatim from the
primary PDF:

> "From the mathematical point of view, **the existence of global smooth solutions for
> arbitrary initial data is an outstanding open problem.** … We prove that the solution
> will stay uniformly smooth **provided that** its mass, energy and entropy densities
> remain bounded, and away from vacuum."

So the flagship conditional-regularity programme at the kinetic level takes the
*macroscopic* (hydrodynamic) bounds as **hypotheses**. There is no reservoir of extra
control at the kinetic level to be spent on the fluid level; the dependency runs the other
way, and Boltzmann's own global regularity is at least as open as NS's.

### 3.4 What about the H-theorem specifically?

The H-theorem gives entropy dissipation, which in the incompressible NS scaling converges
to **exactly the Leray energy inequality** — this is the standard statement of the
Bardos–Golse–Levermore programme, and Golse–Saint-Raymond's paper is its completion. So
the kinetic-entropy structure delivers, in the limit, the one estimate NS already has,
and the one S56 §1.4(c) already established is insufficient (Palasek's shell model
satisfies the energy budget and blows up). Nothing in the moment hierarchy survives the
limit except the conserved moments and the entropy — i.e., mass, momentum, energy.

### 3.5 Angle 3 — verdict

**Clean negative, and it is the strongest of the three.** The kinetic route produces
Leray solutions; its own authors say the regularity of Leray solutions is not known and
that the smooth branch is confined to local/small data; the 2025 state of the art assumes
an NS solution as input; the kinetic equation's own global regularity is open with
conditional results that *hypothesise* the fluid bounds; and the H-theorem's limit is the
Leray energy inequality, an estimate already closed by citation at S56. **No regularity
criterion for NS has ever been derived through the kinetic description, and there is a
structural reason: the kinetic level is strictly downstream in the information flow.**

---

## 4. Recency check on what S56 already closed

Short, targeted, as instructed. **Nothing found dated after S56 (7–8 Sep 2026) that
changes the picture.**

- **BKM/CZ logarithm.** Targeted searches return only logarithmically *improved* Serrin
  criteria and the same 2020 Kanamaru optimality statement S56 recorded. No removal, no
  weakening. Several 2026 keyword hits are of the unprovenanced class S56 §1.7 flagged
  ("Global Regularity for Navier-Stokes on T3 via Bounded Vorticity-Response
  Functionals", the engrXiv "Equilibrium Depletion and Universal Frequency Envelopes",
  the Medium "Recursive Framework"); **not used, and no claim made about them**.
- **Ancient-Euler Liouville.** Nothing new found.
- **Grujić.** No version dated after S56's v3. (One index rendered the paper as v2/13 Jul
  2026, i.e. *staler* than S56's record; S57 read both papers in full from primary
  sources, so this is immaterial.)
- **Tao's blog, posts after 7 Sep 2026.** Four in this area: "Stable singularity of the
  Euler equations on $\mathbb{R}^3$" (10 Sep, guest post on Ganeshram–Duruisseaux–
  Anandkumar — **Euler**, **unforced**, **blowup**, PINN-discovered self-similar profile
  with scaling exponent converging to $0.5$, "certified … with high enough precision to
  carry out the stability analysis"); "Crowdsourcing a list of general resources on AI and
  mathematics" (10 Sep); "A Severe Misalignment of AI in Mathematics" (11 Sep); plus
  Hodge/sofic/Andrews–Curtis items (11 Sep) outside this area. The fluid item is the one
  S56 already recorded via the 7 Sep post's edit. **The AI asymmetry S56 identified holds:
  every effort targets blowup.**

**Methodological caveat, stated honestly.** The arXiv API was persistently rate-limited
through this session's proxy ("Rate exceeded" on every attempt across ~40 minutes and two
backoff strategies), so I could **not** run the metadata sweeps S56 used. The recency
check above rests on targeted web search plus direct primary-source fetches. It is
therefore weaker evidence than S56's sweeps, and a paper indexed under vocabulary my
searches missed could exist. I flag this rather than presenting the check as equivalent.

---

## 5. Recommendation

### 5.1 The recommendation, in one sentence

**No new work package. All three angles reduce to Wall 1 or to "a different toolkit that
meets the same wall in its own vocabulary", and I decline to manufacture a target from
any of them.**

This is the **sixth** independent confirmation of Wall 1 after S48, S50, S52, S53 and
S57 — and the first that reaches the programme's *own founding thesis* rather than an
external candidate. That is worth stating plainly: the audit was aimed inward as much as
outward, and the inward result is the same as the outward one.

### 5.2 What each angle would have had to show, and did not

| Angle | What a positive would have looked like | What was found |
|---|---|---|
| 1. Scalar-$T$ thesis | A regularity argument for R1/R2 that never takes $\|\nabla u\|_\infty$ | Wall 1 transfers verbatim by a *kinematic* argument; the max principle's input **is** $\|\nabla u\|_\infty$; the scalar's a priori content **is** the energy; both limits **are** the Prize problem |
| 2. Stochastic Lagrangian | An estimate on $\nabla u$ or $\omega$ bypassing the $L^\infty$ CZ recovery | $\mathbb{P}$ is still in the formula; the circulation theorem is a pathwise Euler identity; authors say criteria are "translated"; global results are small-data; the exchange step is provably false |
| 2′. Rough paths / reg. structures | A global a priori bound for deterministic 3D NS | Founding hypothesis is *local subcriticality* (necessary, by Lemma 8.10) — the negation of NS's difficulty; output is local; all NS-regularity applications modify the equation |
| 3. Kinetic / hydrodynamic limits | Kinetic structure supplying control on $\nabla u$ the PDE lacks | Limit is **Leray**; authors state NS regularity is the *obstacle*; 2025 SOTA **assumes** an NS solution; Boltzmann's own regularity is conditional on the fluid bounds |

### 5.3 One thing worth recording, if anyone writes to `proofs/` again

Nothing here is a new mathematical object worth a section. But two items are cheap,
internal, and would prevent a future session from re-walking angle 1:

- **A one-paragraph Wall entry** ("Wall 11 — the thermodynamic reformulation inherits
  Wall 1 by construction"), whose content is §1.2 Steps 1–3 only: Wall 1 is a property of
  $\operatorname{div}u=0$; `prop:s51-poloidal-sharp` is a static divergence-free field;
  therefore it is a counterexample for any reformulation that keeps
  $\operatorname{div}u=0$ exactly, including Route 1, Route 2, and the Constantin–Iyer
  representation. That is three sentences and one existing reference, with no new
  estimate.
- **A `rem:` correcting PRD §3.5(c) and §3.6**, since both are wrong in the repository's
  authoritative specification and a future session reading the PRD forward would be
  misled. §3.5(c) has the sign of $d\!\int\!T/dt$ backwards; §3.6 proves too much.

**I am not recommending these as a session.** They are housekeeping, they are optional,
and by the S53/S55/S57 standard "documentation-sharpening only" is not progress toward
(A)/(B). If the owner wants them, they are an hour, not a work package.

### 5.4 What this implies for the programme's next move

S56 §3.3 put a different research question on the table and explicitly declined to
recommend it. S57 returned F2 (the expected negative). This audit returns F2 three more
times, and additionally closes the programme's own origin. The honest summary:

> **The cheapest remaining moves against (A)/(B) have all been tried. The internal
> toolkit is exhausted (S52, S53, S56 §2). The best external candidate was priced and
> refuted (S57). The two largest untried external methodologies — stochastic/probabilistic
> and kinetic — are now assessed and are structurally unable to deliver a universally
> quantified a priori estimate. And the programme's founding thermodynamic thesis, never
> before pressure-tested against Wall 1, inherits it by construction.**

That leaves exactly three honest options, in the S46/S56 tradition of stating them so the
owner decides rather than drifting:

1. **Stop adding to (A)/(B) and consolidate.** The document, the walls ledger and the
   book are a genuine and unusual artefact: a fully documented record of *why* a
   well-defined line of attack does not close, with six independent confirmations. That
   has value as a finished object. This is the option the evidence most supports.
2. **Accept S56 §3.3's different question** (how much forcing a blowup construction
   actually requires, and whether the requirement is removable), with its stated caveat
   that it is work toward (C)/(D), not (A)/(B). Unchanged from S56; still not recommended
   by me, still available.
3. **Wait for an external breakthrough on Wall 1 itself.** The wall's own Consequence
   clause names the only permitted class ("any reformulation in which $\nabla u$ is not
   recovered by an $L^\infty$ Calderón–Zygmund bound"). S57 found and verified exactly one
   technique in that class that works — the commutator/Coifman–Rochberg–Weiss step — and
   showed the smallness it consumes is not available. If a second such technique appears,
   the programme is well positioned to price it in one session; that has now been done
   twice and the machinery is proven. Standing, zero-cost.

**Named failure modes, stated in advance, for whichever is chosen.** (Included because the
handover asked for them and because a "no" recommendation still has failure modes.)

- **G1 — the negative is read as an impossibility proof.** Nothing here proves Wall 1
  cannot be removed. §1 proves that *these three routes* do not remove it. A future
  session must not cite this document as evidence that the logarithm is irremovable;
  `subsec:wall1-cz`(c) already says the same about the wall itself.
- **G2 — angle 1 is re-opened numerically.** Route 1 and Route 2 solvers exist in
  `layer2/`, `layer3/`. Running them proves nothing about §1: the obstruction is an a
  priori estimate, and Wall 6 already classifies this programme's numerics as unable to
  serve as evidence for an estimate. **A Route 1 or Route 2 numerical campaign is
  forbidden in advance as a response to this audit.**
- **G3 — §1.8's scaling is upgraded silently.** The Conjecture 3.5 balance is a formal
  self-consistency calculation, explicitly not a theorem, and it lands in Wall 3 territory
  even if made rigorous. It must not be written into `proofs/` as anything else.
- **G4 — the 2026 stochastic preprint is cited as a source.** `arXiv:2608.16915` is an
  unrefereed 12-page preprint by an unknown single author. I used it **only** for its
  self-declared negatives. Nothing should be built on its positive content without
  independent verification (F5 of S56 applies verbatim).
- **G5 — "exhausted" is heard as "the problem is closed".** It is not. It means this
  programme's reachable moves are used up, which is a statement about the programme, not
  about the mathematics.

---

## 6. Verification record and non-claims

**Read as primary text this session:** `CLAUDE.md` in full; `HANDOVER-S58-AUDIT2.md` in
full; `HANDOVER-S56-AUDIT-FINDINGS.md` in full; `PROGRESS.md` S56 and S57 KEY FINDING
blocks in full plus the S56/S57 rows of the session-log table; `lem:local-enstrophy-typeI`,
`lem:BS-split`, `prop:first-rung`, `rem:routine-list`, `target:K-absorption` (lines
5680–5789); `subsec:wall1-cz` in full (lines 15141–15213); `prop:s51-poloidal-sharp` and
Step 1 of its proof (lines 11676–11705); PRD v0.4 §§3.1–3.7 text-extracted from
`tfirst_prd_v4.docx`; Constantin–Iyer `math/0511067v4` **in full from the primary PDF**
(Theorems 1.1, 2.2, 2.6, Propositions 2.1, 2.9 and the §2 remark on BKM); Golse–Saint-
Raymond, *J. Math. Pures Appl.* 91 (2009), **from the primary PDF** (abstract, §1
introduction, Appendix C); Imbert–Silvestre, *EMS Surv. Math. Sci.* 7 (2020), **from the
primary PDF** (title page and abstract); Hairer, *A theory of regularity structures*
(2015 revision), **from the primary PDF** (§1.1 local subcriticality, Remark 1.1,
Remark 1.18, Lemma 8.10). Abstracts fetched and quoted for `math/0702506`,
`arXiv:2608.16915`, `arXiv:2503.12046`, `arXiv:2004.07528`, `arXiv:2607.18939`.
Tao's September 2026 blog archive index fetched, and the 10 Sep post fetched in full.

**Derived here, not taken from any source:**
the kinematic transfer argument of §1.2 (Steps 1–3), including the observation that
`prop:s51-poloidal-sharp` is a counterexample for Route 1/Route 2 *verbatim*;
the Route 1 vorticity equation (R1-vort) with its three extra terms (E1)–(E3) and their
$\nabla T$/$\nabla^2T$ requirements; the parabolic $L^\infty$ threshold computation and
the chain $\nabla u\in L^{5}\Rightarrow u\in L^5_tL^{15/2}_x\Rightarrow$ LPS;
the exact conservation $\tfrac12\|u\|_{L^2}^2+\int T$ and the correction to PRD §3.5(c);
the identity $\int\theta(t)=\tfrac12(\|u_0\|^2-\|u(t)\|^2)$;
the observation that Boccardo–Gallouët's $q<5/3$ is exactly $u\in L^{10/3}(Q)$ translated;
the scaling table of §1.4; the §1.6 "proves too much" refutation of PRD §3.6;
the §1.7 $\delta$-family argument that both limits *are* the Prize problem;
the §1.8 self-consistent balance $a=4/(4+3\delta)$;
the observation that the Leray projection in Constantin–Iyer's formula is the same CZ
operator as Wall 1's, and that Proposition 2.9's pathwise character makes it an
Euler-level identity.

**Explicit non-claims.**
1. I have **not** proved that Wall 1 cannot be removed. §1 shows Route 1 and Route 2 do
   not remove it; that is a statement about those two systems.
2. §1.8's scaling balance is **formal**, not a theorem. Its only load-bearing use is item
   3 of that subsection, which does not depend on the balance.
3. The parabolic thresholds ($q>(n+2)/2$ for $L^\infty$; Boccardo–Gallouët $q<(n+2)/n$)
   are quoted as standard and justified here by scaling; I did **not** fetch the LSU or
   Boccardo–Gallouët primary texts this session.
4. I read **abstracts only** for Bulíček–Feireisl–Málek 2009 and Bulíček–Málek–Shilkin
   2014 (both behind ScienceDirect paywalls; the citation data came from the authors'
   own institutional publication pages at `karlin.mff.cuni.cz`). I did not read either
   proof.
5. I read **abstracts and rendered theorem statements only** for `arXiv:2608.16915`; I
   did not verify any of its propositions. It is an unrefereed single-author preprint and
   this audit takes no position on whether its results are correct.
6. I did **not** read Flandoli–Luo (2021) or Flandoli–Hofmanová–Luo–Nilssen (2022) beyond
   their abstracts and journal data.
7. **The arXiv API was rate-limited throughout this session**, so the systematic metadata
   sweeps S56 ran were not possible. §4's recency check rests on targeted web search and
   primary fetches and is correspondingly weaker. Absence of evidence, not evidence of
   absence.
8. I did not read, open, or modify anything in `book/`, `layer3/`, `layer4/`, `results/`
   or `proofs/`, other than reading the four `proofs/` passages listed above.
