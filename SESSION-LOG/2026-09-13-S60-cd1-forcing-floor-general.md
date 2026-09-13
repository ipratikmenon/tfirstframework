# S60 — CD-1: reconnaissance for the (C)/(D) line — generalising the forcing floor, and the vanishing-viscosity question for the OpenAI Euler construction

**Date:** 2026-09-13
**Branch / worktree:** `claude/cd1-forcing-floor-general` at `/home/user/cd1-worktree`,
based on `732a3cd` (`[S59] plan: pivot active line to (C)/(D) -- how much forcing, is it
removable`).
**Brief:** `HANDOVER-S60-CD1.md`, read in full and executed as written.
**Scope:** `proofs/claim_a_3d_proof_attempt.tex`, `PROGRESS.md`, `SESSION-LOG/` only.
No numerics written or run. `book/`, `layer3/`, `layer4/`, `results/` untouched.
Not committed, not pushed.

**Session tag: verified twice, at the start and again immediately before writing
`PROGRESS.md` and this log.** `git log --oneline -1` = `732a3cd` both times;
`ls SESSION-LOG/` highest = **S59** both times; `S58` appears in commit `26ba1e1`'s
message but has **no** session-log file (flagged by S59, still true). **S60 is
correct.** No other agent branch or worktree present. **Label prefix is `cd1`, not
`s60`**, following S53/S54/S55/S57: a work-package prefix cannot collide.

**Live file verified before touching anything:** `proofs/claim_a_3d_proof_attempt.tex`
was 16,193 lines at the tip, matching the live working copy exactly.

**Build:** baseline **rebuilt from the tip this session** rather than taken from any
previous log: 0 errors, 0 LaTeX warnings, 109 hyperref bookmark warnings, 65 overfull /
2 underfull hboxes, 215 pp, `check_proof.py` exit 0 (756 labels). *(Note: S55 recorded
64 overfull; the document has grown since. The baseline used here is the measured one,
not the recorded one.)* After the edit: **0 errors, 0 LaTeX warnings, 109 bookmark
warnings, 65 overfull / 2 underfull — every count identical to baseline, and the
overfull *set* is identical line for line.** 221 pp. `check_proof.py` exit 0 (786
labels, 551 `\ref` targets, 61 `\cite` keys, all resolved). Diff **purely additive**:
+447 lines, **0 deletions**.

---

## 0. One line per deliverable

| Deliverable | Answer |
|---|---|
| **Thread 1** — generalise the forcing floor from one profile to a class | **A genuine general theorem exists and is now written.** `thm:cd1-general-floor` + `cor:cd1-rate`, valid for **ν ≥ 0** (Euler included). Its sharpest form needs **no geometric hypothesis at all**. `prop:wpt3-floor` is the special case `(a,b) = (h,1/2)`, **with its constant recovered exactly**. |
| **Thread 1, the correction inside it** | **S55's own record was too generous to the OpenAI profile.** Exact self-similarity is load-bearing for **neither the exponent nor the constant** — only for time-localisation. |
| **Thread 1, bonus** | `prop:cd1-Lq`: a **positive-measure** (`L^q`) floor, answering in one specific sense the uniform-version question S55 left open. Does **not** apply to the OpenAI leading profile, and this is stated. |
| **Thread 2** — does the OpenAI *Euler* mechanism survive viscosity? | **Structural negative at any fixed ν > 0, derived from that paper's own displayed conditions**, corroborated three independent ways. **But the ν → 0 question itself is Constantin's open "gap between the two ranges of viscosities" and is NOT tractable here.** Nothing written to `proofs/`. |
| **Recommendation** | **One scoped work package, CD-2, on Thread 1 only.** Thread 2 is closed as a survey; pursuing it further is **not** recommended. |

---

## 1. Thread 1 — the general theorem

### 1.1 Re-derivation first, per standing discipline

Everything carried over was re-read directly in the proof document at its current line
numbers and re-derived before use.

| Carried-over claim | Location verified | Verdict |
|---|---|---|
| `lem:wpt3-forced-swirl`, statement and proof | lines 13628–13680 | **Correct.** Re-derived; the three ingredients (pressure, advection, vector Laplacian) each check out under axisymmetry. |
| `lem:wpt3-forced-max` (H1)–(H3) and (i)–(v) | lines 13687–13864 | **Correct.** The Danskin argument, both bounds, is sound as written; (H2) is genuinely load-bearing for the upper bound. |
| `prop:wpt3-floor` Steps 1–6 | lines 13914–14015 | **Correct.** `Γ⁽⁰⁾ = q^{-h}H` re-derived; `M(t) = c_*τ^{-h}` exact; `r(x_m) = c_rτ^{1/2}`. |
| `rem:wpt3-pointwise` (a)–(c) | lines 13866–13890 | **Correct**, and still applies verbatim to the new theorem's parts (i)–(ii). |
| S55 log: "this exactness … is what makes the floor exact rather than heuristic" | S55 log §3 Step 2 | **Overstated.** See §1.4. |
| Handover: "identify which hypotheses were profile-specific" | handover | Done; four-way audit at `rem:cd1-audit`. |
| Handover's guess that this is S60 | handover | **Correct**, verified twice. |
| Baseline 109 bookmark / 64 overfull (S55) | S55 log | **Bookmark count correct; overfull is now 65.** Baseline rebuilt from the tip rather than trusted. |

**One thing S55 stated for `ν > 0` that is true for `ν ≥ 0`.** Neither lemma's proof
uses positivity: in `lem:wpt3-forced-swirl` the viscous term is deleted from the
definition of `f` and from the swirl equation simultaneously; in
`lem:wpt3-forced-max` the only appearance of `ν` is the inequality
`ν(∂_r²Γ + ∂_z²Γ) ≤ 0` at an interior maximum, which at `ν = 0` is `0 ≤ 0`. Recorded as
`rem:cd1-euler`. **This matters for the (C)/(D) line specifically**, because the
forced-blowup constructions the field is now removing forcing from — Córdoba–
Martínez-Zoroa, Alpöge–Buckmaster — include forced 3D *Euler*.

### 1.2 The theorem

`thm:cd1-general-floor`. `ν ≥ 0`, `(u,p)` axisymmetric and smooth, `div u = 0`,
`f := ∂_t u + (u·∇)u − νΔu + ∇p`, hypotheses (H1)–(H3) of `lem:wpt3-forced-max`,
`Γ = r u_θ`, `M(t) = max Γ > 0`, `Σ(t)` the maximiser set,
`R(t) = max_{Σ(t)} r`, `F(t) = max_{Σ(t)}|f_θ|`.

```
(i)   RATIO FORM — no geometric hypothesis at all.
      At some x_* in Sigma(t):     f_theta(x_*)/u_theta(x_*)  >=  (d+/dt) log M(t)

(ii)  RADIUS FORM                  F(t)  >=  (1/R(t)) * (d+/dt) M(t)

(iii) INTEGRATED FORM — no differentiability of M anywhere
      int_s^t R(sigma) F(sigma) dsigma  >=  M(t) - M(s)
```

Form (i) is the one worth keeping. It has no core radius, no growth rate, no
self-similarity in it, and it reads:

> **at a swirl maximum, the azimuthal residual measured against the local swirl is at
> least the logarithmic growth rate of the swirl maximum.**

Under `M ≍ τ^{-a}` the right side is `a/τ`, so `|f_θ| ≥ (a/τ)|u_θ|` at that point: the
residual is at least as large as `∂_t u_θ` would be at that growth rate. **That is the
general form of `rem:wpt3-calibration`**, which S55 could only establish for one profile
by comparing four exponents of that profile's own θ-momentum equation. Here it follows
from the swirl maximum principle alone.

### 1.3 The rate form, and that it contains S55 exactly

`cor:cd1-rate`. Add
**(G1)** `M(t) ≥ c₁τ^{-a}`, `a > 0` — a prescribed (Type-II) swirl growth rate; and
**(G2)** `R(σ) ≤ R₀ τ^b` — a prescribed core-radius geometry. Then
`|f_θ| ≤ Kτ^{-p}` at the maximisers forces

```
p >= 1 + a + b ,     and at equality     K >= a*c1/R0 .
```

Proof is three cases on `b − p` against `−1` applied to the integrated form; the
critical case gives `c₁ − M(s₀)τ^a ≤ R₀K τ^{a−β}/β` with `β = p − b − 1`, and letting
`τ ↓ 0` yields both statements.

`cor:cd1-recovers`: for the OpenAI leading field, `eq:wpt3-Mexact` is (G1) with
`a = h`, `c₁ = c_*`, and `eq:wpt3-rm` is (G2) with `b = 1/2`, `R₀ = c_r`. Hence
`p ≥ 1 + h + 1/2 = A + 1` and `K ≥ h c_*/c_r` — **exactly `eq:wpt3-floorB`, exponent and
constant.** The general theorem contains the S55 result rather than paralleling it.

### 1.4 The precise correction to S55's record — the finding the handover asked for

S55's log recorded, of `M(t) = c_*τ^{-h}` being an exact power rather than `≍`:
*"This exactness — not `≍` — is what makes the floor exact rather than a scaling
heuristic, and it is available only because the leading field is exactly self-similar."*

**That is too generous.** `rem:cd1-selfsimilar` splits it three ways:

- **(a) The exponent does not need it.** `p ≥ 1 + a + b` follows from the one-sided
  bounds (G1)–(G2) alone, and because the proof routes through the *integrated* form,
  it needs no differentiability of `M` at any point.
- **(b) The constant does not need it either.** `K ≥ a c₁ / R₀` is recovered from the
  asymptotic constants alone.
- **(c) What exactness actually buys is time-localisation, and only that.** With it, the
  floor holds *at every `t`*, at an identified point. Without it,
  `cor:cd1-rate` excludes a *uniform* power bound on the interval, hence yields only
  "for every `s < T` there is a `σ ∈ (s,T)` at which the floor holds", with `σ` not
  identified.

**The one genuinely profile-specific step of S55 was Step 5** — the derivation of the
maximiser radius `c_rτ^{1/2}` as an *output* of the exact profile. In the general
theorem that radius is an *input*, (G2). And form (i) removes even that, by measuring
the residual against the local swirl rather than in absolute units.

So: **exact self-similarity was a convenience of presentation, not a hypothesis of the
mathematics.** This is a real correction to how S55 characterised its own result, and it
is the answer the handover asked for on "which feature of the OpenAI profile was
load-bearing".

### 1.5 The hypothesis audit (`rem:cd1-audit`)

- **Generic, doing no work:** (H1) smoothness, (H3) orientation convention.
- **Generic, doing all the work:** the pressure-free swirl equation, and the
  annihilation of *both* the drift and the singular `−2ν/r` term at an interior
  maximum. Neither uses anything beyond axisymmetry and `div u = 0`.
- **The one hypothesis that must be checked per construction:** **(H2)**, uniform
  compactness of the maximiser set. It is what makes Danskin's upper bound available,
  and it fails when the swirl maximum escapes to spatial infinity. For the S55 profile
  it needed Step 3; for a compactly supported construction it is automatic.
- **The class boundary is axisymmetry and it is absolute.** The only lever is
  `(∇p)_θ = r^{-1}∂_θ p = 0`. Outside the axisymmetric class there is no scalar swirl
  equation and **no analogue is proved or claimed**. This is why the theorem, like
  §`sec:wpt3-forcing-floor` before it, says nothing whatever about Theorem 1.1 of
  `\cite{OpenAI2026}`, whose field is non-axisymmetric. **Wall 9 governs and is
  unchanged.**

### 1.6 The positive-measure form (`prop:cd1-Lq`)

S55 recorded, and deliberately did not pursue, the question of whether a uniform
(positive-measure) version of the floor is provable. One is, in a specific sense, and it
costs four lines:

```
d/dt ||Gamma(.,t)||_{L^q}  <=  || (r f_theta)(.,t) ||_{L^q} ,      2 <= q < infinity, nu >= 0
```

for compactly supported `Γ`. Testing the swirl equation with `|Γ|^{q−2}Γ`: the drift
integrates away by `div u = 0`; the diffusion gives `−q(q−1)ν∫|Γ|^{q−2}|∇Γ|² ≤ 0`; and
**the singular `−2ν/r` term vanishes identically** because `dx = r\,dr\,dθ\,dz` cancels
the `r^{-1}` and `Γ` vanishes on the axis — again no ε-exhaustion, the same structural
reason as `lem:wpt3-forced-max`(i). Hölder finishes it.

**Three honest caveats, all recorded in `rem:cd1-Lq-scope`:**
1. It bounds `‖r f_θ‖_{L^q}`, not `‖f_θ‖_{L^q}`; converting needs a radius bound on the
   *support of the residual*, strictly stronger than (G2).
2. It does **not** contradict `rem:wpt3-pointwise`(b) — that remark says no `L^q` bound
   follows *from the pointwise statement*, and the spike counterexample still applies.
   This is a different argument from a different hypothesis (growth of `‖Γ‖_{L^q}`, not
   of `max Γ`), and the two are independent.
3. **It does not apply to the OpenAI leading field.** There
   `Γ⁽⁰⁾ = q^{-h}H` with `H ~ √2 c_∞ X^{-h}` and `X = r²/(2q)`, so
   `Γ⁽⁰⁾ → √2 c_∞ 2^h r^{-2h}` at large `X`, *independently of `q`* — not in `L^{q'}(ℝ³)`
   for any finite `q'`. This is a tool for compactly supported constructions.

---

## 2. Thread 2 — the OpenAI Euler paper and vanishing viscosity

**All sources below were fetched and read this session. Nothing is recalled.**
Nothing from this thread is written to `proofs/`.

### 2.1 The mechanism, from the primary source

Fetched: `https://cdn.openai.com/pdf/315b36cd-ec98-4023-8342-93345194ece1/euler.pdf`
(*Finite Time Blowup for the Euler Equation*, OpenAI, 56 pp), read in full for §1, §2,
§3.2, §5.2, §5.3, §5.5, and the reference list. Theorem 1.1 is as S56 recorded, verbatim
here:

> "There exists `u₀ ∈ C^∞_{c,σ}(ℝ³)` such that `0 < T_*(u₀) < ∞`. Its smooth Euler
> solution satisfies `limsup_{t↑T_*}‖∇u(t)‖_{L^∞} = ∞`,
> `∫₀^{T_*}‖curl u(t)‖_{L^∞} dt = ∞`."

**It is not a self-similar construction.** It is an iterated-oscillation scheme in the
Córdoba–Martínez-Zoroa lineage (its own references [10], [11], [13]): at stage `j` a
localized "packet" of localization length `ℓ_j` and frequency `k_j` is added over a
parent flow `U_{j−1}` whose gradient contains a rank-one shear of size `h_{j−1}`; the
packet is amplified by that shear and its terminal gradient becomes the next parent's
shear. From (2.5) the leading velocity increment is `(ℓα/k)χ₁ v f_δ(k m₀·y)` and the
leading gradient is `αχ₁ v⊗m f_δ'`, so — the paper's own words — *"Increasing `k` reduces
the velocity by a factor `k^{-1}`, while leaving this leading gradient unchanged."*
The physical wavenumber of the packet is `N_j = k_j/ℓ_j`.

Scales, (5.3) and (5.11) verbatim:

```
x_{J-1} = x0 ,   x_j = j^2 x_{j-1}
log h_j = x_{j-1}/j^5 ,  log k_j = x_{j-1}/j^2 ,  log l_j^{-1} = x_{j-1}/j^{7/2}
(5.11):  log k_{j-1} = x_{j-1}/(j-1)^4 ,   log h_{j-1} = x_{j-1}/(j-1)^7
```

*(I re-derived (5.11) from (5.3) independently before using it: `x_{j−2} = x_{j−1}/(j−1)²`
gives both identities exactly. They agree.)*

### 2.2 The structural obstruction, from the paper's own displayed conditions

The packet proposition (Proposition 3.1) is available only under

> **(3.12)** `P^Q ≤ k^{ϑ/100}`, `ϑ = 10^{−6}`

with `Q = Q(c₀)` fixed, and (5.10) lists `h_j` and `h_{j−1}` among the entries of the
size parameter `P_j`; (5.12) is the verification that `log P_j / log k_j < ϑ/(100Q)`.
Therefore, **from the paper itself**,

```
k_j  >  P_j^{100Q/theta}  >=  h_j^{1e8 * Q} .
```

**The frequency must exceed a fixed enormous power of the gradient the packet
produces.** This is not a numerical convenience: (3.12) is what makes the
inverse-frequency expansion of the nonlinear corrections converge, i.e. it is a
requirement of the method.

Against this, viscous marginality at viscosity `ν > 0` — the packet's dissipation rate
must not exceed the stretching rate that amplifies it — requires

```
nu * (k_j / l_j)^2  <~  h_{j-1} ,     i.e.    k_j  <~  l_j * (h_{j-1}/nu)^{1/2}  <=  (h_{j-1}/nu)^{1/2} .
```

The two are **incompatible in the exponent by a factor `≈ 2×10⁸ Q`** for every fixed
`ν > 0`, once `h_j → ∞`. Substituting the actual scales makes the size of the failure
explicit:

```
nu (k_j/l_j)^2 / h_{j-1}  =  nu * exp( x_{j-1} [ 2/j^2 + 2/j^{7/2} - 1/(j-1)^7 ] )
```

whose exponent is `≈ 2x_{j−1}/j²` — positive and divergent at every normal stage
`j ≥ J`, with `J` chosen large. The viscous damping factor over the amplification
interval is therefore doubly-exponentially small while the amplification it must beat is
only single-exponentially large.

**The same statement in Tao's amplitude/frequency language** (IPAM 2017 lecture notes,
fetched and read: *"if `A ≪ N` we expect dissipation to dominate … if `A ≫ N` we expect
nonlinear behaviour"*, and the energy constraint `A ≲ N^{3/2}`): the packet at stage `j`
has `A_j = h_j/N_j`, so `A_j/N_j = h_j/N_j² → 0` — deep in the dissipation-dominated
regime — and `A_j ≪ N_j^{3/2}` by an enormous margin, i.e. nowhere near the
energy-saturating cascade Tao identifies as the only viable NS blowup shape.

**And that is not an accident of the scale choice.** The property that makes
`A_j/N_j → 0` is the same property that gives
`Σ_j ‖U_j(0) − U_{j−1}(0)‖_{H^m} < ∞ for every m` (2.1) and hence
`u₀ ∈ C^∞_{c,σ}`: the velocity increment carries a factor `k^{-1}` while the gradient
does not. **The smoothness of the limiting datum and the viscosity-fragility of the
mechanism are the same fact.** That is the sharpest form of the Thread-2 finding, and it
is a statement about the *method*, not about these particular exponents.

### 2.3 Three independent corroborations

**(a) The lineage has been pushed against dissipation, and stopped at ~5% of it.**
Córdoba, Martínez-Zoroa and Zheng, `arXiv:2407.06776` (v2, 5 Aug 2024, 33 pp; to appear
*Arch. Ration. Mech. Anal.*), abstract verbatim:

> "we establish the formation of singularities of classical solutions with finite energy
> of the **forced fractional** Navier Stokes equations where the dissipative term is
> given by `|∇|^α` for any `α ∈ [0, α₀)` (`α₀ = (22 − 8√7)/9 > 0`)…"

`α₀ = (22 − 8√7)/9 = 0.09266…`, computed here, **out of `α = 2` for real viscosity** —
and still *with* forcing. Their §1.2.4 states the obstruction in exactly the form
derived above, verbatim:

> "First, the stretching caused by the outer vortex layers to `ω_n` is of magnitude
> `A_{n−1}`. For the blow-up to happen, **the stretching should beat the dissipation,
> which is of magnitude `M_n^α`**, so we have our first constraint: `M_n^α ≤ A_{n−1}`."

and they take `log M_n = R^n` (doubly-exponential frequencies), reaching
`s ≤ s(α) := (1+α)/2 − √(3α)`, positive only for small `α`. *(Back-of-envelope threshold
`s(α) > 0` is `α < 5 − 2√6 ≈ 0.1010`; the theorem's `α₀ ≈ 0.0927` is slightly below it,
consistent with their own note that "in reality there are other sources of error".)*
**This is the field's own measurement of how much dissipation this mechanism survives,
and the answer is: almost none.** Note also that OpenAI's Euler frequencies grow *faster*
than CMZ's — `log log k_j ~ 2j log j`, versus CMZ's `log log M_n = n log R` — so there is
even less room.

**(b) OpenAI's own Navier–Stokes paper does not use this geometry.** Re-fetched
`https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf` and
read §2.1 and §7.2. Verbatim:

> "`ℓ_r ≍ τ^{1/2}`, `ℓ_z ≍ τ^{1/2−h}`" … "`Re_θ := |u_θ|ℓ_r/ν ≍ τ^{−h} → ∞`,
> `Re_r := |u_r|ℓ_r/ν = O(1)`" … "**The radial Reynolds number remains bounded, so
> viscosity continues to compete with radial inflow.**"

and for the oscillatory pulses:

> "**Viscosity remains a leading-order term in the pulse equation.** … Thus `εk² ≍ 1` in
> the pulse equation throughout the correction process."
> "Background shear amplifies the waves; **viscous dissipation makes them decay**."
> "**Viscous damping strengthens and eventually exceeds the amplification**, so the pulse
> grows and then decays."

So when the same authors needed viscosity to be survivable, they abandoned the
iterated-packet geometry entirely and used a self-similar background whose radial scale
sits **exactly at the viscous scale** (`ℓ_r² ≍ ντ`, `Re_r = O(1)`) with amplitude and
wavenumber tuned to `εk² ≍ 1` — the marginal case — **and still needed forcing**. That
is the strongest available internal evidence for the Thread-2 conclusion.

**(c) The Euler paper's own use of viscosity is not physical.** Its only `ν` is an
*auxiliary* regularisation used to solve for the correction on a lifted phase space and
then removed (`§3.6`, "Existence and removal of viscosity … solving (3.51) for `ν > 0`
and then passing to `ν ↓ 0`"). It carries no information about the question here, and a
future session must not read it as if it did.

### 2.4 What is genuinely known, and what is genuinely open

**The only proved implication runs the other way.** Peter Constantin, *On the Euler
equations of incompressible fluids* (survey; fetched and read), verbatim:

> "It is known that if there are no singularities in the solution of the Euler equations
> with initial data `u₀` on the time interval `[0,T]`, then there can be no
> singularities in the Navier-Stokes solution with the same initial data and small
> enough viscosity ([36]). The regularity for large enough viscosities is also known.
> **Unfortunately, there is a gap between the two ranges of viscosities, and it is not
> clear how to close it.**"

(`[36]` = P. Constantin, *Note on loss of regularity for solutions of the 3-D
incompressible Euler and related equations*; the same reference is cited later in the
survey for "the limit holds for as long as the Euler solution is smooth".)

So: for smooth data the inviscid limit holds up to, but not through, the Euler singular
time; **Euler blowup ⟹ NS blowup at small `ν` is not known and is exactly Constantin's
gap.** Nothing found this session narrows it.

**The Onsager / anomalous-dissipation literature has no purchase on this construction
either**, and the reason is structural rather than a gap in the search: the construction
transfers essentially *no* energy to high frequencies — that is precisely what
`Σ_j ‖U_j(0) − U_{j−1}(0)‖_{H^m} < ∞` for every `m` means, and the packet energy
`≈ A_j² ℓ_j³` is summably tiny. A dissipation anomaly requires the opposite. So this is
exactly the kind of Euler singularity that is *not* an anomalous-dissipation mechanism,
and the mature `ν → 0` theory (Onsager's threshold `C^{1/3}`, Constantin–E–Titi, Isett)
neither constrains it nor is constrained by it.

### 2.5 Thread 2 verdict

- **At fixed `ν > 0`: the mechanism does not survive**, and the obstruction is
  structural — the method's own frequency-dominance condition (3.12) is incompatible in
  the exponent with viscous marginality. Corroborated by the lineage's hypodissipative
  ceiling `α₀ ≈ 0.093` and by the same authors' choice of a completely different, viscous-scale
  geometry when they needed viscosity to be survivable.
- **In the limit `ν → 0`: genuinely open, and not a question this programme can move.**
  It is Constantin's gap. Answering it would require exactly the kind of uniform-in-`ν`
  a priori estimate the whole subject lacks.
- **What is *not* claimed.** No NS statement follows: nothing here says no NS singularity
  exists, that some other mechanism cannot work, that the Euler theorem is wrong, or that
  the `ν → 0` limit of NS with this `u₀` is regular. `\cite{OpenAI2026}`'s Theorem 1.1 and
  the Euler paper's Theorem 1.1 are both untouched, and neither is independently
  verifiable from outside in one session.

---

## 3. Recommendation

### 3.1 In one sentence

**Take one scoped work package on Thread 1 — CD-2, below — and close Thread 2 as a
completed survey; do not open it as a session.**

### 3.2 CD-2, scoped

**Premise.** `thm:cd1-general-floor` is a general lower bound but it is currently a
theorem with no second instance. Its value is entirely in whether it applies to
constructions other than the one it was abstracted from. That is a finite, checkable
question and it is the natural next step.

**Deliverables.**
1. **Find the class's second and third members.** Apply `cor:cd1-rate` to the
   *axisymmetric* forced constructions that actually exist — Córdoba–Martínez-Zoroa
   `arXiv:2309.08495` (forced 3D Euler, `C^{1,1/2−ε} ∩ L²` force) and
   Córdoba–Martínez-Zoroa–Zheng `arXiv:2407.06776` (forced hypodissipative NS) — reading
   their `(a, b)` off the primary sources. **First check whether these constructions are
   axisymmetric at all**; if they are not, the answer is "the class has one member" and
   the session should say so and stop.
2. **Extend the floor to fractional dissipation.** `lem:wpt3-forced-swirl` and
   `lem:wpt3-forced-max` rest on the sign of the diffusion at an interior maximum.
   `(−Δ)^{α/2}` has a maximum principle too (the Córdoba–Córdoba pointwise inequality),
   so the lemma plausibly extends to `|∇|^α`, `0 < α ≤ 2`, with `ν(Δ − 2r^{-1}∂_r)Γ`
   replaced by a nonlocal operator whose sign at a maximum must be established, not
   assumed. **This is the one step that could fail and it should be attempted first.**
3. **Only if 1 and 2 both succeed:** state the floor once for `|∇|^α` and read off what
   it says about how much forcing the hypodissipative constructions require as
   `α ↑ α₀`.

**Named failure modes, stated in advance.**
- **H1 — the class has one member.** Every extant forced construction outside
  `\cite{OpenAI2026}`'s axisymmetric background is non-axisymmetric, so the general
  theorem never applies to anything else. *This is the most likely single outcome.* It is
  a legitimate result — it would say the theorem is a generalisation in form only — and
  the session should end there rather than manufacture an application.
- **H2 — the fractional swirl equation has no clean form.** `(−Δ)^{α/2}` applied to
  `u_θ e_θ` need not reduce to a scalar operator on `Γ` at all; the `−2/r ∂_r` structure
  is special to the local Laplacian. If so, say precisely where it breaks and stop.
- **H3 — rediscovery.** The nonlocal maximum principle for `Γ` in the axisymmetric class
  may already exist in the fractional-NS literature. *Mitigation: search first, derive
  second, attribute.*
- **H4 — drifting into (C)/(D) advocacy.** The floor is a lower bound on a residual. It
  is not evidence for or against any blowup claim, and CD-2 must not be written as if it
  were. `rem:cd1-scope` is the template.
- **H5 — the forcing prohibition being blurred.** Everything remains a statement about
  *given* fields with `f` defined as the residual. **No forcing term enters any argument
  of this document toward (A)/(B), ever.**

**Out of scope for CD-2:** numerics of any kind; `book/`; any statement about whether
any external theorem is correct; any non-axisymmetric extension.

**Definition of done:** either a second worked instance of `cor:cd1-rate` with its
`(a,b)` read off a primary source, or a precise statement of which of H1/H2 fired.

### 3.3 Why Thread 2 is not recommended as a session

The tractable part is already done and is in §2 of this log: the fixed-`ν` question is
answered structurally, from primary sources, and needs no further work. The remaining
part is Constantin's gap between the two viscosity ranges, which is a named open problem
of the subject, is not made easier by this construction's geometry, and would require a
uniform-in-`ν` estimate — i.e. the same class of object Wall 1 obstructs. **Opening it
would be manufacturing a target.** If anything from §2 is ever wanted in the document,
it belongs in a Wall entry of two paragraphs, not a section, and it would need a new
bibliography entry for the Euler paper (currently absent — `\cite{OpenAI2026}` is the
NS paper only).

---

## 4. Verification record and non-claims

**Derived here, not taken from any source:** `thm:cd1-general-floor` (i)–(iii) and the
observation that (i) needs no geometric hypothesis; `cor:cd1-rate` and its three-case
proof; `cor:cd1-recovers`; the extension of both S55 lemmas to `ν ≥ 0`; the
hypothesis audit of `rem:cd1-audit`; the correction to S55's self-similarity claim;
`prop:cd1-Lq` including the vanishing of the `−2ν/r` term against `|Γ|^{q−2}Γ`; the
observation that `Γ⁽⁰⁾ → √2c_∞2^h r^{−2h}` and so lies in no `L^q`; the re-derivation of
(5.11) from (5.3); the combination of (3.12) with (5.10)/(5.12) into
`k_j > h_j^{10⁸Q}`; the comparison of that with viscous marginality; the computation
`ν(k_j/ℓ_j)²/h_{j−1} = ν exp(x_{j−1}[2/j² + 2/j^{7/2} − 1/(j−1)⁷])`; the identification
of `A_j/N_j → 0` in Tao's language; the observation that `H^m`-summability of the initial
increments and viscosity-fragility are the same property; the arithmetic
`α₀ = (22−8√7)/9 = 0.0927…` and `5 − 2√6 = 0.1010…`.

**Read as primary text this session:** `CLAUDE.md` in full; `HANDOVER-S60-CD1.md` in
full; session logs S54 and S55 in full; `HANDOVER-S56-AUDIT-FINDINGS.md` and
`HANDOVER-S59-AUDIT2-FINDINGS.md` in full; `sec:wpt3-forcing-floor` in the proof document
at lines 13601–14314 (`lem:wpt3-forced-swirl`, `lem:wpt3-forced-max` with proof,
`rem:wpt3-pointwise`, `prop:wpt3-floor` with proof, `rem:wpt3-calibration`,
`cor:wpt3-annulus`, `rem:wpt3-scope`); the OpenAI **Euler** PDF in full for §1, §2,
§3.2 ((3.11), (3.12), Prop. 3.1), §5.2 ((5.3)–(5.7)), §5.3, §5.5 ((5.10)–(5.12)), §3.6
and the reference list; the OpenAI **Navier–Stokes** PDF §1 (related work), §2.1 (the
scales and Reynolds numbers) and §7.2 (`εk² ≍ 1`, and Figure 5's caption); Córdoba–
Martínez-Zoroa–Zheng `arXiv:2407.06776` — abstract and metadata from the arXiv listing
page, and §§1.2.3–1.2.4 in full from the arXiv HTML rendering; Constantin, *On the Euler
equations of incompressible fluids*, §2.1 and the reference list, text-extracted from the
primary PDF; Tao, *Finite time blowup constructions for modifications of the
Navier-Stokes and Euler equations* (IPAM 2017), the amplitude/frequency slides, from the
primary PDF.

**Explicit non-claims.**
1. `thm:cd1-general-floor` is **elementary** — the swirl maximum principle, Danskin, and
   one division. No credit is claimed for machinery; the content is the generality, the
   `ν = 0` extension, and the hypothesis audit.
2. It says **nothing** about non-axisymmetric fields, hence nothing about Theorem 1.1 of
   either OpenAI paper. No error in either is asserted, and none was found.
3. Parts (i)–(ii) of the theorem are **pointwise at a moving point**;
   `rem:wpt3-pointwise` applies verbatim. `prop:cd1-Lq` is the only positive-measure
   statement here.
4. §2's conclusion is that the OpenAI Euler *construction* is incompatible with fixed
   `ν > 0`. It is **not** a claim about whether Navier–Stokes blows up, nor about the
   `ν → 0` limit, nor a claim that no variant of the method could work.
5. I did **not** read Alpöge–Buckmaster (their preprints are still not on arXiv) and
   made no use of them beyond Tao's post as quoted at S56.
6. I did **not** verify any proof in `arXiv:2407.06776`; it is used only for its
   abstract, its stated threshold, and its own §1.2.4 statement of the obstruction, all
   quoted.
7. The claim that (H2) "must be checked per construction" is a statement about the
   lemma's hypotheses, not a claim that it fails for any particular construction.
8. I did not read, open or modify anything in `book/`, `layer3/`, `layer4/`, `results/`,
   or `layer1/`, `layer2/`.

---

## 5. Files changed

```
proofs/claim_a_3d_proof_attempt.tex   +447 / -0   (16,193 -> 16,640 lines)
      new sec:cd1-general-floor with rem:cd1-euler, thm:cd1-general-floor,
      rem:cd1-read, cor:cd1-rate, cor:cd1-recovers, rem:cd1-selfsimilar,
      rem:cd1-audit, prop:cd1-Lq, rem:cd1-Lq-scope, rem:cd1-scope.
      No bibliography entries added, no existing text altered. Purely additive.
PROGRESS.md                           KEY FINDING S60 block, session-log row,
                                      Active milestone
SESSION-LOG/2026-09-13-S60-cd1-forcing-floor-general.md   (this file)
```

Not committed, not pushed. Worktree `/home/user/cd1-worktree`, branch
`claude/cd1-forcing-floor-general`, for independent verification before merge.
**Re-check the session tag and the tip at merge time.**
