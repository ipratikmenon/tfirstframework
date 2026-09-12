# S53 — WP-T1: Type-II reconnaissance

**Date:** 2026-09-12
**Branch / worktree:** `claude/wpt1-typeii-recon` at `/home/user/wpt1-worktree`, based on
`2259b4d` ([S52] plan: Type-II reconnaissance handover)
**Scope:** `proofs/claim_a_3d_proof_attempt.tex`, `PROGRESS.md`, `SESSION-LOG/` only.
No numerics written or run. `book/`, `layer3/`, `layer4/`, `results/` untouched.
**Deliverable:** `\S sec:typeII-wpt1` (§39 in the rendered PDF), 12 new bibliography
entries, three new rows (3′, 3″, 3‴) in the S47 walls table.
**Build:** `check_proof.py` exit 0 (657 labels, 456 `\ref` targets, 53 `\cite` keys, all
resolved). `pdflatex` ×3 clean: 0 errors, 0 LaTeX warnings, **109 hyperref bookmark
warnings — exactly the pre-existing S50/S52 baseline**, checked rather than assumed.
188 pp (was 178). Diff **purely additive**: 666 insertions, 0 deletions.

**Session tag:** verified before writing anything. `git log --oneline -1` = `2259b4d`,
`ls SESSION-LOG/` highest = S52, so S53 is correct. Re-checked at the end: unchanged, no
other agent branch present. **Label prefix is `wpt1`, not `s53`** — three consecutive
sessions have had their tag claimed concurrently, and `sec:axisym-s51` already carries
`s51` labels for what the log calls S52. A work-package prefix cannot collide, and if this
turns out to be S54 only the section heading needs changing.

---

## 0. What was asked, and what the answer is in one line each

| Deliverable | Answer |
|---|---|
| **D1** — what is known about Type-II for true 3D NS | One author's programme (Seregin, four papers), **entirely scenario-conditional**. No unconditional exclusion exists. Type-II is **open** and is **not** "expected" — the only unforced numerical candidate is Type-I-shaped, the only Type-II profile is forced. |
| **D2** — is the vacuity fundamental or specific? | **Neither, as posed.** It is fundamental to *NS-preserving* rescaling and specific to nothing else; a rate-adapted rescaling exists and is published; its price is exact and unavoidable — **ν_eff = ν·A/L, so the limit object loses the viscosity and solves Euler.** |
| **D3** — verdict | **(iii)**, referral to a named branch of mathematics, justified by **(ii)**, a derived negative. **Not (i).** No technical target is proposed. |

---

## 1. Deliverable 1 — what is actually known

### 1.1 A third meaning of "Type I" was found, and it is the one the results use

`rem:s51-handover` (S52) records two conventions: **temporal** (this document, KNSS:
`M(t) ≤ C_I(T−t)^{-1}`, `sup|u| ≤ C_I(T−t)^{-1/2}`) and **spatial** (axisymmetric
literature: `|u| ≤ C/r`). The Type-II literature uses a **third**: Caffarelli–Kohn–%
Nirenberg's scale-invariant energies

```
A(v,r) = sup_{−r²<t<0} (1/r)∫_{B(r)}|v|²dx
E(v,r) = (1/r)∬_{Q(r)}|∇v|²dz
C(v,r) = (1/r²)∬_{Q(r)}|v|³dz
```

with Type-I meaning these stay bounded. Two *variants by the same author*:
`g = max{limsup A, limsup E, limsup C} < ∞` (arXiv:2304.04045) and
`g₀ = min{liminf A, liminf E, liminf C} < ∞` (arXiv:2606.29468). Since `g₀ ≤ g`, Type-II
in the `g₀` sense is the *stronger* requirement — all three must diverge.

**Relation to this document, established rather than guessed.** Temporal Type-I **implies**
the CKN one: the pointwise bound `√(T−t)|u| ≤ M'` implies the generalised Type-I bound
(the `A`-quantity uniformly over centres). This is *not formal* — it is a theorem of
Seregin–Zajączkowski 2006, and Barker–Prange (arXiv:2211.16215) state it in footnote 8 in
exactly those terms: *"It turns out to be true, see [78] and the review article [77, pages
844–849]. This makes (4.1) a good notion."* The converse is available from no source read.

**Consequence (`cor:wpt1-inclusion`):** Seregin's Type-II class is **strictly smaller** than
this document's. So the space of putative singularities is three regions, not two:

- **(R1)** temporal Type-I — the pincer's target; closed in the axisymmetric class by
  Seregin–Šverák 2009 (`fact:s51-ss`).
- **(R2)** temporal Type-II but CKN Type-I — **addressed by nobody**. This document needs
  the pointwise rate; Seregin's papers assume `g₀ = ∞`. Whether (R2) is empty is unknown.
  Seregin himself, on the corresponding Type-I question: *"The question about Type I blowups
  is whether the boundedness of g allows blowups or not. It is still open."*
- **(R3)** CKN Type-II — Seregin's programme.

Any future claim to have "covered Type-II" must say what it does about (R2). I do not know
whether (R2) is populated and I make no claim that it is; the point is that the two
literatures do not tile the space between them.

### 1.2 The Type-II literature: shape, and exactly what it excludes

Four papers, all Seregin, all read this session:

| Paper | Status |
|---|---|
| arXiv:2304.04045v1 (8 Apr 2023) = CPAA 23(10):1389–1406, 2024 | Prop. 1.2: the zoom construction |
| arXiv:2402.13229v3 (8 Aug 2026) | axisymmetric case |
| arXiv:2507.08733 (11 Jul 2025, rev. 3 Jan 2026) | Thm 5.1: exclusion under an added LPS hypothesis |
| arXiv:2606.29468v1 (28 Jun 2026) | Thm 2.1: exclusion for an `(f,g)` exponent relation |

Uniform shape: **(a)** fix a *scenario* — a rate function quantifying how fast
`M^{s,l}_κ(v,r)` diverges, plus weighted energy bounds; nothing is assumed about a
pointwise rate. **(b)** zoom by the **Euler scaling**
`v^{λ,α}(y,τ) = λ^α v(λy, λ^{α+1}τ)`, `α > 1`, extracting a nontrivial ancient solution of
the **Euler** equations. **(c)** contradict nontriviality with a Liouville theorem for
ancient Euler flows in the weighted class the scenario defines.

What is actually excluded is narrow and Seregin says so: the limit "cannot be irrotational";
`m < 1/2` is impossible so `m ≥ 1/2` (equivalently `m₀ ≥ 4/5` suffices);
arXiv:2507.08733 Thm 5.1's exclusion carries an *added Ladyzhenskaya–Prodi–Serrin
hypothesis* (`v ∈ L^{s₂,l₂}(Q)` with `3/s₂ + (α+1)/l₂ = α`), and Remark 5.2 says as much.
**No unconditional exclusion of Type-II for 3D NS exists in this literature, and none is
claimed.**

### 1.3 The rate-free branch — real, applicable, and far too weak to separate the types

Leray's lower bounds; Escauriaza–Seregin–Šverák 2003 (no `L^∞_t L³_x` blowup — no rate
assumed, so it applies to Type-II verbatim); Tao arXiv:1908.04958 (critical `L³` exceeds
`(log log log 1/(T−t))^c` along a sequence); Barker–Prange Result D (a slightly
supercritical Orlicz norm blows up, Leray–Hopf, no Type-I hypothesis). None of these
excludes Type-II — they exclude *bounded critical norm*, which any Type-II blowup violates
anyway. A triple logarithm is nowhere near a rate that distinguishes `(T−t)^{-1}` from
`(T−t)^{-1-δ}`.

Barker–Prange also state the boundary of the conditional theory, which is
Fact `fact:s51-ss` seen from the other side: *"the regularity in the Type I case is proved
only for axisymmetric flows"*.

### 1.4 Open ≠ expected, kept apart deliberately

- **Open.** Nothing read here excludes Type-II, unconditionally or under a hypothesis the
  equations discharge.
- **Not expected.** No source read asserts Type-II occurs. The strongest unforced numerical
  candidate is **Type-I-shaped**: Hou arXiv:2107.06509v2 / FoCM 2022, *"Both R(t) and Z(t)
  seem to scale like O((T−t)^{1/2}). The maximum vorticity grows like O((T−t)^{−1})."* Hou
  adds the honest caveat that the solution *"cannot be asymptotically self-similar"* so the
  rate *"is most likely not exactly equal to 1/2"*, and that the uniform upper bound *"would
  be almost impossible to verify numerically"* — i.e. the data cannot exclude a mild Type-II
  correction, but supply no positive indication of one. Seregin's four papers treat Type-II
  throughout as a possibility to be *ruled out*, never as an expectation.
- **What is excluded on the other side.** Backward self-similar blowup with finite energy is
  dead (Nečas–Růžička–Šverák; Tsai). Neither type has positive evidence; one of them has a
  specific mechanism excluded.

### 1.5 The forced Type-II profile, verified from the primary source

The handover asserts the OpenAI construction produced a Type-II profile. **Confirmed, from
the PDF rather than from the book chapter**, and with more structure than was on record.
Theorem 1.1 (for every `ν > 0` there is `f ∈ C_c^∞` and smooth `u` on `ℝ³×[0,1)` with
`u(·,0)=0`, `sup‖u‖_{L²} < ∞`, `limsup‖u‖_{L^∞} = ∞`), and §2 gives, with `τ = 1−t`:

```
ℓ_r ≍ τ^{1/2},   ℓ_z ≍ τ^{1/2−h},   0 < h < 1/100
|u_θ^{(0)}|, |u_z^{(0)}| ≍ τ^{−1/2−h},   |u_r^{(0)}| = O(τ^{−1/2})
```

Type-II in the temporal sense, yes. But the sharper point, not previously recorded here:
**the two length scales differ**, `ℓ_r/ℓ_z ≍ τ^h → 0`. A single-parameter parabolic
rescaling has *one* length; it cannot resolve a two-scale core at any amplitude. So the
obstruction of `rem:s39-distance`(1) is not only a rate mismatch — it is a *geometry*
mismatch, and the geometry mismatch is the one exhibited by the only construction in
existence. (Forced. Wall 9 governs. Quoted as a profile, not as evidence.)

---

## 2. Deliverable 2 — the vacuity, re-derived, and corrected

### 2.1 Confirmed as far as it goes

With `u_λ(y,s) = λu(x*+λy, T+λ²s)` and the Type-I rates:
`|u_λ| = λ|u| ≤ λC_I(λ²|s|)^{−1/2} = C_I|s|^{−1/2}` and
`|ω_λ| = λ²|ω| ≤ C_I|s|^{−1}`, both λ-independent — `eq:rescaled-bounds` is right and
Type-I is used exactly as the document says. With `M(t) ≤ C(T−t)^{−1−δ}`, `δ>0`, the
vorticity bound degrades by `λ^{−2δ} → ∞`. The uniform bound really is lost.

### 2.2 The derivation the document does not have (`lem:wpt1-nueff`)

For `v(y,s) = A·u(x*+Ly, T+τs)`, `p̃ = (Aτ/L)p(...)`:

```
∂_s v + (τ/AL)(v·∇)v − ν(τ/L²)Δv + ∇p̃ = 0
```

Normalising the nonlinearity by `τ = AL` gives NS with

```
ν_eff = ν · A/L .
```

Checks: `A=L=λ` gives `τ=λ²`, `ν_eff=ν` — that is `eq:rescaling`. Seregin's Euler scaling
`A=λ^α, L=λ` gives `τ=λ^{α+1}` and `ν_eff = νλ^{α−1} → 0` for `α>1`. His `f`-form
`v^λ = λf(λ)v(λy, λ²f(λ)τ)` gives `ν_eff = νf(λ) → 0`. All three consistent.

### 2.3 The trichotomy (`cor:wpt1-trichotomy`)

Requiring `L^∞` normalisation of a rate `sup|u| ≍ (T−t)^{−β}` gives `A τ^{−β} ≍ 1`, hence
`A = L^{β/(1−β)}` and `ν_eff/ν = A/L = L^{(2β−1)/(1−β)}`. As `L↓0`:

| β | ν_eff | limit equation |
|---|---|---|
| `< 1/2` | `→ ∞` | Stokes/heat — no blowup mechanism survives |
| `= 1/2` | `= ν` | **NS. `A=L`. This is `eq:rescaling`, and it is the UNIQUE member of the family that preserves NS.** |
| `(1/2, 1)` | `→ 0` | **Euler.** Seregin's zoom; `β = α/(α+1)`, and his `α ∈ (1,3/2]` is `β ∈ (1/2,3/5]` |
| `≥ 1` | — | **no member of the family works.** At `β=1` normalisation forces `L=1` (no spatial zoom); for `β>1` it forces `τ = AL → ∞`, contradicting `τ ≤ T`. |

The last row is a **scope limit on the Euler-zoom method itself** and is stated in no source
read here. Velocity blowup at `(T−t)^{-1}` or faster is outside *both* methods.

### 2.4 The correction to `rem:s39-distance`(1) — reported prominently, as the handover asked

The handover explicitly asked: if the vacuity is less total than stated, say so loudly.
**It is less total, in two independent ways, and the true statement is worse for the
programme, not better.**

- **(a) Within the same rescaling.** A limit object needs uniform bounds in *some* topology,
  not in `L^∞`. Under the strictly weaker hypothesis `g < ∞`, the same λ-family is bounded
  in `L^∞_t L²_loc ∩ L²_t Ḣ¹_loc` on every fixed cylinder (`A(u_λ,a) = A(u,λa) ≤ g`, ditto
  `E`, `C`), and with the pressure controlled the classical CKN/Lin compactness gives an
  ancient **suitable weak NS** solution — viscosity retained. **The document already knows
  this**: the remark following `prop:profile-extraction` says the `C^∞_loc` topology is
  *"far more than the weak-\* L^∞ / weak L²_loc convergence one gets from energy bounds
  alone"*. `rem:s39-distance`(1) simply does not carry the qualification. This is region
  (R2).
- **(b) Outside it.** Even when `g = ∞`, the Euler zoom produces a nontrivial ancient limit.

**Correct statement, for a future session to quote instead:**

> A rescaling that normalises a supercritical blowup rate in `L^∞` necessarily has
> `ν_eff = νA/L → 0`, so its limit object solves **Euler**. The vacuity is not that nothing
> survives; it is that what survives is **inviscid**, and every tool of §33–§44 is parabolic
> — Moser iteration, the magnitude Caccioppoli inequality, De Giorgi level iteration,
> `lem:local-enstrophy-typeI`, CKN ε-regularity, ESS backward uniqueness. None of them says
> anything about a solution of the Euler equations.

`rem:s39-distance`(1) remains correct about the *conclusion* and imprecise about the
*reason*.

### 2.5 The one loophole, closed by dimensions (`cor:wpt1-swirl`)

The obvious hope after §2.3 is: find an Euler zoom that still carries something this
programme controls. In the axisymmetric class the natural candidate is the swirl invariant
`Γ = r·u_θ`, which is free by the maximum principle (`lem:s51-gamma`). It does not work,
and the reason is dimensional:

```
Γ_new(y,s) = |y'| v_θ = (|x'|/L)(A u_θ) = (A/L) Γ = (ν_eff/ν) Γ
```

`Γ` and `ν` both have dimension `L²/T`, so they rescale by the *same* factor and `Γ/ν` is
invariant across the whole family. Hence:

1. Any zoom with `ν_eff → 0` has `Γ → 0`: **the Euler-zoom limit of an axisymmetric Type-II
   blowup is swirl-free.** This reproduces, from dimensions alone, the computation of
   arXiv:2402.13229 (*"The latter means that our limiting Euler equations have no swirl,
   i.e., u_ϑ vanishes"*). No novelty is claimed for the conclusion, only for the reading.
2. Conversely, **any rescaling that retains the swirl retains the viscosity** — i.e. is
   `eq:rescaling`, i.e. requires temporal Type-I. "An Euler zoom that keeps the swirl" is
   *dimensionally impossible*, not merely unknown.

Recorded alongside, with an explicit non-claim (`rem:wpt1-s52-coincidence`): S52's
`prop:s51-aligned` found that alignment forces `∂_zΓ ≡ 0`; the Euler zoom forces `Γ ≡ 0`.
These are different in strength and arise for unrelated reasons. **No connection is
asserted.** They are logged together only so a future session does not conflate them.

### 2.6 Why the other-PDE techniques do not transfer

The handover asked specifically about Giga–Kohn / Merle–Zaag-type modulation.

- *Semilinear heat.* Giga–Kohn similarity variables characterise blowup under a Type-I
  bound; below the Sobolev exponent all blowup *is* Type-I, and Type-II appears only at and
  above criticality. Its construction stabilises an explicit **bubble** — Collot–Merle–%
  Raphaël (arXiv:1709.04941) construct *anisotropic* Type-II blowup for `∂_t u = Δu + u^p`,
  `p ≥ 3`, `d ≥ 14`, by "construction and stabilization of type II bubbles in the parabolic
  setting". The bubble is the stationary ground state.
- *Geometric flows.* Type-II MCF singularities are handled by Hamilton point selection,
  whose limit is an **eternal** solution identified as a translating soliton via the
  equality case of a differential Harnack inequality.

Both need something NS does not supply: an explicit stationary/self-similar profile to
linearise around, or a differential Harnack. The candidate NS profile is a backward
self-similar solution, and those are excluded with finite energy (Nečas–Růžička–Šverák).
**This is a statement about what is missing, not a proof that no substitute exists** — and
Seregin's programme is precisely the admission of it: it gives up on keeping the limit in
the same equation and pays with an Euler limit.

---

## 3. Deliverable 3 — verdict

**(iii), justified by (ii). Not (i).**

1. **Negative, derived.** Type-II is not reachable by any adaptation of the machinery of
   §33–§38. `cor:wpt1-trichotomy` is the obstruction: `ν_eff = νA/L`, and `L^∞`-normalising
   a supercritical rate forces `ν_eff → 0`. Adapting the rescaling does not adapt the tools;
   it deletes the term they are built on. `cor:wpt1-swirl` closes the axisymmetric loophole.
2. **The branch that owns the question, named.** *Liouville theorems for ancient solutions
   of the 3D Euler equations in weighted classes defined by a blowup scenario.* Read, in
   order: arXiv:2304.04045, arXiv:2402.13229, arXiv:2507.08733, arXiv:2606.29468. The
   axisymmetric one is closest to `sec:axisym-s51`, and there the residual problem is a
   Liouville theorem for **swirl-free** ancient axisymmetric Euler. This programme has no
   foothold: it has never written an Euler estimate, and `rem:s51-torus` already established
   that entering the axisymmetric class requires first moving the standing framework off 𝕋³.
3. **Entry cautions, so a future session prices it correctly.**
   (a) Producing another scenario-conditional exclusion is producing another conditional
   statement, not progress on the Clay problem — the same trap as `prop:s51-boundedness`.
   (b) The Euler-zoom family stops at `β = 1`.
   (c) Region (R2) is outside both methods.
4. **What does not change.** Wall 3 keeps its classification `open`. Nothing here is
   progress on Type-II, Type-I, Wall 1, or any Prize or Claim A conclusion.

**No WP-T2 is proposed.** The handover's own instruction applies: do not manufacture a
target to have something to hand off.

---

## 4. Verification record

**Derived from scratch, taken from no source:**
- `ν_eff = ν A/L` and the four-way trichotomy, including the `β ≥ 1` scope limit. Checked
  against `eq:rescaling` (`A=L` ⟹ `ν_eff=ν`), against Seregin's `λ^α` form
  (`ν_eff = νλ^{α−1}`) and against his `f(λ)` form (`ν_eff = νf(λ)`). All three agree.
- `Γ_new = (A/L)Γ = (ν_eff/ν)Γ`, and the dimensional reason (`[Γ] = [ν] = L²/T`). Checked
  against the displayed identity in arXiv:2402.13229 — agrees.
- The vacuity itself, from `eq:rescaling` and the Type-I rates, both directions
  (uniform for `β=1/2`, `λ^{−2δ}` blow-up for `β>1/2`).

**Read as primary text, not as summaries:** Seregin arXiv:2304.04045 (full PDF),
arXiv:2402.13229v3, arXiv:2507.08733, arXiv:2606.29468 (full PDF); Barker–Prange
arXiv:2211.16215 (full PDF, incl. footnote 8 and §7); Hou arXiv:2107.06509v2 (full PDF, incl.
§3.4 and Concluding Remarks); the OpenAI manuscript PDF (Thm 1.1 and §2). Abstracts fetched
for Albritton–Barker arXiv:1811.00502, Collot–Merle–Raphaël arXiv:1709.04941, Hou
arXiv:2405.10916, Cheskidov–Dai–Palasek arXiv:2511.09556.

**A search-engine summary that is wrong, recorded so it is not repeated.** One search
synthesis asserted that Nečas–Růžička–Šverák and Tsai "essentially rules out any reasonable
Type II blowup scenario". That is backwards: NRS/Tsai exclude *backward self-similar*
(i.e. the cleanest **Type-I**) blowup. Nothing in this session was taken from search
syntheses.

**Checked against the handover, per standing discipline:**

| Handover claim | Verdict |
|---|---|
| `rem:s39-distance` at line 7703 | **Correct** |
| `subsec:wall3-typeII` at line 12668 | **Correct** (`\subsection` on 12668, `\label` on 12669) |
| Type-I/Type-II defined by `limsup (T−t)M(t)` | **Correct for this document** (`line 5565` carries both the vorticity and velocity rates), but the handover's instruction to check the sources was the right one: the Type-II literature uses a *third* convention (§1.1), and the handover does not anticipate that. |
| The OpenAI profile is Type-II with `u_θ ~ τ^{−(1/2+h)}`, `h>0` | **Correct**, verified from the primary PDF, and it is also **anisotropic** (`ℓ_r ≍ τ^{1/2}` vs `ℓ_z ≍ τ^{1/2−h}`), which the handover does not mention and which is the sharper obstruction. |
| "No session has ever attempted … even a literature survey" of Type-II | **Correct** |
| "rescale to nothing" is a document assertion, not re-derived | **Correct, and it turns out to be too strong** — see §2.4. |
| S53 is the next session tag | **Correct**, verified twice. |

**Defect found in the document, in passing.** An automated audit of every `(line~N)`
citation whose target bears a resolvable label:

| drift | count | examples |
|---|---|---|
| `+33` | 12 | `rem:s39-distance` (cited 7670, actual 7703 — this is the Wall 3 citation), `prop:s41-nocontraction` (8156/8189), `prop:s41-vacuous` (8191/8224), `prop:s41-coarea` (8270/8303), `rem:s41-shortfall` (8395/8428), `prop:s44-linfty-route` (9209/9242) |
| `−655` | 1 | `lem:typeI-scaleinv-bounds` (cited 7384, actual 6729) |
| `−10` to `−18` | 4 | `sec:moser-s33` (5565/5547), `lem:logtime` (7408/7391), `thm:apriori-sigma` (5883/5873) |
| `+11` | 1 | `obs:s48-wall1` (10408/10419) |
| exact | 8 | |

Every `\ref` resolves, so nothing is broken; the line numbers are decoration that has gone
stale. **Not repaired here** — and the new section was inserted *before* `sec:walls-s47`
precisely so that all of them keep the values they have. A session that repairs them should
do so as its own change, with the audit script re-run afterwards.

**No error was found in my own derivations** on re-check, but two things were softened
after first drafting: the claim that region (R2) is non-empty (downgraded to "whether it is
empty is unknown"), and the claim that the dimensional reading of `cor:wpt1-swirl` is new
(downgraded to "no novelty is claimed for the conclusion, only for the reading").

---

## 5. Files changed

```
proofs/claim_a_3d_proof_attempt.tex   +666 / −0   (sec:typeII-wpt1, 12 bib entries, walls rows 3'/3''/3''')
PROGRESS.md                           KEY FINDING S53 block, session-log row, Active milestone
SESSION-LOG/2026-09-12-S53-typeII-reconnaissance.md   (this file)
```

Not committed, not pushed. Worktree `/home/user/wpt1-worktree`, branch
`claude/wpt1-typeii-recon`, for independent verification before merge. Re-check the session
tag and the tip at merge time — S52's log records that this exact collision has now happened
three times running.
