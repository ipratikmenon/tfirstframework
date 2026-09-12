# S54 — WP-T2: does the OpenAI profile survive as an ancient Euler solution?

**Date:** 2026-09-12
**Branch / worktree:** `claude/wpt2-openai-profile` at `/home/user/wpt2-worktree`, based on
`0a12970` ([S53] plan: WP-T2 handover).
**Scope:** `proofs/claim_a_3d_proof_attempt.tex`, `PROGRESS.md`, `SESSION-LOG/` only.
No numerics written or run (a symbolic re-derivation was done in the scratchpad and is
not committed). `book/`, `layer3/`, `layer4/`, `results/` untouched.
**Deliverable:** `\S sec:wpt2-openai` (§40 in the rendered PDF), 2 new bibliography entries.
**Build:** `check_proof.py` exit 0 (686 labels, 475 `\ref` targets, 55 `\cite` keys, all
resolved). `pdflatex` ×3: **0 errors, 0 LaTeX warnings, 109 hyperref bookmark warnings,
64 overfull / 2 underfull hboxes — every count identical to baseline.** The baseline was
**rebuilt from the tip this session** rather than taken from S53's log; it agreed with
S53's recorded 109. 196 pp (was 188). Diff purely additive.

**Session tag:** verified before writing anything. `git log --oneline -1` = `0a12970`,
`ls SESSION-LOG/` highest = S53, so S54 is correct. Re-checked before writing
PROGRESS.md and the log: unchanged, no other agent branch present. **Label prefix is
`wpt2`, not `s54`**, following S53's reasoning — a work-package prefix cannot collide.

---

## 0. One line per deliverable

| Deliverable | Answer |
|---|---|
| **D0** — get the real profile from a primary source | Source found, fetched, read. The screenshot transcription this programme carried is **correct in every symbol**. The governing profile equations **are stated explicitly** — the favourable branch. One fact was missing and it turned out to govern everything: **the leading balance retains radial viscosity.** |
| **D1** — derive the anisotropic rescaling | A well-**defined**, nontrivial, ancient limit exists — and it is **not Euler**. `ν_eff` splits into two; the radial one is *exactly* `ν`. The limit is a **Prandtl-type anisotropic system**. Not claimed well-**posed**. |
| **D2** — test against ancient-solution exclusion | Outcome **(i), excluded** — but by the **swirl maximum principle**, not by any Liouville theorem, and the Liouville route provably could not have worked. |
| **Verdict** | A narrow positive. One profile family, axisymmetric class, unforced equations. **Not** a step toward (A)/(B). **Not** a statement about the OpenAI theorem. |

---

## 1. Deliverable 0 — the primary source

**Found and fetched.** `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`
— *Finite Time Blowup for Navier–Stokes*, OpenAI, 165 pp + appendices. The S53
bibliography entry carried no locator, so the document could not be re-checked from it;
that is now repaired in `\bibitem{OpenAI2026}`.

**Theorem 1.1, verbatim**, is as S53 recorded: for every `ν>0` there are
`f ∈ C_c^∞(ℝ³×(0,∞);ℝ³)`, a compact `K`, and smooth `u,p` on `ℝ³×[0,1)` solving forced NS
with `u(·,0)=0`, supports in `K`, `sup‖u(t)‖_{L²} < ∞`, `limsup‖u(t)‖_{L^∞} = ∞`. It claims
Fefferman alternative (C), and (D) via Corollary 10.6.

### 1.1 The transcription was right

The screenshot-derived formulas in the handover are confirmed **symbol for symbol** by the
source's (4.1) and (4.3):

```
A = 1/2 + h,  D = 1/2 - h,  0 < h < 1/100
tau = q(1 - eta^2),  z = q^D eta,  s = r^2/2,  X = s/q = r^2/(2q),  q>0, -1<eta<1
u_theta^(0) = q^-A E,  u_z^(0) = q^-A U,  r u_r^(0) = V_0,  p^(0) = q^-2A Pi
E = C^-1 sqrt(2X) phi
```

Profile domain: the closed rectangle `[0,X_c]×[-1,1]` (elsewhere `[0,X_an]×[-1,1]`), with
smoothness imposed on `phi` (equivalently `F = phi/C`), and `V_0 = X v_0`. **Nothing built
on the transcription needs revision.** This was the specific thing the handover flagged as
unverified, and it survived.

### 1.2 The profile equations are explicit — the favourable branch

The handover asked whether the source states governing profile equations or only their
asymptotic consequence. **It states them.** Proposition 4.2 gives

```
-2L (X phi_XX + 2 phi_X)/phi = S_q ,     -2L (X U_XX + U_X) = S_n        (4.13)
```

closed by the incompressibility/centrifugal identities

```
V_0 = (X/L)(2 eta U - 2 D eta A_X(U) - d d_eta A_X(U)) ,   Pi_X = E^2/(2X)  (4.7)
```

with `d = 1-eta^2`, `L = 1-2h eta^2`, `A_X(f) = X^-1 int_0^X f dx`, and the explicit
first-order sources `S_q, S_n` of (4.9). So Deliverable 2 was a test of an explicit object,
not the harder reconstruction task the handover warned about.

Also recorded: the construction runs at **`ν = 1`** throughout ("We use ν = 1 for the rest
of this outline"), with general `ν` recovered by `u_ν(x,t) = √ν u(x/√ν, t)`, singular time
unchanged.

### 1.3 The fact that was missing, and that governs everything after it

**The leading balance retains radial viscosity and omits axial viscosity.** The source says
this three separate ways:

- §2.1 dimensional analysis: `|u_r|/ℓ_r = O(τ^-1)`, `|u_z|/ℓ_z ≍ τ^-1`, `ν/ℓ_r² ≍ τ^-1`,
  and `(ν/ℓ_z²)/(ν/ℓ_r²) = ℓ_r²/ℓ_z² ≍ τ^{2h} → 0`;
- `Re_θ = |u_θ|ℓ_r/ν ≍ τ^{-h} → ∞` but **`Re_r = |u_r|ℓ_r/ν = O(1)`**, of which it says
  "viscosity continues to compete with radial inflow";
- (4.12) defines the leading tangential residual by *adding back* `∂_z² u_j^(0)`, and the
  notation guide describes the resulting stress as "retaining radial viscosity and omitting
  axial viscosity".

S53 did not have this, and it is the reason Deliverable 1 does not come out the way the
handover expected.

---

## 2. Deliverable 1 — the anisotropic rescaling, derived

### 2.1 The general lemma (`lem:wpt2-aniso`)

Derived by direct substitution into the **cylindrical** axisymmetric equations, not by
analogy to `lem:wpt1-nueff`, and independently re-checked by computer algebra. With
`r = L_r ρ`, `z = L_z ζ`, `t = T + T_s s`, `v_r = A_r u_r`, `v_θ = A_θ u_θ`, `v_z = A_z u_z`,
`q̃ = A_p p`:

Continuity becomes `(1/ρ)∂_ρ(ρ v_r) + (A_r L_r)/(A_z L_z) ∂_ζ v_z = 0`, so

```
      the limit is incompressible  <=>  A_r L_r = A_z L_z =: Lambda
```

and normalising the advection by `T_s = Lambda` gives **two** effective viscosities:

```
      nu_r = nu * A_r / L_r ,      nu_z = nu * A_z / L_z
```

with the radial pressure gradient carrying `A_r²/A_p`, the centrifugal term `A_r²/A_θ²`,
and the axial pressure gradient `A_z²/A_p`.

**Reduction check.** `A_r=A_θ=A_z=A`, `L_r=L_z=L` makes the constraint vacuous, gives
`T_s = AL` (S53's `τ=AL`), collapses both viscosities to `ν A/L`, and with `A_p = A²`
(matching S53's `p̃ = (Aτ/L)p = A²p`) puts the pressure and centrifugal coefficients at 1.
**`lem:wpt1-nueff` is contained in `lem:wpt2-aniso`.**

### 2.2 The limit for the OpenAI exponents (`prop:wpt2-limit`)

`L_r=τ^{1/2}`, `L_z=τ^D`, `A_r=τ^{1/2}`, `A_θ=A_z=τ^A`, `A_p=τ^{2A}`.

1. **Incompressibility is automatic, and it is the source's own identity.**
   `A_r L_r = τ` and `A_z L_z = τ^{A+D} = τ` — because **`A + D = 1`**, which the source
   uses explicitly to derive the first half of (4.7). The profile's amplitude anisotropy is
   exactly conjugate to its length anisotropy. `T_s = τ`: the natural time unit is the
   remaining time.
2. **`ν_r = ν τ^{1/2}/τ^{1/2} = ν` — exactly preserved, independent of τ.**
   **`ν_z = ν τ^A/τ^D = ν τ^{2h} → 0`.** This reproduces the source's own dimensional
   statement from the rescaling alone.
3. **Radial momentum degenerates.** Its pressure and centrifugal coefficients are both
   `τ^{1-2A} = τ^{-2h} → ∞` while the axial pressure coefficient is `1`. Multiplying by
   `τ^{2h}` kills every inertial and viscous term and leaves the **cyclostrophic balance**
   `∂_ρ q̃ = v_θ²/ρ` — which is precisely `Π_X = E²/(2X)`, the second half of (4.7).
4. **The limiting system** is that constraint plus swirl and axial transport with
   radial-only diffusion plus continuity. It is closed (`v_r` from continuity, `q̃` from the
   constraint up to a free `q̃_0(ζ,s)`) and nontrivial (the rescaled leading field converges
   to the profile itself, `v_θ = |s|^{-A}E` with `X → ρ²/(2|s|)`, `η → ζ|s|^{-D}`).
5. **So the limit does not solve Euler**, and `cor:wpt1-trichotomy`(iii) does not transfer.

### 2.3 What the limit is, and the honest answer on well-posedness

Structurally it is a **boundary-layer (Prandtl-type) reduction**: thin-direction momentum
degenerated to a constraint; thin-direction velocity recovered from continuity; diffusion
only in the thin direction; no diffusion in the advected direction. The thin direction here
is the **radial** one (`ℓ_r/ℓ_z ≍ τ^h → 0`), and the free outer function is visible in the
source as the axis pressure datum `Π(0,η)` of its (4.31), entering `Π = Π(0,η) + C_p`
exactly as an outer pressure gradient enters Prandtl's system.

Two cautions, both recorded in the document:

- **No source read here writes this system.** The classical "quasi-cylindrical"
  approximation for slender vortices is an *Euler* (Bragg–Hawthorne) reduction and is not
  this — I fetched González et al. (arXiv:1203.2787) and checked. The closest verified
  statement is that the inner viscous core of a slender vortex is treated by Hall (1961) and
  Stewartson–Hall (1963) "by making scaling assumptions of the usual boundary-layer type",
  as reported in Mayer, *Long's vortex revisited* (arXiv:1303.1212), which I did read. The
  identification is therefore **structural**, made from the derived equations, and is offered
  as such.
- **Well-defined ≠ well-posed.** Gérard-Varet–Dormy (JAMS 23:591–609, 2010; arXiv:0904.0434,
  abstract fetched and quoted) prove Prandtl is **linearly ill-posed in Sobolev spaces**,
  "the strong instability is due to viscosity". Their result is *not* transported to this
  system and no claim is made that it applies; it is recorded because the features their
  proof exploits — viscosity in the thin variable, none in the advected one — are exactly
  what `prop:wpt2-limit` produces. **The handover asked whether the limit is well-posed; the
  answer is "well-defined, not known to be well-posed, with specific reason to doubt it".**

### 2.4 The Euler zoom fails both ways (`cor:wpt2-isotropic-dichotomy`)

The Liouville theorems need an *isotropic* Euler zoom. Since `sup|u| ≍ τ^{-A}` with
`A = 1/2+h ∈ (1/2,1)`, `cor:wpt1-trichotomy`(iii) formally applies. But one length cannot
match two:

- **`L ≍ ℓ_r = τ^{1/2}`**: `ν_eff = ν τ^h → 0` (Euler, yes), rescaled radial extent `O(1)`,
  but rescaled **axial extent `ℓ_z/L = τ^{-h} → ∞`**, so the limit is ζ-independent;
  `A|u_r| ≍ τ^h → 0`, so the radial component is annihilated; and `T_s/τ = τ^h → 0`, so the
  limit is **steady**. The limit object is a steady ζ-independent columnar flow
  `v_θ(ρ)e_θ + v_z(ρ)e_z` — an exact steady Euler solution for *any* profiles, nontrivial but
  **of infinite energy** and in no weighted class any of the four Seregin theorems uses.
- **`L ≍ ℓ_z = τ^D`**: rescaled radial extent `ℓ_r/L = τ^h → 0`, the core collapses onto the
  axis, the limit is **trivial**, and nontriviality — the hypothesis every Liouville argument
  contradicts — fails outright.

---

## 3. Deliverable 2 — the exclusion, which needed none of §2

### 3.1 The result (`prop:wpt2-swirl-exclusion`)

From (4.3), with `r = √(2qX)` and `E = C^{-1}√(2X)φ`:

```
Gamma^(0) = r u_theta^(0) = sqrt(2qX) * q^-A * C^-1 sqrt(2X) phi
          = (2 X phi / C) q^(1/2 - A) = q^-h H(X,eta),     H = sqrt(2X) E = 2 X phi / C
```

**This is the source's own identity**, stated verbatim in its proof of Proposition 4.2:
"The angular equation is simplest in angular momentum `r u_θ^(0) = q^{-h} H`." My derivation
reproduces it independently. Since `q ≍ τ` on the core and `H ≍ 1` at the core edge
`X ≍ X_b`, `sup|Γ^(0)| ≳ τ^{-h} → ∞`. Confirmed a second, independent way from the source's
exterior swirl `K = r^{-1-2h}H_ext(τ/r²)`: `rK = r^{-2h}H_ext(τ/r²) ≍ τ^{-h}` at
`r ≍ τ^{1/2}`, matching at the join.

But `lem:s51-gamma` (S52, classical) gives `‖Γ(t)‖_∞ ≤ ‖Γ(0)‖_∞` for smooth axisymmetric
solutions of the **unforced** equations, and Schwartz data have
`‖r u_θ(·,0)‖_∞ ≤ ‖|x|u_0‖_∞ < ∞`. **Contradiction.** The profile cannot be the blowup
asymptotics of an unforced axisymmetric Navier–Stokes solution.

### 3.2 The geometric form (`cor:wpt2-radius-cap`)

The maximum principle gives `|u_θ| ≤ Γ_0/r` — i.e. the swirl component *automatically*
satisfies the **spatial** Type-I bound (S53's convention **(S)**). Hence on a core of radius
`ℓ_r`, `|u_θ| ≲ Γ_0/ℓ_r`. So:

> **An unforced axisymmetric core with the Type-I radial geometry `ℓ_r ≍ τ^{1/2}` can only
> carry Type-I swirl `|u_θ| ≲ τ^{-1/2}`. Type-II swirl `|u_θ| ≍ τ^{-1/2-h}` requires the
> strictly thinner core `ℓ_r ≲ τ^{1/2+h}`.**

The OpenAI profile has `ℓ_r ≍ τ^{1/2}` *together with* `|u_θ| ≍ τ^{-1/2-h}`, and violates the
cap by exactly the factor `τ^{-h}` that makes it Type-II.

### 3.3 Exactly what is and is not excluded

- **Axisymmetric class only** — `lem:s51-gamma` is. Nothing is said about a
  non-axisymmetric flow with the same scales.
- **Unforced only.** That is the point of the exercise, but it also means this is **not a
  statement about `\cite{OpenAI2026}`**: their flow is forced *and* **non-axisymmetric** (the
  oscillatory pulses carry nonzero azimuthal harmonics — "their cylindrical velocity
  components have zero angular average"), so `lem:s51-gamma` never applied to it. **No
  contradiction with their Theorem 1.1 is asserted or implied.** Wall 9 governs.
- **One profile family, not Type-II.** Any axisymmetric core with `ℓ_r ≲ τ^{1/2+h}` escapes
  `cor:wpt2-radius-cap` entirely. Region (R2), Wall 1 and Wall 3 are untouched.
- **Elementary.** The swirl maximum principle is classical and was already in this document.
  The content is the *check*. No credit is claimed for machinery.

---

## 4. Verdict, in the narrow terms the handover required

**Deliverable-2 outcome (i): the profile is excluded.** But:

1. It is excluded in the **axisymmetric** class, for the **unforced** equations, and it is
   **one profile family**.
2. It is **not** a step toward alternative (A) or (B). The chain to global regularity would
   still need Wall 1 closed on the Type-I side, *plus* a general Type-II profile-extraction
   theorem, *plus* a general unconditional Liouville theorem — three open problems WP-T2 does
   not touch. S53's finding that Seregin's programme is scenario-conditional throughout is
   unchanged.
3. It is **not** a statement about the OpenAI construction.
4. **The transferable lesson, and the only thing here worth reusing:** the Liouville
   programme was the wrong tool for this question, and two independent computations say so —
   the correct rescaling does not produce an Euler solution, and the wrong one produces either
   nothing or an infinite-energy shear flow. What settled it was a conservation law the
   axisymmetric class supplies for free. **Before importing profile-extraction machinery
   against a candidate profile, check it against `cor:wpt2-radius-cap` first.**

**No WP-T3 is proposed.** The out-of-scope Burgers-vortex/comparison-principle idea was not
attempted and is not merged in.

---

## 5. Verification record

**Re-derived from scratch, trusted from nothing:**

- `lem:wpt1-nueff` (S53). Every substitution re-checked: `∂_s v = Aτ ∂_t u`,
  `(v·∇_y)v = A²L(u·∇)u`, `Δ_y v = AL²Δ_x u`, `∇_y p̃ = Aτ∇_x p`; hence
  `ν τ/L² = νA/L` at `τ=AL`. **Correct.**
- `cor:wpt1-trichotomy` (S53). `Aτ^{-β} ≍ 1` with `τ=AL` gives `A^{1-β}=L^β`,
  `A = L^{β/(1-β)}`, `ν_eff/ν = L^{(2β-1)/(1-β)}`; all four cases including the `β≥1` scope
  limit (`β=1 ⟹ L=1`; `β>1 ⟹ τ = L^{1/(1-β)} → ∞`). **Correct.**
- `lem:wpt2-aniso` and `prop:wpt2-limit` — new, derived by hand from the cylindrical
  equations and then **independently re-derived with sympy** (scratchpad only, not
  committed). Hand and symbolic results agree term for term.
- `Γ^(0) = q^{-h}H` — derived from (4.3) and found to match the source's stated identity
  exactly; cross-checked against the exterior swirl formula (3.5).

**Read as primary text, not as summaries:** the OpenAI manuscript PDF in full for §1, §2,
§3.1, §4.1–4.2 (incl. (4.1), (4.2), (4.3), (4.5), (4.6), (4.7), (4.8), (4.9), (4.11), (4.12),
(4.13), (4.15), (4.30), (4.31), (5.1), (5.2), (3.5)) and its notation guide; Gérard-Varet–Dormy
arXiv:0904.0434 (abstract); Mayer arXiv:1303.1212 (full text, for the Hall/Stewartson–Hall
boundary-layer statement); González et al. arXiv:1203.2787 (full text — checked and
**rejected** as a citation, its quasi-cylindrical approximation is inviscid).

**Checked against the handover and S53, per standing discipline:**

| Claim checked | Verdict |
|---|---|
| Handover's screenshot transcription of the profile | **Correct in every symbol** — confirmed against (4.1)/(4.3). |
| S53's `lem:wpt1-nueff`, `cor:wpt1-trichotomy` | **Correct**, re-derived independently. |
| S53's `fact:wpt1-openai` scaling data | **Correct**, re-verified against §2 of the source. |
| Handover: "if the source states only asymptotics…" | **Does not apply** — (4.13) is explicit. |
| Handover's expectation that the limit would be Euler-like | **Wrong, and it warned it might be.** `ν_r = ν` exactly; the limit is Prandtl-type. |
| S53 §2.5 / `cor:wpt1-swirl`(i) as phrased in the S53 log | **Needs its unforced hypothesis attached.** See below. |
| S53's "the geometry mismatch is the sharper obstruction" | **True for rescaling, irrelevant to exclusion.** See below. |
| Handover's guess that this is S54 | **Correct**, verified twice (start and before writing). |
| Baseline of 109 bookmark warnings | **Correct**, but rebuilt from the tip rather than trusted. |

**Three refinements to S53, recorded in `rem:wpt2-s53-corrections`:**

1. **`cor:wpt1-swirl`(i) has a load-bearing hypothesis, and this profile witnesses it.** The
   corollary concludes "the Euler-zoom limit of an axisymmetric Type-II blowup is swirl-free"
   from `Γ_new = (ν_eff/ν)Γ` *together with* `‖Γ(t)‖_∞ ≤ ‖Γ(0)‖_∞`. The second ingredient is
   `lem:s51-gamma`, which requires the **unforced** equations. For this profile
   `Γ ≍ τ^{-h} → ∞`, so `Γ_new = τ^h · τ^{-h} = O(1) ≠ 0`: the Euler-zoom limit of *this*
   object **retains its swirl**. The corollary as written in the proof document is correct
   (it cites `lem:s51-gamma`); the **S53 session log's §2.5 states it without the
   qualification**, and that phrasing should not be reused.
2. **`cor:wpt1-trichotomy`(iii) is stated without an anisotropy caveat** and could be read as
   more general than it is. `prop:wpt2-limit` shows the conclusion changes qualitatively for
   two lengths and two amplitudes.
3. **The anisotropy is the obstruction to *rescaling*, not to *exclusion*.**
   `fact:wpt1-openai` calls `ℓ_r/ℓ_z → 0` "the sharper point" and is right —
   `cor:wpt2-isotropic-dichotomy` is its quantitative form. But it is not what kills the
   profile: `cor:wpt2-radius-cap` is single-scale and would exclude an *isotropic*
   `ℓ_r = ℓ_z = τ^{1/2}`, `|u_θ| ≍ τ^{-1/2-h}` profile just as fast.

**Nothing stale found in the handover** beyond the two items above; its session-tag guess
(S54) held, and its warning that the transcription was unverified was the right instruction
even though the transcription turned out to be right.

**Softened after first drafting, deliberately:** the identification of the limit with a named
system (downgraded from "this is the quasi-cylindrical system" to "structurally a
Prandtl-type reduction; no source read writes this system"), and the well-posedness claim
(downgraded from "well-posed" to "well-defined, not known to be well-posed").

---

## 6. Files changed

```
proofs/claim_a_3d_proof_attempt.tex   +632 / -5   (14,058 -> 14,685 lines)
      the 5 deletions are the OpenAI2026 bibliography annotation, replaced by a
      longer one carrying the URL locator and the expanded reading record;
      everything else is purely additive: sec:wpt2-openai + 2 new bib entries.
PROGRESS.md                           KEY FINDING S54 block, session-log row, Active milestone
SESSION-LOG/2026-09-12-S54-wpt2-openai-profile.md   (this file)
```

Not committed, not pushed. Worktree `/home/user/wpt2-worktree`, branch
`claude/wpt2-openai-profile`, for independent verification before merge. Re-check the session
tag and the tip at merge time.
