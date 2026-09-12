# S55 — WP-T3: a quantitative floor on the forcing in Wall 9

**Date:** 2026-09-12
**Branch / worktree:** `claude/wpt3-forcing-floor` at `/home/user/wpt3-worktree`, based on
`1aa483e` ([S54] plan: WP-T3 handover).
**Scope:** `proofs/claim_a_3d_proof_attempt.tex`, `PROGRESS.md`, `SESSION-LOG/` only.
No numerics written or run. `book/`, `layer3/`, `layer4/`, `results/` untouched.
**Deliverable:** `\S sec:wpt3-forcing-floor` (§41 in the rendered PDF), plus a quantitative
addendum to Wall 9, a row 9′ in the status table, and an extended `OpenAI2026` reading
record. No new bibliography entries.
**Build:** `check_proof.py` exit 0 (714 labels, 500 `\ref` targets, 55 `\cite` keys, all
resolved). `pdflatex` ×3: **0 errors, 0 LaTeX warnings, 109 hyperref bookmark warnings,
64 overfull / 2 underfull hboxes — every count identical to baseline.** The baseline was
**rebuilt from the tip this session** rather than taken from S54's log; it agreed with
S54's recorded 109 / 64 / 2 / 196 pp exactly. 206 pp (was 196). Diff purely additive apart
from the `OpenAI2026` annotation.

**Session tag:** verified before writing anything and re-checked before writing
`PROGRESS.md` and this log. `git log --oneline -1` = `1aa483e` both times;
`ls SESSION-LOG/` highest = S54 both times (note S49 and S50 have no log files; S48 is
followed directly by S51). **S55 is correct — the handover's guess held this time.**
No other agent branch or worktree present. **Label prefix is `wpt3`, not `s55`**,
following S53/S54: a work-package prefix cannot collide.

**Live file verified before touching anything:** `proofs/claim_a_3d_proof_attempt.tex` was
14,685 lines at the tip, matching both the live working copy and S54's recorded figure.

---

## 0. One line per deliverable

| Deliverable | Answer |
|---|---|
| **D1** — the forced maximum principle as a standalone lemma | **Proved**, as two lemmas. The forced swirl equation gains **exactly** `r·f_θ` and nothing else. The maximiser is automatically off-axis, so **no ε-exhaustion of the axis is needed** — cleaner than `lem:s51-gamma`'s own proof. Danskin one-sided derivatives written out in full; an integrated form needing no differentiability anywhere. **The handover's own statement of the target was wrong (sign).** |
| **D2** — apply it, get a number | `\|f_θ\| ≥ (h·c_*/c_r)·τ^-(A+1) = τ^-(3/2+h)` pointwise at the running swirl maximum; `\|r·f_θ\| ≥ h·c_*·τ^-(A+1/2)`. Both constants explicit. Full leading order in the θ-momentum equation. |
| **D3** — check against the primary source | **Clean confirmation, exponent-exact, no slack.** The source states `T = q^{-A-1/2}T_0`; converting gives `\|f_θ\| ≍ τ^-(A+1)` — the same two exponents, by a genuinely different route. Constants not comparable from the source. |
| **Verdict** | **Documentation-sharpening of Wall 9, and nothing else.** Not a step toward (A)/(B). Not a challenge to (C)/(D). No error in the OpenAI construction asserted, and none found. |

---

## 1. What was re-derived rather than trusted, and what turned out to be stale

Per standing discipline, every carried-over exponent, formula and line number was
re-checked before use. Results:

| Carried-over claim | Source | Verdict |
|---|---|---|
| `lem:s51-gamma` at line 11446, statement and proof | proof doc | **Correct.** Read directly; the equation `∂_tΓ + u_r∂_rΓ + u_z∂_zΓ = ν(Δ − (2/r)∂_r)Γ` and the max-principle proof are as the handover described. |
| `lem:wpt2-aniso` at line 13159 (S54) | proof doc | **Correct**, read directly; not used in this session's derivations. |
| `cor:wpt2-radius-cap` at line 13451 (S54) | proof doc | **Correct**, read directly. |
| Derivation of `Γ⁽⁰⁾ = q^{-h}H` (S54, `eq:wpt2-gamma`) | proof doc + source (4.3) | **Correct.** Re-derived from scratch: `√(2qX)·q^{-A}·C^{-1}√(2X)φ = (2Xφ/C)q^{1/2−A}` and `1/2 − A = −h`. |
| S54's `sup\|Γ⁽⁰⁾\| ≳ τ^{-h}` | S54 log | **Correct, and sharpened.** It is not merely `≍` — the leading field is exactly self-similar, so `sup Γ⁽⁰⁾ = c_*τ^{-h}` *exactly*, which is what makes the floor exact rather than heuristic. |
| S54's `ℓ_r ≍ τ^{1/2}` as an input to the floor | handover | **Not needed as an input.** The maximiser's radius comes out as `c_r τ^{1/2}` from the profile itself (Step 5 of `prop:wpt3-floor`). Recorded because the handover proposed feeding it in. |
| Handover: "adding `f` adds exactly `r f_θ` … re-derive this yourself" | handover | **Correct**, and re-derived. The three places it could fail (pressure, advection, vector Laplacian) each check out under axisymmetry. |
| Handover's statement of the D1 target | handover | **WRONG.** See §2.3. |
| Handover: "the technical sections almost certainly give a quantitative bound … the residual-absorption piece is the one to check" | handover | **Correct in substance, wrong about location.** The quantitative statement is in **§4** (the stress exponent), not in the pulse/residual-absorption sections §7–§9. |
| Handover's guess that this is S55 | handover | **Correct**, verified twice. |
| Baseline 109 bookmark warnings / 64 overfull / 2 underfull / 196 pp | S54 log | **Correct**, but rebuilt from the tip rather than trusted. |

**Nothing else in S52/S53/S54 was found stale.** The one substantive handover error is the
sign error in D1's statement, recorded below and in `rem:wpt3-handover-sign`.

---

## 2. Deliverable 1 — the lemma

### 2.1 `lem:wpt3-forced-swirl`

For axisymmetric `(u,p)` with `div u = 0` and `f := ∂_t u + (u·∇)u − νΔu + ∇p`:

```
∂_tΓ + u_r ∂_rΓ + u_z ∂_zΓ  =  ν(Δ − (2/r)∂_r)Γ  +  r·f_θ ,      Γ = r u_θ
```

`f_r` and `f_z` never appear, and neither does pressure. Re-derived from the momentum
equations; the three ingredients are `(∇p)_θ = r^{-1}∂_θ p = 0`, the vanishing of the
`∂_θ` term in `((u·∇)u)_θ`, and `(Δu)_θ = Δu_θ − u_θ/r² + 2r^{-2}∂_θ u_r = Δu_θ − u_θ/r²`.
Each is exactly a place a forced version could go wrong, and each is checked in the proof.

### 2.2 `lem:wpt3-forced-max` — the hypotheses it actually needs

The lemma needs **three** hypotheses and no more:

- **(H1) regularity** — `u` is `C¹` in `t`, `C²` in `x`; `p` is `C¹` in `x`. Nothing is
  required of `u_r, u_z` beyond this, because the drift term is *annihilated* at the
  maximiser, not estimated.
- **(H2) attainment with uniform compactness** — for each compact `J` of times there is a
  fixed compact `K_J` containing every maximiser. Satisfied if `Γ(·,t)` is compactly
  supported in a fixed compact, or decays to zero uniformly on `J`. **This is the
  hypothesis that does real work**: without the *uniform* compactness the Danskin upper
  bound in (iii) fails, because there is nowhere for the maximisers `x_k` to converge.
- **(H3) nondegeneracy** — `M(t) := max_x Γ > 0`. Needed only to force the maximiser off
  the axis.

Notably **not** needed: any bound on `f`, any decay rate, any smallness, any hypothesis
about the drift, any differentiability of `M`, and — because of (i) below — any special
treatment of the coordinate singularity at `r = 0`.

Conclusions:

1. **(i) Every maximiser is off the axis.** `Γ = r u_θ = O(r²)` at the axis by smoothness
   of `u`, so `Γ ≡ 0` there, and (H3) puts every maximiser at `r_m > 0`. **Consequence: the
   singular coefficient `−2ν/r` is never evaluated at `r = 0`, and the ε-exhaustion
   `{r ≥ ε}` used in `lem:s51-gamma`'s own proof is not needed here.** This is the one place
   where the forced statement is *cleaner* than the unforced one it generalises.
2. **(ii) `∂_tΓ(x_m,t) ≤ (r f_θ)(x_m,t)`.** At an interior max of the `(r,z)`-function
   `Γ`: `∂_rΓ = ∂_zΓ = 0` and the Hessian is negative semidefinite. So the drift vanishes
   *whatever `u_r,u_z` are*, and
   `ν(∂_r²Γ − r_m^{-1}∂_rΓ + ∂_z²Γ) = ν(∂_r²Γ + ∂_z²Γ) ≤ 0` — **the singular first-order
   term is annihilated by `∂_rΓ = 0`, not bounded.** This is the entire content of why the
   unforced maximum principle holds, and it is untouched by adding forcing.
3. **(iii) Danskin.** `M` is locally Lipschitz, and the one-sided derivatives exist
   everywhere with `d⁺M/dt = max_{Σ(t)} ∂_tΓ` and `d⁻M/dt = min_{Σ(t)} ∂_tΓ`. Both
   directions are written out: the lower bound from `M(t+s) ≥ Γ(x_m,t+s)` for a fixed
   maximiser; the upper bound by extracting `x_k ∈ Σ(t+s_k)`, using (H2) to get
   `x_k → x̄`, joint continuity plus continuity of `M` to get `x̄ ∈ Σ(t)`, and the mean
   value theorem in `t`. `M` is then differentiable a.e. by Rademacher.
4. **(iv) The forced maximum principle**, in Dini form and in an **integrated form that
   assumes no differentiability of `M` at any point**:
   `M(t) − M(s) ≤ ∫_s^t max_{Σ(σ)} (r f_θ) dσ`. At `f ≡ 0` this returns
   `lem:s51-gamma` exactly, so it is a genuine generalisation, not a parallel statement.
5. **(v) The `‖Γ‖_∞` form, with a sign factor** — see below.

### 2.3 The handover's statement of D1 is false, and had to be corrected

The handover asked for `d/dt ‖Γ‖_∞ ≤ (r f_θ)(x_m(t),t)` with
`M_Γ := ‖Γ‖_{L^∞}`. **As written for the absolute-value norm this is false.** If the norm
is attained at a point where `Γ < 0`, that point is a global *minimum* of `Γ`; there the
Hessian is positive semidefinite, the diffusion term is `≥ 0`, and the correct inequality
is `d/dt M_Γ ≤ −(r f_θ)(x_m,t)` — the opposite sign. The correct statements are:

- for the **signed** maximum `M = max_x Γ`: the handover's inequality, verbatim;
- for the **norm**: `d⁺M_Γ/dt ≤ max_{Σ₊∪Σ₋} sgn(Γ)·(r f_θ) ≤ max_{Σ₊∪Σ₋} |r f_θ|`.

Nothing downstream changes — D2 uses the signed maximum and reports a floor on `|r f_θ|`,
so the sign factor is absorbed — but the lemma had to be stated correctly to be usable,
and the erroneous form should not be propagated. Recorded in `rem:wpt3-handover-sign`.

### 2.4 Pointwise vs uniform, stated explicitly (`rem:wpt3-pointwise`)

The handover flagged this as the most likely place for a marginal error, so it is separated
out rather than left implicit. What is proved is a statement **at one moving point**, a set
of measure zero in space-time. Therefore:

- **(a)** No positive-measure statement follows. A forcing consisting of a single tall
  spike tracking the maximiser, supported on a set of measure ε, satisfies the conclusion.
- **(b)** No `L^q` bound for `q < ∞` follows: that spike has `‖f‖_{L^q} = O(ε^{1/q}) → 0`.
- **(c)** The `L^∞` statement *does* follow, trivially, and is the **weaker** one — it is
  what you get by throwing the location information away. The direction matters and is
  easy to get backwards.

The uniform statement one might want — "`|f_θ|` is bounded below throughout the core" — is
a genuinely different and stronger claim and is **not proved**. Everything in D2 and D3 is
stated pointwise at `x_m(t)` and should be quoted that way.

---

## 3. Deliverable 2 — the floor, derived in full

**Step 1.** `Γ⁽⁰⁾ = r u_θ⁽⁰⁾ = √(2qX)·q^{-A}·C^{-1}√(2X)φ = (2Xφ/C)q^{1/2−A} = q^{-h}H(X,η)`
with `H = 2Xφ/C`, using `1/2 − A = −h`. Re-derived independently of S54; it is also the
source's own identity.

**Step 2.** `τ = q(1−η²)` gives `q = τ/(1−η²)`, so `Γ⁽⁰⁾ = τ^{-h}(1−η²)^h H(X,η)`. At fixed
`t` the pair `(X,η)` sweeps the whole profile domain as `(r,z)` sweeps the half-plane
(`η, q` are fixed by `(z,τ)`, then `X = r²/2q` sweeps `(0,∞)` with `r`). Hence

```
M(t) := sup_x Γ⁽⁰⁾ = c_* τ^{-h}   EXACTLY,     c_* := sup_{(X,η)} (1−η²)^h H(X,η)
```

with `c_*` a constant independent of `t`. **This exactness — not `≍` — is what makes the
floor exact rather than a scaling heuristic**, and it is available only because the leading
field is exactly self-similar.

**Step 3 (attainment).** `H(0,η) = 0`; `H → 0` as `X → ∞` because the source's exterior
power is `H_pow = √(2X)c_∞X^{-A} = √2 c_∞ X^{-h}`; and `(1−η²)^h → 0` at `η = ±1` where `H`
stays bounded (smoothness on the closed rectangle). So the sup is attained at some
`(X_*, η_*)` with `0 < X_* < ∞` and `|η_*| < 1`, and `c_* ∈ (0,∞)`. (H1)–(H3) hold.

**Step 4.** `τ = 1 − t`, `dτ/dt = −1`, so
`dM/dt = c_*(−h)τ^{-h-1}·(−1) = h·c_*·τ^{-(1+h)} > 0`.

**Step 5 (the radius is an output).** `r(x_m) = √(2q_*X_*)` with `q_* = τ/(1−η_*²)`, so
`r(x_m) = c_r τ^{1/2}` with `c_r := (2X_*/(1−η_*²))^{1/2} > 0`. **`ℓ_r ≍ τ^{1/2}` is not
assumed** — the maximiser tracks the core radius as a consequence of Steps 1–3.

**Step 6.** Lemma D1 gives `h·c_*·τ^{-(1+h)} = dM/dt ≤ max_{Σ(t)} (r f_θ⁽⁰⁾)`, so

```
FLOOR A:   |r f_θ⁽⁰⁾| ≥ h·c_*·τ^{-(1+h)}       = h·c_*·τ^{-(A+1/2)}
FLOOR B:   |f_θ⁽⁰⁾|   ≥ (h·c_*/c_r)·τ^{-(3/2+h)} = (h·c_*/c_r)·τ^{-(A+1)}
```

pointwise at the running swirl maximum, at radius `c_r τ^{1/2}`. The second forms use
`A = 1/2 + h`, so `A + 1/2 = 1 + h` and `A + 1 = 3/2 + h`.

### 3.1 Calibration — why `τ^{-(A+1)}` is the meaningful number

Every individual term of the leading field's own θ-momentum equation is the same size:

```
|∂_t u_θ⁽⁰⁾|          ≍ τ^{-A}/τ                    = τ^{-(A+1)}
|u_r⁽⁰⁾ ∂_r u_θ⁽⁰⁾|   ≍ τ^{-1/2}·τ^{-A}/τ^{1/2}     = τ^{-(A+1)}
|∂_r² u_θ⁽⁰⁾|         ≍ τ^{-A}/τ                    = τ^{-(A+1)}
|∂_z² u_θ⁽⁰⁾|         ≍ τ^{-A}·τ^{-2D}              = τ^{-(A+1)+2h}   ← smaller by τ^{2h}
```

So the floor says the azimuthal residual is **of full leading order in its own momentum
equation** — not a perturbation, not lower order, and not made so by any choice of the free
profile data. That is the quantitative content of Wall 9(b)'s "load-bearing, not
decorative". As a by-product the same computation reproduces, from the floor's side, the
source's own structural claim that radial viscosity is retained in the leading balance and
axial viscosity is not.

### 3.2 A derived consequence the source does not state (`cor:wpt3-annulus`)

From `lem:wpt3-forced-swirl` at `ν = 1`, and the source's definition
`R_θ⁽⁰⁾ := R(u⁽⁰⁾,p⁽⁰⁾)·e_θ + ∂_z²u_θ⁽⁰⁾` (its (4.12)):

```
r·R_θ⁽⁰⁾ = D_t Γ⁽⁰⁾ − (∂_r² − r^{-1}∂_r) Γ⁽⁰⁾
```

so `R_θ⁽⁰⁾ = 0` on an open set means the material derivative there is **radial diffusion
only**. If the swirl maximiser sat in such a set, then at it `∇Γ = 0`, the drift vanishes,
the first-order radial term is annihilated, and `dM/dt = ∂_r²Γ ≤ 0` — contradicting
`dM/dt = h c_* τ^{-(1+h)} > 0`. Since the source makes `T_0 = 0` for `X ≤ X_a` and
`X ≥ X_b` (Thm 4.6(ii)) and `R_θ⁽⁰⁾` a cylindrical divergence of `q^{-A-1/2}T_0` (Prop 4.2):

> **`X_* ∈ [X_a, X_b]` — the swirl maximum of the leading profile is forced into the active
> annulus where the pulses act. This is a consequence of the swirl maximum principle, not a
> design choice.**

### 3.3 An independent cross-check that validates all of §2 (`rem:wpt3-414-check`)

The forced swirl equation was tested against the source's (4.14), a formula stated there
for an unrelated reason. On `X ≤ X_a` the profiles satisfy (4.13), and (4.14)+(4.13) give

```
D_t(q^{-h}H) = −(q^{-h-1}/L)·H·S_q = 2q^{-h-1}(H/φ)(Xφ_XX + 2φ_X) = (4X/C) q^{-h-1}(Xφ_XX + 2φ_X)
```

Independently, `r R_θ⁽⁰⁾ = 0` there, so the same material derivative equals
`(∂_r² − r^{-1}∂_r)Γ⁽⁰⁾`; computing that from `Γ⁽⁰⁾ = q^{-h}H(X)`, `X = r²/2q`,
`∂_r = (r/q)∂_X` at fixed `(z,t)`:

```
(∂_r² − r^{-1}∂_r)Γ⁽⁰⁾ = 2X q^{-h-1} H_XX = (4X/C) q^{-h-1}(Xφ_XX + 2φ_X)
```

**Identical.** This simultaneously checks the source's (4.14), this session's forced swirl
equation, and the identity `cor:wpt3-annulus` rests on. It is the strongest single piece of
verification in this session.

---

## 4. Deliverable 3 — the comparison

**Primary source re-fetched this session** from
`https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`
(165 pp + appendices), and read for the quantitative forcing statements. Theorem 1.1 is as
S53/S54 recorded.

### 4.1 The handover was right that a bound exists, and wrong about where

The handover expected the quantitative statement in the pulse-seeding / residual-absorption
sections (§7–§9). **It is in §4**, as a stress exponent, and it is stated cleanly.

### 4.2 What the source states, verbatim

**Definition of the force as the residual** (§2):

> "For any incompressible flow u and pressure p, we can always define the external force f
> to be the residual in (1.1)."

**Qualitatively — no rate attached** (§2, and again §2.2):

> "For the background alone, the momentum residual becomes unbounded as t ↑ 1, so it cannot
> supply the smooth force required by the theorem."

> "…needed to sustain the background becomes unbounded as t ↑ 1, so it cannot serve as the
> smooth external force required by the construction."

**Quantitatively** — after its (4.11), and in Proposition 4.2:

> "In the notation of Section 3, the physical stress representation to be proved is
> T = q^{−A−1/2} T₀."

> "The residuals R_θ^{(0)}, R_z^{(0)} in (4.12) are minus the cylindrical divergences
> ∂_r + 2/r and ∂_r + 1/r of q^{−A−1/2} T₀."

with nondegeneracy of `T₀` on the annulus, Theorem 4.6(iv):

> "|T₀| ≥ cζ,  |∂^α T₀| ≤ C_α ζ δ^{−m_α}"

(`ζ > 0` exactly on `X_a < X < X_b`), and, in the proof of Proposition 4.2, the azimuthal
case written out in the same variable this session's lemma uses:

> "The angular equation is simplest in angular momentum r u_θ^{(0)} = q^{−h} H."

with material derivative `−q^{−h−1} H S_q / L`, and

> "The axial material derivative plus the pressure gradient is similarly −q^{−A−1} S_n/L."

The same exponents are carried to the **full axisymmetric background** by Proposition 5.5:

> "R(u_B, p_B) = −(∂_r + 2/r) T_phys,θ e_θ − (∂_r + 1/r) T_phys,z e_z + E_B"   (5.41)

with `|∂^α E_B| ≤ C_{α,M} q^M` for every `M`, and (5.43) fixing
`T̂ = q^{A+1/2} T_phys` with `|D^I(T̂ − T₀)| ≤ C_I q^{2h} ζ δ^{−N_I}`.

**The final force, by contrast, is flat** — Theorem 3.1(iii):

> "|∂_x^α ∂_t^b R(u,p)| ≤ C_{α,b,N,X₁} q^N   (0 ≤ X ≤ X₁, q ↓ 0)"

> "For X ≥ X_ext the residual is identically zero."

### 4.3 Converting and comparing

The floor applies to any axisymmetric divergence-free field whose swirl maximum grows like
`τ^{-h}`. Two objects in the source qualify: `u⁽⁰⁾`, and the **axisymmetric** background
`u_B` of Prop 5.5 (its (5.42) gives `q^A u_θ,B = E₀ + O(q^{2h})`, so
`M_B(t) = c_*τ^{-h}(1+o(1))`). Neither is the object of Theorem 1.1.

For `u_B`: `f_θ,B = −(∂_r + 2/r)T_phys,θ + (E_B)_θ` with
`T_phys,θ = q^{-A-1/2}(T_{0,θ} + O(q^{2h}))` and `E_B` flat. Since `q` depends on `(z,t)`
only, `∂_r = (r/q)∂_X` at fixed `(z,t)` and `r² = 2qX`, so

```
r f_θ,B = −2 q^{-A-1/2} ∂_X(X T_{0,θ}) + O(q^{-A-1/2+2h}) + flat
        ≍ q^{-A-1/2} = q^{-(1+h)} ≍ τ^{-(1+h)}          on the active annulus
|f_θ,B| ≍ τ^{-(A+1)} = τ^{-(3/2+h)}
```

### 4.4 Verdict on the comparison

**Clean confirmation, exponent-exact, with no slack.** Floor A demanded
`τ^{-(A+1/2)}`; the source states `q^{-A-1/2}`. Floor B demanded `τ^{-(A+1)}`; the source
states `−q^{-A-1}S_n/L` for the axial material derivative and the same power follows for
the azimuthal one. **Two genuinely different routes** — the swirl maximum principle plus
the single exponent `Γ⁽⁰⁾ ≍ τ^{-h}` on one side, radial integration of the residual on the
other — **give the same two exponents.** It is not circular: neither derivation uses the
other's input.

**Constants are not compared.** The floor's constant is `h·c_*` and `h < 1/100`, so the
floor is weak at constant level; the source states no numerical value for `T₀` at any
particular `X`. **The comparison is exponent-level only** and must be quoted that way.

### 4.5 The one derived requirement the source does not state — a fact to check, not a discrepancy

`cor:wpt3-annulus` forces `X_* ∈ [X_a, X_b]`. Partial independent evidence from the source,
in both directions: `H(0,η) = 0` and `H ∼ √2 c_∞ X^{-h}` at large `X`, so `X_*` is
interior; but its (4.30) gives `E = c_patch(1+η²)^{-1}X^{-1/2-λ}` with **λ > 0** (its
Theorem 4.6), on the reserved intervals `I_pos, I_mean ⊂ (X_a, X_v)`, so
`H = √(2X)E ∝ X^{-λ}` is already *decreasing* there, placing `X_* ≤ inf I_pos`. The two
together pin `X_*` to the narrow window `[X_a, inf I_pos]` near the inner edge of the
annulus.

**Nothing read this session confirms or contradicts `X_* ≥ X_a`.** Settling it needs the
inner-profile construction of the source's Appendix B, which was not re-derived here and is
outside one session. It is recorded as a requirement the leading profile must satisfy for
its own stated swirl growth to be possible — checkable in one line by anyone holding the
profile — and **explicitly not as a defect**: the construction's whole design puts the
momentum transport in that annulus, so the requirement is the expected one. Per the
handover's instruction, it is reported as a fact to check, not as a refutation.

### 4.6 No tension with Theorem 1.1, and none is implied

The final force is flat (§4.2 above), and the floor says nothing against that, for exactly
the reason Wall 9 already records: **the full field is not axisymmetric** — the oscillatory
pulses' "cylindrical velocity components have zero angular average" but their momentum
fluxes do not — so `lem:wpt3-forced-max` does not apply to it at all. What the floor bounds
is the residual of the *axisymmetric* background, an intermediate object whose residual the
construction never claims is small; it says the opposite and supplies it internally through
the pulses' averaged quadratic products rather than through `f`.

Read the other way round, this is the useful direction: **any mechanism carrying the
axisymmetric background's swirl from `O(1)` to `τ^{-h}` must inject azimuthal momentum at
full leading order.** That is why breaking axisymmetry is not an optional convenience of
that construction, and it is the sharpest form of Wall 9 available.

---

## 5. Verdict, in the narrow terms the handover required

**This is documentation-sharpening, not a step toward (A)/(B) and not a challenge to
(C)/(D).** Specifically:

1. **Documentation-sharpening.** Wall 9's "load-bearing, not decorative" now carries a
   number: `τ^{-(A+1)}` pointwise at the swirl maximum, full leading order in the
   θ-momentum equation, and the source's own stated exponent agrees. Wall 9's
   classification is **unchanged**: `proved impossible` as an import route.
2. **Not a step toward (A) or (B).** Nothing here bears on global regularity for the
   unforced equations. It proves a *lower* bound on a forcing — the opposite direction of
   travel.
3. **Not a challenge to (C) or (D).** No error in `\cite{OpenAI2026}` is asserted and none
   was found. The one derived requirement the source does not state is reported as a fact
   to check (§4.5). That construction is not independently re-verifiable from outside in
   one session and this section does not attempt it.
4. **Wall 1, Wall 3, Region (R2) and `prop:wpt2-swirl-exclusion` are untouched.** The class
   is axisymmetric and the statement is pointwise.
5. **What is genuinely new is small and worth naming exactly:** the forced generalisation
   of `lem:s51-gamma` (elementary, but it did not exist in this document, and the
   handover's form of it was wrong); the exponent `τ^{-(A+1)}`; and `cor:wpt3-annulus`,
   which is the only statement here that says something the source does not already say
   about itself.

**No WP-T4 is proposed.** The natural next questions — whether `X_* ≥ X_a` holds for the
source's profile, and whether a uniform (positive-measure) version of the floor is
provable — are recorded above and deliberately not manufactured into a work package.

---

## 6. Files changed

```
proofs/claim_a_3d_proof_attempt.tex   +737   (14,685 -> 15,422 lines)
      new sec:wpt3-forcing-floor (§41); a quantitative addendum item (e) to Wall 9;
      row 9' in the status table; an extended OpenAI2026 reading record.
      No new bibliography entries. Otherwise purely additive.
PROGRESS.md                           KEY FINDING S55 block, session-log row, Active milestone
SESSION-LOG/2026-09-12-S55-wpt3-forcing-floor.md   (this file)
```

Not committed, not pushed. Worktree `/home/user/wpt3-worktree`, branch
`claude/wpt3-forcing-floor`, for independent verification before merge. **Re-check the
session tag and the tip at merge time.**
