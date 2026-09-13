# S57 — WP-A: pricing Grujić's log-weighted bmo hypothesis against σ*-decay

**Date:** 2026-09-13
**Branch:** `claude/shared-conversation-jgoavo`.
**Tip at start and at finish:** `19f4d14` ([S56] plan: WP-A handover). Highest session
log at start = S56. **S57 is correct**, verified twice (start and before writing this
log). Working tree clean at start; no other agent branch present.
**Proof-document work done in an isolated detached worktree at `19f4d14`:**
`/tmp/claude-0/-home-user-tfirstframework/5dc09b30-469f-5a15-85ac-3ac7e0d30549/scratchpad/wpa-s57`
The live `proofs/claim_a_3d_proof_attempt.tex` is **untouched** (15,422 lines, unchanged).
**Nothing committed, nothing pushed.** `book/`, `layer3/`, `layer4/`, `results/`
untouched. `PROGRESS.md` deliberately not edited (S52/S56 precedent); ready-to-paste
block in §7.

---

## 0. One line per deliverable

| Deliverable | Answer |
|---|---|
| **D1** — read both papers, classify | Done from the primary HTML sources in full. Part I's chain is **proved from two assumed hypotheses**; Part II **does not propagate the hypothesis** — it proves a conditional transfer theorem under four further hypotheses and then exhibits, in its own §7.5, two natural critical cores on which they fail. |
| **D2** — the dictionary | **Refuted as stated.** The S56 Poincaré translation is not an equivalence and not a two-way implication. One one-way implication survives, for an unrestricted scale family. Three explicit divergence-free counterexamples. |
| **D3** — commutator/CRW | **Holds.** Verified link by link, including the CRW application, which reduces to CRW for Riesz transforms with no general CZ version needed. One cosmetic error in the source ("reflexive Lorentz spaces"). |
| **D4** — verdict | Handover's **line 2**: the translation cannot be completed, for three specific stated reasons. The hypothesis is **logically independent** of σ*-decay — not the fifth avatar S56 expected — but it **is** a fourth member of the (AX)/saturation no-concentration family, and it is not cheaper. |

---

## 1. Deliverable 1 — what is proved, what is assumed, what is asserted

Both papers fetched and read in full (arXiv HTML v3 / v2), not from abstracts. Recorded
verbatim-faithfully in the new `Fact` of §`sec:wpa-s57`.

**Part I (`arXiv:2607.08866v3`), the chain:**

1. **Assumed** — Definition 2.1, the "critical point singularity": `ω = Φ(x,t)|x|^{-2}`,
   `Φ` bounded, `|∇Φ| ≲ |x|^{-1}`, `Φ` scale-invariant or log-periodic, with the
   consequence that `A_λ(t) = {ω>λ}` is **contained in** `B_R`, `R ≤ Cλ^{-1/2}`.
2. **Assumed** — `ξ ∈ L^∞_t bmo_{1/|log r|}(ℝ³)`.
3. **Proved** (and not new) — the unidirectional cancellation. It is
   Constantin–Fefferman's `det(ŷ, ξ(x+y), ξ(x))` written as a commutator:
   `α = ξ·[𝒯,ξ](ω)·ξ`.
4. **Proved from 1+2** — Theorem 4.1, `‖α(·,t)‖_{L^{3/2,∞}(B_R)} ≤ C₀/|log R|`.
5. **Proved from 4** — §5–§6: truncated-enstrophy/Lorentz-interpolation gives
   `|{ω>λ}| ≤ Cλ^{-3/2}(log λ)^{-3/2}`, then O'Neil gives
   `|{|u|>λ}| ≤ Cλ^{-3}(log λ)^{-3}`.
6. **Proved from 5** — Theorem 7.4: the sparseness scale
   `r_s ∼ ‖u‖_∞^{-1}(log‖u‖_∞)^{-1}` drops below the analyticity radius
   `ρ_s ∼ ν‖u‖_∞^{-1}`; escape-time + harmonic measure gives the contradiction.

Every internal step I checked in §4–§7 is arithmetically correct (the `L^{3,1}`
rearrangement `‖|y|^{-3}χ_{|y|>2R}‖ ≍ R^{-2}`; the `L²→L^{3/2,∞}` and `L⁴→L^{3,1}`
local embeddings; the telescoping of ball means; the reversed double sum; the Young
absorption; the coercive damping ODE; the O'Neil integrations; the Solynin/Ransford
parameter locking). **I take no position on whether Theorem 7.4 is true.** A referee
report was not attempted and is not a session's work.

**Part II (`arXiv:2609.05720v2`):** its §2 direction PDE (3) is `eq:dir-eq` **term for
term** — S56's claim verified against the primary source, including
`P_{ξ⊥}((ξ·∇)u) = P_{ξ⊥}(Sξ)`. But:

- **Theorem 5.7 / Corollary 5.8 are conditional on four hypotheses**, all assumed:
  (H*k*), an inflow smallness condition `sup_A |x|(u_r)_- + 2νC_Φ ≤ (3−k)ν`;
  a tangential strain **three logs** below critical,
  `|P_{ξ⊥}Sξ| ≤ Λ|x|^{-2}|log|x||^{-3}`; a cap condition; and a core-history
  modulus `ω_core(σ;t) ≤ Λ₀/|log σ|`.
- **§7.5 is the paper's own negative result and it is decisive.** Proposition 7.10: a
  scale-invariant `e_*`-axisymmetric tube core has `|F_tan| = Φ₀|x|^{-2}|sin 2ϑ|` —
  exactly critical, so the strain hypothesis fails at order one and "the direction winds
  logarithmically as ℓ(t)→0". Proposition 7.12 + Remark 7.13: a near-Beltrami core has
  critical twist and `‖ξ‖_{bmo_{1/|log r|}} ≳ |log r|`.
  **The propagation question is open, and the companion says so itself.**

**Two defects found in the sources and recorded (not used):**

1. **The two papers do not meet.** Theorem 7.4 and its restatement (Part II, Thm 7.1)
   both demand `ξ ∈ bmo_{1/|log r|}(ℝ³)` — the sup in the norm is over **all** `x ∈ ℝ³`
   and `0<r<1/2`. Part II's Corollary 5.8 and Theorem 7.9 deliver only
   `‖ξ‖_{bmo_{1/|log r|}(B_{ρ₁})} ≤ C`. The gap is bridgeable (Theorem 4.1's proof uses
   the hypothesis only on balls inside `B_{2R^{1/2}}`) but the local Theorem 7.4 that
   Theorem 7.9 invokes is **stated nowhere**.
2. **The global form is close to vacuous.** By the counterexample below, `ξ` fails the
   bmo condition at *every* scale near any transversal zero of `ω`. So the global
   hypothesis also silently asserts `ω` has no such zero anywhere in `ℝ³`.

---

## 2. Deliverable 2 — the dictionary, settled

### 2.1 What is true (proved)

Exact scale identification first: under Definition 2.1 with `c₀ ≤ Φ ≤ C₀`, the source's
inner scale `ℓ(t)` and this document's `r* = M^{-1/2}` agree up to fixed constants, so
`ω ∼ |x|^{-2}` **is** `M·(r*)² = 1`. And Part II's critical direction bound
`|∇ξ| ≤ C₁/max(|x|,ℓ)` is exactly "**σ\* bounded**" — the borderline of
`rem:marginality` and the Regime-B stalemate of `prop:regime-B`. It is not σ*-decay, and
the source explicitly says it carries no logarithm.

Then, with `σ̃_ρ(y) := ρ^{-1}∫_{B_ρ(y)}|∇ê|²` (**no indicator**), Poincaré–Wirtinger gives

```
⨍_{B_ρ(y)} |ê − ê_{B_ρ}| ≤ C_P (3/4π)^{1/2} σ̃_ρ^{1/2} ≤ ¼ σ̃_ρ^{1/2},
```

`C_P = μ₁(B₁)^{-1/2} ≤ ½`. Hence the **only** surviving implication:

> **(Σ_log)** `σ̃_ρ(y) ≤ K|log ρ|^{-2}` for every `B_ρ(y) ⊆ B_{r*}(x*)`
> **⟹** the local `bmo_{1/|log r|}` seminorm of `ê` on `B_{r*}(x*)` is `≤ ¼K^{1/2}`,

and it is **strict**. At the top scale `ρ = r*` this reads `σ̃_{r*} ≤ 4K/(log M)²` —
which is the handover's sketch, confirmed **in form only**, as a *sufficient* condition,
for the *unrestricted* deficit, at *one* ball.

### 2.2 What is false (three counterexamples, all smooth and divergence-free)

**(C1) bmo ⇏ σ\*.** `ω^(δ) = M(cos g_δ(x₃), sin g_δ(x₃), 0)`,
`g_δ(s) = sin(s/δ)/(2|log δ|)`. Divergence-free, `|ω| ≡ M` so the indicator is **inert**.
Then `‖ê‖_{bmo_φ} ≤ 2` uniformly in δ, while
`σ* ≥ (π/9)(r*)²/(δ² log²(1/δ)) → ∞`. **The hypothesis carries no `L²` information
about `∇ê` at all.**

**(C2) The indicator is not removable.** `ω^(δ) = M cos(x₃/δ) e₁`. On
`{|ω| ≥ M/4}` the direction is locally constant, so **`σ* = 0` exactly**; but on any
ball centred on a zero plane `⨍|ê − ê_B| = 1` at *every* radius, so
`‖ê‖_{bmo_φ} = +∞` and `σ̃_ρ = +∞`. This is generic, not pathological: it is what
happens at any transversal zero of `ω`. **The `1_{|ω|≥M/4}` in `def:sigma-star` is not a
convenience — it is what keeps σ\* finite**, and any unrestricted-oscillation hypothesis
imported here needs a zero-free core supplied as a separate assumption.

**(C3) One scale does not give the family.** A localized defect
`ω^(δ) = e₃ + ∇×(ψ_δ e₁)` with `ψ_δ(x) = δχ₀(|x−y₀|/δ)` has `σ* ≤ Cδ/r* → 0` while its
mean oscillation on `B_δ(y₀)` is a **fixed** `c₀>0` (the field on `B_δ` is an exact
dilate of a fixed one). bmo_φ demands `c₀ ≤ K/|log δ|`: fails by a factor `|log δ|`.

### 2.3 The two named subtleties, handled

- **The `{|ω| ≥ M/4}` restriction.** It cannot be dropped, by (C2), and dropping it is
  not a technicality: the unrestricted deficit is *generically infinite*. Even where `ω`
  is zero-free, the honest restricted version costs a relative-Poincaré constant for the
  set `E = B_{r*} ∩ {|ω|≥M/4}` (not scale-invariantly controlled; `E` need not be
  connected or uniform, so the Jones extension that Theorem 4.1 needs is unavailable on
  `E`) and a volume-fraction penalty `(|B_{r*}|/|E|)^{1/2}`. Grujić's hypothesis is
  therefore **strictly stronger in support**.
- **Every scale, not just `r*`.** Verified directly from Theorem 4.1's proof: the
  hypothesis is consumed on every sub-ball of `B_{2R}` (near field, through
  `‖ξ‖_{BMO(B_{2R})} ≤ Cφ(2R)` and the Jones extension), on `B_R` (John–Nirenberg), and
  on the dyadic balls `B_{2^j R}` up to radius `R^{1/2}` (oscillation tail) — for every
  small `R`, with `R` slaved to the level by `R = Cλ^{-1/2}`. So the condition is needed
  at all scales in `(0, R^{1/2}]`, which straddles `r*` in both directions. It does
  **not** collapse to the single critical scale, and by (C3) the collapse is false.

### 2.4 Implications, stated in both directions

- **vs σ\*-decay:** neither implies the other. `(Σ_log) ⟹ bmo_φ` strictly; nothing runs
  back. σ*-decay is `L²`-gradient / one ball / one scale / level-restricted; the
  hypothesis is `L¹`-oscillation / all balls / all scales / unrestricted. Weaker in
  norm, stronger in scale, stronger in support. **Incomparable.**
- **vs (AX) (`prop:s51-ax`):** (AX) is a pointwise bound on a vorticity *magnitude*
  component inside the axisymmetric class; the hypothesis is an integrated bound on the
  *direction* in no symmetry class. No implication either way; no common refinement.
- **vs saturation (`thm:s50-equivalence`):** spatial vs temporal, no implication. But the
  *currency* comparison is exact and is the structurally interesting statement:
  `thm:s50-equivalence` says removing one logarithm from `‖∇u‖_∞ ≲ M𝓛` is equivalent to
  log-free saturation; Theorem 4.1 buys **exactly one logarithm** — but in
  `‖α‖_{L^{3/2,∞}(B_R)}`, not in `‖∇u‖_∞`. That is a Wall 1(d)-permitted reformulation
  paying the same single-logarithm deficit in a weaker norm, which is precisely why
  `prop:s51-poloidal-sharp` does not contradict it: that counterexample saturates an
  `L^∞` bound and says nothing about a restricted Lorentz norm of `ê·Sê`.
- **Family verdict.** `rem:s51-verdict-log`'s pattern gets a fourth instance:
  scale-invariant, imposed not proved, and — by the companion's own §7.5 — not
  established for NS solutions. **It has not been shown cheaper than what it replaces.**
  What is new is that this one is *not* a relabelling of σ*-decay.

---

## 3. Deliverable 3 — the commutator/CRW step: **it holds**

Checked in isolation, link by link:

1. **The cancellation is Constantin–Fefferman's**, correctly restated. The fully
   symmetric contraction of the strain kernel against `e⊗e⊗e` vanishes — the
   `det(ŷ, ξ(x+y), ξ(x))` factor — so `α = ξ·[𝒯,ξ](ω)·ξ` and `|α| ≤ |[𝒯,ξ](ω)|`.
2. **Smallness of the symbol is correct.** `φ(r)=1/|log r|` is *increasing* on `(0,1)`,
   so sub-balls of `B_{2R}` of radius `r ≤ 2R` have oscillation `≤ ‖ξ‖_{bmo_φ}φ(2R)`,
   giving `‖ξ‖_{BMO(B_{2R})} ≤ Cφ(2R)`. The monotonicity is used in the right direction
   here and again at the dyadic truncation `2^{N+1}R ≈ R^{1/2}`.
3. **Jones extension is legitimate.** A ball is a uniform domain with absolute
   scale-invariant constants. The substitution
   `[𝒯,ξ](ωχ_{B_{2R}})(x) = [𝒯,ξ̃](ωχ_{B_{2R}})(x)` for `x ∈ B_R` is valid because the
   source is supported where `ξ̃ = ξ` and `x ∈ B_R ⊂ B_{2R}`.
4. **CRW applies, and more cleanly than the source claims.** The operator is a
   composition of two Riesz transforms, and `[R_iR_j,b] = R_i[R_j,b] + [R_i,b]R_j`, so
   the original CRW bound `‖[b,R_j]‖_{L^p→L^p} ≤ C_p‖b‖_{BMO}` suffices — **no general
   Calderón–Zygmund version (Janson/Uchiyama) is needed**. Linear in `‖b‖_{BMO}`.
   Vector-valued bookkeeping is componentwise.
5. **The passage to `L^{3/2,∞}` is correct.** Real (Marcinkiewicz–Hunt) interpolation
   between `L^{p₀}` and `L^{p₁}`, `1<p₀<3/2<p₁<∞`, gives
   `(L^{p₀},L^{p₁})_{θ,∞} = L^{3/2,∞}` with norm `≲ ‖b‖_{BMO}^{1−θ}‖b‖_{BMO}^{θ}`, i.e.
   still linear. **One cosmetic error in the source:** it calls these "the reflexive
   Lorentz spaces". `L^{3/2,∞}` is not reflexive (and not normed in the displayed
   quasi-norm, though normable via `f**` for `p>1`). Conclusion unaffected.
6. **The far-field bookkeeping is correct** — rearrangement, embeddings, telescoping,
   reversed sum — all re-derived independently.

**Consequence for this programme.** The technique stands on its own: *the stretching
eigenvalue is a commutator whose symbol is the direction field, so any smallness of the
direction's BMO norm at small scales is a linear gain on α in a restricted critical
Lorentz norm.* It is the first estimate on record here that extracts a gain from the
direction without an `L^∞` Calderón–Zygmund bound, i.e. inside Wall 1(d)'s permitted
class. **It is not by itself a route around Wall 1**: the smallness it consumes is not
available from σ\* (Deliverable 2) and is not available from the companion paper
(Deliverable 1).

---

## 4. Deliverable 4 — verdict

**Handover's line 2: the translation cannot be completed, for three specific reasons**
((C1), (C2), (C3)) — with a positive residue, the strict one-way implication
`(Σ_log) ⟹ bmo_φ`, and with the F2-flavoured family conclusion intact but sharpened:
the hypothesis is a **fourth member of the (AX)/saturation no-concentration family**,
**not** a relabelling of σ*-decay, and **not proved cheaper**.

Not claimed: that Theorem 7.4 is true or false; that the hypothesis can be propagated;
that the counterexamples say anything about the two conditions *along a blow-up
trajectory* (they refute implications between the conditions as stated, which is the
level S50 set and the level at which cheapness claims are made).

**No wall is reclassified. Wall 1, Wall 3, Wall 9, (AX), saturation all unchanged.**

---

## 5. Stale / imprecise / wrong, checked against primary sources

| Claim | Source of claim | Verdict |
|---|---|---|
| Grujić's direction PDE is `eq:dir-eq` term for term, incl. `P_{ê⊥}((ê·∇)u)=P_{ê⊥}(Sê)` | S56 | **Correct** — verified against Part II §2, eq. (3). |
| `ω ∼ \|x\|^{-2}` / `L^{3/2,∞}` is our `M(r*)²=1` | S56 | **Correct**, and now proved as an exact scale identification (`ℓ ≍ r*`). |
| `lem:direction-eq` at line 5349, `def:sigma-star` at 5393 | handover / S56 | **Correct**, line numbers had not drifted. |
| Theorem 7.4's hypothesis is `ω ∈ L^∞_t L^{3/2,∞}` | S56's rendering | **Imprecise, and it matters.** The hypothesis is Definition 2.1; `L^{3/2,∞}` is only "in particular". Definition 2.1 additionally gives the **containment** `A_λ ⊂ B_{Cλ^{-1/2}}`, which `L^{3/2,∞}` does **not** give (that bounds the *measure*, not the *location*) and which §5 uses essentially. |
| The Albritton–Bradshaw prior is "load-bearing and negative" for `2607.08866`'s endgame | S56 | **Wrong on the point that matters.** Read in full: their negative half concerns the *a priori* sparseness classes `Z^{(k)}_{ᾱ_k}` ("not stronger than `u ∈ L^∞_tL²_x`"; the identification `Z^{(k)}_α ∼ W^{k−1/α,∞}` "should not be taken seriously"). Their *first* half is **positive** — they give "a simple proof that sufficiently sparse Navier–Stokes solutions do not develop singularities", explicitly "an alternative to the approach of [Gru13], which is based on analyticity and the 'harmonic measure maximum principle'" (Thm 1.2, Cor 1.3). `2607.08866` §7 uses only that *conditional* direction, i.e. the half AB reprove and endorse. The prior remains a fair caution about *interpretation*; it does not attach to the step S56 attached it to. |
| "(AX) — S53's structural hypothesis on the axis" | WP-A brief | **Wrong session.** `prop:s51-ax` lives in `sec:axisym-s51`, headed "(Session S51)". S56 attributed it to S52. The document's own labels/headings are `s51`, but `SESSION-LOG/` shows the axisymmetric pivot was session **S52** (`2026-09-11-S52-axisymmetric-swirl-pivot.md`) — a pre-existing off-by-one between the document's section headings and the session logs, of the same kind `rem:s47-handover-corrigenda` already warns about for §20–§22. **Not repaired here**; recorded so it is not rediscovered. |
| "the S39 σ\*-decay gap" | WP-A brief | **Imprecise.** §S39 (`sec:rigidity-s39`) converts σ*-decay into a Liouville problem; the gap itself is unified at S37 (`cor:unification`) and proved *necessary* at S38 (`prop:regime-A-impossible`, and the remark after it). |
| Handover's guess that the session tag is S57 | handover | **Correct this time** — tip `19f4d14`, highest log S56, verified twice. |

---

## 6. Files changed

```
SESSION-LOG/2026-09-13-S57-wpa-grujic-bmo-pricing.md         (this file, new, uncommitted)
```

**In the isolated worktree only** (`…/scratchpad/wpa-s57`, detached at `19f4d14`):

```
proofs/claim_a_3d_proof_attempt.tex   15,422 → 16,193 lines
  + \section{...}\label{sec:wpa-s57}  inserted before \section{Routes and Walls...}
  + 6 bibliography entries: Grujic2026Log, Grujic2026Decay, AlbrittonBradshaw2022,
    CRW1976, Jones1980, Hunt1966
```

New labels (no collisions; the string `s57` did not occur in the document before):
`sec:wpa-s57`, `subsec:s57-scope`, `rem:s57-standing`, `subsec:s57-sources`,
`fact:s57-source` (+ items `it:s57-def21/cancel/thm41/degiorgi/thm74/companion/obstructions`),
`rem:s57-defects`, `subsec:s57-scale`, `lem:s57-scale`, `subsec:s57-poincare`,
`def:s57-sigma-tilde`, `lem:s57-poincare`, `rem:s57-sketch`, `subsec:s57-counter`,
`prop:s57-no-gradient`, `prop:s57-indicator`, `prop:s57-onescale`,
`subsec:s57-verdict`, `thm:s57-independence` (+ `it:s57-A/B`), `cor:s57-refuted`,
`rem:s57-correct-translation`, `subsec:s57-family`, `rem:s57-family`,
`subsec:s57-crw`, `rem:s57-crw`, `rem:s57-ab`, `subsec:s57-ledger`,
`rem:s57-notclaimed`, plus equations `eq:s57-bmophi/tilde/sigmalog/osc/implies/survives`.

### Verification (worktree), with the baseline measured this session, not assumed

| Check | Baseline (`19f4d14`, unmodified) | After |
|---|---|---|
| `python3 scripts/check_proof.py` | exit 0, all 8 checks pass | **exit 0**, all 8 pass (756 labels, 529 `\ref` targets, 61 `\cite` keys, all resolved) |
| `pdflatex` × 3, `-halt-on-error` | exit 0, 0, 0 | **exit 0, 0, 0** |
| TeX errors (`^! `) | 0 | **0** |
| `LaTeX Warning` | 0 | **0** |
| hyperref "Token not allowed in a PDF string" (the bookmark warnings) | **109** | **109** — unchanged (new section title uses `\texorpdfstring`) |
| Overfull `\hbox` | 64 | 65 (one new, 2.4pt; the document already carries several of 40–46pt) |
| Underfull `\hbox` | 2 | 2 |

**Not committed, not pushed.** Worktree path for independent verification before merge:
`/tmp/claude-0/-home-user-tfirstframework/5dc09b30-469f-5a15-85ac-3ac7e0d30549/scratchpad/wpa-s57`
(detached HEAD at `19f4d14`; `git worktree remove` it after merging or discarding).

---

## 7. Ready-to-paste `PROGRESS.md` entry

*(Not applied — paste as a new `## 🔑 KEY FINDING — S57` block after the S56 block, if
S56's block has been pasted, and update "Active milestone".)*

```
## 🔑 KEY FINDING — S57 (WP-A): Grujić's log-weighted bmo hypothesis is NOT
## σ*-decay in other clothes — the two are logically independent — but it is a
## fourth member of the (AX)/saturation family and is not cheaper; the
## commutator/CRW step is CORRECT and is a real new tool

**Verdict: the S56 dictionary is REFUTED as stated.** Proof-document content
written in an isolated worktree at `19f4d14`; nothing committed.

**(I) The dictionary, settled.** Poincaré runs gradients → oscillation and does
not run back. What survives is one strict one-way implication, for an
UNRESTRICTED SCALE FAMILY: if `ρ^{-1}∫_{B_ρ(y)}|∇ê|² ≤ K|log ρ|^{-2}` for every
`B_ρ(y) ⊆ B_{r*}`, then `‖ê‖_{bmo_{1/|log r|}(B_{r*})} ≤ ¼K^{1/2}`. At the top
scale this is S56's `σ* ≲ 1/(log M)²` — but only as a SUFFICIENT condition, only
for the unrestricted deficit, only at one ball. Three smooth divergence-free
counterexamples kill everything else: (C1) `bmo_φ ≤ 2` with `σ* → ∞`; (C2)
`σ* = 0` with `‖ê‖_{bmo_φ} = ∞`; (C3) `σ* → 0` with bmo_φ failing by `|log δ|`
at scale δ.

**(II) The `1_{|ω|≥M/4}` indicator is load-bearing, not cosmetic.** The
unrestricted deficit is GENERICALLY INFINITE — at any transversal zero of ω the
direction has an order-one jump at every scale. The indicator is what makes σ*
finite. Any unrestricted-oscillation hypothesis imported here needs a zero-free
core supplied separately.

**(III) Grujić's OWN critical direction bound is `σ* = O(1)`.** His
`|∇ξ| ≤ C₁/max(|x|,ℓ)` is exactly this document's borderline
(`rem:marginality`, Regime B of `prop:regime-B`), with `ℓ(t) ≍ r*` proved. He
assumes the borderline and buys the gain elsewhere.

**(IV) The commutator/CRW step HOLDS.** Verified link by link: the cancellation
is Constantin–Fefferman's in commutator dress; `φ` increasing gives
`‖ξ‖_{BMO(B_{2R})} ≤ Cφ(2R)`; Jones extension on a ball is legitimate; CRW for
Riesz transforms plus `[R_iR_j,b] = R_i[R_j,b] + [R_i,b]R_j` suffices (no general
CZ version needed), linear in `‖b‖_BMO`; Hunt interpolation reaches `L^{3/2,∞}`.
One cosmetic source error ("reflexive Lorentz spaces" — `L^{3/2,∞}` is not
reflexive). **This is the first estimate on record here extracting a gain from
the direction field without an `L^∞` CZ bound — inside Wall 1(d)'s permitted
class — but the smallness it consumes is not available from σ*.**

**(V) The companion paper does NOT propagate the hypothesis**, and says so: its
Thm 5.7/Cor 5.8 are conditional on four assumed hypotheses (an inflow smallness
condition (Hk); tangential strain THREE logs below critical; a cap condition; a
core-history modulus), and its own §7.5 exhibits two natural critical cores where
they fail — a point-concentrated tube core (critical self-strain, Prop 7.10) and
a Beltrami core (critical twist, `‖ξ‖_{bmo_φ} ≳ |log r|`, Prop 7.12).

**(VI) TWO LEDGER CORRECTIONS.** (a) S56 rendered Thm 7.4's hypothesis as
`ω ∈ L^∞_tL^{3/2,∞}`; it is Definition 2.1, which additionally supplies the
CONTAINMENT `A_λ ⊂ B_{Cλ^{-1/2}}` that `L^{3/2,∞}` does not give and §5 needs.
(b) **The Albritton–Bradshaw prior does NOT attach where S56 put it.** Their
negative verdict is about the *a priori* sparseness classes; their *first* half
positively reproves the conditional "sufficiently sparse ⟹ no singularity"
criterion, which is the only half `2607.08866` §7 uses.

**Wall 1, Wall 3, Wall 9, (AX) and saturation all UNCHANGED in classification.**
Full record: `SESSION-LOG/2026-09-13-S57-wpa-grujic-bmo-pricing.md` and
`\S sec:wpa-s57` (worktree, unmerged).
```
