# S52 — WP3: the axisymmetric-with-swirl pivot

**Date:** 2026-09-11
**Branch / worktree:** `claude/wp3-axisym-s51` at `/home/user/wp3-worktree`, based on `9770c61` ([S50])
**Scope:** `proofs/claim_a_3d_proof_attempt.tex` only. No numerics. `layer3/`, `layer4/`,
`results/`, `book/` untouched (WP4 owns those concurrently).
**Deliverable:** `\S sec:axisym-s51` (§38 in the rendered PDF, pp. 152–164), plus 11 new
bibliography entries and three rows in the S47 walls table.
**Build:** `check_proof.py` exit 0 (634 labels, 436 `\ref` targets, 41 `\cite` keys, all
resolved). `pdflatex` ×3 clean: 0 errors, 0 LaTeX warnings, 109 hyperref bookmark warnings
— exactly the pre-existing S50 baseline, so the two new math-bearing headings
(`\texorpdfstring`-guarded) added none. 178 pp (was 164). Diff is **purely additive**:
1147 insertions, 0 deletions.

---

## 1. The central question, answered

**Question (WP3 deliverable 1).** In the axisymmetric class, ∇u is recovered from
(ω_θ, Γ) by 2D-type elliptic problems in (r,z), not by the 3D Biot–Savart operator. Is the
L^∞ bound on ∇u still costing 𝓛 = log(e+‖u‖_{H³}/M), costing less, or log-free?

**Answer: still 𝓛, and provably sharp. The pivot does not evade Wall 1.**

∇u has 8 non-zero entries in the cylindrical frame (`lem:s51-frame`), splitting into:

- **Swirl block** (3 entries: ∂_r u_θ, ∂_z u_θ, −u_θ/r). **Log-free, absolute constant,
  unconditional** (`lem:s51-swirl-block`):

      |∂_z u_θ| = |ω_r| ≤ M,   |u_θ/r| ≤ ‖ω_z‖_∞/2,   |∂_r u_θ| ≤ (3/2)M.

  The middle bound comes from ∂_r(r²·u_θ/r) = r ω_z, integrated from the axis — a
  **first-order ODE in r alone**, whose Green's kernel s ds/r² is locally integrable.
  This is a genuine gain and is strictly better than `lem:s48-biotsavart`, which is
  log-free but returns the wrong exponent (3/5); here both the logarithm and the exponent
  are right.

- **Poloidal block** (5 entries). Two relations (ω_θ and div u = 0) for five unknowns ⟹
  a genuinely 2D elliptic problem, no ODE to exploit. **Pays the full 𝓛, and this is
  sharp** (`prop:s51-poloidal-sharp`): an explicit smooth compactly supported
  axisymmetric divergence-free field, built by an **exact** stream-function embedding of
  the classical dyadic 2D Calderón–Zygmund counterexample, with

      ‖∇u‖_∞ ≥ c N ‖ω‖_∞,   ‖u‖_{H³} ≤ C(r*) 2^{2N},

  **supported inside B_{r*} around the axis point** — i.e. inside the critical cylinder
  where the programme needs the bound. No perturbation of the elliptic operator is
  involved; the embedding is exact, so nothing is lost to error terms.

- **On the axis only** (`lem:s51-axis-logfree`), the poloidal recovery *is* log-free:
  writing L_θ = r·m, m solves Δ₅m = −Ω where Δ₅ = ∂_r² + (3/r)∂_r + ∂_z² is the
  5-dimensional radial Laplacian, and on {r=0} the poloidal block equals 2|∂_z m|, an
  **order-(−1)** operator applied to Ω = ω_θ/r. Order-(−1) kernels (|X|^{-4} in R⁵) are
  locally integrable; order-0 ones are not. That is exactly where the logarithm lives.
  **But the price is one power of the axis weight:** the bound is against ‖ω_θ/r‖_∞, not M.

**Quantitatively, the trade is bad.** At the critical scale r* = M^{-1/2} the unconditional
bound r*‖Ω‖_{L^∞(B_{r*})} can be as large as M^{3/2}, versus M𝓛 for the logarithmic route.
Since M^{1/2} ≫ 𝓛, **the log-free route is strictly the worse of the two** unless the
scale-invariant structural hypothesis

    (AX)   |ω_θ| ≤ K M r/r*  on B_{2r*}    ⟺   ‖Ω‖_{L^∞(B_{2r*})} ≤ K M/r*

is supplied, under which ‖∇u‖_∞ ≤ C(K+c_I)M log-free **on the axis** (`prop:s51-ax`).
(AX) is scale-invariant with the same K, and needs no constant to beat 1 — so it is not a
marginal target in the `rem:s39-distance`(3) sense.

**But (AX) is not free, and that is stated in the document, not softened.** For smooth
fields (AX) is to leading order ‖∇ω‖_∞ ≲ M/r* = M^{3/2} near the axis, and the only
control of ∇ω this document has (`lem:local-enstrophy-typeI`) is an L² bound that itself
carries 𝓛. **(AX) is a no-concentration hypothesis of the same family as saturation, and
S50's lesson applies verbatim: it has not been shown to be cheaper than what it replaces,
and no such claim is made.**

**The recurring shape, named explicitly** (`rem:s51-verdict-log`): this is the *third*
time the programme has produced a log-free bound against the wrong right-hand side.
`obs:s48-wall1` gave ‖∇u‖_∞ ≤ C(T−t)^{-1}, log-free but not against M.
`cor:s50-restatement` gave back exactly the estimate already held.
`lem:s51-axis-logfree` gives log-free against ‖Ω‖_∞ ≥ M/r* ≫ M. In all three the
logarithm is neither paid nor avoided — it is converted into a right-hand-side defect of
equal or greater size. Wall 1 is structural, and this class does not evade it.

---

## 2. The second finding: this class's Type-I branch is already closed, by others

**`Fact fact:s51-ss` — Seregin–Šverák 2009 (CPDE 34, 171–201; arXiv:0804.1803),
Theorem 1.1, quoted:** if v ∈ L³(Q(z₀,R)) is an axisymmetric weak solution with pressure
q ∈ L^{3/2}(Q(z₀,R)) and

    ess sup_{Q(z₀,R)} √(t₀−t) |v̄(x,t)| < +∞     (v̄ = **poloidal** projection only)

then z₀ is a regular point. **No swirl-smallness. No far-field decay. Local.** The standard
Type-I bound |v| ≤ C(T−t)^{-1/2} implies the hypothesis trivially.

**Consequence (`prop:s51-h1h2-moot`):** the pincer of §sec:pincer-s32, with H1 and H2
proved in full, concludes "no Type-I singularity" — which per `rem:s39-distance` is its
ceiling. In this class that conclusion is a **2009 theorem**, proved for weak solutions,
locally, with no (N), no (R), no (A-up) and no σ*. **Proving H1 and H2 here would re-prove
a known theorem under strictly more hypotheses.** The pivot moves the programme onto ground
where its own goal has already been reached without it.

---

## 3. What the session actually contributes (and it is a dictionary, not a theorem)

**`lem:s51-w-identity` — exact decomposition of the coherence deficit:**

    w = |∇ê|² = |∇_{(r,z)}ê|² + (1 − ê_z²)/r² = |∇_{(r,z)}ê|² + (J² + Ω²)/|ω|²

with J = ω_r/r, Ω = ω_θ/r. The curvature floor (1−ê_z²)/r² involves **no derivative of ê
at all** — it comes from the rotation of the frame, (1/r)∂_θ ê = (1/r)(ê_r e_θ − ê_θ e_r).

Consequences:
- `cor:s51-floor`: at an on-axis singularity (CKN), σ* ≥ (r*)^{-3}∫(1−ê_z²). Small σ*
  means the vorticity is, on average, aligned with the symmetry axis.
- `cor:s51-e0`: the constant direction e₀ of the aligned horn is **forced** to be ±e_z.
- **`prop:s51-aligned`: alignment kills the swirl source exactly.** ê ≡ ±e_z on an open
  set ⟹ ω_r ≡ 0 ⟹ ∂_zΓ = −r ω_r ≡ 0 ⟹ ∂_z(Γ²)/r⁴ ≡ 0. The programme's aligned horn,
  built for unrelated reasons at S32–S39, coincides in this class with the exact condition
  under which the axisymmetric problem degenerates to its solved case
  (Ladyzhenskaya / Ukhovskii–Yudovich).
- **Failure-mode-3 check, passed and recorded** (`rem:s51-aligned-scope`): alignment forces
  ∂_zΓ ≡ 0, **not** Γ ≡ 0 — Γ = Γ(r) is admissible and genuinely swirling. The aligned
  horn sits exactly *on* the boundary between the swirling and swirl-free problems without
  falling into the latter.
- `rem:s51-gram`: `lem:direction-eq`, `prop:sigma-star-invariant` and
  `lem:gram-negativity` are Cartesian algebraic identities, so they transfer **verbatim**
  with nothing to re-derive and **nothing gained**. Only the floor is new.

---

## 4. The rediscovery check (failure mode 1) — confirmed to be a real risk

`rem:s51-rediscovery`: by the identity above, the curvature part of σ* **is** the
normalised pair (J,Ω) — the two scalars on which Chen–Fang–Zhang 2017 and Lei–Zhang 2017
build their entire analysis. The literature controls (J,Ω) **dynamically**, by energy
estimates on the closed (J,Ω,ω_z) system, in L^∞_t L²_loc ∩ L²_t Ḣ¹_loc. The σ* machinery
controls the same objects only geometrically and only after normalisation by |ω|², which
destroys the homogeneity that makes their estimate close. **On this class the literature's
handle on the same objects is strictly stronger.**

---

## 5. The intended target is not a target (deliverable 3)

Side-by-side, with exact hypotheses (`fact:s51-gamma-criteria`):

| Source | Hypothesis on the swirl | Conclusion |
|---|---|---|
| Ladyzhenskaya 1968; Ukhovskii–Yudovich 1968 | Γ ≡ 0 | global regularity |
| Chen–Fang–Zhang 2017 (DCDS 37, 1923–1939) | \|u_θ\| ≤ C r^{−1+ε} | global regularity |
| Lei–Zhang 2017 (PJM 289, 169–187) | Γ in the scale-invariant (δ*,C*) form-boundedness class; sufficient: \|Γ\| ≤ C₁\|ln r\|^{−2} near the axis; v₀∈H^{1/2}, Γ₀∈L^∞ | global regularity |
| Wei 2016 (JMAA 435, 402–413) | \|Γ\| ≤ C₁\|ln r\|^{−3/2} near the axis; v₀∈H², Γ₀∈L^∞ | global regularity |
| **free (max principle)** | ‖Γ(t)‖_∞ ≤ ‖Γ(0)‖_∞ | **—** |

**`prop:s51-boundedness`:** any criterion of the form "Γ bounded ⟹ regular", with no decay
toward the axis, is **equivalent to global regularity for the whole class**, because its
hypothesis is discharged for free by the maximum principle. So the brief's stated standard
for "a real advance" (replacing smallness by boundedness) is, in this class, the entire
open problem — not a target.

**No criterion is claimed, and none is produced.** The residual gap is exactly
|ln r|^{−3/2} ⇝ O(1): purely logarithmic, and narrowing by exponent (−2 → −3/2), not by a
change of method.

---

## 6. Type-II ledger (`rem:s51-ledger`, written in the `rem:s39-distance` style)

1. **Type-I branch closed, by others** (`fact:s51-ss`). The programme's ceiling here is
   already occupied.
2. **Type-II axisymmetric is entirely open, and nothing in §33–§39 reaches it.**
   `rem:s39-distance`(1) applies verbatim and without weakening: profile extraction, the
   rescaled family and r* = M^{-1/2} all use the Type-I rates essentially. **Transferring
   the S32–S39 machinery to this class has zero value against the only open branch.**
3. **The live tools against Type-II are rate-free and are not geometric.** The Γ criteria
   assume no blow-up rate at all — that is their real advantage over the Liouville
   literature, and why they are the state of the art. Entering this line means competing
   with Lei–Zhang and Wei on their own ground with their own tools.
4. **The residual is logarithmic in a *second*, independent sense.** ‖Γ‖_∞ is free;
   regularity is known under |Γ| ≲ |ln r|^{−3/2}. **The axisymmetric class does not have
   one fewer logarithm than the general problem — it has one more.** What it *does* have
   is that this second logarithm is **not known to be sharp** (no counterexample shows
   Γ ∈ L^∞ cannot suffice), whereas `prop:s51-poloidal-sharp` shows the first one *is*
   sharp. A provably immovable obstruction traded for an open one is the only honest
   argument for working here, and it is an argument about prospects, not progress.
5. **Cost of entry** (`rem:s51-torus`): 𝕋³ has no continuous rotational symmetry about a
   line, so this class does not exist there. The programme's standing assumptions
   (line 5563) are on 𝕋³; the whole axisymmetric literature is an R³ theory. Moving the
   framework to R³ with finite energy is a prerequisite, not a detail.

---

## 7. Verification record — what was re-derived, what was fetched, what was found wrong

**Re-derived from scratch (nothing taken from the handover):**
- The Γ equation ∂_tΓ + u_r∂_rΓ + u_z∂_zΓ = ν(Δ − (2/r)∂_r)Γ and its maximum principle —
  **handover statement confirmed correct** (`lem:s51-gamma`). Note the exact cancellation
  of the two r^{-2} terms (u_r∂_r(Γ/r) against the Coriolis term u_r u_θ/r), which is
  what makes Γ transported.
- The ω_θ equation (`eq:s51-omtheta`), derived from ∂_z(r-momentum) − ∂_r(z-momentum):
  the centrifugal source enters with a **plus** sign; the convective commutator gives
  (∂_r u_r + ∂_z u_z)ω_θ = −(u_r/r)ω_θ by incompressibility; the viscous part uses
  [Δ,∂_z] = 0 but Δ∂_r f = ∂_rΔf + ∂_r f/r².
- The Ω = ω_θ/r equation, giving source ∂_z(Γ²)/r⁴ — **handover statement confirmed
  correct** (`lem:s51-omega`), and cross-checked against Lei–Zhang's form −2(u_θ/r)J,
  which is the same quantity.
- The three vorticity identities of the setting — all confirmed.

**No error was found in the handover's statement of the setting.** One thing the handover
does *not* flag, and which does mislead, is recorded as `rem:s51-handover`: **"Type I"
carries two incompatible meanings in this subject.** This programme and KNSS use the
*temporal* bound |u| ≤ C(T−t)^{-1/2}; much of the axisymmetric literature (and survey
treatments of it) call the *spatial* bound |u| ≤ C/r "Type I", because that is what makes
the Γ equation critical. Neither implies the other. Every statement in §38 names which.

**Sources fetched, with what they actually say vs. what a summary would suggest:**
- **CSTY 2008** (IMRN 2008, rnn016 — full text read). Theorem 1.1 assumes
  |v| ≤ C*(r²−t)^{-1/2}. Since (r²−t)^{-1/2} ≤ (−t)^{-1/2}, this is **stronger** than
  temporal Type-I: it implies it *and* forces decay in r. So CSTY alone does **not**
  exclude all Type-I.
- **KNSS 2009** (Acta Math. 203 — full text read). Theorem 5.3 needs |u| ≤ C/r globally.
  Theorem 6.2 needs **both** |u| ≤ C(T−t)^{-1/2} **and** a far-field |u| ≤ C/r for r ≥ R₀,
  and the authors state explicitly that it **fails** without the second (take u = b(t)).
  **So the handover's "Type-I excluded with swirl: CSTY / KNSS" is, on those two sources
  alone, an overstatement.** The document's existing use of `KNSS2009` (`thm:knss`, the
  2D bounded-ancient Liouville theorem, and the note that KNSS also treats the 3D
  swirl-free case) was checked and **is accurate**.
- **Seregin–Šverák 2009** (arXiv:0804.1803 — full text read). **This is the correct
  primary citation for Type-I exclusion in this class**, and it is cleaner than either of
  the above: Theorem 1.1 is local, needs no decay, and imposes the Type-I bound on the
  poloidal part alone. Already in the bibliography; now used.
- **Lei–Zhang 2017** and **Q. S. Zhang's survey arXiv:2101.04905** (full text read) for the
  FBC, the |ln r|^{−2} corollary, the closed (J,Ω,ω_z) system, and the criticality
  observation (scaling gap 0 for ASNS — the one genuine structural upside of this class,
  recorded in `fact:s51-criticality`).
- **Wei 2016** (|ln r|^{−3/2}), **Chen–Fang–Zhang 2017** (r^{−1+ε}, source of J),
  **Ladyzhenskaya 1968**, **Ukhovskii–Yudovich 1968**, **Neustupa–Pokorný 2000**,
  **Chae–Lee 2002**, **Kukavica–Ziane 2006**, **Cao–Titi 2011** — exact statements
  confirmed and cited (`rem:s51-onecomponent` for the last four).

**One error found in my own work and repaired before finalising:** the first draft of
`prop:s51-poloidal-sharp` embedded the 2D stream function φ_N directly and claimed the
resulting field was compactly supported. It is not — the 2D Newtonian potential grows
logarithmically. Repaired by cutting off with χ ≡ 1 on B₁ ⊃ supp ω_N, which leaves every
derivative at the origin **exactly** unchanged and adds only an N-independent O(1) to
‖ω‖_∞ (because ‖ω_N‖_{L¹} + ‖ω_N‖_{L^∞} ≤ C uniformly). Geometry adjusted to
r₀ = r*/2, ρ = r*/64 so the support still lies inside B_{r*}.

---

## 8. Ready-to-paste PROGRESS.md entry

*(Not applied here — PROGRESS.md is a merge hazard with the concurrent WP4 session. Paste
this as a new `## 🔑 KEY FINDING — S52` block after the S50 block, and update
"Active milestone".)*

```
## 🔑 KEY FINDING — S52 (WP3): the axisymmetric pivot keeps Wall 1, and its
## Type-I branch is already someone else's theorem

**Verdict: negative on both counts the pivot was commissioned to test.**
`sec:axisym-s51` (proof doc 12,245 → 13,392 lines, 178 pp, compiles clean,
diff purely additive).

**(I) The CZ logarithm is NOT removed, and not weakened.** ∇u splits into a
swirl block (3 of 8 entries) — log-free, ≤ (3/2)M, absolute constant,
unconditional, by a first-order ODE in r (`lem:s51-swirl-block`, strictly
better than `lem:s48-biotsavart`, which is log-free but has the wrong
exponent) — and a poloidal block (5 entries), which pays the full 𝓛 and is
**sharp**: `prop:s51-poloidal-sharp` gives an explicit field, supported
INSIDE B_{r*} around the axis, with ‖∇u‖_∞ ≥ cN‖ω‖_∞ and ‖u‖_{H³} ≲ 2^{2N},
by an exact stream-function embedding of the 2D CZ counterexample. On the
axis alone the poloidal block IS log-free (order-(−1) kernel in the 5D
picture, `lem:s51-axis-logfree`) — but against ‖ω_θ/r‖_∞, which at r* can
reach M^{3/2} vs M𝓛. **The log-free route is strictly the worse of the two**
unless the scale-invariant structural hypothesis (AX) |ω_θ| ≤ KMr/r* holds
(`prop:s51-ax`), and (AX) is a no-concentration hypothesis of the saturation
family — NOT shown cheaper, and not claimed to be. **Third occurrence of the
same pattern** (obs:s48-wall1, cor:s50-restatement, lem:s51-axis-logfree):
log-free bound, wrong right-hand side. Wall 1 is structural.

**(II) Type-I in this class is already closed by Seregin–Šverák 2009 Thm 1.1**
(`fact:s51-ss`): local, weak solutions, q ∈ L^{3/2}, Type-I bound on the
POLOIDAL part alone, no swirl-smallness, no far-field decay. The pincer's
ceiling (rem:s39-distance) is "no Type-I singularity". So H1 + H2, proved in
full here, would re-prove a 2009 theorem under more hypotheses
(`prop:s51-h1h2-moot`). **Correction to the handover:** CSTY 2008 assumes
|v| ≤ C*(r²−t)^{-1/2}, STRONGER than Type-I; KNSS Thm 6.2 needs an extra
far-field |u| ≤ C/r and the authors say it fails without it. Neither alone
excludes Type-I. SS09 does.

**(III) No regularity criterion, and the target is not a target.**
`prop:s51-boundedness`: any "Γ bounded ⟹ regular" statement IS global
regularity for the class, since ‖Γ‖_∞ is free by the maximum principle. Gap
in the literature is |ln r|^{-3/2} → O(1) (Wei 2016; Lei–Zhang 2017 at −2),
purely logarithmic. **Rediscovery risk confirmed real** (`rem:s51-rediscovery`):
by `lem:s51-w-identity`, w = |∇_{(r,z)}ê|² + (J²+Ω²)/|ω|² — the curvature part
of σ* IS the normalised Chen–Fang–Zhang/Lei–Zhang pair, which they control
dynamically and far better.

**What IS new and correct (a dictionary, not a theorem):** the exact
decomposition `lem:s51-w-identity` with its derivative-free curvature floor
(1−ê_z²)/r²; `cor:s51-e0` (the aligned horn's constant direction is FORCED to
be the axis); and `prop:s51-aligned` — **alignment kills the swirl source
exactly**: ê ≡ ±e_z ⟹ ω_r ≡ 0 ⟹ ∂_zΓ ≡ 0 ⟹ ∂_z(Γ²)/r⁴ ≡ 0. Failure-mode-3
check passed and recorded: this forces ∂_zΓ ≡ 0, NOT Γ ≡ 0.

**Type-II ledger (`rem:s51-ledger`).** Type-II axisymmetric fully open;
§33–§39 use Type-I rates essentially so transfer value is ZERO against it;
the live tools (Γ criteria) are rate-free energy estimates on (J,Ω,ω_z), not
geometric. **The class does not have one fewer logarithm than the general
problem — it has one more**; its only honest advantage is that the second one
is not KNOWN to be sharp while the first provably is. Cost of entry:
𝕋³ has no axis, so the standing framework must move to R³ first.

**Verified independently before merge:** re-derived the Γ equation, the ω_θ
equation and the Ω equation from the momentum equations (handover's statement
of all three CONFIRMED correct — no error found in the setting); read CSTY,
KNSS, Seregin–Šverák, Lei–Zhang, Wei, CFZ and the Zhang survey in full rather
than in summary; found and repaired one error in this session's own
sharpness proof (the 2D stream function is not compactly supported — fixed by
a cutoff that leaves derivatives at the origin exactly unchanged).
check_proof.py exit 0; three pdflatex passes, 0 errors, 0 LaTeX warnings,
bookmark warnings at the 109 baseline.

**Wall 1 classification UNCHANGED: `open`, structural.** New walls-table rows
10, 10′, 10″.
```

---

## 9. Recommendation, recorded so it is not softened later

(`rem:s51-verdict`.) The axisymmetric pivot does not do what WP3 was commissioned to test:
it does not evade Wall 1 structurally, and its residual open branch (Type-II) is
unreachable by the machinery the pivot was meant to carry across. A session that wishes to
continue here should do so **on the literature's terms** — energy estimates on the closed
(J,Ω,ω_z) system, aimed at the |ln r|^{−3/2} gap — and should expect to be competing
directly with Lei–Zhang and Wei rather than extending this programme.
