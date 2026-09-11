# Session S51 — 2026-09-11 — WP4: numerics redirect (axisymmetric + helicity)

Executed WP4a, WP4b and WP4c of `HANDOVER-S46-OPUS.md` in an isolated worktree
(`/home/user/wp4-worktree`, branch `claude/wp4-numerics-redirect`, base
`38247c1`). WP4d (true `(r,z)` solver, Hou-style search) not started — it was
gated on 4a, and 4a's answer is recorded below.

No `proofs/` or `book/` files touched (a concurrent agent holds
`claim_a_3d_proof_attempt.tex`). Not committed; for the owner to verify and merge.

## The standing obligation

Wall 6 (S46): the c_K numerics were never discriminating, because a regular
solver cannot exhibit the σ*-decay that S38 proved necessary for blowup, so it
measures c_K ≈ 1 — what *both* hypotheses predict on any flow it can produce.
Four sessions confirmed the non-information. Every experiment below therefore
carries a one-sentence statement, **written before the run**, of what a regular
flow could show that would distinguish the alternatives. All are recorded in
`results.db` in the `notes` column of each row.

---

## WP4a — axisymmetric diagnostics (`layer4/axisym_diagnostics.py`)

Γ = r·u_θ, ω_θ/r, and the source ∂_z(Γ²)/r⁴ localised in (r,z), on the
**existing** Route-2 JAX solver (`layer3/route2_3D_jax.py`) with an
axisymmetric swirling IC. float64 throughout: the Γ test is a cancellation
test and is vacuous in float32.

**Discriminating statement (MP).** Γ obeys a drift–diffusion equation with no
zeroth-order term, so for axisymmetric flow max|Γ| is non-increasing; a regular
flow has a definite falsifiable outcome, and any rise above the measured
axisymmetry defect and the float64 roundoff floor is a **solver bug**.

**Discriminating statement (source).** For smooth axisymmetric flow u_θ ~ c·r
near the axis, so Γ²/r⁴ = (u_θ/r)² is finite there; a regular flow shows
whether the literal r⁻⁴ form is evaluable at achievable resolution or is
swamped by 1/r cancellation noise — which gates whether a true (r,z) solver
(WP4d) would be measuring the term or its own roundoff.

### Calibration first (house rule)

| calibration | result |
|---|---|
| Γ vs exact Burgers closed form | rel err **1.05e-15** — PASS |
| ω_θ, source vs Burgers | identically 0 (exact) — PASS, but **degenerate**: a routine returning 0 would also pass |
| source vs **planted non-degenerate** field, `d_z(η²) = c²e^{-2r²/s²}sin 2z` | rel err **3.0e-14** — PASS |
| solver holds an exact **periodic** solution (ABC/Beltrami, E(t)/E(0)=e^{-2νt}) | energy rel err **2.1e-14**, shape **1.1e-14** — PASS |

**Honest caveat, stated rather than buried:** the Burgers vortex is *not*
periodic (its strain u_r = −ar/2, u_z = az is unbounded), so it is **not** a
steady state of this periodic solver and the solver cannot be expected to hold
it. The handover assumed otherwise. It is used here to calibrate the
*diagnostics* against closed forms, which is what the house rule needs; the
solver's "holds an exact solution" test is done on the ABC flow, which **is**
an exact periodic NS solution.

### Γ maximum-principle result — it fired, and it was not a solver bug

First run (N=64): max|Γ| rose **5.6%**, against an axisymmetry defect of 2.6%.
Chased as instructed. Localisation showed the rise sits in the **core**
(r ≈ 0.8–1.0), not the far field, and tracks the spectral tail fraction. The
resolution series settles it:

| N | max rel rise of \|Γ\|_∞ | D_axi | tail frac | verdict |
|---|---|---|---|---|
| 48 | 1.860e-1 | 4.93e-2 | 5.80e-3 | FAIL |
| 64 | 5.553e-2 | 2.64e-2 | 2.02e-3 | FAIL |
| 96 | 5.193e-3 | 7.52e-3 | 3.36e-4 | PASS |
| **128** | **1.332e-3** | 3.80e-3 | 6.64e-5 | **PASS** |

The violation falls **140×** over a 3× refinement while the spectral tail falls
87×. **No convergence exponent is fitted: 4 points is below this module's
`MIN_FIT_POINTS = 8`** (S49 rule — `results.db` carries a legacy PASS on a
2-point fit with r² = 1.0000).

**Conclusion: no solver bug. The Γ maximum principle is a working free
correctness test, and what it detects is loss of resolution** — it fails exactly
when the run is under-resolved and passes once it is not. It is now a standing
resolution monitor, with `tail_fraction` recorded on every snapshot.

Two genuine defects *were* found and fixed en route:
1. My first `axisymmetry_defect` binned in (r,z) and so conflated within-bin
   **radial** variation with azimuthal variation, reporting D ≈ 0.19 on an
   analytically axisymmetric field. Replaced by an azimuthal **Fourier**
   decomposition on a polar resampling: now 6.9e-4 on that field and 0.117 on
   a deliberately perturbed control.
2. `twisted_ring_IC` — see WP4b.

### Source near-axis (the WP4d gate)

Both branches are implemented: the literal `(1/r⁴)∂_z(Γ²)` and the
axis-regular `∂_z(η²)`, η = u_θ/r (an exact identity, since r is z-independent).
Their relative discrepancy on an **evolved** flow is **1.4e-14** at N=128, with
`r_loss_radius = 0.0` at every resolution — i.e. there is **no** near-axis
evaluability problem at these resolutions. **The WP4d gate is open on this
criterion**; the term can be measured rather than roundoff.

---

## WP4b — twisted ring and helicity diagnostics

`layer3/twisted_ring_IC.py` (ω = f(ρ)[ê_φ + (ι ρ/R₀) ê_χ], Biot–Savart in
float64) and `layer4/helicity_diagnostics.py` (H, h = u·ω, twist density
τ = ê·(∇×ê), stretching α = ê·S·ê, the pointwise correlation on {|ω| ≥ M/4}).

**Discriminating statement.** On a regular flow both |τ| and α are finite O(1)
fields on {|ω| ≥ M/4}; if twist locally depletes stretching then corr(|τ|,α) is
negative there, post-transient and pre-reconnection, with a sign that persists
as |ι| grows and across two resolutions — and a regular flow is equally capable
of showing a positive or null correlation, so either outcome informs.

### Named failure modes, discharged as tests

1. **Zero net helicity by symmetry** — H = 1.7e-18 at ι = 0, exactly linear and
   sign-flipping in ι, H ∝ Γ². Asserted before any run.
3. **Ring core under 4 cells (S40)** — `check_resolution` raises. It also
   revealed that N = 32 cannot resolve a ring *at all* under the joint
   constraints, so 48 is the smallest grid used.
4. **Bespoke finite-difference solver** — none; everything runs through the
   existing, unmodified `route2_3D_jax`.

### The parity trap (not in the handover; found here)

τ is a **pseudoscalar**, α a true scalar. If IC(−ι) is the mirror of IC(+ι) —
it is — then corr(τ,α) is **forced odd in ι** and corr(|τ|,α) **forced even**.
So "a correlation that flips sign with ι is not a mechanism" would, for the
signed correlation, be reading the mirror symmetry of the initial data.
Measured oddness/evenness residual: **3.9e-14** — machine precision. The
parity-correct measurement is corr(|τ|,α); the signed version is kept only as
this correctness check.

### Two real bugs found by controls

**(i) `twist_density` was contaminated at the signal level.** The first
implementation FFT-differentiated ê, which is a *unit* vector up to the edge of
the vorticity support and 0 outside — a jump of size 1, whose spectral
derivative rings across the box. On an **untwisted** ring, where τ must vanish
identically, it returned ±0.74 against a true-signal scale of ~1.5. Replaced by
the exact identity

  ê·(∇×ê) ≡ ω·(∇×ω)/|ω|²,

which differentiates ω (smooth, no truncation). Residual at ι = 0 dropped from
0.74 to **3e-4**. Caught by `test_untwisted_ring_has_near_zero_twist_density`.

**(ii) Self-intersecting ring geometry.** The first parameter choice had
ρ_cut = 3a > R₀, putting tube material on the axis where the toroidal frame is
singular. A third guard now rejects ρ_cut ≥ R₀. Defaults are R₀ = 1.35,
a = 0.55, ρ_cut = 2a — the joint constraints (a/dx ≥ 4, ρ_cut < R₀,
R₀ + ρ_cut < π) are tight and admit no solution with ρ_cut = 3a at N ≤ 64.

### The control failure, diagnosed rather than averaged past

The ι = 0 run gave corr = **+0.15 at N=64 and −0.42 at N=96** — a sign flip
under refinement. This was treated as a stop sign.

Diagnosis: **mean|τ| at ι = 0 is 2.1e-3 (N=64) → 9.1e-4 (N=96) — converging to
zero**, while the ι ≥ 0.5 signal is 0.72–1.95 and agrees between those grids to
**four significant figures**. A Pearson correlation is **scale-invariant**, so
it cannot report "there is no signal": handed a pure-noise τ it returns an
arbitrary number. **ι = 0 was never a valid null control** — that was a design
error in the control, not a defect in the ι > 0 measurements, and it does not
impugn them. Two fixes:

- `MIN_TWIST_SIGNAL = 1e-2`: a correlation is **refused** when mean|τ| on the
  mask falls below it. The ι = 0 rows in `results.db` are corrected to
  `INSUFFICIENT-DEGENERATE`.
- A **valid** null: `surrogate_null_correlation` — phase randomisation, which
  preserves |τ|'s power spectrum (hence its spatial autocorrelation) while
  destroying its pointwise alignment with α. A plain permutation will not do:
  it destroys the autocorrelation too and reports a null far too narrow for
  smooth fields (asserted as a test).

### Result, with the valid null (200 surrogates, t = 2.5)

| N | ι | mean\|τ\| | observed | null 95% band | z | verdict |
|---|---|---|---|---|---|---|
| 64 | 0.0 | 0.0021 | +0.206 | [−0.022,+0.026] | +16.5 | **DEGENERATE — refused** |
| 96 | 0.0 | 0.0009 | −0.420 | [−0.014,+0.015] | −57.1 | **DEGENERATE — refused** |
| 64 | 0.5 | 0.721 | **−0.308** | [−0.132,+0.135] | −4.77 | SIGNIFICANT |
| 96 | 0.5 | 0.721 | **−0.306** | [−0.115,+0.130] | −4.93 | SIGNIFICANT |
| 64 | 1.0 | 1.318 | **−0.181** | [−0.126,+0.144] | −2.76 | SIGNIFICANT |
| 96 | 1.0 | 1.318 | **−0.190** | [−0.133,+0.148] | −2.95 | SIGNIFICANT |
| 64 | 2.0 | 1.950 | −0.025 | [−0.057,+0.054] | −0.79 | NULL |
| 96 | 2.0 | 1.953 | −0.027 | [−0.061,+0.075] | −0.83 | NULL |

Note the ι = 0 rows are called "significant" by the surrogate test with huge z —
because the null for a noise field is extremely narrow. **A significance test
alone is not enough; a signal-magnitude guard is also required.** That is why
both exist.

### Verdict on the lemma candidate — PARTIAL, and stated plainly

**Twist does locally deplete stretching at moderate twist, and the measurement
is trustworthy: resolution-converged (1–5% between N=64 and N=96), significant
against a valid null, sign-stable (negative wherever the twist signal is
non-degenerate).**

**But the ι-dependence runs the wrong way for a mechanism.** A local
twist-depletion lemma requires depletion to *strengthen* with twist. It does the
opposite: −0.31 at ι = 0.5, −0.19 at ι = 1.0, **statistically indistinguishable
from zero at ι = 2.0**, at both resolutions. The strongest depletion is at the
weakest twist tested, and it is gone by the strongest.

So this is **not** the lemma the pivot wanted, and it should not be built on.
It is, however, **informative** — unlike c_K, this measurement could have come
out the other way on a regular flow, and the ι-scaling is a real, reproducible
negative constraint on the local-depletion picture. Effect sizes are small
throughout (|c| ≤ 0.31).

---

## WP4c — reconnection detection

Three independent flags in `detect_reconnection`: connected-component count of
{|ω| ≥ M/2} under **periodic** connectivity (union-find across wrapped faces),
a local minimum of the minimum inter-component distance, and a helicity jump.
No reconnection was detected in any twisted-ring run in the ι sweep, so the
pre/post split did not bite here — but the machinery is in place and tested, and
the sweep's common analysis window is cut at the earliest flagged time when one
occurs.

One bug found and fixed: the relative helicity-jump test divided by H(0), which
is **0 by symmetry** for the untwisted ring, so every ι = 0 run was flagged as a
reconnection. It now falls back to the field's own helicity capacity 2·E·M when
H(0) is negligible against it.

## Methodological fix: a common analysis window

Each run's own post-transient time differed across ι (0.0 for |ι| = 2, 0.8
elsewhere), so the first sweep compared correlations taken over **different**
time ranges. `iota_sweep` now recomputes every run's primary window on one
window common to the whole sweep; the per-run adaptive window is kept as
`primary_own_window`. HANDOVER §WP4b failure mode 2 (S43 early-transient
misreading) is handled by this plus explicit early/middle/late thirds.

## WP4d

**Not started** (out of scope, gated). 4a's gate criterion — near-axis
evaluability of the source term — came out **open**: the r⁻⁴ form agrees with
the regular form to 1.4e-14 with no loss radius.

## Files

- `layer4/axisym_diagnostics.py`, `layer4/test_axisym_diagnostics.py` (35 tests)
- `layer3/twisted_ring_IC.py`, `layer3/test_twisted_ring_IC.py` (24 tests)
- `layer4/helicity_diagnostics.py`, `layer4/test_helicity_diagnostics.py` (45 tests)
- `results/figures/wp4a_gamma_maxprinciple.png`,
  `results/figures/wp4b_twist_stretching.png`
- `results/results.db` — new tables `layer4_axisym_experiments`,
  `layer4_axisym_timeseries`, `layer4_axisym_calibration`,
  `layer4_helicity_experiments`, `layer4_helicity_timeseries`,
  `layer4_helicity_sweep`, `layer4_helicity_null_control`

## Tests

104 new (35 + 24 + 45). Full suite: **layer4 547 passed, layer1+layer2 318
passed, 0 failures.**

## claim_ids written

`wp4a-gamma-maxprinciple`, `wp4a-diag-calibration`,
`wp4b-twist-depletes-stretching`, `wp4b-twist-alpha-null-control`.
Every row carries a claim_id and a key_metric; the count of rows lacking a
claim_id was 0 before this session's writes and is 0 after.
