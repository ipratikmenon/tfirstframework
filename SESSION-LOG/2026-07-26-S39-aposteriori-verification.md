# Session S39 — 2026-07-26

## Summary

Added a new layer4 numerical diagnostic, `layer4/aposteriori_verification.py`,
implementing an a posteriori regularity **verification-style diagnostic** in
the shape of the technique introduced by Chernyshenko–Constantin–Robinson–Titi
(2007) and developed by Morosi–Pizzocchero. This work is independent of the
S31–S38 analytical proof program (Bridge Lemma / σ* / H1 / H2 / Regime A-B) —
it does not touch `proofs/claim_a_3d_proof_attempt.tex` or any of the modules
owned by that thread this session (`layer3/*`, `alignment_bridge.py`,
`k_correlation.py`, `depletion_diagnostic.py`, `enstrophy_exponent.py`).

## Scope — stated explicitly, restated here for the record

This module **does not produce a proof**. A genuinely rigorous
computer-assisted a posteriori proof (CCRT / Morosi–Pizzocchero) requires
rigorous interval-arithmetic / validated-numerics bounds on the residual and
every norm entering the Gronwall estimate. Everything computed here is
ordinary IEEE-754 floating point. Three modeling choices are NOT rigorously
justified and are flagged at every call site in the module:

1. Gronwall amplification rate `A(t) = C·‖∇u_a(t)‖_{L²}` — defensible
   functional form, but `C` is not derived rigorously (`C=1.0` default).
2. `default_threshold` — a chosen ν-proportional scale, not a rigorously
   derived local-existence radius.
3. The residual's `∂_t u_a` is a time-differenced estimate between
   consecutive solver states (the honest residual of the discrete
   trajectory, including its own time-integration error).

Consequently the module never returns a "proof_certificate" — only a
`verification_margin` — and `run_verification_experiment`'s verdict means
"diagnostic satisfied", never "regularity proved".

## What was built

- `layer4/aposteriori_verification.py` (599 lines):
  - `ns_residual` — momentum-equation residual; spatial operator reused
    verbatim from `layer3.route2_3D._rhs_velocity` (nonlinear term, Leray
    projection, spectral viscous term — none reimplemented); `∂_t u_a` from
    time-differencing.
  - `estimate_dudt_time_diff` — the time-differencing helper (primary
    definition (a) from the task spec).
  - `residual_norms` — spectral L² and H^{-1} ((1+|k|²)^{-1} weight) norms.
  - `amplification_rate`, `default_threshold` — the two unrigorous modeling
    constants (SCOPE notes 1–2).
  - `verification_margin` — exact frozen-coefficient (exponential-Euler)
    Gronwall recursion; reduces EXACTLY to closed-form solutions for
    constant A, r (the correctness anchor).
  - `run_verification_experiment` — runs `layer3/route2_3D.py`'s own
    solver (`make_solver`, `taylor_green_ic`/`shear_layer_ic`, `set_ic`,
    `compute_cfl_dt`, `solver_step`, all imported) and tracks the
    diagnostic every `record_every` steps.
  - `log_result` / `_write_result` / `_ensure_db` — new tables
    `layer4_aposteriori_experiments` and `layer4_aposteriori_timeseries` in
    `results/results.db`, with retry-on-locked exponential backoff
    (matching `layer4/depletion_diagnostic.py`'s pattern).

- `layer4/test_aposteriori_verification.py` (28 tests, all PASS):
  - Residual near-zero on an exact NS solution (2D Taylor–Green vortex
    decay embedded in 3D, `u = e^{-2νt}(sin x cos y, -cos x sin y, 0)`,
    which is an exact solution of the full nonlinear incompressible NS
    equations, not merely the linear Stokes problem).
  - Residual grows under a random perturbation.
  - `residual_norms`: H^{-1} ≤ L² ordering on a mid-frequency field; both
    zero on the zero field; bad-dx rejection; immutability.
  - `verification_margin`: exact quadrature for A≡0 (`E(T)=r0·T`); exact
    Gronwall closed form for A≡a>0 (`E(T)=r0(e^{aT}-1)/a`); monotonic
    doubling; degenerate/mismatched/nonpositive-threshold/non-increasing
    ValueError paths; immutability.
  - Smoke tests for `run_verification_experiment` (both `tg` and `shear`
    ICs, N=16).
  - DB logging roundtrip + retry-on-locked backoff (monkeypatched flaky
    write, matching `test_depletion_diagnostic.py`'s pattern).

- Full `layer4/` suite: **347/347 PASS, zero regressions** (28 new).

## Production runs

Both run foreground/blocking (each individually exceeded the tool's 120s
default timeout and was auto-moved to a background task slot by the
harness — not something requested; I did not poll, and waited for the
completion notifications before proceeding), N=64, ν=1e-3, seed=42,
record_every=10, n_steps=500 (chosen to comfortably fit within a single
call; each run's actual wall time was ~155s, well under any budget concern).

| exp_id | ic | n_steps | wall_time_s | residual_L2 range | residual_Hm1 range | A_max | E_final | margin_final | verdict |
|---|---|---|---|---|---|---|---|---|---|
| EXP-L4-VERIFY-TG-001 | tg | 500 | 155.6 | [0, 0.0766] | [0, 4.04e-3] | 3.542 | 5.507e16 | 5.507e19 | PARTIAL |
| EXP-L4-VERIFY-SH-001 | shear | 500 | 153.7 | [0, 1.898] | [0, 0.0735] | 3.131 | 2.017e11 | 2.017e14 | PARTIAL |

Both logged to `results/results.db` (`layer4_aposteriori_experiments`,
52 timeseries rows each in `layer4_aposteriori_timeseries`).

**Honest reading of the numbers:** in both runs the actual residual stays
small and well-behaved throughout (TG: Hm1 peaks at ~4e-3 mid-run, decays
back down by t≈22.6; shear: Hm1 stays under 0.08 throughout, monotonically
decreasing after an initial transient) — there is no sign of the
underlying trajectory being poorly resolved or blowing up. The margin
nonetheless explodes to 1e14–1e19 because the naive `exp(∫A dt)` Gronwall
factor compounds `A(t) ~ O(1–3.5)` over a ~20-time-unit integration window
(`exp(A·T)` with `A·T ~ 30–80`), against a `threshold = ν = 1e-3` that is
tiny by construction. **This is exactly the kind of artifact the module's
SCOPE section warns about** — the margin is extremely sensitive to the two
unrigorous modeling constants (`C_gronwall`, `C_thresh`), and a PARTIAL
verdict here should be read as "this particular unrigorous choice of
constants doesn't close the loop over this run's time window," not as any
evidence bearing on regularity, blow-up, or the separate S31–S38 analytical
proof program's open questions. Shorter windows or a larger threshold would
likely flip the verdict to PASS without changing anything about the
underlying (well-behaved) trajectory — precisely illustrating why this
diagnostic is a precursor to a rigorous proof, not a proof itself.

## Files changed

- `layer4/aposteriori_verification.py` (new)
- `layer4/test_aposteriori_verification.py` (new)
- `results/results.db` (2 new tables, 2 experiment rows, 104 timeseries rows)
- `PROGRESS.md` (session log row + last-updated pointer)
- `SESSION-LOG/2026-07-26-S39-aposteriori-verification.md` (this file)

No existing file was modified apart from the above; `layer3/*`,
`proofs/claim_a_3d_proof_attempt.tex`, `alignment_bridge.py`,
`k_correlation.py`, `depletion_diagnostic.py`, and `enstrophy_exponent.py`
were left untouched, per instructions.

---

## ADDENDUM (Fable) — windowed correction of a coordinator spec error

**The original spec was wrong, and the -001 runs were vacuous.** Accumulating
exp(∫₀ᵀ A ds) over a ~20-time-unit window with A ~ O(1–3.5) gives ~e⁷⁰ ~ 1e30, so
the margin cannot be satisfied by *any* residual, however small. The PARTIAL verdicts
of EXP-L4-VERIFY-TG-001 / SH-001 (margins 5.5e19, 2.0e14) therefore carried **zero
information**. This was an error in Fable's specification of the amplification window,
not in the agent's implementation — recorded here because the program's convention is
that spec errors get named as explicitly as mathematical ones. Real CCRT-style
verifications are applied over short windows for exactly this reason.

**Fix implemented (S39, same session):** sliding/restarting short windows with the
error reset to zero at each window start, post-processed from the same recorded
trajectory at three window lengths. Module now 799 lines; 11 new tests (39 in file);
full layer4 suite **358/358, zero regressions**. New table
`layer4_aposteriori_windows` (+87 rows); `window_length` added as a column beyond the
original spec, since without it rows from the same exp_id at different W collide on
window_index.

**Results (N=64, ν=1e-3, seed=42, n_steps=300, ~117 s wall each):**

| IC | W | n_windows | fraction_satisfied | max Gronwall | median margin | max margin |
|---|---|---|---|---|---|---|
| tg (t_final 12.60) | 0.5 | 26 | 0.423 | 5.11 | 1.58 | 4.81 |
| tg | 1.0 | 13 | 0.154 | 28.6 | 5.93 | 32.0 |
| tg | 2.0 | 7 | 0.000 | 920 | 40.9 | 1022 |
| shear (t_final 11.46) | 0.5 | 23 | **0.957** | 3.37 | 0.0014 | 6.83 |
| shear | 1.0 | 12 | 0.917 | 9.81 | 0.0046 | 19.9 |
| shear | 2.0 | 6 | 0.833 | 79.2 | 0.031 | 161 |

Long-window reference margins: TG 9.13e12, shear 6.58e8 — vacuous, as predicted.
Verdicts: **EXP-L4-VERIFY-TG-002 PARTIAL**, **EXP-L4-VERIFY-SH-002 PARTIAL** (rule:
PASS iff some W has fraction_satisfied == 1.0).

**Honest interpretation (Fable):**
1. Windowing removed the vacuity — Gronwall factor fell from ~1e12 to O(3–5) — so the
   PARTIAL verdicts are now *informative comparisons* rather than foregone conclusions.
   But neither trajectory reaches a clean PASS at any tested W.
2. The shear-vs-TG difference is real but must not be over-read: it reflects the size
   of **our numerical residual relative to an arbitrary threshold** (C_thresh·ν, a
   modeling choice), not either flow's proximity to a singularity. Cross-IC comparison
   is meaningful; the absolute verdict is not.
3. **Suspected artifact, flagged not fixed:** shear's *very first* window (t_start = 0)
   fails at all three W. A one-sided time difference at step 0 gives a poor ∂ₜu_a
   estimate, so this is most likely the residual definition misbehaving at the temporal
   boundary rather than a physical feature. Worth a one-sided/second-order-startup fix
   before any of these numbers are cited.
4. Scope unchanged and still binding: floating point, unrigorous constant C in the
   amplification rate, threshold a modeling surrogate. This is a **verification-style
   diagnostic, not a proof of smoothness**, and no interval arithmetic is involved.

**Bearing on the analytical program: none.** This track is fully independent of the
S31–S38 Bridge Lemma / σ* / H1 / H2 line and of S39's profile rigidity. It does not
support or undermine any of them.
