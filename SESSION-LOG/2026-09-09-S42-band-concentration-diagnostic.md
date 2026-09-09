# Session S42 — 2026-09-09 — Band-Concentration Diagnostic for Conjecture target:band-absorption

Directly continues the S41 thread (per user direction: "let's focus on the core term
short for closure and look for next steps"). Does not touch
proofs/claim_a_3d_proof_attempt.tex — this is a numerical stress-test of the
conjecture S41 left open, not an analytical attempt to prove it.

## Motivation

S41 executed the De Giorgi level iteration sketched in S36 for the annulus residual
`\mathcal{A}_{\mathrm{ann}}` and proved it does NOT close — the recursion has no
contraction for any Young parameter choice, and a co-area identity shows no level
selection whatsoever helps. What survives is an explicit "core" (two families, S1+S2)
with an unconditional crude bound exactly one logarithm short of closure, restated as
`Conjecture target:band-absorption`: the coupling density on the transition band
`{M/8 <= |ω| <= M/4}` is not concentrated there relative to its value on the broader
region `{|ω| >= M/4}`. The proof doc states explicitly this is "the same analytic
type as Conjecture conj:K-refined" — the `c_K` decorrelation conjecture, which this
program has never proved, only numerically stress-tested (`layer4/k_correlation.py`,
S34-S35: `c_K ≈ 1`, no growth across a Reynolds sweep, read as "supported, not
proved"). This session builds the direct analog for the band-localized conjecture.

There was a live reason to suspect the conjecture could be FALSE rather than just
hard to prove: S41's structural finding was that the residual is fundamentally a
*level density* of a measure, with nothing in the argument preventing an adversarial
construction from concentrating structure exactly on a magnitude threshold. Numerical
testing was the natural way to get evidence either way before investing further
analytical effort.

## What was built

New `layer4/band_concentration.py`. Reuses (does not reimplement) `route2_3D`'s
solver (`make_solver`, `set_ic`, `taylor_green_ic`, `shear_layer_ic`,
`compute_cfl_dt`, `solver_step` — exact unforced Prize NS, f=0 throughout, same
solver used by every other layer4 diagnostic), `alignment_bridge`'s
`compute_direction_field`, `grad_ehat_sq`, `_validate_omega_grid`,
`_vorticity_from_velocity`, `adv_vortex_tubes_ic`, and `k_correlation`'s
`grad_scalar_sq`.

Two new spectral primitives (siblings of the existing gradient functions):
- `grad2_scalar_sq(field, kx, ky, kz)` — `|∇²f|² = Σ_{j,k}(∂_j∂_k f)²`, the
  Hessian-Frobenius-norm-squared, via `d_j d_k f = ifft(-k_j k_k · fft(f))`.
- `grad2_ehat_sq(ex, ey, ez, kx, ky, kz)` — `Σ_i |∇²ê_i|²`.

Core diagnostic, `band_concentration_diagnostics(wx, wy, wz, dx, theta_minus=1/8,
theta_plus=1/4)`, computing two ratios (default band matches Definition
`def:annulus-residual`'s `Ann = {M/8 ≤ |ω| ≤ M/4}` exactly; `E = {|ω| ≥ M/4}` matches
`k_correlation.py`'s `E` exactly):

- `conc_Y = mean_Band(|∇m|²w) / mean_E(|∇m|²w)` — tests `eq:s41-bandY`, the
  integrand of family (S1), the term S41 proved cannot be absorbed by any Young
  parameter choice in the level iteration.
- `conc_D = mean_Band(|∇²ê|²) / mean_E(|∇²ê|²)` — tests `eq:s41-bandD`, the
  angular-dissipation density family (S2)'s Young split reduces to.

Same degenerate-covariance convention as `c_K` (`conc := 0` when the denominator is
degenerately zero, not NaN/inf). Runner `run_band_experiment(ic, N, n_steps, ...)`
mirrors `alignment_bridge.run_bridge_experiment`'s structure (own loop, not a call to
it, to avoid running the solver twice) and logs to a new `band_concentration_experiments`
table in `results/results.db`. Verdict convention identical to `c_K`'s: PASS iff every
recorded value is finite and `median(conc_Y) <= 10` and `median(conc_D) <= 10` — the
same "stays O(1) or below" bar, not a proved threshold.

No existing file was modified.

## Tests

`layer4/test_band_concentration.py` — 24/24 PASS. Coverage: `grad2_scalar_sq` exact
on a single-mode analytic field (`f=cos(x)` → `|∇²f|²=cos²(x)` spectrally exact),
zero on a constant field, immutability; `grad2_ehat_sq` zero for a constant-direction
field, shape/non-negativity, immutability; band/E mask disjointness (overlap only at
the measure-zero boundary `m=θ₊M`); degenerate-zero-vorticity raises; invalid
`theta_minus >= theta_plus` raises; `conc := 0` convention verified on a
constant-direction field (denominator degenerately zero); required-keys and
immutability checks; default thetas match `Definition def:annulus-residual` (0.125,
0.25) exactly; solver-run smoke tests for all three ICs (tg/shear/adv); invalid IC
raises; exp_id/claim_id format; DB logging round-trip.

Full `layer4/` suite: 424 passed (the new 24 included). Four pre-existing failures
in `layer4/test_gamma_adversarial_search.py` — `AttributeError: module 'numpy' has
no attribute 'trapz'` (removed in this environment's NumPy 2.x; unrelated to this
session, not touched, out of scope).

## Production experiments (N=64, ν=1e-3, seed=42, n_steps=400 — same parameters as
the original S34 `c_K` sweep, for direct comparability)

Logged to `results/results.db`, table `band_concentration_experiments`:

| exp_id | ic | conc_Y median | conc_Y max | conc_D median | conc_D max | verdict |
|---|---|---|---|---|---|---|
| EXP-L4-BAND-TG-001 | tg | 0.460 | 26.57 | 1.308 | 63.15 | PASS |
| EXP-L4-BAND-SH-001 | shear | 0.061 | 2.301 | 6.505 | 192.77 | PASS |
| EXP-L4-BAND-ADV-001 | adv | 1.526 | 9.810 | 2.890 | 5.184 | PASS |

**Result: 3/3 PASS on the median convention.** All medians for both ratios stay
comfortably below the O(1)-to-10 range used throughout this program's "support, not
proof" diagnostics.

## Honest reading of the result — the part that matters

The medians support the conjecture. The **maxima do not tell a uniformly clean
story**, and this is reported directly rather than smoothed over:

- **TG**: `conc_D` spikes to 63.1 at its worst point — over 6x the PASS threshold,
  even though the median (1.31) is fine.
- **Shear**: `conc_D` spikes to **192.8** — nearly 20x the threshold — again with a
  comfortable median (6.5, still under 10, but closer to the line).
- **Adversarial**: the calmest by far. `conc_Y` max is 9.81 (just under threshold),
  `conc_D` max is 5.18. No large excursions. The per-step trace (recorded in
  `results.db`'s associated timeseries, not reproduced here) shows smooth convergence
  toward a stable band around 1.3-1.7 (conc_Y) / 2.8-3.2 (conc_D) at late time, rather
  than a single outlier event.

This program has an established, directly relevant precedent for exactly this pattern:
S35 found `c_K` at N=64 for TG was **resolution-suppressed by roughly 2x** relative to
the resolved N=128 value (confirmed in S37's resolution sweep) — TG and shear runs at
N=64 have a documented history of being under-resolved in this program (`tail_fraction`
warnings recorded in multiple prior sessions). `conc_D` is built from a
*second*-derivative spectral quantity (`|∇²ê|²`), which is intrinsically more sensitive
to marginal grid resolution than the first-derivative quantities `c_K` and `conc_Y` use.
The specific pattern here — large, transient spikes in the second-derivative quantity,
concentrated in exactly the two IC families with a known resolution history, while the
adversarial run (different numerical character, less prone to this particular
pathology) stays calm — is consistent with a resolution artifact. It is **not
confirmed** to be one; that would require rerunning TG/shear at N=128 and checking
whether the spikes shrink or persist, exactly mirroring S37's `c_K` resolution sweep.
That check was not run this session (time/compute budget).

**Bottom line:** mild numerical support for `Conjecture target:band-absorption`,
weighted most heavily by the adversarial run (the most relevant and least
resolution-suspect of the three), with an open, flagged, unresolved question about
whether the TG/shear spikes are real or numerical. This is not proof either way, and
should not be read as stronger than that.

## Next step (not done this session)

An N=128 TG/shear resolution check on `conc_D`, mirroring `kcorr_reynolds_sweep.py`'s
resolution-sweep methodology from S37, to determine whether the observed spikes shrink
(supporting the resolution-artifact reading) or persist (a genuine finding requiring
further analytical attention on `target:band-absorption`).

## Files changed

- `layer4/band_concentration.py` (new)
- `layer4/test_band_concentration.py` (new)
- `results/results.db` (3 new rows, `band_concentration_experiments` table, new table
  created)
- `PROGRESS.md` (new KEY FINDING S42 block, Last/Second-to-last session pointers,
  Session Log table row)
- `SESSION-LOG/2026-09-09-S42-band-concentration-diagnostic.md` (this file)

No file in `proofs/`, `layer1/`, `layer2/`, or `layer3/` touched.
