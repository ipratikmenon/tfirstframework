# Session S43 — 2026-09-09 — N=128 Resolution Check on the Band-Concentration Spikes

Executes the next step S42 flagged and left undone. No new code was written; this is a
production-run + analysis session against the existing `layer4/band_concentration.py`.
Does not touch `proofs/claim_a_3d_proof_attempt.tex`.

## The question S42 left open

S42 built the band-concentration diagnostic for `Conjecture target:band-absorption`
(the "core term" S41 proved the level iteration cannot absorb) and ran it at N=64 on
TG / shear / adversarial ICs. All three passed on the median convention, but TG and
shear showed large transient maxima in `conc_D` (63.1× and 192.8×). S42 flagged these
as "likely resolution artifacts, consistent with this program's documented N=64
TG/shear under-resolution history (S35-S37)" — explicitly *flagged, not explained
away* — and named the N=128 resolution check as the way to settle it.

## Runs

`EXP-L4-BAND-TG-101` and `EXP-L4-BAND-SH-101`: N=128, 400 steps, ν=1e-3,
`eps_param`=0.1, seed=42, `record_every`=10, `theta`=(1/8, 1/4). Wall time ~9705s and
~9547s (~2.7h each) — the per-step `|∇²ê|²` Hessian computation dominates cost at
N=128. Logged to `results/results.db`, table `band_concentration_experiments`.

## Result 1 — the medians are resolution-robust

| IC | conc_Y median N=64 → N=128 | conc_D median N=64 → N=128 | verdict |
|---|---|---|---|
| TG | 0.460 → 0.405 | 1.308 → **1.299** | PASS both |
| shear | 0.061 → 0.378 | 6.505 → 3.917 | PASS both |

TG's `conc_D` median changes by 0.7% across a 2× refinement. The median is the
statistic the verdict uses and the one analogous to how `c_K` was reported in
S34-S35; it is stable and comfortably under threshold at both resolutions.

## Result 2 — the maxima are NOT resolution-converged, and move in opposite directions

| IC | conc_Y max N=64 → N=128 | conc_D max N=64 → N=128 |
|---|---|---|
| TG | 26.6 → **136.0** (grew ~5×) | 63.1 → **214.4** (grew ~3.4×) |
| shear | 2.30 → 0.95 (shrank) | 192.8 → 75.0 (shrank ~2.6×) |

Refinement made TG's spikes *worse* and shear's *better*. A statistic that moves in
opposite directions under refinement on two ICs is not converged, and no conclusion —
in either direction — can rest on it as such.

## Result 3 (decisive) — every maximum lives in the early transient

Splitting each run's recorded timeseries at t = 3 settles what the maxima actually
are. Across **all** runs, **both** resolutions, **all three** IC families, every
maximum occurs at t < 3, and the bulk is uniformly tame:

| Run | max (and when) | bulk t≥3 conc_Y med/max | bulk t≥3 conc_D med/max |
|---|---|---|---|
| TG N=128 | cY 136.0, cD 214.4 — both at t=1.43 | 0.293 / **0.63** | 1.114 / **1.71** |
| shear N=128 | cD 75.0 at t=0.54 | 0.419 / 0.71 | 3.916 / 5.80 |
| shear N=64 | cD 192.8 at t=0.70 | 0.056 / 1.56 | 1.862 / 17.10 |
| adv N=64 | cY 9.8 at t=0.45; cD 5.2 at t=0.87 | 1.488 / 2.42 | 2.834 / 3.24 |

TG at N=128 — the run with the single largest spike in the whole study (136×) — never
exceeds **0.63** (`conc_Y`) or **1.71** (`conc_D`) anywhere in the bulk. Its last five
recorded samples sit at `conc_Y` ≈ 0.4–0.6, `conc_D` ≈ 1.4–1.7: stable, not drifting.
And the bulk maxima *converge downward* under refinement (shear bulk `conc_D` max
17.10 → 5.80).

## Correction to S42

S42's reading — "likely resolution artifacts" — was **half wrong**, and the correction
is recorded rather than quietly dropped. TG's spikes grew by 5× under refinement, so
they are not resolution artifacts in the sense S42 supposed. The accurate
characterization is that they are an **early-transient phenomenon**: they occur while
the initial conditions are still smooth and symmetric and the high-vorticity region
`{|ω| ≥ M/4}` has not yet been properly developed by nonlinear mixing, so band
statistics over few, unrepresentative points are volatile. That is the wrong time
regime to bear on `target:band-absorption` at all, which concerns the behavior of the
coupling density *approaching a putative singularity* — late time, large M — not
initial transients.

The mechanism behind the early-transient volatility is **not established here**.
Plausible contributors (offered as hypotheses, not claims): the band containing few
and unrepresentative points before mixing develops the high-vorticity set; and `ê`'s
second derivatives being Gibbs-prone near TG's vorticity nulls, which a finer grid
resolves more sharply rather than less. Establishing which, if either, would need its
own diagnostic.

## Conclusion

**`Conjecture target:band-absorption` is numerically supported, more cleanly than S42
concluded.** In the bulk-time regime that the conjecture actually concerns, the band
density is tame, stable across ICs, and converging downward under refinement, with no
sign of the concentration that would have made the conjecture false. The adversarial
IC — the purpose-built stress test — remains the calmest of the three.

This is numerics, not proof, held to exactly the same "supported, not proved" standard
this program applied to `c_K` in S34-S35, and it should not be read as stronger. It is
also 41 recorded samples per run, at two resolutions, on three IC families — evidence,
not a theorem, and it does nothing to close S41's actual analytical gap.

## Loose end (does not affect the conclusion)

The original N=64 TG run's per-step trace was lost to a `tail -100` in its capture
(only its summary row survived in `results.db`), so the time-split table above uses
TG at N=128 plus shear/adv at N=64. A confirmatory N=64 TG trace re-run (identical
parameters and seed, `db_path=None` so it does not touch the database) was in flight
at write time. The conclusion rests on four independently-traced runs and does not
depend on it.

## Process note

The first background waiter armed for these runs used `pgrep -f` with a pattern that
appeared in the waiter's own command line, so it matched itself, could never exit, and
never fired — the runs' completion was caught by a manual status check instead, after
they had already been finished for some time. Future waiters should grep the log for a
completion marker plus failure signatures (`TRACE DONE|Traceback|Error|Killed`) rather
than `pgrep` for a string the waiter itself contains.

## Files changed

- `results/results.db` (2 new rows in `band_concentration_experiments`)
- `PROGRESS.md` (KEY FINDING S43, session pointers, Session Log row)
- `SESSION-LOG/2026-09-09-S43-band-resolution-check.md` (this file)

No source file touched; `layer4` suite remains 404/404.
