# Session S45 — 2026-09-10 — Scale-r* Diagnostics (S44 follow-up) + Jet Strain Sweep

Two threads closed in one job, per user direction ("dump the magnetism idea and go
forward with jet strain sweep and s44"). Numerical; no proof-document changes.

## What S44 asked for, and why S42-S43 didn't answer it

`rem:s44-numerics` named three quantities bearing on `conj:lc-s44`,
`conj:K-refined` and `target:band-absorption`, and explicitly warned that S42-S43's
`band_concentration.py` volume fractions are taken **over the whole periodic box**,
not over `Q(r*)`, and are therefore "motivation and not evidence":

1. `|E|/|Q(r*)|` and `|Band|/|Q(r*)|` — the volume hypotheses (NB-E) / (NB).
2. `c_K^band` — which `prop:s44-linfty-route` requires to satisfy a **decaying**
   bound.
3. The **𝓛-dependence** of both correlation constants — the sharp point, since
   S34-S35 measured `c_K ≈ 1` and explicitly could not distinguish `O(1)` from the
   stronger `≤ C/𝓛` the closure actually needs.

## What was built

- **`layer4/scale_diagnostics.py`** — all three quantities computed inside
  `B_{r*}(x*)` (the spatial section of `Q(r*)`), reusing `_ball_mask`,
  `compute_direction_field`, `grad_ehat_sq`, `grad_scalar_sq`,
  `_vorticity_from_velocity` and `bridge_diagnostics`. Adds `sobolev_H3_norm` and
  `L_factor` (`𝓛 = log(e + ‖u‖_{H³}/M)`, definition at proof-doc line ~5671), and
  `l_dependence_summary`, which regresses `log c_K` on `log 𝓛` over pooled samples.
- **`layer3/quasar_jet_IC.py`** — bipolar swirling-jet IC and strain sweep.

Two design decisions worth recording:

- **Energy is held fixed; the swept parameter is the strain-to-swirl *ratio*.**
  Sweeping raw strain amplitude would confound "more strain" with "more energy" and
  move every diagnostic for the trivial reason.
- **The tilt perturbation is mandatory, not cosmetic.** A perfectly aligned
  axisymmetric jet has `ê` constant, so `w = |∇ê|² ≡ 0` and every correlation
  diagnostic collapses to its zero-by-convention branch.

**Instrument calibration before use:** `l_dependence_summary` was tested against
planted synthetic power laws and recovers slope `−1.0` and `−0.5` exactly (r²=1);
`sobolev_H3_norm` matches its analytic value on a single Fourier mode. 39/39 tests
pass; full `layer4` suite 443/443.

## Result 1 — the 𝓛-decay test. Negative, and this is the session's headline.

Regression of `log c_K` on `log 𝓛`, pooled across samples:

| pool | c_K slope | r² | c_K^band slope | r² |
|---|---|---|---|---|
| Part 1 (TG/shear/adv × 2 viscosities) | **−0.196** | 0.050 | **−0.059** | 0.017 |
| Jet sweep only | +0.431 | 0.144 | +0.747 | 0.303 |
| All 12 runs combined | **+0.088** | 0.020 | **+0.031** | 0.004 |

**The slope is nowhere near −1.** The low r² is not a failure of the fit — it *is*
the finding: across `𝓛 ∈ [1.35, 4.68]`, a 3.5× span over which a clean `c_K ∼ C/𝓛`
law would force a 3.5× **drop**, `c_K` instead stays pinned at 0.92–1.21 and
`c_K^band` at 0.93–1.46 — roughly 20% variation. The measured constants are
**flat in 𝓛**.

So: the data is consistent with **boundedness** and inconsistent with a clean
`1/𝓛` law over the accessible range. That is a genuine negative signal for
`prop:s44-linfty-route` — the one closure route S44 identified as requiring *no*
volume information.

**Standing caveat (S37), at full force.** These are globally-regular *decaying*
flows, not the approach to a singularity. The conjecture concerns `𝓛 → ∞`, and the
accessible range is small. This does **not** refute the conjecture. What it records
is that no onset of the required decay is visible anywhere we can currently look —
and that the program has now been looking at `c_K ≈ 1` for four sessions (S34, S35,
S42-43, S45) without ever seeing the decay the closure needs.

## Result 2 — the volume hypotheses, first measurement at scale r*

| hypothesis | measured | verdict |
|---|---|---|
| **(NB-E)** `\|E\| ≥ κ₂\|Q(r*)\|/𝓛` | `fE` = 0.45–0.99, `𝓛` ≈ 1.7–3.5 → `fE·𝓛` ≈ 1.4–2.6 | comfortably satisfied |
| **(NB)** `\|Band\| ≥ κ₁\|Q(r*)\|/𝓛` | `fB` = 0.00–0.18 → `fB·𝓛` ≈ 0–0.6 | tight, and **fails outright in one run** |

The failure case is structurally interesting and was not previously recorded:
**shear at ν=1e-3 returned `fB = 0.0000` exactly — the transition band
`{M/8 ≤ |ω| ≤ M/4}` was empty inside `B_{r*}`**, all vorticity there lying above
`M/4`. Where that happens, the band-absorption question is locally *vacuous* while
(NB) fails, so the volume route is unavailable precisely there. `c_K^band` is
correspondingly degenerate-by-convention (0.000) in that run.

Net on the two routes S44 laid out: the `‖w‖_∞` route needs a decay we cannot
detect (Result 1); the volume route needs (NB), which is measurably fragile.

## Result 3 — jet strain sweep: the aligned corner probed directly, no breakdown

| strain ratio | `fE` | `fB` | `c_K` | `c_K^band` | σ* |
|---|---|---|---|---|---|
| 0.00 | 0.632 | 0.141 | 1.075 | 1.007 | 1.17e+01 |
| 0.25 | 0.609 | 0.159 | 1.091 | 0.999 | 1.49e+01 |
| 0.50 | 0.553 | 0.170 | 1.072 | 0.985 | 8.05e+00 |
| 1.00 | 0.459 | 0.159 | 1.026 | 0.970 | 6.71e+00 |
| 2.00 | 0.300 | 0.125 | 0.982 | 0.948 | **4.31e+00** |
| 4.00 | 0.405 | 0.136 | 0.918 | 1.457 | 2.15e+01 |

σ* falls from 11.7 (pure swirl) to a **minimum of 4.31 at strain ratio 2**, then
jumps to 21.5 at ratio 4. The descent confirms the physical picture that motivated
the experiment: increasing collimation drives the flow toward alignment, exactly the
`σ* → 0` corner that S30's Theorem E and S38 identify as the one a singularity must
approach. The jump at ratio 4 is consistent with the jet destabilising once strain
dominates swirl.

**And `c_K` stays 0.92–1.09 across the entire sweep, including at the σ* minimum.**
`c_K^band` stays ≈0.95–1.01 except at ratio 4 (1.457). So the aligned corner was
probed deliberately and the correlation constants did not degrade there — mild
positive evidence for the constants' robustness, and a null result for the hope that
this geometry would expose adversarial structure.

## Resolution caveat — flagged, not buried

The ball `B_{r*}` at N=64 contains only **81 grid cells** for Taylor-Green, so TG's
scale-r* fractions and correlations are coarse and should carry the least weight.
Shear has 147–179; the adversarial runs ≈1650–1750; the jet runs 1300–4000. The
better-resolved runs tell the same story, which is why the conclusions above are
stated at all — but a TG-only reading of these numbers would not be trustworthy.

## Experiments logged

12 rows in the new `scale_diagnostics_experiments` table of `results/results.db`:
`EXP-L4-SCALE-{TG,SH,ADV}-{001,002}` (ν = 1e-3, 5e-4) and
`EXP-L4-JET{000,025,050,100,200,400}-064`.

## Where this leaves the band/c_K numerical programme

It has returned what it can. The constants are O(1) and stable everywhere tested —
across three IC families, two viscosities, and a deliberate sweep into the aligned
corner. What is *not* in evidence is the 𝓛-decay that `prop:s44-linfty-route`
requires, and the alternative volume route rests on an (NB) that can fail outright.
Further runs of this kind seem unlikely to change the picture; what would change it
is either a genuinely different regime (which regular decaying flows cannot reach)
or an analytical argument.

## Files changed

- `layer4/scale_diagnostics.py` (new), `layer4/test_scale_diagnostics.py` (new)
- `layer3/quasar_jet_IC.py` (new)
- `results/results.db` (new table, 12 rows)
- `PROGRESS.md`, this session log

No proof-document changes; no existing module modified.
