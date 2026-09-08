# Session S36 — 2026-07-19 — C1/C2 Verification, Depletion Numerics

## Purpose

Verify the two couplings (C1, C2) flagged in S35's K-self-absorption identity —
decisive for H1 — and build the first H2-feeding diagnostic (depletion ratio).

## Analytical results — §verification-s36 (proof doc, now 6275 lines, 502/502 balanced)

**C2 verified.** α-weighted terms close via the a priori bound (S34);
𝓛²-degraded, allowed class.

**C1 verified modulo one classical item.** Gram negativity (Lemma
gram-negativity): the drift-curvature terms contract ∇m⊗∇m against the PSD Gram
tensor G_jk = ∂ⱼê·∂ₖê with a NEGATIVE sign — the most dangerous couplings are
good-signed, not merely bounded. Stretching/commutator couplings close via the
existing Biot–Savart split. The one remaining item: transition-annulus terms
(supported on {M/8 ≤ |ω| ≤ M/4}) require the standard De Giorgi nested-level
iteration — classical bookkeeping, recorded as "open-mechanical," distinct in kind
from the genuinely open items (H2, σ*-decay).

**Correction to S35 (self-audit again catching own error):** the "bonus" second
good term ν∬m²w²η² claimed in S35's K-self-absorption identity CANCELS EXACTLY
against the νw² feedback inside D_tw (Lemma w2-cancel: the +νwê term of the
direction equation contributes exactly +νw² to ½D_tw). The K-control conclusion
stands; the reinforcement claim was wrong and is amended. Critical quadratic is
handled solely by the second rung (S35), not twice.

**Net effect on H1:** K-self-absorption holds (absolute constant, no log, no
correlation hypothesis) modulo the annulus residual. Routine items (i) transport,
(ii) commutator, (iv) time-slice all discharged with short proofs. H1's estimate
layer is now: proved, except for (a) σ*-smallness itself and (b) the annulus
level-iteration (mechanical).

## Numerical results

**Depletion diagnostic** (new layer4/depletion_diagnostic.py, 499 lines, 18 tests,
program 301/301 zero regressions) — measures α(x*)/M, the quantity H2's depletion
claim needs to be small and falling:

| exp_id | median/max/final ratio | M range |
|---|---|---|
| EXP-L4-DEPLETION-TG-001 | 0.231 / 0.538 / 0.271 | [1.40, 20.9] |
| EXP-L4-DEPLETION-ADV-001 | 0.016 / 0.032 / 0.021 | [1.75, 1.96] |

TG timeseries (own analysis of raw rows): depletion_ratio anti-correlates with M
— 0.538 at M=1.40 (step 30, near M-minimum) down to 0.085 at M=18.4 (step 110,
near M-peak). theta_eff (defined for M>10) ranges −0.35 to +0.02, i.e. the
effective scaling exponent of α_sup_high in M is only ~0.15–0.52 — dramatically
below the naive exponent 1. Strong evidence for genuine depletion (θ ≪ 1/2 in the
naive-vs-depleted convention).

**Reynolds sweep resolution fix (N=128 TG) — completed during S37, addendum:**
TG-101/102/103 (ν=1e-3/5e-4/2.5e-4) all landed clean (tail_fraction 1.3e-5–2.6e-5,
vs 0.07–0.59 at N=64). Regression slope +0.0028 — flat, confirming the qualitative
N=64 adversarial trend. But **c_K_median ≈ 2.24 at N=128 vs ≈1.02–1.06 at N=64 for
the same ν** — under-resolution had suppressed TG's measured c_K by roughly 2×, not
noise but systematic spectral pile-up at the grid cutoff. TG now sits noticeably
above the adversarial series (≈1.12–1.16), not near it as the under-resolved
comparison suggested. **New open caveat:** a stray longer probe (TG-104, N=128,
full t=8.02) showed tail_fraction=0.0082 — even N=128 becomes mildly under-resolved
late in the TG cascade; resolution adequacy is time-window-dependent, not just an
(N,ν) property. Flagged, not resolved — a candidate for an N=256 check if the exact
magnitude of c_K becomes analytically load-bearing.

## Process note (important for future sessions)

Three separate instances this session of agents detaching long-running numerical
experiments into background shells/watchers and then ending their own turn — the
detached process sometimes survives (harness-level auto-backgrounding on a 120s
Bash timeout) and sometimes does not (agent explicitly used run_in_background/
nohup and the process died with the session). Mitigation used successfully: querying
results.db directly rather than trusting an agent's self-report, and — when a
process is confirmed alive via `ps aux` — wrapping a `while kill -0 $PID; do sleep;
done` blocking wait in the orchestrator's own Bash (not the sub-agent's) to get a
reliable completion notification. Recommend stating "no detaching, ever" even more
bluntly in future long-run prompts.

## Actions

- proofs/claim_a_3d_proof_attempt.tex: §verification-s36 + amendment remark
  inserted (502/502 balanced, all refs verified).
- layer4/depletion_diagnostic.py + tests new; 2 experiments logged.
- layer4/kcorr_reynolds_sweep.py extended for N=128 TG runs (in progress at
  session close; will complete and log independently).
- PROGRESS.md update dispatched.
- No commit — not a git repository.

## Immediately followed by S37 (same session, per user's "next step already")

See SESSION-LOG/2026-07-19-S37-*.md: the H2 localized-CF derivation, which found
that H2's remaining gap and the standalone σ*-decay gap are structurally the SAME
open problem — a real reduction in the program's open-item count.
