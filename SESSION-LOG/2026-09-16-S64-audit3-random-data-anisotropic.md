# S64 — third strategic audit: random-data well-posedness and anisotropic BKM criteria (both negative; Cao-Titi found CZ-free but bottoms out one level up)

Date: 2026-09-16. Third audit for (A)/(B), following S56 (six angles) and
S59 (three angles); launched in parallel with S63 (book chapter, `book/`
only, no interaction with this session). Owner's direction, verbatim: "Yes
I would recommend one more audit" — explicit advance acceptance that a
third negative result is a legitimate outcome. Executed by a subagent under
full orchestrator supervision; every load-bearing claim below was
independently re-verified against primary sources (see Verification).

## Angle 1 — random-data / probabilistic well-posedness: clean negative, on a different axis from anything S56/S59 found

Three structurally distinct strands were checked from primary sources.

**Strand A** (Nahmod–Pavlović–Staffilani, arXiv:1204.5444, and a 2025
fractional-NS continuation, arXiv:2509.17171): Wiener-randomization of
initial data gives almost-sure **global weak solution existence** for data
*below* `L²` regularity (negative-Sobolev, "super-critical" on the data
axis). Independently confirmed via direct abstract fetch (title, authors,
and abstract text match exactly). **This never engages the strong/classical
solution question at all** — Leray weak solutions from `L²` data already
exist unconditionally and deterministically since 1934; these papers only
push how rough the *data* can be while keeping weak existence, which is
orthogonal to (A)/(B) (continuing an already-existing smooth solution
smoothly for all time). No `L^∞(∇u)`-type bound is ever needed, so this
never comes close enough to Wall 1 to price against it.

**Strand B**: the specific "random tensor / propagation of randomness"
machinery associated with the Deng–Nahmod–Yue lineage was searched for
directly; every application found is to dispersive equations (NLS, cubic
wave). No NS/Euler application exists in the literature currently
searchable. The subagent's own structural explanation (the mechanism spends
a probabilistic gain against a *dispersion relation*, which NS's parabolic
smoothing and Euler's absence of dispersion do not supply) is flagged
explicitly as its own inference rather than a sourced claim — noted here
and not treated as established.

**Strand C** (Földes–Sy, arXiv:2401.00332, "Almost sure global
well-posedness for 3D Euler equation..."): the closest thing found to an
actual probabilistic global-regularity theorem for a genuine vortex-stretching
fluid PDE. Independently re-fetched and read in full via `pdftotext`
extraction of the primary PDF (not just the abstract). Confirmed verbatim,
at the exact location the subagent cited: *"A central question we are
interested in, is whether the blow up is generic, or it occurs only for
well prepared initial data,"* immediately followed by the paper's own
stated mechanism — a fluctuation-dissipation regularization (stochastic
forcing `~√ε` + dissipation `~ε`, invariant measure, `ε→0` limit) explicitly
introduced, in the paper's own words, to "randomize the dynamics and
**presumably make the solutions generic**." This confirms the subagent's
reading precisely: the construction is deliberately about statistically
generic behavior; its own framing declines to say anything about
non-generic ("well prepared") data — exactly the kind of data a
hypothetical structured blow-up profile would be. It is also confirmed for
**3D Euler**, not Navier-Stokes.

**Verdict, Angle 1**: clean negative, and precisely characterized as
distinct in kind from S59's finding (which found the negation of NS's
difficulty baked into a hypothesis) — here the mismatch is an axis
mismatch: data-roughness (Strands A/B) or generic-vs-arbitrary (Strand C),
neither of which is the axis (A)/(B) asks about.

## Angle 2 — anisotropic Littlewood-Paley criteria: the most interesting single finding of three audits, precisely qualified

**Cao–Titi, arXiv:1005.4463** (published *Arch. Ration. Mech. Anal.* 202
(2011)), Theorem 1: a sufficient condition for global regularity stated
entirely in terms of **one entry** of the velocity gradient tensor,
`∫₀ᵀ‖∂u_j/∂x_k‖_α^β ds ≤ M` under an explicit `(α,β)` scaling condition.
Independently re-verified from the primary PDF (extracted via `pdftotext`):
title, authors, and the exact theorem statement (conditions (11)/(12))
match the subagent's report verbatim. **Independently confirmed the proof
mechanism directly**: grepped the extracted text for every
Riesz/Calderón/singular-integral/Fourier-multiplier keyword — zero hits
anywhere in the paper — and read the actual proof of Theorem 1 (lines
138-300 of the extracted text): it proceeds by taking the inner product
with `−Δ_h u`, integrating by parts, using the incompressibility condition
`∂₁u₁+∂₂u₂+∂₃u₃=0` algebraically to collapse the nonlinear term to
`≤ C∫|u₃||∇u||∇_h∇u|`, then closing via two elementary physical-space
anisotropic inequalities (Lemmas 2-3, proved from scratch by 1D-slice
Hölder/Cauchy-Schwarz, no operators beyond elementary calculus). **This is
a genuinely Calderón-Zygmund-free proof of a genuine global regularity
sufficient condition — a real, novel-to-this-programme technical fact,
independently confirmed by direct reading rather than by title.**

**But this does not supply a route to (A)/(B).** Condition (11)/(12) is
itself a Serrin-critical, BKM-strength space-time integrability hypothesis
on a single gradient entry — strictly stronger than what the Leray-Hopf
energy inequality alone gives, and verifying it unconditionally from smooth
data is exactly as open as the original problem. The subagent's citation of
a 2026 continuation (Guo–Wang–Xiong, arXiv:2609.03877, described as calling
the base one-component condition "still unresolved" as of 9 days before
this audit) was **not independently re-verified by the orchestrator** in
this pass — flagged here as an unverified supporting detail, not load-bearing
for the session's central conclusion (which rests on the independently
confirmed Cao-Titi proof mechanism and theorem statement alone).

**Verdict, Angle 2**: sharper than "bottoms out in CZ in the remaining
directions" (the a priori expectation stated in the handover, and the S57
Grujić pattern). Here it genuinely does not bottom out in CZ — it bottoms
out in requiring above-energy-class control by *any* means, which is a
restatement of the Millennium problem in one-component coordinates, not a
CZ artifact. Recorded as a kept technical fact for the ledger (in the S57
commutator-technique tradition), explicitly not a work package: there is no
known dictionary connecting this programme's own machinery to the Cao-Titi
hypothesis the way S57 built one for Grujić's bmo condition, so there is no
well-scoped pricing exercise on offer.

**Hygiene note, independently spot-checked**: the subagent flagged and
excluded arXiv:2606.11720 (Runlong Yu) as belonging to the same
unprovenanced, single-author preprint pattern S56 §1.7 already declined to
use (≥10 similarly-styled NS preprints from the same author in one month).
Not independently re-checked by the orchestrator, but this is a correct
application of an already-established programme policy, not a new claim
needing separate verification.

## Angle 3 — bounded recency rescan: nothing new, one item cross-confirmed

Independently re-verified the single substantive item: fetched the live
`navier-stokes.org/navier-stokes-problem-solved/` tracker page directly and
confirmed, verbatim, both halves of the subagent's quote — OpenAI's Sept 8
2026 forced-blowup announcement and Clay's Sept 11 response ("apparently
been settled"; evaluation "deliberately unhurried"), **and** the page's own
explicit statement that *"The unforced 3D Navier–Stokes global-regularity
question remains open: the announced forced construction does not settle
it,"* plus a second sentence found on independent re-fetch not quoted by
the subagent: *"The announced smooth-data blowup does not disprove Clay
alternatives A/B."* This is exactly the (C)/(D) claim already priced in
full at S60-S62 (now closed) and bears on (C)/(D), not (A)/(B); it does not
touch Wall 1.

## Verification performed by the orchestrator (independent of the subagent)

- Fetched arXiv:1005.4463's abstract page directly: title and abstract
  confirmed verbatim.
- Downloaded and `pdftotext`-extracted the full Cao-Titi PDF; grepped for
  Riesz/Calderón/singular-integral/Fourier-multiplier terms (zero hits) and
  read the actual theorem statement and proof (lines 138-300) directly —
  confirms the CZ-free claim independently, not merely re-stating the
  subagent's assertion.
- Downloaded and `pdftotext`-extracted the full Földes–Sy PDF; located and
  read the exact "generic... well prepared" sentence and its surrounding
  paragraph, confirming both the quote and its context are represented
  fairly (the paper's own stated purpose is to make solutions generic via
  randomization, not to address arbitrary/non-generic data).
- Fetched arXiv:1204.5444's abstract directly: title, authors, and abstract
  confirmed verbatim.
- Fetched the live navier-stokes.org tracker page directly: confirmed the
  OpenAI/Clay quotes and the "remains open" / "does not disprove A/B"
  statements verbatim.
- Did **not** independently re-verify: Guo–Wang–Xiong (arXiv:2609.03877)
  beyond taking the subagent's characterization on trust; the
  Deng-Nahmod-Yue dispersive-only claim (Strand B) beyond the subagent's
  own search description; the Chemin-Zhang continuation literature; the
  arXiv:2606.11720 author-pattern exclusion. None of these gaps affects the
  session's central, independently-confirmed findings (Cao-Titi's mechanism,
  Földes-Sy's genericity caveat, Nahmod-Pavlović-Staffilani's scope, and the
  tracker page's content).

No arithmetic or attribution error was found in the subagent's report on
the parts re-examined this session.

## Recommendation

**No new work package. Third consecutive audit (after S56's six angles and
S59's three) returns a negative on every angle tried, and this is the
outcome the owner explicitly anticipated and pre-approved.** Angle 1 is a
clean negative for reasons distinct in kind from S59's. Angle 2 is the most
technically interesting finding across all three audits — a genuinely
Calderón-Zygmund-free sufficient condition, independently confirmed by
direct proof-reading — but it does not close the gap, for a reason now
precisely stated rather than assumed. Angle 3 confirms nothing new has
appeared since S59/S60.

Combined with the internal toolkit exhaustion (S48-S53, S56 §2) and the
now-separately-closed (C)/(D) forcing-floor line (S60-S62), the honest
statement is: **the programme's productive search space for (A)/(B), at
its current scale and tools, is exhausted.** Standing options, unchanged
from S59: consolidate the existing analytical/numerical record as a
finished negative-result document (the book already underway does this);
treat (C)/(D) as a separate, now-closed research direction; or wait, at
zero ongoing cost, for an external breakthrough on Wall 1 itself — this
programme's pricing machinery (exercised at S57 on Grujić, and here on
Cao-Titi) is ready to evaluate one if it appears.

**Named failure modes for any future session, carried forward:**
- **H1** — angle 2's CZ-free finding is misread as progress toward (A)/(B).
  It is not; the Cao-Titi hypothesis is the unresolved problem in different
  coordinates.
- **H2** — Földes-Sy's growth bound is cited without its genericity caveat.
  Forbidden; the caveat is load-bearing and independently confirmed here.
- **H3** — this session's bounded recency check is treated as exhaustive.
  It is not (targeted search, not a systematic metadata sweep).
- **H4** — "exhausted" is read as "the mathematics is settled." It is not;
  it is a statement about this programme's own reachable moves.

## Standing discipline confirmed

- Session tag: pre-assigned S64 (concurrent with S63, the book-chapter
  session) to avoid the collision problem that has recurred six times;
  verified at session start (`git log --oneline -1` showed tip `5e436e9`,
  highest `SESSION-LOG` tag S62) — no collision.
- `book/`, `layer3/`, `layer4/`, `results/` untouched; no `proofs/` edit
  made; no forcing introduced into this programme's own arguments (out of
  scope for this audit entirely, and not touched).
- No worktree needed (no proof-document edit); nothing to merge.
