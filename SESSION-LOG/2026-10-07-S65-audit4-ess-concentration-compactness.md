# S65 — fourth strategic audit: ESS/backward-uniqueness and concentration-compactness both restate the problem in a different critical currency

Date: 2026-10-07. Fourth audit for (A)/(B), following S56 (six angles), S59
(three angles), S64 (three angles) — twelve angles total, all negative.
Executed by a subagent under full orchestrator supervision; every
load-bearing claim below was independently re-verified against primary
sources (see Verification).

## Angle 1 — Escauriaza–Seregin–Šverák (ESS) and backward uniqueness: clean negative, a genuinely new category

Primary sources read in full: Kenig–Koch, arXiv:0908.3349 (the cleanest
arXiv-available rendering of ESS's machinery, since ESS 2003 itself
predates arXiv); Gallagher–Koch–Planchon, arXiv:1012.0145; Tao,
arXiv:1908.04958.

**Kenig–Koch's Theorem 0.1** (critical-norm boundedness in `Ḣ^{1/2}`
implies global smoothness) is confirmed, independently fetched and quoted,
to carry the authors' own parenthetical: *"(Theorem 0.1 of course follows
from the result in [14])"* — [14] being ESS 2003. This is an alternative
*proof* of an already-known theorem, via Kenig–Merle-style
concentration-compactness/rigidity (the first such application to a
parabolic equation), not new mathematical content.

**The three-step rigidity architecture** (verified by the subagent against
the primary text and spot-checked here): (1) profile decomposition
(Gallagher) extracts a critical element — pure functional-analytic
compactness, no CZ sup-norm step; (2) the critical element is compact in
`L³` — uses only the *assumed* a priori critical-norm bound, not an
`L^∞(∇u)` bound; (3) rigidity forces the element to vanish, via a CKN-type
local smallness criterion, a backward-uniqueness theorem (Carleman-based),
and unique continuation. None of these three is a Biot-Savart recovery of
`∇u` from `ω` — confirmed by the subagent's step-by-step reading and not
contradicted by anything the orchestrator checked.

**Why this still doesn't help (A)/(B) — the key structural point,
independently confirmed via the Kenig-Koch and Tao abstracts.** ESS's
theorem is an *iff*: blow-up at `T*` `⟺` `limsup‖u(t)‖_{L³}=∞`. So (A)/(B)
is **logically equivalent**, not merely implied, to "`‖u(t)‖_{L³}` stays
bounded for every smooth datum, for all time." The entire
backward-uniqueness apparatus exists to prove that equivalence, not to
decide the equivalent statement — and that statement is exactly as hard as
the original, because `L³` is the critical Lebesgue space and the only
known global a priori control is the strictly *subcritical* Leray energy.
**This is the identical abstract shape S59 found for the programme's own
founding thesis**: a subcritical global quantity, a supercritical or
critical target, no known bridge.

**Quantitative form, independently confirmed.** Tao's arXiv:1908.04958
abstract, fetched directly, confirms verbatim: *"the critical L³_x norm
must blow up at a rate `(log log log 1/(T*−t))^c` or faster"* — and this is
already on this programme's own ledger (S56 §1.4(b), "Tao's triple-log").
The quantitative ESS line is not new territory for this programme; it is a
previously-logged item now correctly attributed to its ESS parentage. It
quantifies an orthogonal statement (the `L³` norm), not Wall 1's `∇u`
sup-norm.

**A genuine internal precedent, independently located and confirmed
verbatim.** `proofs/claim_a_3d_proof_attempt.tex` lines 7215-7219 already
record a prior attempt to invoke ESS backward uniqueness directly, as
candidate mechanism (β) for the open gap (A-up): *"($\beta$) Backward
uniqueness in the sense of Escauriaza–Seregin–Šverák
\cite{EscauriazaSereginSverak2003}, applied to the difference between
$\omega_U$ and its projection onto $e_0$; the parabolic operator is right
but the required decay at spatial infinity is not available for a
profile."* The orchestrator independently located and read this passage
directly (not taking the subagent's citation on faith) and confirms it
verbatim. This is a third, independent confirmation — from inside the
programme's own prior work — that this specific tool does not transfer
cleanly into this programme's machinery, worth recording so no future
session rediscovers it as new.

**Verdict, angle 1**: clean negative, in a precisely new category — not
"Wall 1 in different vocabulary" (the technique genuinely never uses CZ
sup-norm recovery) and not "cannot reach (A)/(B) for a structural reason"
(it is a genuine iff for the right quantity), but *"a different, equally
hard a priori-bound problem, logically equivalent to (A)/(B) itself for the
critical `L³` norm rather than Wall 1's `∇u` sup-norm, with the same
subcritical-energy-vs-critical-target scaling gap as the cause in both
cases."*

## Angle 2 — concentration-compactness / critical element: same technology as angle 1, same bucket

Kenig–Koch's own paper *is* the concentration-compactness paper the
handover asked about; GKP is its direct `L³`-critical sequel. Neither
paper reaches unconditional regularity — both reach ESS's *conditional*
criterion by an alternative route. GKP's genuine byproduct, independently
confirmed verbatim: *"assuming a singularity-producing initial datum for
Navier-Stokes exists in a critical Lebesgue or Besov space, we show there
is one with minimal norm, generalizing a result of Rusin and Sverak."* —
exactly the handover's anticipated shape ("if blow-up occurs, a minimal
blow-up solution with property X exists," with existence itself left open).

**The structural parallel to this programme's own S48/S50 work** — checked
and confirmed with one honest qualification. S48 proved a genuine, clean
escape from Wall 1's specific CZ move on the blow-up-profile branch
(`thm:s48-R`), but S50 then proved the residual quantity ("saturation") is
*exactly equivalent* to Wall 1 (`thm:s50-equivalence`), closing the loop.
The Kenig–Merle program exhibits the same two-stage shape: the "rule out
the compact limit" step is the part that gets cleanly proved (step 3 above,
paralleling S48's clean restart), and the genuinely unresolved part is the
step *before* it — the a priori compactness/boundedness hypothesis needed
to extract the limit object at all (critical-norm boundedness, paralleling
"saturation"). **The orchestrator's independent check of this claim's
sharpest sub-point**: the subagent reported that, unlike the S48/S50
profile-route case, no proof exists closing this gap back into Wall 1 — no
known dictionary between `‖u‖_{L³}` bounded and `‖∇u‖_∞≤KM`. This is
reported as an *absence of a proof*, correctly flagged by the subagent
itself as weaker than S50's proven equivalence and explicitly not to be
assumed either way by a future session. The orchestrator did not attempt
to construct such a dictionary either; this remains an open question about
whether the parallel to S48/S50 is exact or only structural.

**Independently re-verified, correcting two unreliable automated
summaries.** The subagent's claim that Jia–Šverák (arXiv:1201.1592) proves
"compactness of the minimal-data set conditional on the same unresolved
existence question" was checked directly by the orchestrator. Two
automated fetches of this paper gave inconsistent, partially wrong
summaries (one said the paper doesn't address singularity formation at
all; both missed the compactness result). The orchestrator downloaded the
actual PDF and extracted text via `pdftotext` directly: confirmed the
paper's own words, *"we also recover the compactness of the set of
'minimal blow up initial data' in `L³(R³)` modulo translations and
scalings,"* and its Theorem 4.1, *"Suppose `ρ_max<∞`. Then `M` is nonempty,
and moreover, `M` is compact,"* where `ρ_max` is defined (line 747 of the
extracted text) as the threshold below which no blow-up occurs — i.e.
`ρ_max<∞` **is** exactly the hypothesis that blow-up occurs for some
critical-norm-bounded datum, the open existence question. The subagent's
characterization is confirmed correct; the automated single-fetch
summaries were simply unreliable here, a useful reminder to fall back to
direct text extraction when a claim is load-bearing.

**Dated literature search** (arXiv API, reported by the subagent): nothing
found for 3D NS specifically via concentration-compactness/critical-element
methods since GKP/Jia–Šverák (2011-2015); not independently re-run by the
orchestrator, but consistent with the general sparsity of hits found
elsewhere this session.

**Verdict, angle 2**: clean negative, same underlying technology as angle 1
— reduces to a different, structurally analogous open rigidity question,
not proven (or disproven) to reduce further to Wall 1 itself.

## Recency rescan (three-week window, 16 Sept – 7 Oct 2026)

Independently re-confirmed the single most load-bearing item: fetched the
live `navier-stokes.org` tracker directly. Confirmed verbatim: *"The
announced construction relies on a carefully chosen smooth external force.
It does not decide whether every smooth flow without external forcing
stays smooth forever (alternatives A and B),"* with the page's most recent
substantive review dated 29 Sept 2026, no retraction, no official
verification, no prize awarded, Clay's "deliberately unhurried" position
unchanged. This is exactly the already-closed (C)/(D) claim; not reopened
per the handover's instruction.

Not independently re-run by the orchestrator: the subagent's dated arXiv-API
negative searches for BKM-logarithm weakening, Grujić/Tao-Barker/Cao-Titi-
lineage updates, and new ESS/concentration-compactness preprints in this
window (all reported negative). The subagent's flagged hygiene note (a
non-arXiv, unaffiliated "Dustyn Stanley" preprint claiming a full
resolution, correctly excluded per the S56 §1.7 precedent) was not
independently re-checked but is a correct application of established
programme policy.

## Recommendation

**No new work package. Fourth consecutive audit (fourteen angles across
four sessions) returns negative on every angle tried.** Both new angles
this session are genuinely different technique from anything tried before
(confirmed, not assumed) and both land in a single identifiable bucket:
each requires an a priori bound on a critical-or-above-scaling quantity
(`L³(u)` here; `∇u` sup-norm for Wall 1) that the only available global
control (subcritical Leray energy) cannot supply — the same gap S59
already found in the programme's own founding thesis, now confirmed to
recur in borrowed external machinery too.

Standing options unchanged from S59/S64: consolidate the existing record,
treat (C)/(D) as separately closed (S60-S62), or remain ready at zero
ongoing cost to price the next external technique if one appears — a task
this programme can now do efficiently, having done it three times (Grujić
at S57, Cao-Titi at S64, ESS/Kenig-Koch here).

## Standing discipline confirmed

- Session tag S65 verified at start (`git log --oneline -1` showed tip
  `7643461`, highest `SESSION-LOG` tag S64) and re-checked unchanged
  through the session.
- `book/`, `layer3/`, `layer4/`, `results/` untouched; no `proofs/` edit
  made (only read for verification); (C)/(D) line not reopened.
- No worktree needed; nothing to merge beyond this log and PROGRESS.md.
