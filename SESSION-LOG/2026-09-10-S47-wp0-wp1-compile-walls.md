# Session S47 — 2026-09-10 — WP0 (Verification Substrate) + WP1 (Walls Document)

First execution session under `HANDOVER-S46-OPUS.md`. Opus subagent in an
isolated git worktree; merged only after independent verification.

## Headline

**The proof document compiles for the first time in the programme's history:
142 pages, zero errors, zero LaTeX warnings.** S39 discovered it had never been
compiled and that all prior "verification" was grep-based brace counting. That
gap is now closed.

## WP0 — verification substrate

**Toolchain.** `texlive-latex-base`, `-latex-extra`, `-science`,
`-fonts-recommended`, **plus `cm-super`**. The last is not optional and is the
kind of thing that wastes a session: without it the first run dies with
`pdfTeX error (font expansion): auto expansion is only possible with scalable
fonts` — microtype has no Type-1 CM fonts to work with — and the toolchain
looks broken when it is not. Now recorded in CLAUDE.md.

**Three defects beyond the handover's enumerated list.** The compile found
undefined control sequences the structural checker cannot see by design:

- `\colonequals` (line ~4584)
- `\fint` (~5254 and ~5912, three error instances)

**These three formulas were not typesetting at all.** Fixed with
`\providecommand`; `\fint` aliases the preamble's existing `\Xint`, which
builds exactly the averaged-integral glyph and is otherwise never used in the
document — the alias line had evidently been lost at some point.
`\colonequals` maps to mathtools' `\coloneqq`. Neither changes any formula's
meaning.

**Process finding worth carrying forward.** The agent's first grep for compile
errors reported zero on a run that had five: `-file-line-error` reformats TeX
errors as `file:line:` with **no leading `!`**, so `grep "^!"` silently misses
everything. Any future error-scraping must not rely on the leading bang.

**Known defects fixed as specified.** `observation` and `openproblem` declared
via `\newtheorem[theorem]` alongside the existing declarations, sharing the
theorem counter so existing prose references number correctly.

**All four dangling references resolved — none guessed.** The handover named
this as the trap (a guessed target injects a false cross-reference into a
mathematical document). Each site now carries a TeX comment recording the
evidence:

| was | resolved to | evidence |
|---|---|---|
| `thm:vorticity-uniform` | `thm:uniform-enstrophy` (4387) | citing line quotes the displayed conclusion verbatim |
| `prop:blowup-alignment` | `thm:blowup-alignment` (4827) | 18 lines below, in the subsection of that name; citing word also corrected Proposition→Theorem |
| `lem:GN-ehat` | `lem:gn-grad-ehat` (4702) | only lemma bounding ‖∇ê‖_∞; citing sentence displays its conclusion |
| `eq:NS-vorticity` | `eq:vorticity-eq` (4223) | the document's only evolution equation for ω |

Machinery for the alternative is in place regardless: `check_proof.py` honours
a `% CHECK-PROOF: ALLOW-DANGLING <label>` marker, so a genuinely
unidentifiable reference can be *declared* rather than guessed.

**`scripts/check_proof.py`.** Brace balance (`\{ \} \\ \%` and comments
handled), `\begin`/`\end` multiset **and nesting order via a stack**,
undeclared environments, duplicate labels, unresolved `\ref`/`\eqref`/`\cite`,
and the malformed `end{...>` typo class found in S44. Exits non-zero on
failure. Its docstring states explicitly that **green is necessary but not
sufficient** — it cannot see undefined macros, which is precisely how the three
defects above survived every previous "verification". Invocation documented in
CLAUDE.md's Test Commands.

## WP1 — the walls document

`\section{Routes and Walls: A Referenced Record}`, `\label{sec:walls-s47}`,
lines 9716–10348, immediately before the bibliography. Ten walls, each on the
fixed four-part template (route and hoped-for conclusion / exact obstruction
with label references / classification / what it forbids or permits), plus a
status table and two closing remarks.

Walls: 0 (S31 audit, all eight grounds), 1 (the CZ logarithm, three
appearances), 2 (the (R) relocation), 3 (the Type-II ceiling), 4 (level
iteration), 5 (c_K/band non-implication + the S44 corrigendum), 6
(non-discriminating numerics), 7 (the S35 retracted term), 8
(`prop:localized-cf` non-transfer plus the companion H1 pull-back
prohibition), 9 (the S40 de-forcing finding).

**Classifications kept deliberately honest**, which was the named failure mode
(narrative drift into self-defence):

- **W1 is OPEN, not proved-impossible.** No proof in the document says the
  logarithm cannot be removed — only that it is paid once and inherited
  everywhere. Recording it as a proved barrier would have overstated it.
- W8 is "proved impossible *as written*", with a localized substitute on ℝ³
  left open. W9 is "proved impossible *as an import*".
- W6's classification rests on an argument about the *form* of the conjectures
  (each conditional on a regime a resolved run cannot represent), stated as
  such rather than implied to be a theorem.
- The S45 𝓛-slope data is fenced in a paragraph marked "Outside this
  document", noting it does not refute the conjecture and that nothing in the
  document rests on it.

## Three overstatements caught — all in my own S46 work

The agent was told to verify every claim against the document rather than
trusting the brief or `PROGRESS.md`, and that where they disagree the document
wins. That instruction paid for itself three times. All three are corrected in
`HANDOVER-S46-OPUS.md`; all three are recorded here because the erroneous
versions are still visible in older `PROGRESS.md` entries.

1. **The S31 audit has eight grounds, not seven.** `sec:audit-s31` lists
   Defects 1–8. `PROGRESS.md`'s S31 block lists seven. The omitted one is D8
   (unproved coercivity in `lem:A-evol`, plus the L^∞ CZ bound) — which is
   **the first appearance of W1**. Dropping it was not neutral: the programme's
   central structural obstruction was named in the audit that started the
   rebuild, and the summary lost it.
2. **`conj:lc-s44` is not "a common weaker target"** relative to
   `conj:K-refined`. By `rem:s44-lc-reading`, *neither implies the other*. It
   is verbatim the c_K hypothesis — same normalisation, same `1/𝓛` decay —
   asserted on the larger region `{m ≥ M/16}`, supplying the *conclusion* both
   targets were designed to give. It is weaker than JLE-1, not than
   `conj:K-refined`.
3. **"c_K and band-absorption proved independent" overstates S44.**
   `rem:s44-scope` is emphatic: an admissible test configuration need not come
   from an NS solution, so the counterexamples do **not** establish
   independence for Navier–Stokes solutions — "that question is not decidable
   by any argument available here". What is proved is that no implication
   follows from (A1)–(A3). The S44 KEY FINDING and its Session Log row carry
   the unqualified phrasing and are annotated accordingly.

Also surfaced by the first-ever compile: **the hand-written numbers in the
§20–§22 prose are off by one from LaTeX's** (prose "Theorem 20.4" is LaTeX's
21.4). Pre-existing, out of scope, not repaired — but a reader following prose
numbers is sent to the wrong section.

## Verification I performed before merging

- **The safety check first.** The worktree was indeed stale (7737 lines, at
  commit `14c4514`, no `sec:relation-s44`) — exactly the hazard the brief
  warned of. The agent detected it and copied the live 9787-line file in before
  editing. I confirmed the merge is safe by diffing the worktree file against
  the **live** file, not against the worktree's own base: **4 lines changed,
  and they are precisely the four bad `\ref` sites**; everything else purely
  additive. All seven major sections (S31, S36, S37, S38, S39, S41, S44)
  verified present.
- Re-derived the three corrections against the document: counted Defects 1–8;
  read `rem:s44-scope` and `rem:s44-lc-reading` in full.
- Spot-checked the highest-risk reference resolution semantically: citing line
  4740 quotes `Z(r,z_0) ≤ D·r^{2/3}` verbatim, exactly what
  `thm:uniform-enstrophy` displays at 4387. Confirmed all four target labels
  exist and are unique.
- Ran `check_proof.py` myself (exit 0) and `pdflatex` myself (two passes exit
  0; the single log warning is the benign "Label(s) may have changed. Rerun",
  cleared by a third pass; the 55 overfull boxes are pre-existing typesetting).

**One correction to the agent's own report:** it states
`lem:local-enstrophy-typeI` is at line 5694 "in the pristine file" and that my
~5664 was wrong. In the pristine file it is at **5660**, with `𝓛` defined at
5671; 5694 is its position *after* the agent's preamble additions shifted line
numbers. Minor, but line numbers are used as evidence in this document, so it
is recorded.

## Build artifacts

`proofs/claim_a_3d_proof_attempt.log` and `.pdf` are tracked, committed in the
original import `14c4514`, and now predate roughly 5000 lines of the document.
They were **not** deleted — they are imported files, not mine to remove as a
side effect — but a stale committed PDF alongside a now-compiling source is a
hazard for any reader. `.gitignore` now covers newly generated artifacts.
**Open decision for the owner:** remove the stale pair, given a current
142-page PDF can be regenerated on demand.

## Status against the handover

WP0 ✅ · WP1 ✅ · **WP2 (settle hypothesis (R)) is next** — the decisive gate,
one session, literature-only. Its named trap: ε-regularity needs *smallness* of
scaled local energy; Type-I supplies only *boundedness*.

## Files changed

- `proofs/claim_a_3d_proof_attempt.tex` (+679 lines net; 4 `\ref` sites
  corrected; preamble declarations and `\providecommand`s; new
  `§sec:walls-s47`)
- `scripts/check_proof.py` (new), `CLAUDE.md` (Test Commands, toolchain and
  ALLOW-DANGLING conventions), `.gitignore` (TeX artifacts)
- `HANDOVER-S46-OPUS.md` (three corrections + the incomplete-defect-list note)
- `PROGRESS.md`, this session log

No code, no numerics, no `results.db` changes.
