# S63 — book chapter: "The Forcing-Floor Class Theorem, and Where It Stops" (S60-S62)

Date: 2026-09-16. Run in parallel with S64 (third strategic audit),
`book/` only, no interaction between the two sessions. Executed by a
subagent under full orchestrator supervision; every claim below was
independently re-verified by the orchestrator before merge (see
Verification).

## What was written

New chapter `book/chapters/17c-forcing-floor.tex` (578 lines), placed in
`\part{The Landscape}` immediately after `17b-pressure-testing`, one
`\include` line added to `book/tfirst-book.tex`. The chapter documents the
S60-S62 arc as a complete, closed three-session story: the pivot from
(A)/(B) to a bounded study of (C)/(D) forcing requirements; S60's general
class theorem (`thm:cd1-general-floor`, three equivalent forms, the rate
corollary, the correction to S55's own record about what exact
self-similarity actually buys, the hypothesis audit pinning axisymmetry as
the theorem's one absolute boundary, the L^q floor and its scope, and
Thread 2 reported and closed as a survey); S61 as a scoping-only handover;
S62's two deliverables (the fractional extension's precise H2 failure, and
the Zhang-instance search landing in the theorem's degenerate corner). One
TikZ figure (the six-step arc, drawn as a straight chain since nothing in
it doubled back) and two tables (the three candidate papers checked; the
rate corollary applied to both known instances).

Two small, additive "Amended at S63: ..." patches were made: `18-walls.tex`
(Wall 9's S55 forcing-number paragraph, noting it is now one instance of
S60's class theorem, with S62's fractional/second-instance result and the
closure decision) and `22-where-it-stands.tex` (the ranked "other live
targets" list, noting the ranking is unchanged but is no longer what the
programme is actively doing, and pointing to the new chapter). Neither
list/ordering itself was touched; both patches are purely additive
paragraphs appended after existing content.

## Verification performed by the orchestrator

- **Read the full 578-line chapter directly.** Content matches the S60 and
  S62 session logs and PROGRESS.md KEY FINDING blocks closely on every
  numbered claim checked (the theorem's three forms, the rate corollary's
  exponent/constant recovery, the fractional counterexample's mechanism,
  the Zhang-construction derivation and its `a=0`/`b=1/2` exponents, the
  three-row candidate table). No claim was found that contradicts or
  overstates the underlying session logs.
- **`\ref` discipline, checked by grep, not by trusting the subagent's own
  claim of having checked it.** Extracted every `\ref{...}` in the new
  chapter (10 occurrences, 7 distinct targets) and confirmed each target's
  `\label` exists in a book chapter file: `ch:pressure-testing` and
  `sec:pt-s55` (in `17b-pressure-testing.tex`), `ch:openai` (in
  `17-openai-forced.tex`), `sec:wall9` (in `18-walls.tex`),
  `sec:ff-s60-audit`, `tab:ff-candidates`, `tab:ff-rate` (all in the new
  chapter itself). **Zero illegal references to proof-document labels** —
  separately grepped for every `thm:cd1-*`, `cor:cd1-*`, `lem:wpt*`,
  `prop:wpt*`, `prop:cd1-*`, `rem:cd1-*`, `rem:wpt*` pattern used with
  `\ref{}` specifically: no matches, confirming every proof-document object
  is cited as `\texttt{...}` text throughout, as required.
- **Diff to `book/tfirst-book.tex`**: read directly, confirmed to be exactly
  the single claimed `\include` line, correctly placed.
- **The two ch.18/ch.22 patches**: read both diffs directly. Both are
  purely additive paragraphs, both correctly use `\Sn{63}` and
  `\ref{ch:forcing-floor}` (confirmed resolvable above), both stop short of
  touching the existing enumerated content they follow.
- **Build claims, independently reproduced from scratch, not taken on
  trust.** Stashed the working-tree changes (`git stash push -u`), built
  the unmodified book fresh (3 `pdflatex` passes): **233 pages, 0 errors,
  0 LaTeX warnings, 355 total overfull+underfull hboxes (37+318, matching
  the subagent's reported baseline split exactly), 0 bookmark "Token not
  allowed" warnings, exactly 1 pre-existing `pdfTeX warning (dest):
  name{Hfootnote.12}...`.** Restored the changes (`git stash pop`) and
  rebuilt (3 passes): **243 pages (+10), 0 errors, 0 LaTeX warnings, 355
  total overfull+underfull (identical to baseline, not merely close), 0
  bookmark warnings, the identical single pre-existing `pdfTeX warning
  (dest)`.** Every number in the subagent's report reproduced exactly by
  independent rebuild.

No error was found in the subagent's work. No changes were requested before
merge.

## Merge

All four files (`book/chapters/17c-forcing-floor.tex` new,
`book/chapters/18-walls.tex`, `book/chapters/22-where-it-stands.tex`,
`book/tfirst-book.tex` modified) staged and committed as-is, alongside this
session log and the corresponding `PROGRESS.md` update.

## Standing discipline confirmed

- Session tag pre-assigned S63 (concurrent with S64); verified at session
  start against tip `5e436e9`, no collision with S64's independent S64 tag.
- `proofs/`, `layer3/`, `layer4/`, `results/` untouched — confirmed via
  `git status` before merge (only the four `book/` files plus this session's
  own new files were ever modified).
