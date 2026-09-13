# Handover — a new book chapter: "Putting the Landscape to the Test"

Written 2026-09-13. Owner's direction: close the S53-S57 thread by writing it
up properly as a new chapter of `book/`, documenting the S53-S57 arc as a
coherent whole.

## Scope: `book/` only. Do not touch `proofs/`, `layer3/`, `layer4/`, `results/`.

A second, independent session is concurrently running a fresh strategic
audit against `proofs/` and the wider literature. Do not coordinate with it,
do not read its output, do not touch anything outside `book/`.

## What this chapter covers, and why it deserves one

Five sessions (S53-S57) form a single coherent story that the existing book
(finished at S49, 27 chapters, 219 pages) does not contain at all:

1. **S53**: proved (`lem:wpt1-nueff`, `cor:wpt1-trichotomy`) that any
   rescaling normalising a supercritical blow-up rate forces the effective
   viscosity to zero, so the limit solves Euler — corrected the document's
   own prior assertion ("rescale to nothing") to something sharper and
   truer. Verdict on Type-II: open, unreachable by this programme's tools,
   belongs to Seregin's ancient-Euler Liouville programme.
2. **S54**: took the one concrete external Type-II candidate in existence
   (OpenAI's forced NS construction) and tested its profile against
   ancient-solution exclusion. Found — via a route shorter than the one
   planned, using S52's own pre-existing swirl maximum principle — that the
   profile is impossible in the unforced axisymmetric class. A narrow,
   verified positive result, explicitly not a statement about OpenAI's
   actual (forced, non-axisymmetric) theorem.
3. **S55**: turned Wall 9's qualitative "forcing is load-bearing" into a
   derived, source-confirmed quantitative floor on the forcing OpenAI's
   construction requires — two independent derivations (a forced maximum
   principle here; the source's own residual-integration formalism there)
   landing on the same exponent.
4. **S56**: a fresh strategic audit. Found the single most relevant external
   candidate available anywhere in the recent literature — Grujić's
   log-weighted BMO depletion mechanism (arXiv:2607.08866v3) — the only
   methodology found that plausibly satisfies Wall 1's own escape clause.
5. **S57**: priced that candidate rigorously and refuted it as a route
   around Wall 1, via three explicit counterexamples, while keeping one
   genuine byproduct (a verified commutator/Coifman-Rochberg-Weiss
   technique) as a small, real contribution independent of the paper's
   overall fate.

**The chapter's spine, stated plainly and to be the chapter's own explicit
thesis, not left implicit**: this arc is the fifth, sixth, and seventh
independent confirmation (after S48 and S50, S52) that Wall 1 is real and
does not move under relabeling, reformulation, or the best available
external candidate this programme could find. That is the "new revelation"
— not a breakthrough on the Prize, but the most thorough stress-test this
programme's central obstruction has ever received, and its survival of that
stress-test.

## Where it goes in the book, and what to touch

**Write it as `book/chapters/17b-pressure-testing.tex`** (a `\chapter{...}`
with its own `\label{ch:pressure-testing}`), and add exactly one line to
`book/tfirst-book.tex`:
```
\include{chapters/17b-pressure-testing}
```
placed immediately after `\include{chapters/17-openai-forced}`, still inside
`\part{The Landscape}` — the existing chapter 17 surveys the external
landscape; this chapter is the direct technical response to it, and belongs
in the same part. **Read `book/chapters/17-openai-forced.tex` in full
first** to match its voice, its citation discipline, and to avoid
duplicating any of its content (it already covers the OpenAI/Alpöge-
Buckmaster background and the S49 two-strands correction — this new chapter
picks up from there, at S53, and should reference it rather than re-explain
it).

**A second, small, explicitly bounded task**: `book/chapters/18-walls.tex`
and `book/chapters/22-where-it-stands.tex` were written before S50-S57
happened and are now stale on specific points. Do **not** rewrite either
chapter. Find the specific stale entries (Wall 1's confirmation count, any
statement about the axisymmetric pivot's status, any final ledger row
touching Wall 9 or the OpenAI construction) and patch them with minimal,
precise amendments in the style ch. 11's table used at S48 ("Amended at
S48: ...") — a small addition noting the update and pointing to the new
chapter, not a silent rewrite. If you are not confident a given passage is
actually stale, leave it alone rather than guessing.

## Sources to read before writing a word

All five session logs in full: `SESSION-LOG/2026-09-12-S53-typeII-
reconnaissance.md`, `2026-09-12-S54-wpt2-openai-profile.md`, `2026-09-12-
S55-wpt3-forcing-floor.md` (dated 2026-09-13 as filed — check the actual
filename), `2026-09-12-S56-strategic-audit-literature.md`,
`2026-09-13-S57-wpa-grujic-bmo-pricing.md`. Cross-check every claim in them
against the current PROGRESS.md KEY FINDING blocks for S53 through S57, and
where finer detail is needed, against the live proof document sections
directly: `sec:typeII-wpt1` (S53), the swirl-exclusion material added at
S54, `sec:wpt3-forced-floor`-region (S55, verify actual label), and
`sec:wpa-s57` (S57) in `proofs/claim_a_3d_proof_attempt.tex`. **Do not read
`proofs/` in order to write into it — you are reading it only as a source
for the book chapter, which is a separate document.**

## Style contract — identical to every existing chapter, restated because it matters most on citation discipline

- **This is a different LaTeX document from `proofs/claim_a_3d_proof_attempt.tex`.**
  Cite proof-document objects as **text**, e.g.
  `\texttt{lem:wpt1-nueff}, \S\,\texttt{sec:typeII-wpt1}` — never as
  `\ref{...}`, which would silently fail to resolve or, worse, coincidentally
  resolve to an unrelated label in the book's own namespace. Only `\ref` a
  label that exists in a book chapter file.
- **Sources, always.** Every quantitative claim traces to a session log, a
  `results.db` row (none apply here — this thread produced no numerics), or
  a labelled object in the proof document, with the object's name given as
  text per above.
- **The proof document governs**, where the two disagree.
- **No narrative drift into self-congratulation.** This chapter documents
  five real, correct, verified results that all failed to find a route
  around Wall 1. Say that plainly. The chapter's honest emotional register
  is "we tested this as hard as we know how, and the wall held" — not
  "look how clever the machinery is."
- Use the master's macros (`\Sn{57}`, `\Lfac`, `\sstar`, `\rstar`, `\ehat`,
  `\norm{}`, `\abs{}`, the status macros `\PROVED \OPEN \IMPOSSIBLE
  \WITHDRAWN \SUPERSEDED \PROVEDMOD{} \CONDITIONAL{}`, the `record`
  environment, `statement`/`defn`/`obsv`). Do not define new macros.
- **Guard math in section headings** with `\texorpdfstring{$x$}{x}` or
  hyperref emits bookmark warnings — check the current baseline yourself
  before and after (`grep -c "Token not allowed" <log>`), do not assume a
  number from a session report.
- **Diagrams and tables are expected**, matching the density of the rest of
  the book. At minimum: one figure showing the five-session arc as a chain
  (candidate found → tested → excluded/quantified/refuted, with what
  survived at each step), and the table comparing the three counterexamples
  from S57 (or citing/adapting ch. 17's existing comparison-table style).
- Write for a mathematician who has read the rest of the book but not this
  thread. Prose, not bullet fragments — this is a book, in the tradition of
  every other chapter.

## Build verification — required before reporting done

```
mkdir -p /tmp/bookbuild-s58/chapters
cd /home/user/tfirstframework/book
pdflatex -interaction=nonstopmode -file-line-error -output-directory=/tmp/bookbuild-s58 tfirst-book.tex
```
Three passes. Target: **0 errors, 0 LaTeX warnings**, bookmark warnings at or
below whatever the current baseline is (measure it yourself on the
unmodified book first, then again after your changes — do not assume a
number). `cm-super` must be installed or microtype aborts.

Do not commit or push. Report: the chapter's line count, the master-file
diff, the final page count, the before/after warning counts, exactly which
(if any) passages in ch. 18/22 you patched and the precise diff for each,
and any place a source disagreed with another and how you resolved it.
