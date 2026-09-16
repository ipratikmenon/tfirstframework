# Handover — a new book chapter: the (C)/(D) forcing-floor arc (S60-S62)

Written 2026-09-16. Owner's direction: launch this alongside a third
strategic audit (`HANDOVER-S64-AUDIT3.md`), in parallel, as two independent
sessions. **You are the book-chapter session. Pre-assigned tag: S63.**
Verify this yourself (`git log --oneline -1`, `ls SESSION-LOG/`) before
writing anything, but the tag is pre-allocated specifically so you and the
concurrent audit session do not collide — use S63, and if your own check
somehow shows S63 already taken, stop and report rather than guessing a
replacement.

## Scope: `book/` only. Do not touch `proofs/`, `layer3/`, `layer4/`, `results/`.

A second, independent session is concurrently running a fresh strategic
audit against the wider literature. Do not coordinate with it, do not read
its output, do not touch anything outside `book/`.

## What this chapter covers, and why it deserves one

Three sessions (S60-S62) form a second coherent story after the S53-S57 arc
already documented in `book/chapters/17b-pressure-testing.tex`, this time
about the programme's pivot from (A)/(B) to (C)/(D):

1. **S60 (CD-1)**: after the internal toolkit (S52/S53), the best external
   candidate (S57), and the programme's own founding thesis (S59) all
   independently confirmed Wall 1 is real, the owner pivoted the active
   research question from "prove (A)/(B)" to "how much forcing do existing
   forced blow-up constructions actually require, and is it removable?" —
   explicitly NOT a step toward (A)/(B), and explicitly never introducing
   forcing into this programme's own arguments (a standing, absolute rule).
   S60 generalised S55's one-profile forcing floor (`lem:wpt3-forced-max`)
   into a genuine class theorem (`thm:cd1-general-floor`, `cor:cd1-rate`),
   valid for the whole class of axisymmetric fields with a prescribed swirl
   growth rate and core-radius geometry, for every `ν≥0` including Euler —
   the first genuinely general theorem this line produced, not just a check
   against one external object. It also found the OpenAI unforced Euler
   construction structurally viscosity-fragile at every fixed `ν>0` (Thread
   2), while flagging that the general vanishing-viscosity question itself
   is a real, open question in the literature (Constantin's own "gap
   between the two ranges of viscosities") and explicitly NOT recommending
   it as a work package.
2. **S61**: a scoping handover only (no session log; `HANDOVER-S61-CD2.md`),
   pivoting Thread 1 forward: extend the swirl maximum principle to
   fractional dissipation, and find a second instance to test the class
   theorem against.
3. **S62 (CD-2)**: Deliverable 1 — the fractional extension **fails**, and
   precisely why is now known: the classical proof needs the vector-Laplacian
   identity `(Δu)_θ=Δu_θ-u_θ/r²` to close a scalar equation on `Γ=ru_θ`,
   which the fractional operator's genuinely nonlocal, sign-changing angular
   weight has no analogue for. A concrete counterexample makes this
   precise rather than asserted. Deliverable 2 — a second genuine instance
   was found (Qi Zhang, arXiv:2311.12306, an axisymmetric forced ASNS
   construction matching this document's own system exactly) but it sits in
   a degenerate corner of the theorem (growth exponent `a=0`), giving only
   the trivial floor. **Recommendation, taken: the (C)/(D) forcing-floor
   line is closed.**

**The chapter's spine, stated plainly as its own thesis, not left implicit**:
this arc is the point where the programme stopped asking "can we prove
(A)/(B)" and started asking a well-posed, answerable question about
existing forced constructions instead — and, having built one genuinely new
piece of general mathematics (the forcing-floor class theorem) and pushed it
to its natural limit (fractional dissipation, a second instance), found that
limit reached within three sessions. This is a complete, small, honest
research arc: a real theorem, a real negative result precisely located, and
a clean stopping point — not a breakthrough on the Prize, and the chapter
should say so as plainly as ch. 17b did for the S53-S57 arc.

## Where it goes in the book, and what to touch

**Write it as `book/chapters/17c-forcing-floor.tex`** (a `\chapter{...}` with
its own `\label{ch:forcing-floor}`), and add exactly one line to
`book/tfirst-book.tex`:
```
\include{chapters/17c-forcing-floor}
```
placed immediately after `\include{chapters/17b-pressure-testing}` (still
inside `\part{The Landscape}`) — **read `book/chapters/17b-pressure-testing.tex`
in full first** to match its voice, its citation discipline, and to avoid
duplicating any of its content (it already covers S53-S57 and the "Wall 1
survives every test" thesis; this chapter picks up exactly where it leaves
off, at the (A)/(B)-to-(C)/(D) pivot, and should reference it rather than
re-explain it). Also skim `book/chapters/17-openai-forced.tex` for voice
consistency, since both this chapter and 17b build on it.

**A second, small, explicitly bounded task**: `book/chapters/18-walls.tex`
and `book/chapters/22-where-it-stands.tex` may now be stale on the specific
point of the programme's *active research question* (they likely still
describe the programme as working toward (A)/(B) directly, without
reflecting the C/D pivot). Do **not** rewrite either chapter. Find the
specific stale entries and patch them with minimal, precise amendments in
the "Amended at S63: ..." style used at S48 and S58. If you are not
confident a given passage is actually stale, leave it alone rather than
guessing.

## Sources to read before writing a word

`SESSION-LOG/2026-09-13-S60-cd1-forcing-floor-general.md`,
`HANDOVER-S61-CD2.md`, `SESSION-LOG/2026-09-16-S62-cd2-fractional-and-second-instance.md`
in full. Cross-check every claim against the PROGRESS.md KEY FINDING blocks
for S60 and S62. Where finer technical detail is needed, go to the live
proof-document section `sec:cd1-general-floor` (S60; note S62 made no
proofs/ edit, so there is nothing new there to read for S62's content —
its findings live only in the session log and PROGRESS.md). For the
fractional-extension counterexample and the Zhang-paper self-similar
reduction specifically, the session log's own derivations are the primary
record; if you want the actual external sources for extra confidence, the
relevant arXiv IDs are 2309.08495, 2407.06776, and 2311.12306 (all already
fetched and quoted in the S62 session log — re-fetching is optional, not
required, since this chapter documents this programme's own work, not the
external papers' content directly).

## Style contract — identical to every existing chapter

- **This is a different LaTeX document from `proofs/claim_a_3d_proof_attempt.tex`.**
  Cite proof-document objects as **text**, e.g. `\texttt{thm:cd1-general-floor}`
  — never as `\ref{...}`. Only `\ref` a label that exists in a book chapter file.
- **Sources, always.** Every quantitative claim traces to a session log or a
  labelled object in the proof document, given as text per above.
- **The proof document and session logs govern**, where sources disagree.
- **No narrative drift into self-congratulation or advocacy for (C)/(D).**
  This chapter documents a real theorem, a real negative result, and an
  honest stopping point. It does not campaign for (C)/(D) as a goal — that
  remains this programme's standing prohibition (H4/H5 in S60-S62's own
  language), and the chapter should make clear the programme studies other
  people's forcing, never its own.
- Use the master's macros (`\Sn{}`, `\Lfac`, `\sstar`, `\rstar`, `\ehat`,
  `\norm{}`, `\abs{}`, the status macros, the `record` environment,
  `statement`/`defn`/`obsv`). Do not define new macros.
- **Guard math in section headings** with `\texorpdfstring{$x$}{x}` — check
  the current baseline yourself (`grep -c "Token not allowed" <log>`) before
  and after, do not assume a number.
- **Diagrams and tables expected**, matching the rest of the book's density.
  At minimum: one figure or table showing the three-session arc (pivot →
  general theorem → fractional attempt/second-instance search → closed),
  and a small table for the two candidate papers checked in S62 (ruled out
  on axisymmetry) plus the one instance found (Zhang, degenerate).
- Write for a mathematician who has read the rest of the book but not this
  thread. Prose, not bullet fragments.

## Build verification — required before reporting done

```
mkdir -p /tmp/bookbuild-s63/chapters
cd /home/user/tfirstframework/book
pdflatex -interaction=nonstopmode -file-line-error -output-directory=/tmp/bookbuild-s63 tfirst-book.tex
```
Three passes. Target: **0 errors, 0 LaTeX warnings**, bookmark warnings at or
below whatever the current baseline is — measure it yourself on the
unmodified book first, then again after your changes.

Do not commit or push. Report: the chapter's line count, the master-file
diff, the final page count, the before/after warning counts, exactly which
(if any) passages in ch. 18/22 you patched and the precise diff for each.

## Standing discipline

- Session tag pre-assigned as S63 (see above) — verify, don't assume, but
  this is the coordination mechanism for this specific parallel launch.
- Fetch/re-read primary programme sources (session logs, PROGRESS.md) directly;
  do not recall this handover's own summary as a substitute for reading them.
- Do not commit or push. Report your findings and the build results directly.
