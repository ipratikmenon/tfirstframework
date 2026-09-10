#!/usr/bin/env python3
"""
check_proof.py — mechanical structural checker for the T-First proof documents.

Motivation (HANDOVER-S46-OPUS.md, WP0 / W7): six sessions (S31, S34, S35, S39,
S41, S44) shipped consequential errors into a LaTeX document that had never been
compiled.  This script catches the *syntactic* class of those errors without
needing a TeX toolchain, so it can run in any container and as a pre-commit gate.

Checks implemented
------------------
  1. brace balance            — comment-stripped, with \\{  \\}  \\\\  \\%  handled
  2. \\begin/\\end multiset    — equal counts per environment name, and correct
                                nesting order (stack discipline)
  3. undeclared environments  — every \\begin{X} must be a LaTeX/package builtin
                                or declared via \\newtheorem / \\newenvironment
                                in the preamble
  4. duplicate labels         — no \\label{k} may appear twice
  5. unresolved references    — every \\ref/\\eqref/\\autoref/\\Cref/\\cref target
                                must have a \\label; every \\cite key must have a
                                \\bibitem (or a .bib entry, if one is present)
  6. malformed-end typo class — the `end{...>` / `\\end{...>` corruption found in
                                S44, plus `\\begin{...>`

Not covered (run pdflatex for these)
------------------------------------
Undefined control sequences.  S47 found three (\\colonequals at line 4584,
\\fint at 5254 and 5912) that no purely structural check could have seen; only
a compile catches that class.  `check_proof.py` green is necessary, not
sufficient.

Exit status: 0 iff every check passes; 1 otherwise.

Usage
-----
    python scripts/check_proof.py                       # default document
    python scripts/check_proof.py path/to/file.tex ...  # explicit targets
    python scripts/check_proof.py --quiet               # errors only
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import Counter, defaultdict

DEFAULT_TARGETS = ["proofs/claim_a_3d_proof_attempt.tex"]

# Environments LaTeX / the loaded packages provide.  Anything not here and not
# declared in the document itself is flagged.
BUILTIN_ENVIRONMENTS = {
    # core
    "document", "abstract", "titlepage", "figure", "figure*", "table", "table*",
    "center", "flushleft", "flushright", "quote", "quotation", "verse",
    "verbatim", "verbatim*", "minipage", "tabular", "tabular*", "tabbing",
    "itemize", "enumerate", "description", "list", "trivlist", "picture",
    "thebibliography", "theindex", "footnotesize", "small", "sloppypar",
    "filecontents", "filecontents*",
    # amsmath / mathtools
    "equation", "equation*", "align", "align*", "alignat", "alignat*",
    "aligned", "alignedat", "gather", "gather*", "gathered", "multline",
    "multline*", "split", "flalign", "flalign*", "eqnarray", "eqnarray*",
    "array", "cases", "dcases", "rcases", "drcases", "subequations",
    "matrix", "pmatrix", "bmatrix", "Bmatrix", "vmatrix", "Vmatrix",
    "smallmatrix", "psmallmatrix", "bsmallmatrix", "displaymath", "math",
    "proof",
    # enumitem
    "enumerate*", "itemize*", "description*",
    # tikz / pgf, harmless to allow
    "tikzpicture", "scope", "axis",
}

# Declaration forms that create an environment name.
DECL_PATTERNS = [
    re.compile(r"\\newtheorem\*?\s*\{([^}]*)\}"),
    re.compile(r"\\newenvironment\*?\s*\{([^}]*)\}"),
    re.compile(r"\\renewenvironment\*?\s*\{([^}]*)\}"),
    re.compile(r"\\DeclareDocumentEnvironment\s*\{([^}]*)\}"),
    re.compile(r"\\NewDocumentEnvironment\s*\{([^}]*)\}"),
]

BEGIN_RE = re.compile(r"\\begin\s*\{([^}]*)\}")
END_RE = re.compile(r"\\end\s*\{([^}]*)\}")
LABEL_RE = re.compile(r"\\label\s*\{([^}]*)\}")
REF_RE = re.compile(r"\\(?:eqref|autoref|ref|Cref|cref|nameref|pageref)\s*\{([^}]*)\}")
CITE_RE = re.compile(
    r"\\(?:cite|citep|citet|citeauthor|citeyear)\s*(?:\[[^\]]*\])*\s*\{([^}]*)\}"
)
BIBITEM_RE = re.compile(r"\\bibitem\s*(?:\[[^\]]*\])?\s*\{([^}]*)\}")
SECTION_RE = re.compile(r"\\(?:sub)*section\*?\s*\{")
# The S44 corruption class: an environment delimiter terminated by `>` instead
# of `}`, e.g.  \end{proof>  or  end{lemma>
MALFORMED_RE = re.compile(r"\\?(?:begin|end)\s*\{[A-Za-z@*]*>")

# Reference targets that are deliberately recorded as dangling.  A dangling
# reference must be *declared* here (and rendered in red in the document) rather
# than silently guessed at — see HANDOVER-S46-OPUS.md, WP0 named failure mode.
DANGLING_MARKER_RE = re.compile(r"%\s*CHECK-PROOF:\s*ALLOW-DANGLING\s+(\S+)")


class Findings:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.notes: list[str] = []

    def error(self, msg: str) -> None:
        self.errors.append(msg)

    def note(self, msg: str) -> None:
        self.notes.append(msg)


def strip_comments(line: str) -> str:
    r"""Return the line with its TeX comment removed, honouring \% and \\."""
    out = []
    i = 0
    n = len(line)
    while i < n:
        c = line[i]
        if c == "\\":
            out.append(line[i:i + 2])
            i += 2
            continue
        if c == "%":
            break
        out.append(c)
        i += 1
    return "".join(out)


def scan_braces(lines: list[str], findings: Findings, path: str) -> None:
    """Brace balance with escapes and comments handled; reports section depths."""
    depth = 0
    section_depths: list[tuple[int, int, str]] = []
    first_negative = None
    for lineno, raw in enumerate(lines, start=1):
        line = strip_comments(raw)
        if SECTION_RE.search(line):
            section_depths.append((lineno, depth, raw.strip()[:70]))
        i = 0
        n = len(line)
        while i < n:
            c = line[i]
            if c == "\\":
                i += 2
                continue
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth < 0 and first_negative is None:
                    first_negative = lineno
            i += 1
    if first_negative is not None:
        findings.error(
            f"{path}: brace depth went NEGATIVE (unmatched '}}') first at "
            f"line {first_negative}"
        )
    if depth != 0:
        findings.error(f"{path}: brace imbalance — final depth {depth:+d} (expected 0)")
        bad = [s for s in section_depths if s[1] != 0]
        if bad:
            findings.error(
                f"{path}: first section heading reached at nonzero depth: "
                f"line {bad[0][0]} (depth {bad[0][1]:+d}) — {bad[0][2]}"
            )
    else:
        findings.note(f"{path}: brace balance OK (final depth 0)")


def scan_environments(lines: list[str], findings: Findings, path: str) -> None:
    declared: set[str] = set()
    for raw in lines:
        line = strip_comments(raw)
        for pat in DECL_PATTERNS:
            for m in pat.finditer(line):
                declared.add(m.group(1).strip())

    begins: Counter = Counter()
    ends: Counter = Counter()
    stack: list[tuple[str, int]] = []
    order_errors = 0

    for lineno, raw in enumerate(lines, start=1):
        line = strip_comments(raw)
        events: list[tuple[int, str, str]] = []
        for m in BEGIN_RE.finditer(line):
            events.append((m.start(), "b", m.group(1).strip()))
        for m in END_RE.finditer(line):
            events.append((m.start(), "e", m.group(1).strip()))
        events.sort()
        for _, kind, name in events:
            if kind == "b":
                begins[name] += 1
                stack.append((name, lineno))
            else:
                ends[name] += 1
                if not stack:
                    findings.error(
                        f"{path}:{lineno}: \\end{{{name}}} with no open environment"
                    )
                    order_errors += 1
                elif stack[-1][0] != name:
                    open_name, open_line = stack[-1]
                    findings.error(
                        f"{path}:{lineno}: \\end{{{name}}} closes "
                        f"\\begin{{{open_name}}} opened at line {open_line} "
                        f"(nesting mismatch)"
                    )
                    order_errors += 1
                    stack.pop()
                else:
                    stack.pop()

    for name, lineno in stack:
        findings.error(f"{path}:{lineno}: \\begin{{{name}}} is never closed")

    names = set(begins) | set(ends)
    mismatched = [n for n in sorted(names) if begins[n] != ends[n]]
    for n in mismatched:
        findings.error(
            f"{path}: environment '{n}' count mismatch — "
            f"{begins[n]} \\begin vs {ends[n]} \\end"
        )
    if not mismatched and not order_errors and not stack:
        findings.note(
            f"{path}: \\begin/\\end multiset OK ({sum(begins.values())}/"
            f"{sum(ends.values())} balanced, {len(names)} distinct environments)"
        )

    undeclared = sorted(
        n for n in names if n and n not in BUILTIN_ENVIRONMENTS and n not in declared
    )
    if undeclared:
        first_use: dict[str, int] = {}
        for lineno, raw in enumerate(lines, start=1):
            for m in BEGIN_RE.finditer(strip_comments(raw)):
                nm = m.group(1).strip()
                if nm in undeclared and nm not in first_use:
                    first_use[nm] = lineno
        for n in undeclared:
            findings.error(
                f"{path}: environment '{n}' used but never declared "
                f"(first use line {first_use.get(n, '?')}) — add a "
                f"\\newtheorem/\\newenvironment in the preamble"
            )
    else:
        findings.note(
            f"{path}: all environments declared ({len(declared)} local declarations)"
        )


def scan_labels_and_refs(lines: list[str], findings: Findings, path: str,
                         bib_keys: set) -> None:
    labels: dict[str, list[int]] = defaultdict(list)
    refs: dict[str, list[int]] = defaultdict(list)
    cites: dict[str, list[int]] = defaultdict(list)
    bibitems: set = set()
    allowed_dangling: set = set()

    for lineno, raw in enumerate(lines, start=1):
        for m in DANGLING_MARKER_RE.finditer(raw):
            allowed_dangling.add(m.group(1).strip())
        line = strip_comments(raw)
        for m in LABEL_RE.finditer(line):
            labels[m.group(1).strip()].append(lineno)
        for m in REF_RE.finditer(line):
            for key in m.group(1).split(","):
                key = key.strip()
                if key:
                    refs[key].append(lineno)
        for m in CITE_RE.finditer(line):
            for key in m.group(1).split(","):
                key = key.strip()
                if key:
                    cites[key].append(lineno)
        for m in BIBITEM_RE.finditer(line):
            bibitems.add(m.group(1).strip())

    dupes = {k: v for k, v in labels.items() if len(v) > 1}
    for k in sorted(dupes):
        findings.error(
            f"{path}: duplicate \\label{{{k}}} at lines "
            + ", ".join(str(x) for x in dupes[k])
        )
    if not dupes:
        findings.note(f"{path}: {len(labels)} labels, no duplicates")

    unresolved = sorted(k for k in refs if k not in labels)
    flagged, hard = [], []
    for k in unresolved:
        (flagged if k in allowed_dangling else hard).append(k)
    for k in hard:
        findings.error(
            f"{path}: unresolved reference '{k}' at line(s) "
            + ", ".join(str(x) for x in refs[k][:5])
            + " — no matching \\label"
        )
    for k in flagged:
        findings.note(
            f"{path}: reference '{k}' is a DECLARED DANGLING reference "
            f"(line(s) {', '.join(str(x) for x in refs[k][:5])}) — intentional, "
            f"rendered in red in the document"
        )
    if not hard:
        findings.note(
            f"{path}: {len(refs)} distinct \\ref targets, all resolved"
            + (f" ({len(flagged)} declared dangling)" if flagged else "")
        )

    known_keys = bibitems | bib_keys
    if known_keys:
        bad_cites = sorted(k for k in cites if k not in known_keys)
        for k in bad_cites:
            findings.error(
                f"{path}: unresolved \\cite{{{k}}} at line(s) "
                + ", ".join(str(x) for x in cites[k][:5])
            )
        if not bad_cites:
            findings.note(f"{path}: {len(cites)} distinct \\cite keys, all resolved")
    elif cites:
        findings.error(
            f"{path}: {len(cites)} \\cite keys but no \\bibitem or .bib entries found"
        )


def scan_malformed(lines: list[str], findings: Findings, path: str) -> None:
    hits = []
    for lineno, raw in enumerate(lines, start=1):
        line = strip_comments(raw)
        for m in MALFORMED_RE.finditer(line):
            hits.append((lineno, m.group(0)))
    for lineno, text in hits:
        findings.error(
            f"{path}:{lineno}: malformed environment delimiter '{text}' "
            f"(the S44 'end{{...>' typo class — '>' should be '}}')"
        )
    if not hits:
        findings.note(f"{path}: no malformed 'end{{...>' delimiters")


def collect_bib_keys(texpath: str) -> set:
    keys: set = set()
    directory = os.path.dirname(os.path.abspath(texpath)) or "."
    for root in {directory, os.getcwd()}:
        try:
            entries = os.listdir(root)
        except OSError:
            continue
        for name in entries:
            if name.endswith(".bib"):
                try:
                    with open(os.path.join(root, name), encoding="utf-8",
                              errors="replace") as fh:
                        for m in re.finditer(r"@\w+\s*\{\s*([^,\s]+)", fh.read()):
                            keys.add(m.group(1).strip())
                except OSError:
                    pass
    return keys


def check_file(path: str, findings: Findings) -> None:
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    findings.note(f"{path}: {len(lines)} lines")
    scan_braces(lines, findings, path)
    scan_environments(lines, findings, path)
    scan_labels_and_refs(lines, findings, path, collect_bib_keys(path))
    scan_malformed(lines, findings, path)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Structural checker for the T-First LaTeX proof documents.")
    ap.add_argument("targets", nargs="*", default=None,
                    help="LaTeX files to check (default: the Claim-A proof document)")
    ap.add_argument("--quiet", action="store_true", help="print failures only")
    args = ap.parse_args(argv)

    targets = args.targets or DEFAULT_TARGETS
    findings = Findings()
    for t in targets:
        if not os.path.isfile(t):
            findings.error(f"{t}: file not found")
        else:
            check_file(t, findings)

    if not args.quiet:
        for n in findings.notes:
            print(f"  ok   {n}")
    for e in findings.errors:
        print(f"FAIL   {e}", file=sys.stderr)

    if findings.errors:
        print(f"\ncheck_proof.py: {len(findings.errors)} FAILURE(S)", file=sys.stderr)
        return 1
    if not args.quiet:
        print(f"\ncheck_proof.py: all checks passed for {', '.join(targets)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
