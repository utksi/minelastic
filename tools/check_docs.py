# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Checks of the documentation text: typography, spelling, prose style and the separation of the tutorials.

1. ASCII typography in docs/, README.md, examples/, src/, tests/ and tools/: no typographic dashes, minus
   signs, ellipses, arrows, comparison or multiplication signs, curly quotes or Unicode super- and
   subscripts. Letters with diacritics and Greek letters (Å, ö, α, β) are allowed.
2. The tutorials (docs/get-started/) do not link into the theory pages, so that they can be followed without
   them.
3. Spelling with codespell, with the intentional names of IGNORE.
4. Prose style with proselint on the prose only (code blocks, inline code, maths, links and HTML comments
   removed), without its rule typography.symbols, which asks for the symbols that rule 1 forbids.

    python tools/check_docs.py          # prints every finding; exit status 1 if there is any
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEXT_DIRS = ["docs", "examples", "src", "tests", "tools"]
FORBIDDEN = re.compile("[\u2010-\u2015\u2026\u2190-\u21ff\u2212\u2264\u2265\u00d7\u2018\u2019\u201c\u201d"
                       "\u00b2\u00b3\u00b9\u2070-\u209f]")
# the scheme name asess (all six single strains), the constant EQUIL, the CP2K keyword OT and VASP's spelling
# "cummulative" in the OUTCAR excerpts of the tests
IGNORE = ["asess", "equil", "ot", "cummulative"]
PROSELINT_OFF = ["typography.symbols"]


def files(exts):
    """Paths (relative to the repository) of the files with the given extensions in TEXT_DIRS and README.md."""
    out = ["README.md"]
    for d in TEXT_DIRS:
        for base, _, names in os.walk(os.path.join(ROOT, d)):
            out += [os.path.relpath(os.path.join(base, n), ROOT) for n in names if n.endswith(exts)]
    return sorted(out)


def typography():
    """Findings 'path:line: character' for every forbidden typographic character."""
    found = []
    for path in files((".md", ".py", ".yml", ".txt", ".json")):
        for i, line in enumerate(open(os.path.join(ROOT, path), encoding="utf-8"), 1):
            for ch in FORBIDDEN.findall(line):
                found.append(f"{path}:{i}: typography: U+{ord(ch):04X} {ch!r}")
    return found


def tutorial_links():
    """Findings for links from docs/get-started/ into docs/theory/."""
    found = []
    for path in files((".md",)):
        if not path.startswith(os.path.join("docs", "get-started")):
            continue
        for i, line in enumerate(open(os.path.join(ROOT, path)), 1):
            if re.search(r"\]\([^)]*theory/", line):
                found.append(f"{path}:{i}: tutorial links into the theory pages")
    return found


def spelling():
    """codespell findings in all text files, except the names in IGNORE."""
    paths = [os.path.join(ROOT, p) for p in files((".md", ".py"))]
    r = subprocess.run(["codespell", "--ignore-words-list", ",".join(IGNORE)] + paths,
                       capture_output=True, text=True)
    return [f"{os.path.relpath(ln.split(':')[0], ROOT)}:{ln.split(':', 1)[1]}" for ln in r.stdout.splitlines()]


def prose(text):
    """The prose of a markdown page: code blocks, inline code, maths, link targets, attribute lists, admonition
    markers and HTML comments are replaced by blanks of the same line structure, so that line numbers are kept."""
    keep_lines = lambda m: "\n" * m.group(0).count("\n")                  # noqa: E731
    text = re.sub(r"(?ms)^(\s*)(```|~~~).*?^\1\2\s*$", keep_lines, text)
    text = re.sub(r"(?s)<!--.*?-->", keep_lines, text)
    text = re.sub(r"(?s)\$\$.*?\$\$", keep_lines, text)
    text = re.sub(r"\$[^$\n]+\$", "x", text)
    text = re.sub(r"`[^`\n]+`", "x", text)
    text = re.sub(r"\]\([^)\n]*\)", "]", text)
    text = re.sub(r"(?m)^:::.*$", "", text)
    text = re.sub(r"\{:?\s*[.#][^}\n]*\}", "", text)                         # attribute lists { .class }
    text = re.sub(r"(?m)^\s*(!!!|\?\?\?\+?) .*$", "", text)                     # admonition markers
    return text


def style():
    """proselint findings on the prose of every markdown page."""
    from proselint.config import default
    from proselint.tools import lint
    config = {**default, "checks": {**default["checks"], **{k: False for k in PROSELINT_OFF}}}
    found = []
    for path in files((".md",)):
        for check, message, line, column, *_ in lint(prose(open(os.path.join(ROOT, path)).read()), config=config):
            found.append(f"{path}:{line + 1}:{column}: {check}: {message}")
    return found


def main():
    """Runs all checks, prints the findings and returns 1 if there are any."""
    found = typography() + tutorial_links() + spelling() + style()
    print("\n".join(found) if found else "no findings")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
