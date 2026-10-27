# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Executable documentation: the code blocks of docs/ and README.md are run and their printed output compared.

Every fenced code block (except BibTeX) carries a marker, an HTML comment on the line before it:

    <!-- run -->                          python: executed, in one namespace per page
    <!-- output -->                       compared with the stdout of the run and shell blocks since the
                                          previous output block
    <!-- shell -->                        console or bash: the minelastic commands are executed; in a console
                                          block, the lines after "$ command" are its expected output
    <!-- output-of: examples/NAME.py -->  compared with the stdout of the example script
    <!-- excerpt: examples/NAME.py -->    every code line of the block occurs in the script
    <!-- file: PATH -->                   compared with the file PATH of the repository
    <!-- input -->                        a fragment of an input file, shown for reading only
    <!-- skip: REASON -->                 not executed (needs a DFT code, a cluster, MACE, ...)

The blocks of a page run in order in one temporary directory that holds copies of the example templates
(vasp_aln with a placeholder POTCAR, qe_aln, abinit_aln). Shell lines that run a DFT code (mpirun, srun,
sbatch) are not executed: `minelastic collect` then reads the stored stresses, forces and dipoles of the real
runs of these templates (tests/data_ionic, tests/data_docs) in place of the output files, after checking that
the generated strain states are those of the stored runs.

Output comparison: whitespace is ignored, a line "..." matches any number of lines and "..." within a line
matches any text. A number matches if it differs from the printed one by at most one unit in its last printed
digit (integers exactly). Example scripts that take minutes (SLOW) run only with MINELASTIC_DOCS_SLOW=1.
"""
import contextlib
import functools
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys

import numpy as np
import pytest

from minelastic import cli
from minelastic.io import abinit, qe, vasp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGES = sorted(os.path.relpath(os.path.join(d, f), ROOT) for d, _, fs in os.walk(os.path.join(ROOT, "docs"))
               for f in fs if f.endswith(".md")) + ["README.md"]
SLOW = {"examples/emt_finite_temperature.py", "examples/bks_quartz_piezo.py"}
DFT_RUN = re.compile(r"\b(mpirun|srun|sbatch|vasp_std|pw\.x|abinit|aims\.x)\b")
FENCE = re.compile(r"^(\s*)(```+|~~~+)\s*([\w+-]*)\s*$")
MARKER = re.compile(r"^\s*<!--\s*(output-of|output|run|shell|excerpt|file|input|skip)\s*:?\s*(.*?)\s*-->\s*$")
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def blocks(page):
    """Fenced code blocks of a markdown file: (line number, language, marker, argument, dedented body lines).

    The marker is the HTML comment on the nearest non-blank line above the fence (None if there is none).
    """
    lines = open(os.path.join(ROOT, page)).read().splitlines()
    out, i = [], 0
    while i < len(lines):
        m = FENCE.match(lines[i])
        if not m:
            i += 1
            continue
        indent, fence, lang = m.groups()
        j = i + 1
        while j < len(lines) and lines[j].strip() != fence:
            j += 1
        body = [ln[len(indent):] if ln.startswith(indent) else ln.strip() for ln in lines[i + 1:j]]
        k = i - 1
        while k >= 0 and not lines[k].strip():
            k -= 1
        mk = MARKER.match(lines[k]) if k >= 0 else None
        out.append((i + 1, lang, mk.group(1) if mk else None, mk.group(2) if mk else "", body))
        i = j + 1
    return out


def _unit(text):
    """One unit in the last printed digit of a number string; 0 for integers (exact match)."""
    mant, _, exp = text.lower().partition("e")
    if "." not in mant:
        return 0.0 if not exp else 10.0 ** int(exp)
    return 10.0 ** (-len(mant.split(".")[1]) + (int(exp) if exp else 0))


def _line_ok(expected, actual):
    """True if actual equals expected up to whitespace runs and to one unit in the last digit of each number;
    "..." inside a line matches any text."""
    nums, pieces = [], []
    for k, part in enumerate(expected.strip().split("...")):
        if k:
            pieces.append(".*?")
        texts = NUMBER.split(part)
        nums += NUMBER.findall(part)
        words = [r"\s*".join(map(re.escape, t.split())) for t in texts]
        pieces.append(r"\s*" + (r"\s*(" + NUMBER.pattern + r")\s*").join(words) + r"\s*")
    m = re.fullmatch("".join(pieces), " ".join(actual.split()))
    if not m:
        return False
    return all(abs(float(a) - float(e)) <= _unit(e) * (1 + 1e-9) for a, e in zip(m.groups(), nums))


def compare(expected, actual):
    """Expected lines (with "..." wildcards) against the actual output; leading and trailing blank lines are
    ignored. Returns None on a match, else a message with both texts."""
    def trim(lines):
        """Lines without trailing whitespace and without leading and trailing blank lines."""
        lines = [ln.rstrip() for ln in lines]
        while lines and not lines[0]:
            lines.pop(0)
        while lines and not lines[-1]:
            lines.pop()
        return lines
    exp, act = trim(expected), trim(actual.splitlines())

    @functools.lru_cache(maxsize=None)
    def match(i, j):
        """True if exp[i:] matches act[j:]; a "..." line consumes any number of actual lines."""
        if i == len(exp):
            return j == len(act)
        if exp[i].strip() == "...":
            return any(match(i + 1, k) for k in range(j, len(act) + 1))
        return j < len(act) and _line_ok(exp[i], act[j]) and match(i + 1, j + 1)
    if match(0, 0):
        return None
    return "expected:\n" + "\n".join(exp) + "\n\nactual:\n" + "\n".join(act)


@functools.lru_cache(maxsize=None)
def script_output(script):
    """stdout of an example script, run with the current interpreter in a temporary copy of examples/."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        env = {**os.environ, "PYTHONPATH": os.path.join(ROOT, "src") + os.pathsep + os.environ.get("PYTHONPATH", "")}
        r = subprocess.run([sys.executable, "-W", "ignore", os.path.join(ROOT, script)], cwd=tmp, env=env,
                           capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, f"{script} failed:\n{r.stderr}"
    return r.stdout


@functools.lru_cache(maxsize=None)
def stored_runs():
    """Stored results of the real runs of the example templates: {code: record with setup, stress, forces,
    dipoles}. VASP: the 12 single points of the --ionic run with LCALCPOL (the 4 strained cells are those of
    the clamped-ion run); QE and ABINIT: the 4 strained cells of the minimal set."""
    with open(os.path.join(ROOT, "tests", "data_ionic", "aln_run_ionic.json")) as f:
        runs = {"vasp": json.load(f)}
    with open(os.path.join(ROOT, "tests", "data_docs", "aln_dft_runs.json")) as f:
        runs.update({k: v for k, v in json.load(f).items() if k in ("qe", "abinit")})
    return runs


def _stored(run_dir, code):
    """Stored record of the state in run_dir, after checking that the generated state (strain, displaced
    atom and vector) and the reference structure equal those of the stored run."""
    sid = os.path.basename(os.path.normpath(run_dir))
    rec = stored_runs()[code]
    with open(os.path.join(os.path.dirname(os.path.normpath(run_dir)), "minelastic.json")) as f:
        setup = json.load(f)
    gen = {st["id"]: st for st in setup["states"]}[sid]
    ref = {st["id"]: st for st in rec["setup"]["states"]}
    assert sid in ref, f"no stored {code} run for state {sid}"
    for key in ("strain", "vector"):
        if key in gen:
            assert np.allclose(gen[key], ref[sid][key], atol=1e-9), f"{sid}: {key} differs from the stored run"
    a, b = setup["input_structure"], rec["setup"]["input_structure"]
    assert a["symbols"] == b["symbols"] and np.allclose(a["cell"], b["cell"], atol=1e-6)
    assert np.allclose(a["scaled_positions"], b["scaled_positions"], atol=1e-6)
    return sid, rec


def _replay_readers(monkeypatch):
    """Code readers that return the stored stress (GPa), forces (eV/A) and dipoles (e A) of the real runs."""
    def reader(code):
        """read_result of one code that returns the stored stress and forces."""
        def read_result(run_dir, output=None):
            """(stress (6,) in GPa, forces (N, 3) in eV/A) of the state in run_dir, from the stored run."""
            sid, rec = _stored(run_dir, code)
            return np.array(rec["stress_GPa"][sid]), np.array(rec["forces_eV_A"][sid])
        return read_result

    def read_dipole(run_dir, output="OUTCAR"):
        """Stored Berry-phase dipole (3,) in e A of the VASP state in run_dir."""
        sid, rec = _stored(run_dir, "vasp")
        return np.array(rec["dipole_eA"][sid])
    monkeypatch.setattr(vasp, "read_result", reader("vasp"))
    monkeypatch.setattr(vasp, "read_dipole", read_dipole)
    monkeypatch.setattr(qe, "read_result", reader("qe"))
    monkeypatch.setattr(abinit, "read_result", reader("abinit"))


def _workspace(tmp):
    """Copies of the example templates in tmp, with a placeholder POTCAR (ENMAX = 400 eV) for vasp_aln, and
    ibrion6/: the vasp_aln structure with the INCAR of a VASP IBRION = 6 run at VASP's defaults."""
    for name in ("vasp_aln", "qe_aln", "abinit_aln"):
        shutil.copytree(os.path.join(ROOT, "examples", name), os.path.join(tmp, name))
    with open(os.path.join(tmp, "vasp_aln", "POTCAR"), "w") as f:
        f.write("placeholder for the docs test\n   ENMAX  =  400.000; ENMIN  =  300.000 eV\n")
    shutil.copytree(os.path.join(tmp, "vasp_aln"), os.path.join(tmp, "ibrion6"))
    with open(os.path.join(tmp, "ibrion6", "INCAR"), "w") as f:
        f.write("IBRION = 6\nISIF = 3\nNFREE = 2\nPOTIM = 0.015\n")


def run_command(line):
    """Output of one shell line: minelastic commands run in-process, cat and ls on the workspace; DFT runs are
    replaced by the stored results (empty output). Anything else is an error in the docs test."""
    argv = shlex.split(line, comments=True)
    if not argv:
        return ""
    if argv[0] not in ("minelastic", "cat", "ls") and DFT_RUN.search(line):
        return ""
    if argv[0] == "cat":
        return "".join(open(f).read() for f in argv[1:])
    if argv[0] == "ls":
        return "\n".join(sorted(os.listdir(argv[1] if len(argv) > 1 else "."))) + "\n"
    assert argv[0] == "minelastic", f"the docs test runs only minelastic, cat and ls, not: {line}"
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            cli.main(argv[1:])
        except SystemExit as e:
            assert not e.code, f"{line} exited with {e.code}:\n{buf.getvalue()}"
    return buf.getvalue()


def shell_commands(body, console):
    """(command, expected output lines) of a console block ("$ " prompts) or of a bash block (one command per
    logical line, no expected output). Backslash continuations are joined."""
    out, cur = [], None
    for ln in body:
        if console and not ln.startswith("$ "):
            if cur is not None:
                cur[1].append(ln)
            continue
        text = ln[2:] if console else ln
        if cur is not None and cur[0].endswith("\\"):
            cur[0] = cur[0][:-1] + " " + text.strip()
            continue
        cur = [text, []]
        out.append(cur)
    return out


@pytest.mark.parametrize("page", PAGES)
def test_page(page, tmp_path, monkeypatch):
    """Runs every block of the page in order and compares the printed output with the page."""
    monkeypatch.chdir(tmp_path)
    _workspace(str(tmp_path))
    _replay_readers(monkeypatch)
    ns, buf = {"__name__": "__docs__"}, io.StringIO()
    for lineno, lang, marker, arg, body in blocks(page):
        where = f"{page}:{lineno}"
        if lang == "bibtex":
            continue
        assert marker, f"{where}: {lang or 'plain'} block without a docs-test marker"
        if marker == "skip":
            assert arg, f"{where}: skip needs a reason"
        elif marker == "run":
            assert lang == "python", where
            with contextlib.redirect_stdout(buf):
                exec(compile("\n".join(body), where, "exec"), ns)
        elif marker == "output":
            msg = compare(body, buf.getvalue())
            assert msg is None, f"{where}: output differs\n{msg}"
            buf = io.StringIO()
        elif marker == "shell":
            assert lang in ("console", "bash"), where
            for cmd, expected in shell_commands(body, lang == "console"):
                text = run_command(cmd)
                if lang == "console":
                    msg = compare(expected, text)
                    assert msg is None, f"{where}: output of `{cmd}` differs\n{msg}"
                else:
                    buf.write(text)
        elif marker == "output-of":
            if arg in SLOW and not os.environ.get("MINELASTIC_DOCS_SLOW"):
                continue
            msg = compare(body, script_output(arg))
            assert msg is None, f"{where}: output of {arg} differs\n{msg}"
        elif marker == "excerpt":
            src = {" ".join(ln.split()) for ln in open(os.path.join(ROOT, arg)).read().splitlines()}
            for ln in body:
                code = " ".join(ln.split("  #")[0].split())
                if code and code not in ("...",) and not code.startswith("#"):
                    assert any(code in s for s in src), f"{where}: line not in {arg}: {ln}"
        elif marker == "file":
            with open(os.path.join(ROOT, arg)) as f:
                msg = compare(body, f.read())
            assert msg is None, f"{where}: differs from {arg}\n{msg}"
