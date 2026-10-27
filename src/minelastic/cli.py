# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Command line interface: minelastic info | generate | collect | vasp-tensor.

info prints the symmetry, tensor form and strain set of a structure; generate and collect call
workflow.generate and workflow.collect; vasp-tensor converts a tensor printed by VASP (IBRION = 6) to GPa in
Voigt order and adds the Born stability check and the Voigt-Reuss-Hill moduli.
"""
import argparse
import json
import sys

import numpy as np

from . import __version__
from .core import ElasticSetup
from .properties import born_stable, vrh
from .workflow import CODES, collect, generate, read_structure


def _common(p):
    """Options of the strain set shared by info and generate."""
    p.add_argument("--scheme", default="minimal", choices=["minimal", "asess"],
                   help="minimal uncoupled set (default) or the six single-component strains")
    p.add_argument("--delta", type=float, default=0.01, help="strain step (default 0.01)")
    p.add_argument("--stencil", default=None, choices=["central", "four-point", "forward"],
                   help="finite differences: four-point (default; -2, -1, 1, 2 steps) or central (+-1 step; the "
                        "default with --md) or forward")
    p.add_argument("--symprec", type=float, default=1e-5, help="spglib tolerance in Angstrom")
    p.add_argument("--frame", default="auto", choices=["auto", "input", "standard"],
                   help="orientation in which the strains are applied (default: input if it is standard)")


def main(argv=None):
    """Parse the command line and run the subcommand; returns the exit status 0."""
    ap = argparse.ArgumentParser(
        prog="minelastic", epilog="Documentation: https://utksi.github.io/minelastic/",
        description="Elastic constants from the minimal uncoupled stress-strain sets. info shows the symmetry and "
                    "the strain set of a structure, generate writes the inputs of the strained cells, collect reads "
                    "the finished calculations and prints the elastic tensor, and vasp-tensor analyses the tensor "
                    "that VASP prints for IBRION = 6.")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("info", help="symmetry, tensor form and strain set of a structure")
    p.add_argument("structure")
    p.add_argument("--code", default="generic", choices=sorted(CODES),
                   help="how to read the structure (template dir for vasp/aims/abinit, pw.in for qe)")
    p.add_argument("--format", default=None, help="ASE format for --code generic")
    _common(p)

    p = sub.add_parser("generate", help="write the strained inputs")
    p.add_argument("code", choices=sorted(CODES))
    p.add_argument("template", help="vasp/aims/abinit: directory with the reference inputs; qe: pw.x input; "
                                    "generic: structure file")
    p.add_argument("-o", "--out", default="elastic", help="output directory (default ./elastic)")
    p.add_argument("--relax-ions", action="store_true", help="relax the ions in each strained cell")
    p.add_argument("--format", default=None, help="generic: ASE format of the structure files")
    p.add_argument("--filename", default=None, help="generic: name of the structure file in each directory")
    p.add_argument("--extra", nargs="*", default=(), help="generic: further files copied into each directory")
    p.add_argument("--potcar", default="copy", choices=["copy", "link"], help="vasp: copy or symlink POTCAR")
    p.add_argument("--fix-incar", action="store_true",
                   help="vasp: apply all recommended settings that the INCAR does not meet (EDIFF, ENCUT, PREC, "
                        "LREAL); without it they are only reported")
    p.add_argument("--ionic", action="store_true",
                   help="also displace the symmetry-irreducible atoms (relaxed-ion tensors, Born charges and "
                        "piezoelectric response from clamped cells, as VASP IBRION=6)")
    p.add_argument("--displacement", type=float, default=0.01, help="--ionic: atomic step in Angstrom")
    p.add_argument("--md", action="store_true",
                   help="finite T: fixed-cell molecular dynamics in each strained cell, as set up in the template "
                        "(vasp: IBRION = 0; qe: calculation = 'md')")
    _common(p)

    p = sub.add_parser("collect", help="elastic constants from the finished calculations")
    p.add_argument("directory", nargs="?", default="elastic")
    p.add_argument("--method", default="uncoupled", choices=["uncoupled", "designated", "lstsq"])
    p.add_argument("--output", default=None, help="name of the output file in each directory")
    p.add_argument("--format", default=None, help="generic: ASE format of the output file")
    p.add_argument("--md", action="store_true", help="average the stress over the ionic steps (finite T)")
    p.add_argument("--skip", type=int, default=0, help="--md: steps discarded for equilibration")
    p.add_argument("--samples", default=None,
                   help="glob of per-snapshot sub-directories in each state directory (e.g. 'conf-*'); "
                        "their stresses are averaged")
    p.add_argument("--kinetic-temperature", type=float, default=None,
                   help="add the ideal-gas term n kB T (stresses without kinetic part, e.g. VASP MD)")
    p.add_argument("--internal-strain", action="store_true", help="print the internal strain tensor")
    p.add_argument("--piezo", action="store_true",
                   help="piezoelectric tensor from the Berry-phase dipoles of the strained cells (VASP LCALCPOL)")
    p.add_argument("--born", default=None, metavar="OUTCAR",
                   help="Born charges and clamped-ion piezoelectric tensor from a VASP LEPSILON/LCALCEPS run of "
                        "the reference structure (with --ionic data)")
    p.add_argument("--stress-measure", default="cauchy", choices=["cauchy", "kirchhoff"],
                   help="differentiate the Cauchy stress (default) or J*sigma (as VASP IBRION=6)")

    p = sub.add_parser("vasp-tensor", help="analyse the elastic tensor printed by VASP IBRION=6")
    p.add_argument("outcar")
    p.add_argument("--which", default="total", choices=["total", "clamped", "ionic"])

    a = ap.parse_args(argv)
    kw = dict(scheme=a.scheme, delta=a.delta, stencil=a.stencil, symprec=a.symprec, frame=a.frame) \
        if a.cmd in ("info", "generate") else {}

    if a.cmd == "info":
        atoms = read_structure(a.code, a.structure, a.format)
        print(ElasticSetup(atoms, **{**kw, "stencil": a.stencil or "four-point"}).describe())
    elif a.cmd == "generate":
        if a.fix_incar and a.code != "vasp":
            ap.error("--fix-incar is for vasp")
        extra = dict(potcar=a.potcar, fix_incar=a.fix_incar) if a.code == "vasp" else {}
        setup, notes = generate(a.code, a.template, a.out, relax_ions=a.relax_ions, fmt=a.format,
                                filename=a.filename, extra_files=a.extra, ionic=a.ionic,
                                displacement=a.displacement, md=a.md, **kw, **extra)
        print(setup.describe())
        print(f"\n{len(setup.states)} calculations written to {a.out}/ (list in {a.out}/STATES)")
        for n in notes:
            print("NOTE:", n)
    elif a.cmd == "collect":
        r = collect(a.directory, method=a.method, output=a.output, fmt=a.format, md=a.md, skip=a.skip,
                    kinetic_temperature=a.kinetic_temperature, samples=a.samples,
                    stress_measure=a.stress_measure, polarization=a.piezo, born_file=a.born)
        print(r.summary())
        if a.internal_strain and r.internal_strain is not None:
            print("internal strain tensor dF/de (eV/A), columns xx yy zz yz xz xy:")
            for n, L in enumerate(r.internal_strain):
                for i, ax in enumerate("xyz"):
                    print(f"  ion {n + 1:3d} {ax} " + " ".join(f"{v:10.4f}" for v in L[i]))
        print(f"\nwritten {a.directory}/elastic.json")
    elif a.cmd == "vasp-tensor":
        from .io.vasp import read_elastic_tensor
        which = {"total": "TOTAL ELASTIC MODULI", "clamped": "SYMMETRIZED ELASTIC MODULI",
                 "ionic": "ELASTIC MODULI CONTR FROM IONIC RELAXATION"}[a.which]
        C = read_elastic_tensor(a.outcar, which)
        print(f"{which} (GPa, Voigt xx yy zz yz xz xy):")
        for row in C:
            print("  " + " ".join(f"{v:9.2f}" for v in row))
        if a.which != "ionic":
            stable, ev = born_stable(C)
            print(f"Born stable: {stable}")
            if stable:
                print(json.dumps({k: round(v, 3) for k, v in vrh(C).items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
