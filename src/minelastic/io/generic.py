# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Any code with an ASE writer for its structure file and an ASE reader for its output.

Each state directory gets the strained structure, written with ase.io.write in the given format (e.g.
'vasp', 'castep-cell', 'cif'), and copies of further input files. The result is the last image of the output
file read with ase.io.read: stress in eV/A^3 x 160.2177 = GPa (ASE sign, tension positive), forces in eV/A.
"""
import os
import shutil

from ase.io import read, write

from ..calculator import EV_PER_A3_IN_GPA


def read_structure(template, fmt=None):
    """Atoms of the template structure file (ASE format fmt, or guessed from the name)."""
    return read(template, format=fmt)


def prepare(setup, out_dir, filename, fmt=None, extra_files=(), relax_ions=False):
    """One directory per state with the strained structure in filename (format fmt) and the extra files."""
    notes = ["generic interface: only the structure is written; make sure the inputs request the stress "
             + ("and relax the ions at fixed cell" if relax_ions else "and keep the ions fixed")]
    for st in setup.states:
        d = os.path.join(out_dir, st["id"])
        os.makedirs(d, exist_ok=True)
        write(os.path.join(d, filename), setup.strained_atoms(st), format=fmt)
        for f in extra_files:
            shutil.copy(f, d)
    return notes


def read_result(run_dir, output, fmt=None):
    """(stress in GPa, forces in eV/A or None) of the last image of run_dir/output read by ASE."""
    at = read(os.path.join(run_dir, output), format=fmt, index=-1)
    try:
        forces = at.get_forces()
    except Exception:
        forces = None
    return at.get_stress(voigt=True) * EV_PER_A3_IN_GPA, forces
