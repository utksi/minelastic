# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""ABINIT: one directory per state with the strained POSCAR and the template input; stress and forces from
the ABINIT output.

The template directory holds the reference structure as POSCAR and one input file (*.abi) that reads it with
structure "poscar:POSCAR" and sets everything else (pseudopotentials in the order of the POSCAR species,
cutoff, k points, ...). ABINIT's POSCAR reader takes the element of each atom from the end of its
coordinate line, so the POSCARs are written with the symbols there. The initial moments for the symmetry
analysis are those of spinat in the input. The output <input stem>.abo is read with ASE: stress in Ha/bohr^3
converted to GPa (ABINIT's sign, tension positive), forces in eV/A.
"""
import glob
import os
import re
import shutil

import numpy as np
from ase.io import read

from ..calculator import EV_PER_A3_IN_GPA
from ..symmetry import SPECIES
from .vasp import poscar_species_blocks


def _input(template):
    """The single *.abi file of the template directory."""
    files = sorted(glob.glob(os.path.join(template, "*.abi")))
    if len(files) != 1:
        raise FileNotFoundError(f"{template} must contain exactly one ABINIT input file *.abi (found {len(files)})")
    return files[0]


def _uncommented(text):
    """Input text in lower case, without comments (text after # or !)."""
    return re.sub(r"[#!].*", "", text).lower()


def _int_var(low, var, default=None):
    """Integer value of variable var in the uncommented input low, or default; dataset-specific (var1, var:,
    var+, var?) or conflicting values are an error."""
    found = re.findall(r"\b%s([0-9:+?]*)\s+(-?\d+)" % var, low)
    if any(suffix for suffix, _ in found) or len({v for _, v in found}) > 1:
        raise ValueError(f"{var}: dataset-specific or conflicting values are not supported")
    return int(found[0][1]) if found else default


def _numbers(low, var, count):
    """The count numbers after variable var in the uncommented input low, or None if var is absent; n*x stands
    for n copies of x, and Fortran exponents (2.0d0) are read."""
    found = re.findall(r"\b%s([0-9:+?]*)\s" % var, low)
    if any(found) or len(found) > 1:
        raise ValueError(f"{var}: dataset-specific or repeated values are not supported")
    if not found:
        return None
    vals = []
    for tok in low[re.search(r"\b%s\s" % var, low).end():].split():
        m = re.fullmatch(r"(?:(\d+)\*)?([-+]?(?:\d+\.?\d*|\.\d+)(?:[ed][-+]?\d+)?)", tok)
        if m is None or len(vals) >= count:
            break
        vals += [float(m.group(2).replace("d", "e"))] * int(m.group(1) or 1)
    if len(vals) != count:
        raise ValueError(f"{var}: expected {count} numbers, read {len(vals)} (plain numbers and n*x only)")
    return np.array(vals)


def _spin_moments(low, natom):
    """Initial moments from spinat (3 numbers per atom, Cartesian), as ABINIT uses them: the z components for
    collinear spins (nspden = 2: nsppol = 2, or the antiferromagnetic mode nsppol = 1), the vectors for
    noncollinear magnetism (nspden = 4), None without spin polarization (nspden = 1) or spinat. Defaults
    nsppol = 1, nspden = nsppol."""
    nspden = _int_var(low, "nspden", _int_var(low, "nsppol", 1))
    s = _numbers(low, "spinat", 3 * natom) if nspden in (2, 4) else None
    if s is None:
        return None
    s = s.reshape(natom, 3)
    return s if nspden == 4 else s[:, 2]


def read_structure(template):
    """Atoms of template/POSCAR, with the index of each atom's POSCAR species block in atoms.arrays[SPECIES]
    and the initial magnetic moments from spinat of the template input (_spin_moments)."""
    atoms = read(os.path.join(template, "POSCAR"), format="vasp")
    blocks = poscar_species_blocks(os.path.join(template, "POSCAR"))
    if blocks is not None and len(blocks) == len(atoms):
        atoms.new_array(SPECIES, blocks)
    moments = _spin_moments(_uncommented(open(_input(template)).read()), len(atoms))
    if moments is not None and np.any(np.abs(moments) > 1e-8):
        atoms.set_initial_magnetic_moments(moments)
    return atoms


def write_poscar(path, atoms, comment):
    """POSCAR in direct coordinates with the element symbol at the end of every coordinate line.

    Species blocks are the runs of equal consecutive symbols, so the order of the atoms is kept.
    """
    symbols = atoms.get_chemical_symbols()
    blocks = []
    for s in symbols:
        if blocks and blocks[-1][0] == s:
            blocks[-1][1] += 1
        else:
            blocks.append([s, 1])
    out = [comment, "1.0"]
    out += ["  " + " ".join(f"{x:21.16f}" for x in v) for v in atoms.cell[:]]
    out.append(" ".join(s for s, _ in blocks))
    out.append(" ".join(str(n) for _, n in blocks))
    out.append("Direct")
    for f, s in zip(atoms.get_scaled_positions(wrap=False), symbols):
        out.append("  " + " ".join(f"{x:19.16f}" for x in f) + f"  {s}")
    open(path, "w").write("\n".join(out) + "\n")


def prepare(setup, template, out_dir, relax_ions=False):
    """One directory per state with the strained POSCAR and a copy of the template input.

    The template must not have two species of one element (ABINIT's POSCAR reader identifies species by
    element), nor dataset-specific or conflicting ionmov/optcell. The cell must stay fixed (no optcell other
    than 0). With relax_ions, 'ionmov 2', 'ntime 100' and
    'tolmxf 1.0d-6' (Ha/bohr) are appended unless ionmov is set; ionmov 0 with relax_ions, or ionmov other
    than 0 without it, is an error (text after # or ! is a comment). If the cell was rotated to its standard
    orientation, Cartesian vectors of the input (noncollinear spinat, efield, dfield, zeemanfield) are
    refused. Returns notes for the user.
    """
    inp = _input(template)
    text = open(inp).read()
    low = _uncommented(text)
    notes = []
    sp = setup.input_atoms.arrays.get(SPECIES)
    if sp is not None and len(set(zip(setup.input_atoms.numbers.tolist(), sp.tolist()))) > len(set(setup.input_atoms.numbers)):
        raise ValueError("ABINIT's POSCAR reader cannot distinguish two species of one element (repeated element "
                         "blocks in POSCAR); for magnetic sublattices use spinat")
    for var in ("ionmov", "optcell"):
        _int_var(low, var)
    if not np.allclose(setup.rotation, np.eye(3)):
        # WHY: Cartesian vectors of the input would have to be rotated with the cell
        vec = [v for v in ("efield", "dfield", "zeemanfield") if re.search(r"\b%s[0-9:+?]*\s" % v, low)]
        if setup.input_atoms.get_initial_magnetic_moments().ndim == 2:
            vec.append("noncollinear spinat")
        if vec:
            raise ValueError("the cell was rotated to its standard orientation (frame='auto'); "
                             f"{vec[0]} is direction-dependent and is not rotated: use frame='input' (--frame input)")
    if not re.search(r'structure\s+"poscar:poscar"', low):
        notes.append(f'{os.path.basename(inp)} should read the structure with structure "poscar:POSCAR"')
    if re.search(r"\boptcell\s+[1-9]", low):
        raise ValueError("the template relaxes the cell (optcell); the strained cells must stay fixed")
    m = re.search(r"\bionmov\s+(\d+)", low)
    if m and int(m.group(1)) != 0 and not relax_ions:
        raise ValueError("the template relaxes the ions (ionmov): remove it for clamped ions, or use relax_ions")
    if m and int(m.group(1)) == 0 and relax_ions:
        raise ValueError("relax_ions with ionmov 0 in the template: remove ionmov or set a relaxation algorithm")
    if relax_ions and not m:
        text += "\nionmov 2\nntime 100\ntolmxf 1.0d-6\n"
    for st in setup.states:
        d = os.path.join(out_dir, st["id"])
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), setup.strained_atoms(st), f"minelastic {st['id']}")
        open(os.path.join(d, os.path.basename(inp)), "w").write(text)
    shutil.copy(inp, os.path.join(out_dir, "template.abi"))
    return notes


def read_result(run_dir, output=None):
    """(stress in GPa, Voigt, tension positive; forces in eV/A or None) of the last step of the ABINIT output.

    The output is run_dir/output, or the most recent *.abo file of run_dir.
    """
    if output is None:
        files = sorted(glob.glob(os.path.join(run_dir, "*.abo*")), key=os.path.getmtime)
        if not files:
            raise FileNotFoundError(f"no ABINIT output *.abo in {run_dir}")
        path = files[-1]
    else:
        path = os.path.join(run_dir, output)
    at = read(path, format="abinit-out")
    try:
        forces = at.get_forces()
    except Exception:
        forces = None
    return at.get_stress(voigt=True) * EV_PER_A3_IN_GPA, forces
