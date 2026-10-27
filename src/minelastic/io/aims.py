# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""FHI-aims: geometry.in per state (written by ASE, fractional coordinates), control.in copied from the
template directory; stress and forces read from aims.out with ASE (eV/A^3 -> GPa, tension positive)."""
import os
import shutil

from ase.io import read, write

from ..calculator import EV_PER_A3_IN_GPA

DEFAULT_OUTPUT = "aims.out"


def read_structure(template):
    """Atoms of template/geometry.in (initial moments are read by ASE)."""
    return read(os.path.join(template, "geometry.in"), format="aims")


def prepare(setup, template, out_dir, relax_ions=False):
    """One directory per state with geometry.in of the strained cell and the template's control.in.

    The active directives of control.in (text after '#' is a comment) must request the stress
    (compute_analytical_stress or compute_numerical_stress .true.) and keep the cell fixed (no relax_unit_cell
    other than none). With relax_ions, 'relax_geometry bfgs 1E-3' is appended unless relax_geometry is set;
    without it, an active relax_geometry (other than none) is an error. Returns notes for the user.
    """
    control = open(os.path.join(template, "control.in")).read()
    notes = []
    active = {}
    for line in control.splitlines():
        f = line.split("#", 1)[0].lower().split()
        if f:
            active[f[0]] = f[1:] or [""]
    if not any(active.get(k, [""])[0].startswith(".t") for k in ("compute_analytical_stress",
                                                                   "compute_numerical_stress")):
        raise ValueError("control.in must request the stress: compute_analytical_stress .true.")
    if active.get("relax_unit_cell", ["none"])[0] != "none":
        raise ValueError("control.in relaxes the unit cell (relax_unit_cell); the strained cells must stay fixed")
    relax = active.get("relax_geometry", ["none"])[0] != "none"
    if relax and not relax_ions:
        raise ValueError("control.in relaxes the ions (relax_geometry): remove it for clamped ions, or use "
                         "relax_ions (--relax-ions)")
    if relax_ions and not relax:
        control += "\nrelax_geometry bfgs 1E-3\n"
    for st in setup.states:
        d = os.path.join(out_dir, st["id"])
        os.makedirs(d, exist_ok=True)
        write(os.path.join(d, "geometry.in"), setup.strained_atoms(st), format="aims", scaled=True)
        open(os.path.join(d, "control.in"), "w").write(control)
    return notes


def read_result(run_dir, output=None):
    """(stress in GPa, Voigt, tension positive; forces in eV/A or None) of the last step in aims.out."""
    at = read(os.path.join(run_dir, output or DEFAULT_OUTPUT), format="aims-output", index=-1)
    try:
        forces = at.get_forces()
    except Exception:
        forces = None
    return at.get_stress(voigt=True) * EV_PER_A3_IN_GPA, forces
