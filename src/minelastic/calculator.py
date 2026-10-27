# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Elastic constants with an ASE calculator: every state of an ElasticSetup computed in memory.

ASE returns the stress in eV/A^3, tension positive, Voigt order xx yy zz yz xz xy; it is converted to GPa
(x 160.2177) and passed to ElasticSetup.compute with the forces of the clamped cells and, if a dipole function
is given, the dipole moment of every cell.
"""
import numpy as np

from .core import ElasticSetup

EV_PER_A3_IN_GPA = 160.21766208


def calculate(atoms, calc, relax_ions=False, method="uncoupled", fmax=1e-4, steps=500, optimizer=None,
              dipole=None, **setup_kwargs):
    """Elastic tensor of atoms from the stresses of every strained cell, computed with an ASE calculator.

    Builds ElasticSetup(atoms, **setup_kwargs), computes the stress (eV/A^3 from calc, converted to GPa,
    tension positive) and the forces (eV/A) of every state with calc, and returns ElasticSetup.compute of them.

    Parameters
    ----------
    atoms : ase.Atoms
        Reference structure, relaxed and stress-free for the textbook constants.
    calc : ase calculator
        Calculator that returns the stress; it is attached to every strained cell.
    relax_ions : bool
        Relax the atoms of each strained cell at fixed cell (optimizer, default BFGS, to fmax eV/A within steps
        steps); the constants are then the relaxed-ion constants and no forces are passed on. Otherwise the
        cells are clamped and their forces give the internal-strain tensor.
    method : str
        Extraction method of ElasticSetup.compute: 'uncoupled', 'designated' or 'lstsq'.
    fmax, steps, optimizer :
        Convergence (eV/A), step limit and ASE optimizer class of the relaxations.
    dipole : callable, optional
        Function of an Atoms object (calculator attached, after any relaxation) that returns the dipole moment
        of the cell in e*Angstrom, physical sign. It gives the piezoelectric tensor e = dP/de (clamped-ion for
        clamped cells, relaxed-ion with relax_ions) and, with ionic=True, the Born charges, the total e and d.
    **setup_kwargs :
        Options of ElasticSetup: scheme, delta, stencil, symprec, frame, ionic, displacement.

    Returns
    -------
    ElasticResult
        With ionic=True, result.C_relaxed is the relaxed-ion stiffness from the clamped cells and the
        displacements.
    """
    setup = ElasticSetup(atoms, **setup_kwargs)
    if relax_ions and setup.ionic:
        raise ValueError("ionic=True gives the relaxed-ion constants from clamped cells; do not combine it with "
                         "relax_ions=True")
    stresses, forces, dipoles = {}, {}, {}
    for st in setup.states:
        at = setup.strained_atoms(st)
        at.calc = calc
        if relax_ions:
            from ase.optimize import BFGS
            (optimizer or BFGS)(at, logfile=None).run(fmax=fmax, steps=steps)
        # ASE stress: eV/A^3, xx yy zz yz xz xy, positive in tension
        stresses[st["id"]] = at.get_stress(voigt=True) * EV_PER_A3_IN_GPA
        forces[st["id"]] = at.get_forces()
        if dipole is not None:
            dipoles[st["id"]] = np.asarray(dipole(at), float)
    return setup.compute(stresses, forces=None if relax_ions else forces, method=method,
                         polarization=dipoles if dipole is not None else None)
