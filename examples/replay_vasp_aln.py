# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Bring your own stresses: the Python API on real VASP data (wurtzite AlN, PBE; under a second).

tests/data_berry/aln_berry_minimal.json holds what VASP computed for the four strained cells of the minimal
set (clamped ions, delta = 0.01, LCALCPOL = .TRUE.): the stress and the Berry-phase dipole of each cell.
The stresses of any code, run in any way, are used like this.
"""
import json
import os

import numpy as np
from ase import Atoms

from minelastic import ElasticSetup

here = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(here, "..", "tests", "data_berry", "aln_berry_minimal.json")))
s = data["setup"]["input_structure"]
atoms = Atoms(s["symbols"], cell=s["cell"], scaled_positions=s["scaled_positions"], pbc=True)

setup = ElasticSetup(atoms, delta=0.01, stencil="central")   # hexagonal minimal set: 2 patterns at +-delta, as run
for st in setup.states:
    # run your code on setup.strained_atoms(st): Cauchy stress in GPa, tension positive, xx yy zz yz xz xy
    print(st["id"], "applied strain (Voigt, engineering shear):", np.round(st["strain"], 4))

res = setup.compute(data["stress_GPa"], polarization=data["dipole_eA"])
print()
print(res.summary())

e, dfpt = res.piezo, data["dfpt_clamped_piezo_C_m2"]    # e: (3, 6), C/m^2, Voigt columns xx yy zz yz xz xy
print("\nclamped-ion piezoelectric tensor (C/m^2): this run (Berry phase, 4 cells) vs DFPT (LEPSILON)")
for name, (i, j) in {"e33": (2, 2), "e31": (2, 0), "e15": (0, 4)}.items():
    print(f"  {name} = {e[i, j]:7.3f}   {dfpt[name]:7.3f}")
