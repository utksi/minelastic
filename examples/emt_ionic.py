# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Relaxed-ion elastic constants of hcp Cu (ASE's EMT potential) in three ways (a few seconds).

The two atoms of the hcp cell relax when the cell is strained, which softens c11, c12 and c66. minelastic
obtains the relaxed-ion constants

1. by relaxing the atoms in every strained cell of the minimal set (relax_ions=True);
2. from clamped single points only (ionic=True): the strained cells give the clamped-ion constants C0 and,
   together with the symmetry-irreducible atomic displacements, the Gamma-point force constants Phi and the
   internal strain tensor Lambda; then C = C0 - Lambda^T Phi^+ Lambda / Omega, as VASP's IBRION = 6 does;
3. as 1, but with the six single-component strains instead of the minimal set (scheme="asess").

Requires ASE >= 3.23 (FrechetCellFilter).
"""
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS

import minelastic

atoms = bulk("Cu", "hcp", a=2.55, c=4.16)
atoms.calc = EMT()
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-7)    # stress-free reference cell
print("reference: a = %.4f A, c = %.4f A" % (atoms.cell.cellpar()[0], atoms.cell.cellpar()[2]))

relaxed = minelastic.calculate(atoms, EMT(), relax_ions=True, fmax=1e-6)
ionic = minelastic.calculate(atoms, EMT(), ionic=True)
six = minelastic.calculate(atoms, EMT(), relax_ions=True, fmax=1e-6, scheme="asess")

n_strained = sum(1 for st in ionic.setup.states if "atom" not in st)
n_displaced = len(ionic.setup.states) - n_strained
print(f"\n{'':6s} {'clamped C0':>11s} {'relaxed (1)':>12s} {'ionic (2)':>10s} {'six (3)':>9s}   GPa")
c_ionic = dict(zip(ionic.setup.form.labels, ionic.setup.form.project(ionic.C_relaxed)))
for k, c0 in ionic.as_dict().items():
    print(f"{k:6s} {c0:11.2f} {relaxed.as_dict()[k]:12.2f} {c_ionic[k]:10.2f} {six.as_dict()[k]:9.2f}")
print(f"\n(1) {len(relaxed.setup.states)} relaxations; (2) {n_strained} strained and {n_displaced} displaced-atom "
      f"single points; (3) {len(six.setup.states)} relaxations")
# (1) and (3) contain finite-strain (anharmonic) terms of order delta^2, (2) is the linear response
print("largest difference (2) - (1): %.3f GPa" % np.abs(ionic.C_relaxed - relaxed.C).max())
print("\nGamma-point force constants of atom 1 on itself (eV/A^2):")
print(ionic.force_constants[0, :, 0, :].round(4))
