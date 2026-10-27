# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Elastic constants of fcc and hcp Cu with ASE's EMT potential (seconds on a laptop)."""
from ase.build import bulk
from ase.calculators.emt import EMT

import minelastic

fcc = minelastic.calculate(bulk("Cu", "fcc", a=3.59, cubic=True), EMT(), delta=0.005)
print(fcc.summary(), "\n")
hcp = bulk("Cu", "hcp", a=2.55, c=4.16)
clamped = minelastic.calculate(hcp, EMT(), delta=0.005)
relaxed = minelastic.calculate(hcp, EMT(), delta=0.005, relax_ions=True, fmax=1e-6)
print(relaxed.summary())
print(f"clamped-ion c11 = {clamped.as_dict()['c11']:.2f} GPa, c12 = {clamped.as_dict()['c12']:.2f} GPa")
print("internal strain of atom 1 (eV/A per unit strain, columns xx yy zz yz xz xy):")
print(clamped.internal_strain[0].round(4))
