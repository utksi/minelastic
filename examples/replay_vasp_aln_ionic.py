# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Elastic, internal-strain, Born-charge and piezoelectric tensors of AlN from stored VASP data (under a second).

tests/data_ionic/aln_run_ionic.json holds the run of docs/how-to/piezoelectricity.md: `minelastic generate vasp vasp_aln
--ionic --displacement 0.015` with LCALCPOL = .TRUE., 4 strained and 8 displaced-atom cells, all clamped single
points with VASP 6.6.1. Stored per cell: the stress, the forces and the Berry-phase dipole. From them:
C0 from the stresses of the strained cells, the internal-strain tensor from their forces (and again from the
stresses of the displaced cells), the force constants Phi and Born charges Z* from the displaced cells, then
C = C0 - L^T Phi^+ L / Omega, e = e0 + Z* Phi^+ L / Omega and d = e C^-1. The references are VASP's own
IBRION = 6, ISIF = 3, LEPSILON (six strains, DFPT) of the same input.
"""
import json
import os

import numpy as np

from minelastic import ElasticSetup, ionic

here = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(here, "..", "tests", "data_ionic", "aln_run_ionic.json")))
setup = ElasticSetup.from_dict(data["setup"])            # the setup that generated the 12 cells
res = setup.compute(data["stress_GPa"], forces=data["forces_eV_A"], polarization=data["dipole_eA"])

ref = {k: np.array(v) for k, v in data["stock"].items()}
e_ref = ref["piezo_clamped_dfpt"] + ref["piezo_ionic_dfpt_born"]
d_ref = ionic.strain_coefficients(e_ref, ref["C_total"])
L = ref["internal_strain"]                               # printed by VASP, eV/A, Voigt columns

print(f"{len(setup.states)} single points; Voigt order xx yy zz yz xz xy\n")
print("                              minelastic   VASP IBRION=6")
rows = [("c11 clamped ions (GPa)", res.C[0, 0], ref["C_clamped"][0, 0]),
        ("c33 clamped ions (GPa)", res.C[2, 2], ref["C_clamped"][2, 2]),
        ("c11 relaxed ions (GPa)", res.C_relaxed[0, 0], ref["C_total"][0, 0]),
        ("c33 relaxed ions (GPa)", res.C_relaxed[2, 2], ref["C_total"][2, 2]),
        ("c44 relaxed ions (GPa)", res.C_relaxed[3, 3], ref["C_total"][3, 3]),
        ("c66 relaxed ions (GPa)", res.C_relaxed[5, 5], ref["C_total"][5, 5]),
        ("Lambda Al, z zz (eV/A)", res.internal_strain[0, 2, 2], L[0, 2, 2]),
        ("Lambda Al, x xy (eV/A)", res.internal_strain[0, 0, 5], L[0, 0, 5]),
        ("Lambda Al, x xz (eV/A)", res.internal_strain[0, 0, 4], L[0, 0, 4]),
        ("Z*_xx Al (|e|)", res.born[0, 0, 0], ref["born_dfpt"][0, 0, 0]),
        ("Z*_zz Al (|e|)", res.born[0, 2, 2], ref["born_dfpt"][0, 2, 2]),
        ("e33 clamped ions (C/m^2)", res.piezo[2, 2], ref["piezo_clamped_dfpt"][2, 2]),
        ("e33 total (C/m^2)", res.piezo_total[2, 2], e_ref[2, 2]),
        ("e31 total (C/m^2)", res.piezo_total[2, 0], e_ref[2, 0]),
        ("e15 total (C/m^2)", res.piezo_total[0, 4], e_ref[0, 4]),
        ("d33 (pC/N)", res.piezo_d[2, 2], d_ref[2, 2]),
        ("d31 (pC/N)", res.piezo_d[2, 0], d_ref[2, 0]),
        ("d15 (pC/N)", res.piezo_d[0, 4], d_ref[0, 4])]
for name, mine, theirs in rows:
    print(f"  {name:26s} {mine:9.3f}   {theirs:9.3f}")
print(f"\ninternal-strain tensor, all atoms, largest difference to VASP: "
      f"{np.abs(res.internal_strain - L).max():.3f} eV/A from the forces of the strained cells, "
      f"{np.abs(res.internal_strain_displaced - L).max():.3f} eV/A from the stresses of the displaced atoms")
print(f"acoustic sum of the Born charges: {np.abs(res.born.sum(axis=0)).max():.4f} |e|")
