# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Fixture for the ionic route (minelastic.ionic) from VASP IBRION = 6 runs of wurtzite AlN.

Raw data come from an IBRION = 6 run with LCALCPOL and the minimal strain set, in which the Berry-phase
dipole was computed in every cell (patch_tests/pz/aln_pol_min):
* per strained cell: the stress;
* per displaced-atom cell: the atom, the displacement, the forces, the stress and the dipole.

The reference values come from STOCK VASP runs of the same input (six standard strains, LEPSILON):
* the ionic relaxation contribution to C and the total C (patch_tests/aln_std/pristine);
* the DFPT Born charges, the clamped-ion piezoelectric tensor and the ionic piezoelectric contribution
  computed by VASP with those Born charges (same run).

    PATCH_TESTS=/path/to/patch_tests python make_ionic_fixture.py
"""
import json
import os

import numpy as np
from ase.io import read

from minelastic.io.vasp import read_born_charges, read_clamped_piezo, read_elastic_tensor

T = os.environ["PATCH_TESTS"]                 # directory of the VASP test runs (not in this repository)
RAW = f"{T}/pz/aln_pol_min"
STOCK = f"{T}/aln_std/pristine"
VASP_TO_VOIGT = [0, 1, 2, 4, 5, 3]
OUT = os.path.join(os.path.dirname(__file__), "data_ionic", "aln_ionic.json")


def piezo_block(outcar, head):
    """Last 3x6 piezoelectric block 'head' of an OUTCAR (C/m^2), columns reordered to Voigt."""
    L = open(outcar, errors="replace").read().splitlines()
    i = [k for k, l in enumerate(L) if head in l][-1]
    while not L[i].strip().startswith("---"):
        i += 1
    return np.array([[float(x) for x in L[i + 1 + r].split()[1:7]] for r in range(3)])[:, VASP_TO_VOIGT]


def main():
    """Collects the per-cell stresses, forces and dipoles of the patched IBRION = 6 run and the reference tensors
    of the stock LEPSILON run; writes tests/data_ionic/aln_ionic.json.
    """
    L = open(f"{RAW}/OUTCAR", errors="replace").read().splitlines()
    atoms = read(f"{RAW}/POSCAR")
    n = len(atoms)
    potim = float(next(l for l in L if l.strip().startswith("POTIM")).split("=")[1].split()[0])
    nstr = int(next(l for l in L if "Strain:" in l and "degrees of freedom" in l).split()[1])
    i = next(k for k, l in enumerate(L) if "pattern   eps_xx" in l)
    pats = []
    for r in range(nstr):
        e = [float(x) for x in L[i + 1 + r].split()[1:7]]
        pats.append([e[0], e[1], e[2], 2 * e[4], 2 * e[5], 2 * e[3]])
    stress = [(-0.1 * np.array([float(x) for x in l.split()[2:8]])[VASP_TO_VOIGT]).tolist()
              for l in L if l.strip().startswith("in kB")]
    pos, frc = [], []
    for k, l in enumerate(L):
        if "TOTAL-FORCE (eV/Angst)" in l:
            blk = np.array([[float(x) for x in L[k + 2 + m].split()[:6]] for m in range(n)])
            pos.append(blk[:, :3]); frc.append(blk[:, 3:])
    pel = [np.array([float(x) for x in l.split("(")[1].split(")")[0].split()]) for l in L if "p[elc]=(" in l]
    pio = [np.array([float(x) for x in l.split("(")[1].split(")")[0].split()]) for l in L if "p[ion]=(" in l]
    dip = [(a + b).tolist() for a, b in zip(pel, pio)]
    nsteps = len(stress)
    assert len(frc) == nsteps == len(dip), (len(frc), nsteps, len(dip))
    strained = [dict(pattern=k, sign=s, stress_GPa=stress[1 + 2 * k + j], dipole_eA=dip[1 + 2 * k + j])
                for k in range(nstr) for j, s in enumerate((1, -1))]
    cell0 = atoms.cell[:]
    displaced = []
    for t in range(1 + 2 * nstr, nsteps):
        d = pos[t] - pos[0]
        f = d @ np.linalg.inv(cell0)
        d = (f - np.round(f)) @ cell0
        m = int(np.argmax(np.linalg.norm(d, axis=1)))
        assert abs(np.linalg.norm(d[m]) - potim) < 1e-4, np.linalg.norm(d[m])
        displaced.append(dict(atom=m, displacement_A=d[m].round(8).tolist(), forces=frc[t].tolist(),
                              stress_GPa=stress[t], dipole_eA=dip[t]))
    ref = dict(
        C_ionic=read_elastic_tensor(f"{STOCK}/OUTCAR", "ELASTIC MODULI CONTR FROM IONIC RELAXATION").tolist(),
        C_total=read_elastic_tensor(f"{STOCK}/OUTCAR", "TOTAL ELASTIC MODULI").tolist(),
        C_clamped=read_elastic_tensor(f"{STOCK}/OUTCAR", "SYMMETRIZED ELASTIC MODULI").tolist(),
        born_dfpt=read_born_charges(f"{STOCK}/OUTCAR").tolist(),
        piezo_clamped_dfpt=read_clamped_piezo(f"{STOCK}/OUTCAR").tolist(),
        piezo_ionic_dfpt_born=piezo_block(f"{STOCK}/OUTCAR", "PIEZOELECTRIC TENSOR IONIC CONTR").tolist())
    d = dict(label="wurtzite AlN (PBE), IBRION=6 minimal set with Berry-phase dipoles in every cell; references "
                   "from stock VASP 6.6.1 (six strains, LEPSILON)", raw=os.path.relpath(RAW, T),
             reference=os.path.relpath(STOCK, T), potim=potim, patterns=pats,
             structure=dict(symbols=atoms.get_chemical_symbols(), cell=cell0.tolist(),
                            scaled_positions=atoms.get_scaled_positions().tolist()),
             reference_stress_GPa=stress[0], reference_dipole_eA=dip[0], reference_forces=frc[0].tolist(),
             strained=strained, displaced=displaced, stock=ref)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(d, open(OUT, "w"), indent=0)
    print(f"{OUT}: {len(strained)} strained, {len(displaced)} displaced cells",
          [(x["atom"], x["displacement_A"]) for x in displaced])


if __name__ == "__main__":
    main()
