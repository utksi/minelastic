# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Fixture of the documented VASP run (docs/how-to/piezoelectricity.md): wurtzite AlN with --ionic and LCALCPOL.

`minelastic generate vasp vasp_aln -o el --ionic --displacement 0.015` with LCALCPOL = .TRUE., all 12 single
points run with stock VASP 6.6.1 (piezo_tests/aln_min_ionic). Stored per cell, as minelastic's VASP readers
return them: the stress (GPa, Voigt, tension positive), the forces (eV/A) and the Berry-phase dipole (e A).
The reference tensors are those of stock VASP 6.6.1 IBRION = 6, ISIF = 3, LEPSILON (six strains) of the same
input, as in data_ionic/aln_ionic.json, plus the internal-strain tensor that run prints.

    PIEZO_TESTS=/path/to/piezo_tests PATCH_TESTS=/path/to/patch_tests python make_aln_run_fixture.py
"""
import json
import os

import numpy as np

from minelastic.io.vasp import read_dipole, read_result

RUN = os.path.join(os.environ["PIEZO_TESTS"], "aln_min_ionic")      # VASP runs (not in this repository)
STOCK = os.path.join(os.environ["PATCH_TESTS"], "aln_std", "pristine", "OUTCAR")
HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "data_ionic", "aln_run_ionic.json")


def internal_strain(outcar, natoms):
    """Blocks 'INTERNAL STRAIN TENSOR FOR ION n' of a VASP OUTCAR (eV/A): (natoms, 3, 6), rows x y z, columns
    reordered from VASP's X Y Z XY YZ ZX to Voigt xx yy zz yz xz xy."""
    lines = open(outcar, errors="replace").read().splitlines()
    out = []
    for i, line in enumerate(lines):
        if "INTERNAL STRAIN TENSOR FOR ION" in line and len(out) < natoms:
            rows = [[float(x) for x in lines[i + 3 + k].split()[1:7]] for k in range(3)]
            out.append(np.array(rows)[:, [0, 1, 2, 4, 5, 3]])
    return np.array(out)


def main():
    """Writes data_ionic/aln_run_ionic.json from the run directory and the stock references."""
    setup = json.load(open(os.path.join(RUN, "minelastic.json")))
    setup["template"] = "examples/vasp_aln with LCALCPOL = .TRUE."      # instead of the local absolute path
    states = open(os.path.join(RUN, "STATES")).read().split()
    stress, forces, dipole = {}, {}, {}
    for sid in states:
        s, f = read_result(os.path.join(RUN, sid))
        stress[sid], forces[sid] = np.round(s, 6).tolist(), np.round(f, 6).tolist()
        dipole[sid] = np.round(read_dipole(os.path.join(RUN, sid)), 8).tolist()
    stock = json.load(open(os.path.join(HERE, "data_ionic", "aln_ionic.json")))["stock"]
    stock["internal_strain"] = internal_strain(STOCK, len(setup["input_structure"]["symbols"])).tolist()
    data = {"label": "wurtzite AlN (PBE): minelastic --ionic --displacement 0.015 with LCALCPOL, 12 single points "
                     "with stock VASP 6.6.1; references from stock VASP 6.6.1 IBRION = 6, ISIF = 3, LEPSILON",
            "source": "piezo_tests/aln_min_ionic; stock: patch_tests/aln_std/pristine",
            "setup": setup, "stress_GPa": stress, "forces_eV_A": forces, "dipole_eA": dipole, "stock": stock}
    json.dump(data, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}: {len(states)} cells")


if __name__ == "__main__":
    main()
