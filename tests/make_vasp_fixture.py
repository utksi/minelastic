# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Compact test fixture from a VASP IBRION=6, ISIF>=3 run (stock VASP or LMINIMAL_ELASTIC).

Stores the reference structure, POTIM, the strain patterns of the run (in VASP's order), the stresses of the
+-POTIM strained cells (raw "in kB" lines; strain DOFs come first, displacement 1 = +POTIM, 2 = -POTIM) and
the tensors printed by VASP. Usage:

    python make_vasp_fixture.py RUN_DIR OUT.json --label "AlN, stock VASP" [--paper c11=114.2,...]
"""
import argparse
import json
import re

import numpy as np
from ase.io import read

# stock IBRION=6 distortions (finite_diff.F, MAKE_DISPLACED_ID): xx, yy, zz and simple shears xy, yz, zx.
# Their symmetric parts as engineering Voigt strains (xx yy zz yz xz xy) per unit step:
STOCK = [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0]]
VASP_TO_VOIGT = [0, 1, 2, 4, 5, 3]
# the distortion matrices themselves (A_new = (1 + step * D) A on column vectors), DATA DISTORTION of VASP
STOCK_D = [np.diag([1., 0, 0]), np.diag([0, 1., 0]), np.diag([0, 0, 1.]),
           np.array([[0, 0, 0], [1., 0, 0], [0, 0, 0]]), np.array([[0, 0, 0], [0, 0, 0], [0, 1., 0]]),
           np.array([[0, 0, 1.], [0, 0, 0], [0, 0, 0]])]


def tensor(lines, head):
    """Last 6x6 block 'head' of an OUTCAR in GPa and Voigt order: 0.1 M[P, P] of the printed kBar matrix."""
    idx = [i for i, l in enumerate(lines) if l.strip().startswith(head)]
    if not idx:
        return None
    i = idx[-1]
    while not lines[i].strip().startswith("XX"):
        i += 1
    M = np.array([[float(x) for x in lines[i + r].split()[1:7]] for r in range(6)]) / 10.0
    return M[np.ix_(VASP_TO_VOIGT, VASP_TO_VOIGT)].round(4).tolist()


def main():
    """Reads POSCAR and OUTCAR of an IBRION = 6 run: POTIM, the strain patterns (stock: the six single strains;
    patched runs: as printed), the raw stress of every +-POTIM strained cell with J = det F, and the printed
    tensors; writes them to the JSON fixture.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--label", required=True)
    ap.add_argument("--paper", default=None, help="published constants, e.g. c11=114.2,c12=59.8")
    ap.add_argument("--paper-upper", default=None,
                    help="published tensor as the 21 upper-triangle values c11 c12 ... c66 (comma separated, 0 where zero)")
    a = ap.parse_args()
    L = open(f"{a.run}/OUTCAR", errors="replace").read().splitlines()
    atoms = read(f"{a.run}/POSCAR", format="vasp")
    potim = float(next(l for l in L if l.strip().startswith("POTIM")).split("=")[1].split()[0])
    nfree = int(next(l for l in L if l.strip().startswith("NFREE")).split("=")[1].split()[0])
    assert nfree == 2, "only NFREE = 2 runs are supported here"
    minimal = any("LMINIMAL_ELASTIC: minimal set" in l or "MINIMAL ELASTIC CONSTANTS METHOD ACTIVE" in l for l in L)
    nstrain = int(next(l for l in L if "Strain:" in l and "degrees of freedom" in l).split()[1])
    if minimal and any("LMINIMAL_ELASTIC: minimal set" in l for l in L):        # v2 prints the patterns
        i = next(k for k, l in enumerate(L) if "pattern   eps_xx" in l)
        pats = []
        for r in range(nstrain):
            e = [float(x) for x in L[i + 1 + r].split()[1:7]]                  # xx yy zz xy yz zx (tensor)
            pats.append([e[0], e[1], e[2], 2 * e[4], 2 * e[5], 2 * e[3]])
    elif minimal and nstrain == 6:                                            # v1 triclinic: VASP order, symmetric
        pats = STOCK
    elif minimal:                                                             # v1: the paper's table, input frame
        from minelastic import ElasticSetup
        pats = ElasticSetup(atoms, frame="input").patterns.tolist()
        assert len(pats) == nstrain
    else:
        pats = STOCK
    from minelastic.voigt import strain_to_tensor
    dist = [STOCK_D[k] if not minimal else strain_to_tensor(pats[k]) for k in range(nstrain)]
    kb = [np.array([float(x) for x in l.split()[2:8]]) for l in L if l.strip().startswith("in kB")]
    # step 0: reference; then DOF 1..nstrain, each +POTIM, -POTIM
    strained = []
    for k in range(nstrain):
        for s, sign in ((1 + 2 * k, 1), (2 + 2 * k, -1)):
            v = kb[s]
            F = np.eye(3) + sign * potim * dist[k]
            strained.append(dict(pattern=k, sign=sign, F=F.tolist(), J=float(np.linalg.det(F)),
                                 stress_GPa=(-0.1 * v[VASP_TO_VOIGT]).round(6).tolist()))   # physical sign
    d = dict(label=a.label, source=a.run, potim=potim, minimal=minimal, nstrain=nstrain, patterns=pats,
             structure=dict(symbols=atoms.get_chemical_symbols(), cell=atoms.cell[:].tolist(),
                            scaled_positions=atoms.get_scaled_positions().tolist()),
             reference_stress_GPa=(-0.1 * kb[0][VASP_TO_VOIGT]).round(6).tolist(), strained=strained,
             vasp_clamped=tensor(L, "SYMMETRIZED ELASTIC MODULI"), vasp_total=tensor(L, "TOTAL ELASTIC MODULI"),
             vasp_ionic=tensor(L, "ELASTIC MODULI CONTR FROM IONIC RELAXATION"))
    if a.paper:
        d["paper"] = {k: float(v) for k, v in (kv.split("=") for kv in a.paper.split(","))}
    if a.paper_upper:
        d["paper_upper"] = [float(x) for x in a.paper_upper.split(",")]
    json.dump(d, open(a.out, "w"), indent=0)
    print(f"{a.out}: {len(strained)} strained cells, minimal={minimal}")


if __name__ == "__main__":
    main()
