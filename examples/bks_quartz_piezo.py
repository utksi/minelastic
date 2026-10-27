# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Piezoelectric, Born-charge and relaxed-ion tensors of alpha-quartz with a rigid-ion model (about 10 s).

The BKS potential of silica (B. W. H. van Beest, G. J. Kramer and R. A. van Santen, Phys. Rev. Lett. 64, 1955
(1990)) has fixed charges q_Si = 2.4, q_O = -1.2 and Buckingham terms A exp(-b r) - C / r^6 for Si-O and O-O.
The calculator below sums every pair term within the cutoff r_c, the Coulomb term by damped shifted-force
summation (C. J. Fennell and J. D. Gezelter, J. Chem. Phys. 124, 234104 (2006)), and shifts each term so that
it and its force vanish at r_c. The dipole of a cell is d = sum_n q_n r_n, so the Born charges are the q_n and
the clamped-ion piezoelectric tensor e0 is zero: the whole piezoelectric response comes from the internal
relaxation.

Two routes give the relaxed-ion tensors:
1. ionic=True: clamped single points of the strained cells and of displaced atoms give C0, the internal-strain
   tensor L, the force constants Phi and Z*, then C = C0 - L^T Phi^+ L / Omega and e = e0 + Z* Phi^+ L / Omega;
2. relax_ions=True: the ions are relaxed in every strained cell, and the stress and the dipole of the relaxed
   cells are differentiated directly.
"""
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.filters import FrechetCellFilter
from ase.neighborlist import neighbor_list
from ase.optimize import BFGS
from ase.spacegroup import crystal
from scipy.special import erfc

import minelastic

KE = 14.399645          # e^2 / (4 pi eps0) in eV A


class BKS(Calculator):
    """BKS silica: E = 1/2 sum_{i != j, r_ij < r_c} v_ij(r_ij), with v_ij the shifted Coulomb and Buckingham
    terms of pair_terms; forces F_i = sum_j v'(r_ij) r_ij / r_ij (r_ij = r_j - r_i) and the stress (eV/A^3,
    tension positive) sigma = 1/(2V) sum_{i != j} v'(r_ij) r_ij (x) r_ij / r_ij."""

    implemented_properties = ["energy", "forces", "stress"]
    charges = {"Si": 2.4, "O": -1.2}
    buckingham = {frozenset(["Si", "O"]): (18003.7572, 4.87318, 133.5381),       # A (eV), b (1/A), C (eV A^6)
                  frozenset(["O"]): (1388.7730, 2.76000, 175.0000)}

    def __init__(self, cutoff=10.0, alpha=0.25, **kwargs):
        """cutoff r_c (A) of all pair terms; alpha (1/A) the damping of the Coulomb term."""
        super().__init__(**kwargs)
        self.cutoff, self.alpha = cutoff, alpha

    def pair_terms(self, r, qq, A, b, C):
        """(v, dv/dr) at distances r < r_c: v = k qq [erfc(a r)/r - erfc(a r_c)/r_c + g(r_c)(r - r_c)] with
        g(x) = erfc(a x)/x^2 + 2a/sqrt(pi) exp(-a^2 x^2)/x, plus B(r) - B(r_c) - B'(r_c)(r - r_c) with
        B(r) = A exp(-b r) - C/r^6."""
        a, rc = self.alpha, self.cutoff
        g = lambda x: erfc(a * x) / x**2 + 2 * a / np.sqrt(np.pi) * np.exp(-(a * x) ** 2) / x      # noqa: E731
        B = lambda x: A * np.exp(-b * x) - C / x**6                                                 # noqa: E731
        dB = lambda x: -A * b * np.exp(-b * x) + 6 * C / x**7                                       # noqa: E731
        v = KE * qq * (erfc(a * r) / r - erfc(a * rc) / rc + g(rc) * (r - rc)) + B(r) - B(rc) - dB(rc) * (r - rc)
        dv = KE * qq * (g(rc) - g(r)) + dB(r) - dB(rc)
        return v, dv

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        """Energy, forces and stress from all pairs within r_c (each pair counted twice by neighbor_list)."""
        super().calculate(atoms, properties, system_changes)
        i, j, r, D = neighbor_list("ijdD", self.atoms, self.cutoff)
        sym = np.array(self.atoms.get_chemical_symbols())
        q = np.array([self.charges[s] for s in sym])
        v, dv = np.zeros(len(r)), np.zeros(len(r))
        kind = np.array([frozenset([s, t]) for s, t in zip(sym[i], sym[j])])
        for key in set(kind):
            m = kind == key
            v[m], dv[m] = self.pair_terms(r[m], q[i[m]] * q[j[m]], *self.buckingham.get(key, (0.0, 1.0, 0.0)))
        forces = np.zeros((len(self.atoms), 3))
        np.add.at(forces, i, (dv / r)[:, None] * D)
        s = 0.5 * np.einsum("p,pa,pb->ab", dv / r, D, D) / self.atoms.get_volume()
        self.results = {"energy": 0.5 * v.sum(), "forces": forces, "stress": s[[0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]]}


def alpha_quartz():
    """alpha-quartz, P3_221: Si 3a (0.4697, 0, 2/3), O 6c (0.4135, 0.2669, 0.7858) (Levien et al. 1980, origin
    of the International Tables)."""
    return crystal(["Si", "O"], [(0.4697, 0.0, 2 / 3), (0.4135, 0.2669, 0.7858)], spacegroup=154,
                   cellpar=[4.916, 4.916, 5.405, 90, 90, 120])


atoms = alpha_quartz()
atoms.calc = BKS()
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-7, steps=4000)       # stress-free reference
q = np.array([BKS.charges[s] for s in atoms.get_chemical_symbols()])


def dipole(a):
    """Dipole of the cell, d = sum_n q_n r_n (e A); the positions are continuous, so there is no branch issue."""
    return q @ a.positions


print(f"BKS alpha-quartz: a = {atoms.cell.cellpar()[0]:.4f} A, c = {atoms.cell.cellpar()[2]:.4f} A\n")
ionic = minelastic.calculate(atoms, BKS(), ionic=True, dipole=dipole, delta=0.005, displacement=0.005)
relaxed = minelastic.calculate(atoms, BKS(), relax_ions=True, fmax=1e-7, steps=4000, dipole=dipole, delta=0.005)

print("relaxed-ion constants (GPa)   (1) ionic=True   (2) relax_ions=True")
for name, (i, j) in {"c11": (0, 0), "c33": (2, 2), "c44": (3, 3), "c66": (5, 5), "c12": (0, 1), "c13": (0, 2),
                     "c14": (0, 3)}.items():
    print(f"  {name}                          {ionic.C_relaxed[i, j]:9.2f}   {relaxed.C[i, j]:9.2f}")
print("piezoelectric tensor (C/m^2)")
for name, (i, j) in {"e11": (0, 0), "e14": (0, 3)}.items():
    print(f"  {name}                          {ionic.piezo_total[i, j]:9.4f}   {relaxed.piezo[i, j]:9.4f}")
print(f"\nBorn charges of Si (diagonal, |e|): {np.diag(ionic.born[0]).round(4)}; clamped-ion e0 = "
      f"{np.abs(ionic.piezo).max():.4f} C/m^2")
print(f"d11 = {ionic.piezo_d[0, 0]:.3f} pC/N, d14 = {ionic.piezo_d[0, 3]:.3f} pC/N")
print(f"largest differences of the routes: {np.abs(ionic.C_relaxed - relaxed.C).max():.3f} GPa, "
      f"{np.abs(ionic.piezo_total - relaxed.piezo).max():.1e} C/m^2")
