# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Internal-strain tensor of hcp Cu with ASE's EMT potential, checked against brute force (seconds).

The internal-strain tensor L[n, a, J] = dF_na/de_J is the force on atom n, component a, per unit strain J
(Voigt, engineering shear) at fixed fractional coordinates. minelastic obtains it twice:
1. from the forces in the strained cells of the minimal set (for hcp the two patterns xx+yz and zz),
   completed for all six strain components by symmetry;
2. with ionic=True, from the stresses of the displaced atoms, through Omega dsigma_J/du_na = -L[n, a, J].
The brute force below applies each of the six single-component strains and differentiates the forces of all
atoms numerically (four-point stencil). The ionic relaxation per unit strain, u = Phi^+ L with the Gamma-point
force constants Phi, is checked by relaxing the atoms in a strained cell.
"""
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize import BFGS

import minelastic
from minelastic.voigt import deformation_gradient

res = minelastic.calculate(bulk("Cu", "hcp", a=2.55, c=4.16), EMT(), ionic=True)
ref = res.setup.atoms                    # the structure in the working frame of the result


def strained(e):
    """Copy of ref with lattice vectors a_i' = F a_i, F = I + eps(e) (Voigt e), fractional coordinates kept."""
    at = ref.copy()
    at.set_cell(np.asarray(ref.cell[:]) @ deformation_gradient(e).T, scale_atoms=True)
    at.calc = EMT()
    return at


h = 0.005
brute = np.zeros((len(ref), 3, 6))
for J in range(6):
    F = {m: strained(m * h * np.eye(6)[J]).get_forces() for m in (-2, -1, 1, 2)}
    brute[:, :, J] = (F[-2] - 8 * F[-1] + 8 * F[1] - F[2]) / (12 * h)

print("internal-strain tensor of atom 1 (eV/A per unit strain), rows x y z, columns xx yy zz yz xz xy:")
print(res.internal_strain[0].round(4))
print(f"largest difference to brute force, all atoms: {np.abs(res.internal_strain - brute).max():.1e} eV/A "
      f"from the strained cells, {np.abs(res.internal_strain_displaced - brute).max():.1e} eV/A from the "
      "displaced atoms")

# u_n = sum_J (Phi^+ L)[n, :, J] e_J; check for e_xx = +-h by relaxing the atoms at fixed cell
u = []
for m in (1, -1):
    at = strained(m * h * np.eye(6)[0])
    x0 = at.get_positions()
    BFGS(at, logfile=None).run(fmax=1e-8, steps=1000)
    du = at.get_positions() - x0
    u.append(du - du.mean(axis=0))          # without a rigid translation
direct = (u[0] - u[1]) / (2 * h)
print("relaxation of atom 1 per unit strain xx (A):", res.ionic_relaxation[0, :, 0].round(4),
      "; relaxing the atoms:", direct[0].round(4))
