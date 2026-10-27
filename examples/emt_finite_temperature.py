# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Elastic constants of fcc Al at 300 K from molecular dynamics with ASE's EMT potential (about 3 minutes).

The finite-temperature protocol of the paper:
1. the lattice parameter at T from a zero-pressure NPT run;
2. fixed-cell NVT runs in the strained cells of the minimal set (cubic: one strain pattern, at +-delta);
3. C_ij(T) from the strain derivative of the ensemble-averaged stress.

The stress includes the kinetic (ideal-gas) part, get_stress(include_ideal_gas=True). For codes whose MD
stress contains only the virial (VASP), pass kinetic_temperature=T to compute() (or --kinetic-temperature to
minelastic collect) instead. Block averages are passed as the samples, so that the standard errors are those
of nearly uncorrelated data. Requires ASE >= 3.23.
"""
import numpy as np
from ase import units
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.filters import FrechetCellFilter
from ase.md.langevin import Langevin
from ase.md.nptberendsen import NPTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import BFGS

import minelastic

T = 300.0                                  # K
DT = 2 * units.fs
EQUIL, PROD, BLOCK = 1000, 6000, 250       # MD steps per strained cell: equilibration, production, block length
rng = np.random.default_rng(7)

# 0 K reference: relaxed cell and static elastic constants
ref = bulk("Al", "fcc", a=4.05, cubic=True)
ref.calc = EMT()
BFGS(FrechetCellFilter(ref), logfile=None).run(fmax=1e-7)
static = minelastic.calculate(ref, EMT())

# 1. lattice parameter at T: isotropic NPT (Berendsen) of a 3x3x3 supercell (108 atoms) at zero pressure
cell = ref.repeat(3)
cell.calc = EMT()
MaxwellBoltzmannDistribution(cell, temperature_K=T, rng=rng)
npt = NPTBerendsen(cell, DT, temperature_K=T, pressure_au=0.0, taut=100 * units.fs, taup=1000 * units.fs,
                   compressibility_au=1 / (75 * units.GPa))
npt.run(2000)
volumes = []
npt.attach(lambda: volumes.append(cell.get_volume()), interval=10)
npt.run(3000)
a_T = (np.mean(volumes) / 27) ** (1 / 3)
print(f"a(0 K) = {ref.cell[0, 0]:.4f} A, a({T:.0f} K) = {a_T:.4f} A")

# 2. NVT in every strained cell of the minimal set; the reference structure is the supercell at a(T)
# central differences: the four-point stencil would amplify the noise of the MD averages
setup = minelastic.ElasticSetup(bulk("Al", "fcc", a=a_T, cubic=True).repeat(3), delta=0.01, stencil="central")
stresses = {}
for st in setup.states:
    atoms = setup.strained_atoms(st)
    atoms.calc = EMT()
    MaxwellBoltzmannDistribution(atoms, temperature_K=T, rng=rng)
    dyn = Langevin(atoms, DT, temperature_K=T, friction=0.01 / units.fs, rng=rng)
    dyn.run(EQUIL)
    s = []
    dyn.attach(lambda: s.append(atoms.get_stress(voigt=True, include_ideal_gas=True) / units.GPa))
    dyn.run(PROD)
    s = np.array(s[:len(s) // BLOCK * BLOCK])   # GPa, tension positive, Voigt xx yy zz yz xz xy
    stresses[st["id"]] = s.reshape(-1, BLOCK, 6).mean(axis=1)
    print(f"{st['id']}: {len(s)} steps, <stress> = {np.round(s.mean(axis=0), 3)} GPa")

# 3. elastic constants with standard errors
res = setup.compute(stresses)
print(f"\n{'':6s} {'0 K':>8s} {T:8.0f} K")
for (k, c0), cT, err in zip(static.as_dict().items(), res.constants, res.errors):
    print(f"{k:6s} {c0:8.1f} {cT:8.1f} +- {err:.1f}  GPa")
