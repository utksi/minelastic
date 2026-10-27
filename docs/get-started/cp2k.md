# First run with CP2K

CP2K is driven through ASE's CP2K calculator, so the whole calculation is one Python script, with no
input templates or output files to handle. Here we compute the clamped- and relaxed-ion elastic constants
of silicon in a 64-atom cubic supercell at the $\Gamma$ point. The output shown below is from this run.

## 1. Connect ASE to CP2K

ASE's calculator talks to CP2K's interactive shell, which current CP2K versions start with
`cp2k.psmp --shell`. ASE requires the string `cp2k_shell` in the command, so the shell is started from a
small script with that name. Every MPI rank reads the shell commands, so srun must pass the standard input
to all tasks, which is its default:

<!-- skip: shell script for CP2K -->
```bash
#!/bin/bash
# file: cp2k_shell
exec srun -u -n "$SLURM_NTASKS" cp2k.psmp --shell "$@"
```

Make it executable (`chmod +x cp2k_shell`) and run the Python script below inside a Slurm job.

## 2. The script

<!-- skip: needs CP2K -->
```python
import os

from ase.build import bulk
from ase.calculators.cp2k import CP2K
from ase.units import Ry

import minelastic

data = "/path/to/cp2k/data"           # the directory with CP2K's BASIS_MOLOPT and GTH_POTENTIALS
inp = """&FORCE_EVAL
  &DFT
    &SCF
      EPS_SCF 1.0E-8
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
        ENERGY_GAP 0.05
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-8
        MAX_SCF 20
      &END OUTER_SCF
    &END SCF
    &MGRID
      REL_CUTOFF 60
    &END MGRID
  &END DFT
&END FORCE_EVAL
"""
calc = CP2K(command=os.path.abspath("cp2k_shell"), xc="PBE", cutoff=500 * Ry, max_scf=40,
            basis_set="DZVP-MOLOPT-SR-GTH", pseudo_potential="GTH-PBE",
            basis_set_file=f"{data}/BASIS_MOLOPT", potential_file=f"{data}/GTH_POTENTIALS",
            stress_tensor=True, inp=inp)

si = bulk("Si", "diamond", a=5.468, cubic=True).repeat(2)          # 64 atoms
result = minelastic.calculate(si, calc, ionic=True, stencil="central")
print(result.summary())
```

The calculator uses PBE with GTH pseudopotentials, the DZVP-MOLOPT-SR-GTH basis, a density cutoff of 500 Ry
and a tight SCF threshold; `stress_tensor=True` makes CP2K compute the stress. This run uses $a$ = 5.468 Å;
relax the lattice parameter with the same calculator first.

!!! info "Convergence"

    Elastic constants are obtained from derivatives of the stress with respect to the strain, which are estimated
    from stress differences of a few GPa. They are therefore much more sensitive to the density cutoff, the basis
    set, the supercell size (or k-point mesh) and the SCF threshold than total energies are. The settings here
    are chosen to show the workflow; they are not the result of a convergence study. Before using the constants
    in a publication, check that they change by less than the accuracy you need when the cutoff is raised and the
    basis and supercell are enlarged, and relax the reference structure with the final settings.


`calculate` finds the symmetry (cubic, one strain pattern), builds every strained cell, computes its stress
with CP2K and returns the result. `ionic=True` adds cells in which one atom is displaced, which give the
relaxed-ion constants without relaxing anything. `stencil="central"` applies each strain and displacement
at two amplitudes instead of the default four, which halves the number of CP2K calculations.

## 3. The result

<!-- file: tests/data_docs/si_cp2k_summary.txt -->
```text
minelastic 0.1.0: Fd-3m (227), cubic, 1 strain pattern(s), method 'uncoupled'
  c11 =     156.50 GPa
  c12 =      60.56 GPa
  c44 =     102.57 GPa
stiffness matrix (GPa, Voigt xx yy zz yz xz xy, working frame):
     156.50     60.56     60.56      0.00      0.00      0.00
      60.56    156.50     60.56      0.00      0.00      0.00
      60.56     60.56    156.50      0.00      0.00      0.00
       0.00      0.00      0.00    102.57      0.00      0.00
       0.00      0.00      0.00      0.00    102.57      0.00
       0.00      0.00      0.00      0.00      0.00    102.57
Born stable: True (eigenvalues 95.9, 95.9, 102.6, 102.6, 102.6, 277.6)
K_VRH = 92.54  G_VRH = 75.60  E = 178.26 GPa  nu = 0.179  A_U = 0.727
ionic relaxation contribution to C (GPa):
      -0.00     -0.00      0.00      0.00      0.00      0.00
      -0.00     -0.00     -0.00     -0.00     -0.00     -0.00
       0.00     -0.00     -0.00     -0.00      0.00     -0.00
       0.00     -0.00     -0.00    -22.85      0.00      0.00
       0.00     -0.00      0.00      0.00    -22.85      0.00
       0.00     -0.00     -0.00      0.00      0.00    -22.85
relaxed-ion stiffness matrix (GPa):
     156.50     60.56     60.56      0.00      0.00      0.00
      60.56    156.50     60.56     -0.00     -0.00     -0.00
      60.56     60.56    156.50     -0.00      0.00     -0.00
       0.00     -0.00     -0.00     79.72      0.00      0.00
       0.00     -0.00      0.00      0.00     79.72      0.00
       0.00     -0.00     -0.00      0.00      0.00     79.72
relaxed ions: Born stable True (eigenvalues 79.7, 79.7, 79.7, 95.9, 95.9, 277.6)
consistency (redundant stress components): largest |computed - fitted| = 0.00 GPa (pattern 1, xz)
...
```

The first block holds the clamped-ion constants, with the atoms kept at their strained positions. In
diamond, a shear strain moves the two atoms of the primitive cell against each other, which lowers c44 from
102.57 GPa (clamped ions) to 79.72 GPa (relaxed ions); c11 and c12 do not change. `result.as_dict()` returns
the clamped-ion constants, `result.C_relaxed` the relaxed-ion matrix.

## 4. Check the result

The same script with `scheme="asess"` in place of `ionic=True` applies the six single strains, 12 cells
instead of 2, and gives

<!-- file: tests/data_docs/si_cp2k_summary.txt -->
```text
...
minelastic 0.1.0: Fd-3m (227), cubic, 6 strain pattern(s), method 'uncoupled'
  c11 =     156.47 GPa
  c12 =      60.53 GPa
  c44 =     102.52 GPa
...
```

which agrees with the single pattern to 0.05 GPa.
