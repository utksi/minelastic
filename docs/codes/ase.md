# ASE calculators

Any ASE calculator that returns the stress can be used directly from Python: EMT and other interatomic
potentials, machine-learned potentials such as MACE, and DFT codes driven through ASE such as CP2K. The
[first tutorial](../get-started/first-tensor.md) uses EMT, and the [CP2K tutorial](../get-started/cp2k.md) a
DFT code.

## `calculate`

<!-- input -->
```python
result = minelastic.calculate(atoms, calc, relax_ions=False, method="uncoupled", fmax=1e-4, steps=500,
                              optimizer=None, dipole=None, **setup_kwargs)
```

builds the strain set of `atoms` (`ElasticSetup(atoms, **setup_kwargs)`), computes every strained and
displaced cell with the ASE calculator `calc` and returns the `ElasticResult`. For each cell it takes

* the stress, `calc.get_stress()` in eV/Å$^3$, converted to GPa (tension positive, as in ASE);
* the forces, for the internal-strain tensor and, with `ionic=True`, the force constants;
* the dipole `dipole(atoms)` in e Å, if a function is given, for the Born charges and the piezoelectric
  tensors.

With `relax_ions=True` the atoms of each strained cell are first relaxed at fixed cell with `optimizer`
(default BFGS) to `fmax` (eV/Å). `setup_kwargs` are the options of `ElasticSetup`: `scheme`, `delta`,
`stencil`, `symprec`, `frame`, `ionic` and `displacement` ([Python API](../reference/api.md)).

The how-to pages show each quantity with an ASE calculator:

| quantity | example | page |
|---|---|---|
| clamped- and relaxed-ion constants | fcc and hcp Cu, EMT | [Your first elastic tensor](../get-started/first-tensor.md) |
| relaxed-ion constants, two routes | hcp Cu, EMT (`examples/emt_ionic.py`) | [Relaxed-ion constants](../how-to/relaxed-ions.md#with-an-ase-calculator) |
| internal-strain tensor, force constants | hcp Cu, EMT (`examples/emt_internal_strain.py`) | [Internal strain](../how-to/internal-strain.md#with-an-ase-calculator) |
| Born charges, piezoelectric e and d | alpha-quartz, BKS rigid-ion model (`examples/bks_quartz_piezo.py`) | [Piezoelectricity](../how-to/piezoelectricity.md#with-an-ase-calculator) |
| finite temperature | fcc Al, EMT molecular dynamics (`examples/emt_finite_temperature.py`) | [Finite temperature](../how-to/finite-temperature.md#with-an-ase-calculator) |

Relax the reference cell with the same calculator first: the textbook constants refer to a stress-free
reference.

## Machine-learned potentials

Any ASE calculator works. Relax the reference cell with the same calculator, and use double precision: the
stresses of the strained cells differ by only about 1 GPa.

<!-- skip: needs MACE -->
```python
from ase.io import read
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from mace.calculators import mace_mp
import minelastic

atoms = read("POSCAR")
calc = mace_mp(model="medium", default_dtype="float64")
atoms.calc = calc
BFGS(FrechetCellFilter(atoms)).run(fmax=1e-4)           # stress-free reference
result = minelastic.calculate(atoms, calc, relax_ions=True)
print(result.summary())
```

For wurtzite AlN with MACE-MPA-0, the minimal set and the six single strains agree to 0.1 GPa, and the two
routes to the relaxed-ion constants to 0.3 GPa ([Validation](../reference/validation.md)).
