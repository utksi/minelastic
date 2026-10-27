# Your own workflow

The command line covers the usual case: write the inputs, run them, read the outputs. The Python API lets
you replace any of these steps, for example to run the calculations with a workflow manager, with a code
that minelastic has no interface for, or with stresses that you already have.

## Bring your own stresses

`ElasticSetup` lists the strain states of a structure and builds their cells, and `compute` turns the
stresses into the elastic tensor. How the stresses are computed is up to you:

<!-- run -->
```python
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.units import GPa
from minelastic import ElasticSetup


def my_code_stress(atoms):
    """Stand-in for your own calculation: Cauchy stress in GPa, tension positive, Voigt xx yy zz yz xz xy."""
    atoms.calc = EMT()
    return atoms.get_stress(voigt=True) / GPa


setup = ElasticSetup(bulk("Cu", "fcc", a=3.5898, cubic=True))
stresses = {}
for state in setup.states:                         # state["id"] = 's1_m2', 's1_m1', ...
    cell = setup.strained_atoms(state)             # ASE Atoms of this strained cell
    stresses[state["id"]] = my_code_stress(cell)
result = setup.compute(stresses)
print({name: round(value, 2) for name, value in result.as_dict().items()})
```

<!-- output -->
```text
{'c11': 172.62, 'c12': 115.44, 'c44': 89.91}
```

Compute the stress of each cell exactly as `strained_atoms` returns it: the cell may have been rotated into
the standard orientation of the crystal ([Conventions](../reference/conventions.md#orientation)). `compute`
also takes the forces of the cells (`forces={state_id: (N, 3) array}` in eV/Å) for the internal-strain
tensor, and their dipoles (`polarization={state_id: (3,) array}` in e Å) for the piezoelectric tensor.

[`examples/replay_vasp_aln.py`](https://github.com/utksi/minelastic/tree/main/examples/replay_vasp_aln.py)
does this with the stresses and Berry-phase dipoles that VASP computed for the 4 strained cells of wurtzite
AlN, stored in `tests/data_berry/`:

<!-- excerpt: examples/replay_vasp_aln.py -->
```python
setup = ElasticSetup(atoms, delta=0.01, stencil="central")
res = setup.compute(data["stress_GPa"], polarization=data["dipole_eA"])
print(res.summary())
```

<!-- output-of: examples/replay_vasp_aln.py -->
```text
s1_m1 applied strain (Voigt, engineering shear): [-0.01 -0.   -0.   -0.01 -0.   -0.  ]
s1_p1 applied strain (Voigt, engineering shear): [0.01 0.   0.   0.01 0.   0.  ]
s2_m1 applied strain (Voigt, engineering shear): [-0.   -0.   -0.01 -0.   -0.   -0.  ]
s2_p1 applied strain (Voigt, engineering shear): [0.   0.   0.01 0.   0.   0.  ]

minelastic 0.1.0: P6_3mc (186), hexagonal, 2 strain pattern(s), method 'uncoupled'
  c11 =     440.26 GPa
...
piezoelectric tensor e_iJ (C/m^2, rows x y z, columns xx yy zz yz xz xy, working frame):
    -0.0000    0.0000    0.0000    0.0000    0.3351   -0.0000
    -0.0000    0.0000    0.0000    0.3351   -0.0000    0.0000
     0.2668    0.2668   -0.4285   -0.0000    0.0000   -0.0000
...

clamped-ion piezoelectric tensor (C/m^2): this run (Berry phase, 4 cells) vs DFPT (LEPSILON)
  e33 =  -0.428    -0.424
  e31 =   0.267     0.270
  e15 =   0.335     0.340
```

## Save and restore a setup

When the calculations run for days, store the setup and rebuild it later. `to_dict` returns a JSON-ready
record, the one that `generate` writes to `minelastic.json`, and `from_dict` rebuilds the same strain states
from it:

<!-- run -->
```python
import json

with open("setup.json", "w") as f:
    json.dump(setup.to_dict(), f)
with open("setup.json") as f:
    restored = ElasticSetup.from_dict(json.load(f))
print(restored.compute(stresses).as_dict() == result.as_dict())
```

<!-- output -->
```text
True
```

## Drive a DFT code from Python

`minelastic.generate` and `minelastic.collect` are the functions behind the two commands, so a VASP, Quantum
ESPRESSO or ABINIT workflow can be scripted without the command line. With the template of the
[VASP tutorial](../get-started/vasp.md):

<!-- run -->
```python
import minelastic

setup, notes = minelastic.generate("vasp", "vasp_aln", "el", stencil="central")
print(len(setup.states), "calculations")
for note in notes:
    print("NOTE:", note)
```

<!-- output -->
```text
4 calculations
NOTE: INCAR: set IBRION = -1, NSW = 0, ISIF = 2
```

After the runs, `collect` reads the results and writes `el/elastic.json`:

<!-- run -->
```python
result = minelastic.collect("el")
print({name: round(value, 2) for name, value in result.as_dict().items()})
```

<!-- output -->
```text
{'c11': 440.26, 'c12': 96.78, 'c13': 62.16, 'c33': 432.32, 'c44': 125.02}
```

The keyword arguments are the options of the commands ([Python API](../reference/api.md)).
