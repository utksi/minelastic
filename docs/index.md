<h1 class="wordmark-title"><img src="assets/wordmark-light.svg#only-light" alt="minelastic"><img src="assets/wordmark-dark.svg#only-dark" alt="minelastic"></h1>

Elastic, internal-strain and piezoelectric tensors from the fewest possible stress calculations, for every
crystal symmetry and with any electronic-structure code.

![The usual approach applies six strain patterns; minelastic applies one to a cubic crystal and reads c11, c12 and c44 from the stress components of the strained cell. Below, the number of strain patterns for each crystal system.](assets/hero-light.svg#only-light)
![The usual approach applies six strain patterns; minelastic applies one to a cubic crystal and reads c11, c12 and c44 from the stress components of the strained cell. Below, the number of strain patterns for each crystal system.](assets/hero-dark.svg#only-dark)

**minelastic** implements the minimal and uncoupled strain sets scheme of S. I. Simak and U. Singh,
[arXiv:2610.11733](https://arxiv.org/abs/2610.11733) (2026). Each independent elastic constant is read off
one stress component of one strained cell, and the number of strain patterns is the smallest that any
stress-strain method can use: one for a cubic crystal and two for tetragonal, trigonal and hexagonal
crystals, where the usual approach applies six. The constants agree with those of the six single strains to
a few hundredths of a GPa.

## Example

<!-- run -->
```python
from ase.build import bulk
from ase.calculators.emt import EMT
import minelastic

result = minelastic.calculate(bulk("Cu", "fcc", a=3.59, cubic=True), EMT())
print(result.summary())
```

<!-- output -->
```text
minelastic 0.1.0: Fm-3m (225), cubic, 1 strain pattern(s), method 'uncoupled'
  c11 =     172.50 GPa
  c12 =     115.35 GPa
  c44 =      89.84 GPa
stiffness matrix (GPa, Voigt xx yy zz yz xz xy, working frame):
     172.50    115.35    115.35      0.00      0.00      0.00
     115.35    172.50    115.35      0.00      0.00      0.00
     115.35    115.35    172.50      0.00      0.00      0.00
       0.00      0.00      0.00     89.84      0.00      0.00
       0.00      0.00      0.00      0.00     89.84      0.00
       0.00      0.00      0.00      0.00      0.00     89.84
Born stable: True (eigenvalues 57.2, 57.2, 89.8, 89.8, 89.8, 403.2)
K_VRH = 134.40  G_VRH = 56.85  E = 149.48 GPa  nu = 0.315  A_U = 1.755
consistency (redundant stress components): largest |computed - fitted| = 0.00 GPa ...
```

Any ASE calculator works the same way. For DFT codes, minelastic writes the input files of the strained
cells and reads the finished calculations:

<!-- skip: needs a VASP template and the VASP runs -->
```bash
minelastic generate vasp ref/ -o el    # ref/: the VASP input of the relaxed structure
# run VASP in every directory listed in el/STATES
minelastic collect el
```

## What it computes

| quantity | calculations beyond the strained cells | guide |
|---|---|---|
| elastic constants, bulk and shear moduli, stability | none | [Elastic constants in practice](how-to/elastic.md) |
| internal-strain tensor | none: the forces of the same cells | [Internal strain](how-to/internal-strain.md) |
| clamped-ion piezoelectric tensor | none: Berry phase in the same cells | [Piezoelectricity](how-to/piezoelectricity.md) |
| relaxed-ion elastic constants | relaxations, or single points with displaced atoms | [Relaxed-ion constants](how-to/relaxed-ions.md) |
| Born charges, total piezoelectric tensors e and d | single points with displaced atoms | [Piezoelectricity](how-to/piezoelectricity.md) |
| elastic constants at finite temperature, with error bars | one MD run per strained cell | [Finite temperature](how-to/finite-temperature.md) |

The symmetry analysis includes magnetic order ([Magnetic materials](how-to/magnetic.md)). Every quantity is
built on the same strained cells, so each inherits the saving of the minimal set
([Cost](theory/cost.md)).

## Supported codes

minelastic writes the inputs for VASP and Quantum ESPRESSO from a template of the relaxed structure, and for
ABINIT, FHI-aims and every structure format that ASE can write (CASTEP, SIESTA, ...). Any ASE calculator,
such as EMT, MACE or CP2K, is used directly from Python. There are tutorials for
[VASP](get-started/vasp.md), [Quantum ESPRESSO](get-started/qe.md), [ABINIT](get-started/abinit.md) and
[CP2K](get-started/cp2k.md), with the output of real calculations.

## Usage with agentic workflows

An AI coding agent such as Claude Code can run minelastic for you. The repository includes an agent skill,
[`skills/minelastic`](https://github.com/utksi/minelastic/tree/main/skills/minelastic), with the commands and
the checks that a correct calculation needs. In Claude Code:

<!-- skip: Claude Code commands -->
```text
/plugin marketplace add utksi/minelastic
/plugin install minelastic@minelastic
```

Then ask in plain words, for example "compute the relaxed-ion elastic constants of this POSCAR with MACE".
Other agents and more examples: [Usage with agentic workflows](how-to/agentic-workflows.md).

## How this documentation is organized

If you are new to minelastic, start with the [installation](get-started/installation.md) and
[your first elastic tensor](get-started/first-tensor.md), which needs no DFT code, and then the tutorial for
your code. The [how-to guides](how-to/elastic.md) treat each quantity and each code in more depth. The
[theory](theory/index.md) section has the method and all formulas; you do not need it to use minelastic. The
[reference](reference/api.md) documents the Python API, the command line, units and conventions, and the
validation against VASP, Quantum ESPRESSO, ABINIT and CP2K.

If you use minelastic in published work, please cite the paper ([Citing](reference/citing.md)).
