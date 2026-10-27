# Internal strain and force constants

The internal-strain tensor is the force on each atom per unit strain, with the atoms held at their strained
fractional positions,

$$
\Lambda_{n\alpha, J} = \frac{\partial F_{n\alpha}}{\partial \varepsilon_J}
$$

for atom $n$, Cartesian direction $\alpha$ and strain component $J$, in eV/Å. It describes how the atoms
move when the crystal is strained, and it enters the relaxed-ion elastic constants and the piezoelectric
tensor ([theory](../theory/internal-strain.md)).

minelastic obtains $\Lambda$ in two ways. The forces in the strained cells give it without further calculations: by
symmetry, the few patterns of the minimal set determine the whole tensor. With `--ionic`, the stresses of
cells with displaced atoms give it a second time, together with the $\Gamma$-point force constants
$\Phi_{m\alpha, n\beta} = -\partial F_{m\alpha}/\partial u_{n\beta}$ in eV/Å$^2$.

## With an ASE calculator

[`examples/emt_internal_strain.py`](https://github.com/utksi/minelastic/tree/main/examples/emt_internal_strain.py)
computes both routes for hcp Cu with EMT. `ionic=True` adds the cells with displaced atoms:

<!-- run -->
```python
from ase.build import bulk
from ase.calculators.emt import EMT
import minelastic

res = minelastic.calculate(bulk("Cu", "hcp", a=2.55, c=4.16), EMT(), ionic=True)
print(res.internal_strain.shape)              # Lambda from the strained cells, (N, 3, 6), eV/A
print(res.internal_strain_displaced.shape)    # Lambda from the displaced atoms
print(res.force_constants.shape)              # Phi, (N, 3, N, 3), eV/A^2
print(res.ionic_relaxation.shape)             # Phi^+ Lambda, (N, 3, 6), A per unit strain
```

<!-- output -->
```text
(2, 3, 6)
(2, 3, 6)
(2, 3, 2, 3)
(2, 3, 6)
```

The script compares both tensors with brute force, the forces of all atoms under each of the six single
strains:

<!-- output-of: examples/emt_internal_strain.py -->
```text
...
internal-strain tensor of atom 1 (eV/A per unit strain), rows x y z, columns xx yy zz yz xz xy:
[[ 0.      0.     -0.      0.      0.     -1.9159]
 [-1.9159  1.9159 -0.     -0.      0.     -0.    ]
 [ 0.     -0.     -0.      0.     -0.     -0.    ]]
largest difference to brute force, all atoms: 3.5e-05 eV/A from the strained cells, 4.2e-06 eV/A from the displaced atoms
relaxation of atom 1 per unit strain xx (A): [ 0.     -0.6114 -0.    ] ; relaxing the atoms: [ 0.     -0.6114 -0.    ]
```

The two patterns of the minimal set determine the whole tensor, and the predicted relaxation of the atoms
per unit strain, $u = \Phi^+ \Lambda$, is the one that relaxing the atoms gives.

## With VASP

The forces are read together with the stresses, so `collect --internal-strain` prints the tensor of any
clamped-ion run, here that of the [VASP tutorial](../get-started/vasp.md):

<!-- shell -->
```bash
minelastic generate vasp vasp_aln -o el --stencil central
for d in $(cat el/STATES); do (cd el/$d && mpirun vasp_std > out); done
```

<!-- shell -->
```console
$ minelastic collect el --internal-strain
...
internal strain tensor dF/de (eV/A), columns xx yy zz yz xz xy:
  ion   1 x    -0.0000     0.0000     0.0000    -0.0000    -4.9075    -8.2521
  ion   1 y    -8.2521     8.2521     0.0000    -4.9075    -0.0000     0.0000
  ion   1 z    -5.0802    -5.0802    11.0955     0.0000    -0.0000    -0.0000
...
```

The tensor has the sign and units of the one that VASP prints for `IBRION = 6` (`INTERNAL STRAIN TENSOR FOR
ION`), which VASP computes from displaced atoms:

| Al (ion 1), eV/Å | z, zz | x, xy | x, xz | z, xx |
|---|---|---|---|---|
| minelastic, 4 strained cells | 11.095 | -8.252 | -4.908 | -5.080 |
| VASP `IBRION = 6` | 11.130 | -8.243 | -4.897 | -5.084 |

The largest difference over all atoms and components is 0.04 eV/Å. With `--ionic`, the tensor from the
displaced atoms agrees with VASP's to 0.06 eV/Å, and the force constants are written to `elastic.json`
(`force_constants_eV_per_A2`).

## With Quantum ESPRESSO and other codes

pw.x forces are read from the same output as the stresses (`generate` switches `tprnfor` on), and the
[Quantum ESPRESSO](../get-started/qe.md) and [ABINIT](../get-started/abinit.md) tutorials end with the
tensor. For wurtzite AlN with pw.x, the tensors from the minimal set and from the six single strains agree
to 0.007 eV/Å. Every interface that returns forces gives the tensor; with your own calculations, pass the
forces as `compute(stresses, forces={state_id: (N, 3) array})` ([Your own workflow](own-workflow.md)).
