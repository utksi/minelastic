# Magnetic materials

Magnetic order can lower the symmetry of a crystal, and with it the form of the elastic tensor and the
strain set. minelastic uses the magnetic space group (from spglib), including the operations combined with
time reversal, which also constrain the internal-strain and piezoelectric tensors. Rocksalt NiO, for
example, is cubic without magnetic order, but its AFM-I order, ferromagnetic (001) planes of Ni with
alternating moments, makes it tetragonal:

<!-- run -->
```python
from ase.build import bulk
import minelastic

nio = bulk("NiO", "rocksalt", a=4.17, cubic=True)
print(minelastic.ElasticSetup(nio).describe().splitlines()[0])

z = nio.get_scaled_positions()[:, 2]
moments = [(2.0 if zi < 0.25 else -2.0) if s == "Ni" else 0.0
           for s, zi in zip(nio.get_chemical_symbols(), z)]
nio.set_initial_magnetic_moments(moments)
setup = minelastic.ElasticSetup(nio)
print(setup.describe().splitlines()[0])
print(len(setup.patterns), "strain patterns:", setup.form.labels)
```

<!-- output -->
```text
space group Fm-3m (225), point group m-3m, elastic type cubic
space group P4/mmm (123), point group 4/mmm, elastic type tetragonal-II
2 strain patterns: ['c11', 'c12', 'c13', 'c33', 'c44', 'c66']
```

There are six constants and two strain patterns instead of three and one. The difference matters: in DFT
calculations of NiO AFM-I with VASP, pw.x and ABINIT, c11 and c33 differ by 33 to 40 GPa, which the cubic
form, used if the moments were ignored, would miss
([Validation](../reference/validation.md#end-to-end-runs-with-real-codes)).

## Where the moments come from

The initial moments are read from the template:

| code | moments used for the symmetry |
|---|---|
| VASP | `MAGMOM` when `ISPIN = 2`, `LNONCOLLINEAR` or `LSORBIT` is set (collinear as scalars, noncollinear as axial vectors) |
| ABINIT | `spinat`: the z components with `nsppol 2` or `nspden 2`, the vectors with `nspden 4` |
| pw.x | the species labels; with `noncolin`, also the directions from `starting_magnetization`, `angle1`, `angle2` |
| FHI-aims, ASE | the initial magnetic moments of the structure (`initial_moment` in geometry.in) |

Check the space group that `minelastic info` prints before running anything.

## Points to watch

A collinear pw.x run distinguishes the magnetic sublattices by their species labels, as pw.x requires, but
then the operations that exchange them together with a spin flip are not used. An antiferromagnet such as
MnF2 then gets the orthorhombic set of 3 patterns instead of the tetragonal set of 2; the constants are the
same, at the cost of more cells.

!!! warning "ASE calculators"

    With an ASE calculator, the moments of the structure must be those the calculator uses. EMT and MACE
    ignore them, and ASE's CP2K calculator has no per-atom moments.

If the cell has to be rotated to its standard orientation, noncollinear inputs are refused, because the
moment directions would have to be rotated too. Use `--frame input` (`frame="input"`) to keep the input
orientation.
