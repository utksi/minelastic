# Python API

<!-- input -->
```python
import minelastic
from minelastic import ElasticSetup, ElasticResult, calculate, generate, collect
```

All quantities follow the [conventions](conventions.md): GPa, Voigt order xx yy zz yz xz xy, engineering
shear strain, stress positive in tension.

| name | what it does |
|---|---|
| [`calculate`](#minelastic.calculate) | the whole calculation with an ASE calculator, in memory |
| [`ElasticSetup`](#minelastic.ElasticSetup) | the strain set of a structure, its strained cells, and the tensor from their stresses |
| [`ElasticResult`](#minelastic.ElasticResult) | the constants, tensors, moduli and checks |
| [`generate`](#minelastic.generate), [`collect`](#minelastic.collect) | the `minelastic generate` and `minelastic collect` commands, for DFT codes |

::: minelastic.calculate

::: minelastic.ElasticSetup
    options:
      members:
        - compute
        - strained_atoms
        - describe
        - relation_table
        - to_dict
        - from_dict

::: minelastic.ElasticResult
    options:
      members:
        - as_dict
        - summary
        - to_dict
        - C_input_frame
        - C_relaxed
        - properties

::: minelastic.generate

::: minelastic.collect

## Code interfaces and utilities

| module | contents |
|---|---|
| `minelastic.io.vasp` | `read_vasprun`, `read_outcar`, `read_result`, `read_elastic_tensor`, `read_dipole`, `read_born_charges`, `read_clamped_piezo`, `read_incar`, `kspacing_grid` |
| `minelastic.io.qe` | `read_structure` (any `ibrav`), `read_output`, `qe_lattice` |
| `minelastic.io.abinit`, `minelastic.io.aims`, `minelastic.io.generic` | `read_structure`, `prepare`, `read_result` |
| `minelastic.properties` | `vrh(C)`, `born_stable(C)` |
| `minelastic.voigt` | `rotate_stiffness(C, R)`, `mandel_rotation(R)`, conversions between strain and stress vectors and tensors |
| `minelastic.symmetry` | `analyze(atoms, symprec)`, `elastic_type(spacegroup_number)`, `rotate_atoms(atoms, R)` |
| `minelastic.tensor_form` | `TensorForm(rotations)`, `invariant_strain_dimension(rotations)` |
| `minelastic.algorithm` | `paper_greedy(form)`, `exact_minimum(form)`, `relations(form, pattern)`, `choose_amplitudes` |
| `minelastic.ionic` | `irreducible_displacements`, `force_constants`, `born_charges`, `ionic_contributions`, `strain_coefficients`, `rotate_born`, `rotate_piezo` |
| `minelastic.internal_strain`, `minelastic.piezo` | symmetry completion of the internal-strain and piezoelectric tensors |
