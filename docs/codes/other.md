# ABINIT, FHI-aims and other codes

## ABINIT

A complete example is the [ABINIT tutorial](../get-started/abinit.md). The template
directory holds the reference structure as `POSCAR` and one ABINIT input file (`*.abi`) that reads it with
`structure "poscar:POSCAR"`. The input sets everything else: pseudopotentials in the order of the POSCAR
species, cutoff, k-point mesh and convergence.

ABINIT's POSCAR reader takes the element of each atom from the end of its coordinate line, so minelastic
writes the strained POSCARs in that form; for the same reason two species of one element cannot be used.
Magnetic sublattices are set with `spinat`, which minelastic also reads for the symmetry
([Magnetic materials](../how-to/magnetic.md)). `collect` reads the stress and forces from the newest
`*.abo` file of each directory. With `--relax-ions`, `ionmov 2`, `ntime 100` and `tolmxf 1.0d-6` are
appended unless the input sets `ionmov`. The cell must stay fixed (`optcell 0`), the ions must stay fixed for
clamped-ion runs (no `ionmov`), and dataset-specific values of `ionmov` and `optcell` are not accepted.

## FHI-aims

The template directory holds `geometry.in` (the reference structure) and `control.in`:

<!-- skip: needs an FHI-aims template and the FHI-aims runs -->
```bash
minelastic generate aims ref/ -o el          # geometry.in per strain state, control.in copied
# run aims in every directory, output in aims.out
minelastic collect el
```

`control.in` must request the stress (`compute_analytical_stress .true.`) and must not relax the cell
(`relax_unit_cell`); otherwise `generate` stops with an error. For clamped ions it must not relax the atoms
either. With `--relax-ions`, `relax_geometry bfgs 1E-3` is added unless it is present. The FHI-aims reader
is tested on an output distributed with ASE; minelastic has not been run with FHI-aims itself.

## Any format ASE can write (CASTEP, SIESTA, ...)

The `generic` interface writes only the strained structure in each directory and copies further input
files. The results are read back with ASE:

<!-- skip: needs CASTEP inputs and the CASTEP runs -->
```bash
minelastic generate generic crystal.cell -o el --format castep-cell --filename crystal.cell \
    --extra crystal.param
# run the code in every directory
minelastic collect el --output crystal.castep --format castep-castep
```

The reference structure (`crystal.cell` here) is read with the same `--format` as the files that are
written. Make sure the inputs request the stress and keep the cell fixed. Keep the ions fixed too, unless
you want relaxed ions. `minelastic info crystal.cell --format castep-cell` shows the strain set first.

## Anything else

The Python API separates the strain set from how the calculations are run, so stresses from any source can
be used ([Your own workflow](../how-to/own-workflow.md)).
