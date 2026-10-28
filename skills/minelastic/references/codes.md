# minelastic with each code

What `minelastic generate` writes for each code and what `minelastic collect` reads. Every output
directory also gets `STATES` (the cell directories, one per line) and `minelastic.json` (the setup, which
`collect` needs).

## VASP

Template: a folder with `POSCAR`, `INCAR`, `POTCAR` and `KPOINTS` (or `KSPACING` in the INCAR), for the
relaxed reference structure. Each cell directory gets the strained or displaced `POSCAR` (with the species
lines and selective-dynamics flags of the template), `INCAR`, `KPOINTS` and `POTCAR` (`--potcar link`
symlinks it).

`generate` copies the INCAR and changes only what the type of calculation requires, and lists every change:

* `IBRION = -1`, `NSW = 0`, `ISIF = 2` for clamped ions (`ISIF = 2` computes the stress);
* `IBRION = 2`, `NSW = 99`, `ISIF = 2` and, unless set, `EDIFFG = -5E-4` with `--relax-ions`;
* with `--md`, the template's `IBRION = 0`, `NSW` and `POTIM` are kept and `ISIF = 2` is set;
* `NFREE`, `POTIM`, `LEPSILON` and `LCALCEPS` are removed, so an `IBRION = 6` input works as it is;
* `KSPACING` is replaced by the explicit k-point mesh of the reference cell, so all cells use one mesh.

It also checks, and reports without changing them unless `--fix-incar` is given: `EDIFF` at most 1E-6
(1E-8 suggested; 1E-6 for MD), `ENCUT` at least 1.3 times the largest `ENMAX` of the POTCAR,
`PREC = Accurate` (or `High`) and `LREAL = .FALSE.`. `--fix-incar` applies all of them or none.

`collect` reads `vasprun.xml`, then `OUTCAR` (`--output` for another name). For piezoelectric tensors, put
`LCALCPOL = .TRUE.` in the template INCAR: `collect --piezo` reads the Berry-phase dipole of every cell.
`--born OUTCAR` takes the Born charges and the clamped-ion piezoelectric tensor of a `LEPSILON` or
`LCALCEPS` run of the reference instead. `minelastic vasp-tensor OUTCAR` converts the tensor of an
`IBRION = 6, ISIF >= 3` run to GPa (`--which total|clamped|ionic`).

Magnetic moments: `MAGMOM` with `ISPIN = 2` (collinear) or with `LNONCOLLINEAR`/`LSORBIT` (vectors).

## Quantum ESPRESSO (pw.x)

Template: the pw.x input file of the relaxed reference. Set `pseudo_dir` to an absolute path, because
pw.x runs in one sub-directory per cell. Every strained input is written with `ibrav = 0`, the strained
`CELL_PARAMETERS` in Å, `ATOMIC_POSITIONS crystal`, the species labels and `if_pos` flags of the
template, and `tstress` and `tprnfor` switched on; all other namelists and cards are copied. Any `ibrav`
is accepted. With `--relax-ions`, `calculation = 'relax'` is set, with `forc_conv_thr = 2.0d-5` and
`etot_conv_thr = 1.0d-8` unless the template sets them (the pw.x defaults stop too early). With `--md`, the
template has `calculation = 'md'` with a thermostat in `&IONS` and a fixed cell.

<!-- skip: needs pw.x -->
```bash
minelastic generate qe pw.in -o el
for d in $(cat el/STATES); do (cd el/$d && mpirun pw.x -in pw.in > pw.out); done
minelastic collect el
```

`collect` reads `pw.out` in each directory (`--output` for another name). pw.x prints the stress in
Ry/bohr^3, positive in compression; minelastic converts it. Magnetic sublattices need separate species
labels (one `starting_magnetization` each), and minelastic treats the labels as distinct species; with
`noncolin = .true.` it also reads `angle1` and `angle2`. A collinear antiferromagnet can then get a larger
set than its magnetic group needs (MnF2: 3 patterns instead of 2); the constants are the same.

## ABINIT

Template: a folder with `POSCAR` and one input `*.abi` that reads it with `structure "poscar:POSCAR"`. The
input sets everything else (pseudopotentials in the order of the POSCAR species, cutoff, k-points,
convergence). The cell must stay fixed (`optcell 0`) and, for clamped ions, the ions too (no `ionmov`).
`--relax-ions` appends `ionmov 2`, `ntime 100` and `tolmxf 1.0d-6` unless `ionmov` is set. `collect` reads
the newest `*.abo` file of each directory. Moments: `spinat` (z components with `nsppol 2` or `nspden 2`,
vectors with `nspden 4`).

## FHI-aims

Template: a folder with `geometry.in` and `control.in`. `control.in` must request the stress
(`compute_analytical_stress .true.`) and must not relax the cell; for clamped ions it must not relax the
atoms. `--relax-ions` adds `relax_geometry bfgs 1E-3`. `collect` reads `aims.out`. Moments: `initial_moment`
in `geometry.in`.

## Any code ASE can write (CASTEP, SIESTA, ...)

`generic` writes only the strained structure and copies further files; ASE reads the results back. The
inputs must request the stress and keep the cell fixed.

<!-- skip: needs CASTEP -->
```bash
minelastic info crystal.cell --format castep-cell
minelastic generate generic crystal.cell -o el --format castep-cell --filename crystal.cell --extra crystal.param
minelastic collect el --output crystal.castep --format castep-castep
```

## CP2K and other codes through ASE

Use `minelastic.calculate` with the ASE calculator (for CP2K, `ase.calculators.cp2k.CP2K` with a
`cp2k_shell` command). ASE's CP2K calculator has no per-atom magnetic moments.

## Running the cells

Any scheduler works; a Slurm job array over the lines of `STATES` is the usual way:

<!-- skip: needs Slurm and VASP -->
```bash
#!/bin/bash
#SBATCH --array=1-8              # one task per line of el/STATES
d=el/$(sed -n "${SLURM_ARRAY_TASK_ID}p" el/STATES)
cd "$d" && srun vasp_std > out
```
