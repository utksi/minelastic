# Quantum ESPRESSO

A complete example is the [Quantum ESPRESSO tutorial](../get-started/qe.md). Below we describe how
minelastic handles the pw.x inputs and outputs.

## Template and commands

The template is the pw.x input of the relaxed reference structure. pw.x runs in one sub-directory per
state, so set `pseudo_dir` to an absolute path.

<!-- skip: needs a pw.x template and the pw.x runs -->
```bash
minelastic generate qe pw.in -o el                  # clamped ions
minelastic generate qe pw.in -o el_relaxed --relax-ions
for d in $(cat el/STATES); do (cd el/$d && mpirun pw.x -in pw.in > pw.out); done
minelastic collect el
```

## The inputs

Every strained input is written with

* `ibrav = 0`, the strained `CELL_PARAMETERS` in Å and `ATOMIC_POSITIONS crystal`;
* the species labels and `if_pos` flags of the template;
* `tstress` and `tprnfor` switched on.

Any `ibrav` is accepted, with `celldm` or with `A, B, C, cosAB, ...`, and the positions may be given in
`alat`, `bohr`, `angstrom` or `crystal`. All other namelists and cards are copied unchanged. For relaxed
ions, `calculation = 'relax'` is set, with tight `forc_conv_thr` and `etot_conv_thr` unless the template
sets them: the defaults of pw.x stop the relaxation before the internal strain is relaxed. `collect` reads
`pw.out` in each directory (`--output` for another name). pw.x prints the stress in Ry/bohr$^3$, positive in
compression; minelastic converts it to GPa, positive in tension.

## Magnetic order and molecular dynamics

pw.x needs a separate species label for every magnetic sublattice (starting magnetization per species), and
minelastic treats the labels as distinct species in the symmetry analysis. With `noncolin = .true.` it also
reads the moment directions from `starting_magnetization`, `angle1` and `angle2`
([Magnetic materials](../how-to/magnetic.md)). For finite temperature, the template has `calculation = 'md'`
with a thermostat in `&IONS` and a fixed cell; `generate --md` keeps it and `collect --md` averages the
stresses of all steps ([Finite temperature](../how-to/finite-temperature.md#with-vasp-or-quantum-espresso)).
