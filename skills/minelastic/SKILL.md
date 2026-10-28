---
name: minelastic
description: Elastic constants, relaxed-ion constants, internal-strain, Born-charge and piezoelectric tensors, and elastic constants at finite temperature, with the minelastic Python package and VASP, Quantum ESPRESSO, ABINIT, FHI-aims, other DFT codes or any ASE calculator (EMT, MACE, CP2K). Use when the user wants elastic constants, a stiffness tensor, bulk, shear or Young's moduli, piezoelectric coefficients or strained DFT inputs, or mentions minelastic.
---

# minelastic

minelastic computes the elastic tensor from the fewest strain patterns the crystal symmetry allows: 1 for
cubic crystals, 2 for hexagonal, tetragonal and trigonal ones, up to 6 for triclinic ones, with each
constant taken from one stress component. With an ASE calculator it runs everything in memory. With a DFT
code it writes the inputs of the strained cells (`generate`), the user runs them, and it reads the outputs
(`collect`). Units are GPa, Voigt order xx yy zz yz xz xy, stress positive in tension. Documentation:
https://utksi.github.io/minelastic/

## Rules for every task

1. **Start from a stress-free reference.** The constants are derivatives of the stress at the reference
   structure. With an ASE calculator, relax the cell and the positions with the same calculator before
   computing anything, even when the user did not ask, and say so. A DFT template must already be relaxed
   with the settings the strained cells will use; minelastic cannot check this, so state the assumption.
2. **One set of settings for every cell.** All cells of a set use the same code, functional, cutoff,
   k-point mesh, smearing and pseudopotentials as each other and as the reference relaxation. The constants
   come from stress differences of a few GPa, and a systematic offset between cells does not cancel. Do not
   set up reruns of part of a set with changed settings to combine with the old cells, even when asked to
   save time: explain why, and offer to regenerate and rerun the whole set (it has only a few cells).
3. **Converged settings.** Elastic constants need tighter settings than total energies. When the cutoff,
   the k-point mesh, the electronic convergence or the smearing look loose, name them and recommend
   converging the constants themselves before they are published. Do not change the user's settings
   silently: `generate vasp` reports its recommendations and `--fix-incar` applies them.
4. **Check the symmetry first.** Run `minelastic info` before `generate`. For a template, run it on the
   template (`info --code CODE TEMPLATE`), not on the bare structure file, so that magnetic moments are
   included: magnetic order can lower the symmetry, and then the larger set is correct. If a structure
   comes out with lower symmetry than it should have (P1 for a relaxed wurtzite, say), it is slightly
   distorted: rerun `info` with `--symprec 1e-3` (Å), check that the expected space group appears,
   pass the same `--symprec` to `generate`, and tell the user.
5. **Do not run DFT codes** unless the user asks and the code is available. Hand over the run step: every
   directory listed in `DIR/STATES` is one calculation, for example one task of a job array.
6. **Report** the constants in GPa with the space group, the number of strain patterns and calculations,
   the consistency check and Born stability that minelastic prints, and every assumption you made.

## Pick the route

| situation | route |
|---|---|
| an ASE calculator: EMT, MACE, CP2K or another code through ASE | Python, `minelastic.calculate(atoms, calc, ...)` |
| input files of a DFT code | shell, `minelastic info`, `generate`, run, `collect` |
| stresses from somewhere else | Python, `ElasticSetup(atoms)` and `setup.compute(stresses)` ([quantities](references/quantities.md#your-own-stresses)) |

With an ASE calculator, call `calculate`; do not rebuild it from `generate generic` and a loop over
directories.

## ASE calculators

<!-- run -->
```python
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
import minelastic

atoms = bulk("Cu", "hcp", a=2.55, c=4.16)               # or ase.io.read("POSCAR")
atoms.calc = EMT()
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-6)      # stress-free reference first

clamped = minelastic.calculate(atoms, EMT())                     # clamped ions
relaxed = minelastic.calculate(atoms, EMT(), relax_ions=True, fmax=1e-6)  # ions relaxed in every cell
ionic = minelastic.calculate(atoms, EMT(), ionic=True)           # relaxed ions from displaced atoms
constants = relaxed.as_dict()                                    # {"c11": ..., "c12": ...}, GPa
relaxed_ionic = dict(zip(ionic.setup.form.labels, ionic.setup.form.project(ionic.C_relaxed)))
print(relaxed.summary())
```

* `result.as_dict()` gives the independent constants of `result.C`. With `ionic=True`, `C` and `as_dict()`
  stay clamped-ion; the relaxed-ion matrix is `result.C_relaxed`, projected onto the constants as above.
* `relax_ions=True` relaxes the atoms in every strained cell; `ionic=True` adds single points with displaced
  atoms instead (the method of VASP's `IBRION = 6`) and also gives the force constants. They agree to second
  order in the strain step. If every atom sits at a centre of inversion (fcc, bcc, rocksalt), C equals the
  clamped-ion C.
* `dipole=f` (a function of the `Atoms` returning the dipole in e Å) adds the Born charges
  (`result.born`) and the piezoelectric tensors (`result.piezo`, `result.piezo_total`, `result.piezo_d`).
* Other results: `result.properties` (Voigt-Reuss-Hill moduli and more), `result.internal_strain`
  (eV/Å), `result.C_input_frame`. Options of the strain set: `delta` (step, default 0.01), `stencil`,
  `symprec`, `scheme`, `frame`, `displacement`.
* Run machine-learned potentials in double precision (MACE: `default_dtype="float64"`): the stresses of
  the strained cells differ by only about 1 GPa.

## DFT codes

<!-- shell -->
```bash
minelastic info --code vasp vasp_aln          # space group, constants, strain set; writes nothing
minelastic generate vasp vasp_aln -o el       # one directory per cell, listed in el/STATES
```

<!-- skip: needs the finished DFT runs -->
```bash
for d in $(cat el/STATES); do (cd el/$d && mpirun vasp_std > out); done   # or one job per line
minelastic collect el                         # prints the constants, writes el/elastic.json
```

| code | template |
|---|---|
| `vasp` | folder with `POSCAR`, `INCAR`, `POTCAR` and `KPOINTS` (or `KSPACING` in the INCAR) |
| `qe` | the pw.x input file; `pseudo_dir` must be an absolute path, because pw.x runs in sub-folders |
| `abinit` | folder with `POSCAR` and one `*.abi` that reads it with `structure "poscar:POSCAR"` |
| `aims` | folder with `geometry.in` and `control.in` (with `compute_analytical_stress .true.`) |
| `generic` | any structure file ASE can write, with `--format`, `--filename` and `--extra` files |

What `generate` changes in each input, and which outputs `collect` reads: [codes](references/codes.md).

`generate` takes `-o DIR`, `--stencil four-point|central|forward` (default four-point: 4 cells per
pattern; central: 2, the default with `--md`), `--delta` (default 0.01), `--symprec`, `--scheme asess`
(the six single strains, for comparison) and, for VASP, `--fix-incar` and `--potcar link`.

| quantity | `generate` | `collect` |
|---|---|---|
| clamped-ion constants | | |
| relaxed-ion constants by relaxations | `--relax-ions` | |
| relaxed-ion constants from displaced atoms, force constants | `--ionic` | |
| internal-strain tensor | | `--internal-strain` |
| clamped-ion piezoelectric tensor (VASP) | `LCALCPOL = .TRUE.` in the template INCAR | `--piezo` |
| relaxed-ion e and d, Born charges (VASP) | `LCALCPOL = .TRUE.` in the template INCAR, `--ionic` | `--piezo` |
| Born charges and clamped-ion e from DFPT (VASP) | `--ionic` | `--born OUTCAR` of a `LEPSILON` run |
| finite temperature from MD (VASP, pw.x) | NVT template at the thermal volume, `--md` | `--md --skip N --kinetic-temperature T` |
| finite temperature from snapshots | single points in `conf-*` sub-folders | `--samples 'conf-*'` |

* `generate` never adds `LCALCPOL`: put it in the template INCAR before `generate`. minelastic reads
  Berry-phase dipoles from VASP outputs only; for other codes see
  [quantities](references/quantities.md#piezoelectricity).
* With `--ionic`, `constants` in `elastic.json` are clamped-ion; the relaxed-ion matrix is `C_relaxed`.
* MD: `--skip N` discards the first N steps as equilibration, and `--kinetic-temperature T` adds the
  kinetic stress n k_B T, which VASP and pw.x leave out. `collect` refuses MD runs without `--md` or
  `--samples`. Give the user the complete `collect` command.

## Check the result

* **Consistency.** `collect` compares every stress component with the fitted tensor. A large deviation
  means a lower symmetry than assumed, a step outside the linear regime, loose convergence, k-point meshes
  that differ between cells, or MD noise. `--method lstsq` shows how much the constants depend on which
  components are used.
* **Stability.** Born stability is printed. A negative eigenvalue points to a reference that is not at
  equilibrium or to an unstable structure.
* **VASP's `IBRION = 6`.** Differences of a few tenths of a GPa are expected (another numerical
  derivative; `--stress-measure kirchhoff` differentiates J sigma as VASP does); several GPa mean other
  settings. `minelastic vasp-tensor OUTCAR` converts VASP's own tensor to GPa and Voigt order.
* **Orientation.** The textbook constants refer to the standard orientation of the crystal; if the input
  cell had another orientation, it was rotated rigidly, and `C_input_frame` is the tensor in the input
  orientation.

## Errors

* `--ionic` stops with an unstable or soft mode: the reference is not at a minimum of the energy; relax the
  atoms tightly.
* Noncollinear inputs are refused when the cell must be rotated: use `--frame input`.
* Atoms of one element in different VASP species blocks or with different pw.x labels count as different
  species, as in the DFT code.

More: [codes](references/codes.md), [quantities](references/quantities.md), and the troubleshooting page of
the documentation.
