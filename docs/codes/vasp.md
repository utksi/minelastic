# VASP

A complete example is the [VASP tutorial](../get-started/vasp.md). Below we describe how minelastic
handles the VASP inputs and outputs.

## Template and commands

The template directory holds the input of the relaxed reference structure: `POSCAR`, `INCAR`, `POTCAR`, and
`KPOINTS` or `KSPACING` in the INCAR.

<!-- skip: needs a template directory ref/ and the VASP runs -->
```bash
minelastic info --code vasp ref/             # symmetry, tensor form and strain set; no files written
minelastic generate vasp ref/ -o el          # one directory per state, listed in el/STATES
minelastic collect el                        # after the runs; writes el/elastic.json
```

Each `el/<state>/` holds `POSCAR` (the strained or displaced cell, with the species lines and the
selective-dynamics flags of the template), `INCAR`, `KPOINTS` and `POTCAR` (`--potcar link` symlinks it).
`el/minelastic.json` records the setup; `collect` needs it.

## The INCAR

`generate` copies the template's INCAR and changes only what the type of calculation requires:

* `IBRION = -1`, `NSW = 0`, `ISIF = 2` for clamped ions (`ISIF = 2` computes the stress);
* `IBRION = 2`, `NSW = 99`, `ISIF = 2` and, unless set, `EDIFFG = -5E-4` with `--relax-ions`;
* with `--md`, the template's `IBRION = 0`, `NSW` and `POTIM` are kept and `ISIF = 2` is set;
* `NFREE`, `POTIM`, `LEPSILON` and `LCALCEPS` are removed, so that an `IBRION = 6` input can be used as it
  is (relaxations then use VASP's default `POTIM`);
* `KSPACING` is replaced by the explicit k-point mesh of the reference cell, so that all cells use the same
  mesh.

Every change is listed. The settings that matter for stress differences of a few GPa are checked, but left as
they are: `EDIFF` at most 1E-6 (1E-8 suggested; 1E-6 for MD), `ENCUT` at least 1.3 times the largest
`ENMAX` of the POTCAR, `PREC = Accurate` (or `High`) and `LREAL = .FALSE.`. The template of the tutorial
meets all of them. The INCAR of an `IBRION = 6` run with VASP's defaults gives:

<!-- shell -->
```console
$ cat ibrion6/INCAR
IBRION = 6
ISIF = 3
NFREE = 2
POTIM = 0.015
$ minelastic generate vasp ibrion6 -o el
...
NOTE: INCAR: set IBRION = -1 (template: 6), NSW = 0, ISIF = 2 (template: 3); removed NFREE, POTIM
NOTE: INCAR: recommended for small stress differences: EDIFF default 1E-4 -> 1E-8, ENCUT default ENMAX = 400 -> 520, PREC default Normal -> Accurate. The template values are kept; --fix-incar (fix_incar=True) applies all of these
```

`--fix-incar` applies all of the recommended values. There is no partial fix: the INCAR is either the
template's or follows all the recommendations.

## Running

Any scheduler works. For example, a Slurm job array over the lines of `STATES`:

<!-- skip: submits a Slurm job -->
```bash
#!/bin/bash
#SBATCH --array=1-8              # one task per line of el/STATES
d=el/$(sed -n "${SLURM_ARRAY_TASK_ID}p" el/STATES)
cd "$d" && srun vasp_std > out
```

## Reading the results

`collect` reads the stress and the forces of each directory from `vasprun.xml`, or from `OUTCAR` if there is
no readable `vasprun.xml`. VASP prints the stress in kB, positive in compression, in the order XX YY ZZ XY YZ
ZX; minelastic converts it to GPa, positive in tension, in Voigt order. With `--md` the stress of every ionic
step is read; with `--piezo` the Berry-phase dipole of each cell (`LCALCPOL`); with `--born OUTCAR` the Born
charges and the clamped-ion piezoelectric tensor of a DFPT (`LEPSILON`) or `LCALCEPS` run.

`minelastic vasp-tensor OUTCAR` reads the tensor that VASP prints for `IBRION = 6` with `ISIF = 3` or more,
converts it to GPa and Voigt order, and adds the Voigt-Reuss-Hill moduli and the Born stability test
([Elastic constants in practice](../how-to/elastic.md#compare-with-vasps-ibrion-6)).

## Notes

* Use the same `ENCUT`, `PREC`, k-point mesh and `EDIFF` for the reference and all strained cells. The
  plane-wave basis changes with the volume of the strained cell, and a raised cutoff (1.3 to 1.4 times
  ENMAX) keeps the resulting Pulay stress small.
* `MAGMOM`, `LDAU*`, `ISPIN` and all other tags are copied unchanged, and `MAGMOM` enters the symmetry
  analysis ([Magnetic materials](../how-to/magnetic.md)). If the cell has to be rotated to its standard
  orientation, direction-dependent inputs (`LNONCOLLINEAR`, `LSORBIT`, `SAXIS`, `EFIELD`, Cartesian k points)
  are refused; use `--frame input`.
* The elastic constants, the internal-strain tensor, the relaxed-ion constants and the piezoelectric tensors
  of the tutorial structure are on the how-to pages:
  [elastic constants](../how-to/elastic.md), [relaxed ions](../how-to/relaxed-ions.md#with-vasp),
  [internal strain](../how-to/internal-strain.md#with-vasp) and
  [piezoelectricity](../how-to/piezoelectricity.md#with-vasp).
