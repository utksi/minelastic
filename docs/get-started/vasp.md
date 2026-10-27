# First run with VASP

Here we compute the five elastic constants of wurtzite AlN with VASP, starting from the input files of
the relaxed structure. This takes four single-point calculations, 42 s in total on
one 32-core node. The output shown below is from these runs.

## 1. Prepare the template

minelastic writes the input of every strained cell from a template: a directory with the VASP input of the
relaxed reference structure. We use
[`examples/vasp_aln`](https://github.com/utksi/minelastic/tree/main/examples/vasp_aln) from the repository.
Copy the directory and add a PBE POTCAR for Al and N, in the order of the POSCAR. The INCAR is

<!-- file: examples/vasp_aln/INCAR -->
```text
PREC = Accurate
ENCUT = 520
EDIFF = 1E-8
ISMEAR = 0
SIGMA = 0.05
LREAL = .FALSE.
LWAVE = .FALSE.
LCHARG = .FALSE.
```

and the KPOINTS file asks for a Gamma-centred 8 x 8 x 5 mesh. The cutoff is 1.3 times the largest ENMAX of
the POTCAR (520 eV for N), the electronic convergence is tight, real-space projection is off, and the
structure in the POSCAR was relaxed with these settings.

!!! info "Convergence"

    Elastic constants are differences of stresses of a few GPa, so they are much more sensitive to the
    plane-wave cutoff, the k-point mesh, the smearing and the electronic convergence than total energies are. The settings here are chosen to show the workflow; they are not the
    result of a convergence study. Before using the constants in a publication, check that they change by
    less than the accuracy you need when the cutoff is raised and the k-point mesh refined, and relax the reference structure with the final
    settings.


## 2. Look at the strain set

`minelastic info` analyses the structure without writing anything:

<!-- shell -->
```console
$ minelastic info --code vasp vasp_aln --stencil central
space group P6_3mc (186), point group 6mm, elastic type hexagonal
5 independent constants: c11, c12, c13, c33, c44
invariant strains (free lattice parameters) = minimum number of strain patterns: 2
stiffness matrix in the working frame:
          c11            c12            c13              0              0              0
          c12            c11            c13              0              0              0
          c13            c13            c33              0              0              0
            0              0              0            c44              0              0
            0              0              0              0            c44              0
            0              0              0              0              0  (c11 - c12)/2
strain set 'minimal': 2 pattern(s), step 0.01, multipliers (-1, 1)
  pattern 1: xx+yz          s1 = 1 c11, s2 = 1 c12, s3 = 1 c13, s4 = 1 c44
  pattern 2: zz             s1 = 1 c13, s2 = 1 c13, s3 = 1 c33
```

The crystal is hexagonal, with five independent elastic constants, and two strain patterns determine all of
them. Pattern 1 stretches the cell along x and shears it in the yz plane at the same time; pattern 2
stretches it along z. The entry `s1 = 1 c11` means that the first stress component of pattern 1,
$\sigma_{xx}$, equals $c_{11}$ times the strain. Every constant appears alone in at least one stress
component, so each is read off directly.

!!! note "Two points per strain"

    By default minelastic applies each pattern at $\pm 0.01$ and $\pm 0.02$, which makes the numerical
    derivative more accurate but doubles the number of calculations. For this first run, `--stencil central`
    keeps only $\pm 0.01$ (the multipliers -1 and 1 in the output above).


## 3. Write the inputs

<!-- shell -->
```console
$ minelastic generate vasp vasp_aln -o el --stencil central
...
4 calculations written to el/ (list in el/STATES)
NOTE: INCAR: set IBRION = -1, NSW = 0, ISIF = 2
$ cat el/STATES
s1_m1
s1_p1
s2_m1
s2_p1
$ ls el/s1_p1
INCAR
KPOINTS
POSCAR
POTCAR
```

`generate` prints the same analysis as `info` (left out here) and writes one directory per strained cell:
`s1_p1` is pattern 1 at $+0.01$, `s2_m1` pattern 2 at $-0.01$. Each holds the strained POSCAR and copies of
the INCAR, KPOINTS and POTCAR. The note lists every change to the INCAR: a single point with fixed ions
(`IBRION = -1`, `NSW = 0`) that computes the stress (`ISIF = 2`). If the template misses one of the settings
recommended above, a second note names it; `--fix-incar` applies all of them.

## 4. Run VASP

The four calculations are independent. Run them one after another,

<!-- skip: runs VASP -->
```bash
for d in $(cat el/STATES); do (cd el/$d && mpirun vasp_std > out); done
```

or at the same time, for example as a Slurm job array:

<!-- skip: submits a Slurm job -->
```bash
#!/bin/bash
#SBATCH --array=1-4
d=el/$(sed -n "${SLURM_ARRAY_TASK_ID}p" el/STATES)
cd "$d" && srun vasp_std > out
```

## 5. Collect the results

<!-- shell -->
```console
$ minelastic collect el
minelastic 0.1.0: P6_3mc (186), hexagonal, 2 strain pattern(s), method 'uncoupled'
  c11 =     440.26 GPa
  c12 =      96.78 GPa
  c13 =      62.16 GPa
  c33 =     432.32 GPa
  c44 =     125.02 GPa
stiffness matrix (GPa, Voigt xx yy zz yz xz xy, working frame):
     440.26     96.78     62.16      0.00      0.00      0.00
      96.78    440.26     62.16      0.00      0.00      0.00
      62.16     62.16    432.32      0.00      0.00      0.00
       0.00      0.00      0.00    125.02      0.00      0.00
       0.00      0.00      0.00      0.00    125.02      0.00
       0.00      0.00      0.00      0.00      0.00    171.74
Born stable: True (eigenvalues 125.0, 125.0, 171.7, 343.5, 382.4, 587.0)
K_VRH = 194.83  G_VRH = 154.66  E = 366.89 GPa  nu = 0.186  A_U = 0.165
consistency (redundant stress components): largest |computed - fitted| = 0.00 GPa (pattern 1, xz)

written el/elastic.json
```

`collect` reads the stress of each cell from `vasprun.xml` and prints:

* the five independent constants and the full stiffness matrix, in GPa, with $C_{66} = (c_{11} - c_{12})/2$;
* the Born stability test: the crystal is mechanically stable if all eigenvalues of the matrix are positive;
* the Voigt-Reuss-Hill bulk modulus K, shear modulus G and Young's modulus E, Poisson's ratio nu and the
  universal anisotropy index A_U;
* a consistency check: the stress components that are not needed for the constants are compared with the
  values that the fitted tensor predicts. A deviation of more than a few tenths of a GPa points to an
  unrelaxed reference structure, loose convergence or a lower symmetry than assumed.

## 6. Use the results

Everything is also written to `el/elastic.json`:

<!-- run -->
```python
import json

with open("el/elastic.json") as f:
    res = json.load(f)
print({name: round(value, 2) for name, value in res["constants"].items()})
print(f"bulk modulus {res['properties']['K_H']:.1f} GPa, shear modulus {res['properties']['G_H']:.1f} GPa")
```

<!-- output -->
```text
{'c11': 440.26, 'c12': 96.78, 'c13': 62.16, 'c33': 432.32, 'c44': 125.02}
bulk modulus 194.8 GPa, shear modulus 154.7 GPa
```

`res["C"]` is the 6 x 6 matrix. The other entries are listed under
[Command line](../reference/cli.md#output-files).

## 7. Check the result

The same calculation with the six single strains (`--scheme asess`, 12 cells) and VASP's own finite
differences (`IBRION = 6`, `ISIF = 3`, 12 cells) give, with the same settings:

| GPa | c11 | c12 | c13 | c33 | c44 | c66 |
|---|---|---|---|---|---|---|
| minelastic (4 cells) | 440.26 | 96.78 | 62.16 | 432.32 | 125.02 | 171.74 |
| six single strains (12 cells) | 440.21 | 96.81 | 62.12 | 432.32 | 124.96 | 171.70 |
| VASP `IBRION = 6` (12 cells) | 440.47 | 96.89 | 62.27 | 432.91 | 125.14 | 171.79 |

The differences of a few tenths of a GPa come from the numerical derivatives: VASP uses a strain of
0.015. With minelastic's default of four strains per pattern, the minimal set and the six single strains
agree to 0.005 GPa.

The forces in the same four cells also give the internal-strain tensor, the force on each atom per unit
strain, without further calculations:

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

The relaxed-ion constants, the Born charges and the piezoelectric tensors need a few more calculations, with
displaced atoms, described in [Relaxed-ion constants](../how-to/relaxed-ions.md#with-vasp) and
[Born charges and piezoelectricity](../how-to/piezoelectricity.md#with-vasp).
