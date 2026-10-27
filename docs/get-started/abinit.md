# First run with ABINIT

Here we compute the five elastic constants of wurtzite AlN with ABINIT, starting from the input files of
the relaxed structure. This takes four single-point calculations, 373 s in total
on one 32-core node. The output shown below is from these runs.

## 1. Prepare the template

minelastic writes the input of every strained cell from a template: a directory with the reference structure
as a `POSCAR` and one ABINIT input file that reads it. We use
[`examples/abinit_aln`](https://github.com/utksi/minelastic/tree/main/examples/abinit_aln) from the
repository. The input file `aln.abi` is

<!-- file: examples/abinit_aln/aln.abi -->
```text
# wurtzite AlN, PBE, PseudoDojo v0.4 standard norm-conserving pseudopotentials (Al.psp8, N.psp8)
# the structure is read from POSCAR; set pp_dirpath to the directory with the pseudopotentials
structure "poscar:POSCAR"
pp_dirpath "./pseudo"
pseudos "Al.psp8, N.psp8"
ecut 48
ngkpt 8 8 5
nshiftk 1
shiftk 0 0 0
occopt 1
nstep 100
tolvrs 1.0d-18
chksymbreak 0
prtwf 0
prtden 0
prteig 0
```

and the POSCAR holds the structure, with the element of each atom at the end of its line, as ABINIT's
POSCAR reader requires:

<!-- file: examples/abinit_aln/POSCAR -->
```text
AlN, wurtzite, PBE (relaxed with VASP; see vasp_aln)
1.0
     3.1284000000000001    0.0000000000000000    0.0000000000000000
    -1.5641999999999994    2.7092738731992383    0.0000000000000000
     0.0000000000000000    0.0000000000000000    5.0163000000000002
Al N
2 2
Direct
   0.3333333333333333  0.6666666666666666  0.0000000000000000  Al
   0.6666666666666667  0.3333333333333334  0.4999999999999999  Al
   0.3333333333333333  0.6666666666666666  0.3817999999999999  N
   0.6666666666666667  0.3333333333333334  0.8817999999999998  N
```

The pseudopotentials are the PBE standard set of PseudoDojo v0.4. Every calculation runs in a directory of
its own, so set `pp_dirpath` to the absolute path of the directory with `Al.psp8` and `N.psp8`. The elastic
constants are differences of stresses of a few GPa, so the cutoff and the convergence threshold are tighter
than for a total-energy calculation.

## 2. Look at the strain set

`minelastic info` analyses the structure without writing anything:

<!-- shell -->
```console
$ minelastic info --code abinit abinit_aln --stencil central
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

Each pattern is applied at the strains $-0.01$ and $+0.01$ (multipliers -1 and 1 of the step). By default
minelastic also uses $\pm 0.02$, which makes the numerical derivative more accurate but doubles the number of
calculations. `--stencil central` keeps the cheaper two-point derivative for this first run.

## 3. Write the inputs

<!-- shell -->
```console
$ minelastic generate abinit abinit_aln -o el --stencil central
...
4 calculations written to el/ (list in el/STATES)
$ cat el/STATES
s1_m1
s1_p1
s2_m1
s2_p1
$ ls el/s1_p1
POSCAR
aln.abi
```

`generate` prints the same analysis as `info` (left out here) and writes one directory per strained cell:
`s1_p1` is pattern 1 at $+0.01$, `s2_m1` pattern 2 at $-0.01$. Each holds the strained POSCAR and a copy of
the input file. The input must keep the cell and the ions fixed; `generate` stops with an error if it relaxes
either.

## 4. Run ABINIT

The four calculations are independent. Run them one after another,

<!-- skip: runs ABINIT -->
```bash
for d in $(cat el/STATES); do (cd el/$d && mpirun abinit aln.abi > log); done
```

or at the same time, for example as a Slurm job array:

<!-- skip: submits a Slurm job -->
```bash
#!/bin/bash
#SBATCH --array=1-4
d=el/$(sed -n "${SLURM_ARRAY_TASK_ID}p" el/STATES)
cd "$d" && srun abinit aln.abi > log
```

## 5. Collect the results

<!-- shell -->
```console
$ minelastic collect el
minelastic 0.1.0: P6_3mc (186), hexagonal, 2 strain pattern(s), method 'uncoupled'
  c11 =     441.69 GPa
  c12 =      94.39 GPa
  c13 =      59.20 GPa
  c33 =     434.52 GPa
  c44 =     125.84 GPa
stiffness matrix (GPa, Voigt xx yy zz yz xz xy, working frame):
     441.69     94.39     59.20      0.00      0.00      0.00
      94.39    441.69     59.20      0.00      0.00      0.00
      59.20     59.20    434.52      0.00      0.00      0.00
       0.00      0.00      0.00    125.84      0.00      0.00
       0.00      0.00      0.00      0.00    125.84      0.00
       0.00      0.00      0.00      0.00      0.00    173.65
Born stable: True (eigenvalues 125.8, 125.8, 173.6, 347.3, 387.4, 583.2)
K_VRH = 193.55  G_VRH = 156.16  E = 369.19 GPa  nu = 0.182  A_U = 0.170
consistency (redundant stress components): largest |computed - fitted| = 0.00 GPa (pattern 1, xz)

written el/elastic.json
```

`collect` reads the stress of each cell from the newest `*.abo` file in its directory and prints:

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
{'c11': 441.69, 'c12': 94.39, 'c13': 59.2, 'c33': 434.52, 'c44': 125.84}
bulk modulus 193.6 GPa, shear modulus 156.2 GPa
```

`res["C"]` is the 6 x 6 matrix. The other entries are listed under
[Command line](../reference/cli.md#output-files).

## 7. Check the result

The same calculation with the six single strains (`--scheme asess`, 12 cells) gives, with the same
settings:

| GPa | c11 | c12 | c13 | c33 | c44 |
|---|---|---|---|---|---|
| minelastic (4 cells) | 441.69 | 94.39 | 59.20 | 434.52 | 125.84 |
| six single strains (12 cells) | 441.65 | 94.41 | 59.15 | 434.52 | 125.79 |

The two agree to 0.05 GPa with a third of the calculations, and both agree with ABINIT's own
density-functional perturbation theory (DFPT) tensor to 0.6 GPa.

The forces in the same four cells also give the internal-strain tensor, the force on each atom per unit
strain, without further calculations:

<!-- shell -->
```console
$ minelastic collect el --internal-strain
...
internal strain tensor dF/de (eV/A), columns xx yy zz yz xz xy:
  ion   1 x    -0.0000     0.0000     0.0000     0.0000    -5.0587    -8.4931
  ion   1 y    -8.4931     8.4931     0.0000    -5.0587     0.0000     0.0000
  ion   1 z    -5.2537    -5.2536    11.4090     0.0000    -0.0000    -0.0000
...
```

The relaxed-ion constants need the ions relaxed in each strained cell, or a few more single points with
displaced atoms; both are described in [Relaxed-ion constants](../how-to/relaxed-ions.md).
