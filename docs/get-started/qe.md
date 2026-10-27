# First run with Quantum ESPRESSO

Here we compute the five elastic constants of wurtzite AlN with pw.x, starting from the input files of
the relaxed structure. This takes four single-point calculations, 9 s in total on
one 32-core node. The output shown below is from these runs.

## 1. Prepare the template

minelastic writes the input of every strained cell from a template: the pw.x input of the relaxed reference
structure. We use
[`examples/qe_aln/pw.in`](https://github.com/utksi/minelastic/tree/main/examples/qe_aln) from the
repository:

<!-- file: examples/qe_aln/pw.in -->
```fortran
&CONTROL
  calculation = 'scf'
  prefix = 'aln'
  pseudo_dir = './pseudo'
  outdir = './tmp'
/
&SYSTEM
  ibrav = 4, celldm(1) = 5.911804, celldm(3) = 1.603472
  nat = 4, ntyp = 2
  ecutwfc = 50, ecutrho = 400
  occupations = 'fixed'
/
&ELECTRONS
  conv_thr = 1.0d-10
/
&IONS
/
ATOMIC_SPECIES
Al 26.982 Al.pbe-n-kjpaw_psl.1.0.0.UPF
N  14.007 N.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
Al 0.3333333333 0.6666666667 0.0000
Al 0.6666666667 0.3333333333 0.5000
N  0.3333333333 0.6666666667 0.3818
N  0.6666666667 0.3333333333 0.8818
K_POINTS automatic
8 8 5 0 0 0
```

It uses the PBE functional and the PAW potentials `Al.pbe-n-kjpaw_psl.1.0.0.UPF` and
`N.pbe-n-kjpaw_psl.1.0.0.UPF` of pslibrary 1.0.0, with tight cutoffs and convergence threshold; the
structure was relaxed with these settings.

!!! note "Pseudopotential path"

    Every calculation runs in a directory of its own, so set `pseudo_dir` to the absolute path of the
    directory with the pseudopotential files.

!!! info "Convergence"

    Elastic constants are obtained from derivatives of the stress with respect to the strain, which are estimated
    from stress differences of a few GPa. They are therefore much more sensitive to the cutoffs, the k-point
    mesh, the smearing and the convergence threshold than total energies are. The settings here are chosen to
    show the workflow; they are not the result of a convergence study. Before using the constants in a
    publication, check that they change by less than the accuracy you need when the cutoffs are raised and the
    k-point mesh refined, and relax the reference structure with the final settings.


## 2. Look at the strain set

`minelastic info` analyses the structure without writing anything:

<!-- shell -->
```console
$ minelastic info --code qe qe_aln/pw.in --stencil central
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
$ minelastic generate qe qe_aln/pw.in -o el --stencil central
...
4 calculations written to el/ (list in el/STATES)
$ cat el/STATES
s1_m1
s1_p1
s2_m1
s2_p1
$ ls el/s1_p1
pw.in
```

`generate` prints the same analysis as `info` (left out here) and writes one directory per strained cell:
`s1_p1` is pattern 1 at $+0.01$, `s2_m1` pattern 2 at $-0.01$. Each `pw.in` is the template with the
strained cell written out (`ibrav = 0`, `CELL_PARAMETERS` in Å, `ATOMIC_POSITIONS crystal`) and with
`tstress` and `tprnfor` switched on; everything else is copied unchanged.

## 4. Run pw.x

The four calculations are independent. Run them one after another,

<!-- skip: runs pw.x -->
```bash
for d in $(cat el/STATES); do (cd el/$d && mpirun pw.x -in pw.in > pw.out); done
```

or at the same time, for example as a Slurm job array:

<!-- skip: submits a Slurm job -->
```bash
#!/bin/bash
#SBATCH --array=1-4
d=el/$(sed -n "${SLURM_ARRAY_TASK_ID}p" el/STATES)
cd "$d" && srun pw.x -in pw.in > pw.out
```

## 5. Collect the results

<!-- shell -->
```console
$ minelastic collect el
minelastic 0.1.0: P6_3mc (186), hexagonal, 2 strain pattern(s), method 'uncoupled'
  c11 =     445.37 GPa
  c12 =     100.00 GPa
  c13 =      64.21 GPa
  c33 =     437.66 GPa
  c44 =     125.02 GPa
stiffness matrix (GPa, Voigt xx yy zz yz xz xy, working frame):
     445.37    100.00     64.21      0.00      0.00      0.00
     100.00    445.37     64.21      0.00      0.00      0.00
      64.21     64.21    437.66      0.00      0.00      0.00
       0.00      0.00      0.00    125.02      0.00      0.00
       0.00      0.00      0.00      0.00    125.02      0.00
       0.00      0.00      0.00      0.00      0.00    172.68
Born stable: True (eigenvalues 125.0, 125.0, 172.7, 345.4, 385.9, 597.1)
K_VRH = 198.18  G_VRH = 155.29  E = 369.38 GPa  nu = 0.189  A_U = 0.171
consistency (redundant stress components): largest |computed - fitted| = 0.00 GPa (pattern 1, xz)

written el/elastic.json
```

`collect` reads the stress of each cell from `pw.out` and prints:

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
{'c11': 445.37, 'c12': 100.0, 'c13': 64.21, 'c33': 437.66, 'c44': 125.02}
bulk modulus 198.2 GPa, shear modulus 155.3 GPa
```

`res["C"]` is the 6 x 6 matrix. The other entries are listed under
[Command line](../reference/cli.md#output-files).

## 7. Check the result

The same calculation with the six single strains (`--scheme asess`, 12 cells) gives, with the same
settings:

| GPa | c11 | c12 | c13 | c33 | c44 |
|---|---|---|---|---|---|
| minelastic (4 cells) | 445.37 | 100.00 | 64.21 | 437.66 | 125.02 |
| six single strains (12 cells) | 445.35 | 100.04 | 64.18 | 437.66 | 124.97 |

The two agree to 0.05 GPa with a third of the calculations. The VASP tutorial gives values a few GPa
lower, because the PAW datasets and cutoffs differ.

The forces in the same four cells also give the internal-strain tensor, the force on each atom per unit
strain, without further calculations:

<!-- shell -->
```console
$ minelastic collect el --internal-strain
...
internal strain tensor dF/de (eV/A), columns xx yy zz yz xz xy:
  ion   1 x     0.0000    -0.0000     0.0000     0.0000    -5.1051    -8.5575
  ion   1 y    -8.5576     8.5575     0.0000    -5.1051    -0.0000     0.0000
  ion   1 z    -5.2838    -5.2838    11.4889    -0.0000     0.0000     0.0000
...
```

The relaxed-ion constants need the ions relaxed in each strained cell, or a few more single points with
displaced atoms; both are described in
[Relaxed-ion constants](../how-to/relaxed-ions.md#with-quantum-espresso).
