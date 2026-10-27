# Born charges and piezoelectricity

The piezoelectric stress tensor is the change of the polarization with strain, and the strain coefficients
follow with the compliance $S = C^{-1}$:

$$
e_{iJ} = \frac{\partial P_i}{\partial \varepsilon_J} \quad (\mathrm{C/m^2}), \qquad
d = e\, S \quad (\mathrm{pC/N}).
$$

The Born effective charges $Z^*_{n,\alpha\beta} = \partial p_\alpha / \partial u_{n\beta}$ are the change of
the dipole $p$ of the cell when atom $n$ moves ([theory](../theory/piezoelectricity.md)).

minelastic takes the polarization of every cell from its dipole, for example from a Berry-phase calculation.
The strained cells alone give the clamped-ion tensor $e^0$, with no further calculations if the dipoles are computed
with the stresses. With `--ionic`, the displaced atoms add the Born charges, and the relaxed-ion (total)
tensor follows as $e = e^0 + Z^* \Phi^+ \Lambda / \Omega$, together with the relaxed-ion elastic constants
and $d$.

## With an ASE calculator

`calculate` takes the dipole as a function of the `Atoms` of each cell, returning the dipole in e Å.
[`examples/bks_quartz_piezo.py`](https://github.com/utksi/minelastic/tree/main/examples/bks_quartz_piezo.py)
uses the BKS rigid-ion model of silica (charges $q_{Si} = 2.4$, $q_O = -1.2$ and Buckingham terms, as a small
ASE calculator in the script) for alpha-quartz, point group 32. In a rigid-ion model the dipole is
$p = \sum_n q_n r_n$:

<!-- excerpt: examples/bks_quartz_piezo.py -->
```python
def dipole(a):
    return q @ a.positions

ionic = minelastic.calculate(atoms, BKS(), ionic=True, dipole=dipole, delta=0.005, displacement=0.005)
relaxed = minelastic.calculate(atoms, BKS(), relax_ions=True, fmax=1e-7, steps=4000, dipole=dipole, delta=0.005)
```

The results are `ionic.born` ($Z^*$), `ionic.piezo` ($e^0$), `ionic.piezo_total` ($e$) and `ionic.piezo_d`
($d$). In a rigid-ion model $Z^*$ equals the charges and $e^0 = 0$: the whole response comes from the
relaxation of the ions. Route (1) obtains it from clamped single points with displaced atoms, route (2) from
the stresses and dipoles of cells with relaxed ions:

<!-- output-of: examples/bks_quartz_piezo.py -->
```text
BKS alpha-quartz: a = 4.9593 A, c = 5.4572 A

relaxed-ion constants (GPa)   (1) ionic=True   (2) relax_ions=True
  c11                              93.64       93.63
  c33                             103.22      103.11
  c44                              49.12       49.01
  c66                              43.07       43.07
  c12                               7.50        7.49
  c13                              14.24       14.16
  c14                              18.01       18.02
piezoelectric tensor (C/m^2)
  e11                             0.1409      0.1411
  e14                            -0.0027     -0.0028

Born charges of Si (diagonal, |e|): [2.4 2.4 2.4]; clamped-ion e0 = 0.0000 C/m^2
d11 = 1.946 pC/N, d14 = -1.483 pC/N
largest differences of the routes: 0.117 GPa, 1.8e-04 C/m^2
```

The two routes agree to 0.1 GPa and 2e-4 C/m$^2$ with a strain step of 0.005; with 0.01, the relaxations at
finite strain add 0.6 GPa to route (2). BKS gives c11 = 93.6 GPa and e11 = 0.141 C/m$^2$, against 86.6 GPa
and 0.171 C/m$^2$ measured for quartz.

## With VASP

### Inputs

Add `LCALCPOL = .TRUE.` to the INCAR of the [tutorial template](../get-started/vasp.md), so that VASP
computes the Berry-phase dipole of every cell, and add the atomic displacements:

<!-- shell -->
```console
$ minelastic generate vasp vasp_aln -o el --ionic --displacement 0.015 --stencil central
...
strain set 'minimal': 2 pattern(s), step 0.01, multipliers (-1, 1)
atomic displacements (symmetry-irreducible), step 0.015 A: atom 1 x, atom 1 z, atom 3 x, atom 3 z
  pattern 1: xx+yz          s1 = 1 c11, s2 = 1 c12, s3 = 1 c13, s4 = 1 c44
  pattern 2: zz             s1 = 1 c13, s2 = 1 c13, s3 = 1 c33

12 calculations written to el/ (list in el/STATES)
NOTE: INCAR: set IBRION = -1, NSW = 0, ISIF = 2
```

Wurtzite has two inequivalent atoms, Al and N, each displaced along x and z (y follows from the 3-fold axis)
by $\pm 0.015$ Å. Together with the 4 strained cells that makes 12 single points. VASP's own `IBRION = 6`
uses the same 8 displaced cells, but 12 strained cells. Like the tutorial, this run uses two points per
strain and displacement, as VASP does with its default `NFREE = 2`; the default four points double every
count.

### Collect

<!-- shell -->
```bash
for d in $(cat el/STATES); do (cd el/$d && mpirun vasp_std > out); done
```

<!-- shell -->
```console
$ minelastic collect el --piezo
...
relaxed-ion stiffness matrix (GPa):
     376.11    128.41     97.71      0.00      0.00     -0.00
     128.41    376.11     97.71      0.00     -0.00     -0.00
      97.71     97.71    354.60     -0.00     -0.00      0.00
       0.00      0.00     -0.00    112.47      0.00     -0.00
       0.00     -0.00     -0.00      0.00    112.47     -0.00
      -0.00     -0.00      0.00     -0.00     -0.00    123.85
relaxed ions: Born stable True (eigenvalues 112.5, 112.5, 123.9, 247.7, 272.3, 586.8)
Born effective charges (|e|, from dipoles of the displaced atoms); acoustic sum 0.0007:
  ion   1    2.5230    0.0000    0.0000
  ion   1    0.0000    2.5230    0.0000
  ion   1    0.0000    0.0000    2.6807
...
total (relaxed-ion) piezoelectric tensor e (C/m^2, columns xx yy zz yz xz xy, working frame):
     0.0000    0.0000   -0.0000    0.0000   -0.3122   -0.0000
     0.0000    0.0000   -0.0000   -0.3122    0.0000   -0.0000
    -0.5924   -0.5924    1.4500    0.0000    0.0000   -0.0000
piezoelectric strain coefficients d = e S (relaxed ion) (pC/N, columns xx yy zz yz xz xy, working frame):
     0.0000   -0.0000   -0.0000    0.0000   -2.7762   -0.0000
     0.0000    0.0000   -0.0000   -2.7762    0.0000   -0.0000
    -2.2012   -2.2012    5.3023    0.0000    0.0000   -0.0000
...
```

The full output also lists the clamped-ion constants, the clamped-ion and ionic piezoelectric tensors, and
the ionic contribution to C. The sum of the Born charges over all atoms, zero for an exact calculation,
shows how well they are converged.

Instead of the Berry-phase dipoles, the Born charges and $e^0$ can come from a density-functional
perturbation theory run (`LEPSILON`) or an `LCALCEPS` run of the reference structure:

<!-- skip: needs the OUTCAR of a VASP LEPSILON run -->
```bash
minelastic collect el --born dfpt/OUTCAR
```

### Comparison with VASP's own tensors

The reference is VASP 6.6.1 with `IBRION = 6`, `ISIF = 3` and `LEPSILON = .TRUE.`, with the same structure,
cutoff and k-point mesh:

| | minelastic, Berry phase | minelastic, `--born` | VASP |
|---|---|---|---|
| strained + displaced cells | 4 + 8 | 4 + 8, and 1 DFPT run | 12 + 8, and DFPT |
| c11, c33 relaxed ions (GPa) | 376.11, 354.60 | 376.11, 354.60 | 376.33, 355.24 |
| c44, c66 relaxed ions (GPa) | 112.47, 123.85 | 112.47, 123.85 | 112.57, 123.89 |
| Z\*_xx, Z\*_zz of Al (\|e\|) | 2.523, 2.681 | 2.520, 2.680 | 2.520, 2.680 |
| e33 clamped ions (C/m$^2$) | -0.428 | -0.424 | -0.424 |
| e33, e31, e15 total (C/m$^2$) | 1.450, -0.592, -0.312 | 1.454, -0.589, -0.306 | 1.454, -0.588, -0.307 |
| d33, d31, d15 (pC/N) | 5.30, -2.20, -2.78 | 5.31, -2.20, -2.72 | 5.30, -2.19, -2.72 |

The elastic constants differ by the effect of the strain step (0.01 against VASP's 0.015, see
[Relaxed-ion constants](relaxed-ions.md#with-vasp)). The Berry-phase and DFPT routes to $Z^*$ and $e^0$
agree to 0.004 \|e\| and 0.005 C/m$^2$.

### Reproduce it without VASP

The stresses, forces and dipoles of these 12 cells are stored in `tests/data_ionic/aln_run_ionic.json`.
[`examples/replay_vasp_aln_ionic.py`](https://github.com/utksi/minelastic/tree/main/examples/replay_vasp_aln_ionic.py)
recomputes everything from them in under a second and compares with VASP:

<!-- output-of: examples/replay_vasp_aln_ionic.py -->
```text
12 single points; Voigt order xx yy zz yz xz xy

                              minelastic   VASP IBRION=6
  c11 clamped ions (GPa)       440.257     440.473
  c33 clamped ions (GPa)       432.317     432.910
  c11 relaxed ions (GPa)       376.110     376.325
  c33 relaxed ions (GPa)       354.597     355.240
  c44 relaxed ions (GPa)       112.467     112.569
  c66 relaxed ions (GPa)       123.852     123.886
  Lambda Al, z zz (eV/A)        11.095      11.130
  Lambda Al, x xy (eV/A)        -8.252      -8.243
  Lambda Al, x xz (eV/A)        -4.908      -4.897
  Z*_xx Al (|e|)                 2.523       2.520
  Z*_zz Al (|e|)                 2.681       2.680
  e33 clamped ions (C/m^2)      -0.428      -0.424
  e33 total (C/m^2)              1.450       1.454
  e31 total (C/m^2)             -0.592      -0.588
  e15 total (C/m^2)             -0.312      -0.307
  d33 (pC/N)                     5.302       5.298
  d31 (pC/N)                    -2.201      -2.191
  d15 (pC/N)                    -2.776      -2.724

internal-strain tensor, all atoms, largest difference to VASP: 0.040 eV/A from the forces of the strained cells, 0.057 eV/A from the stresses of the displaced atoms
acoustic sum of the Born charges: 0.0007 |e|
```

Lambda is the internal-strain tensor of an atom in eV/Å per unit strain, obtained once from the forces in
the strained cells and once from the stresses of the displaced atoms.

### Clamped-ion tensor from the strained cells alone

Without `--ionic`, `collect --piezo` on clamped cells computed with `LCALCPOL` gives the clamped-ion tensor
$e^0$ from the same 4 cells as the elastic constants: the polarization of the few strained cells determines
the whole tensor by symmetry. For this AlN example, e33 = -0.428 C/m$^2$, against -0.424 C/m$^2$ from
DFPT; [`examples/replay_vasp_aln.py`](own-workflow.md#bring-your-own-stresses) reproduces it from the stored
stresses and dipoles.

!!! note "Berry phase and relaxations"

    VASP evaluates the Berry phase only for the first ionic step. In a relaxation the dipole would belong to
    the unrelaxed structure, so minelastic refuses dipoles from relaxed-ion cells; the total tensor comes from
    `--ionic`, which needs single points only.

## With other codes

The elastic part of `--ionic` (force constants, internal strain, relaxed-ion constants) works with every
interface. minelastic reads Berry-phase dipoles from VASP output only. With other codes, pass the dipoles
(in e Å, physical sign) to `ElasticSetup.compute(stresses, forces=..., polarization={state_id: dipole})`,
or $Z^*$ and $e^0$ from a perturbation-theory run with `born=` and `clamped_piezo=`
([Your own workflow](own-workflow.md)).
