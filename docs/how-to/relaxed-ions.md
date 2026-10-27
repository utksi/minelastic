# Relaxed-ion constants

When a crystal is strained, atoms that are not fixed by symmetry move inside the cell. The constants with
the atoms held at their strained fractional positions are the clamped-ion constants $C^0$. Experiments
measure the relaxed-ion constants, because at acoustic frequencies the atoms follow the strain:

$$
C = C^0 - \frac{1}{\Omega}\, \Lambda^{\mathsf T} \Phi^{+} \Lambda ,
$$

lowered by the relaxation of the atoms. Here $\Lambda$ is the internal-strain tensor (force on each atom
per unit strain), $\Phi$ the $\Gamma$-point force constants and $\Omega$ the cell volume
([theory](../theory/relaxed-ions.md)). If every atom sits at a centre of inversion, as in fcc, bcc,
rocksalt or cubic perovskites, nothing relaxes and $C = C^0$.

minelastic computes $C$ in two ways:

| | relax the ions | displace the atoms |
|---|---|---|
| option | `--relax-ions`, `relax_ions=True` | `--ionic`, `ionic=True` |
| calculations | one relaxation per strained cell | one single point per strained cell, plus one per symmetry-irreducible atomic displacement and step |
| also gives | | the force constants and, with dipoles, the Born charges and piezoelectric tensors |
| watch out for | relaxations must be converged tightly: a strain of 0.01 induces forces of only about 0.05 eV/Å | an unstable or soft optical mode makes the relaxed-ion response undefined; minelastic stops with an error |

The second route is the linear response itself, the method of VASP's `IBRION = 6`. It needs no
relaxations, which makes it cheaper and more robust in most cases. The two routes differ by terms of second
order in the strain step.

## With an ASE calculator

[`examples/emt_ionic.py`](https://github.com/utksi/minelastic/tree/main/examples/emt_ionic.py) relaxes hcp
Cu with EMT and computes the relaxed-ion constants by both routes, and by relaxing the ions in the six single
strains:

<!-- excerpt: examples/emt_ionic.py -->
```python
relaxed = minelastic.calculate(atoms, EMT(), relax_ions=True, fmax=1e-6)          # route 1
ionic = minelastic.calculate(atoms, EMT(), ionic=True)                             # route 2
six = minelastic.calculate(atoms, EMT(), relax_ions=True, fmax=1e-6, scheme="asess")
```

<!-- output-of: examples/emt_ionic.py -->
```text
reference: a = 2.5386 A, c = 4.1430 A

        clamped C0  relaxed (1)  ionic (2)   six (3)   GPa
c11         233.60       216.34     216.34    216.35
c12          94.81       112.07     112.07    112.07
c13          74.77        74.77      74.77     74.78
c33         254.00       254.00     254.00    254.00
c44          49.30        49.30      49.30     49.30

(1) 8 relaxations; (2) 8 strained and 8 displaced-atom single points; (3) 24 relaxations
largest difference (2) - (1): 0.003 GPa

Gamma-point force constants of atom 1 on itself (eV/A^2):
[[ 1.6453  0.     -0.    ]
 [ 0.      1.6453 -0.    ]
 [-0.     -0.      8.4942]]
```

The routes agree to 0.003 GPa. hcp has no optical mode that couples linearly to the yz and xz shears, so c44
is the same with clamped and relaxed ions.

## With VASP

<!-- shell -->
```bash
minelastic generate vasp vasp_aln -o el_relaxed --relax-ions --stencil central     # route 1
minelastic generate vasp vasp_aln -o el_ionic --ionic --displacement 0.015 --stencil central   # route 2
```

With `--relax-ions`, `generate` sets `IBRION = 2`, `NSW = 99`, `ISIF = 2` and, unless the template sets it,
`EDIFFG = -5E-4`. With `--ionic`, wurtzite AlN needs its two inequivalent atoms displaced along x and z: with
two points per displacement, 8 displaced cells in addition to the 4 strained cells. VASP's own `IBRION = 6`
uses the same 8 displaced cells but 12 strained cells, because it always applies the six single strains.
For the template of the [VASP tutorial](../get-started/vasp.md), with two points per strain:

| GPa | c11 | c12 | c13 | c33 | c44 | c66 |
|---|---|---|---|---|---|---|
| clamped ions $C^0$ (4 cells) | 440.26 | 96.78 | 62.16 | 432.32 | 125.02 | 171.74 |
| relaxed ions, route 1 (4 relaxations) | 375.66 | 127.98 | 97.63 | 355.22 | 112.26 | 123.84 |
| relaxed ions, route 2 (4 + 8 single points) | 376.11 | 128.41 | 97.71 | 354.60 | 112.47 | 123.85 |
| VASP `IBRION = 6`, `ISIF = 3` (12 + 8 single points) | 376.33 | 128.55 | 97.79 | 355.24 | 112.57 | 123.89 |

The differences of a few tenths of a GPa come from the numerical derivatives: VASP strains the cell by 0.015
and differentiates $J\sigma$, minelastic used 0.01 and the stress itself. The run of route 2 also gives the
Born charges and the piezoelectric tensors ([Born charges and piezoelectricity](piezoelectricity.md)).

## With Quantum ESPRESSO

<!-- shell -->
```bash
minelastic generate qe qe_aln/pw.in -o el_relaxed --relax-ions --stencil central
```

sets `calculation = 'relax'` with `forc_conv_thr = 2.0d-5` and `etot_conv_thr = 1.0d-8` unless the template
sets them; the pw.x defaults stop the relaxation before the internal strain is relaxed. For the template of
the [Quantum ESPRESSO tutorial](../get-started/qe.md):

| GPa | c11 | c12 | c13 | c33 | c44 |
|---|---|---|---|---|---|
| minimal set (4 relaxations) | 377.67 | 132.00 | 101.72 | 356.99 | 111.51 |
| six single strains (12 relaxations) | 377.62 | 132.09 | 101.70 | 356.99 | 111.45 |

Route 2 with pw.x (`--ionic`, 12 single points) agrees with the relaxations to 1.1 GPa.

## With other codes

Both routes work with every interface. With ABINIT, `--relax-ions` appends `ionmov 2`, `ntime 100` and
`tolmxf 1.0d-6` unless the input sets `ionmov`; with FHI-aims, it adds `relax_geometry bfgs 1E-3`. With CP2K
through ASE, `ionic=True` lowers c44 of diamond Si from the clamped-ion 102.6 GPa to 79.7 GPa
([CP2K tutorial](../get-started/cp2k.md)).
