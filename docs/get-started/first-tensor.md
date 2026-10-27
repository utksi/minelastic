# Your first elastic tensor

Here we compute the elastic constants of copper with EMT, a simple interatomic potential that comes with
ASE. The calculation takes a few seconds on a laptop and needs nothing besides minelastic, yet the steps are
those of a DFT calculation: a relaxed reference structure, the strain set its symmetry requires, the stresses
of the strained cells, and the constants.

## 1. Relax the reference structure

The elastic constants are derivatives of the stress at the equilibrium structure, so the cell is relaxed
first, with the same calculator:

<!-- run -->
```python
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS

atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
atoms.calc = EMT()
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-6)
print(f"a = {atoms.cell[0, 0]:.4f} A")
```

<!-- output -->
```text
a = 3.5898 A
```

## 2. Look at the strain set

`ElasticSetup` finds the symmetry of the structure and the strains that it needs:

<!-- run -->
```python
import minelastic

setup = minelastic.ElasticSetup(atoms)
print(setup.describe())
```

<!-- output -->
```text
space group Fm-3m (225), point group m-3m, elastic type cubic
3 independent constants: c11, c12, c44
invariant strains (free lattice parameters) = minimum number of strain patterns: 1
stiffness matrix in the working frame:
c11  c12  c12    0    0    0
c12  c11  c12    0    0    0
c12  c12  c11    0    0    0
  0    0    0  c44    0    0
  0    0    0    0  c44    0
  0    0    0    0    0  c44
strain set 'minimal': 1 pattern(s), step 0.01, multipliers (-2, -1, 1, 2)
  pattern 1: xx+yz          s1 = 1 c11, s2 = 1 c12, s3 = 1 c12, s4 = 1 c44
```

A cubic crystal has three independent elastic constants, and one strain pattern determines all of them. The
pattern `xx+yz` stretches the cell along x and shears it in the yz plane at the same time. The relations on
the last line say which stress component gives which constant: the stress $\sigma_{xx}$ (s1) equals $c_{11}$
times the strain, $\sigma_{yy}$ and $\sigma_{zz}$ (s2, s3) give $c_{12}$, and the shear stress
$\sigma_{yz}$ (s4) gives $c_{44}$. The usual approach applies six separate strains instead.

The pattern is applied at four strains, $-0.02$, $-0.01$, $+0.01$ and $+0.02$, and the stresses are
differentiated numerically. These are the four strained cells:

<!-- run -->
```python
for state in setup.states:
    print(state["id"], state["strain"])
```

<!-- output -->
```text
s1_m2 [-0.02, -0.0, -0.0, -0.02, -0.0, -0.0]
s1_m1 [-0.01, -0.0, -0.0, -0.01, -0.0, -0.0]
s1_p1 [0.01, 0.0, 0.0, 0.01, 0.0, 0.0]
s1_p2 [0.02, 0.0, 0.0, 0.02, 0.0, 0.0]
```

The strains are written as six components in the order xx, yy, zz, yz, xz, xy, with the shear components
as engineering strains (twice the tensor component).

## 3. Compute the constants

`calculate` builds the strained cells, computes their stresses with the calculator and returns the result:

<!-- run -->
```python
result = minelastic.calculate(atoms, EMT())
print(result.summary())
```

<!-- output -->
```text
minelastic 0.1.0: Fm-3m (225), cubic, 1 strain pattern(s), method 'uncoupled'
  c11 =     172.61 GPa
  c12 =     115.43 GPa
  c44 =      89.90 GPa
stiffness matrix (GPa, Voigt xx yy zz yz xz xy, working frame):
     172.61    115.43    115.43      0.00      0.00      0.00
     115.43    172.61    115.43      0.00      0.00      0.00
     115.43    115.43    172.61      0.00      0.00      0.00
       0.00      0.00      0.00     89.90      0.00      0.00
       0.00      0.00      0.00      0.00     89.90      0.00
       0.00      0.00      0.00      0.00      0.00     89.90
Born stable: True (eigenvalues 57.2, 57.2, 89.9, 89.9, 89.9, 403.5)
K_VRH = 134.49  G_VRH = 56.88  E = 149.57 GPa  nu = 0.315  A_U = 1.755
consistency (redundant stress components): largest |computed - fitted| = 0.00 GPa ...
```

The summary lists

* the independent constants and the full stiffness matrix, in GPa;
* the Born stability test: the crystal is mechanically stable if all eigenvalues of the matrix are positive;
* the Voigt-Reuss-Hill bulk modulus K, shear modulus G and Young's modulus E, Poisson's ratio nu and the
  universal anisotropy index A_U;
* a consistency check: the stress components that are not needed for the constants (here $\sigma_{xz}$ and
  $\sigma_{xy}$, which must vanish) are compared with the values that the fitted tensor predicts.

## 4. Use the result

The result holds the numbers behind the summary:

<!-- run -->
```python
print({name: round(value, 2) for name, value in result.as_dict().items()})
print(result.C.round(1))
p = result.properties
print(f"K = {p['K_H']:.1f} GPa, G = {p['G_H']:.1f} GPa, E = {p['E_H']:.1f} GPa, nu = {p['nu_H']:.3f}")
```

<!-- output -->
```text
{'c11': 172.61, 'c12': 115.43, 'c44': 89.9}
[[172.6 115.4 115.4   0.    0.    0. ]
 [115.4 172.6 115.4   0.    0.    0. ]
 [115.4 115.4 172.6   0.    0.    0. ]
 [  0.    0.    0.   89.9   0.    0. ]
 [  0.    0.    0.    0.   89.9   0. ]
 [  0.    0.    0.    0.    0.   89.9]]
K = 134.5 GPa, G = 56.9 GPa, E = 149.6 GPa, nu = 0.315
```

`result.to_dict()` returns everything as a dictionary that can be saved as JSON:

<!-- run -->
```python
import json

with open("cu_elastic.json", "w") as f:
    json.dump(result.to_dict(), f, indent=1)
```

## 5. A crystal with internal relaxation

In hcp copper the atoms do not sit at centres of inversion, so a strain also moves them relative to the
cell. The constants computed with the atoms held at their strained positions are the clamped-ion constants;
the physical, relaxed-ion constants need the atoms relaxed in every strained cell, which `relax_ions=True`
does:

<!-- run -->
```python
hcp = bulk("Cu", "hcp", a=2.55, c=4.16)
hcp.calc = EMT()
BFGS(FrechetCellFilter(hcp), logfile=None).run(fmax=1e-6)

clamped = minelastic.calculate(hcp, EMT())
relaxed = minelastic.calculate(hcp, EMT(), relax_ions=True, fmax=1e-6)
for name in clamped.labels:
    print(f"{name}   clamped {clamped.as_dict()[name]:7.2f}   relaxed {relaxed.as_dict()[name]:7.2f}  GPa")
```

<!-- output -->
```text
c11   clamped  233.60   relaxed  216.34  GPa
c12   clamped   94.81   relaxed  112.07  GPa
c13   clamped   74.77   relaxed   74.77  GPa
c33   clamped  254.00   relaxed  254.00  GPa
c44   clamped   49.30   relaxed   49.30  GPa
```

The hexagonal crystal has five constants and needs two strain patterns, 8 strained cells in place of the 24
of the six single strains. The relaxation changes c11 and c12 by 17 GPa and leaves the other constants
unchanged.

The same steps with a DFT code are in the tutorials for [VASP](vasp.md), [Quantum ESPRESSO](qe.md),
[ABINIT](abinit.md) and [CP2K](cp2k.md).
