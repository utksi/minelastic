# Cost of the minimal sets

The number of strain patterns in the minimal set equals the number of lattice parameters that the symmetry
leaves free, and no stress-strain scheme can use fewer
([proof sketch](index.md#why-fewer-cells-are-impossible)). The usual approach, used by VASP's
`IBRION = 6` and by most workflows, applies the six single-component strains. Both apply every pattern at
the same step sizes: four by default in minelastic ($\pm\delta$, $\pm 2\delta$, as VASP's `NFREE = 4`) and
two with `--stencil central` (as `NFREE = 2`). The saving is therefore the ratio of the pattern counts.

| elastic type | strain patterns | single points (four-point / central) | six single strains (four-point / central) | fewer by |
|---|---|---|---|---|
| triclinic | 6 | 24 / 12 | 24 / 12 | 1 |
| monoclinic | 4 | 16 / 8 | 24 / 12 | 1.5 |
| orthorhombic | 3 | 12 / 6 | 24 / 12 | 2 |
| tetragonal I and II | 2 | 8 / 4 | 24 / 12 | 3 |
| trigonal I and II | 2 | 8 / 4 | 24 / 12 | 3 |
| hexagonal | 2 | 8 / 4 | 24 / 12 | 3 |
| cubic | 1 | 4 / 2 | 24 / 12 | 6 |

Most materials of interest are in the last four rows. The paper also compares the minimal sets with earlier
tabulated schemes, which need more configurations for several symmetry types.

## What each quantity costs

Everything that minelastic computes is built on the same strained cells, so every quantity inherits the
factor above. $n$ is the number of patterns, $p$ the number of steps per pattern (4 or 2) and $D$ the number
of symmetry-irreducible atomic displacements.

| quantity | minimal set | six single strains |
|---|---|---|
| clamped-ion constants | $n p$ single points | $6 p$ single points |
| internal-strain tensor | no extra cost: the forces of the same cells | $D p$ displaced-atom single points (VASP) |
| clamped-ion piezoelectric tensor $e^0$ | no extra cost: Berry phase in the same cells | a DFPT run, or $6 p$ Berry-phase cells |
| relaxed-ion constants by relaxing the ions | $n p$ relaxations | $6 p$ relaxations |
| relaxed-ion constants, Born charges, total $e$ and $d$ (`--ionic`) | $n p + D p$ single points | $6 p + D p$ single points (VASP `IBRION = 6`) |
| finite-temperature constants | $n p$ MD runs ($p = 2$) | $6 p$ MD runs |

For wurtzite AlN ($n = 2$, $D = 4$) and two points per strain, `--ionic` needs 4 + 8 = 12 single points
where VASP's `IBRION = 6` needs 12 + 8 = 20. For fcc Al at finite temperature it is 2 MD runs instead of
12, the largest saving in absolute terms, since every cell is a full molecular-dynamics run.

## Measured

Wall times of the clamped-ion runs of wurtzite AlN in the tutorials ([VASP](../get-started/vasp.md),
[Quantum ESPRESSO](../get-started/qe.md), [ABINIT](../get-started/abinit.md)), with the same
settings for both sets (two points per strain, one 32-core node):

| code | minimal set (4 cells) | six single strains (12 cells) | ratio |
|---|---|---|---|
| VASP 6.6.1 | 42 s | 163 s | 3.9 |
| Quantum ESPRESSO 7.6 | 9 s | 27 s | 3.0 |
| Quantum ESPRESSO 7.6, relaxed ions | 22 s | 65 s | 3.0 |
| ABINIT 10.6.7 | 373 s | 1220 s | 3.3 |

The constants of the two sets agree to 0.05 to 0.1 GPa ([Validation](../reference/validation.md)). In VASP the ratio
exceeds the ratio of the cell counts because the cells strained along zz keep the hexagonal symmetry (30
irreducible k points against 93 to 102 for the sheared cells, 6 s against 15 to 18 s per cell). They are half
of the minimal set but only a sixth of the six-strain set.
