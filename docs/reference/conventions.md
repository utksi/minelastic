# Conventions

## Units and notation

| quantity | unit | convention |
|---|---|---|
| stress $\sigma$ | GPa | Cauchy stress, positive in tension; Voigt order xx, yy, zz, yz, xz, xy |
| strain $\varepsilon$ | dimensionless | Voigt order as above, engineering shear components ($\varepsilon_4 = 2\varepsilon_{yz}$) |
| elastic constants $C$ | GPa | 6x6 Voigt matrix, $\sigma_i = \sum_j C_{ij} \varepsilon_j$ |
| internal-strain tensor $\Lambda$ | eV/Å | $\Lambda[n, \alpha, J] = \partial F_{n\alpha} / \partial \varepsilon_J$ (per unit engineering strain) |
| force constants $\Phi$ | eV/Å$^2$ | $\Phi[m, \alpha, n, \beta] = -\partial F_{m\alpha} / \partial u_{n\beta}$ |
| Born charges $Z^*$ | \|e\| | $Z[n, \alpha, j] = \partial p_\alpha / \partial u_{nj}$ |
| piezoelectric stress tensor $e$ | C/m$^2$ | $e_{iJ} = \partial P_i / \partial \varepsilon_J$, rows x, y, z, Voigt columns |
| piezoelectric strain tensor $d$ | pC/N | $d = e S$, $S = C^{-1}$ the relaxed-ion compliance |
| dipoles (input to `compute`) | e Å | physical sign (positive charge) |

The code readers convert on reading. VASP prints the stress in kB, pressure-like (positive in compression),
in the order XX YY ZZ XY YZ ZX. pw.x prints Ry/bohr$^3$, also pressure-like. ASE uses eV/Å$^3$, positive in
tension. The VASP elastic tensors (`vasp-tensor`) are converted from kBar to GPa and reordered to Voigt
order.

## Strains

A strain is applied as the symmetric deformation gradient $F = I + \varepsilon$ (no rigid rotation) to the lattice
vectors, with the fractional coordinates fixed. Each pattern is applied at multiples of the step $\delta$ (default
0.01) given by the `stencil` option:

* `four-point` (default): $\pm\delta$, $\pm 2\delta$, weights $(1, -8, 8, -1)/12\delta$, error $O(\delta^4)$;
* `central` (the default for MD inputs, `generate --md`): $\pm\delta$, error $O(\delta^2)$, half the
  calculations;
* `forward`: $+\delta$ and the unstrained reference, error $O(\delta)$.

The four-point formula makes the result independent of the strain set to within $O(\delta^4)$, so the
minimal set and the six single strains agree closely. For averages of noisy stresses (MD) it amplifies
the noise; central differences are then the better choice.

The stress is differentiated as the Cauchy stress by default. `stress_measure="kirchhoff"`
(`--stress-measure kirchhoff`) differentiates $J\sigma$ instead, with $J = V/V_0$, as VASP does for `IBRION = 6`. For
a stress-free reference both give the same constants. They differ at $O(\delta^2)$ in the finite differences.

## Orientation

The tensor form, and with it the strain set, depends on how the crystal axes lie in the Cartesian frame.
minelastic derives the form from the point-group rotations of the given cell, so any orientation is handled
correctly. The textbook constants (c11, c12, c13, c33, c44, ...) refer to the standard orientation, though.
Examples of non-standard orientations are a hexagonal c axis along x, or a rhombohedral cell with its 3-fold
axis along [111]. With `frame="auto"` (the default), minelastic keeps the input orientation if the stiffness
matrix has its standard form there. Otherwise it rotates the cell rigidly into the spglib standard
orientation, the working frame.

* The constants, `result.C` and the strained cells refer to the working frame.
* `result.C_input_frame` (and `C_input_frame` in `elastic.json`) is the full tensor rotated back to the
  input orientation.
* `frame="input"` or `"standard"` forces one or the other.

For the trigonal space groups of type 312 (P312, P31m, P-31m and the others), the standard setting with a parallel to x has c15
instead of c14. minelastic handles both forms.

## Symmetry tolerance

The symmetry is found by spglib with `symprec = 1e-5` Å by default, the value VASP uses. Increase it
(`--symprec 1e-3`) for structures relaxed with loose tolerances. The space group, the elastic type and the
number of strain patterns are printed by `minelastic info`. If they are not what you expect, fix the
structure (or the tolerance) before running anything.
