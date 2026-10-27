# Elastic constants in practice

The tutorials go through the basic workflow. Here we discuss the choices that determine the accuracy and
the cost of the constants; the formulas are in the theory section, under
[Elastic constants](../theory/elastic.md).

## Prepare the reference structure

The elastic constants are derivatives of the stress at the stress-free equilibrium structure. Relax the
reference structure with the same settings (functional, cutoff, k-point mesh, smearing) that the strained
cells will use. The stresses of the strained cells differ from it by about $\delta C$, 1 to 5 GPa for the
default step $\delta = 0.01$, so a residual stress of 0.1 GPa is harmless; it cancels in the difference to
first order.

The same reasoning sets the accuracy the calculations need. In VASP, raise `ENCUT` to 1.3 to 1.4 times the
largest `ENMAX` of the POTCAR, so that the Pulay stress changes little with the cell volume, and converge the
electronic structure tightly (`EDIFF` of 1E-8). `minelastic generate vasp` reports the settings that the
template does not meet ([VASP](../codes/vasp.md#the-incar)). All strained cells use the k-point mesh of the
template; a `KSPACING` in the INCAR is replaced by the explicit mesh of the reference cell.

!!! warning "One set of settings"

    Never combine cells computed with different codes or settings. The constants are differences of
    stresses, and a systematic difference between two cells does not cancel.

## Check the symmetry before running

`minelastic info` prints the space group, the elastic type and the number of strain patterns, without writing
anything:

<!-- shell -->
```console
$ minelastic info --code vasp vasp_aln
space group P6_3mc (186), point group 6mm, elastic type hexagonal
...
strain set 'minimal': 2 pattern(s), step 0.01, multipliers (-2, -1, 1, 2)
...
```

If the space group is lower than expected, the structure is slightly distorted (symmetrize it, or pass a
larger tolerance with `--symprec`, default 1e-5 Å), or magnetic order or different species labels lower the
symmetry ([Magnetic materials](magnetic.md)). The constants are still correct, but more cells are computed
than necessary. A supercell whose lattice is less symmetric than the crystal also gets the strain set of the
lower symmetry; start from the primitive or conventional cell.

## Choose the strain step and the number of points

Each strain pattern is applied at several multiples of the step $\delta$, and the stresses are
differentiated numerically. minelastic offers three finite-difference formulas, chosen with `--stencil`
(`stencil=` in Python):

| option | strains | error | calculations per pattern |
|---|---|---|---|
| `four-point` (default) | $\pm\delta$, $\pm 2\delta$ | of order $\delta^4$ | 4 |
| `central` (default with `--md`) | $\pm\delta$ | of order $\delta^2$ | 2 |
| `forward` | $+\delta$ and the unstrained cell | of order $\delta$ | 1 plus the reference |

The four-point formula is the best choice for static calculations, where the stresses are precise and the
error of the numerical derivative dominates. For hcp Cu (EMT), the two-point formula is off by 0.49 GPa at
$\delta = 0.01$. For wurtzite AlN in VASP, the minimal set and the six single strains differ by 0.14 GPa with
two points and by 0.005 GPa with four. Averages over molecular dynamics are noisy, and the four-point weights
$(1, -8, 8, -1)/12$ amplify that noise, so MD inputs use two points
([Finite temperature](finite-temperature.md)).

The step `--delta` (default 0.01) should stay in the linear regime: 0.005 to 0.015 usually is. Smaller steps
amplify numerical noise, larger ones add anharmonic contributions. For soft or strongly anharmonic materials,
compare two steps.

## Check the result

The last line of the summary compares the stress components that are not needed for the constants with the
values that the fitted tensor predicts:

<!-- input -->
```text
consistency (redundant stress components): largest |computed - fitted| = 0.00 GPa (pattern 1, xz)
```

A deviation of more than a few tenths of a GPa points to a reference structure that is not relaxed, a
strain step outside the linear regime, loose convergence, a k-point mesh that differs between cells, or a
lower symmetry than assumed.

Two further checks cost more calculations. `--scheme asess` applies the six single strains, the usual
approach, and should give the same constants within the error of the numerical derivative. `collect --method
lstsq` fits all stress components of all patterns by least squares instead of reading each constant from the
components that contain it alone (the default, `uncoupled`); a large change shows that some components are
inconsistent.

## Orientation of the result

The textbook constants (c11, c12, c13, c33, c44, ...) refer to the standard orientation of the crystal, for
example the hexagonal c axis along z. If the stiffness matrix does not have its standard form in the input
orientation, minelastic rotates the cell rigidly into the standard orientation before straining it, and
reports the constants there. The full tensor in the input orientation is `C_input_frame`, in `elastic.json`
and as `result.C_input_frame`. `--frame input` keeps the input orientation and `--frame standard` always
rotates ([Conventions](../reference/conventions.md#orientation)).

## Compare with VASP's `IBRION = 6`

`minelastic vasp-tensor OUTCAR` reads the elastic tensor that VASP prints for `IBRION = 6` and `ISIF = 3`,
converts it to GPa and Voigt order and adds the stability test and the moduli:

<!-- skip: needs an OUTCAR of an IBRION = 6 run -->
```bash
minelastic vasp-tensor OUTCAR                     # TOTAL ELASTIC MODULI (relaxed ions)
minelastic vasp-tensor OUTCAR --which clamped     # SYMMETRIZED ELASTIC MODULI (clamped ions)
```

Its numbers differ from minelastic's by a few tenths of a GPa. VASP strains the cell by 0.015 (`POTIM`),
applies simple shears instead of symmetric strains and differentiates $J\sigma$, the stress times the volume
ratio, instead of the stress itself. `--delta 0.015 --stencil central --stress-measure kirchhoff --scheme
asess` comes closest to what VASP does with its default `NFREE = 2`.

## Use other strain sets

Any set of patterns can be passed to `ElasticSetup` as `scheme=[[...6 numbers...], ...]`, for example the
strain sets tabulated by Le Page and Saxe. The constants are then fitted by least squares
(`method="lstsq"`), since such sets are not uncoupled in general.
