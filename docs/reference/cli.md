# Command line

<!-- input -->
```text
minelastic info | generate | collect | vasp-tensor   [--help for the options of each]
```

## `minelastic info`

Shows the symmetry, the tensor form and the strain set of a structure, without writing anything.

<!-- skip: usage examples; need a structure file and a template -->
```bash
minelastic info POSCAR                      # any file ASE can read (--format to force a format)
minelastic info --code vasp ref/            # template directory (vasp, aims) or pw.x input (qe)
```

## `minelastic generate CODE TEMPLATE`

Writes one input directory per strain state, plus `STATES` and `minelastic.json`.

| CODE | TEMPLATE |
|---|---|
| `vasp` | directory with `POSCAR`, `INCAR`, `POTCAR` and `KPOINTS` (or `KSPACING` in the INCAR) |
| `qe` | pw.x input file |
| `aims` | directory with `geometry.in` and `control.in` |
| `abinit` | directory with `POSCAR` and one input `*.abi` that reads it with `structure "poscar:POSCAR"` |
| `generic` | structure file; needs `--filename` (and `--format` unless ASE can tell the format from the file name) |

| option | meaning |
|---|---|
| `-o, --out DIR` | output directory (default `elastic`) |
| `--relax-ions` | relax the ions in each strained cell (relaxed-ion constants) |
| `--ionic` | also displace the symmetry-irreducible atoms ([relaxed-ion and piezoelectric response from clamped cells](../how-to/relaxed-ions.md)) |
| `--displacement A` | `--ionic`: atomic step in Å (default 0.01) |
| `--md` | fixed-cell molecular dynamics in each strained cell, as set up in the template ([finite temperature](../how-to/finite-temperature.md)); vasp and qe |
| `--scheme minimal\|asess` | minimal uncoupled set (default) or the six single-component strains |
| `--delta D` | strain step (default 0.01) |
| `--stencil four-point\|central\|forward` | finite-difference formula (default four-point, $\pm\delta$ and $\pm 2\delta$; central with `--md`) |
| `--symprec S` | spglib tolerance in Å (default 1e-5) |
| `--frame auto\|input\|standard` | orientation of the applied strains ([Conventions](conventions.md#orientation)) |
| `--potcar copy\|link` | vasp: copy or symlink the POTCAR |
| `--fix-incar` | vasp: apply all recommended settings the INCAR does not meet (`EDIFF`, `ENCUT`, `PREC`, `LREAL`); without it they are only reported ([The INCAR](../codes/vasp.md#the-incar)) |
| `--format FMT`, `--filename NAME`, `--extra FILES` | generic: ASE format, structure file name, further files to copy |

## `minelastic collect [DIR]`

Reads the finished calculations in `DIR` (default `elastic`), prints the results and writes
`DIR/elastic.json`.

| option | meaning |
|---|---|
| `--method uncoupled\|designated\|lstsq` | extraction method ([Method](../theory/index.md#extraction-and-consistency-checks)) |
| `--output NAME` | output file in each directory (defaults: `vasprun.xml`, then `OUTCAR`; `pw.out`; `aims.out`; newest `*.abo`) |
| `--format FMT` | generic: ASE format of the output file |
| `--internal-strain` | print the internal-strain tensor (clamped ions) |
| `--piezo` | piezoelectric tensor (and Born charges with `--ionic`) from the Berry-phase dipoles (VASP `LCALCPOL`) |
| `--born OUTCAR` | Born charges and clamped-ion piezoelectric tensor from a VASP `LEPSILON`/`LCALCEPS` run |
| `--md`, `--skip N` | average the stress over the MD steps after the first N |
| `--samples GLOB` | average the stresses of per-snapshot sub-directories (e.g. `'conf-*'`) |
| `--kinetic-temperature T` | add the ideal-gas term n k_B T (stresses without the kinetic part) |
| `--stress-measure cauchy\|kirchhoff` | differentiate $\sigma$ (default) or $J\sigma$ (as VASP `IBRION = 6`) |

## `minelastic vasp-tensor OUTCAR`

Converts the tensor printed by VASP `IBRION = 6, ISIF >= 3` to GPa and Voigt order. `--which
total|clamped|ionic` selects TOTAL ELASTIC MODULI (default), SYMMETRIZED ELASTIC MODULI or the ionic
relaxation contribution. For the first two it also prints the Born stability and, for a stable tensor, the
Voigt-Reuss-Hill moduli.

## Output files

| file | content |
|---|---|
| `STATES` | the strain-state directory names, one per line (for job arrays) |
| `minelastic.json` | the setup: program name, version and URL, reference structure, frame rotation, patterns, step, finite-difference multipliers, states, code, template; `collect` needs it |
| `elastic.json` | the results, described below |

`elastic.json` holds:

* always: `program`, `version` and `url` of minelastic, `constants` (GPa), the matrices `C` (working frame) and `C_input_frame`, `checks` (every
  stress component with its computed and fitted value), `properties` (VRH moduli, Young's modulus,
  Poisson's ratio, anisotropy index; `null` if C is not positive definite), `spacegroup`, `elastic_type`
  and the collect options;
* `errors`: standard errors, for MD or snapshot averages with at least two samples per state;
* `internal_strain_eV_per_A` (N, 3, 6): with the forces of clamped cells;
* `piezoelectric_C_per_m2` (3, 6): with the dipoles of the strained cells;
* with `--ionic`: `force_constants_eV_per_A2`, `internal_strain_displaced_eV_per_A`, `C_ionic`,
  `C_relaxed`, `ionic_relaxation_A_per_strain`, and `properties_relaxed` if the relaxed-ion C is positive
  definite. With dipoles or `--born` also `born_charges`, `born_source`, `piezoelectric_ionic_C_per_m2`,
  `piezoelectric_total_C_per_m2` and `piezoelectric_d_pC_per_N`.

The same tensors rotated to the orientation of the input structure (identical to the working frame when the
cell was not rotated) are `C_input_frame`, `C_relaxed_input_frame`, `C_ionic_input_frame`,
`internal_strain_input_frame`, `internal_strain_displaced_input_frame`, `force_constants_input_frame`,
`ionic_relaxation_input_frame`, `born_charges_input_frame`, `piezoelectric_input_frame`,
`piezoelectric_ionic_input_frame`, `piezoelectric_total_input_frame` and `piezoelectric_d_input_frame`.
