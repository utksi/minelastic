# Quantities beyond the clamped-ion constants

## Relaxed-ion constants

Atoms that are not fixed by symmetry move when the crystal is strained, and the measured (relaxed-ion)
constants are lower than the clamped-ion ones:

    C = C0 - Lambda^T Phi^+ Lambda / Omega

with the internal-strain tensor Lambda (force per unit strain, eV/Å), the Gamma-point force constants
Phi (eV/Å^2) and the cell volume Omega. Two routes:

| | relax the ions | displace the atoms |
|---|---|---|
| option | `--relax-ions`, `relax_ions=True` | `--ionic`, `ionic=True` |
| calculations | one relaxation per strained cell | one single point per strained cell, plus one per symmetry-irreducible atomic displacement and step |
| also gives | | force constants and, with dipoles, Born charges and piezoelectric tensors |
| watch out for | tight relaxations: a strain of 0.01 induces forces of only about 0.05 eV/Å | an unstable or soft optical mode makes the response undefined; minelastic stops |

The second route is the method of VASP's `IBRION = 6`, needs no relaxations, and is usually cheaper and
more robust. With `--ionic`, `--displacement` is the atomic step (default 0.01 Å; VASP uses 0.015).
The relaxed-ion results are `C_relaxed` (`elastic.json`, `result.C_relaxed`); `constants` and
`result.as_dict()` stay clamped-ion.

## Internal strain

The forces in the strained cells give the internal-strain tensor without further calculations (shape
(N, 3, 6), eV/Å): `collect --internal-strain`, `result.internal_strain`, `internal_strain_eV_per_A` in
`elastic.json`. With `--ionic`, the displaced cells give it a second time (`internal_strain_displaced`)
together with the force constants (`force_constants`, (N, 3, N, 3)).

## Piezoelectricity

    e = e0 + Z* Phi^+ Lambda / Omega,    d = e S,    S = C^-1

e0 is the clamped-ion tensor (C/m^2) from the dipoles of the strained cells, Z* the Born charges from the
dipoles of the displaced cells, and d (pC/N) uses the relaxed-ion compliance.

* VASP: add `LCALCPOL = .TRUE.` to the template INCAR (generate does not add it), generate with `--ionic`
  for Z*, e and d (without it, only e0), and `collect --piezo`. `elastic.json` then holds `born_charges`,
  `piezoelectric_C_per_m2` (e0), `piezoelectric_ionic_C_per_m2`, `piezoelectric_total_C_per_m2` (e) and
  `piezoelectric_d_pC_per_N` (d). The acoustic sum of the Born charges, printed, shows how well they are
  converged. Alternatively `--born OUTCAR` takes Z* and e0 from a `LEPSILON` or `LCALCEPS` run.
* ASE calculators: `calculate(..., ionic=True, dipole=f)` with `f(atoms)` the dipole in e Å (physical
  sign); results `result.born`, `result.piezo` (e0), `result.piezo_total` (e), `result.piezo_d` (d).
* Other codes: minelastic reads Berry-phase dipoles from VASP only. Pass the dipoles to
  `ElasticSetup.compute(stresses, forces=..., polarization={state_id: dipole})`, or Z* and e0 from a
  perturbation-theory run with `born=` and `clamped_piezo=`.

## Finite temperature

The constants at temperature T are averages over MD in each strained cell, at fixed cell, at the volume of
temperature T (from an NPT run or the thermal expansion).

* VASP or pw.x: the template is the NVT run of the thermally expanded cell (VASP: `IBRION = 0`, `NSW`,
  `POTIM`, a thermostat such as `MDALGO = 3` with `LANGEVIN_GAMMA` per species, `TEBEG = TEEND = T`).
  `generate CODE TEMPLATE -o el_T --md` keeps the MD settings and fixes the cell (default stencil central).
  After the runs: `minelastic collect el_T --md --skip N --kinetic-temperature T`, where N is the number of
  equilibration steps to discard (for example the first 1 to 2 ps), and the kinetic term n k_B T is added
  because VASP and pw.x print the stress without it (about 0.25 GPa for Al at 300 K).
* Snapshots: single points on MD snapshots (for example from a machine-learned potential), in sub-directories
  `conf-*` of each cell: `collect el_T --samples 'conf-*'`.
* ASE: run the MD yourself in each cell of `ElasticSetup(atoms, stencil="central")` and pass the averaged
  stresses to `setup.compute`, or see `examples/emt_finite_temperature.py`.
* Standard errors (`errors`) assume uncorrelated steps; thin or block-average the trajectory for reliable
  error bars. A smaller strain step makes the constants noisier.

## Magnetic materials

minelastic uses the magnetic space group (spglib), including the operations combined with time reversal.
Magnetic order can lower the symmetry (AFM-I NiO is tetragonal, 2 patterns and 6 constants, instead of
cubic). The moments come from `MAGMOM` (VASP), `spinat` (ABINIT), the species labels and `angle1`/`angle2`
(pw.x), `initial_moment` (FHI-aims) or `atoms.get_initial_magnetic_moments()` (ASE). Check the space group
that `minelastic info --code CODE TEMPLATE` prints. With an ASE calculator the moments must be the ones the
calculator uses; EMT and MACE ignore them.

## Your own stresses

<!-- run -->
```python
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.units import GPa
from minelastic import ElasticSetup

setup = ElasticSetup(bulk("Cu", "fcc", a=3.5898, cubic=True))
stresses = {}
for state in setup.states:                       # state["id"] = 's1_m2', 's1_m1', ...
    cell = setup.strained_atoms(state)           # ASE Atoms of this cell, possibly rotated
    cell.calc = EMT()                            # stand-in for any calculation
    stresses[state["id"]] = cell.get_stress(voigt=True) / GPa   # GPa, tension positive, Voigt
result = setup.compute(stresses)
print(result.as_dict())
```

Compute the stress of each cell exactly as `strained_atoms` returns it. `compute` also takes
`forces={state_id: (N, 3) array}` (eV/Å) and `polarization={state_id: (3,) array}` (e Å).
`setup.to_dict()` and `ElasticSetup.from_dict(d)` save and restore a setup.
