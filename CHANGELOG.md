# Changelog

## 0.1.0 (2026-10-22)

First public release, accompanying S. I. Simak and U. Singh, arXiv:2610.11733.

* Minimal uncoupled strain sets for all nine elastic symmetry types, in any orientation; extraction
  methods `uncoupled`, `designated` and `lstsq`; consistency checks; VRH moduli and Born stability.
* Four-point finite differences by default (central differences for MD inputs).
* `generate`/`collect` from the command line or from Python (`minelastic.generate`, `minelastic.collect`).
  For VASP, every INCAR change is listed and the settings that matter for small stress differences are
  checked; `--fix-incar` applies them.
* `calculate(..., dipole=f)`: Born charges and piezoelectric tensors with any ASE calculator.
* Interfaces: VASP, Quantum ESPRESSO (pw.x, any `ibrav`), ABINIT, FHI-aims, any ASE format, any ASE
  calculator.
* Clamped- and relaxed-ion constants; internal-strain tensor completed by symmetry.
* `--ionic`: Gamma-point force constants, Born charges, relaxed-ion stiffness and the piezoelectric e and
  d tensors from clamped single points.
* Berry-phase piezoelectric tensor from the strained cells (VASP `LCALCPOL`).
* Finite temperature: MD averages (`generate --md`, `collect --md`) and snapshot averages (`--samples`),
  with standard errors and the kinetic term.
* Symmetry with magnetic order (magnetic space group, including time-reversed operations) from VASP
  `MAGMOM`, ABINIT `spinat`, pw.x species labels and noncollinear directions, and ASE initial moments; VASP
  species blocks and pw.x species labels.
* Documentation with tutorials for ASE, VASP, Quantum ESPRESSO, ABINIT and CP2K (outputs of real runs),
  how-to guides per quantity and per code, the theory, and the API reference generated from the
  docstrings. The tests run every code block of the documentation and compare the printed output.
