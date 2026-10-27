# Troubleshooting

## `info` shows more strain patterns than expected

The number of patterns follows from the symmetry that spglib finds with the tolerance `symprec` (default
1e-5 Å). A slightly distorted structure is found with a lower symmetry: symmetrize it, or pass a larger
`--symprec`. Magnetic order can also lower the symmetry ([Magnetic materials](magnetic.md)), as can atoms of
one element in different VASP species blocks or with different pw.x species labels, which minelastic treats
as different species, as the DFT code does. spglib lists only the operations that map the lattice of the
given cell onto itself, so a supercell whose lattice is less symmetric than the crystal gets the set of that
lower symmetry. In all these cases the constants are still correct, but more cells are computed than
necessary.

## The consistency check reports a large deviation

The stress components that are not used for the constants disagree with the fitted tensor. Possible causes
are a structure of lower symmetry than assumed (see above), a strain step outside the linear regime, too loose
electronic or ionic convergence, a k-point mesh that differs between cells, or the noise of MD averages.
`collect --method lstsq` shows how much the constants depend on which stress components are used.

## The constants differ from VASP's `IBRION = 6`

By a few tenths of a GPa, because the numerical derivatives differ: VASP strains the cell by 0.015, applies
simple shears and differentiates the stress times the volume ratio
([Elastic constants in practice](elastic.md#compare-with-vasps-ibrion-6)). Differences of several GPa point to
different settings: VASP's `IBRION = 6` run may use another cutoff, k-point mesh or `PREC` than the template.

## `collect` refuses to read MD runs

Without `--md`, minelastic reads the stress of the last ionic step only. So that the last step of an MD run
is never taken by mistake, `collect` stops if the inputs were written with `generate --md`; collect them with
`--md --skip N` (or `--samples` for snapshots).

## `--ionic` stops with an unstable or soft mode

The force constants have a negative or nearly zero eigenvalue besides the three uniform translations. The
relaxed-ion response is then undefined: the reference structure is not at a minimum of the energy, or a
phonon at $\Gamma$ is unstable. Relax the atoms of the reference structure tightly, and check the structure.

## Noncollinear inputs are refused

The cell had to be rotated into the standard orientation of the crystal, and the moment directions (and
other direction-dependent inputs such as `SAXIS`, `EFIELD` or Cartesian k points) would have to be rotated
with it. Use `--frame input` to keep the input orientation.

## Can I combine cells from different codes or settings?

No. All strained cells of a set must be computed with identical settings. The constants are estimated from
stress differences, and systematic differences between cells do not cancel.
