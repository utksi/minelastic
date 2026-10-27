# Validation

`pytest tests` runs the test suite in about 20 seconds. Nearly all of it runs on data shipped in `tests/`,
so no DFT code is needed.

## Algorithm

* The paper's strain sets and extraction formulas are reproduced for all nine elastic types. Their
  minimality is confirmed by exhaustive search.
* Random symmetry-adapted tensors are recovered exactly from synthetic stresses with nonlinear terms. This
  holds in standard and random orientations, for every extraction method, with the central and
  four-point formulas.
* With ASE's EMT potential, the minimal and six-strain sets agree for clamped ions, relaxed ions and the
  internal-strain tensor.
* The symmetry completion of the internal-strain, piezoelectric and force-constant tensors is exact.
  The completed EMT force constants equal those from displacing every atom (1e-8 eV/Å$^2$).
* The internal-strain tensor of hcp Cu (EMT) from the two patterns of the minimal set equals brute force
  (all six single strains, all atoms) to 4e-5 eV/Å, and the tensor from the displaced atoms to 4e-6 eV/Å.
  The predicted relaxation per unit strain equals that of relaxing the atoms (`examples/emt_internal_strain.py`).
* For alpha-quartz with the BKS rigid-ion model, the relaxed-ion constants and the piezoelectric tensor from
  clamped single points (`ionic=True`) equal those from relaxing the ions in the strained cells to 0.1 GPa and
  2e-4 C/m$^2$, and the Born charges equal the model charges (`examples/bks_quartz_piezo.py`).
* The VASP and Quantum ESPRESSO writers and readers are tested. The pw.x lattice construction was also
  checked against pw.x's own output for every `ibrav` (1, 2, 3, -3, 4, 5, -5, 6 to 11, -9, 91, 12, -12, 13, -13 and 14, with
  `celldm` and with `A, B, C`).

## Real VASP data

From the raw $\pm\delta$ stresses of `IBRION = 6` runs (`tests/data`), minelastic reproduces the tensors that VASP
printed to 5e-5 GPa. This covers unmodified VASP and patched VASP builds that apply the minimal sets.
The minimal and six-strain sets agree for AlN (0.14 GPa), Se in two orientations (up to 0.07 GPa) and Al
(1.05 GPa with a 31x31x31 k-mesh).

## The paper's benchmark runs

`tests/data_paper` holds the stresses of the runs behind Table IV of the paper, one material per elastic
type, plus CaWO4, Bi and BaTiO3:

| elastic type | material |
|---|---|
| triclinic | albite |
| monoclinic | β-Ga2O3 |
| orthorhombic | TiSi2 |
| tetragonal I | CaMoO4 |
| tetragonal II | rutile TiO2 |
| trigonal I | dolomite |
| trigonal II | α-quartz |
| hexagonal | Mg |
| cubic | Al |

Each material has a run with VASP's six strains and a run with the paper's minimal set (a patched VASP). From the raw stresses, minelastic
reproduces the clamped-ion tensor VASP printed. The six-strain runs agree exactly (below 1e-3 GPa). The minimal runs
agree within 0.2 GPa, because the paper's implementation uses a different finite-difference estimator. With
VASP's ionic relaxation term added, the results reproduce the published Table IV:

* the eigenvalues agree within 0.2 GPa;
* the individual constants agree within 0.16 GPa where the run frame is the paper's. For α-quartz the
  agreement is up to the sign of c14, which depends on how the a axes are placed.

The minimal and six-strain sets agree within 0.6 GPa (Al: 1.05 GPa). The β-Ga2O3 runs are sibling runs of the
published ones, within about 2%.

## Ionic route and piezoelectricity

The reference is VASP's own `IBRION = 6` with `LEPSILON`, for wurtzite AlN (`tests/data_ionic`,
`tests/data_berry`):

* the ionic contribution to C agrees to 3e-4 GPa (of 77.7);
* with the DFPT Born charges, the ionic piezoelectric term agrees to 1e-5 C/m$^2$;
* the Berry-phase Born charges agree with DFPT to 0.004 |e|;
* the all-Berry-phase total e and d agree with the DFPT route to 0.006 C/m$^2$ (e33 = 1.454) and 0.05 pC/N
  (d33 = 5.30);
* the clamped-ion piezoelectric tensor from the Berry-phase dipoles of the 4 minimal strained cells agrees
  with DFPT to 0.005 C/m$^2$ (e33 = -0.424 C/m$^2$).

## End-to-end runs with real codes

Each code below was run on the inputs written by `minelastic generate`, and `minelastic collect` (or
`minelastic.calculate`) read the results. Agreement is with an independent reference computed by the same
code with the same settings.

| code | system | check | agreement |
|---|---|---|---|
| VASP 6.6.1 | wurtzite AlN | clamped ions, 4 vs 12 cells; relaxed ions, 4 relaxations vs `--ionic` vs VASP's `IBRION = 6` | 0.06 GPa; 0.7 GPa |
| VASP 6.6.1 | wurtzite AlN | `--ionic` with Berry phase: Born charges, e, d vs VASP's `IBRION = 6` + `LEPSILON` | 0.004 \|e\|, 0.006 C/m$^2$, 0.05 pC/N |
| VASP 6.6.1 | wurtzite AlN | internal-strain tensor from the forces of the 4 strained cells, and from the stresses of the displaced atoms, vs the tensor printed by VASP's `IBRION = 6` (largest entry 11.1 eV/Å) | 0.04 eV/Å; 0.06 eV/Å |
| VASP 6.6.1 | α-quartz (P3_221, 9 atoms) | `--ionic` with Berry phase vs VASP's `IBRION = 6` + `LEPSILON`: Born charges (not symmetric), C_ion, ionic piezoelectric term (with the DFPT Born charges), clamped-ion C | 0.011 \|e\|; 0.3 GPa of 230; 0.0002 C/m$^2$; 2.3 GPa of 310 (steps 0.01 vs 0.015) |
| VASP 6.6.1 | NiO, AFM-II (`MAGMOM`) | R-3m in the rhombohedral AFM-II cell: 4 cells vs VASP's `IBRION = 6` (12 cells), including c14 | 1.2 GPa (0.4 %) |
| VASP 6.6.1 | NiO, AFM-I, cubic 8-atom cell (`MAGMOM`) | only the moments lower Fm-3m to P4/mmm: 4 vs 12 cells | 1.6 GPa |
| VASP 6.6.1 | fcc Al, 32 atoms | `KSPACING` grid of minelastic vs VASP's own grid | identical |
| VASP 6.6.1 | fcc Al, 32 atoms, 300 K | `generate --md`, Langevin MD, `collect --md` vs an independent average of `vasprun.xml` (ASE) | identical |
| Quantum ESPRESSO 7.6 | wurtzite AlN | clamped ions, 2 vs 6 patterns; relaxed ions; internal strain | 0.05 GPa; 0.1 GPa; 0.007 eV/Å |
| Quantum ESPRESSO 7.6 | fcc Al, 32 atoms, 300 K | `generate --md`, Berendsen MD, `collect --md` vs ASE's reader of the pw.x output | identical |
| Quantum ESPRESSO 7.6 | wurtzite AlN | `--ionic` (12 single points) vs relaxing the ions | 1.1 GPa (0.3 %) |
| Quantum ESPRESSO 7.6 | NiO, AFM-I, cubic 8-atom cell (species labels Ni1/Ni2) | only the moments lower Fm-3m to P4/mmm: 4 vs 12 cells | 0.85 GPa |
| ABINIT 10.6.7 | wurtzite AlN | clamped ions, 4 vs 12 cells | 0.05 GPa |
| ABINIT 10.6.7 | wurtzite AlN | clamped ions and `--ionic` relaxed ions vs ABINIT's DFPT tensors (anaddb) | 0.6 GPa (0.15 %) |
| ABINIT 10.6.7 | NiO, AFM-I, cubic 8-atom cell (`spinat`) | only the moments lower Fm-3m to P4/mmm: 4 vs 12 cells | 1.5 GPa |
| MACE-MPA-0 (ASE) | wurtzite AlN | clamped ions, 2 vs 6 patterns; `ionic=True` vs `relax_ions=True` | 0.1 GPa; 0.3 GPa |
| CP2K 2026.2 (ASE) | diamond Si, 64 atoms, $\Gamma$ point | clamped ions, 1 vs 6 patterns; `ionic=True` gives c44 = 79.7 GPa against the clamped 102.6 GPa | 0.05 GPa |

The two MD runs check the finite-temperature pipeline only: their k-meshes (3x3x3 and 2x2x2 for 32 atoms)
are far too coarse for a metal, so the constants themselves are not converged; the static lattice at the
same settings gives the same values.
In the NiO AFM-I runs the cell is kept cubic (only its volume is relaxed), so a cubic form would be chosen
if the moments were not read; the tetragonal constants differ by up to 40 GPa (c11 - c33 = -33 GPa in VASP
and pw.x, -40 GPa in ABINIT). The 4- and 12-cell results differ by about the non-hydrostatic residual stress
of the cubic cell (up to 0.9 GPa).
The FHI-aims reader is tested on an FHI-aims output distributed with ASE; FHI-aims itself was not run. The
`generic` interface is tested with synthetic files only.
