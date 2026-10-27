# Examples

Runnable scripts (seconds to minutes on a laptop, with ASE's EMT potential, a rigid-ion model or stored data):

| script | what it shows |
|---|---|
| `emt_cu.py` | fcc and hcp Cu: clamped and relaxed ions, internal-strain tensor |
| `emt_ionic.py` | relaxed-ion constants from clamped cells (`ionic=True`) vs relaxing the ions, hcp Cu |
| `emt_internal_strain.py` | hcp Cu: internal-strain tensor from the minimal set and from displaced atoms, checked against brute force; ionic relaxation per unit strain |
| `bks_quartz_piezo.py` | alpha-quartz with the BKS rigid-ion model (a small ASE calculator in the script): Born charges, piezoelectric e and d, relaxed-ion constants, by two routes |
| `emt_finite_temperature.py` | fcc Al at 300 K: NPT for a(T), NVT in the strained cells, constants with errors (about 3 min) |
| `replay_vasp_aln.py` | the Python API on stored VASP data: elastic and Berry-phase piezoelectric tensors of AlN |
| `replay_vasp_aln_ionic.py` | the `--ionic` run of the docs from stored VASP data: relaxed-ion C, internal strain, Born charges, e and d of AlN against stock VASP |

Templates for DFT codes (wurtzite AlN, PBE):

| directory | code | notes |
|---|---|---|
| `vasp_aln/` | VASP | add a PBE POTCAR (Al, N); add `LCALCPOL = .TRUE.` for the piezoelectric route |
| `qe_aln/` | Quantum ESPRESSO | pslibrary 1.0.0 PAW potentials; set `pseudo_dir` to an absolute path |
| `abinit_aln/` | ABINIT | PseudoDojo v0.4 norm-conserving potentials; set `pp_dirpath` |

The documentation walks through them with the full output: https://utksi.github.io/minelastic/
