# Citing

If you use minelastic in published work, please cite the paper that introduces the method:

> S. I. Simak and U. Singh, *A general algorithm for minimal uncoupled stress-strain calculations of
> elastic constants*, [arXiv:2610.11733](https://arxiv.org/abs/2610.11733) (2026).

```bibtex
@misc{simak2026minimal,
  title         = {A General Algorithm for Minimal Uncoupled Stress--Strain Calculations of Elastic Constants},
  author        = {Simak, Sergei I. and Singh, Utkarsh},
  year          = {2026},
  eprint        = {2610.11733},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.mtrl-sci},
  doi           = {10.48550/arXiv.2610.11733},
  url           = {https://arxiv.org/abs/2610.11733}
}
```

The repository also contains a [`CITATION.cff`](https://github.com/utksi/minelastic/blob/main/CITATION.cff)
file. GitHub shows it as "Cite this repository".

If you use the relaxed-ion or piezoelectric route (`--ionic`), the Berry-phase piezoelectric tensor, or the
finite-temperature averaging, please also cite the original methods:

* X. Wu, D. Vanderbilt and D. R. Hamann, Phys. Rev. B **72**, 035105 (2005): the ionic contributions from
  force constants, internal strain and Born charges;
* D. Vanderbilt, J. Phys. Chem. Solids **61**, 147 (2000): the proper piezoelectric response.

minelastic uses [spglib](https://spglib.readthedocs.io) for the symmetry analysis and
[ASE](https://wiki.fysik.dtu.dk/ase/) for structures and calculators.
