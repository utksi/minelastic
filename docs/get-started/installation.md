# Installation

minelastic requires Python 3.9 or newer. Its dependencies, numpy, [spglib](https://spglib.readthedocs.io)
and [ASE](https://wiki.fysik.dtu.dk/ase/), are installed with it:

<!-- skip: installs the package -->
```bash
pip install git+https://github.com/utksi/minelastic.git
```

The command-line tool is called `minelastic`. To check the installation:

<!-- shell -->
```console
$ minelastic --version
0.1.0
```

`minelastic --help` lists the commands, and `minelastic generate --help` the options of one of them.

## From a clone

The repository also contains the example scripts, the input templates of the tutorials and the tests:

<!-- skip: clones the repository and installs the package -->
```bash
git clone https://github.com/utksi/minelastic.git
cd minelastic
pip install -e ".[test]"
pytest tests
```

The tests take about a minute and need no DFT code: they use stored stresses of real calculations.

## Optional packages

The DFT codes are not needed to install minelastic, to write the inputs or to read the results; they are
needed only to run the calculations. Machine-learned potentials are used through their ASE calculators, for
example `pip install mace-torch` for [MACE](https://github.com/ACEsuit/mace). The tutorial with ASE's EMT
potential relaxes the cell with `ase.filters.FrechetCellFilter`, which needs ASE 3.23 or newer. The
documentation is built with `pip install -e ".[docs]"` followed by `mkdocs serve`.
