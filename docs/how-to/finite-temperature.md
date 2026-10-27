# Finite temperature

At a temperature $T$ the elastic constants are the strain derivatives of the thermally averaged stress,

$$
C_{ij}(T) = \frac{\partial \langle \sigma_i \rangle_T}{\partial \varepsilon_j} ,
$$

taken at the thermally expanded lattice ([theory](../theory/finite-temperature.md)). The strain sets are the
same as at 0 K; only the stress of each strained cell becomes an average over a molecular-dynamics run. The
saving is largest here, because every strained cell is a full MD run: a cubic crystal needs 2 runs instead
of 12.

The protocol has three steps:

1. an NPT run at zero pressure gives the lattice parameters at $T$; most of the softening of the constants
   with temperature comes from thermal expansion;
2. an NVT run with the cell held fixed in every strained cell of the expanded lattice (NPT would undo the
   strain);
3. the constants from the averaged stresses, with standard errors from the fluctuations.

The stress fluctuations scale as $1/\sqrt{N}$ with the number of atoms. A larger supercell, a longer run or
a larger strain step (within the linear regime) reduce the error. MD inputs use two points per pattern,
$\pm\delta$, because the four-point formula amplifies the noise of the averages: at equal total MD time its
result is about twice as noisy. For metals, use an electronic smearing that matches the temperature (VASP
`ISMEAR = -1` with `SIGMA` = $k_B T$).

## With an ASE calculator

[`examples/emt_finite_temperature.py`](https://github.com/utksi/minelastic/tree/main/examples/emt_finite_temperature.py)
follows the three steps for fcc Al with ASE's EMT potential in a 108-atom supercell. Berendsen NPT gives
$a(T)$. Langevin NVT then runs for 7000 steps of 2 fs (1000 of them for equilibration) in each of the two
strained cells. The stresses include the kinetic part (`include_ideal_gas=True`) and are passed to
`ElasticSetup.compute` as block averages, an array of samples per strain state:

<!-- excerpt: examples/emt_finite_temperature.py -->
```python
setup = minelastic.ElasticSetup(bulk("Al", "fcc", a=a_T, cubic=True).repeat(3), delta=0.01, stencil="central")
stresses = {}
for st in setup.states:
    atoms = setup.strained_atoms(st)
    ...
    stresses[st["id"]] = s.reshape(-1, BLOCK, 6).mean(axis=1)
res = setup.compute(stresses)
```

<!-- output-of: examples/emt_finite_temperature.py -->
```text
a(0 K) = 3.9943 A, a(300 K) = 4.0269 A
s1_m1: 6000 steps, <stress> = [-0.509 -0.375 -0.372 -0.303 -0.003 -0.003] GPa
s1_p1: 6000 steps, <stress> = [ 0.425  0.259  0.253  0.303  0.002 -0.009] GPa

            0 K      300 K
c11        53.3     46.7 +- 0.4  GPa
c12        32.9     31.5 +- 0.5  GPa
c44        36.2     30.3 +- 0.4  GPa
```

The run takes about 3 minutes on one core. EMT is a toy potential for Al, so only the workflow carries over,
not the numbers.

## With VASP or Quantum ESPRESSO

The template is the NVT run of the thermally expanded supercell. `generate --md` keeps its MD settings and
fixes the cell. For VASP, for example Langevin NVT at 600 K:

<!-- input -->
```text
IBRION = 0 ; NSW = 10000 ; POTIM = 2
MDALGO = 3 ; LANGEVIN_GAMMA = 10 ; TEBEG = 600 ; TEEND = 600   # one LANGEVIN_GAMMA per species
ISMEAR = -1 ; SIGMA = 0.0517              # k_B T at 600 K (eV), for a metal
```

<!-- skip: needs a VASP MD template and the MD runs -->
```bash
minelastic generate vasp ref_600K -o el_600K --md      # sets ISIF = 2 (fixed cell), keeps IBRION, NSW, POTIM
# run the MD in every directory of el_600K/STATES
minelastic collect el_600K --md --skip 1000 --kinetic-temperature 600
```

`--skip` discards the equilibration steps. VASP and pw.x print the stress without its kinetic part, and
`--kinetic-temperature` adds it ($n k_B T$ in $C_{11}$, $C_{12}$, ..., about 0.25 GPa for Al at 300 K). The
standard errors assume uncorrelated steps, so thin or block-average the trajectory if you need reliable
error bars. For pw.x, the template has `calculation = 'md'` with a thermostat in `&IONS` and a fixed cell,
and the rest is the same. minelastic refuses to collect MD runs without `--md` (or `--samples`), so that
the last step is never taken by mistake.

A cheaper route runs the MD with a machine-learned force field and computes accurate single points on
selected snapshots, in sub-directories `conf-*` of each strain state. Collect them with

<!-- skip: needs the snapshot calculations -->
```bash
minelastic collect el_600K --samples 'conf-*'
```

This works with every interface. The paper uses this protocol for fcc Al from 300 to 900 K: ab initio MD in
a 108-atom supercell, accelerated by on-the-fly machine-learned force fields.

The constants from fixed-cell NVT runs are isothermal. Ultrasonic experiments measure adiabatic constants,
which differ in c11 and c12 of a cubic crystal (through the bulk modulus) but not in c44
([theory](../theory/finite-temperature.md)).
