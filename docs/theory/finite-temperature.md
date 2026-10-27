# Finite temperature

At temperature $T$ the elastic constants are the strain derivatives of the thermally averaged stress,

$$
C_{ij}(T) = \left.\frac{\partial \langle \sigma_i \rangle_T}{\partial \varepsilon_j}\right|_{\varepsilon = 0},
$$

taken at the thermally expanded reference cell, with the cell held fixed in every strained state. These are
isothermal constants. The minimal strain sets apply unchanged; only the stress of each strained cell becomes
an average. The protocol of the paper is:

1. **Lattice at $T$.** An NPT run at zero pressure gives the thermally expanded cell. Most of the softening of
    the elastic constants with temperature comes from thermal expansion, so the expanded cell, not the 0 K
    one, is strained.
2. **Strained cells at $T$.** NVT runs with the cell held fixed (NPT would undo the strain), one in each
    strained cell of the minimal set.
3. **Constants.** For the states $x$ of pattern $k$ (amplitudes $m_x\delta$), with $n_x$ samples of the stress
    (MD steps after equilibration, block averages, or snapshots), the mean stresses $\bar\sigma_x$ give

    $$
    r_k = \frac{1}{\delta} \sum_x g_x\, \bar\sigma_x, \qquad
    \mathrm{Cov}(r_k) = \frac{1}{\delta^2} \sum_x g_x^2\, \frac{S_x}{n_x},
    $$

    with the finite-difference weights $g_x$ and the sample covariance matrix $S_x$ of the six stress components. The
    constants follow from the $r_k$ as at 0 K ([Elastic constants](elastic.md)), and their standard errors from
    $\mathrm{Cov}(r_k)$, correlations between the components of one pattern included. The formula assumes
    uncorrelated samples, so block-average or thin long trajectories.

MD inputs use the two-point formula ($\pm\delta$) by default. The four-point weights $(1, -8, 8, -1)/12$ amplify the noise of
the averages: at equal total MD time the result is about twice as noisy.

**Kinetic stress.** The Cauchy stress of an MD run contains the kinetic term $-n k_B T\, \delta_{ij}$ with
$n = N/V$. Under a normal strain $V$ changes by the factor $1 + \varepsilon_j$, so this term adds

$$
\frac{\partial (-n k_B T)}{\partial \varepsilon_j} = n k_B T \qquad (i, j \le 3)
$$

to $C_{ij}$ (about 0.25 GPa for Al at 300 K). VASP and pw.x print the stress of the instantaneous
configuration without it; `--kinetic-temperature T` (`kinetic_temperature=T`) adds the term. With
`stress_measure="kirchhoff"` the kinetic part of $J\sigma$ is $-N k_B T/V_0$, independent of the strain, and
nothing is added.

**Isothermal and adiabatic.** Fixed-cell NVT runs give isothermal constants; ultrasonic experiments measure
adiabatic ones, $C^S_{ij} = C^T_{ij} + T V \lambda_i \lambda_j / C_V$ with the thermal stress coefficients
$\lambda_i = \partial \sigma_i / \partial T|_\varepsilon$ and the heat capacity $C_V$ of the cell. For cubic crystals
$\lambda_1 = \lambda_2 = \lambda_3$, so only c11 and c12 differ (through the bulk modulus), not c44 or
(c11 - c12)/2.
