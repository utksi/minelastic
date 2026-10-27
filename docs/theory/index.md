# The method

Here we explain the method of the paper ([arXiv:2610.11733](https://arxiv.org/abs/2610.11733)) and give
the formulas for every quantity that minelastic computes. None of this is needed to use minelastic; the
tutorials and how-to guides can be followed without it.

We first show why a few strained cells suffice, how the strain patterns are chosen and why no stress-strain
scheme can use fewer. Definitions and formulas for each quantity are on separate pages:
[elastic constants](elastic.md), [relaxed-ion constants](relaxed-ions.md),
[internal-strain and force-constant tensors](internal-strain.md), [piezoelectric tensors](piezoelectricity.md),
[finite-temperature constants](finite-temperature.md) and [cost](cost.md).

## Stress-strain finite differences

A homogeneous strain $\varepsilon$ (Voigt vector with engineering shear components) is applied to the
lattice vectors of the relaxed reference cell, and the stress of the strained cell is computed. To first
order, $\sigma_i = \sum_j C_{ij}\,\varepsilon_j$. A strain *pattern* is a direction $u$ in strain space,
applied at $\pm\delta$. The central difference

$$
\frac{\sigma(+\delta u) - \sigma(-\delta u)}{2\delta} = C\,u + O(\delta^2)
$$

gives six linear relations between the stress components and the elastic constants. The default
four-point formula, $[\sigma(-2\delta u) - 8\sigma(-\delta u) + 8\sigma(\delta u) - \sigma(2\delta u)]/12\delta$,
gives the same relations with an error of order $\delta^4$. The point group
reduces $C$ to $n$ independent constants, from 21 (triclinic) to 3 (cubic). The question is how few
patterns determine all of them.

## Symmetry-adapted form

minelastic computes the form of $C$ numerically from the Cartesian point-group rotations of the given cell.
It group-averages the 21 symmetric unit matrices in Mandel notation, so the form holds in any orientation.
Each independent constant is labelled by its leading matrix entry (c11, c12, ...), and dependent entries
are written out, e.g. $C_{66} = (c_{11}-c_{12})/2$ for hexagonal crystals (`minelastic info`).

## Uncoupled patterns and minimal sets

Let the pattern $u$ have nonzero amplitudes on the Voigt components in its support $S$. There are
$2^6 - 1 = 63$ possible supports. The stress component $\sigma_i$ of the pattern is uncoupled if, for
generic amplitudes on $S$, it depends on a single independent constant. That constant then follows from
one stress component of one strained cell. For a cubic crystal, the single pattern $u = (1,0,0,1,0,0)$
(`xx + yz`) gives

$$
\sigma_{xx} = c_{11}\,\delta, \qquad \sigma_{yy} = \sigma_{zz} = c_{12}\,\delta, \qquad
\sigma_{yz} = c_{44}\,\delta ,
$$

so all three constants come from one strained cell.

The paper constructs the patterns from the sparsity of the symmetry-reduced stiffness matrix, represented
by an auxiliary binary matrix. A greedy procedure selects supports until every constant appears in an
uncoupled stress component. Since there are only 63 supports, the minimum number can be verified by
exhaustive search. minelastic accepts the greedy set for a given crystal only if it reaches the lower bound
derived [below](#why-fewer-cells-are-impossible), so it is provably minimal. Otherwise it searches all
combinations of supports. The amplitude signs are chosen such that no uncoupled relation cancels. The
resulting sets, for the standard settings:

| elastic type | example SG | constants | patterns | minimal set |
|---|---|---|---|---|
| triclinic | 2 | 21 | 6 | xx ; yy ; zz ; yz ; xz ; xy |
| monoclinic (unique axis b) | 12 | 13 | 4 | xx + yz ; yy + xy ; zz ; xz |
| orthorhombic | 70 | 9 | 3 | xx + yz + xz + xy ; yy ; zz |
| tetragonal I | 88 | 7 | 2 | xx + yz ; zz + xy |
| tetragonal II | 136 | 6 | 2 | xx + yz + xy ; zz |
| trigonal I | 148 | 7 | 2 | xx ; zz + yz |
| trigonal II | 154 | 6 | 2 | xx ; zz + yz |
| hexagonal | 194 | 5 | 2 | xx + yz ; zz |
| cubic | 225 | 3 | 1 | xx + yz |

## Extraction and consistency checks

* `uncoupled` (default): each constant is the average of all stress components that contain it alone. For
  cubic crystals, for example, c12 is the average from $\sigma_{yy}$ and $\sigma_{zz}$.
* `designated`: only the first such component, as in the tables of the paper.
* `lstsq`: least squares over all stress components of all patterns.

The stress components that contain several constants are not used for the extraction. They are compared
with the fitted tensor and reported as `checks` (summary line "consistency"). A large mismatch points to
nonlinearity, noise, or a reference structure that is not relaxed.

## Why fewer cells are impossible

The strains invariant under the point group, i.e. the free lattice parameters, form a subspace $V_0$ of
dimension 6, 4, 3, 2, 2, 2, 2, 2, 1 for the nine types above. Every rotation maps the $V_0$-component of a
strain onto itself. Suppose fewer strained cells were used than $\dim V_0$. Then there is an invariant strain
$w \neq 0$ orthogonal to the $V_0$-components of all applied strains. Adding $t\,w w^{\mathsf T}$ to $C$
(in Mandel notation) leaves every computed stress unchanged, yet it is an admissible elastic tensor of the
same symmetry. So $C$ would be undetermined. The minimal uncoupled sets reach this bound in every class. No
scheme, coupled or uncoupled, can use fewer strained cells.

## Internal-strain and piezoelectric tensors

Let a space-group operation $\{R|t\}$ map atom $n$ onto atom $g(n)$. Then the forces satisfy
$F_{g(n)}(R\varepsilon R^{\mathsf T}) = R\,F_n(\varepsilon)$. The forces of the few strained cells
therefore give the response of every atom to the whole symmetry orbit of the applied strains. That orbit
spans all six strain components, since otherwise the elastic tensor would be undetermined by the argument
above. The internal-strain tensor $\Lambda_n = \partial F_n/\partial\varepsilon$ of every atom follows
exactly, by least squares over the orbit (`collect --internal-strain`).

The polarization is a polar vector, $P(R\varepsilon R^{\mathsf T}) = R\,P(\varepsilon)$, so the
polarizations of the same cells determine the whole piezoelectric tensor $e_{iJ} = \partial
P_i/\partial\varepsilon_J$ in the same way. Berry-phase dipoles are defined only modulo a quantum.
minelastic writes the dipole of each cell in lattice coordinates, $p = \sum_a \varphi_a \mathbf{a}_a$, and aligns them
with the reference. This gives the *proper* piezoelectric response,
$e = \Omega_0^{-1} A_0\, \partial\varphi/\partial\varepsilon$ (D. Vanderbilt, J. Phys. Chem. Solids
**61**, 147 (2000)). Clamped-ion cells give the clamped-ion tensor $e^0$.

## Relaxed-ion response from clamped calculations

With atomic displacements $u$ and strain $\varepsilon$, the energy per cell is, to second order,

$$
E(u, \varepsilon) = \frac{\Omega}{2}\,\varepsilon^{\mathsf T} C^0 \varepsilon
- u^{\mathsf T}\Lambda\,\varepsilon + \frac12\, u^{\mathsf T}\Phi\, u ,
\qquad
P = e^0\varepsilon + \frac{1}{\Omega}\, Z^* u .
$$

Here $\Phi$ are the $\Gamma$-point force constants, $\Lambda$ the internal-strain tensor, and $Z^*$ the Born
effective charges. Relaxing the ions, $u = \Phi^{+}\Lambda\varepsilon$, gives

$$
C = C^0 - \frac{1}{\Omega}\,\Lambda^{\mathsf T}\Phi^{+}\Lambda , \qquad
e = e^0 + \frac{1}{\Omega}\, Z^*\Phi^{+}\Lambda , \qquad d = e\,C^{-1} ,
$$

where $\Phi^{+}$ is the inverse of $\Phi$ on the complement of the three uniform translations (X. Wu,
D. Vanderbilt and D. R. Hamann, Phys. Rev. B **72**, 035105 (2005)). With `--ionic`, minelastic displaces
only the symmetry-irreducible atoms. Each one moves along x, then y, then z, as long as the new direction
enlarges the span of its images under the site symmetry. $\Phi$ and $Z^*$ are completed over the symmetry
orbit of the displacements. $\Lambda$ is taken from the stresses of the displaced cells (the Maxwell
relation $\partial\sigma/\partial u = -\Lambda^{\mathsf T}/\Omega$, as in VASP). The forces in the
strained cells give it a second time, as a cross-check.

## Finite temperature

At finite temperature, $C_{ij}(T) = \partial\langle\sigma_i\rangle/\partial\varepsilon_j$, with the
stresses averaged over NVT trajectories in the strained cells of the thermally expanded lattice. The strain
sets are the same. See [Finite temperature](finite-temperature.md) for the protocol, the kinetic
term and the error estimates.
