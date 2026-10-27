# Elastic constants

Stress and strain are 6-vectors in Voigt order xx, yy, zz, yz, xz, xy, with engineering shear strains:

$$
\sigma = (\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{yz}, \sigma_{xz}, \sigma_{xy}), \qquad
\varepsilon = (\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, 2\varepsilon_{yz}, 2\varepsilon_{xz}, 2\varepsilon_{xy}).
$$

$\sigma$ is the Cauchy stress in GPa, positive in tension. The elastic constants are its derivatives at the
stress-free reference structure,

$$
C_{ij} = \left.\frac{\partial \sigma_i}{\partial \varepsilon_j}\right|_{\varepsilon = 0}, \qquad i, j = 1, \dots, 6 .
$$

With the fractional coordinates of the atoms held fixed they are the clamped-ion constants $C^0$;
[Relaxed ions](relaxed-ions.md) adds the relaxation of the atoms.

**Strained cells.** A strain pattern is a 6-vector $u_k$. Pattern $k$ at the amplitude $m\delta$ deforms
the lattice vectors $a_i$ of the reference cell as

$$
a_i' = F a_i, \qquad F = I + \varepsilon(m\delta u_k),
$$

where $\varepsilon(\cdot)$ is the symmetric strain tensor of a Voigt vector (no rigid rotation). The
fractional coordinates of the atoms are kept.

**Finite differences.** The response of pattern $k$ is

$$
r_k = \frac{1}{\delta} \sum_m g_m\, \sigma(m\delta u_k) = C u_k + O(\delta^4),
$$

with $m = (-2, -1, 1, 2)$ and $g = (1, -8, 8, -1)/12$ for the default four-point formula, or $m = (-1, 1)$,
$g = (-1/2, 1/2)$ and an error of order $\delta^2$ for the two-point (central) formula. The default step is
$\delta = 0.01$.

**Symmetry.** The point group of the crystal (with magnetic order, if any) restricts $C$ to

$$
C = \sum_{p=1}^{n} c_p B_p ,
$$

with $n$ independent constants, from 21 for triclinic to 3 for cubic crystals. Each constant $c_p$ is
named after the leading entry of its basis matrix $B_p$ (c11, c12, ...), and the basis is computed from the
point-group rotations of the given cell ([Method](index.md)). Each response gives six linear relations

$$
r_{ki} = \sum_p a_{kip}\, c_p, \qquad a_{kip} = (B_p u_k)_i .
$$

**Extraction.** In a minimal uncoupled set every constant appears alone in at least one relation, so
$c_p = r_{ki}/a_{kip}$. The default method `uncoupled` averages $c_p$ over all such relations, `designated`
takes the first one, and `lstsq` fits all relations by weighted least squares. The relations that contain
more than one constant are not needed; minelastic reports the largest residual
$|r_{ki} - \sum_p a_{kip} c_p|$ as a consistency check. For hexagonal crystals the two patterns
$u_1$ = xx + yz and $u_2$ = zz give

$$
r_1 = (c_{11},\ c_{12},\ c_{13},\ c_{44},\ 0,\ 0), \qquad r_2 = (c_{13},\ c_{13},\ c_{33},\ 0,\ 0,\ 0),
$$

so c13 is obtained twice.

**Moduli and stability.** With $S = C^{-1}$, the Voigt and Reuss bounds of the bulk and shear moduli are

$$
K_V = \frac{C_{11} + C_{22} + C_{33} + 2(C_{12} + C_{23} + C_{13})}{9}, \qquad
G_V = \frac{C_{11} + C_{22} + C_{33} - (C_{12} + C_{23} + C_{13}) + 3(C_{44} + C_{55} + C_{66})}{15},
$$

$$
\frac{1}{K_R} = S_{11} + S_{22} + S_{33} + 2(S_{12} + S_{23} + S_{13}), \qquad
\frac{15}{G_R} = 4(S_{11} + S_{22} + S_{33}) - 4(S_{12} + S_{23} + S_{13}) + 3(S_{44} + S_{55} + S_{66}),
$$

with the Hill averages $K = (K_V + K_R)/2$ and $G = (G_V + G_R)/2$, Young's modulus $E = 9KG/(3K + G)$,
Poisson's ratio $\nu = (3K - 2G)/(2(3K + G))$ and the universal anisotropy index
$A_U = 5 G_V/G_R + K_V/K_R - 6$. The crystal is Born stable if all eigenvalues of $C$ are positive.

**Stress measure.** VASP's `IBRION = 6` differentiates $J\sigma$ with $J = \det F$;
`stress_measure="kirchhoff"` (`--stress-measure kirchhoff`) does the same. For a stress-free reference both
give the same constants.
