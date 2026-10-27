# Piezoelectricity

The **piezoelectric stress tensor** is the change of the polarization $P = d/\Omega$ (dipole moment $d$ of
the cell per volume) with strain,

$$
e_{iJ} = \frac{\partial P_i}{\partial \varepsilon_J}
$$

in C/m$^2$, for Cartesian direction $i$ and Voigt strain component $J$. The **strain coefficients** follow
with the relaxed-ion compliance $S = C^{-1}$,

$$
d = e\, S
$$

in pC/N, and the **Born effective charges** are the dipole change per atomic displacement,

$$
Z^*_{n, \alpha\beta} = \frac{\partial d_\alpha}{\partial u_{n\beta}}
$$

in units of the elementary charge. Charge neutrality requires $\sum_n Z^*_n = 0$ (acoustic sum rule);
minelastic reports the sum as a convergence check.

**Proper response.** A Berry-phase polarization is defined modulo a quantum. minelastic writes the dipole of
each cell in lattice coordinates, $d = \sum_a \phi_a a_a$, brings the $\phi_a$ of all cells of a pattern to one
branch ($\phi \to \phi - q\, \mathrm{round}((\phi - \phi_{\mathrm{prev}})/q)$, in order of increasing strain),
and differentiates

$$
e_{iJ} = \frac{1}{\Omega_0} \sum_a a^0_{a,i}\, \frac{\partial \phi_a}{\partial \varepsilon_J},
$$

with the lattice vectors $a^0_a$ and the volume $\Omega_0$ of the reference cell (Vanderbilt, J. Phys. Chem.
Solids 61, 147 (2000)). Unlike the derivative of $P$ itself (the improper tensor), this proper tensor does
not depend on the branch of the polarization.

**Clamped ions.** The dipoles of the strained cells give the clamped-ion tensor $e^0$ with the same
finite-difference formula as the stresses. The polarization is a polar vector,
$P(R\varepsilon R^T) = R\, P(\varepsilon)$ for every point-group rotation $R$, so the responses of the
computed patterns determine the whole tensor by least squares over their symmetry orbits, as for the
[internal-strain tensor](internal-strain.md).

**Relaxed ions.** The displaced atoms give $Z^*$ from their dipoles. With the relaxation per unit strain
$u = \Phi^+ \Lambda$ ([Relaxed ions](relaxed-ions.md)), the polarization $P = e^0 \varepsilon + Z^* u/\Omega$
gives the total (relaxed-ion) tensor

$$
e = e^0 + \frac{1}{\Omega}\, Z^* \Phi^+ \Lambda .
$$

With $Z^*$ in units of $e$, $\Phi^+\Lambda$ in Å and $\Omega$ in Å$^3$, the second term is in $e$/Å$^2$
(1 $e$/Å$^2$ = 16.02 C/m$^2$).

The dipoles can come from a Berry-phase calculation (read from VASP with `LCALCPOL`), from any function of
the atoms (`calculate(..., dipole=f)`), or be passed to `ElasticSetup.compute(..., polarization=...)`.
Instead, $Z^*$ and $e^0$ can be taken from a DFPT run of the reference structure.
