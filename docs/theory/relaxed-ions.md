# Relaxed ions

A strain moves the atoms that are not fixed by symmetry. To second order in the atomic displacements $u$
(3N components) and the strain $\varepsilon$, the energy of a cell of volume $\Omega$ is

$$
E(u, \varepsilon) = E_0 + \frac{\Omega}{2}\, \varepsilon \cdot C^0 \varepsilon
  - u \cdot \Lambda \varepsilon + \frac{1}{2}\, u \cdot \Phi u ,
$$

with the clamped-ion constants $C^0$, the internal-strain tensor $\Lambda = \partial F/\partial \varepsilon$
and the $\Gamma$-point force constants $\Phi = -\partial F/\partial u$ ([Internal strain](internal-strain.md)).
The forces $F = -\partial E/\partial u = \Lambda\varepsilon - \Phi u$ vanish for

$$
u = \Phi^+ \Lambda\, \varepsilon ,
$$

where $\Phi^+$ is the inverse of $\Phi$ on the complement of the three uniform translations. The stress
$\sigma = \Omega^{-1}\, \partial E/\partial \varepsilon$ at these positions gives the relaxed-ion constants

$$
C = C^0 - \frac{1}{\Omega}\, \Lambda^T \Phi^+ \Lambda .
$$

With $\Lambda$ in eV/Å, $\Phi$ in eV/Å$^2$ and $\Omega$ in Å$^3$, the second term is in eV/Å$^3$
(1 eV/Å$^3$ = 160.2 GPa). minelastic computes $C$ in one of two ways.

1. **Relaxing the ions** (`relax_ions=True`, `--relax-ions`). The atoms of every strained cell are relaxed
    at fixed cell, and the stresses are differentiated as for $C^0$. The relaxation at finite strain also
    contains terms beyond linear order; the four-point formula removes them to order $\delta^4$.
2. **Clamped single points** (`ionic=True`, `--ionic`), the route of VASP's `IBRION = 6`. The strained cells
    give $C^0$. The symmetry-irreducible atoms, displaced by $\pm h$ and $\pm 2h$ (default $h$ = 0.01 Å), give
    $\Phi$ from their forces and $\Lambda$ from their stresses. $C$ follows from the formula above, as the exact
    linear response. If $\Phi$ has an unstable or soft optical mode, the relaxed-ion response is undefined and
    minelastic stops with an error.

The second route needs no relaxations. Relaxations are costly, and they must be converged tightly: the forces
that a strain of 0.01 induces are only about 0.05 eV/Å.
