# Internal strain and force constants

The **internal-strain tensor** is the force on each atom per unit strain at fixed fractional coordinates,

$$
\Lambda_{n\alpha, J} = \left.\frac{\partial F_{n\alpha}}{\partial \varepsilon_J}\right|_{\varepsilon = 0},
$$

for atom $n$, Cartesian direction $\alpha$ and Voigt strain component $J$ (engineering shear), in eV/Å per
unit strain. The net force vanishes, so $\sum_n \Lambda_n = 0$. minelastic obtains it in two ways.

1. **From the strained cells**, without further calculations. The force responses of each pattern,
    $f_k = \delta^{-1} \sum_m g_m F(m\delta u_k)$, use the same finite-difference formula as the stresses. If a space-group
    operation $\{R|t\}$ maps atom $n$ onto atom $g(n)$, the forces obey

    $$
    F_{g(n)}(R \varepsilon R^T) = R\, F_n(\varepsilon) .
    $$

    The computed patterns therefore give the response of every atom to the whole orbit of strains
    $\{R^T \varepsilon_k R\}$. This orbit spans all six strain components whenever the patterns determine $C$,
    so the minimal set gives the complete tensor; $\Lambda_n$ follows by least squares over the orbit.

2. **From displaced atoms** (`ionic=True`). The mixed second derivative of the energy taken in the other
    order gives

    $$
    \Lambda_{n\beta, J} = -\Omega\, \frac{\partial \sigma_J}{\partial u_{n\beta}} ,
    $$

    from the stresses of the cells in which atom $n$ is displaced along $\beta$ (the route of VASP's
    `IBRION = 6`).

The displaced atoms also give the **$\Gamma$-point force constants**

$$
\Phi_{m\alpha, n\beta} = -\frac{\partial F_{m\alpha}}{\partial u_{n\beta}}
$$

in eV/Å$^2$. Only the symmetry-irreducible displacements are computed: for the first atom of each orbit,
the Cartesian directions x, y, z are kept as long as they add a new direction to the images under its site
group, as VASP chooses them. The other atoms and directions follow from $F_{g(m)} = R F_m$, and $\Phi$ is
symmetrized as a $3N \times 3N$ matrix. The **relaxation of the atoms per unit strain** is

$$
u = \Phi^+ \Lambda
$$

in Å, with $\Phi^+$ the inverse of $\Phi$ on the complement of the uniform translations. It enters the
relaxed-ion constants ([Relaxed ions](relaxed-ions.md)) and the piezoelectric tensor
([Born charges and piezoelectricity](piezoelectricity.md)).
