# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Relaxed-ion (ionic) contributions to the elastic and piezoelectric tensors from atomic displacements.

To second order in the atomic displacements u (3N) and the engineering strain e (6), per cell of volume
Omega (Wu, Vanderbilt and Hamann, Phys. Rev. B 72, 035105 (2005)):

    E(u, e) = Omega/2 e.C0.e - u.Lambda.e + 1/2 u.Phi.u,     P = e0 e + (1/Omega) Z* u,

with Phi the Gamma-point force constants, Lambda = dF/de the internal-strain tensor and Z* the Born charges
(F = -dE/du). At fixed strain the ions relax to u = Phi^+ Lambda e, which gives

    C = C0 - (1/Omega) Lambda^T Phi^+ Lambda,     e = e0 + (1/Omega) Z* Phi^+ Lambda,     d = e C^-1,

where Phi^+ is the inverse of Phi on the complement of the three uniform translations.

Only symmetry-irreducible displacements are computed. If {R|t} maps atom n onto g(n), displacing atom
g(a) by R d instead of atom a by d gives the forces F_g(m) -> R F_m, the dipole change R dp and the stress
R sigma R^T. Each atom's response to the whole orbit {R d} follows, and every response matrix is a least
squares fit over that orbit.
"""
import numpy as np

from .voigt import W, stress_to_tensor, tensor_to_mandel, tensor_to_stress

EV_PER_A3_IN_GPA = 160.21766208
E_PER_A2_IN_C_PER_M2 = 16.021766208


def irreducible_displacements(operations, natoms):
    """[(atom, unit vector)] whose symmetry orbits give the response of every atom to every direction.

    For the first atom n of each orbit of atoms, with site group S_n = {R : g(n) = n}, the Cartesian unit
    vectors x, y, z are taken in this order, each one kept if it raises rank{R d : R in S_n, d kept}, until
    the rank is 3 (as VASP chooses the displacements for IBRION = 6).
    """
    done, out = set(), []
    for n in range(natoms):
        if n in done:
            continue
        done |= {int(perm[n]) for _, perm in operations}
        site = [R for R, perm in operations if perm[n] == n]
        images, rank = np.zeros((0, 3)), 0
        for d in np.eye(3):
            trial = np.vstack([images] + [R @ d for R in site])
            r = np.linalg.matrix_rank(trial, tol=1e-6)
            if r > rank:
                out.append((n, d.copy()))
                images, rank = trial, r
            if rank == 3:
                break
    return out


def _orbit_fit(operations, natoms, displacements, transform):
    """M (natoms, n_out, 3) with response = M[j] @ u for a displacement u of atom j, by least squares.

    transform(R, perm, k) is the response (flattened, n_out) of the cell in which atom g(a_k) is displaced by
    R d_k. For each atom j, the pairs (u, y) = (R d_k, transform(R, perm, k)) over all operations with
    g(a_k) = j give M_j = Y^T U (U^T U)^-1. Raises ValueError if U^T U is singular to 1e-8 of its trace.
    """
    U = [[] for _ in range(natoms)]
    Y = [[] for _ in range(natoms)]
    for k, (a, d) in enumerate(displacements):
        for R, perm in operations:
            j = int(perm[a])
            U[j].append(R @ d)
            Y[j].append(transform(R, perm, k))
    out = []
    for j in range(natoms):
        Uj, Yj = np.array(U[j]), np.array(Y[j])
        M = Uj.T @ Uj
        if np.linalg.eigvalsh(M)[0] < 1e-8 * np.trace(M):
            raise ValueError(f"the displacements do not determine the response of atom {j + 1} in all directions")
        out.append(Yj.T @ Uj @ np.linalg.inv(M))
    return np.array(out)


def force_constants(operations, natoms, displacements, force_responses, symmetrize=True):
    """Gamma-point force constants Phi[m, a, n, b] = -dF_ma / du_nb (eV/A^2), shape (natoms, 3, natoms, 3).

    force_responses : (K, natoms, 3) dF/du of the K displacements (atom a_k along the unit vector d_k)
    The orbit images use F_g(m) = R F_m. With symmetrize, Phi is replaced by (Phi + Phi^T) / 2 as a
    3 natoms x 3 natoms matrix.
    """
    fr = np.asarray(force_responses, float)

    def transform(R, perm, k):
        """Forces when atom g(a_k) is displaced by R d_k: G[g(m)] = R F_m (flattened)."""
        G = np.empty_like(fr[k])
        G[perm] = fr[k] @ R.T                        # atom g(m) feels R F_m
        return G.ravel()

    M = _orbit_fit(operations, natoms, displacements, transform)      # (j, 3 natoms, 3): dF/du_j
    Phi = -M.reshape(natoms, natoms, 3, 3).transpose(1, 2, 0, 3)       # [m, a, j, b]
    if symmetrize:
        P2 = Phi.reshape(3 * natoms, 3 * natoms)
        Phi = (0.5 * (P2 + P2.T)).reshape(natoms, 3, natoms, 3)
    return Phi


def born_charges(operations, natoms, displacements, dipole_responses):
    """Born effective charges Z[n, a, b] = dp_a / du_nb (|e|), shape (natoms, 3, 3).

    dipole_responses : (K, 3) dp/du of the cell dipole (e Angstrom per Angstrom, physical sign) for the K
    displacements; the orbit images use dp -> R dp.
    """
    dr = np.asarray(dipole_responses, float)
    return _orbit_fit(operations, natoms, displacements, lambda R, perm, k: R @ dr[k])


def internal_strain_from_stress(operations, natoms, displacements, stress_responses, volume):
    """Internal-strain tensor from the stresses of the displaced cells: Lambda[n, b, J] = -Omega dsigma_J/du_nb.

    This is the mixed derivative -d2E/(du de) taken in the other order (sigma_J = (1/Omega) dE/de_J), so it
    equals the Lambda from the forces in the strained cells. Units: eV/Angstrom per unit engineering strain
    (Omega in Angstrom^3, the stress converted from GPa). stress_responses : (K, 6) dsigma/du in GPa/Angstrom,
    Voigt, tension positive; the orbit images use sigma -> R sigma R^T.
    """
    sr = np.asarray(stress_responses, float)
    M = _orbit_fit(operations, natoms, displacements,
                   lambda R, perm, k: tensor_to_mandel(R @ stress_to_tensor(sr[k]) @ R.T))   # (n, 6, 3), Mandel
    # Mandel -> Voigt stress components, then Lambda = -Omega dsigma/du in eV/A
    return -volume / EV_PER_A3_IN_GPA * (M / W[None, :, None]).transpose(0, 2, 1)


def acoustic_projector(natoms):
    """(P, T): T = [I3; I3; ...] / sqrt(natoms) (3 natoms x 3) spans the uniform translations, P = I - T T^T."""
    T = np.kron(np.ones((natoms, 1)), np.eye(3)) / np.sqrt(natoms)
    return np.eye(3 * natoms) - T @ T.T, T


def phi_pinv(Phi, tol=1e-6):
    """Phi^+ = sum_k v_k v_k^T / w_k over the eigenpairs (w_k, v_k) of P Phi P without the three translations.

    The three eigenvectors with the largest overlap with the translations are excluded (for one atom per cell
    nothing remains and Phi^+ = 0). Raises ValueError if
    a remaining eigenvalue is below tol times the largest |w| (an unstable or soft optical mode at Gamma: the
    relaxed-ion response is then undefined).
    """
    n = Phi.shape[0]
    if n == 1:
        return np.zeros((3, 3))                       # one atom per cell: no optical modes, no relaxation
    P2 = Phi.reshape(3 * n, 3 * n)
    Pr, T = acoustic_projector(n)
    A = Pr @ P2 @ Pr
    w, V = np.linalg.eigh(A)
    keep = np.ones(len(w), bool)
    keep[np.argsort(np.abs(V.T @ T).sum(axis=1))[-3:]] = False       # the three translations
    wk = w[keep]
    if wk.min() <= tol * np.abs(wk).max():
        raise ValueError(f"force-constant matrix has an unstable or soft optical mode at Gamma "
                         f"(eigenvalue {wk.min():.3g} eV/A^2); the relaxed-ion response is undefined")
    Vk = V[:, keep]
    return Vk @ np.diag(1.0 / wk) @ Vk.T


def ionic_contributions(Lambda, Phi, volume, born=None):
    """(C_ion, e_ion, A) with A = Phi^+ P Lambda, C_ion = -(P Lambda)^T A / Omega and e_ion = Z* A / Omega.

    A (natoms, 3, 6): relaxation of the ions per unit engineering strain (Angstrom), u = A e. C_ion (6, 6) in
    GPa, symmetrized; e_ion (3, 6) in C/m^2, None without born. Lambda (natoms, 3, 6) in eV/A, Phi in eV/A^2,
    born (natoms, 3, 3) in |e|, volume Omega in A^3. P removes the net force (projector above).
    """
    n = Phi.shape[0]
    L = np.asarray(Lambda, float).reshape(3 * n, 6)
    Pr, _ = acoustic_projector(n)
    L = Pr @ L                                        # WHY: a net force under strain only shifts the crystal
    A = phi_pinv(Phi) @ L                             # u per unit strain (A)
    C_ion = -(L.T @ A) / volume * EV_PER_A3_IN_GPA
    C_ion = 0.5 * (C_ion + C_ion.T)
    e_ion = None
    if born is not None:
        Z = np.asarray(born, float).transpose(1, 0, 2).reshape(3, 3 * n)   # Z[alpha, (n, j)]
        e_ion = Z @ A / volume * E_PER_A2_IN_C_PER_M2
    return C_ion, e_ion, A.reshape(n, 3, 6)


def acoustic_sum(born):
    """sum_n Z*_n (3, 3); zero for an exact calculation (acoustic sum rule, charge neutrality)."""
    return np.asarray(born, float).sum(axis=0)


def strain_coefficients(e, C):
    """d = e C^-1 (3, 6) in pC/N from e (C/m^2) and C (GPa): the strain e_J = d_iJ E_i per unit field
    (d_i4, d_i5, d_i6 refer to engineering shear strains)."""
    # WHY: C/m^2 times 1/GPa is 1e-9 C/N = 1000 pC/N
    return 1000.0 * np.asarray(e, float) @ np.linalg.inv(np.asarray(C, float))


def rotate_born(born, R):
    """Z'_n = R Z_n R^T for every atom n: the Born charges in the frame x' = R x."""
    return np.einsum("ia,nab,jb->nij", R, np.asarray(born, float), R)


def rotate_piezo(e, R):
    """e'_ijk = R_ia R_jb R_kc e_abc: the piezoelectric tensor (3, 6, Voigt columns) in the frame x' = R x."""
    # WHY: P_i = e_ijk eps_jk = sum_J e_iJ e_J with engineering shear e_J, so e_iJ = e_ijk (no factor 2):
    #      the rows transform like stress vectors
    full = np.array([stress_to_tensor(row) for row in np.asarray(e, float)])       # e[i, j, k]
    rot = np.einsum("ia,jb,kc,abc->ijk", R, R, R, full)
    return np.array([tensor_to_stress(t) for t in rot])
