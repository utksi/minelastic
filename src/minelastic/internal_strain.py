# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Internal-strain tensor Lambda_n = dF_n/d eps of every atom from the forces in the strained cells.

For a space-group operation {R|t} that maps atom n onto atom g(n), the forces in a strained cell satisfy
F_g(n)(R eps R^T) = R F_n(eps). Hence the force responses f_k = dF/d(step) of the K computed patterns give
F_n(R^T eps_k R) = R^T F_g(n)(eps_k) for every operation: the response of every atom to the whole orbit of
strains {R^T eps_k R}. When the strain set determines C, this orbit spans all six strain components
(otherwise C would be undetermined), and Lambda_n follows by linear least squares over the orbit.
"""
import numpy as np

from .voigt import W, strain_to_tensor, tensor_to_mandel


def internal_strain(operations, patterns, force_responses, rcond=1e-8):
    """Lambda (N, 3, 6), eV/Angstrom per unit engineering strain: F_n = Lambda[n] @ e to first order.

    operations      : [(R, perm)] space-group operations of the reference cell, perm[n] = g(n)
    patterns        : (K, 6) strain patterns u_k per unit step (engineering Voigt)
    force_responses : (K, N, 3) f_k = dF/d(step) of each pattern

    With m = Mandel vector of R^T eps(u_k) R and f = R^T f_k[g(n)] for all (R, k), Lambda_n in Mandel form
    is Z_n Y^-1 with Y = sum m m^T and Z_n = sum f m^T (normal equations); the Voigt form divides column a
    by w_a. Raises ValueError if the smallest eigenvalue of Y is below rcond times the largest (the orbit
    does not span the strain space).
    """
    patterns = np.asarray(patterns, float)
    fr = np.asarray(force_responses, float)
    K, N, _ = fr.shape
    Y = np.zeros((6, 6))
    Z = np.zeros((N, 3, 6))
    for R, perm in operations:
        for k in range(K):
            m = tensor_to_mandel(R.T @ strain_to_tensor(patterns[k]) @ R)
            Y += np.outer(m, m)
            f = fr[k][perm] @ R                      # row n: R^T F_{g(n)}
            Z += f[:, :, None] * m[None, None, :]
    ev = np.linalg.eigvalsh(Y)
    if ev[0] <= rcond * ev[-1]:
        raise ValueError("the symmetry orbit of the strain patterns does not span all six strain components")
    return (Z @ np.linalg.inv(Y)) / W[None, None, :]
