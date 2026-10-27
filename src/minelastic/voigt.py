# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Voigt and Mandel notation for strain, stress and stiffness, and their rotation.

Conventions of the whole package:
* Voigt index a = 0..5 stands for the pair (i_a, j_a) = xx, yy, zz, yz, xz, xy.
* Strain vectors e use engineering shear: e_a = eps_ii for a < 3, e_a = 2 eps_ij for a >= 3.
  Stress vectors s do not: s_a = sigma_ij for all a.
* sigma is the Cauchy stress in GPa, positive in tension. The stiffness C (GPa) maps e to s: s = C e.
* Mandel form uses the weights w = (1, 1, 1, sqrt2, sqrt2, sqrt2): m_a(t) = w_a t_{i_a j_a} for a symmetric
  tensor t, so that eps_M = e / w and sigma_M = w s. In this form a rotation t -> R t R^T acts as an
  orthogonal 6x6 matrix Q(R), and the stiffness becomes M = W C W with W = diag(w).
"""
import numpy as np

PAIRS = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))
LABELS = ("xx", "yy", "zz", "yz", "xz", "xy")
W = np.array([1.0, 1.0, 1.0, np.sqrt(2.0), np.sqrt(2.0), np.sqrt(2.0)])


def strain_to_tensor(e):
    """Symmetric 3x3 strain eps from the engineering Voigt vector e: eps_ii = e_a (a < 3),
    eps_ij = eps_ji = e_a / 2 (a >= 3)."""
    e = np.asarray(e, float)
    t = np.empty((3, 3))
    for k, (i, j) in enumerate(PAIRS):
        t[i, j] = t[j, i] = e[k] if k < 3 else 0.5 * e[k]
    return t


def tensor_to_strain(t):
    """Engineering Voigt vector of a 3x3 strain: e_a = t_ii (a < 3), e_a = t_ij + t_ji (a >= 3)."""
    t = np.asarray(t, float)
    return np.array([t[i, j] if k < 3 else t[i, j] + t[j, i] for k, (i, j) in enumerate(PAIRS)])


def stress_to_tensor(s):
    """Symmetric 3x3 stress from the Voigt vector s: sigma_ij = sigma_ji = s_a."""
    s = np.asarray(s, float)
    t = np.empty((3, 3))
    for k, (i, j) in enumerate(PAIRS):
        t[i, j] = t[j, i] = s[k]
    return t


def tensor_to_stress(t):
    """Voigt vector of a 3x3 stress: s_a = (t_ij + t_ji) / 2."""
    t = np.asarray(t, float)
    return np.array([0.5 * (t[i, j] + t[j, i]) for i, j in PAIRS])


def mandel_rotation(R):
    """Orthogonal 6x6 matrix Q with m(R t R^T) = Q m(t) for every symmetric t (m = Mandel vector).

    Column b is m(R B_b R^T), where B_b is the symmetric tensor with m(B_b) = unit vector b, i.e.
    (B_b)_kl = (B_b)_lk = 1 / w_b for (k, l) = pair b.
    """
    R = np.asarray(R, float)
    Q = np.empty((6, 6))
    for b, (k, l) in enumerate(PAIRS):
        t = np.zeros((3, 3))
        t[k, l] = t[l, k] = 1.0 / W[b]          # Mandel unit vector b as a tensor
        r = R @ t @ R.T
        Q[:, b] = [r[i, j] * W[a] for a, (i, j) in enumerate(PAIRS)]
    return Q


def stiffness_to_mandel(C):
    """Mandel stiffness M_ab = w_a C_ab w_b (so that sigma_M = M eps_M)."""
    return W[:, None] * np.asarray(C, float) * W[None, :]


def stiffness_from_mandel(M):
    """Voigt stiffness C_ab = M_ab / (w_a w_b)."""
    return np.asarray(M, float) / (W[:, None] * W[None, :])


def rotate_stiffness(C, R):
    """Voigt stiffness of the crystal rotated rigidly by R (x' = R x):
    C'_ijkl = R_ip R_jq R_kr R_ls C_pqrs, evaluated as C' = W^-1 Q M Q^T W^-1 with M = W C W, Q = Q(R)."""
    Q = mandel_rotation(R)
    return stiffness_from_mandel(Q @ stiffness_to_mandel(C) @ Q.T)


def deformation_gradient(e):
    """Deformation gradient F = I + eps(e) of a pure strain (symmetric F, no rigid rotation)."""
    return np.eye(3) + strain_to_tensor(e)


def tensor_to_mandel(t):
    """Mandel vector of the symmetric part of t: m_a = w_a (t_ij + t_ji) / 2."""
    t = np.asarray(t, float)
    return np.array([0.5 * (t[i, j] + t[j, i]) * W[a] for a, (i, j) in enumerate(PAIRS)])
