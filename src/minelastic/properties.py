# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Moduli derived from the stiffness matrix: Voigt, Reuss and Hill averages; Born stability.

With C in GPa (Voigt notation) and S = C^-1:
    K_V = [C11 + C22 + C33 + 2 (C12 + C23 + C13)] / 9
    G_V = [C11 + C22 + C33 - (C12 + C23 + C13) + 3 (C44 + C55 + C66)] / 15
    1/K_R = S11 + S22 + S33 + 2 (S12 + S23 + S13)
    15/G_R = 4 (S11 + S22 + S33) - 4 (S12 + S23 + S13) + 3 (S44 + S55 + S66)
    K_H = (K_V + K_R)/2,  G_H = (G_V + G_R)/2
    E_H = 9 K_H G_H / (3 K_H + G_H),  nu_H = (3 K_H - 2 G_H) / (2 (3 K_H + G_H))
    A_U = 5 G_V/G_R + K_V/K_R - 6   (universal anisotropy index; 0 for an isotropic material)
"""
import numpy as np


def vrh(C):
    """Dict K_V, K_R, K_H, G_V, G_R, G_H, E_H, nu_H, A_U by the formulas above (GPa; nu_H, A_U dimensionless)."""
    C = np.asarray(C, float)
    S = np.linalg.inv(C)
    KV = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * (C[0, 1] + C[1, 2] + C[0, 2])) / 9
    GV = (C[0, 0] + C[1, 1] + C[2, 2] - (C[0, 1] + C[1, 2] + C[0, 2]) + 3 * (C[3, 3] + C[4, 4] + C[5, 5])) / 15
    KR = 1 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * (S[0, 1] + S[1, 2] + S[0, 2]))
    GR = 15 / (4 * (S[0, 0] + S[1, 1] + S[2, 2]) - 4 * (S[0, 1] + S[1, 2] + S[0, 2])
               + 3 * (S[3, 3] + S[4, 4] + S[5, 5]))
    KH, GH = (KV + KR) / 2, (GV + GR) / 2
    E = 9 * KH * GH / (3 * KH + GH)
    nu = (3 * KH - 2 * GH) / (2 * (3 * KH + GH))
    AU = 5 * GV / GR + KV / KR - 6
    return dict(K_V=KV, K_R=KR, K_H=KH, G_V=GV, G_R=GR, G_H=GH, E_H=E, nu_H=nu, A_U=AU)


def born_stable(C):
    """(stable, eigenvalues): an unstressed crystal is stable if C is positive definite (Born criterion).

    Tested with the eigenvalues of the symmetric part of the Voigt matrix; positive definiteness of C and of
    the Mandel matrix W C W are equivalent.
    """
    ev = np.linalg.eigvalsh(0.5 * (np.asarray(C) + np.asarray(C).T))
    return bool(ev.min() > 0), ev
