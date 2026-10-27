# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Piezoelectric stress tensor e_iJ = dP_i/d eps_J from the dipoles of the strained cells.

The polarization is a polar vector: P(R eps R^T) = R P(eps) for every rotation R of the point group. The
responses of the K computed patterns therefore give the response to the whole orbit {R eps_k R^T}, which
spans the strain space when the strain set determines C; e follows by least squares over the orbit.

Berry-phase dipoles are defined modulo a quantum. With the dipole of a cell written in lattice coordinates,
d = sum_a phi_a a_a, the phases of the cells at +-delta are brought to the same branch (within half a
quantum), and the proper piezoelectric response is (Vanderbilt, J. Phys. Chem. Solids 61, 147 (2000))

    e_iJ = (1 / Omega_0) sum_a a0_{a,i} dphi_a / d eps_J,

with a0_a and Omega_0 the lattice vectors and volume of the reference cell. Clamped ions give the
clamped-ion tensor e0; ions relaxed in every strained cell give the total tensor.
"""
import numpy as np

from .voigt import W, strain_to_tensor, tensor_to_mandel

E_PER_A2_IN_C_PER_M2 = 16.021766208


def reduced_dipole(cell, dipole):
    """Lattice coordinates phi of a Cartesian dipole: d = A^T phi (rows of A = lattice vectors)."""
    return np.linalg.solve(np.asarray(cell, float).T, np.asarray(dipole, float))


def align_branch(phis, quantum=1.0, x=None):
    """Reduced dipoles brought to one Berry-phase branch, componentwise.

    Taken in the order of increasing x (the strain or displacement of each point; given order if None),
    each phi_k is moved by whole quanta q to within q/2 of its predecessor:
    phi_k -> phi_k - q round((phi_k - phi_prev) / q).
    """
    phis = np.array(phis, float)
    order = np.argsort(x, kind="stable") if x is not None else np.arange(len(phis))
    for a, b in zip(order[:-1], order[1:]):
        phis[b] -= quantum * np.round((phis[b] - phis[a]) / quantum)
    return phis


def orbit_complete(rotations, patterns, responses):
    """e (3, 6) in Voigt columns (per unit engineering strain) from the responses r_k = dP/d(step) of the
    patterns u_k.

    With m = Mandel vector of R eps(u_k) R^T and the response R r_k for all (R, k): e in Mandel form is
    Z Y^-1 with Y = sum m m^T, Z = sum (R r_k) m^T; the Voigt form divides column a by w_a. Raises ValueError
    if Y is singular to 1e-8 relative (the orbit does not span the strain space).
    """
    Y = np.zeros((6, 6))
    Z = np.zeros((3, 6))
    for R in rotations:
        for u, r in zip(np.asarray(patterns, float), np.asarray(responses, float)):
            m = tensor_to_mandel(R @ strain_to_tensor(u) @ R.T)
            Y += np.outer(m, m)
            Z += np.outer(R @ r, m)
    ev = np.linalg.eigvalsh(Y)
    if ev[0] <= 1e-8 * ev[-1]:
        raise ValueError("the symmetry orbit of the strain patterns does not span all six strain components")
    return (Z @ np.linalg.inv(Y)) / W[None, :]
