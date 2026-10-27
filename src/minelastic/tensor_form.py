# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Symmetry-adapted form of the stiffness matrix for a set of Cartesian point-group rotations.

The admissible stiffness matrices are the symmetric 6x6 matrices invariant under M -> Q(R) M Q(R)^T (Mandel
form, Q from voigt.mandel_rotation) for every rotation R of the point group G. The group average
P(M) = (1/|G|) sum_R Q(R) M Q(R)^T projects onto them. Projecting the 21 Voigt unit matrices spans the
admissible subspace. The reduced row echelon form of a basis of that span, with the entries ordered
c11, c12, ..., c16, c22, ..., c66, gives one basis matrix B_p per independent constant, scaled so that the
constant equals its pivot entry (the first entry in that order that is independent of the earlier ones):

    C = sum_p c_p B_p,   c_p = C[pivot_p].

This reproduces the usual labels (hexagonal: c11, c12, c13, c33, c44; C66 = (c11 - c12)/2 is dependent).
Because the form is computed from the rotations, it holds for any orientation of the cell.
"""
from fractions import Fraction

import numpy as np

from .voigt import W, mandel_rotation, stiffness_from_mandel

UPPER = [(i, j) for i in range(6) for j in range(i, 6)]          # c11 c12 ... c16 c22 ... c66
_TOL = 1e-8


def _rref(A, tol=_TOL):
    """Reduced row echelon form of A by Gauss-Jordan elimination with partial pivoting.

    A column is skipped when its largest remaining entry is below tol * max(|A|, 1). Returns the nonzero
    rows and the list of pivot columns.
    """
    A = np.array(A, float)
    rows, cols = A.shape
    pivots, r = [], 0
    scale = max(np.abs(A).max(), 1.0)
    for c in range(cols):
        if r == rows:
            break
        p = r + int(np.argmax(np.abs(A[r:, c])))
        if abs(A[p, c]) <= tol * scale:
            continue
        A[[r, p]] = A[[p, r]]
        A[r] /= A[r, c]
        for i in range(rows):
            if i != r:
                A[i] -= A[i, c] * A[r]
        pivots.append(c)
        r += 1
    return A[:r], pivots


def _snap(x):
    """Entries below 1e-10 in magnitude set to 0, entries within 1e-9 of a fraction with denominator <= 12
    replaced by that fraction (the symmetry relations between stiffness entries are rational)."""
    x = np.where(np.abs(x) < 1e-10, 0.0, x)
    for idx, v in np.ndenumerate(x):
        f = Fraction(v).limit_denominator(12)
        if abs(float(f) - v) < 1e-9:
            x[idx] = float(f)
    return x


def invariant_strain_dimension(rotations):
    """dim V0 = (1/|G|) sum_R tr Q(R): the dimension of the space of strains with R eps R^T = eps for all R.

    This is the multiplicity of the identity representation in the strain representation, i.e. the number of
    free lattice parameters. It is the lower bound on the number of strain patterns from which C can be
    determined by stress-strain differences (docs/theory/index.md).
    """
    Q = np.mean([mandel_rotation(R) for R in rotations], axis=0)
    return int(round(np.trace(Q)))


class TensorForm:
    """Basis B_p (m, 6, 6) of the admissible stiffness matrices: C = sum_p c_p B_p.

    Attributes: basis (m, 6, 6), Voigt; labels ['c11', ...]; pivots [(i, j)] (0-based) with c_p = C[pivots[p]];
    coef (m, 21), the rows of the echelon form over the upper triangle; rotations.
    """

    def __init__(self, rotations, tol=1e-6):
        """Projects the 21 unit matrices E_ij = E_ji = 1 (Voigt, upper triangle) with P, takes the column space
        of the resulting 21 x 21 matrix (SVD, singular values above tol * the largest), and brings a basis of
        it into reduced row echelon form."""
        Qs = [mandel_rotation(R) for R in rotations]
        cols = []
        for i, j in UPPER:
            E = np.zeros((6, 6))
            E[i, j] = E[j, i] = W[i] * W[j]                         # Voigt unit entry, in Mandel form
            P = np.mean([Q @ E @ Q.T for Q in Qs], axis=0)           # group average (Mandel)
            V = stiffness_from_mandel(P)
            cols.append([V[a, b] for a, b in UPPER])
        M = np.array(cols).T
        U, s, _ = np.linalg.svd(M)
        rank = int(np.sum(s > tol * s[0]))
        R, piv = _rref(U[:, :rank].T)
        self.coef = _snap(R)                                         # (m, 21): entries per unit constant
        self.pivots = [UPPER[c] for c in piv]
        self.labels = [f"c{i + 1}{j + 1}" for i, j in self.pivots]
        self.rotations = [np.asarray(R_) for R_ in rotations]
        B = np.zeros((len(piv), 6, 6))
        for p in range(len(piv)):
            for c, (i, j) in enumerate(UPPER):
                B[p, i, j] = B[p, j, i] = self.coef[p, c]
        self.basis = B                                               # (m, 6, 6) Voigt

    @property
    def n(self):
        """Number m of independent constants."""
        return len(self.labels)

    @property
    def nonzero(self):
        """6x6 boolean mask of the entries that are not identically zero for this symmetry."""
        return np.any(np.abs(self.basis) > _TOL, axis=0)

    def entry(self, i, j):
        """Coefficients a_p of entry C_ij (0-based): C_ij = sum_p a_p c_p."""
        return self.basis[:, i, j]

    def matrix(self, c):
        """C = sum_p c_p B_p."""
        return np.tensordot(np.asarray(c, float), self.basis, axes=1)

    def project(self, C):
        """Constants c minimising || W (C - sum_p c_p B_p) W ||_F (least squares): the orthogonal projection of
        C onto the admissible matrices in the Mandel metric."""
        WW = W[:, None] * W[None, :]
        A = np.array([(Bp * WW).ravel() for Bp in self.basis]).T
        c, *_ = np.linalg.lstsq(A, (np.asarray(C, float) * WW).ravel(), rcond=None)
        return c

    def expression(self, i, j):
        """Entry C_ij as a formula in the constants, e.g. '(c11 - c12)/2'. Coefficients that are not fractions
        with denominator <= 12 are printed as decimals."""
        a = self.entry(i, j)
        if not np.any(a):
            return "0"
        fr = [Fraction(float(v)).limit_denominator(12) for v in a]
        if any(abs(float(f) - v) > 1e-9 for f, v in zip(fr, a)):       # not a simple rational relation
            return " ".join(f"{v:+.4g}*{lab}" for v, lab in zip(a, self.labels) if v != 0)
        terms = [(f, self.labels[p]) for p, f in enumerate(fr) if f != 0]
        if len(terms) == 1:
            f, lab = terms[0]
            return {1: lab, -1: "-" + lab}.get(f, f"{f}*{lab}")
        den = np.lcm.reduce([f.denominator for f, _ in terms])
        s = ""
        for f, lab in terms:
            k = int(f * den)
            s += ("+ " if k > 0 and s else "- " if k < 0 and s else "-" if k < 0 else "") + \
                 (f"{abs(k)}{lab} " if abs(k) != 1 else f"{lab} ")
        s = s.strip()
        return f"({s})/{den}" if den != 1 else s

    def describe(self):
        """The 6x6 matrix of expression(i, j), as text."""
        width = max(len(self.expression(i, j)) for i in range(6) for j in range(6))
        return "\n".join("  ".join(self.expression(i, j).rjust(width) for j in range(6)) for i in range(6))
