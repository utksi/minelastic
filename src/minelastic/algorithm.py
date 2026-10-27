# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Minimal uncoupled strain sets (Simak and Singh): the paper's construction and an exhaustive search.

A strain pattern is a vector u with support S (the Voigt components j with u_j != 0). Per unit step it
produces the stresses s_i = sum_j C_ij u_j = sum_p c_p sum_{j in S} (B_p)_ij u_j (C = sum_p c_p B_p, see
tensor_form). For generic amplitudes, s_i depends on constant p exactly when (B_p)_ij != 0 for some j in S.
s_i is uncoupled if it depends on one constant only; that constant then follows from s_i alone. A set of
patterns is complete if every constant is uncoupled in at least one of its stress components.

paper_greedy is the iterative construction of the paper (Sec. III); exact_minimum searches all
combinations of the 63 supports for the smallest complete set.
"""
from itertools import combinations, product

import numpy as np

_TOL = 1e-8


def single_constant_rows(form, support):
    """[(i, p)]: the stress components i that, for generic amplitudes on the support S, depend on exactly
    one constant p, i.e. (B_q)_ij != 0 for some j in S only for q = p."""
    out = []
    for i in range(6):
        involved = [p for p in range(form.n) if any(abs(form.basis[p, i, j]) > _TOL for j in support)]
        if len(involved) == 1:
            out.append((i, involved[0]))
    return out


def extractable(form, support):
    """The set of constants p that are uncoupled in some stress component of a pattern on the support."""
    return {p for _, p in single_constant_rows(form, support)}


def reduced_auxiliary(form):
    """Reduced auxiliary matrix of the paper (Sec. III B).

    A_ij = 1 where C_ij is not identically zero, else 0. Row i is set to zero when the diagonal entry C_ii is
    not +-1 times a single constant (a dependent diagonal such as C66 = (c11 - c12)/2), or when its constant
    already appears on an earlier diagonal (an equivalent row, such as C22 = c11).
    """
    A = form.nonzero.astype(int)
    seen = set()
    for i in range(6):
        d = form.entry(i, i)
        nz = np.flatnonzero(np.abs(d) > _TOL)
        if len(nz) != 1 or abs(abs(d[nz[0]]) - 1) > _TOL or nz[0] in seen:
            A[i, :] = 0              # dependent diagonal (e.g. c66 = (c11-c12)/2) or equivalent to an earlier row
        else:
            seen.add(nz[0])
    return A


def paper_greedy(form):
    """Strain supports from the iterative construction of the paper (Sec. III C).

    Start from the reduced auxiliary matrix A and the candidate support b (b_j = 1 for the nonzero rows of A).
    While some component of v = A b exceeds 1 (row i would respond to more than one strained component),
    remove the highest such index j from b and zero row j. The remaining b is a pattern; its rows are removed
    from A and from the candidates, and the procedure repeats until A is zero. Stops early, without a
    complete set, when no component survives (possible in a non-standard frame); the caller then uses
    exact_minimum.
    """
    A = reduced_auxiliary(form)
    b_in = (A.sum(axis=1) > 0).astype(int)
    supports = []
    while A.any():
        At, b = A.copy(), b_in.copy()
        while True:
            v = At @ b
            bad = np.flatnonzero(v > 1)
            if bad.size == 0:
                break
            j = bad.max()
            b[j] = 0
            At[j, :] = 0
        if not b.any():
            break                    # no progress possible (non-standard frame); caller falls back
        supports.append(tuple(int(j) for j in np.flatnonzero(b)))
        A[b == 1, :] = 0
        b_in[b == 1] = 0
    return supports


def covers(form, supports):
    """True if the union of extractable(S) over the supports contains every constant."""
    got = set()
    for s in supports:
        got |= extractable(form, s)
    return got == set(range(form.n))


def exact_minimum(form, max_patterns=6):
    """A complete set of supports of minimum size, by exhaustive search.

    Each of the 63 nonempty supports S is assigned X(S) = extractable(S). A support is dropped if another
    support S2 has X(S) a subset of X(S2) and no more components (one of equal ones is kept). For
    k = 1, 2, ..., max_patterns all
    k-combinations are tested for complete coverage; among those of the smallest k, the one with the
    fewest strained components in total (least symmetry lowering, then lexicographic) is returned.
    """
    target = set(range(form.n))
    X = {}
    for mask in range(1, 64):
        s = tuple(j for j in range(6) if mask >> j & 1)
        e = frozenset(extractable(form, s))
        if e:
            X[s] = e
    # WHY: a support is never needed if another one extracts at least the same constants with no more strained
    #      components (dropping only such supports keeps both the count and the tie-break exact)
    order = sorted(X, key=lambda s: (len(s), s))
    cands = [s for n, s in enumerate(order)
             if not any(X[s] <= X[s2] and len(s2) <= len(s) and (X[s] < X[s2] or len(s2) < len(s) or m < n)
                        for m, s2 in enumerate(order) if s2 != s)]
    for k in range(1, max_patterns + 1):
        best = None
        for combo in combinations(cands, k):
            if set().union(*(X[s] for s in combo)) == target:
                key = (sum(len(s) for s in combo), combo)
                if best is None or key < best[0]:
                    best = (key, combo)
        if best is not None:
            return list(best[1])
    raise RuntimeError("no complete set of uncoupled strain patterns exists")


def choose_amplitudes(form, support):
    """Pattern u on the support, first component 1, with every single-constant coefficient (B_p)_i . u nonzero.

    For generic amplitudes each single-constant row (i, p) has a nonzero coefficient, but special amplitudes
    can cancel it (e.g. c16 (u_1 - u_2) for u_1 = u_2). Amplitudes +-1 are tried first (all positive first),
    then +-1 and +-2. RuntimeError if none keeps all relations.
    """
    rows = single_constant_rows(form, support)
    for values in ((1.0, -1.0), (1.0, -1.0, 2.0, -2.0)):
        for amps in product(values, repeat=len(support) - 1):
            u = np.zeros(6)
            u[list(support)] = (1.0,) + amps
            if all(abs(form.basis[p, i] @ u) > _TOL for i, p in rows):
                return u
    raise RuntimeError(f"no amplitudes on support {support} keep all single-constant relations")


def relations(form, pattern):
    """The six relations s_i = sum_p a_ip c_p of a pattern u (per unit step), a_ip = (B_p)_i . u.

    One dict per row i: row, coef = (a_ip)_p, and kind 'zero' (all a_ip = 0), 'single' (one nonzero a_ip;
    then also constant = p and factor = a_ip, so that s_i = factor * c_p) or 'coupled'.
    """
    u = np.asarray(pattern, float)
    out = []
    for i in range(6):
        a = form.basis[:, i, :] @ u
        nz = np.flatnonzero(np.abs(a) > _TOL)
        if nz.size == 0:
            out.append(dict(row=i, kind="zero", coef=a))
        elif nz.size == 1:
            out.append(dict(row=i, kind="single", coef=a, constant=int(nz[0]), factor=float(a[nz[0]])))
        else:
            out.append(dict(row=i, kind="coupled", coef=a))
    return out
