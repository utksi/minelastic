# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Space-group operations of a crystal as Cartesian rotations and atom permutations (spglib).

An operation {W|w} acts on fractional coordinates as x -> W x + w. With the lattice vectors as the rows of A,
its Cartesian rotation is R = A^T W A^-T, and it maps atom n onto the atom m = g(n) at W x_n + w (modulo
lattice vectors). The elastic symmetry type (one of nine, as in the paper) follows from the space-group
number.

With initial magnetic moments, the operations are those of the magnetic space group (spglib), with and
without time reversal: stress, strain, forces and polarization are even under time reversal, so the spatial
parts of all magnetic operations constrain C and the other tensors. Collinear moments are scalars, vector
moments are axial vectors. Atoms of one element that belong to different species of the DFT code (array
SPECIES of the Atoms object, set by the readers) are distinguished.
"""

SPECIES = "minelastic_species"
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import spglib

# elastic symmetry types of the paper, by space-group number
ELASTIC_TYPES = (
    (1, 2, "triclinic"),
    (3, 15, "monoclinic"),
    (16, 74, "orthorhombic"),
    (75, 88, "tetragonal-I"),
    (89, 142, "tetragonal-II"),
    (143, 148, "trigonal-I"),
    (149, 167, "trigonal-II"),
    (168, 194, "hexagonal"),
    (195, 230, "cubic"),
)


def elastic_type(spacegroup_number: int) -> str:
    """Elastic symmetry type of space group n: triclinic 1-2, monoclinic 3-15, orthorhombic 16-74,
    tetragonal-I 75-88, tetragonal-II 89-142, trigonal-I 143-148, trigonal-II 149-167, hexagonal 168-194,
    cubic 195-230."""
    for lo, hi, name in ELASTIC_TYPES:
        if lo <= spacegroup_number <= hi:
            return name
    raise ValueError(f"invalid space-group number {spacegroup_number}")


def _cell(atoms):
    """(A, x, numbers) for spglib: lattice vectors as rows, fractional positions in [0, 1), species numbers.

    Atoms of one element with different labels in atoms.arrays[SPECIES] (separate VASP species blocks, pw.x
    species labels) get different numbers, as they are different species in the DFT code.
    """
    numbers = atoms.numbers
    if SPECIES in atoms.arrays:
        keys = list(zip(numbers.tolist(), np.asarray(atoms.arrays[SPECIES]).tolist()))
        ids = {k: i + 1 for i, k in enumerate(dict.fromkeys(keys))}
        numbers = np.array([ids[k] for k in keys])
    return np.asarray(atoms.cell[:], float), atoms.get_scaled_positions(wrap=True), numbers


def _standard_rotation(lattice, P):
    """Rotation to spglib's standard orientation, x_std = R x. With the standardized basis
    (a_s b_s c_s) = (a b c) P^-1, the rows of R are e1 = a_s/|a_s|, e3 = a_s x b_s/|a_s x b_s|, e2 = e3 x e1:
    a_s along x, b_s in the xy-plane. For a non-magnetic dataset this equals spglib's std_rotation_matrix; the
    magnetic dataset leaves its standardized lattice in the input orientation, so R is built here."""
    Ls = lattice.T @ np.linalg.inv(P)               # columns a_s, b_s, c_s
    e1 = Ls[:, 0] / np.linalg.norm(Ls[:, 0])
    e3 = np.cross(Ls[:, 0], Ls[:, 1])
    e3 /= np.linalg.norm(e3)
    return np.array([e1, np.cross(e3, e1), e3])


@dataclass
class SymmetryInfo:
    """Symmetry of a crystal in its own Cartesian frame.

    number, international, pointgroup: space group and point group (spglib); rotations: the distinct
    Cartesian rotations R of the point group; operations: all space-group operations as (R, perm) with
    perm[n] = g(n), the atom onto which atom n is mapped; std_rotation: spglib's rotation from the input
    orientation to its standard orientation (x_std = std_rotation x).
    """

    number: int
    international: str
    pointgroup: str
    rotations: List[np.ndarray]
    operations: List[tuple] = field(repr=False)
    std_rotation: Optional[np.ndarray] = field(default=None, repr=False)

    @property
    def elastic_type(self) -> str:
        """Elastic symmetry type of the space group."""
        return elastic_type(self.number)


def analyze(atoms, symprec: float = 1e-5, angle_tolerance: float = -1.0) -> SymmetryInfo:
    """SymmetryInfo of atoms from spglib (tolerance symprec, in Angstrom).

    Each operation {W|w} gives R = A^T W A^-T and the permutation g with x_g(n) = W x_n + w (mod 1). The
    structure must be mapped onto itself within max(10 symprec, 1e-4) Angstrom, else RuntimeError. With
    nonzero initial magnetic moments, the operations of the magnetic space group are used (all of them, with
    and without time reversal) and the space-group number is that of the group spglib reports for it (the
    family space group of a type I-III group, the unitary subgroup of a type IV group), which has the point
    group of all spatial parts. The standard orientation is that of spglib's standardized cell
    (_standard_rotation).
    """
    lattice, frac, numbers = _cell(atoms)
    moments = atoms.get_initial_magnetic_moments()
    # spglib >= 2.5 returns objects, older versions dicts
    getter = lambda d: (lambda k: getattr(d, k)) if not isinstance(d, dict) else d.__getitem__  # noqa: E731
    if np.any(np.abs(moments) > 1e-8):
        ds = spglib.get_magnetic_symmetry_dataset((lattice, frac, numbers, moments), is_axial=moments.ndim == 2,
                                                  symprec=symprec, angle_tolerance=angle_tolerance)
        if ds is None:
            raise RuntimeError("spglib could not determine the magnetic symmetry")
        get = getter(ds)
        sg = getter(spglib.get_spacegroup_type(int(get("hall_number"))))
        number, international, pointgroup = sg("number"), sg("international_short"), sg("pointgroup_international")
        std_rot = _standard_rotation(lattice, np.asarray(get("transformation_matrix"), float))
    else:
        ds = spglib.get_symmetry_dataset((lattice, frac, numbers), symprec=symprec,
                                         angle_tolerance=angle_tolerance)
        if ds is None:
            raise RuntimeError("spglib could not determine the symmetry")
        get = getter(ds)
        number, international, pointgroup = get("number"), get("international"), get("pointgroup")
        std_rot = np.asarray(get("std_rotation_matrix"), float)
    LT = lattice.T                                  # columns = lattice vectors
    LTinv = np.linalg.inv(LT)
    ops, rots = [], []
    nat = len(numbers)
    tol = max(10 * symprec, 1e-4)
    for Wm, w in zip(get("rotations"), get("translations")):
        R = LT @ Wm @ LTinv                          # Cartesian rotation
        new = frac @ Wm.T + w
        perm = np.empty(nat, int)
        for n in range(nat):
            d = frac - new[n]
            d -= np.round(d)
            dist = np.linalg.norm(d @ lattice, axis=1)
            m = int(np.argmin(dist))
            if dist[m] > tol or numbers[m] != numbers[n]:          # tol in Angstrom, like symprec
                raise RuntimeError("symmetry operation does not map the structure onto itself; "
                                   "check symprec")
            perm[n] = m
        ops.append((R, perm))
        if not any(np.allclose(R, U, atol=1e-6) for U in rots):
            rots.append(R)
    return SymmetryInfo(number=int(number), international=str(international), pointgroup=str(pointgroup),
                        rotations=rots, operations=ops, std_rotation=std_rot)


def rotate_atoms(atoms, R):
    """Copy of atoms rotated rigidly, x' = R x: lattice vectors a_i' = R a_i, fractional coordinates kept, vector
    (noncollinear) initial moments m' = det(R) R m (axial vectors); collinear moments are unchanged."""
    R = np.asarray(R, float)
    new = atoms.copy()
    new.set_cell(np.asarray(atoms.cell[:]) @ R.T, scale_atoms=True)
    m = atoms.get_initial_magnetic_moments()
    if m.ndim == 2:
        new.set_initial_magnetic_moments(np.linalg.det(R) * m @ R.T)
    return new
