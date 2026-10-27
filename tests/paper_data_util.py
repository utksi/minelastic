# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Helpers for the regression against the paper's benchmark runs (tests/data_paper)."""
import json
import os

import numpy as np
from ase import Atoms

from minelastic import ElasticSetup

DATA = os.path.join(os.path.dirname(__file__), "data_paper")
MATERIALS = ["albite", "ga2o3", "tisi2", "camoo4", "tio2", "dolomite", "sio2", "mg", "al", "cawo4", "bi", "batio3"]
W = np.array([1, 1, 1, np.sqrt(2), np.sqrt(2), np.sqrt(2)])


def load(name):
    """Fixture tests/data_paper/<name>.json."""
    return json.load(open(os.path.join(DATA, name + ".json")))


def setup_of(d):
    """ElasticSetup with the run's own patterns, step POTIM, in the input frame."""
    s = d["structure"]
    atoms = Atoms(s["symbols"], cell=s["cell"], scaled_positions=s["scaled_positions"], pbc=True)
    return ElasticSetup(atoms, scheme=d["patterns"], delta=d["potim"], stencil="central", frame="input")


def clamped(d, method=None):
    """Clamped-ion C (GPa, Voigt) from the raw +-POTIM stresses (J sigma, as VASP)."""
    method = method or ("uncoupled" if d["minimal"] else "lstsq")
    setup = setup_of(d)
    st = {f"s{c['pattern'] + 1}_{'p' if c['sign'] > 0 else 'm'}1": np.array(c["stress_GPa"]) * c["J"]
          for c in d["strained"]}
    return setup.compute(st, method=method), setup


def eig(C):
    """Sorted eigenvalues of the Mandel matrix W C W (frame-independent, compared across frames)."""
    return np.sort(np.linalg.eigvalsh(W[:, None] * np.asarray(C) * W[None, :]))


def paper_tensor(d, setup):
    """Paper tensor: from the 21 upper-triangle values, or from the labelled constants in the run's form."""
    if "paper_upper" in d:
        P = np.zeros((6, 6))
        k = 0
        for i in range(6):
            for j in range(i, 6):
                P[i, j] = P[j, i] = d["paper_upper"][k]
                k += 1
        return P
    if "paper" in d:
        return setup.form.matrix([d["paper"][lab] for lab in setup.form.labels])
    return None
