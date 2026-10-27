# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Regression against real VASP IBRION=6 runs (fixtures in tests/data, made with make_vasp_fixture.py).

For every run the clamped-ion tensor is recomputed by minelastic from the raw stresses of the +-POTIM strained
cells and compared with the tensor VASP printed. VASP differentiates J*sigma (SIF/Omega0), so the stresses are
multiplied by J = det F of the deformation VASP applied. Runs of the minimal strain set are compared with the
stock six-strain runs of the same input, and the paper's runs with its Table IV.
"""
import glob
import json
import os

import numpy as np
import pytest
from ase import Atoms

from minelastic import ElasticSetup
from minelastic.io.vasp import read_elastic_tensor  # noqa: F401  (same parser as the fixtures)
from minelastic.voigt import rotate_stiffness

DATA = os.path.join(os.path.dirname(__file__), "data")
FIXTURES = sorted(glob.glob(os.path.join(DATA, "*.json")))


def load(name):
    """Fixture tests/data/<name>."""
    return json.load(open(os.path.join(DATA, name)))


def recompute(d, method):
    """C from the raw stresses of a VASP run: the run's own patterns with step POTIM in the input frame, and the
    stresses multiplied by J = det F (VASP differentiates J sigma).
    """
    s = d["structure"]
    atoms = Atoms(s["symbols"], cell=s["cell"], scaled_positions=s["scaled_positions"], pbc=True)
    setup = ElasticSetup(atoms, scheme=d["patterns"], delta=d["potim"], stencil="central", frame="input")
    stresses = {f"s{c['pattern'] + 1}_{'p' if c['sign'] > 0 else 'm'}1": np.array(c["stress_GPa"]) * c["J"]
                for c in d["strained"]}
    return setup.compute(stresses, method=method)


@pytest.mark.parametrize("path", FIXTURES, ids=[os.path.basename(p)[:-5] for p in FIXTURES])
def test_reproduces_vasp_clamped_tensor(path):
    """The recomputed clamped-ion C equals VASP's SYMMETRIZED ELASTIC MODULI: stock runs to 2e-3 relative (least
    squares is VASP's projection), minimal-set runs to 0.05 relative.
    """
    d = json.load(open(path))
    ref = np.array(d["vasp_clamped"])
    tol = 2e-3 if not d["minimal"] else 0.05          # stock: the same projection as TSYM4
    for method in (("lstsq",) if not d["minimal"] else ("uncoupled", "lstsq")):
        C = recompute(d, method).C
        assert np.abs(C - ref).max() < tol * max(1.0, np.abs(ref).max() / 100), (method, np.abs(C - ref).max())


PAIRS = [("al_paper_minimal.json", "al_paper_stock.json", 1.5),   # metal: k-point error differs per pattern
         ("aln_std_v2.json", "aln_std_stock.json", 0.3),
         ("aln_cx_v2.json", "aln_cx_stock.json", 0.3),
         ("se_std_v2.json", "se_std_stock.json", 0.15),
         ("se_rot90_v2.json", "se_rot90_stock.json", 0.15)]


@pytest.mark.parametrize("minimal,stock,tol", PAIRS, ids=[p[0][:-5] for p in PAIRS])
def test_minimal_matches_stock(minimal, stock, tol):
    """Minimal set and six single strains of the same input agree within tol GPa."""
    a, b = recompute(load(minimal), "uncoupled").C, recompute(load(stock), "lstsq").C
    assert np.abs(a - b).max() < tol, np.abs(a - b).max()


@pytest.mark.parametrize("rot,std,R", [
    ("aln_cx_v2.json", "aln_std_v2.json", np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)),
    ("se_rot90_v2.json", "se_std_v2.json", np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float))])
def test_rotated_cells_give_rotated_tensor(rot, std, R):
    """The c||x AlN and the 90-degree rotated Se cell give the rotated tensor of the standard cell."""
    a = recompute(load(rot), "uncoupled").C
    b = rotate_stiffness(recompute(load(std), "uncoupled").C, R)
    assert np.abs(a - b).max() < 0.1, np.abs(a - b).max()


@pytest.mark.parametrize("name", ["al_paper_stock.json", "al_paper_minimal.json"])
def test_paper_table(name):
    """The Al runs of the paper: VASP's total tensor and minelastic's recomputation both round to Table IV."""
    d = load(name)
    C = np.array(d["vasp_total"])
    got = dict(c11=C[0, 0], c12=C[0, 1], c13=C[0, 2], c33=C[2, 2], c44=C[3, 3], c66=C[5, 5])
    for k, v in d["paper"].items():
        assert abs(got[k] - v) < 0.051, (k, got[k], v)
    # and minelastic's own extraction from the raw stresses rounds to the published numbers
    C2 = recompute(d, "uncoupled" if d["minimal"] else "lstsq").C
    got2 = dict(c11=C2[0, 0], c12=C2[0, 1], c44=C2[3, 3])
    for k, v in d["paper"].items():
        assert abs(got2[k] - v) < 0.051, (k, got2[k], v)


def test_berry_phase_piezo_matches_dfpt():
    """Clamped-ion piezo tensor of AlN from the Berry-phase dipoles of the 4 minimal strained cells (VASP
    LCALCPOL) agrees with the DFPT (LEPSILON) tensor of the same setup to ~1%."""
    d = json.load(open(os.path.join(os.path.dirname(__file__), "data_berry", "aln_berry_minimal.json")))
    setup = ElasticSetup.from_dict({**d["setup"], "program": "minelastic", "version": "x",
                                    "constants": [], "spacegroup": "", "spacegroup_number": 186,
                                    "elastic_type": "hexagonal"})
    r = setup.compute(d["stress_GPa"], polarization=d["dipole_eA"])
    e, ref = r.piezo, d["dfpt_clamped_piezo_C_m2"]
    assert abs(e[2, 2] - ref["e33"]) < 0.015 and abs(e[2, 0] - ref["e31"]) < 0.015
    assert abs(e[0, 4] - ref["e15"]) < 0.015 and abs(e[1, 3] - ref["e15"]) < 0.015
    assert np.abs(e[:, :]).max() < 0.5 and abs(e[2, 0] - e[2, 1]) < 1e-6      # 6mm form
