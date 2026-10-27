# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Regression against the benchmark runs of the paper (Table IV), one material per elastic symmetry type.

The fixtures (tests/data_paper, made by make_paper_fixtures.sh) hold the raw stresses of the +-POTIM strained
cells of each IBRION = 6 run, the tensors VASP printed, and the published constants. e1 = stock VASP (Table IV
"Standard"), e2 = LMINIMAL_ELASTIC (v1 patch, "This work"); VASP 6.5.1 on Dardel.
"""
import os

import numpy as np
import pytest

from minelastic.symmetry import ELASTIC_TYPES
from paper_data_util import DATA, MATERIALS, clamped, eig, load, paper_tensor

if not all(os.path.exists(os.path.join(DATA, f"{m}_{e}.json")) for m in MATERIALS for e in ("e1", "e2")):
    pytest.skip("the paper's benchmark data (tests/data_paper) are not included", allow_module_level=True)

CASES = [(m, e) for m in MATERIALS for e in ("e1", "e2")]
# constants whose sign depends on the in-plane orientation of the axes (a 180 degree rotation about z flips it)
ODD = {"c14", "c15", "c16"}


@pytest.mark.parametrize("material,run", CASES)
def test_clamped_tensor_from_raw_stresses(material, run):
    """minelastic reproduces the clamped-ion tensor VASP printed from the raw stresses of the strained cells.

    Stock runs: the same projection as VASP's TSYM4, so exact. Minimal runs: v1 completed the tensor with its
    per-class tables and TSYM4, minelastic averages the single-constant relations; both are estimators of the
    same tensor and differ by the noise of the redundant stress components (< 0.2 GPa).
    """
    d = load(f"{material}_{run}")
    res, setup = clamped(d)
    tol = 1e-3 if not d["minimal"] else 0.2
    assert np.abs(res.C - np.array(d["vasp_clamped"])).max() < tol
    if d["minimal"]:
        assert len(setup.patterns) == setup.n_min           # the minimal number of strained cells


@pytest.mark.parametrize("material,run", [c for c in CASES if "paper" in load(f"{c[0]}_{c[1]}")
                                          or "paper_upper" in load(f"{c[0]}_{c[1]}")])
def test_paper_table(material, run):
    """Clamped (minelastic, from raw stresses) + ionic relaxation (VASP) reproduces Table IV."""
    d = load(f"{material}_{run}")
    res, setup = clamped(d)
    total = res.C + np.array(d["vasp_ionic"])
    P = paper_tensor(d, setup)
    assert np.abs(eig(total) - eig(P)).max() < 0.2           # rotation invariant (frames may differ)
    if "paper" in d:                                           # same frame: constant by constant
        for p, lab in enumerate(setup.form.labels):
            got = setup.form.project(total)[p]
            want = d["paper"][lab]
            if lab in ODD:
                got, want = abs(got), abs(want)
            assert abs(got - want) < 0.16, (lab, got, want)


STOCK_VS_MINIMAL_TOL = {"al": 1.2}                              # metal: k-point error differs per strain pattern


@pytest.mark.parametrize("material", MATERIALS)
def test_minimal_matches_stock(material):
    """The minimal set and the six stock strains give the same clamped-ion tensor (finite-difference level)."""
    a, b = clamped(load(f"{material}_e2"))[0].C, clamped(load(f"{material}_e1"))[0].C
    assert np.abs(a - b).max() < STOCK_VS_MINIMAL_TOL.get(material, 0.6)


def test_all_nine_types_covered():
    """The benchmark materials span all nine elastic symmetry types."""
    types = {clamped(load(f"{m}_e1"))[1].symmetry.elastic_type for m in MATERIALS}
    assert types == {name for _, _, name in ELASTIC_TYPES}
