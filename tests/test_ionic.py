# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Ionic route (minelastic.ionic): force constants, Born charges, ionic contributions, total e and d.

* EMT hcp Cu: the symmetry-completed force constants equal those of all displacements, and the relaxed-ion
  constants from clamped cells plus displacements equal those from relaxing the ions in every strained cell.
* VASP, wurtzite AlN (tests/data_ionic; make_ionic_fixture.py, make_aln_run_fixture.py): raw stresses, forces
  and Berry-phase dipoles of the strained and displaced cells give the ionic contribution to C, the Born
  charges, the ionic and total piezoelectric tensors and d as stock VASP (IBRION = 6, LEPSILON) does.
"""
import json
import os

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.io import read, write

import minelastic
from minelastic import ElasticSetup, ionic
from minelastic.io import vasp
from minelastic.symmetry import analyze
from minelastic.workflow import collect, generate
from test_minelastic import _fake_vasprun

DATA = os.path.join(os.path.dirname(__file__), "data_ionic", "aln_ionic.json")
GPA = 160.21766208


def hcp_cu():
    """hcp Cu (two atoms per cell, internal relaxation under xx - yy and xy strains)."""
    return bulk("Cu", "hcp", a=2.55, c=4.16)


def emt_forces(atoms, n, d, h=1e-3):
    """dF/du of atom n along d by central differences with EMT: [F(+h d) - F(-h d)] / (2h)."""
    out = []
    for s in (1, -1):
        a = atoms.copy()
        a.calc = EMT()
        a.positions[n] += s * h * d
        out.append(a.get_forces())
    return (out[0] - out[1]) / (2 * h)


def test_force_constants_completion():
    """hcp Cu: 2 irreducible displacements (atom 1, x and z) completed by symmetry give the force constants of
    displacing every atom along x, y and z (symmetrized), to 1e-8 eV/A^2.
    """
    atoms = hcp_cu()
    sym = analyze(atoms)
    disps = ionic.irreducible_displacements(sym.operations, len(atoms))
    assert len(disps) == 2                                   # one atom, x and z (as VASP)
    Phi = ionic.force_constants(sym.operations, len(atoms), disps, [emt_forces(atoms, n, d) for n, d in disps])
    N = len(atoms)
    full = np.zeros((N, 3, N, 3))
    for n in range(N):
        for j in range(3):
            full[:, :, n, j] = -emt_forces(atoms, n, np.eye(3)[j])
    full = 0.5 * (full.reshape(3 * N, 3 * N) + full.reshape(3 * N, 3 * N).T).reshape(N, 3, N, 3)
    assert np.abs(Phi - full).max() < 1e-8


def test_relaxed_ion_constants_from_displacements():
    """EMT hcp Cu: C0 + C_ion from clamped cells equals C from relaxing the ions in every strained cell
    (0.05 GPa); Lambda from the stresses of the displaced cells equals Lambda from the forces in the strained
    cells (1e-3 eV/A); ionic and relax_ions together are refused.
    """
    atoms = hcp_cu()
    r = minelastic.calculate(atoms, EMT(), ionic=True, delta=0.005, displacement=0.005, method="lstsq")
    rr = minelastic.calculate(atoms, EMT(), relax_ions=True, delta=0.005, fmax=1e-5, method="lstsq")
    assert np.abs(r.C_ionic).max() > 10                      # hcp: a real internal relaxation
    assert np.abs(r.C_relaxed - rr.C).max() < 0.05
    # the internal strain tensor from the stresses of the displaced cells equals the one from the forces
    assert np.abs(r.internal_strain - r.internal_strain_displaced).max() < 1e-3
    with pytest.raises(ValueError):
        minelastic.calculate(atoms, EMT(), ionic=True, relax_ions=True)


def aln_case():
    """Setup and per-state stresses, forces and dipoles of the stored VASP AlN run (IBRION = 6 strains and
    displacements, POTIM = 0.015), with the stock VASP reference tensors.
    """
    d = json.load(open(DATA))
    s = d["structure"]
    atoms = Atoms(s["symbols"], cell=s["cell"], scaled_positions=s["scaled_positions"], pbc=True)
    setup = ElasticSetup(atoms, scheme=d["patterns"], delta=d["potim"], stencil="central", frame="input",
                         ionic=True, displacement=d["potim"])
    S, F, P = {}, {}, {}
    for c in d["strained"]:
        sid = f"s{c['pattern'] + 1}_{'p' if c['sign'] > 0 else 'm'}1"
        S[sid], P[sid], F[sid] = c["stress_GPa"], c["dipole_eA"], d["reference_forces"]
    for c in d["displaced"]:
        u = np.array(c["displacement_A"])
        k = next(i for i, (n, v) in enumerate(setup.displacements) if n == c["atom"] and abs(abs(u @ v) - d["potim"]) < 1e-6)
        sid = f"d{k + 1}_{'p' if u @ setup.displacements[k][1] > 0 else 'm'}1"
        S[sid], F[sid], P[sid] = c["stress_GPa"], c["forces"], c["dipole_eA"]
    return setup, S, F, P, {k: np.array(v) for k, v in d["stock"].items()}


def test_vasp_aln_ionic_route():
    """VASP AlN: minelastic chooses VASP's displacements; C_ion agrees with VASP to 0.005 GPa, the Berry-phase
    Born charges with DFPT to 0.01 |e| (acoustic sum < 0.002), e and d with the DFPT route to 0.01 C/m^2 and
    0.1 pC/N; with the DFPT Born charges, e_ion equals VASP's ionic piezoelectric term to 1e-4 C/m^2.
    """
    setup, S, F, P, ref = aln_case()
    # the displacements VASP chose are the irreducible ones of minelastic
    assert [(n, tuple(v)) for n, v in setup.displacements] == [(0, (1, 0, 0)), (0, (0, 0, 1)), (2, (1, 0, 0)), (2, (0, 0, 1))]
    r = setup.compute(S, forces=F, method="lstsq", stress_measure="kirchhoff", polarization=P)
    assert np.abs(r.C_ionic - ref["C_ionic"]).max() < 0.005          # VASP ELASTIC_ION: 77.7 GPa
    assert np.abs(r.born - ref["born_dfpt"]).max() < 0.01            # Berry phase vs DFPT (2.68 |e|)
    assert np.abs(r.born.sum(axis=0)).max() < 0.002                   # acoustic sum rule
    e_ref = ref["piezo_clamped_dfpt"] + ref["piezo_ionic_dfpt_born"]
    d_ref = ionic.strain_coefficients(e_ref, ref["C_total"])
    assert np.abs(r.piezo_total - e_ref).max() < 0.01                # e33 = 1.454 C/m^2
    assert np.abs(r.piezo_d - d_ref).max() < 0.1                     # d33 = 5.30 pC/N
    # with the DFPT Born charges, the ionic piezoelectric term is VASP's PIEZO_ION
    rd = setup.compute(S, forces=F, method="lstsq", stress_measure="kirchhoff", born=ref["born_dfpt"],
                       clamped_piezo=ref["piezo_clamped_dfpt"])
    assert np.abs(rd.piezo_ionic - ref["piezo_ionic_dfpt_born"]).max() < 1e-4
    assert np.abs(rd.piezo_d - d_ref).max() < 0.01
    assert rd.piezo_clamped_source == "external" and r.piezo_clamped_source == "strained cells"
    txt = r.summary()
    assert "relaxed-ion stiffness" in txt and "pC/N" in txt and "Born effective charges" in txt


def test_documented_vasp_run():
    """The --ionic run of docs/how-to/piezoelectricity.md (VASP 6.6.1, data_ionic/aln_run_ionic.json): the numbers
    the docs print, and VASP's own IBRION = 6 + LEPSILON within the step-size differences: C to 0.7 GPa (C_ion to 0.06),
    the internal strain to 0.04 eV/A from the strained cells and 0.06 eV/A from the displaced atoms, Z* to
    0.004 |e|, e to 0.006 C/m^2, d to 0.06 pC/N."""
    d = json.load(open(os.path.join(os.path.dirname(__file__), "data_ionic", "aln_run_ionic.json")))
    r = ElasticSetup.from_dict(d["setup"]).compute(d["stress_GPa"], forces=d["forces_eV_A"],
                                                   polarization=d["dipole_eA"])
    ref = {k: np.array(v) for k, v in d["stock"].items()}
    e_ref = ref["piezo_clamped_dfpt"] + ref["piezo_ionic_dfpt_born"]
    printed = [r.C_relaxed[0, 0], r.C_relaxed[2, 2], r.C_relaxed[3, 3], r.piezo_total[2, 2], r.piezo_d[2, 2]]
    assert np.allclose(printed, [376.11, 354.60, 112.47, 1.450, 5.302], atol=5e-3)
    assert np.abs(r.C_relaxed - ref["C_total"]).max() < 0.7 and np.abs(r.C_ionic - ref["C_ionic"]).max() < 0.06
    assert np.abs(r.internal_strain - ref["internal_strain"]).max() < 0.04
    assert np.abs(r.internal_strain_displaced - ref["internal_strain"]).max() < 0.06
    assert np.abs(r.born - ref["born_dfpt"]).max() < 0.004
    assert np.abs(r.piezo_total - e_ref).max() < 0.006
    assert np.abs(r.piezo_d - ionic.strain_coefficients(e_ref, ref["C_total"])).max() < 0.06


def test_piezo_and_born_rotations():
    """e is invariant under the point-group rotations; rotate_piezo and rotate_born are inverted by R^T."""
    setup, S, F, P, ref = aln_case()
    r = setup.compute(S, forces=F, method="lstsq", stress_measure="kirchhoff", polarization=P)
    for R in setup.symmetry.rotations:                                # invariant under the point group
        assert np.abs(ionic.rotate_piezo(r.piezo_total, R) - r.piezo_total).max() < 1e-8
    R = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)
    back = ionic.rotate_piezo(ionic.rotate_piezo(r.piezo_total, R), R.T)
    assert np.abs(back - r.piezo_total).max() < 1e-12
    assert np.abs(ionic.rotate_born(ionic.rotate_born(r.born, R), R.T) - r.born).max() < 1e-12


def test_vasp_readers_born_and_piezo(tmp_path):
    """read_born_charges and read_clamped_piezo on an OUTCAR excerpt (columns XY YZ ZX to Voigt yz xz xy)."""
    out = tmp_path / "OUTCAR"
    out.write_text(""" BORN EFFECTIVE CHARGES (including local field effects) (in |e|, cummulative output)
 ---------------------------------------------------------------------------------
 ion    1
    1     2.51965     0.00000     0.00000
    2     0.00000     2.51965     0.00000
    3     0.00000     0.00000     2.67975
 ion    2
    1    -2.51965     0.00000     0.00000
    2     0.00000    -2.51965     0.00000
    3     0.00000     0.00000    -2.67975

 PIEZOELECTRIC TENSOR (including local field effects)  for field in x, y, z        (C/m^2)
               XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------
  x     0.00000     0.00000     0.00000     0.00000     0.00000     0.34039
  y     0.00000     0.00000     0.00000     0.00000     0.34039     0.00000
  z     0.27039     0.27039    -0.42396     0.00000     0.00000     0.00000
""")
    Z = vasp.read_born_charges(out)
    assert Z.shape == (2, 3, 3) and Z[0, 2, 2] == 2.67975
    e = vasp.read_clamped_piezo(out)                                  # Voigt xx yy zz yz xz xy
    assert e[0, 4] == 0.34039 and e[1, 3] == 0.34039 and e[2, 2] == -0.42396


def test_ionic_workflow_vasp(tmp_path):
    """generate --ionic writes the displaced-atom cells; collect gives the relaxed-ion constants."""
    tpl = tmp_path / "ref"
    tpl.mkdir()
    atoms = hcp_cu()
    write(tpl / "POSCAR", atoms, format="vasp", direct=True)
    (tpl / "INCAR").write_text("ENCUT = 400\n")
    (tpl / "POTCAR").write_text("dummy\n")
    (tpl / "KPOINTS").write_text("auto\n0\nGamma\n4 4 4\n")
    setup, _ = generate("vasp", str(tpl), str(tmp_path / "el"), delta=0.005, ionic=True, displacement=0.005)
    ids = [st["id"] for st in setup.states]
    assert "d1_p1" in ids and "d2_m1" in ids
    back = read(tmp_path / "el" / "d1_p1" / "POSCAR", format="vasp")
    assert np.allclose(back.positions, setup.strained_atoms("d1_p1").positions, atol=1e-10)
    for st in setup.states:
        a = setup.strained_atoms(st)
        a.calc = EMT()
        _fake_vasprun(tmp_path / "el" / st["id"] / "vasprun.xml", a.get_stress(voigt=True) * GPA, a.get_forces())
    r = collect(str(tmp_path / "el"), method="lstsq")
    rr = minelastic.calculate(atoms, EMT(), relax_ions=True, delta=0.005, fmax=1e-5, method="lstsq")
    assert np.abs(r.C_relaxed - rr.C).max() < 0.05
    d = json.load(open(tmp_path / "el" / "elastic.json"))
    assert "C_relaxed" in d and "force_constants_eV_per_A2" in d
