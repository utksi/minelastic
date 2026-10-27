# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Tests of the strain sets, the extraction, the frames, the workflow and the code interfaces.

Exact tests use synthetic stresses sigma = C e + T(e, e) with a random symmetry-adapted C and a random
quadratic term T, which central differences cancel exactly; the extracted C must equal the input to 1e-8
relative. EMT tests compare the minimal set with the six single strains for a real (anharmonic) potential.
"""
import json
import os

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.io import read, write
from ase.spacegroup import crystal

import minelastic
from minelastic import ElasticSetup, algorithm
from minelastic.io import qe, vasp
from minelastic.symmetry import rotate_atoms
from minelastic.tensor_form import TensorForm, invariant_strain_dimension
from minelastic.voigt import mandel_rotation, rotate_stiffness, stiffness_from_mandel
from minelastic.workflow import collect, generate

BASIS = [(0.1123, 0.2234, 0.3345), (0.4111, 0.0789, 0.2222)]
# representative space groups of the nine elastic symmetry types, the paper's strain supports
# (Voigt 0..5 = xx yy zz yz xz xy) and its extraction formulas (row -> constant) per pattern
PAPER = {
    "triclinic": (2, [4.1, 4.7, 5.3, 81, 97, 103], [(0,), (1,), (2,), (3,), (4,), (5,)],
                  [{0: "c11", 1: "c12", 2: "c13", 3: "c14", 4: "c15", 5: "c16"},
                   {1: "c22", 2: "c23", 3: "c24", 4: "c25", 5: "c26"}, {2: "c33", 3: "c34", 4: "c35", 5: "c36"},
                   {3: "c44", 4: "c45", 5: "c46"}, {4: "c55", 5: "c56"}, {5: "c66"}]),
    "monoclinic": (12, [12.2, 3.0, 5.8, 90, 103.8, 90], [(0, 3), (1, 5), (2,), (4,)],
                   [{0: "c11", 1: "c12", 2: "c13", 3: "c44", 4: "c15", 5: "c46"},
                    {1: "c22", 2: "c23", 4: "c25", 5: "c66"}, {2: "c33", 4: "c35"}, {4: "c55"}]),
    "orthorhombic": (70, [8.3, 4.8, 8.6, 90, 90, 90], [(0, 3, 4, 5), (1,), (2,)],
                     [{0: "c11", 1: "c12", 2: "c13", 3: "c44", 4: "c55", 5: "c66"}, {1: "c22", 2: "c23"},
                      {2: "c33"}]),
    "tetragonal-I": (88, [5.2, 5.2, 11.4, 90, 90, 90], [(0, 3), (2, 5)],
                     [{0: "c11", 1: "c12", 2: "c13", 3: "c44", 5: "c16"}, {2: "c33", 5: "c66"}]),
    "tetragonal-II": (136, [4.6, 4.6, 2.96, 90, 90, 90], [(0, 3, 5), (2,)],
                      [{0: "c11", 1: "c12", 2: "c13", 3: "c44", 5: "c66"}, {2: "c33"}]),
    "trigonal-I": (148, [4.8, 4.8, 16.0, 90, 90, 120], [(0,), (2, 3)],
                   [{0: "c11", 1: "c12", 2: "c13", 3: "c14", 4: "c15"}, {2: "c33", 3: "c44"}]),
    "trigonal-II": (154, [4.9, 4.9, 5.4, 90, 90, 120], [(0,), (2, 3)],
                    [{0: "c11", 1: "c12", 2: "c13", 3: "c14"}, {2: "c33", 3: "c44"}]),
    "hexagonal": (194, [3.2, 3.2, 5.2, 90, 90, 120], [(0, 3), (2,)],
                  [{0: "c11", 1: "c12", 2: "c13", 3: "c44"}, {2: "c33"}]),
    "cubic": (225, [4.04, 4.04, 4.04, 90, 90, 90], [(0, 3)], [{0: "c11", 1: "c12", 3: "c44"}]),
}


def structure(kind):
    """Si-O test crystal of the representative space group of an elastic type (PAPER table)."""
    sg, cp = PAPER[kind][:2]
    return crystal(["Si", "O"], basis=BASIS, spacegroup=sg, cellpar=cp)


def random_rotation(seed):
    """Random proper rotation (det = +1) from the QR decomposition of a Gaussian matrix."""
    q = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))[0]
    return q * np.sign(np.linalg.det(q))


def random_stiffness(form, seed=0):
    """Random positive-definite stiffness matrix with the symmetry of ``form``."""
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(6, 6))
    M = A @ A.T + 6 * np.eye(6)                                      # Mandel, positive definite
    Qs = [mandel_rotation(R) for R in form.rotations]
    M = np.mean([Q @ M @ Q.T for Q in Qs], axis=0) * 50              # group average stays positive definite
    c = form.project(stiffness_from_mandel(M))
    return c, form.matrix(c)


def synthetic_stresses(setup, C, seed=1):
    """sigma = C e + quadratic term (cancels in central differences), in the working frame."""
    rng = np.random.default_rng(seed)
    T = rng.normal(size=(6, 6, 6)) * 300
    out = {}
    for st in setup.states:
        e = np.array(st["strain"])
        out[st["id"]] = C @ e + np.einsum("ijk,j,k->i", T, e, e)
    return out


@pytest.mark.parametrize("kind", list(PAPER))
def test_paper_tables(kind):
    """The minimal set has the paper's supports and single-constant relations for every elastic type, its size
    equals dim V0, and the exhaustive search finds no smaller complete set.
    """
    setup = ElasticSetup(structure(kind))
    assert setup.symmetry.elastic_type == kind
    sg, cp, supports, formulas = PAPER[kind]
    got = [tuple(np.flatnonzero(u)) for u in setup.patterns]
    assert got == supports
    assert len(supports) == setup.n_min == invariant_strain_dimension(setup.symmetry.rotations)
    for rel, formula in zip(setup.relation_table(), formulas):
        singles = {r["row"]: setup.form.labels[r["constant"]] for r in rel if r["kind"] == "single"}
        for row, lab in formula.items():
            assert singles.get(row) == lab, (kind, row, lab, singles)
    # the greedy construction is minimal: exhaustive search finds no smaller set
    assert len(algorithm.exact_minimum(setup.form)) == len(supports)


@pytest.mark.parametrize("kind", list(PAPER))
@pytest.mark.parametrize("method", ["uncoupled", "designated", "lstsq"])
def test_exact_recovery(kind, method):
    """Synthetic stresses with quadratic terms give C exactly (central differences), for every extraction method."""
    setup = ElasticSetup(structure(kind), delta=0.01)
    c, C = random_stiffness(setup.form)
    r = setup.compute(synthetic_stresses(setup, C), method=method)
    assert np.allclose(r.constants, c, atol=1e-8 * np.abs(c).max())
    assert np.allclose(r.C, C, atol=1e-8 * np.abs(C).max())


@pytest.mark.parametrize("kind", list(PAPER))
def test_rotated_input(kind):
    """A randomly oriented cell is strained in the standard frame; the input-frame tensor is consistent."""
    ref = ElasticSetup(structure(kind))
    c, C = random_stiffness(ref.form)
    R = random_rotation(7)
    setup = ElasticSetup(rotate_atoms(structure(kind), R))
    # triclinic: the form is the same in every frame, so the input orientation is kept
    assert np.allclose(setup.rotation, np.eye(3)) == (kind == "triclinic")
    assert len(setup.patterns) == setup.n_min
    # the true tensor in the working frame of this setup: rotate from the reference frame
    R_tot = setup.rotation @ R @ ref.rotation.T
    C_work = rotate_stiffness(C, R_tot)
    r = setup.compute(synthetic_stresses(setup, C_work))
    assert np.allclose(r.C, C_work, atol=1e-8 * np.abs(C).max())
    assert np.allclose(r.C_input_frame, rotate_stiffness(C, R @ ref.rotation.T), atol=1e-8 * np.abs(C).max())


def test_trigonal_312_orientation():
    """P312: 2-fold axes perpendicular to a, so in the a || x setting the form has c15 instead of c14."""
    at = crystal(["Si", "O"], basis=BASIS, spacegroup=149, cellpar=[5.0, 5.0, 4.5, 90, 90, 120])
    setup = ElasticSetup(at)
    assert "c15" in setup.form.labels and "c14" not in setup.form.labels
    assert len(setup.patterns) == 2
    c, C = random_stiffness(setup.form)
    assert np.allclose(setup.compute(synthetic_stresses(setup, C)).C, C)


def test_asess_scheme():
    """The six single strains with the four-point stencil (24 cells) recover C exactly."""
    setup = ElasticSetup(structure("hexagonal"), scheme="asess", stencil="four-point")
    assert len(setup.patterns) == 6 and len(setup.states) == 24
    c, C = random_stiffness(setup.form)
    assert np.allclose(setup.compute(synthetic_stresses(setup, C)).C, C)


def test_md_samples_errors():
    """Noisy samples (400 per state): the constants agree with the true ones within five standard errors."""
    setup = ElasticSetup(bulk("Al", "fcc", a=4.04, cubic=True))
    c, C = random_stiffness(setup.form)
    rng = np.random.default_rng(3)
    clean = synthetic_stresses(setup, C)
    noisy = {k: v + rng.normal(scale=0.05, size=(400, 6)) for k, v in clean.items()}
    r = setup.compute(noisy)
    assert r.errors is not None
    assert np.all(np.abs(r.constants - c) < 5 * r.errors + 1e-12)


def test_emt_internal_strain_and_relaxed():
    """EMT hcp Cu: minimal set and six single strains agree for clamped and relaxed ions (0.2 GPa) and for the
    internal-strain tensor (1e-3 eV/A); relaxing the ions lowers c11.
    """
    hcp = bulk("Cu", "hcp", a=2.55, c=4.16)
    r_min = minelastic.calculate(hcp, EMT(), delta=0.005)
    r_all = minelastic.calculate(hcp, EMT(), delta=0.005, scheme="asess")
    assert np.allclose(r_min.constants, r_all.constants, atol=0.2)
    assert np.abs(r_all.internal_strain).max() > 1.0
    assert np.allclose(r_min.internal_strain, r_all.internal_strain, atol=1e-3)
    rr_min = minelastic.calculate(hcp, EMT(), delta=0.005, relax_ions=True, fmax=1e-6)
    rr_all = minelastic.calculate(hcp, EMT(), delta=0.005, relax_ions=True, fmax=1e-6, scheme="asess")
    assert np.allclose(rr_min.constants, rr_all.constants, atol=0.2)
    assert rr_min.as_dict()["c11"] < r_min.as_dict()["c11"] - 5        # ionic relaxation softens c11


def test_json_roundtrip(tmp_path):
    """to_dict / from_dict restores the frame rotation, the patterns and the state ids (rotated tetragonal cell)."""
    setup = ElasticSetup(rotate_atoms(structure("tetragonal-I"), random_rotation(3)), stencil="four-point")
    d = json.loads(json.dumps(setup.to_dict()))
    s2 = ElasticSetup.from_dict(d)
    assert np.allclose(s2.patterns, setup.patterns) and np.allclose(s2.rotation, setup.rotation)
    assert [s["id"] for s in s2.states] == [s["id"] for s in setup.states]


# --------------------------------------------------------------------------------------------- VASP
def _fake_vasprun(path, sigma_gpa, forces):
    """Minimal vasprun.xml with one step: stress written as VASP does, -10 sigma in kB (XX YY ZZ / XY YZ ZX
    as a 3x3 block), and the forces.
    """
    s = -10 * np.array([[sigma_gpa[0], sigma_gpa[5], sigma_gpa[4]], [sigma_gpa[5], sigma_gpa[1], sigma_gpa[3]],
                        [sigma_gpa[4], sigma_gpa[3], sigma_gpa[2]]])          # kB, pressure-like
    rows = "\n".join(f"   <v> {a:16.8f} {b:16.8f} {c:16.8f} </v>" for a, b, c in s)
    frows = "\n".join(f"   <v> {a:16.8f} {b:16.8f} {c:16.8f} </v>" for a, b, c in forces)
    open(path, "w").write(f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
 <calculation>
  <varray name="forces" >
{frows}
  </varray>
  <varray name="stress" >
{rows}
  </varray>
 </calculation>
</modeling>
""")


def test_vasp_workflow(tmp_path):
    """generate vasp from an IBRION = 6 template (rotated hexagonal cell, KSPACING): INCAR tags replaced, POSCAR
    equal to strained_atoms; collect of synthetic vasprun.xml files returns C (to the 8 written decimals).
    """
    tpl = tmp_path / "ref"
    tpl.mkdir()
    atoms = rotate_atoms(structure("hexagonal"), random_rotation(5))
    write(tpl / "POSCAR", atoms, format="vasp", direct=True)
    (tpl / "INCAR").write_text("ENCUT = 520 ! cutoff\nIBRION = 6 ; ISIF = 3\nPOTIM = 0.015\nLMINIMAL_ELASTIC = .TRUE.\n"
                               "KSPACING = 0.2\nMAGMOM = 4*0\n")
    (tpl / "POTCAR").write_text("dummy\n")
    setup, notes = generate("vasp", str(tpl), str(tmp_path / "el"))
    assert any("KSPACING" in n for n in notes)
    tags = dict(vasp.read_incar(tmp_path / "el" / "s1_p1" / "INCAR"))
    assert tags["IBRION"] == "-1" and tags["ISIF"] == "2" and "POTIM" not in tags and "KSPACING" not in tags
    assert tags["MAGMOM"] == "4*0" and tags["ENCUT"] == "520"
    from ase.io import read
    back = read(tmp_path / "el" / "s2_m1" / "POSCAR", format="vasp")
    ref = setup.strained_atoms("s2_m1")
    assert back.get_chemical_symbols() == ref.get_chemical_symbols()
    assert np.allclose(back.cell[:], ref.cell[:], atol=1e-12) and np.allclose(back.positions, ref.positions, atol=1e-10)
    c, C = random_stiffness(setup.form)
    sig = synthetic_stresses(setup, C)
    for st in setup.states:
        _fake_vasprun(tmp_path / "el" / st["id"] / "vasprun.xml", sig[st["id"]], np.zeros((len(atoms), 3)))
    r = collect(str(tmp_path / "el"))
    assert np.allclose(r.C, C, atol=1e-5)                 # 8 decimals in kB
    assert os.path.exists(tmp_path / "el" / "elastic.json")


def test_vasp_elastic_tensor_block(tmp_path):
    """read_elastic_tensor: kBar to GPa and VASP order XX YY ZZ XY YZ ZX to Voigt (C44 = YZ, C66 = XY)."""
    block = """ TOTAL ELASTIC MODULI (kBar)
 Direction    XX          YY          ZZ          XY          YZ          ZX
 --------------------------------------------------------------------------------
 XX        3763.2529   1285.5238    977.9161      0.0000      0.0000     -0.0000
 YY        1285.5238   3763.2529    977.9161      0.0000      0.0000      0.0000
 ZZ         977.9161    977.9161   3552.4042      0.0000      0.0000     -0.0000
 XY           0.0000      0.0000      0.0000   1238.8646     -0.0000     -0.0000
 YZ          -0.0000      0.0000      0.0000      0.0000   1125.6858      0.0000
 ZX          -0.0000      0.0000     -0.0000     -0.0000      0.0000   1125.6858
"""
    p = tmp_path / "OUTCAR"
    p.write_text(block)
    C = vasp.read_elastic_tensor(str(p))
    assert np.isclose(C[0, 0], 376.32529) and np.isclose(C[3, 3], 112.56858) and np.isclose(C[5, 5], 123.88646)


# ----------------------------------------------------------------------------------------------- QE
QE_TEMPLATE = """&control
  calculation='scf', prefix='x', pseudo_dir='./', outdir='./tmp'
/
&system
  ibrav = 4, celldm(1) = 5.911804, celldm(3) = 1.603472, nat = 4, ntyp = 2,
  ecutwfc = 50, ecutrho = 400
/
&electrons
  conv_thr = 1.0d-10
/
ATOMIC_SPECIES
Al1 26.982 Al.UPF
N  14.007 N.UPF
ATOMIC_POSITIONS crystal
Al1 0.3333333333 0.6666666667 0.0000 0 0 1
Al1 0.6666666667 0.3333333333 0.5000
N  0.3333333333 0.6666666667 0.3818
N  0.6666666667 0.3333333333 0.8818
K_POINTS automatic
8 8 5 0 0 0
"""


def _fake_pwout(path, sigma_gpa, nat):
    """Minimal pw.x output: forces block and 'total   stress' block in Ry/bohr^3, pressure-like (-sigma)."""
    s = -np.array([[sigma_gpa[0], sigma_gpa[5], sigma_gpa[4]], [sigma_gpa[5], sigma_gpa[1], sigma_gpa[3]],
                   [sigma_gpa[4], sigma_gpa[3], sigma_gpa[2]]]) / qe.RY_PER_BOHR3_IN_GPA
    lines = [f"     number of atoms/cell      =            {nat}", "     Forces acting on atoms (cartesian axes, Ry/au):", ""]
    lines += [f"     atom    {n + 1} type  1   force =     0.00000000    0.00000000    0.00000000" for n in range(nat)]
    lines += ["", "     total   stress  (Ry/bohr**3)                   (kbar)     P=        0.00"]
    lines += ["  " + " ".join(f"{x:16.12f}" for x in row) + "   0.00 0.00 0.00" for row in s]
    open(path, "w").write("\n".join(lines) + "\n")


def test_qe_workflow(tmp_path):
    """generate qe from an ibrav = 4 template: ibrav = 0, the cell equals strained_atoms, labels and if_pos kept;
    collect of synthetic pw.x outputs returns C.
    """
    tpl = tmp_path / "pw.in"
    tpl.write_text(QE_TEMPLATE)
    at = qe.read_structure(str(tpl))
    assert np.isclose(at.cell.cellpar()[0], 5.911804 * qe.BOHR) and np.isclose(at.cell.cellpar()[5], 120)
    setup, _ = generate("qe", str(tpl), str(tmp_path / "el"))
    txt = (tmp_path / "el" / "s1_p1" / "pw.in").read_text()
    assert "ibrav = 0" in txt and "celldm" not in txt and "tstress = .true." in txt
    assert "Al1" in txt and "0 0 1" in txt
    gen = qe.read_structure(str(tmp_path / "el" / "s1_p1" / "pw.in"))
    assert np.allclose(gen.cell[:], setup.strained_atoms("s1_p1").cell[:], atol=1e-9)
    c, C = random_stiffness(setup.form)
    sig = synthetic_stresses(setup, C)
    for st in setup.states:
        _fake_pwout(tmp_path / "el" / st["id"] / "pw.out", sig[st["id"]], 4)
    r = collect(str(tmp_path / "el"))
    assert np.allclose(r.C, C, atol=1e-6)


@pytest.mark.parametrize("ibrav,cd,expect", [
    (2, {"celldm(1)": "7.0"}, np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]]) * 3.5),
    (3, {"celldm(1)": "7.0"}, np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]]) * 3.5),
    (4, {"celldm(1)": "6.0", "celldm(3)": "1.6"}, np.array([[6, 0, 0], [-3, 3 * np.sqrt(3), 0], [0, 0, 9.6]])),
])
def test_qe_latgen(ibrav, cd, expect):
    """qe_lattice reproduces the pw.x lattice vectors for ibrav = 2, 3 and 4."""
    system = [("ibrav", str(ibrav))] + list(cd.items())
    assert np.allclose(qe.qe_lattice(system), expect * qe.BOHR)


def test_snapshot_samples(tmp_path):
    """Finite-T protocol: several single points (snapshots) per strain state, averaged with errors."""
    tpl = tmp_path / "ref"
    tpl.mkdir()
    write(tpl / "POSCAR", bulk("Al", "fcc", a=4.04, cubic=True), format="vasp", direct=True)
    (tpl / "INCAR").write_text("ENCUT = 300\n")
    (tpl / "KPOINTS").write_text("a\n0\nGamma\n4 4 4\n0 0 0\n")
    (tpl / "POTCAR").write_text("dummy\n")
    setup, _ = generate("vasp", str(tpl), str(tmp_path / "el"))
    c, C = random_stiffness(setup.form)
    clean = synthetic_stresses(setup, C)
    rng = np.random.default_rng(0)
    for st in setup.states:
        for i in range(40):
            d = tmp_path / "el" / st["id"] / f"conf-{i:03d}"
            d.mkdir()
            _fake_vasprun(d / "vasprun.xml", clean[st["id"]] + rng.normal(scale=0.02, size=6), np.zeros((4, 3)))
    r = collect(str(tmp_path / "el"), samples="conf-*")
    assert r.errors is not None and np.all(r.errors > 0)
    assert np.all(np.abs(r.constants - c) < 5 * r.errors)


def _afm_fcc():
    """fcc Ni, conventional cell, layered (AFM-I) moments along z: P4/mmm instead of Fm-3m."""
    at = bulk("Ni", "fcc", a=3.52, cubic=True)[[0, 3, 1, 2]]        # z = 0 layer first, then z = 1/2
    return at


def test_magnetic_order_lowers_symmetry(tmp_path):
    """Layered AFM moments on fcc Ni give P4/mmm (tetragonal-II) instead of Fm-3m: from initial moments, from a
    spin-polarized VASP MAGMOM, from two POSCAR species blocks and from two pw.x species labels; MAGMOM with
    ISPIN = 1 is ignored; the moments survive to_dict / from_dict.
    """
    at = _afm_fcc()
    assert ElasticSetup(at).symmetry.elastic_type == "cubic"
    at.set_initial_magnetic_moments([1.0, 1.0, -1.0, -1.0])
    setup = ElasticSetup(at)
    assert setup.symmetry.number == 123 and setup.symmetry.elastic_type == "tetragonal-II"
    back = ElasticSetup.from_dict(json.loads(json.dumps(setup.to_dict())))      # as collect() restores it
    assert back.symmetry.number == 123 and np.allclose(back.patterns, setup.patterns)
    # VASP: MAGMOM of a spin-polarized INCAR, or the same element in two species blocks
    for name, incar, species in [("magmom", "ISPIN = 2\nMAGMOM = 2*1.0 2*-1.0\n", "Ni"),
                                 ("blocks", "ISPIN = 1\n", "Ni Ni"),
                                 ("nonmagnetic", "ISPIN = 1\nMAGMOM = 2*1.0 2*-1.0\n", "Ni")]:
        tpl = tmp_path / name
        tpl.mkdir()
        cell = "\n".join(" ".join(f"{x:.10f}" for x in v) for v in at.cell[:])
        pos = "\n".join(" ".join(f"{x:.10f}" for x in p) for p in at.get_scaled_positions())
        counts = "2 2" if species == "Ni Ni" else "4"
        (tpl / "POSCAR").write_text(f"afm\n1.0\n{cell}\n{species}\n{counts}\nDirect\n{pos}\n")
        (tpl / "INCAR").write_text("ENCUT = 400\n" + incar)
        (tpl / "KPOINTS").write_text("a\n0\nGamma\n8 8 8\n0 0 0\n")
        (tpl / "POTCAR").write_text("dummy\n")
        setup, _ = generate("vasp", str(tpl), str(tmp_path / f"el_{name}"))
        assert setup.symmetry.elastic_type == ("cubic" if name == "nonmagnetic" else "tetragonal-II"), name
    # pw.x: two species labels of one element
    tpl = tmp_path / "pw.in"
    rows = "\n".join(f"Ni{1 if n < 2 else 2} " + " ".join(f"{x:.10f}" for x in p)
                     for n, p in enumerate(at.get_scaled_positions()))
    tpl.write_text("&control\n calculation='scf'\n/\n&system\n ibrav = 0, nat = 4, ntyp = 2, ecutwfc = 40\n/\n"
                   "&electrons\n/\nATOMIC_SPECIES\nNi1 58.69 Ni.UPF\nNi2 58.69 Ni.UPF\n"
                   f"CELL_PARAMETERS angstrom\n{cell}\nATOMIC_POSITIONS crystal\n{rows}\n"
                   "K_POINTS automatic\n8 8 8 0 0 0\n")
    assert qe.read_structure(str(tpl)).get_chemical_symbols() == ["Ni"] * 4
    setup, _ = generate("qe", str(tpl), str(tmp_path / "el_qe"))
    assert setup.symmetry.elastic_type == "tetragonal-II"


def test_abinit_interface(tmp_path):
    """generate abinit writes POSCARs with the element at the end of each coordinate line (as ABINIT's POSCAR
    reader needs) and copies the input; read_result returns the stress of a real ABINIT output in GPa with
    -tr(sigma)/3 equal to the pressure ABINIT printed."""
    from minelastic.io import abinit
    tpl = tmp_path / "ref"
    tpl.mkdir()
    atoms = structure("hexagonal")
    abinit.write_poscar(tpl / "POSCAR", atoms, "reference")
    (tpl / "run.abi").write_text('structure "poscar:POSCAR"\npseudos "O.psp8, Si.psp8"\necut 40\n')
    setup, notes = generate("abinit", str(tpl), str(tmp_path / "el"), relax_ions=True)
    assert notes == []
    lines = (tmp_path / "el" / "s1_p1" / "POSCAR").read_text().splitlines()
    assert lines[8].split()[-1] == atoms.get_chemical_symbols()[0]
    back = read(tmp_path / "el" / "s1_p1" / "POSCAR", format="vasp")
    assert np.allclose(back.cell[:], setup.strained_atoms("s1_p1").cell[:], atol=1e-12)
    assert "ionmov 2" in (tmp_path / "el" / "s1_p1" / "run.abi").read_text()
    out = os.path.join(os.path.dirname(__file__), "data_abinit")
    s, f = abinit.read_result(out, "aln_ecut10.abo")
    printed = next(l for l in open(os.path.join(out, "aln_ecut10.abo")) if "Pressure=" in l)
    assert abs(-s[:3].mean() - float(printed.split("Pressure=")[1].split()[0])) < 1e-3 and f.shape == (4, 3)


def test_aims_reader(tmp_path):
    """FHI-aims output shipped with ASE (numerical stress): GPa, and -tr(sigma)/3 equals the printed pressure."""
    import ase
    from minelastic.io import aims
    src = os.path.join(os.path.dirname(ase.__file__), "test", "testdata", "aims", "numerical_stress.out")
    if not os.path.exists(src):
        pytest.skip("ASE test data not installed")
    (tmp_path / "aims.out").write_text(open(src).read())
    s, _ = aims.read_result(str(tmp_path))
    printed = next(l for l in open(src) if "Pressure:" in l)
    assert abs(-s[:3].mean() - float(printed.split()[2]) * 160.21766208) < 1e-6


def test_vasp_md_workflow(tmp_path):
    """Finite-T protocol: fixed-cell MD in every strained cell, stress averaged over the steps after --skip."""
    tpl = tmp_path / "ref"
    tpl.mkdir()
    write(tpl / "POSCAR", bulk("Al", "fcc", a=4.04, cubic=True), format="vasp", direct=True)
    (tpl / "INCAR").write_text("ENCUT = 300\nIBRION = 0\nNSW = 5000\nPOTIM = 2\nMDALGO = 3\nISIF = 3\nTEBEG = 600\n")
    (tpl / "KPOINTS").write_text("a\n0\nGamma\n4 4 4\n0 0 0\n")
    (tpl / "POTCAR").write_text("dummy\n")
    with pytest.raises(ValueError):
        generate("vasp", str(tpl), str(tmp_path / "bad"), md=True, relax_ions=True)
    setup, _ = generate("vasp", str(tpl), str(tmp_path / "el"), md=True)
    tags = dict(vasp.read_incar(tmp_path / "el" / "s1_p1" / "INCAR"))
    assert tags["IBRION"] == "0" and tags["NSW"] == "5000" and tags["POTIM"] == "2" and tags["ISIF"] == "2"
    assert tags["MDALGO"] == "3" and tags["TEBEG"] == "600"
    c, C = random_stiffness(setup.form)
    clean = synthetic_stresses(setup, C)
    rng = np.random.default_rng(1)
    for st in setup.states:
        # 100 equilibration steps far off, then 400 noisy steps around the true stress
        steps = [clean[st["id"]] + 5.0] * 100 + [clean[st["id"]] + rng.normal(scale=0.02, size=6) for _ in range(400)]
        blocks = []
        for s in steps:
            t = -10 * np.array([[s[0], s[5], s[4]], [s[5], s[1], s[3]], [s[4], s[3], s[2]]])
            rows = "\n".join(f"   <v> {a:16.8f} {b:16.8f} {c_:16.8f} </v>" for a, b, c_ in t)
            blocks.append(f' <calculation>\n  <varray name="stress" >\n{rows}\n  </varray>\n </calculation>')
        (tmp_path / "el" / st["id"] / "vasprun.xml").write_text(
            '<?xml version="1.0" encoding="ISO-8859-1"?>\n<modeling>\n' + "\n".join(blocks) + "\n</modeling>\n")
    with pytest.raises(ValueError):
        collect(str(tmp_path / "el"))                       # MD runs need md=True or samples=
    r = collect(str(tmp_path / "el"), md=True, skip=100)
    assert r.errors is not None and np.all(np.abs(r.constants - c) < 5 * r.errors)


def test_qe_md_inputs(tmp_path):
    """generate qe with md: an scf template is refused; an md template keeps calculation = 'md' and &IONS."""
    tpl = tmp_path / "pw.in"
    tpl.write_text(QE_TEMPLATE)
    with pytest.raises(ValueError):
        generate("qe", str(tpl), str(tmp_path / "bad"), md=True)        # template is an scf run
    tpl.write_text(QE_TEMPLATE.replace("calculation='scf'", "calculation='md', nstep=2000, dt=40")
                   + "&ions\n  ion_temperature='andersen', tempw=600\n/\n")
    generate("qe", str(tpl), str(tmp_path / "el"), md=True)
    txt = (tmp_path / "el" / "s1_p1" / "pw.in").read_text()
    assert "calculation = 'md'" in txt and "tstress = .true." in txt and "tempw" in txt


def _random_piezo(rotations, seed=0):
    """Random 3x6 (Voigt, tensor shear components) piezoelectric tensor with the point-group symmetry."""
    T = np.random.default_rng(seed).normal(size=(3, 3, 3))
    T = 0.5 * (T + T.transpose(0, 2, 1))
    T = np.mean([np.einsum("ia,jb,kc,abc->ijk", R, R, R, T) for R in rotations], axis=0)
    return np.array([[T[i, a, b] for a, b in ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))] for i in range(3)])


@pytest.mark.parametrize("rotate", [False, True])
def test_piezo_from_polarization(rotate):
    """Proper piezo tensor from Berry-phase dipoles of the minimal strained cells, with branch jumps."""
    from minelastic.piezo import E_PER_A2_IN_C_PER_M2
    aln = crystal(["Al", "N"], basis=[(1 / 3, 2 / 3, 0), (1 / 3, 2 / 3, 0.3818)], spacegroup=186,
                  cellpar=[3.13, 3.13, 5.02, 90, 90, 120])
    atoms = rotate_atoms(aln, random_rotation(11)) if rotate else aln
    setup = ElasticSetup(atoms)
    assert len(setup.patterns) == 2
    e = _random_piezo(setup.symmetry.rotations)                    # in the working frame
    assert abs(e[2, 2]) > 0.1 and np.allclose(e[2, 3:], 0, atol=1e-10)   # 6mm form: e33 != 0, e3(4,5,6) = 0
    c, C = random_stiffness(setup.form)
    A0, V0 = np.asarray(setup.atoms.cell[:]), setup.atoms.get_volume()
    phi0 = np.array([0.37, -0.81, 0.49])
    rng = np.random.default_rng(2)
    dip = {}
    for st in setup.states:
        P = e @ np.array(st["strain"])                              # C/m^2 (proper, linear)
        phi = phi0 + np.linalg.solve(A0.T, P * V0 / E_PER_A2_IN_C_PER_M2) + rng.integers(-2, 3, 3)   # branch jumps
        dip[st["id"]] = np.asarray(setup.strained_atoms(st).cell[:]).T @ phi
    r = setup.compute(synthetic_stresses(setup, C), polarization=dip)
    assert np.allclose(r.piezo, e, atol=1e-8), np.abs(r.piezo - e).max()
