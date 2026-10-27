# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Edge cases: magnetic symmetry, kinetic term, frames, Berry branches, parsers, input checks, setup records."""
import json

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.io import write
from ase.spacegroup import crystal

import minelastic
from minelastic import ElasticSetup
from minelastic.io import abinit, aims, qe, vasp
from minelastic.piezo import align_branch
from minelastic.symmetry import analyze, rotate_atoms
from minelastic.voigt import rotate_stiffness
from minelastic.workflow import generate
from test_minelastic import random_rotation, structure


def test_noncollinear_moment_lowers_symmetry():
    """A moment along z (axial vector) makes a simple cubic cell tetragonal (P4/mmm), along [111] trigonal;
    rotate_atoms rotates vector moments with the cell, so the symmetry type is kept."""
    sc = bulk("Fe", "sc", a=2.8)
    sc.set_initial_magnetic_moments([[0, 0, 1.0]])
    assert analyze(sc).number == 123 and analyze(sc).elastic_type == "tetragonal-II"
    sc.set_initial_magnetic_moments([[1.0, 1.0, 1.0]])
    assert analyze(sc).elastic_type == "trigonal-II"
    assert analyze(rotate_atoms(sc, random_rotation(2))).elastic_type == "trigonal-II"


def test_antiferromagnet_keeps_time_reversed_operations():
    """MnF2: the two Mn are exchanged by operations combined with time reversal; the elastic tensor keeps the
    tetragonal form of P4_2/mnm (the moments as plain species labels would give an orthorhombic group)."""
    at = crystal(["Mn", "F"], basis=[(0, 0, 0), (0.305, 0.305, 0)], spacegroup=136,
                 cellpar=[4.87, 4.87, 3.31, 90, 90, 90])
    at.set_initial_magnetic_moments([5.0 if n == 0 else -5.0 if n == 1 else 0.0 for n in range(len(at))])
    s = analyze(at)
    assert s.number == 136 and s.elastic_type == "tetragonal-II" and len(s.rotations) == 16


def test_vasp_lsorbit_magmom(tmp_path):
    """LSORBIT alone makes the run noncollinear: MAGMOM with 3N values is read as vectors."""
    tpl = tmp_path / "ref"
    tpl.mkdir()
    write(tpl / "POSCAR", bulk("Fe", "sc", a=2.8), format="vasp", direct=True)
    (tpl / "INCAR").write_text("LSORBIT = .TRUE.\nMAGMOM = 0 0 2\n")
    at = vasp.read_structure(str(tpl))
    assert at.get_initial_magnetic_moments().shape == (1, 3) and analyze(at).elastic_type == "tetragonal-II"


def test_magnetic_standard_frame():
    """The standard orientation of a magnetic cell is that of its magnetic group: NiO with AFM-II order in the
    rhombohedral cell with the 3-fold axis along [111] is rotated to the textbook trigonal form (axis along z,
    c11 c12 c13 c14 c33 c44, 2 patterns); a randomly rotated hexagonal cell with moments along c is rotated
    back to the hexagonal form of the unrotated cell."""
    from ase import Atoms
    a = 4.17
    nio = Atoms("NiNiOO", cell=a * np.array([[1, .5, .5], [.5, 1, .5], [.5, .5, 1]]), pbc=True,
                positions=a * np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, 0], [1, .5, 0]]))
    nio.set_initial_magnetic_moments([2.0, -2.0, 0, 0])
    setup = ElasticSetup(nio)
    assert setup.symmetry.number == 166 and setup.n_min == 2 and len(setup.patterns) == 2
    assert set(setup.form.labels) == {"c11", "c12", "c13", "c14", "c33", "c44"}
    assert np.allclose(np.abs(setup.rotation @ np.ones(3) / np.sqrt(3)), [0, 0, 1])
    hexa = structure("hexagonal")
    hexa.set_initial_magnetic_moments([[0, 0, 1.0]] * len(hexa))
    rotated = ElasticSetup(rotate_atoms(hexa, random_rotation(5)))
    assert rotated.form.labels == ElasticSetup(hexa).form.labels and len(rotated.patterns) == 2


def test_abinit_spinat(tmp_path):
    """ABINIT: spinat gives the moments, as z components for nsppol 2 or nspden 2 (layered AFM fcc Ni: P4/mmm),
    as vectors for nspden 4 (moment along z: tetragonal, along [111]: trigonal), and is ignored without spin
    polarization; n*x repeats x. A cell rotated to its standard orientation is refused with noncollinear
    spinat."""
    from test_minelastic import _afm_fcc
    tpl = tmp_path / "afm"
    tpl.mkdir()
    abinit.write_poscar(tpl / "POSCAR", _afm_fcc(), "afm")
    for spin, number in [("nsppol 2\nspinat 0 0 1.0d0  0 0 1  0 0 -1  0 0 -1\n", 123),
                         ("nspden 2 # nsppol 1\nspinat 2*0 1 0 0 1 2*0 -1 0 0 -1\n", 123),
                         ("spinat 2*0 1 0 0 1 2*0 -1 0 0 -1\n", 225)]:
        (tpl / "run.abi").write_text('structure "poscar:POSCAR"\n' + spin)
        assert analyze(abinit.read_structure(str(tpl))).number == number, spin
    (tpl / "run.abi").write_text("nsppol 2\nspinat1 0 0 1 0 0 1 0 0 -1 0 0 -1\n")
    with pytest.raises(ValueError, match="dataset-specific"):
        abinit.read_structure(str(tpl))
    (tpl / "run.abi").write_text("nsppol 2\nspinat 0 0 1\n")
    with pytest.raises(ValueError, match="expected 12 numbers"):
        abinit.read_structure(str(tpl))
    tpl = tmp_path / "ncl"
    tpl.mkdir()
    abinit.write_poscar(tpl / "POSCAR", bulk("Fe", "sc", a=2.8), "fe")
    for spinat, kind in [("0 0 1", "tetragonal-II"), ("1 1 1", "trigonal-II")]:
        (tpl / "run.abi").write_text(f"nspinor 2\nnspden 4\nspinat {spinat}\n")
        assert analyze(abinit.read_structure(str(tpl))).elastic_type == kind
    hexa = structure("hexagonal")
    hexa.set_initial_magnetic_moments([[0, 0, 1.0]] * len(hexa))          # along c: hexagonal family group
    atoms = rotate_atoms(hexa, random_rotation(5))
    tpl = tmp_path / "rot"
    tpl.mkdir()
    abinit.write_poscar(tpl / "POSCAR", atoms, "rotated")
    m = " ".join(f"{x:.12f}" for x in atoms.get_initial_magnetic_moments().ravel())
    (tpl / "run.abi").write_text(f'structure "poscar:POSCAR"\nnspinor 2\nnspden 4\nspinat {m}\n')
    with pytest.raises(ValueError, match="noncollinear spinat"):
        generate("abinit", str(tpl), str(tmp_path / "el"))
    setup, _ = generate("abinit", str(tpl), str(tmp_path / "el2"), frame="input")
    assert setup.symmetry.elastic_type == "hexagonal"


def test_qe_noncollinear_moments(tmp_path):
    """pw.x with noncolin: starting_magnetization with angle1/angle2 gives vector moments (along z: tetragonal,
    along [111]: trigonal); a collinear ferromagnet keeps the cubic symmetry."""
    def pw(system):
        """Elastic type of simple cubic Fe from a pw.x input with the given &system settings."""
        p = tmp_path / "pw.in"
        p.write_text("&control\n/\n&system\n ibrav = 1, celldm(1) = 5.3, nat = 1, ntyp = 1, ecutwfc = 40,\n"
                     f" {system}\n/\n&electrons\n/\nATOMIC_SPECIES\nFe 55.85 Fe.UPF\n"
                     "ATOMIC_POSITIONS crystal\nFe 0 0 0\nK_POINTS automatic\n4 4 4 0 0 0\n")
        return analyze(qe.read_structure(str(p))).elastic_type
    assert pw("nspin = 2, starting_magnetization(1) = 0.5") == "cubic"
    assert pw("noncolin = .true., lspinorb = .true., starting_magnetization(1) = 0.5") == "tetragonal-II"
    assert pw("noncolin = .true., starting_magnetization(1) = 0.5, angle1(1) = 54.7356103, angle2(1) = 45") \
        == "trigonal-II"
    assert pw("noncolin = .true., lspinorb = .true.") == "cubic"


def test_kinetic_term_only_for_cauchy_stress():
    """Zero virial response: d sigma/d eps = n kB T for the Cauchy stress, and 0 for J sigma (the kinetic part
    of J sigma, -N kB T / V0, does not depend on the strain)."""
    setup = ElasticSetup(bulk("Al", "fcc", a=4.04, cubic=True))
    zero = {st["id"]: np.zeros(6) for st in setup.states}
    nkT = 4 / setup.atoms.get_volume() * 8.617333262e-5 * 300 * 160.21766208
    c = setup.compute(zero, kinetic_temperature=300).C
    k = setup.compute(zero, kinetic_temperature=300, stress_measure="kirchhoff").C
    assert np.allclose(c[:3, :3], nkT) and np.allclose(c[3:, 3:], 0)
    assert np.allclose(k, 0)


def test_branch_alignment_along_the_stencil():
    """Phases changing by 0.2 per point stay on one branch when unwrapped along increasing strain."""
    phis = np.array([[-0.4], [-0.2], [0.2], [0.4]])
    assert np.allclose(align_branch(phis, 1.0, x=[-2, -1, 1, 2]), phis)
    assert np.allclose(align_branch([[0.49], [-0.49]], 1.0, x=[-1, 1]), [[0.49], [0.51]])


def test_one_atom_cell_ionic_route():
    """One atom per cell: no optical modes, C_ion = 0."""
    r = minelastic.calculate(bulk("Cu", "fcc", a=3.59), EMT(), ionic=True)
    assert np.allclose(r.C_ionic, 0) and np.allclose(r.C_relaxed, r.C)


def test_frame_rotation_refuses_direction_dependent_inputs(tmp_path):
    """A cell rotated to the standard orientation is refused with LSORBIT, SAXIS or Cartesian k points."""
    atoms = rotate_atoms(structure("hexagonal"), random_rotation(5))
    for incar, kpts in [("LSORBIT = .TRUE.\n", "a\n0\nGamma\n4 4 4\n"), ("SAXIS = 0 0 1\n", "a\n0\nGamma\n4 4 4\n"),
                        ("ENCUT = 400\n", "a\n1\nCartesian\n0 0 0 1\n")]:
        tpl = tmp_path / f"t{len(list(tmp_path.iterdir()))}"
        tpl.mkdir()
        write(tpl / "POSCAR", atoms, format="vasp", direct=True)
        (tpl / "INCAR").write_text(incar)
        (tpl / "KPOINTS").write_text(kpts)
        (tpl / "POTCAR").write_text("dummy\n")
        with pytest.raises(ValueError, match="frame='input'"):
            generate("vasp", str(tpl), str(tpl / "el"))
    generate("vasp", str(tpl), str(tpl / "el2"), frame="input")          # allowed in the input orientation


def test_qe_parser():
    """Single-line namelists end at '/' outside quotes; species labels are case-insensitive."""
    nl, cards = qe.parse_input("&CONTROL calculation='scf', pseudo_dir='/a/b/' /\n&SYSTEM ibrav=1, nat=1, ntyp=1 /\n"
                               "ATOMIC_SPECIES\nFE1 55.8 Fe.UPF\nK_POINTS gamma\n")
    assert dict(nl["CONTROL"])["pseudo_dir"] == "'/a/b/'" and dict(nl["SYSTEM"])["ntyp"] == "1"
    assert [c[0] for c in cards] == ["ATOMIC_SPECIES", "K_POINTS gamma"]
    assert [qe._element(x) for x in ("Fe1", "FE1", "fe_up", "C1", "Ca2")] == ["Fe", "Fe", "Fe", "C", "Ca"]


def test_vasp_fixed_width_and_blocks(tmp_path):
    """Fields that run together are split by their decimals; the piezoelectric block is found without a
    dashed separator; unequal numbers of stress and force blocks are refused."""
    assert vasp._fixed("  in kB  -123456.78901-234567.89012     0.93457", 5) == [-123456.78901, -234567.89012, 0.93457]
    out = tmp_path / "OUTCAR"
    out.write_text(" PIEZOELECTRIC TENSOR (including local field effects)  for field in x, y, z        (C/m^2)\n"
                   "          XX          YY          ZZ          XY          YZ          ZX\n"
                   "  x     0.00000     0.00000     0.00000     0.00000     0.00000     0.34039\n"
                   "  y     0.00000     0.00000     0.00000     0.00000     0.34039     0.00000\n"
                   "  z     0.27039     0.27039    -0.42396     0.00000     0.00000     0.00000\n")
    e = vasp.read_clamped_piezo(str(out))
    assert e[0, 4] == 0.34039 and e[2, 2] == -0.42396
    out.write_text("   number of dos      NEDOS =    301   number of ions     NIONS =      1\n"
                   "  in kB       1.00000     1.00000     1.00000     0.00000     0.00000     0.00000\n"
                   " POSITION                                       TOTAL-FORCE (eV/Angst)\n -----\n"
                   "      0.00000      0.00000      0.00000         0.000000      0.000000      0.000000\n"
                   "  in kB       2.00000     2.00000     2.00000     0.00000     0.00000     0.00000\n")
    with pytest.raises(ValueError, match="incomplete"):
        vasp.read_outcar(str(out))


def test_aims_and_abinit_relaxation_checks(tmp_path):
    """Clamped-ion runs refuse templates that relax the ions; commented directives are ignored."""
    tpl = tmp_path / "aims"
    tpl.mkdir()
    write(tpl / "geometry.in", bulk("Cu", "fcc", a=3.6), format="aims")
    (tpl / "control.in").write_text("compute_analytical_stress .true.\n# relax_geometry bfgs 1E-3\n")
    generate("aims", str(tpl), str(tmp_path / "a1"))
    assert not any(l.startswith("relax_geometry") for l in (tmp_path / "a1" / "s1_p1" / "control.in").read_text().splitlines())
    generate("aims", str(tpl), str(tmp_path / "a2"), relax_ions=True)
    assert (tmp_path / "a2" / "s1_p1" / "control.in").read_text().count("\nrelax_geometry bfgs 1E-3") == 1
    (tpl / "control.in").write_text("compute_analytical_stress .true.\nrelax_geometry bfgs 1E-3\n")
    with pytest.raises(ValueError, match="relaxes the ions"):
        generate("aims", str(tpl), str(tmp_path / "a3"))
    (tpl / "control.in").write_text("compute_analytical_stress .false.\n")
    with pytest.raises(ValueError, match="stress"):
        generate("aims", str(tpl), str(tmp_path / "a4"))
    tpl = tmp_path / "abi"
    tpl.mkdir()
    abinit.write_poscar(tpl / "POSCAR", bulk("Cu", "fcc", a=3.6), "Cu")
    (tpl / "run.abi").write_text('structure "poscar:POSCAR"\nionmov 2\n')
    with pytest.raises(ValueError, match="relaxes the ions"):
        generate("abinit", str(tpl), str(tmp_path / "b1"))
    (tpl / "run.abi").write_text('structure "poscar:POSCAR"\nionmov 0\n')
    with pytest.raises(ValueError, match="ionmov 0"):
        generate("abinit", str(tpl), str(tmp_path / "b2"), relax_ions=True)


def test_setup_record_is_validated():
    """from_dict refuses a record whose stored strains differ from those its delta and patterns give."""
    d = json.loads(json.dumps(ElasticSetup(bulk("Al", "fcc", a=4.04, cubic=True)).to_dict()))
    ElasticSetup.from_dict(d)
    d["delta"] = 0.02
    with pytest.raises(ValueError, match="differ"):
        ElasticSetup.from_dict(d)


def test_input_frame_outputs_and_singular_summary():
    """Rotated cell with the ionic route: every frame-dependent tensor is also given in the input frame; a
    singular C is reported as unstable instead of failing."""
    at = rotate_atoms(bulk("Cu", "hcp", a=2.55, c=4.16), random_rotation(4))
    r = minelastic.calculate(at, EMT(), ionic=True)
    d = r.to_dict()
    Rt = r.setup.rotation.T
    assert np.allclose(d["C_relaxed_input_frame"], rotate_stiffness(r.C_relaxed, Rt))
    assert np.allclose(d["C_input_frame"], rotate_stiffness(r.C, Rt))
    setup = ElasticSetup(bulk("Al", "fcc", a=4.04, cubic=True))
    res = setup.compute({st["id"]: np.zeros(6) for st in setup.states})
    assert "Born stable: False" in res.summary()


def test_lstsq_errors_with_samples():
    """Least squares with sample errors: inverse-variance weights, errors consistent with the scatter."""
    from test_minelastic import random_stiffness, synthetic_stresses
    setup = ElasticSetup(structure("hexagonal"))
    c, C = random_stiffness(setup.form)
    rng = np.random.default_rng(5)
    noisy = {k: v + rng.normal(scale=0.05, size=(400, 6)) for k, v in synthetic_stresses(setup, C).items()}
    r = setup.compute(noisy, method="lstsq")
    assert np.all(np.abs(r.constants - c) < 5 * r.errors)


def test_errors_include_correlations():
    """Samples with sigma_yy = sigma_zz exactly: c12 from both components is no better than from one, so its
    error equals that of a single relation (independent errors would make it smaller by sqrt 2)."""
    from test_minelastic import random_stiffness, synthetic_stresses
    setup = ElasticSetup(bulk("Al", "fcc", a=4.04, cubic=True))
    c, C = random_stiffness(setup.form)
    rng = np.random.default_rng(9)
    noisy = {}
    for k, v in synthetic_stresses(setup, C).items():
        x = v + rng.normal(scale=0.05, size=(400, 6))
        x[:, 2] = x[:, 1] - v[1] + v[2]                      # sigma_zz noise = sigma_yy noise
        noisy[k] = x
    both = setup.compute(noisy).errors[1]
    one = setup.compute(noisy, method="designated").errors[1]
    assert abs(both / one - 1) < 1e-8
    gls = setup.compute(noisy, method="lstsq")
    assert np.all(np.abs(gls.constants - c) < 5 * gls.errors)


def test_options_are_validated():
    """Unknown schemes, frames, stencils and internal_strain_from values, duplicate or zero multipliers and a
    non-positive step are refused."""
    at = bulk("Al", "fcc", a=4.04, cubic=True)
    for kw in (dict(scheme="asesss"), dict(frame="std"), dict(stencil="five-point"), dict(stencil=(1, 1)),
               dict(stencil=(-1,)), dict(stencil=(0, 1)), dict(delta=0.0), dict(ionic=True, displacement=-0.01)):
        with pytest.raises(ValueError):
            ElasticSetup(at, **kw)
    setup = ElasticSetup(at)
    with pytest.raises(ValueError):
        setup.compute({st["id"]: np.zeros(6) for st in setup.states}, internal_strain_from="strain")


def test_internal_strain_in_the_input_frame():
    """A rotated hcp cell, strained in the standard frame: Lambda rotated back to the input frame equals Lambda
    computed directly in the input frame."""
    at = rotate_atoms(bulk("Cu", "hcp", a=2.55, c=4.16), random_rotation(6))
    std = minelastic.calculate(at, EMT(), delta=0.005)
    inp = minelastic.calculate(at, EMT(), delta=0.005, frame="input", method="lstsq")
    assert not np.allclose(std.setup.rotation, np.eye(3))
    assert np.abs(np.array(std.to_dict()["internal_strain_input_frame"]) - inp.internal_strain).max() < 1e-3


def test_more_input_checks(tmp_path):
    """pw.x: comments only outside quotes, incomplete force blocks refused, Cartesian constraints refused for a
    rotated cell; VASP: Cartesian selective dynamics and a non-default SAXIS refused; ABINIT: two species of one
    element and conflicting or dataset-specific ionmov refused."""
    nl, _ = qe.parse_input("&CONTROL pseudo_dir='/a!b#c' ! comment\n/\n")
    assert dict(nl["CONTROL"])["pseudo_dir"] == "'/a!b#c'"
    out = tmp_path / "pw.out"
    out.write_text("     number of atoms/cell      =            2\n     Forces acting on atoms (cartesian axes, Ry/au):\n\n"
                   "     atom    1 type  1   force =     0.0 0.0 0.0\n\n     total   stress  (Ry/bohr**3)   (kbar)  P= 0\n"
                   "  0.0 0.0 0.0  0 0 0\n  0.0 0.0 0.0  0 0 0\n  0.0 0.0 0.0  0 0 0\n")
    with pytest.raises(ValueError, match="incomplete force block"):
        qe.read_output(str(out))
    from test_minelastic import QE_TEMPLATE
    tpl = tmp_path / "pw.in"
    tpl.write_text(QE_TEMPLATE)
    generate("qe", str(tpl), str(tmp_path / "q0"))                    # standard orientation: if_pos accepted
    at = rotate_atoms(qe.read_structure(str(tpl)), random_rotation(8))   # same crystal, random orientation
    cell = "\n".join(" ".join(f"{x:.10f}" for x in v) for v in at.cell[:])
    pos = "\n".join(("Al1 " if n < 2 else "N ") + " ".join(f"{x:.10f}" for x in f) + (" 0 0 1" if n == 0 else "")
                    for n, f in enumerate(at.get_scaled_positions()))
    head = QE_TEMPLATE.split("ATOMIC_SPECIES")[0].replace("ibrav = 4, celldm(1) = 5.911804, celldm(3) = 1.603472,",
                                                           "ibrav = 0,")
    tpl.write_text(head + "ATOMIC_SPECIES\nAl1 26.982 Al.UPF\nN  14.007 N.UPF\nCELL_PARAMETERS angstrom\n" + cell
                   + "\nATOMIC_POSITIONS crystal\n" + pos + "\nK_POINTS automatic\n8 8 5 0 0 0\n")
    with pytest.raises(ValueError, match="if_pos"):
        generate("qe", str(tpl), str(tmp_path / "q1"))
    tv = tmp_path / "v"
    tv.mkdir()
    (tv / "POSCAR").write_text("Cu\n3.6\n0 0.5 0.5\n0.5 0 0.5\n0.5 0.5 0\nCu\n1\nSelective dynamics\nCartesian\n"
                               "0 0 0 F T T\n")
    (tv / "INCAR").write_text("ENCUT = 400\n")
    (tv / "KPOINTS").write_text("a\n0\nGamma\n8 8 8\n")
    (tv / "POTCAR").write_text("dummy\n")
    with pytest.raises(ValueError, match="Cartesian"):
        generate("vasp", str(tv), str(tmp_path / "v1"))
    (tv / "POSCAR").write_text("Cu\n3.6\n0 0.5 0.5\n0.5 0 0.5\n0.5 0.5 0\nCu\n1\nDirect\n0 0 0\n")
    (tv / "INCAR").write_text("LNONCOLLINEAR = .TRUE.\nMAGMOM = 0 0 1\nSAXIS = 1 0 0\n")
    with pytest.raises(ValueError, match="SAXIS"):
        vasp.read_structure(str(tv))
    ta = tmp_path / "abi"
    ta.mkdir()
    (ta / "POSCAR").write_text("Fe2\n1.0\n2.8 0 0\n0 2.8 0\n0 0 2.8\nFe Fe\n1 1\nDirect\n0 0 0 Fe\n0.5 0.5 0.5 Fe\n")
    (ta / "run.abi").write_text('structure "poscar:POSCAR"\n')
    with pytest.raises(ValueError, match="two species of one element"):
        generate("abinit", str(ta), str(tmp_path / "a1"))
    abinit.write_poscar(ta / "POSCAR", bulk("Cu", "fcc", a=3.6), "Cu")
    for text in ('structure "poscar:POSCAR"\nionmov 0\nionmov 2\n', 'structure "poscar:POSCAR"\nionmov1 2\n'):
        (ta / "run.abi").write_text(text)
        with pytest.raises(ValueError, match="ionmov"):
            generate("abinit", str(ta), str(tmp_path / "a2"))


def test_default_stencils(tmp_path):
    """Four-point stencil (-2, -1, 1, 2) by default, also for generate; central (+-1) for MD inputs."""
    al = bulk("Al", "fcc", a=4.04, cubic=True)
    assert ElasticSetup(al).multipliers == (-2, -1, 1, 2) and len(ElasticSetup(al).states) == 4
    tpl = tmp_path / "ref"
    tpl.mkdir()
    write(tpl / "POSCAR", al, format="vasp", direct=True)
    (tpl / "KPOINTS").write_text("a\n0\nGamma\n4 4 4\n0 0 0\n")
    (tpl / "POTCAR").write_text("dummy\n")
    (tpl / "INCAR").write_text("ENCUT = 400\nIBRION = 0\nNSW = 10\nPOTIM = 2\n")
    assert generate("vasp", str(tpl), str(tmp_path / "md"), md=True)[0].multipliers == (-1, 1)
    (tpl / "INCAR").write_text("ENCUT = 400\n")
    assert generate("vasp", str(tpl), str(tmp_path / "st"))[0].multipliers == (-2, -1, 1, 2)


def test_vasp_incar_report_and_fix(tmp_path):
    """generate vasp lists the tags it sets and removes; EDIFF, ENCUT (1.3 ENMAX of the POTCAR), PREC and LREAL
    are only reported unless fix_incar applies all of them; MD runs get the looser EDIFF suggestion."""
    tpl = tmp_path / "ref"
    tpl.mkdir()
    write(tpl / "POSCAR", bulk("Al", "fcc", a=4.04), format="vasp", direct=True)
    (tpl / "KPOINTS").write_text("a\n0\nGamma\n4 4 4\n0 0 0\n")
    (tpl / "POTCAR").write_text("  PAW_PBE Al 04Jan2001\n   ENMAX  =  240.300; ENMIN  =  180.225 eV\n")
    (tpl / "INCAR").write_text("PREC = Normal\nIBRION = 6\nISIF = 3\nNFREE = 2\nPOTIM = 0.015\nLEPSILON = .TRUE.\n"
                               "LREAL = Auto\n")
    _, notes = generate("vasp", str(tpl), str(tmp_path / "a"), stencil="central")
    text = " ".join(notes)
    assert "set IBRION = -1 (template: 6), NSW = 0, ISIF = 2 (template: 3)" in text
    assert "removed NFREE, POTIM, LEPSILON" in text
    assert "EDIFF default 1E-4 -> 1E-8" in text and "ENCUT default ENMAX = 240.3 -> 320" in text
    assert "PREC Normal -> Accurate" in text and "LREAL Auto -> .FALSE." in text and "--fix-incar" in text
    assert dict(vasp.read_incar(tmp_path / "a" / "s1_p1" / "INCAR"))["PREC"] == "Normal"   # kept without the fix
    _, notes = generate("vasp", str(tpl), str(tmp_path / "b"), stencil="central", fix_incar=True)
    incar = dict(vasp.read_incar(tmp_path / "b" / "s1_p1" / "INCAR"))
    assert {k: incar[k] for k in ("EDIFF", "ENCUT", "PREC", "LREAL")} == \
        {"EDIFF": "1E-8", "ENCUT": "320", "PREC": "Accurate", "LREAL": ".FALSE."}
    assert any(n.startswith("INCAR: applied the recommended settings") for n in notes)
    (tpl / "INCAR").write_text("PREC = Accurate\nENCUT = 400\nEDIFF = 1E-8\n")
    assert not any("recommended" in n for n in generate("vasp", str(tpl), str(tmp_path / "c"))[1])
    (tpl / "INCAR").write_text("PREC = Accurate\nENCUT = 400\nEDIFF = 1E-5\nIBRION = 0\nNSW = 10\nPOTIM = 2\n")
    assert not any("recommended" in n for n in generate("vasp", str(tpl), str(tmp_path / "d"), md=True)[1])
