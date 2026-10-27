# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Quantum ESPRESSO (pw.x): input files from a template pw.in, stresses and forces from the pw.x output.

The template is the input of the reference calculation. For every state it is rewritten with ibrav = 0,
the strained CELL_PARAMETERS (Angstrom) and ATOMIC_POSITIONS crystal (species labels and if_pos flags of the
template kept); all other namelist variables and cards are copied. The lattice of any ibrav is built as
pw.x's latgen does. pw.x prints the stress pressure-like in Ry/bohr^3; the reader returns the Cauchy stress
in GPa, tension positive: s = -14710.5 (m_xx, m_yy, m_zz, m_yz, m_xz, m_xy). Forces: Ry/bohr x 25.711 = eV/A.
"""
import os
import re

import numpy as np

from ..symmetry import SPECIES

DEFAULT_OUTPUT = "pw.out"
_CARDS = ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS",
          "ATOMIC_VELOCITIES", "ATOMIC_FORCES", "ADDITIONAL_K_POINTS", "SOLVENTS", "HUBBARD")
_CELL_KEYS = re.compile(r"^(celldm\(\d\)|a|b|c|cosab|cosac|cosbc)$", re.I)
RY_PER_BOHR3_IN_GPA = 14710.507848260711


def _split_assignments(text):
    """Namelist assignments split at commas and newlines outside quotes: "a = 1, b = 'x,y'" ->
    [('a', '1'), ('b', "'x,y'")]."""
    out, cur, q = [], "", None
    for ch in text:
        if q:
            cur += ch
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            cur += ch
        elif ch in ",\n":
            if cur.strip():
                out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    pairs = []
    for a in out:
        if "=" in a:
            k, v = a.split("=", 1)
            pairs.append((k.strip(), v.strip()))
    return pairs


def _strip_comment(line):
    """line without its comment: text from the first ! or # outside quotes."""
    q = None
    for k, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch in "!#":
            return line[:k]
    return line


def _namelist_end(text):
    """Index of the first '/' outside quotes in text, or -1."""
    q = None
    for k, ch in enumerate(text):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "/":
            return k
    return -1


def parse_input(text):
    """(namelists, cards) of a pw.x input: {NAME: [(key, value)]} in input order, and [(header, [lines])].
    Comments after ! or # are removed; a namelist runs from &NAME to the first '/' outside quotes."""
    lines = [_strip_comment(l).rstrip() for l in text.splitlines()]
    namelists, cards, i = {}, [], 0
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("&"):
            name = s[1:].split()[0].upper()
            body = s[len(name) + 1:]
            # WHY: a namelist ends at the first '/' outside quotes, also on the same line ('&CONTROL ... /')
            while _namelist_end(body) < 0 and i + 1 < len(lines):
                i += 1
                body += "\n" + lines[i]
            end = _namelist_end(body)
            namelists[name] = _split_assignments(body[:end] if end >= 0 else body)
        elif s and s.split()[0].upper() in _CARDS:
            cards.append((s, []))
        elif s and cards:
            cards[-1][1].append(lines[i])
        i += 1
    return namelists, cards


def _set(nl, key, value):
    """Set key (case-insensitive) in the namelist list nl, replacing the value or appending the pair."""
    for n, (k, _) in enumerate(nl):
        if k.lower() == key.lower():
            nl[n] = (k, value)
            return
    nl.append((key, value))


BOHR = 0.529177210903          # Angstrom (CODATA 2018, as ase.units.Bohr)


def _num(v):
    """Fortran number to float ('1.0d-8' -> 1e-8)."""
    return float(str(v).strip().lower().replace("d", "e"))


def _true(v):
    """Fortran logical ('.true.', '.t.', 'T', any case) to bool."""
    return str(v).strip().strip(".").lower().startswith("t")


def _get(pairs, key, default=None):
    """Value of key in a namelist list (case- and space-insensitive), or default."""
    for k, v in pairs:
        if k.replace(" ", "").lower() == key.lower():
            return v
    return default


def qe_lattice(system):
    """Lattice vectors (rows, Angstrom) for ibrav != 0, as pw.x's latgen builds them (INPUT_PW conventions).

    The parameters are celldm(1..6), or A, B, C, cosAB, cosAC, cosBC converted to celldm as QE's abc2celldm
    does. Supported: ibrav 1, 2, +-3, 4, +-5, 6, 7, 8, +-9, 91, 10, 11, +-12, +-13, 14.
    """
    ibrav = int(_get(system, "ibrav"))
    if _get(system, "celldm(1)") is not None:
        cd = [0.0] + [_num(_get(system, f"celldm({i})", 0)) for i in range(1, 7)]
    else:                                        # A, B, C, cosAB, ... -> celldm as abc2celldm of QE
        A = _num(_get(system, "A"))
        cd = [0.0, A / BOHR, _num(_get(system, "B", 0)) / A, _num(_get(system, "C", 0)) / A, 0.0, 0.0, 0.0]
        cosab, cosac, cosbc = (_num(_get(system, k, 0)) for k in ("cosAB", "cosAC", "cosBC"))
        if ibrav == 14:
            cd[4:7] = [cosbc, cosac, cosab]
        elif ibrav in (-12, -13):
            cd[5] = cosac
        else:
            cd[4] = cosab
    a = cd[1] * BOHR
    b, c = cd[2] * a, cd[3] * a
    cos5 = cosab = cd[4]                         # 5, -5: cos(gamma); 12, 13: cos(ab)
    cosac = cd[5]                                # -12, -13
    if ibrav == 14:
        cosbc, cosac, cosab = cd[4], cd[5], cd[6]
    r3 = np.sqrt(3.0)
    if ibrav == 1:
        v = [[a, 0, 0], [0, a, 0], [0, 0, a]]
    elif ibrav == 2:
        v = np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]]) * a / 2
    elif ibrav == 3:
        v = np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]]) * a / 2
    elif ibrav == -3:
        v = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]) * a / 2
    elif ibrav == 4:
        v = [[a, 0, 0], [-a / 2, a * r3 / 2, 0], [0, 0, c]]
    elif ibrav in (5, -5):
        tx, ty, tz = np.sqrt((1 - cos5) / 2), np.sqrt((1 - cos5) / 6), np.sqrt((1 + 2 * cos5) / 3)
        if ibrav == 5:
            v = np.array([[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]]) * a
        else:
            u, w = tz - 2 * np.sqrt(2) * ty, tz + np.sqrt(2) * ty
            v = np.array([[u, w, w], [w, u, w], [w, w, u]]) * a / r3
    elif ibrav == 6:
        v = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    elif ibrav == 7:
        v = np.array([[a, -a, c], [a, a, c], [-a, -a, c]]) / 2
    elif ibrav == 8:
        v = [[a, 0, 0], [0, b, 0], [0, 0, c]]
    elif ibrav == 9:
        v = [[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [0, 0, c]]
    elif ibrav == -9:
        v = [[a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c]]
    elif ibrav == 91:
        v = [[a, 0, 0], [0, b / 2, -c / 2], [0, b / 2, c / 2]]
    elif ibrav == 10:
        v = [[a / 2, 0, c / 2], [a / 2, b / 2, 0], [0, b / 2, c / 2]]
    elif ibrav == 11:
        v = [[a / 2, b / 2, c / 2], [-a / 2, b / 2, c / 2], [-a / 2, -b / 2, c / 2]]
    elif ibrav == 12:
        sg = np.sqrt(1 - cosab ** 2)
        v = [[a, 0, 0], [b * cosab, b * sg, 0], [0, 0, c]]
    elif ibrav == -12:
        sb = np.sqrt(1 - cosac ** 2)
        v = [[a, 0, 0], [0, b, 0], [c * cosac, 0, c * sb]]
    elif ibrav == 13:
        sg = np.sqrt(1 - cosab ** 2)
        v = [[a / 2, 0, -c / 2], [b * cosab, b * sg, 0], [a / 2, 0, c / 2]]
    elif ibrav == -13:
        sb = np.sqrt(1 - cosac ** 2)
        v = [[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [c * cosac, 0, c * sb]]
    elif ibrav == 14:
        sg = np.sqrt(1 - cosab ** 2)
        v = [[a, 0, 0], [b * cosab, b * sg, 0],
             [c * cosac, c * (cosbc - cosac * cosab) / sg,
              c * np.sqrt(1 + 2 * cosbc * cosac * cosab - cosbc ** 2 - cosac ** 2 - cosab ** 2) / sg]]
    else:
        raise ValueError(f"ibrav = {ibrav} is not supported")
    return np.array(v, float)


def _alat(system):
    """Lattice parameter alat in Angstrom (celldm(1) in bohr or A in Angstrom), or None."""
    if _get(system, "celldm(1)") is not None:
        return _num(_get(system, "celldm(1)")) * BOHR
    if _get(system, "A") is not None:
        return _num(_get(system, "A"))
    return None


def _element(label):
    """Chemical symbol of a species label, case-insensitive, as pw.x reads it: the first two letters if they
    form an element symbol, else the first letter ('Fe1', 'FE1', 'fe_up' -> 'Fe'; 'C1' -> 'C')."""
    from ase.data import chemical_symbols
    s = label.strip()
    if len(s) >= 2 and s[:2].isalpha() and s[:2].capitalize() in chemical_symbols:
        return s[:2].capitalize()
    if s[:1].upper() in chemical_symbols:
        return s[:1].upper()
    raise ValueError(f"cannot determine the element of species label {label!r}")


def read_structure(template):
    """Atoms of a pw.x input: lattice from CELL_PARAMETERS (alat, bohr or angstrom) or ibrav, positions in
    alat, bohr, angstrom or crystal units; atoms.arrays[SPECIES] is the index of each atom's species label.
    With noncolin, the initial moment of an atom of species i is the vector
    starting_magnetization(i) (sin a1 cos a2, sin a1 sin a2, cos a1), a1 = angle1(i), a2 = angle2(i) (degrees)."""
    from ase import Atoms
    namelists, cards = parse_input(open(template).read())
    system = namelists.get("SYSTEM", [])
    ibrav = int(_get(system, "ibrav", 0))
    alat = _alat(system)
    card = {h.split()[0].upper(): (h, body) for h, body in cards}
    if ibrav == 0:
        h, body = card["CELL_PARAMETERS"]
        unit = (h.split()[1] if len(h.split()) > 1 else ("alat" if alat else "bohr")).strip("{}()").lower()
        cell = np.array([[_num(x) for x in l.split()[:3]] for l in body if l.strip()][:3])
        cell *= {"angstrom": 1.0, "bohr": BOHR, "alat": alat or 0.0}[unit]
        if alat is None and unit == "alat":
            raise ValueError("CELL_PARAMETERS alat needs celldm(1) or A")
        alat = alat or np.linalg.norm(cell[0])
    else:
        cell = qe_lattice(system)
    species = {l.split()[0]: _element(l.split()[0]) for l in card["ATOMIC_SPECIES"][1] if l.strip()}
    h, body = card["ATOMIC_POSITIONS"]
    unit = (h.split()[1] if len(h.split()) > 1 else "alat").strip("{}()").lower()
    rows = [l.split() for l in body if l.strip()]
    xyz = np.array([[_num(x) for x in r[1:4]] for r in rows])
    symbols = [species[r[0]] for r in rows]
    if unit == "crystal":
        atoms = Atoms(symbols, cell=cell, scaled_positions=xyz, pbc=True)
    else:
        atoms = Atoms(symbols, cell=cell, positions=xyz * {"angstrom": 1.0, "bohr": BOHR, "alat": alat}[unit], pbc=True)
    # WHY: pw.x treats every species label as its own type, so atoms with different labels are not equivalent
    atoms.new_array(SPECIES, np.array([list(species).index(r[0]) for r in rows]))
    # WHY: the labels do not carry the moment directions, which lower the symmetry with spin-orbit coupling
    if _true(_get(system, "noncolin", ".false.")):
        per_species = []
        for i in range(1, len(species) + 1):
            m = _num(_get(system, f"starting_magnetization({i})", 0))
            a1, a2 = (np.radians(_num(_get(system, f"angle{j}({i})", 0))) for j in (1, 2))
            per_species.append(m * np.array([np.sin(a1) * np.cos(a2), np.sin(a1) * np.sin(a2), np.cos(a1)]))
        moments = np.array([per_species[n] for n in atoms.arrays[SPECIES]])
        if np.any(np.abs(moments) > 1e-8):
            atoms.set_initial_magnetic_moments(moments)
    return atoms


def render(namelists, cards):
    """pw.x input text from namelists (in dict order) and cards."""
    out = []
    for name, pairs in namelists.items():
        out.append(f"&{name}")
        out += [f"  {k} = {v}" for k, v in pairs]
        out.append("/")
    for head, body in cards:
        out.append(head)
        out += body
    return "\n".join(out) + "\n"


def prepare(setup, template, out_dir, relax_ions=False, input_name="pw.in", md=False):
    """One directory per state with a pw.x input: ibrav = 0, strained CELL_PARAMETERS (angstrom),
    ATOMIC_POSITIONS crystal, tstress and tprnfor on; all other variables and cards from the template.

    calculation = 'scf', or 'relax' with relax_ions (forc_conv_thr = 2e-5, etot_conv_thr = 1e-8 unless set);
    with md the template must have calculation = 'md', which is kept. If the cell was rotated to its standard
    orientation, direction-dependent inputs (K_POINTS tpiba, noncolin, lspinorb, angle1/angle2, electric
    fields, ATOMIC_VELOCITIES, ATOMIC_FORCES, CONSTRAINTS, if_pos flags other than 1 1 1) are refused.
    Returns notes for the user.
    """
    namelists, cards = parse_input(open(template).read())
    for need in ("CONTROL", "SYSTEM", "ELECTRONS"):
        namelists.setdefault(need, [])
    ctrl, system = namelists["CONTROL"], namelists["SYSTEM"]
    if md:
        if (_get(ctrl, "calculation") or "").strip("'\"").lower() != "md":
            raise ValueError("md: the template must set calculation = 'md' (fixed-cell molecular dynamics)")
    else:
        _set(ctrl, "calculation", "'relax'" if relax_ions else "'scf'")
    _set(ctrl, "tstress", ".true.")
    _set(ctrl, "tprnfor", ".true.")
    namelists["SYSTEM"] = [(k, v) for k, v in system if not _CELL_KEYS.match(k)]
    _set(namelists["SYSTEM"], "ibrav", "0")
    notes = []
    if relax_ions:
        namelists.setdefault("IONS", [])
        # WHY: the forces induced by the strain are ~ delta * |dF/de| ~ 0.05 eV/A for delta = 0.01; the pw.x
        #      defaults (forc_conv_thr 1e-3 Ry/bohr = 0.026 eV/A, etot_conv_thr 1e-4 Ry) stop the relaxation
        #      before it has relaxed the internal strain (seen for shear strains of AlN)
        for key, val, lim in (("forc_conv_thr", "2.0d-5", 1e-4), ("etot_conv_thr", "1.0d-8", 1e-6)):
            cur = _get(ctrl, key)
            if cur is None:
                _set(ctrl, key, val)
            elif _num(cur) > lim:
                notes.append(f"{key} = {cur} is loose for relaxed-ion constants (suggest {val})")
    namelists.pop("CELL", None)
    # reorder: CONTROL SYSTEM ELECTRONS IONS (others after)
    order = ["CONTROL", "SYSTEM", "ELECTRONS", "IONS"]
    namelists = {**{k: namelists[k] for k in order if k in namelists},
                 **{k: v for k, v in namelists.items() if k not in order}}
    pos_card = next((c for c in cards if c[0].split()[0].upper() == "ATOMIC_POSITIONS"), None)
    if pos_card is None:
        raise ValueError("template has no ATOMIC_POSITIONS card")
    rows = [l.split() for l in pos_card[1] if l.strip()]
    if len(rows) != len(setup.input_atoms):
        raise ValueError("ATOMIC_POSITIONS does not match the structure read by ASE")
    kcard = next((c for c in cards if c[0].split()[0].upper() == "K_POINTS"), None)
    tpiba = kcard is not None and kcard[0].split()[-1].lower().strip("{}()") in ("tpiba", "tpiba_b", "tpiba_c")
    if not np.allclose(setup.rotation, np.eye(3)):
        # WHY: Cartesian directions in the inputs would have to be rotated with the cell
        vec = [k for k, v in namelists["SYSTEM"] if k.lower().startswith(("angle1", "angle2", "lforcet")) or
               (k.lower() in ("noncolin", "lspinorb") and _true(v))]
        vec += [k for k, v in namelists["CONTROL"] if k.lower() in ("lelfield", "tefield", "lberry", "lorbm") and _true(v)]
        vec += [c[0].split()[0].upper() for c in cards
                if c[0].split()[0].upper() in ("ATOMIC_VELOCITIES", "ATOMIC_FORCES", "CONSTRAINTS")]
        if any(len(r) > 4 and r[4:7] != ["1", "1", "1"] for r in rows):
            vec.append("if_pos flags")
        if tpiba or vec:
            raise ValueError("the cell was rotated to its standard orientation (frame='auto'); "
                             f"{'K_POINTS tpiba' if tpiba else vec[0]} is direction-dependent and is not rotated: "
                             "use frame='input' (--frame input)")
    if tpiba:
        notes.append("K_POINTS in tpiba units are Cartesian: the strained cells use slightly shifted k points")
    other = [c for c in cards if c[0].split()[0].upper() not in ("ATOMIC_POSITIONS", "CELL_PARAMETERS")]
    os.makedirs(out_dir, exist_ok=True)
    for st in setup.states:
        at = setup.strained_atoms(st)
        frac = at.get_scaled_positions(wrap=False)
        cell = ["CELL_PARAMETERS angstrom"] + ["  " + " ".join(f"{x:16.10f}" for x in v) for v in at.cell[:]]
        pos = ["ATOMIC_POSITIONS crystal"] + [
            f"  {r[0]:4s} " + " ".join(f"{x:16.10f}" for x in f) + ("  " + " ".join(r[4:]) if len(r) > 4 else "")
            for r, f in zip(rows, frac)]
        d = os.path.join(out_dir, st["id"])
        os.makedirs(d, exist_ok=True)
        text = render(namelists, other)
        text += "\n".join(cell) + "\n" + "\n".join(pos) + "\n"
        open(os.path.join(d, input_name), "w").write(f"! minelastic strain state {st['id']}\n" + text)
    return notes


def read_output(path, all_steps=False):
    """(stress, forces) of the last step of a pw.x output, or all steps with all_steps (stress array, list of
    force arrays); ValueError if the numbers of stress and force blocks differ.

    Every 'total   stress' block (Ry/bohr^3, pressure-like) gives s in GPa, tension positive, Voigt; every
    'Forces acting on atoms' block gives the first nat force lines (Ry/bohr -> eV/Angstrom).
    """
    lines = open(path, errors="replace").read().splitlines()
    stresses, forces = [], []
    nat = None
    for i, l in enumerate(lines):
        if nat is None and "number of atoms/cell" in l:
            nat = int(l.split("=")[1])
        if "total   stress" in l:
            m = np.array([[float(x) for x in lines[i + r].split()[:3]] for r in (1, 2, 3)])
            # Ry/bohr^3, pressure-like -> physical stress in GPa
            stresses.append(-RY_PER_BOHR3_IN_GPA * np.array([m[0, 0], m[1, 1], m[2, 2], m[1, 2], m[0, 2], m[0, 1]]))
        if l.strip().startswith("Forces acting on atoms") and nat:
            f, j = [], i + 1
            while len(f) < nat and j < len(lines):
                if "force =" in lines[j] and "atom" in lines[j]:
                    f.append([float(x) for x in lines[j].split("=")[1].split()[:3]])
                j += 1
            if len(f) != nat:
                raise ValueError(f"{path}: incomplete force block ({len(f)} of {nat} atoms)")
            forces.append(np.array(f) * 25.71104309541616)       # Ry/bohr -> eV/A
    if not stresses:
        raise ValueError(f"no stress in {path} (tstress = .true. needed)")
    if forces and len(forces) != len(stresses):
        raise ValueError(f"{path}: {len(stresses)} stress and {len(forces)} force blocks (incomplete output?)")
    if all_steps:
        return np.array(stresses), forces
    return stresses[-1], (forces[-1] if forces else None)


def read_result(run_dir, output=None):
    """(stress, forces) of a finished run directory (pw.out unless output is given)."""
    return read_output(os.path.join(run_dir, output or DEFAULT_OUTPUT))
