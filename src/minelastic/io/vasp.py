# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""VASP: input directories from a template, stresses and forces from vasprun.xml/OUTCAR, IBRION = 6 tensors.

VASP writes the stress pressure-like (positive in compression) in kB, in the order XX YY ZZ XY YZ ZX. The
readers return the Cauchy stress in GPa, positive in tension, in Voigt order: s = -0.1 (XX, YY, ZZ, YZ, ZX, XY).
Elastic tensors are printed in kBar in the same order and are converted by C = 0.1 M[P, P] (P = Voigt
permutation); dipoles and Born charges are read in e Angstrom and |e|.
"""
import math
import os
import re
import shutil
import xml.etree.ElementTree as ET

import numpy as np
from ase.io import read

from ..symmetry import SPECIES

KB_TO_GPA = 0.1
# VASP order XX YY ZZ XY YZ ZX -> Voigt xx yy zz yz xz xy
VASP_TO_VOIGT = [0, 1, 2, 4, 5, 3]
# tags removed from the template INCAR: set by minelastic or specific to IBRION=5-8 / the VASP patch
_DROP = {"IBRION", "NSW", "ISIF", "NFREE", "POTIM", "LMINIMAL_ELASTIC", "LEPSILON", "LCALCEPS",
         "KSPACING", "KGAMMA"}

DEFAULT_OUTPUT = "vasprun.xml"


# ----------------------------------------------------------------------------------------- INCAR
def read_incar(path):
    """Ordered list of (TAG, value) of an INCAR (comments after # or ! dropped; ';' separates tags)."""
    entries = []
    for line in open(path):
        body = re.split(r"[#!]", line, maxsplit=1)[0]
        for part in body.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                entries.append((k.strip().upper(), v.strip()))
    return entries


def write_incar(path, entries, header="minelastic"):
    """INCAR with a comment line and one 'TAG = value' line per entry, in the given order."""
    with open(path, "w") as f:
        f.write(f"# {header}\n")
        for k, v in entries:
            f.write(f"{k} = {v}\n")


def potcar_enmax(path):
    """Largest ENMAX (eV) of the POTCAR at path, or None if the file is missing or has no ENMAX."""
    if not os.path.exists(path):
        return None
    vals = re.findall(r"ENMAX\s*=\s*([0-9.]+)", open(path, errors="replace").read())
    return max(float(v) for v in vals) if vals else None


def recommended_settings(tags, enmax, md=False):
    """{TAG: (template value or None, suggested value)} for the INCAR settings that stress differences of a few
    GPa need and the template does not meet: EDIFF <= 1E-6 (suggested 1E-8; with md <= 1E-5, suggested 1E-6;
    VASP's default is 1E-4), ENCUT >= 1.3 ENMAX (Pulay stress of the strained cells; VASP's default is ENMAX,
    suggested 1.3 ENMAX rounded up to 10 eV), PREC = Accurate or High (default Normal) and LREAL = .FALSE.
    (real-space projection adds noise to forces and stresses)."""
    num = lambda v: float(v.split()[0].lower().replace("d", "e"))          # noqa: E731
    out = {}
    limit, suggest = (1e-5, "1E-6") if md else (1e-6, "1E-8")
    if "EDIFF" not in tags or num(tags["EDIFF"]) > limit:
        out["EDIFF"] = (tags.get("EDIFF"), suggest)
    if enmax is not None and ("ENCUT" not in tags or num(tags["ENCUT"]) < 1.3 * enmax - 1e-6):
        out["ENCUT"] = (tags.get("ENCUT"), str(math.ceil(1.3 * enmax / 10 - 1e-9) * 10))
    if tags.get("PREC", "Normal").strip().lower()[:1] not in ("a", "h"):
        out["PREC"] = (tags.get("PREC"), "Accurate")
    if "LREAL" in tags and tags["LREAL"].strip().strip(".").lower() not in ("false", "f"):
        out["LREAL"] = (tags["LREAL"], ".FALSE.")
    return out


def kspacing_grid(atoms, kspacing):
    """The Gamma-centred grid VASP uses for KSPACING: N_i = max(1, ceil(|b_i| / KSPACING)), with the reciprocal
    lattice vectors b_i = 2 pi (A^-1)^T rows (including the factor 2 pi), in 1/Angstrom."""
    b = 2 * math.pi * np.linalg.inv(np.asarray(atoms.cell[:])).T
    return [max(1, math.ceil(np.linalg.norm(bi) / kspacing - 1e-8)) for bi in b]


def poscar_species_blocks(path):
    """Index of the species block of every atom of a POSCAR (from the counts line), or None if unreadable."""
    lines = open(path).read().splitlines()
    k = 6 if not lines[5].split()[0].lstrip("+-").isdigit() else 5          # counts line (VASP 5 / 4)
    counts = []
    for x in lines[k].split():
        if not x.isdigit():
            break
        counts.append(int(x))
    return np.repeat(np.arange(len(counts)), counts) if counts else None


def _true(value):
    """INCAR logical: .TRUE., T, True -> True."""
    return value.split()[0].lstrip(".").upper().startswith("T")


def _kpoints_cartesian(path):
    """True if a KPOINTS file gives k points or generating vectors in Cartesian coordinates (mode 'C' or 'K')."""
    lines = [l.strip() for l in open(path).read().splitlines()]
    mode = lines[2][:1].lower() if len(lines) > 2 and lines[2] else ""
    if mode == "l":                                  # line mode: the coordinate type is on the next line
        mode = lines[3][:1].lower() if len(lines) > 3 and lines[3] else ""
    return mode in ("c", "k")


def read_structure(template):
    """Reference structure of a template directory with the symmetry information VASP uses: the index of the
    POSCAR species block of every atom (atoms.arrays[SPECIES]) and, for ISPIN = 2, LNONCOLLINEAR or LSORBIT,
    MAGMOM as the initial magnetic moments (N values, or 3N Cartesian components for noncollinear runs with
    the default SAXIS; 'n*value' expanded)."""
    atoms = read(os.path.join(template, "POSCAR"), format="vasp")
    blocks = poscar_species_blocks(os.path.join(template, "POSCAR"))
    if blocks is not None and len(blocks) == len(atoms):
        atoms.new_array(SPECIES, blocks)
    incar = os.path.join(template, "INCAR")
    tags = dict(read_incar(incar)) if os.path.exists(incar) else {}
    spin = tags.get("ISPIN", "1").split()[0] == "2"
    noncollinear = _true(tags.get("LNONCOLLINEAR", "F")) or _true(tags.get("LSORBIT", "F"))
    # WHY: noncollinear MAGMOM components refer to the spin quantization axis SAXIS; only the default
    #      SAXIS = (0, 0, 1), with Cartesian components, is read
    if noncollinear and "MAGMOM" in tags and "SAXIS" in tags and \
            not np.allclose([float(x) for x in tags["SAXIS"].split()[:3]], [0, 0, 1]):
        raise ValueError("noncollinear MAGMOM with SAXIS other than (0, 0, 1): give MAGMOM in Cartesian "
                         "components with the default SAXIS")
    if "MAGMOM" in tags and (spin or noncollinear):
        m = []
        for x in tags["MAGMOM"].split():
            n, _, v = x.rpartition("*")
            m += [float(v)] * (int(n) if n else 1)
        if len(m) == len(atoms) and not noncollinear:
            atoms.set_initial_magnetic_moments(m)
        elif len(m) == 3 * len(atoms) and noncollinear:
            atoms.set_initial_magnetic_moments(np.reshape(m, (len(atoms), 3)))
    return atoms


def write_poscar(path, template_poscar, atoms, comment):
    """POSCAR of ``atoms`` with the species/count lines and selective-dynamics flags of the template.

    WHY: keeps the species blocks exactly as in the template (and hence as in its POTCAR), including
    repeated species, independently of the ASE version.
    """
    t = open(template_poscar).read().splitlines()
    i = 5
    species = t[5] if not t[5].split()[0].lstrip("+-").isdigit() else None
    if species is not None:
        i = 6
    counts = t[i]
    n = sum(int(x) for x in counts.split())
    sel = t[i + 1].strip()[:1].lower() == "s"
    first = i + 2 + int(sel)
    flags = [l.split()[3:6] if sel else [] for l in t[first:first + n]]
    # WHY: the flags refer to the coordinates as given; the cells are written in direct coordinates
    if sel and t[i + 2].strip()[:1].lower() in ("c", "k") and any("F" in f.upper() for fl in flags for f in fl):
        raise ValueError("template POSCAR: selective dynamics with Cartesian coordinates; give direct coordinates")
    if n != len(atoms):
        raise ValueError("template POSCAR and structure differ in the number of atoms")
    out = [comment, "   1.0"]
    out += ["  " + " ".join(f"{x:21.16f}" for x in v) for v in atoms.cell[:]]
    if species is not None:
        out.append(species)
    out.append(counts)
    if sel:
        out.append("Selective dynamics")
    out.append("Direct")
    for f, fl in zip(atoms.get_scaled_positions(wrap=False), flags):
        out.append("  " + " ".join(f"{x:19.16f}" for x in f) + ("  " + " ".join(fl) if fl else ""))
    open(path, "w").write("\n".join(out) + "\n")


def prepare(setup, template, out_dir, relax_ions=False, nsw=99, potcar="copy", md=False, fix_incar=False):
    """One VASP directory per state in out_dir: POSCAR of the strained cell, INCAR, KPOINTS and POTCAR.

    The template holds INCAR, POTCAR, POSCAR (reference structure) and KPOINTS or KSPACING. The INCAR is
    copied without IBRION, NSW, ISIF, NFREE, POTIM, LEPSILON, LCALCEPS, KSPACING and KGAMMA, and gets
    IBRION = -1, NSW = 0, ISIF = 2 (single point with stress), or IBRION = 2, NSW = nsw, ISIF = 2 and
    EDIFFG = -5E-4 unless set (ions relaxed in the fixed cell), or, with md, the template's IBRION = 0, NSW
    and POTIM with ISIF = 2 (fixed-cell MD). KSPACING becomes the fixed Gamma-centred grid of the reference
    cell, so that all cells use the same k points. The notes list every tag set or removed, and the settings
    of recommended_settings that the template does not meet; fix_incar applies all of those. If the cell
    was rotated to its standard orientation, direction-dependent inputs (LNONCOLLINEAR, LSORBIT, SAXIS,
    EFIELD, EFIELD_PEAD, Cartesian KPOINTS) are refused, as they are not rotated. Returns the notes.
    """
    entries = read_incar(os.path.join(template, "INCAR"))
    tags = dict(entries)
    notes = []
    kpoints = os.path.join(template, "KPOINTS")
    if setup.rotation is not None and not np.allclose(setup.rotation, np.eye(3)):
        # WHY: Cartesian directions in the inputs would have to be rotated with the cell
        for t in ("LNONCOLLINEAR", "LSORBIT", "SAXIS", "EFIELD", "EFIELD_PEAD"):
            if t in tags and (t in ("SAXIS", "EFIELD", "EFIELD_PEAD") or _true(tags[t])):
                raise ValueError(f"the cell was rotated to its standard orientation (frame='auto'); {t} is "
                                 "direction-dependent and is not rotated: use frame='input' (--frame input)")
        if os.path.exists(kpoints) and _kpoints_cartesian(kpoints):
            raise ValueError("the cell was rotated to its standard orientation (frame='auto'); KPOINTS has Cartesian "
                             "k points, which are not rotated: use frame='input' (--frame input)")
    kgrid = None
    if not os.path.exists(kpoints):
        if "KSPACING" not in tags:
            raise FileNotFoundError("template has neither KPOINTS nor KSPACING")
        kgrid = kspacing_grid(setup.atoms, float(tags["KSPACING"].split()[0]))
        notes.append(f"KSPACING replaced by the fixed grid {kgrid} of the reference cell")
    md_kept = ("IBRION", "NSW", "POTIM") if md else ()
    keep = [(k, v) for k, v in entries if k not in _DROP or k in md_kept]
    if md:
        # WHY: finite T: fixed-cell MD in each strained cell; ISIF = 2 keeps the cell (ISIF = 3 would be NPT
        #      with the Langevin thermostat) and still gives the stress at every step
        if tags.get("IBRION", "").split()[:1] != ["0"]:
            raise ValueError("md: the template INCAR must set up molecular dynamics (IBRION = 0)")
        added = [("ISIF", "2")]
    elif relax_ions:
        added = [("IBRION", "2"), ("NSW", str(nsw)), ("ISIF", "2")]
        # WHY: residual forces must be small against the strain-induced forces (~0.05 eV/A at delta = 0.01)
        if "EDIFFG" not in tags:
            added.append(("EDIFFG", "-5E-4"))
        elif not (float(tags["EDIFFG"].split()[0]) < 0 and abs(float(tags["EDIFFG"].split()[0])) <= 1e-3):
            notes.append(f"EDIFFG = {tags['EDIFFG']} may be too loose for relaxed-ion constants (suggest -5E-4)")
    else:
        added = [("IBRION", "-1"), ("NSW", "0"), ("ISIF", "2")]
    keep += added
    # report the tags that minelastic sets or removes (KSPACING and KGAMMA are reported with the k grid)
    new = dict(added)
    changes = [f"{k} = {v}" + (f" (template: {tags[k]})" if k in tags else "") for k, v in added
               if tags.get(k, "").split()[:1] != [v]]
    removed = [k for k, _ in entries if k in _DROP and k not in new and k not in md_kept
               and k not in ("KSPACING", "KGAMMA")]
    if changes or removed:
        notes.append("INCAR: " + "; ".join(([f"set {', '.join(changes)}"] if changes else [])
                                            + ([f"removed {', '.join(removed)}"] if removed else [])))
    # settings recommended for small stress differences: all applied with fix_incar, else none
    enmax = potcar_enmax(os.path.join(template, "POTCAR"))
    rec = recommended_settings(tags, enmax, md=md)
    if rec:
        default = {"EDIFF": "default 1E-4", "ENCUT": f"default ENMAX = {enmax:g}" if enmax else "",
                   "PREC": "default Normal"}
        text = ", ".join(f"{k} {old if old is not None else default[k]} -> {val}" for k, (old, val) in rec.items())
        if fix_incar:
            keep = [(k, rec[k][1] if k in rec else v) for k, v in keep]
            keep += [(k, val) for k, (old, val) in rec.items() if old is None]
            notes.append(f"INCAR: applied the recommended settings: {text}")
        else:
            notes.append(f"INCAR: recommended for small stress differences: {text}. The template values are kept; "
                         "--fix-incar (fix_incar=True) applies all of these")
    ref = read_structure(template)
    os.makedirs(out_dir, exist_ok=True)
    for st in setup.states:
        d = os.path.join(out_dir, st["id"])
        os.makedirs(d, exist_ok=True)
        at = setup.strained_atoms(st)
        write_poscar(os.path.join(d, "POSCAR"), os.path.join(template, "POSCAR"), at,
                     f"minelastic {st['id']} strain {np.round(st['strain'], 8).tolist()}")
        write_incar(os.path.join(d, "INCAR"), keep, header=f"minelastic strain state {st['id']}")
        if kgrid is None:
            shutil.copy(kpoints, os.path.join(d, "KPOINTS"))
        else:
            open(os.path.join(d, "KPOINTS"), "w").write("minelastic (from KSPACING)\n0\nGamma\n"
                                                          + " ".join(map(str, kgrid)) + "\n0 0 0\n")
        dst = os.path.join(d, "POTCAR")
        if os.path.lexists(dst):
            os.remove(dst)
        if potcar == "link":
            os.symlink(os.path.abspath(os.path.join(template, "POTCAR")), dst)
        else:
            shutil.copy(os.path.join(template, "POTCAR"), dst)
    assert len(ref) == len(setup.input_atoms)
    return notes


# --------------------------------------------------------------------------------------- outputs
def _fixed(line, decimals):
    """Numbers of a fixed-format line with the given number of decimals, also where adjacent fields run
    together ('-123456.78901-234567.89012', '1234567.8901234567.8901')."""
    return [float(x) for x in re.findall(r"[-+]?\d*\.\d{%d}" % decimals, line)]


def _voigt_from_vasp3x3(m):
    """Voigt stress in GPa, tension positive, from VASP's 3x3 stress in kB (pressure-like): s = -0.1 (m_xx,
    m_yy, m_zz, m_yz, m_xz, m_xy)."""
    m = np.asarray(m, float)
    return -KB_TO_GPA * np.array([m[0, 0], m[1, 1], m[2, 2], m[1, 2], m[0, 2], m[0, 1]])


def read_vasprun(path, all_steps=False):
    """(stress, forces) of the last ionic step of vasprun.xml, or arrays over all steps with all_steps.

    Stress in GPa (tension positive, Voigt), forces in eV/Angstrom (None when not written).
    """
    stresses, forces = [], []
    for _, el in ET.iterparse(path, events=("end",)):
        if el.tag == "calculation":
            s = el.find("./varray[@name='stress']")
            f = el.find("./varray[@name='forces']")
            if s is not None:
                stresses.append(_voigt_from_vasp3x3([[float(x) for x in v.text.split()] for v in s]))
                forces.append(None if f is None else np.array([[float(x) for x in v.text.split()] for v in f]))
            el.clear()
    if not stresses:
        raise ValueError(f"no stress tensor in {path} (ISIF >= 1 needed)")
    if all_steps:
        return np.array(stresses), forces
    return stresses[-1], forces[-1]


def read_outcar(path, all_steps=False):
    """(stress, forces) of the last ionic step of an OUTCAR ('in kB' line and TOTAL-FORCE block), or all
    steps with all_steps (stress array, list of force arrays); units and signs as read_vasprun. ValueError if
    the numbers of stress and force blocks differ."""
    stresses, forces = [], []
    lines = open(path, errors="replace").read().splitlines()
    nions = None
    for i, line in enumerate(lines):
        if nions is None and "NIONS =" in line:
            nions = int(line.split("NIONS =")[1].split()[0])
        if line.strip().startswith("in kB"):
            v = np.array(_fixed(line, 5)[:6])                    # 6F12.5
            stresses.append(-KB_TO_GPA * v[VASP_TO_VOIGT])
        if "TOTAL-FORCE (eV/Angst)" in line and nions:
            forces.append(np.array([[float(x) for x in lines[i + 2 + n].split()[3:6]] for n in range(nions)]))
    if not stresses:
        raise ValueError(f"no stress tensor in {path}")
    if forces and len(forces) != len(stresses):
        raise ValueError(f"{path}: {len(stresses)} stress and {len(forces)} force blocks (incomplete output?)")
    if all_steps:
        return np.array(stresses), forces
    return stresses[-1], (forces[-1] if forces else None)


def read_result(run_dir, output=None):
    """(stress, forces) of a finished run directory: vasprun.xml if present, else OUTCAR."""
    for name in ([output] if output else ["vasprun.xml", "OUTCAR"]):
        p = os.path.join(run_dir, name)
        if os.path.exists(p):
            try:
                return (read_vasprun if name.endswith(".xml") else read_outcar)(p)
            except (ET.ParseError, ValueError):
                if output:
                    raise
    raise FileNotFoundError(f"no readable vasprun.xml or OUTCAR in {run_dir}")


def read_elastic_tensor(outcar, which="TOTAL ELASTIC MODULI"):
    """6x6 elastic tensor printed by IBRION = 6, ISIF >= 3, in GPa and Voigt order (last block of that name).

    which: 'TOTAL ELASTIC MODULI' (relaxed ions), 'SYMMETRIZED ELASTIC MODULI' (clamped ions) or
    'ELASTIC MODULI CONTR FROM IONIC RELAXATION'. C = 0.1 M[P, P] with the printed kBar matrix M.
    """
    lines = open(outcar, errors="replace").read().splitlines()
    idx = [i for i, l in enumerate(lines) if l.strip().startswith(which)]
    if not idx:
        raise ValueError(f"'{which}' not found in {outcar}")
    i = idx[-1]
    while not lines[i].strip().startswith("XX"):
        i += 1
    M = np.array([_fixed(lines[i + r], 4)[:6] for r in range(6)])       # rows 'XX  6F12.4'
    return KB_TO_GPA * M[np.ix_(VASP_TO_VOIGT, VASP_TO_VOIGT)]


# WHY: the dipoles in the OUTCAR have the physical sign: the clamped-ion piezoelectric tensor of AlN from the
#      Berry-phase dipoles of strained cells agrees in sign and to 1% with the LEPSILON (DFPT) tensor
VASP_DIPOLE_SIGN = 1.0


def read_dipole(run_dir, output="OUTCAR"):
    """Dipole of the cell d = p[elc] + p[ion] (e Angstrom, physical sign) from an LCALCPOL = .TRUE. run.

    VASP evaluates the Berry phase for the first ionic step only. If the run has more ionic steps than the
    one with the dipole, the dipole belongs to the unrelaxed structure and ValueError is raised.
    """
    pelc = pion = None
    nsteps = steps_at_dipole = 0
    for line in open(os.path.join(run_dir, output), errors="replace"):
        if "TOTAL-FORCE (eV/Angst)" in line:
            nsteps += 1
        elif "p[elc]=(" in line:
            pelc = np.array([float(x) for x in line.split("(")[1].split(")")[0].split()])
            steps_at_dipole = nsteps
        elif "p[ion]=(" in line:
            pion = np.array([float(x) for x in line.split("(")[1].split(")")[0].split()])
    if pelc is None or pion is None:
        raise ValueError(f"no Berry-phase dipole in {run_dir}/{output} (LCALCPOL = .TRUE. needed)")
    if nsteps > 1 and steps_at_dipole < nsteps:
        raise ValueError(f"{run_dir}: the Berry-phase dipole was computed for ionic step {steps_at_dipole + 1} of "
                         f"{nsteps} (VASP computes it for the first step only); run a single point on the final "
                         "structure")
    return VASP_DIPOLE_SIGN * (pelc + pion)


def read_born_charges(outcar):
    """Born charges Z[n, a, b] = dp_a / du_nb (N, 3, 3) in |e| from a LEPSILON or LCALCEPS run (the last
    block 'including local field effects'), in the Cartesian frame of that run.

    VASP prints the polarization (field) direction a as the row and the displacement direction b as the
    column, so the block is read without transposing. Checked against Berry-phase charges dp_a/du_b for
    alpha-quartz, whose Born tensors are not symmetric (0.01 |e|; the transpose differs by 0.48 |e|).
    """
    lines = open(outcar, errors="replace").read().splitlines()
    idx = [i for i, l in enumerate(lines) if "BORN EFFECTIVE CHARGES" in l and "including local field" in l]
    if not idx:
        raise ValueError(f"no Born effective charges in {outcar} (LEPSILON or LCALCEPS needed)")
    i, out = idx[-1], []
    while True:
        j = next((k for k in range(i + 1, min(i + 8, len(lines))) if lines[k].strip().startswith("ion")), None)
        if j is None:
            break
        out.append([[float(x) for x in lines[j + 1 + r].split()[1:4]] for r in range(3)])
        i = j + 3
    return np.array(out)


def read_clamped_piezo(outcar):
    """Clamped-ion piezoelectric tensor e0 (3, 6) in C/m^2 from a LEPSILON or LCALCEPS run (rows x y z,
    columns reordered from XX YY ZZ XY YZ ZX to Voigt), in the frame of that run; None if not printed."""
    lines = open(outcar, errors="replace").read().splitlines()
    idx = [i for i, l in enumerate(lines)
           if l.strip().startswith("PIEZOELECTRIC TENSOR (including local field") and "C/m^2" in l]
    if not idx:
        return None
    rows = {}
    for l in lines[idx[-1] + 1:idx[-1] + 8]:                 # header XX YY ..., optional dashes, rows x y z
        f = l.split()
        if f and f[0] in ("x", "y", "z") and f[0] not in rows:
            rows[f[0]] = _fixed(l, 5)[:6]
    if len(rows) != 3:
        raise ValueError(f"{outcar}: incomplete piezoelectric tensor block")
    return np.array([rows["x"], rows["y"], rows["z"]])[:, VASP_TO_VOIGT]
