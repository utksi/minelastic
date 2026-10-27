# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""File-based workflow: generate() writes one input directory per state, collect() reads the results.

generate() records the setup in minelastic.json and the state names in STATES; collect() rebuilds the
setup from minelastic.json, reads the stress (and forces, dipoles) of every state with the code's reader,
and calls ElasticSetup.compute. For MD runs the stress of each state is the set of samples of all steps
after skip; with samples, it is the set of results of the sub-directories (e.g. single points on MD
snapshots).
"""
import glob
import json
import os

import numpy as np

from .core import ElasticSetup
from .io import abinit, aims, generic, qe, vasp

CODES = {"vasp": vasp, "qe": qe, "aims": aims, "abinit": abinit, "generic": generic}
META = "minelastic.json"
RESULT = "elastic.json"


def read_structure(code, template, fmt=None):
    """Reference structure of a template: a directory (vasp, aims, abinit), a pw.x input (qe) or a file (generic)."""
    mod = CODES[code]
    return mod.read_structure(template, fmt) if code == "generic" else mod.read_structure(template)


def generate(code, template, out_dir, scheme="minimal", delta=0.01, stencil=None, symprec=1e-5,
             frame="auto", relax_ions=False, fmt=None, filename=None, extra_files=(), ionic=False,
             displacement=0.01, md=False, **code_kwargs):
    """Writes the inputs of every strained (and displaced-atom) cell from a template, as minelastic generate.

    One directory per state of ElasticSetup(reference structure of the template, ...) is written into out_dir,
    plus minelastic.json (the setup record that collect needs) and STATES (the state names, one per line).

    Parameters
    ----------
    code : str
        'vasp', 'qe', 'abinit', 'aims' or 'generic'.
    template : str
        vasp, abinit, aims: directory with the inputs of the relaxed reference structure; qe: pw.x input;
        generic: structure file.
    out_dir : str
        Output directory.
    scheme, delta, symprec, frame, displacement :
        As for ElasticSetup.
    stencil : str, optional
        As for ElasticSetup; by default 'four-point', and 'central' with md.
    relax_ions : bool
        Relax the atoms in each strained cell (relaxed-ion constants).
    ionic : bool
        Also the symmetry-irreducible atomic displacements; the clamped single points then give the
        relaxed-ion tensors, the Born charges and the piezoelectric response.
    md : bool
        Fixed-cell molecular dynamics in every strained cell as set up in the template (VASP: IBRION = 0 with
        NSW, POTIM and the thermostat; QE: calculation = 'md'); collect with md=True.
    fmt, filename, extra_files :
        generic interface: ASE format, structure file name, files copied into every directory.
    **code_kwargs :
        Options of the code's prepare (vasp: potcar='copy' or 'link', nsw, fix_incar).

    Returns
    -------
    (ElasticSetup, list of str)
        The setup and the notes: the changes made to the template and the recommended settings it misses.
    """
    if code not in CODES:
        raise ValueError(f"unknown code {code!r}; choose from {sorted(CODES)}")
    if ionic and relax_ions:
        raise ValueError("ionic=True obtains the relaxed-ion response from clamped cells; do not combine it with "
                         "relax_ions")
    if md and (relax_ions or ionic):
        raise ValueError("md writes molecular dynamics inputs; do not combine it with relax_ions or ionic")
    if md and code not in ("vasp", "qe"):
        raise ValueError("md inputs are written for VASP and Quantum ESPRESSO; for other codes run the MD (or "
                         "single points on its snapshots) yourself and collect them with samples=")
    if stencil is None:
        # WHY: the four-point stencil removes the O(delta^2) error of static stresses, but amplifies the noise of
        #      MD averages (weights (1, -8, 8, -1)/12); for MD, central differences are the better estimator
        stencil = "central" if md else "four-point"
    atoms = read_structure(code, template, fmt)
    setup = ElasticSetup(atoms, scheme=scheme, delta=delta, stencil=stencil, symprec=symprec, frame=frame,
                         ionic=ionic, displacement=displacement)
    os.makedirs(out_dir, exist_ok=True)
    if code == "generic":
        if not filename:
            raise ValueError("the generic interface needs the structure file name (filename=...)")
        notes = generic.prepare(setup, out_dir, filename, fmt=fmt, extra_files=extra_files, relax_ions=relax_ions)
    else:
        notes = CODES[code].prepare(setup, template, out_dir, relax_ions=relax_ions,
                                    **(dict(md=True) if md else {}), **code_kwargs)
    meta = setup.to_dict()
    meta.update(code=code, template=os.path.abspath(template), relax_ions=bool(relax_ions), md=bool(md),
                structure_format=fmt, structure_file=filename, notes=notes)
    with open(os.path.join(out_dir, META), "w") as f:
        json.dump(meta, f, indent=1)
    with open(os.path.join(out_dir, "STATES"), "w") as f:
        f.write("\n".join(st["id"] for st in setup.states) + "\n")
    return setup, notes


def collect(out_dir, method="uncoupled", output=None, fmt=None, md=False, skip=0, kinetic_temperature=None,
            samples=None, stress_measure="cauchy", polarization=False, born_file=None, save=True):
    """Elastic tensor from the finished calculations in out_dir (as written by generate), as minelastic collect.

    Reads the stress (and forces, and dipoles) of every state and returns ElasticSetup.compute of them.

    Parameters
    ----------
    out_dir : str
        Directory written by generate.
    method : str
        'uncoupled', 'designated' or 'lstsq' (see ElasticSetup.compute).
    output, fmt :
        Output file name in each directory (default per code) and, for generic, its ASE format.
    md : bool
        The stress of each state is the set of all MD steps after the first skip (VASP, QE); the standard
        errors assume uncorrelated steps, so thin or block-average long trajectories.
    skip : int
        Number of initial MD steps discarded.
    kinetic_temperature, stress_measure :
        Passed to ElasticSetup.compute.
    samples : str, optional
        Glob of sub-directories of each state (e.g. 'conf-*'), whose stresses are the samples.
    polarization : bool
        Read the Berry-phase dipoles (VASP LCALCPOL) of all cells: e0 from the strained cells and, with
        displaced atoms, the Born charges.
    born_file : str, optional
        OUTCAR of a VASP LEPSILON or LCALCEPS run of the reference structure; its Born charges and e0, rotated
        to the working frame, are used for the ionic terms.
    save : bool
        Also write the result to out_dir/elastic.json.

    Returns
    -------
    ElasticResult
    """
    meta = json.load(open(os.path.join(out_dir, META)))
    setup = ElasticSetup.from_dict(meta)
    born = clamped = None
    if born_file:
        from .ionic import rotate_born, rotate_piezo
        R = setup.rotation                               # x_work = R x_input
        born = rotate_born(vasp.read_born_charges(born_file), R)
        e0 = vasp.read_clamped_piezo(born_file)
        clamped = None if e0 is None else rotate_piezo(e0, R)
        if len(born) != len(setup.atoms):
            raise ValueError(f"{born_file}: {len(born)} Born charges for {len(setup.atoms)} atoms")
    code = meta["code"]
    if meta.get("md") and not (md or samples):
        raise ValueError(f"{out_dir} holds molecular dynamics runs: average them with md=True (--md --skip N) "
                         "or use samples=")
    if (md or samples) and setup.ionic:
        raise ValueError("the atomic displacements (ionic) are for clamped single-point calculations, not MD")
    stresses, forces, dipoles = {}, {}, {}
    for st in setup.states:
        d = os.path.join(out_dir, st["id"])
        if samples:
            subs = sorted(p for p in glob.glob(os.path.join(d, samples)) if os.path.isdir(p))
            if not subs:
                raise FileNotFoundError(f"no sample directories {samples!r} in {d}")
            if code == "generic":
                stresses[st["id"]] = np.array([generic.read_result(p, output, fmt=fmt)[0] for p in subs])
            else:
                stresses[st["id"]] = np.array([CODES[code].read_result(p, output)[0] for p in subs])
        elif md:
            if code == "vasp":
                p = os.path.join(d, output or "vasprun.xml")
                s, _ = (vasp.read_vasprun if p.endswith(".xml") else vasp.read_outcar)(p, all_steps=True)
            elif code == "qe":
                s, _ = qe.read_output(os.path.join(d, output or qe.DEFAULT_OUTPUT), all_steps=True)
            else:
                raise ValueError("md averaging is implemented for VASP and Quantum ESPRESSO outputs")
            if len(s) <= skip:
                raise ValueError(f"{d}: only {len(s)} steps, skip = {skip}")
            stresses[st["id"]] = np.asarray(s)[skip:]
        else:
            if code == "generic":
                s, f = generic.read_result(d, output or meta.get("output") or "OUTCAR", fmt=fmt)
            else:
                s, f = CODES[code].read_result(d, output)
            stresses[st["id"]] = s
            forces[st["id"]] = f
            if polarization:
                if code != "vasp":
                    raise ValueError("reading polarizations is implemented for VASP (LCALCPOL); for other codes "
                                     "pass polarization= to ElasticSetup.compute")
                dipoles[st["id"]] = vasp.read_dipole(d)
    use_forces = (not md and not samples and not meta.get("relax_ions") and forces
                  and all(v is not None for v in forces.values()))
    result = setup.compute(stresses, forces=forces if use_forces else None, method=method,
                           kinetic_temperature=kinetic_temperature, stress_measure=stress_measure,
                           polarization=dipoles if polarization else None, born=born, clamped_piezo=clamped)
    if save:
        d = result.to_dict()
        d.update(spacegroup=meta["spacegroup"], elastic_type=meta["elastic_type"], code=code,
                 relax_ions=meta.get("relax_ions"), md=md, skip=skip, samples=samples,
                 kinetic_temperature=kinetic_temperature)
        with open(os.path.join(out_dir, RESULT), "w") as f:
            json.dump(d, f, indent=1)
    return result
