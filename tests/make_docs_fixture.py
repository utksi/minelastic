# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Fixture of the documented Quantum ESPRESSO and ABINIT runs (docs/get-started/qe.md, abinit.md): wurtzite AlN.

`minelastic generate qe pw.in -o el --stencil central` and `minelastic generate abinit ref/ -o el --stencil
central` with the templates of examples/qe_aln and examples/abinit_aln, run with pw.x 7.6 and ABINIT 10.6.7.
Stored per cell, as minelastic's readers return them: the Cauchy stress (GPa, Voigt xx yy zz yz xz xy,
tension positive) and the forces (eV/A). The docs test replays `minelastic collect` from these numbers.

    QE_RUN=/path/to/qe_tests/aln_minimal ABINIT_RUN=/path/to/code_tests/abinit_aln/minimal python make_docs_fixture.py
"""
import json
import os

import numpy as np

from minelastic.workflow import CODES

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data_docs", "aln_dft_runs.json")
RUNS = {"qe": (os.environ.get("QE_RUN"), "examples/qe_aln/pw.in"),
        "abinit": (os.environ.get("ABINIT_RUN"), "examples/abinit_aln")}


def read_run(run, template):
    """Setup record and the stress (GPa, Voigt) and forces (eV/A) of every state of one run directory.

    The stress is rounded to 1e-6 GPa and the forces to 1e-6 eV/A, the precision the codes print.
    """
    setup = json.load(open(os.path.join(run, "minelastic.json")))
    setup["template"] = template                        # instead of the local absolute path
    stress, forces = {}, {}
    for sid in open(os.path.join(run, "STATES")).read().split():
        s, f = CODES[setup["code"]].read_result(os.path.join(run, sid), None)
        stress[sid], forces[sid] = np.round(s, 6).tolist(), np.round(f, 6).tolist()
    return {"setup": setup, "stress_GPa": stress, "forces_eV_A": forces}


def main():
    """Writes data_docs/aln_dft_runs.json with one entry per code."""
    data = {"label": "wurtzite AlN (PBE), minimal set, central stencil, delta = 0.01: pw.x 7.6 (pslibrary 1.0.0 "
                     "PAW, ecutwfc 50 Ry) and ABINIT 10.6.7 (PseudoDojo v0.4, ecut 48 Ha)",
            "source": "qe_tests/aln_minimal; code_tests/abinit_aln/minimal"}
    for code, (run, template) in RUNS.items():
        data[code] = read_run(run, template)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(data, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
