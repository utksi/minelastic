# Notes for coding agents

To compute elastic constants or other tensors with minelastic, follow
[`skills/minelastic/SKILL.md`](skills/minelastic/SKILL.md).

To change minelastic itself:

* `pip install -e ".[test,docs]"`, then `pytest tests` runs the tests, including every code block of
  `docs/`, `README.md` and `skills/` (`tests/test_docs.py`). No DFT code is needed: the tests use stored
  stresses of real calculations.
* `mkdocs build --strict` builds the documentation, and `python tools/check_docs.py` checks its typography,
  spelling and prose.
* Every Python file starts with the SPDX header and a docstring that states its purpose and the
  mathematics it implements; every function has a docstring with its formula, inputs, outputs and units.
* Text and code use ASCII punctuation: no typographic dashes, curly quotes, arrows or Unicode minus signs.
