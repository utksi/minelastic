# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""minelastic: second-order elastic constants from minimal uncoupled stress-strain calculations.

The elastic constants C_ij = dsigma_i/de_j are obtained from finite differences of the stress of a few
strained cells, with the minimal uncoupled strain sets of S. I. Simak and U. Singh, arXiv:2610.11733 (2026),
for any crystal symmetry. Inputs are generated and outputs read for VASP, Quantum ESPRESSO, ABINIT, FHI-aims
and any ASE format (generate, collect); ASE calculators are used directly (calculate).
"""
__version__ = "0.1.0"
__url__ = "https://github.com/utksi/minelastic"

from .core import ElasticResult, ElasticSetup  # noqa: E402
from .calculator import calculate  # noqa: E402
from .workflow import collect, generate  # noqa: E402

__all__ = ["ElasticSetup", "ElasticResult", "calculate", "generate", "collect", "__version__"]
