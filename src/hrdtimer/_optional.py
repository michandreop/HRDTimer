"""Helpers for optional dependencies that are not installable from PyPI.

The signature-fitting stages depend on MuSiCal (GitHub only) and on
SigProfiler reference-genome data, neither of which can be expressed as
standard PyPI metadata. These helpers turn a missing dependency into a clear,
actionable error instead of a bare ``ModuleNotFoundError``.
"""

from __future__ import annotations

import importlib
from types import ModuleType

_HINTS = {
    "musical": (
        "MuSiCal is required for signature fitting and is not on PyPI. Clone and\n"
        "install it EDITABLE (a plain 'pip install git+...' drops its data files):\n"
        "  git clone --branch v1.0.0 https://github.com/parklab/MuSiCal.git\n"
        "  pip install -e MuSiCal\n"
        "or run setup_environment.sh."
    ),
    "SigProfilerAssignment": (
        "SigProfilerAssignment is required for signature assignment. Install the "
        "signature extras with:\n  pip install 'hrdtimer[signatures]'"
    ),
    "SigProfilerMatrixGenerator": (
        "SigProfilerMatrixGenerator is required for mutation-matrix generation. "
        "Install the signature extras with:\n  pip install 'hrdtimer[signatures]'"
    ),
}


def require(module: str) -> ModuleType:
    """Import ``module`` or raise ImportError with an install hint."""
    try:
        return importlib.import_module(module)
    except ImportError as exc:
        hint = _HINTS.get(module, f"Optional dependency '{module}' is not installed.")
        raise ImportError(hint) from exc
