from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("scalefree")
except PackageNotFoundError:
    __version__ = "0.0.0"

from .vmoments import vprofile, ScaleFreeRunner, ScaleFreeResult
from .mock import mock

# Export the module as scalefree.hermite
from . import hermite as hermite

# Optional: also export the callable fitter directly under a non-colliding name
from .hermite import hermite as hermite_fit

__all__ = [
    "vprofile",
    "ScaleFreeRunner",
    "ScaleFreeResult",
    "mock",
    "hermite",  # module
    "hermite_fit",  # function (optional but useful)
    "__version__",
]
