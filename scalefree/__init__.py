"""
scalefree package

Public API:
- vprofile: main user-facing function to compute moments/VPs
via the Fortran backend
- ScaleFreeRunner: advanced use
(reuse a runner instance, custom workdir, etc.)
- ScaleFreeResult: structured return container
"""

from .vmoments import vprofile, ScaleFreeRunner, ScaleFreeResult

__all__ = ["vprofile", "ScaleFreeRunner", "ScaleFreeResult"]
