"""
scalefree.vmoments

Prompt-driven Python wrapper for the ScaleFree Fortran executable.

Behavior regarding the Fortran backend
-------------------------------------
- If exe_path is provided: use it; if missing and gfortran exists,
we can build it there.
- If exe_path is not provided:
  1) use SCALEFREE_EXE env var if set
  2) else use a cached executable in a user cache directory
  3) else, if gfortran exists, auto-compile from packaged
  fortran_src/scalefree.f
  4) else raise a clear, actionable error instructing how to install gfortran

Behavior regarding output (Option A)
-----------------------------------
- Default: do NOT leave output files behind.
- We still answer the Fortran "Output file" prompt with a temporary filename
  (so the interactive program never blocks), but we prefer parsing the
  structured
  "# kind=..." blocks from STDOUT. Any temporary file is deleted at the end.
- If the caller provides output_path, we treat it as an explicit
request to write
  a persistent file and we parse that file (more stable for regression tests).

Note on pip install messages
----------------------------
Reliable messages at *pip install time* are not guaranteed for wheels.
We therefore warn at import-time (non-fatal) and error at runtime
(fatal if missing).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Union, List, Tuple

import os
import re
import shutil
import subprocess
import warnings

import numpy as np


# ---------------------------------------------------------------------
# Fortran number parsing helpers
# ---------------------------------------------------------------------

_FORTRAN_EXP_RE = re.compile(
    r"""^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d{2,4})$""", re.VERBOSE
)


def _to_float(tok: Any) -> float:
    """
    Parse a token emitted by Fortran into a Python float.

    Handles:
    - native numeric types
    - Fortran D exponents (1.0D-10)
    - rare "mantissa-EXP" tokens without E (0.12-322 -> 0.12e-322)
    """
    if isinstance(tok, (int, float, np.integer, np.floating)):
        return float(tok)

    t = str(tok).strip()
    t = t.replace("D", "E").replace("d", "e")

    m = _FORTRAN_EXP_RE.match(t)
    if m:
        t = f"{m.group(1)}e{m.group(2)}"

    return float(t)


def _fmt(x: Union[float, int, bool]) -> str:
    """Format numeric scalars robustly for Fortran stdin."""
    if isinstance(x, bool):
        return "1" if x else "0"
    if isinstance(x, int):
        return str(x)
    return format(_to_float(x), ".17g")  # round-trip safe for IEEE-754 double


def _potential_code(potential: Any) -> int:
    """
    Fortran prompt: Kepler (1) or Logarithmic (2).

    Accepts:
      - int already in {1,2}
      - string: "kepler"/"logarithmic"/"log"
      - callable returning int
      - object with .ipot/.code/.fortran_id
    """
    if isinstance(potential, int):
        return int(potential)

    if isinstance(potential, str):
        key = potential.strip().lower()
        mapping = {"kepler": 1, "k": 1, "logarithmic": 2, "log": 2}
        if key not in mapping:
            raise ValueError(
                f"Unknown potential='{potential}'. "
                "Use 'kepler'/'logarithmic' or an int 1/2."
            )
        return mapping[key]

    for attr in ("ipot", "code", "fortran_id"):
        if hasattr(potential, attr):
            return int(getattr(potential, attr))

    if callable(potential):
        v = potential()
        if isinstance(v, (int, np.integer)):
            return int(v)

    raise TypeError(
        "Could not interpret 'potential'. "
        "Provide int 1/2, string, callable->int, "
        "or an object with ipot/code/fortran_id."
    )


# ---------------------------------------------------------------------
# Backend resolution / compilation helpers
# ---------------------------------------------------------------------


class ScaleFreeBackendError(RuntimeError):
    """Raised when the Fortran backend cannot be found or built."""


def _have_gfortran() -> bool:
    return shutil.which("gfortran") is not None


def _site_packages_root() -> Path:
    # scalefree/vmoments.py -> scalefree/ -> (site-packages or repo root)
    return Path(__file__).resolve().parent.parent


def _packaged_fortran_source() -> Path:
    """
    Locate the packaged Fortran source.

    In wheels, this should be installed as:
      <site-packages>/fortran_src/scalefree.f

    In editable/repo contexts, this typically is:
      <repo-root>/fortran_src/scalefree.f

    Both are covered because _site_packages_root()
    is the parent of 'scalefree/'.
    """
    return _site_packages_root() / "fortran_src" / "scalefree.f"


def _user_cache_dir() -> Path:
    """
    Cross-platform-ish cache directory without extra dependencies.
    - Linux: $XDG_CACHE_HOME/scalefree or ~/.cache/scalefree
    - macOS: ~/Library/Caches/scalefree
    - Windows: %LOCALAPPDATA%\\scalefree\\Cache
    """
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or str(
            Path.home() / "AppData" / "Local",
        )
        return Path(base) / "scalefree" / "Cache"

    # macOS
    sysname = ""
    if hasattr(os, "uname"):
        try:
            sysname = os.uname().sysname.lower()
        except Exception:
            sysname = ""
    if sysname == "darwin":
        return Path.home() / "Library" / "Caches" / "scalefree"

    # Linux/other POSIX
    xdg = os.environ.get("XDG_CACHE_HOME")
    if xdg:
        return Path(xdg) / "scalefree"
    return Path.home() / ".cache" / "scalefree"


def _default_cached_exe() -> Path:
    exe_name = "scalefree.e" if os.name != "nt" else "scalefree.exe"
    return _user_cache_dir() / exe_name


def _compile_backend(*, exe: Path, src: Path) -> None:
    """
    Compile the Fortran backend. Raises ScaleFreeBackendError on failure.
    """
    if not _have_gfortran():
        raise ScaleFreeBackendError(_missing_gfortran_message())

    if not src.exists():
        raise ScaleFreeBackendError(
            "ScaleFree Fortran source file was not "
            "found inside the installation.\n\n"
            f"Expected: {src}\n"
            "This likely means the wheel/sdist was built without including "
            "fortran_src/scalefree.f."
        )

    exe.parent.mkdir(parents=True, exist_ok=True)

    cmd = ["gfortran", "-O2", "-std=legacy", "-o", str(exe), str(src)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise ScaleFreeBackendError(
            "gfortran was found, "
            "but compilation of the ScaleFree backend failed.\n\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"STDOUT:\n{e.stdout}\n\nSTDERR:\n{e.stderr}\n"
        ) from e

    # Ensure executable bit on POSIX
    if os.name != "nt":
        try:
            mode = exe.stat().st_mode
            exe.chmod(mode | 0o111)
        except Exception:
            pass


def _missing_gfortran_message() -> str:
    return (
        "ScaleFree Fortran backend is not available.\n\n"
        "This package requires a Fortran compiler "
        "(gfortran) to build the backend "
        "executable.\n"
        "Please install gfortran and re-run.\n\n"
        "Typical install commands:\n"
        "  Debian/Ubuntu:  sudo apt-get install gfortran\n"
        "  Fedora:         sudo dnf install gcc-gfortran\n"
        "  macOS (brew):   brew install gcc    # provides gfortran\n"
        "  Windows:        use WSL or MSYS2 to install gfortran\n\n"
        "Alternatively, if you already built the executable elsewhere, set:\n"
        "  export SCALEFREE_EXE=/path/to/scalefree.e\n"
        "or pass exe_path=... explicitly."
    )


def _resolve_executable(
    exe_path: Optional[Union[str, Path]],
    workdir: Optional[Union[str, Path]],
) -> Tuple[Path, Path]:
    """
    Resolve and (if needed) build the backend executable.

    Returns (exe, resolved_workdir).
    """
    # Workdir: if caller provides, respect it; else pick a safe default
    if workdir is not None:
        wd = Path(workdir).expanduser().resolve()
        wd.mkdir(parents=True, exist_ok=True)
    else:
        wd = _user_cache_dir()
        wd.mkdir(parents=True, exist_ok=True)

    # 1) explicit exe_path
    if exe_path is not None:
        exe = Path(exe_path).expanduser().resolve()
        if exe.exists():
            return exe, wd
        # If they provided an exe_path that doesn't exist,
        # build there if possible
        src = _packaged_fortran_source()
        _compile_backend(exe=exe, src=src)
        return exe, wd

    # 2) env var
    env = os.environ.get("SCALEFREE_EXE")
    if env:
        exe = Path(env).expanduser().resolve()
        if exe.exists():
            return exe, wd
        raise ScaleFreeBackendError(
            f"SCALEFREE_EXE is set but does not exist: {exe}\n"
            "Either unset SCALEFREE_EXE "
            "or point it to a valid compiled executable."
        )

    # 3) cached exe
    cached = _default_cached_exe()
    if cached.exists():
        return cached, wd

    # 4) build into cache if gfortran exists; else fail with clear instructions
    src = _packaged_fortran_source()
    _compile_backend(exe=cached, src=src)
    return cached, wd


# Best-effort import-time warning (non-fatal).
try:
    _cached = _default_cached_exe()
    _src = _packaged_fortran_source()
    if (not _cached.exists()) and _src.exists() and (not _have_gfortran()):
        warnings.warn(
            "scalefree: gfortran not found. "
            "The Fortran backend will not be usable "
            "until you install gfortran. "
            "See scalefree.vmoments.ScaleFreeBackendError "
            "for install instructions.",
            RuntimeWarning,
            stacklevel=2,
        )
except Exception:
    pass


# ---------------------------------------------------------------------
# Output container
# ---------------------------------------------------------------------


@dataclass
class ScaleFreeResult:
    blocks: Dict[str, Any]
    raw_text: str
    output_path: Optional[Path]
    stdout: str
    stderr: str


# ---------------------------------------------------------------------
# Parser for structured output
# ---------------------------------------------------------------------


def parse_scalefree_output(text: str) -> Dict[str, Any]:
    """
    Parse the structured ASCII output produced by the modified Fortran code.

    Recognizes:
      - "# kind=XYZ" blocks with optional "# columns: ..." line
      - "# kind=vp" and "# kind=vp_intrinsic" blocks (with embedded vp_table sub-blocks)
      - "# vp_table iproj X" and "# vp_table icomp X" blocks with optional "# columns: ..." line
    """
    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    blocks: Dict[str, Any] = {}
    i = 0

    # The modified Fortran sometimes emits multi-line column headers.
    # Example (vp block):
    #   # columns: iproj true_gam true_V true_sig
    #   #          gauss_gam gauss_V gauss_sig
    #   #         h0 h1 h2 h3 h4 h5 h6
    #
    # The upstream file may only parse the first line. We extend it to
    # consume continuation lines that contain additional identifier tokens.
    _COLNAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

    def parse_columns(ln: str) -> List[str]:
        _, rhs = ln.split(":", 1)
        return rhs.strip().split()

    def parse_columns_continuation(ln: str) -> List[str]:
        s = ln.strip()
        if not s.startswith("#"):
            return []
        # Strip leading '#', keep the remainder.
        tail = s[1:].strip()
        # If it looks like a new header or marker, do not treat as continuation.
        if (not tail) or (":" in tail) or tail.startswith("kind=") or tail.startswith("vp_table"):
            return []
        toks = tail.split()
        if toks and all(_COLNAME_RE.match(t) for t in toks):
            return toks
        return []

    while i < len(lines):
        s = lines[i].strip()
        if not (s.startswith("# kind=") or s.startswith("# vp_table")):
            i += 1
            continue

        # vp_table iproj X
        if s.startswith("# vp_table"):
            parts = s.split()
            iproj = int(parts[-1])
            i += 1
            cols = None
            while i < len(lines) and lines[i].strip().startswith("#"):
                if lines[i].strip().startswith("# columns:"):
                    cols = parse_columns(lines[i].strip())
                elif cols is not None:
                    cols += parse_columns_continuation(lines[i])
                i += 1

            data = []
            while i < len(lines):
                row = lines[i].strip()
                if row == "" or row.startswith("#"):
                    break
                data.append([_to_float(x) for x in row.split()])
                i += 1
            label = parts[2].lower() if len(parts) >= 3 else "iproj"
            vp_table_key = "vp_table_intrinsic" if label == "icomp" else "vp_table"
            blocks.setdefault(vp_table_key, {})[iproj] = {
                "columns": cols if cols else ["v", "vp"],
                "data": np.array(data, dtype=float),
            }
            continue

        # kind=...
        kind = s.replace("# kind=", "").strip()
        i += 1
        cols = None
        while i < len(lines) and lines[i].strip().startswith("#"):
            if lines[i].strip().startswith("# columns:"):
                cols = parse_columns(lines[i].strip())
            elif cols is not None:
                cols += parse_columns_continuation(lines[i])
            i += 1

        # Special handling for the VP summary block: the modified Fortran output
        # interleaves additional VP summary rows with vp_table sub-blocks.
        # Only the first row is preceded by '# kind=vp'; subsequent rows (iproj=2,3,...)
        # appear later as numeric lines with many columns (while vp_table rows have only 2).
        if kind in ("vp", "vp_intrinsic"):
            data = []
            while i < len(lines):
                row = lines[i].strip()

                # End of vp block when a new '# kind=' begins (different kind)
                if row.startswith("# kind=") and not row.startswith(f"# kind={kind}"):
                    break

                # Skip blanks
                if row == "":
                    i += 1
                    continue

                # Parse embedded vp_table blocks
                if row.startswith("# vp_table"):
                    parts = row.split()
                    iproj = int(parts[-1])
                    i += 1
                    tcols = None
                    while i < len(lines) and lines[i].strip().startswith("#"):
                        if lines[i].strip().startswith("# columns:"):
                            tcols = parse_columns(lines[i].strip())
                        elif tcols is not None:
                            tcols += parse_columns_continuation(lines[i])
                        i += 1

                    tdata = []
                    while i < len(lines):
                        r2 = lines[i].strip()
                        if r2 == "" or r2.startswith("#"):
                            break
                        tdata.append([_to_float(x) for x in r2.split()])
                        i += 1

                    vp_table_key = "vp_table" if kind == "vp" else "vp_table_intrinsic"
                    blocks.setdefault(vp_table_key, {})[iproj] = {
                        "columns": tcols if tcols else ["v", "vp"],
                        "data": np.array(tdata, dtype=float),
                    }
                    continue

                # Skip other comment lines
                if row.startswith("#"):
                    i += 1
                    continue

                toks = row.split()
                # vp summary rows are wide (>2); vp_table rows are 2 columns and are
                # handled above after their '# vp_table' marker.
                if len(toks) > 2:
                    data.append([_to_float(x) for x in toks])
                i += 1

            arr = (
                np.array(data, dtype=float)
                if data
                else np.empty((0, 0), dtype=float)
            )

            # Defensive: ensure the header width matches the numeric table width.
            if cols is None:
                cols_norm: List[str] = []
            else:
                cols_norm = list(cols)

            if arr.ndim == 2 and arr.size > 0:
                n_data_cols = int(arr.shape[1])
                if len(cols_norm) > n_data_cols:
                    cols_norm = cols_norm[:n_data_cols]
                elif len(cols_norm) < n_data_cols:
                    cols_norm = cols_norm + [
                        f"col{j+1}" for j in range(len(cols_norm), n_data_cols)
                    ]

            block: Dict[str, Any] = {"columns": cols_norm, "data": arr}

            if cols_norm and arr.shape[0] > 0:
                first = cols_norm[0].lower()
                if first in ("iproj", "icomp"):
                    by_id: Dict[int, Dict[str, float]] = {}
                    for r in arr:
                        idx = int(r[0])
                        by_id[idx] = {
                            cols_norm[j]: r[j]
                            for j in range(min(len(cols_norm), len(r)))
                        }
                    if first == "iproj":
                        block["by_iproj"] = by_id
                    else:
                        block["by_icomp"] = by_id

            blocks[kind] = block
            continue

        data = []
        while i < len(lines):
            row = lines[i].strip()
            if row == "" or row.startswith("#"):
                break
            data.append([_to_float(x) for x in row.split()])
            i += 1

        arr = (
            np.array(data, dtype=float)
            if data
            else np.empty(
                (0, 0),
                dtype=float,
            )
        )

        # Defensive: ensure the header width matches the numeric table width.
        # The modified Fortran code has emitted different VP header variants
        # across versions (e.g. with/without h6). We retain the parsed header
        # tokens, but normalize to the actual data width so downstream code can
        # safely align columns with values.
        if cols is None:
            cols_norm: List[str] = []
        else:
            cols_norm = list(cols)

        if arr.ndim == 2 and arr.size > 0:
            n_data_cols = int(arr.shape[1])
            if len(cols_norm) > n_data_cols:
                cols_norm = cols_norm[:n_data_cols]
            elif len(cols_norm) < n_data_cols:
                cols_norm = cols_norm + [f"col{j+1}" for j in range(len(cols_norm), n_data_cols)]

        block: Dict[str, Any] = {"columns": cols_norm, "data": arr}

        # Convenience indexing for tables where first column is iproj
        if cols_norm and cols_norm[0].lower() == "iproj" and arr.shape[0] > 0:
            by_iproj: Dict[int, Dict[str, float]] = {}
            for r in arr:
                ip = int(r[0])
                by_iproj[ip] = {
                    cols_norm[j]: r[j]
                    for j in range(
                        min(len(cols_norm), len(r)),
                    )
                }
            block["by_iproj"] = by_iproj

        blocks[kind] = block

    return blocks


# ---------------------------------------------------------------------
# STDOUT structured extraction (Option A)
# ---------------------------------------------------------------------


_NUMERIC_ROW_RE = re.compile(
    r"""^\s*[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:(?:[EeDd])[+-]?\d+)?(?:\s+|$)"""
)


def _extract_structured_from_stdout(stdout_text: str) -> str:
    """
    Extract only the structured blocks from STDOUT.

    We keep:
      - lines starting with '#'
      - numeric rows ONLY when we are "inside" a block (after a '# kind=' or
        '# vp_table' marker), stopping when we hit a non-numeric, non-# line.
    """
    keep: List[str] = []
    in_block = False

    for ln in stdout_text.splitlines():
        s = ln.strip()
        if not s:
            # do not force-close blocks on blank lines; just skip
            continue

        if s.startswith("#"):
            keep.append(ln)
            in_block = (
                s.startswith(
                    "# kind=",
                )
                or s.startswith("# vp_table")
                or in_block
            )
            continue

        if in_block and _NUMERIC_ROW_RE.match(ln):
            keep.append(ln)
            continue

        # Any other non-comment, non-numeric line breaks a block context
        in_block = False

    return "\n".join(keep)


# ---------------------------------------------------------------------
# Intrinsic moments extraction (text-only in current Fortran)
# ---------------------------------------------------------------------

_INTRINSIC_SHELL_MARK = "Mass-weighted average spherical shell"
_INTRINSIC_HDR_RE = re.compile(r"\brho\b.*<v_ph>.*<v_r\^2>.*<v_th\^2>.*<v_ph\^2>", re.IGNORECASE)
_FLOAT_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?|[+-]?(?:\d+(?:\.\d*)?|\.\d+)[+-]\d{2,4}")


def _extract_last_intrinsic_block_from_stdout(stdout_text: str) -> str:
    """
    ScaleFree's intrinsic results are printed as human-readable text (WRITE(*))
    rather than '# kind=...' blocks. This helper parses the *last* intrinsic
    moment table in STDOUT and re-encodes it as a synthetic structured block
    so downstream parsing remains uniform.

    Returns an empty string if no intrinsic table is detected.
    """
    lines = stdout_text.splitlines()
    # We scan for the last occurrence of an intrinsic table header, then parse
    # the numeric row that follows.
    last_idx = None
    shell_avg = False

    for i, ln in enumerate(lines):
        if "Intrinsic velocity moments" in ln:
            # reset; we will decide shell_avg when we see the next marker
            last_idx = i
            shell_avg = False
        if last_idx is not None and _INTRINSIC_SHELL_MARK in ln:
            shell_avg = True

    if last_idx is None:
        return ""

    # From last_idx onward, find the header line, then the numeric row.
    hdr_i = None
    for i in range(last_idx, len(lines)):
        if _INTRINSIC_HDR_RE.search(lines[i]):
            hdr_i = i
            break
    if hdr_i is None:
        return ""

    # Find first non-empty line after header that contains at least 5 floats.
    data_i = None
    for i in range(hdr_i + 1, min(hdr_i + 6, len(lines))):
        s = lines[i].strip()
        if not s:
            continue
        toks = _FLOAT_RE.findall(s.replace("D", "E").replace("d", "e"))
        if len(toks) >= 5:
            data_i = i
            break
    if data_i is None:
        return ""

    toks = _FLOAT_RE.findall(lines[data_i].replace("D", "E").replace("d", "e"))
    vals = [_to_float(t) for t in toks]

    if shell_avg:
        # Expected: rho, <v_ph>, <v_r^2>, <v_th^2>, <v_ph^2>, beta
        if len(vals) < 6:
            return ""
        cols = ["rho", "vphi", "vr2", "vth2", "vphi2", "beta"]
        kind = "intrinsic_shell_average"
        vals = vals[:6]
    else:
        # Expected: rho, <v_ph>, <v_r^2>, <v_th^2>, <v_ph^2>
        cols = ["rho", "vphi", "vr2", "vth2", "vphi2"]
        kind = "intrinsic_point"
        vals = vals[:5]

    # Synthetic structured block understood by parse_scalefree_output()
    out = [f"# kind={kind}", "# columns: " + " ".join(cols), "  " + "  ".join(_fmt(v) for v in vals)]
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------
# Prompt-driven runner
# ---------------------------------------------------------------------


class ScaleFreeRunner:
    """
    Runs the ScaleFree Fortran executable and parses structured output.

    Defaults to "no persistent files":
    - if output_path is None, we still answer the Fortran
    "Output file" prompt with a
      temporary filename, but we delete it after the run and
      return output_path=None.
    """

    def __init__(
        self,
        exe_path: Optional[Union[str, Path]] = None,
        workdir: Optional[Union[str, Path]] = None,
    ):
        exe, wd = _resolve_executable(exe_path, workdir)
        self.exe_path = exe
        self.workdir = wd

        if not self.exe_path.exists():
            raise FileNotFoundError(f"Executable not found: {self.exe_path}")

    def vprofile(
        self,
        *,
        potential: Any,
        gamma: float,
        q: float,
        beta: float,
        s: float,
        t: float,
        inclination: float,
        xi: float,
        theta: float,
        df: int = 1,
        integration: int = 1,
        # 0 Romberg, 1 Gauss-Legendre
        ngl_or_eps: float = 0.0,
        # eps if Romberg; nGL if Gauss-Legendre (0 -> default)
        algorithm: int = 1,
        # VP algorithm (1 default)
        vp_reg_param: float = 1.0,
        # Only used when algorithm==2 (fixed regularization strength)
        vp_smooth_eps: float = 0.0,
        # Only used when algorithm==3 (VP smoothness factor; 0.0 => Fortran default)
        maxmom: int = 4,
        average: bool = False,
        kinematics: Union[str, int] = "projected",
        usevp: bool = False,
        verbose_vp: int = 0,
        output_path: Optional[Union[str, Path]] = None,
        timeout_s: int = 120,
        parse_stdout_fallback: bool = False,
        # legacy flag; kept for compatibility
        debug_prompts: bool = False,
    ) -> ScaleFreeResult:
        ipot = _potential_code(potential)

        # ---------------------------------------------------------
        # Output handling
        # ---------------------------------------------------------
        persist_file = output_path is not None
        delete_after = False

        if output_path is None:
            # Create a temp filename to satisfy Fortran prompt,
            # but delete afterwards.
            outname = (
                "scalefree_"
                + "tmp_"
                + f"{os.getpid()}_{id(self) % 10_000_000}"
                + ".txt"
            )
            delete_after = True
        else:
            outname = str(output_path)
            outname = (
                Path(
                    outname,
                ).name
                if Path(outname).is_absolute()
                else outname
            )

        out_path = self.workdir / outname

        # ---------------------------------------------------------
        # Prompt answers
        # ---------------------------------------------------------
        # NOTE: Fortran asks extra questions when algorithm != 1.
        # We keep the public API stable by providing sensible defaults
        # that reproduce the interactive "default" behavior:
        #  - algorithm==2: regularization parameter (>0), default 1.0
        #  - algorithm==3: smoothness eps (0.0 yields default), default 0.0
        _vp_reg_param = float(vp_reg_param) if int(algorithm) == 2 else 1.0
        if _vp_reg_param <= 0.0:
            # Fortran requires >0. Keep it minimally safe.
            _vp_reg_param = 1.0
        _vp_smooth_eps = float(vp_smooth_eps) if int(algorithm) == 3 else 0.0

        answers: Dict[str, str] = {
            "Kepler (1) or Logarithmic (2)": str(ipot),
            "Power-law slope gamma": _fmt(gamma),
            "Intrinsic axial ratio q": _fmt(q),
            "Case I (1) or Case II (2) DF": str(int(df)),
            "Case I (1) or Case II (2)": str(int(df)),
            "Anisotropy parameter beta": _fmt(beta),
            "Odd part parameters s and t": f"{_fmt(s)} {_fmt(t)}",
            "Viewing inclination i": _fmt(inclination),
            "Use Romberg (0) or Gauss-Legendre (1)": str(int(integration)),
            "Give the fractional accuracy epsilon": _fmt(ngl_or_eps),
            "Give number of quadrature points": (
                str(int(ngl_or_eps)) if float(ngl_or_eps).is_integer() else "0"
            ),
            "Choose 1 for default.": str(int(algorithm)),
            # Extra prompts for VP algorithms:
            "Give the regularization parameter": _fmt(_vp_reg_param),
            "Give smoothness factor eps": _fmt(_vp_smooth_eps),
            "Give the maximum number of projected moments": str(int(maxmom)),
            "Give the number of projected moments": str(int(maxmom)),
            # Always answer output file prompt to avoid blocking.
            "Output file": outname,
        }

        # Choose iwhat values based on average flag and requested kinematics
        def _norm_kinematics(k: Union[str, int]) -> Tuple[str, Optional[int]]:
            """
            Returns (mode, iwhat_single)
              - mode == "single": run exactly one calculation with iwhat_single in {0,1,2,3}
              - mode == "both"  : run intrinsic then projected (legacy behavior)
            """
            if isinstance(k, (int, np.integer)):
                iw = int(k)
                if iw not in (0, 1, 2, 3):
                    raise ValueError("kinematics int must be one of {0,1,2,3}.")
                return ("single", iw)

            key = str(k).strip().lower()
            if key in ("projected", "proj", "p"):
                return ("single", 3 if average else 1)
            if key in ("intrinsic", "intr", "i"):
                return ("single", 2 if average else 0)
            if key in ("both", "all"):
                return ("both", None)

            raise ValueError(
                "kinematics must be 'projected', 'intrinsic', 'both', or an int in {0,1,2,3}."
            )

        mode, iwhat_single = _norm_kinematics(kinematics)

        if average:
            iwhat_intr = 2
            iwhat_proj = 3
        else:
            iwhat_intr = 0
            iwhat_proj = 1

        phase = {"step": 0}  # 0 intrinsic, 1 projected (only when mode == "both")

        def respond(line: str) -> Optional[str]:
            for key, val in answers.items():
                if key in line:
                    return val

            # iwhat prompt
            if "Calculate intrinsic (0) or projected (1)" in line:
                if mode == "both":
                    return str(iwhat_intr) if phase["step"] == 0 else str(iwhat_proj)
                assert iwhat_single is not None
                return str(iwhat_single)

            # theta prompt (intrinsic point only)
            if "Give angle theta in the meridional plane" in line:
                return _fmt(theta)

            # xi prompt (projected only)
            if "Give angle on the projected plane" in line:
                return _fmt(xi)

            # verbose prompt (only for projected modes)
            if "Give verbose output of intermediate steps" in line:
                return str(int(verbose_vp))

            # VP prompt variants
            if (
                ("Calculate VPs" in line)
                or ("Use VPs" in line)
                or ("VP" in line and "?" in line)
            ):
                return "1" if usevp else "0"

            # continue? ("Calculate something else for this model?")
            if "Calculate something else for this model" in line:
                if mode == "both" and phase["step"] == 0:
                    phase["step"] = 1
                    return "1"
                return "0"

            return None
# ---------------------------------------------------------
        # Run Fortran interactively
        # ---------------------------------------------------------
        p = subprocess.Popen(
            [str(self.exe_path)],
            cwd=str(self.workdir),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        stdout_lines: List[str] = []
        stderr_text = ""

        try:
            assert p.stdout is not None and p.stdin is not None

            for line in p.stdout:
                stdout_lines.append(line)
                if debug_prompts:
                    print(line, end="")

                ans = respond(line)
                if ans is not None:
                    p.stdin.write(ans + "\n")
                    p.stdin.flush()

            stderr_text = p.stderr.read() if p.stderr else ""
            rc = p.wait(timeout=timeout_s)

        except subprocess.TimeoutExpired:
            p.kill()
            raise RuntimeError(f"Fortran run timed out after {timeout_s}s.")

        stdout_text = "".join(stdout_lines)

        try:
            if rc != 0 or "STOP Wrong answer" in stderr_text:
                raise RuntimeError(
                    "Fortran execution failed.\n\n"
                    f"Return code: {rc}\n"
                    f"STDERR:\n{stderr_text}\n\n"
                    f"STDOUT (first 2000 chars):\n{stdout_text[:2000]}\n"
                )

            def _finalize_blocks(blocks_in: Dict[str, Any], raw_in: str) -> Tuple[Dict[str, Any], str]:
                """
                Post-process parsed blocks according to requested kinematics:
                  - projected-only: keep projected blocks
                  - intrinsic-only: return only intrinsic blocks parsed from STDOUT text tables
                  - both: merge intrinsic blocks (from STDOUT) with projected blocks (from structured output)
                """
                want_intr = False
                want_proj = False
                if mode == "both":
                    want_intr = True
                    want_proj = True
                else:
                    assert iwhat_single is not None
                    want_intr = iwhat_single in (0, 2)
                    want_proj = iwhat_single in (1, 3)

                blocks = dict(blocks_in) if blocks_in else {}
                raw = raw_in or ""

                intrinsic_raw = ""
                if want_intr:
                    if ('intrinsic_point' not in blocks) and ('intrinsic_shell_average' not in blocks):
                        intrinsic_raw = _extract_last_intrinsic_block_from_stdout(stdout_text)
                    else:
                        intrinsic_raw = ""

                    if intrinsic_raw.strip():
                        iblocks = parse_scalefree_output(intrinsic_raw)
                        # Should be exactly one kind, but we merge defensively
                        for k, v in iblocks.items():
                            blocks[k] = v
                        # Preserve provenance in raw_text for debugging/regression tests
                        if raw and not raw.endswith("\n"):
                            raw += "\n"
                        raw += intrinsic_raw

                if not want_proj:
                    # Intrinsic-only: strip projected outputs if present
                    blocks = {
                        k: v
                        for k, v in blocks.items()
                        if k.startswith("intrinsic") or k in ("vp_intrinsic", "vp_table_intrinsic")
                    }
                    raw = intrinsic_raw if intrinsic_raw.strip() else raw

                if not want_intr:
                    # Projected-only: strip intrinsic outputs if present
                    blocks = {
                        k: v
                        for k, v in blocks.items()
                        if (not k.startswith("intrinsic")) and k not in ("vp_intrinsic", "vp_table_intrinsic")
                    }

                if not blocks:
                    raise RuntimeError(
                        "Fortran returned success but no parsable output was detected. "
                        "If you requested intrinsic quantities, note that intrinsic tables "
                        "are printed to STDOUT and must match the expected header format."
                    )

                return blocks, raw

            # ---------------------------------------------------------
            # Parse output: stable preference depends on whether
            # caller requested a file
            # ---------------------------------------------------------

            # (A) If caller explicitly requested a file,
            # prefer parsing that file
            if persist_file and out_path.exists():
                raw = out_path.read_text(encoding="utf-8", errors="replace")
                if debug_prompts:
                    print("\n--- Fortran structured output (from file) ---")
                    print(raw, end="" if raw.endswith("\n") else "\n")
                blocks = parse_scalefree_output(raw)
                blocks, raw = _finalize_blocks(blocks, raw)
                return ScaleFreeResult(
                    blocks=blocks,
                    raw_text=raw,
                    output_path=out_path,
                    stdout=stdout_text,
                    stderr=stderr_text,
                )

            # (B) Parse structured output from the file first (authoritative)
            raw_file = None
            blocks_file: Dict[str, Any] = {}
            if out_path.exists():
                raw_file = out_path.read_text(encoding="utf-8", errors="replace")
                if debug_prompts:
                    print("\n--- Fortran structured output (from file) ---")
                    print(raw_file, end="" if raw_file.endswith("\n") else "\n")
                blocks_file = parse_scalefree_output(raw_file)

            # (C) Also parse any structured blocks present in STDOUT (may be partial)
            raw_stdout_struct = _extract_structured_from_stdout(stdout_text)
            blocks_stdout: Dict[str, Any] = {}
            if raw_stdout_struct.strip():
                if debug_prompts:
                    print("\n--- Fortran structured output (from stdout structured blocks) ---")
                    print(raw_stdout_struct, end="" if raw_stdout_struct.endswith("\n") else "\n")
                blocks_stdout = parse_scalefree_output(raw_stdout_struct)

            # Merge (prefer file blocks; fill gaps with stdout blocks)
            blocks = dict(blocks_file)
            for k, v in blocks_stdout.items():
                blocks.setdefault(k, v)

            base_raw = raw_file if raw_file is not None else raw_stdout_struct
            blocks, base_raw = _finalize_blocks(blocks, base_raw)

            if blocks:
                return ScaleFreeResult(
                    blocks=blocks,
                    raw_text=base_raw,
                    output_path=out_path if persist_file else None,
                    stdout=stdout_text,
                    stderr=stderr_text,
                )

            # (D) Legacy fallback flag (kept for backward compatibility)
            if parse_stdout_fallback and raw_stdout_struct.strip():
                blocks = parse_scalefree_output(raw_stdout_struct)
                blocks, raw_stdout_struct = _finalize_blocks(blocks, raw_stdout_struct)
                return ScaleFreeResult(
                    blocks=blocks,
                    raw_text=raw_stdout_struct,
                    output_path=None,
                    stdout=stdout_text,
                    stderr=stderr_text,
                )

            raise RuntimeError(
                "Fortran returned success but no structured"
                "output was detected.\n"
                "Expected either structured '# kind=...' "
                "blocks in STDOUT or an output file.\n"
                f"STDOUT (first 2000 chars):\n{stdout_text[:2000]}\n"
            )

        finally:
            # Clean up temporary file if we created it
            # purely to satisfy Fortran prompts
            if delete_after:
                try:
                    if out_path.exists():
                        out_path.unlink()
                except Exception:
                    # Non-fatal cleanup failure
                    pass


# ---------------------------------------------------------------------
# Simple functional API (what most users will call)
# ---------------------------------------------------------------------


def vprofile(
    *,
    exe_path: Optional[Union[str, Path]] = None,
    potential: Any,
    gamma: float,
    q: float,
    beta: float,
    s: float,
    t: float,
    inclination: float,
    xi: float,
    theta: float,
    df: int = 1,
    integration: int = 1,
    ngl_or_eps: float = 0.0,
    algorithm: int = 1,
    maxmom: int = 4,
    average: bool = False,
    kinematics: Union[str, int] = "projected",
    usevp: bool = False,
    verbose_vp: int = 0,
    output_path: Optional[Union[str, Path]] = None,
    timeout_s: int = 120,
    parse_stdout_fallback: bool = False,
    debug_prompts: bool = False,
    workdir: Optional[Union[str, Path]] = None,
) -> ScaleFreeResult:
    """
    Convenience function that instantiates a runner and executes vprofile.

    Users can omit exe_path; the backend will be resolved/built automatically.

    Output behavior:
    - output_path=None (default): do not leave output files behind
    - output_path="something.txt": write and keep that file
    (useful for regression)
    """
    runner = ScaleFreeRunner(exe_path=exe_path, workdir=workdir)
    return runner.vprofile(
        potential=potential,
        gamma=gamma,
        q=q,
        beta=beta,
        s=s,
        t=t,
        inclination=inclination,
        xi=xi,
        theta=theta,
        df=df,
        integration=integration,
        ngl_or_eps=ngl_or_eps,
        algorithm=algorithm,
        maxmom=maxmom,
        average=average,
        kinematics=kinematics,
        usevp=usevp,
        verbose_vp=verbose_vp,
        output_path=output_path,
        timeout_s=timeout_s,
        parse_stdout_fallback=parse_stdout_fallback,
        debug_prompts=debug_prompts,
    )
