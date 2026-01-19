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
      - "# vp_table iproj X" blocks with optional "# columns: ..." line
    """
    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    blocks: Dict[str, Any] = {}
    i = 0

    def parse_columns(ln: str) -> List[str]:
        _, rhs = ln.split(":", 1)
        return rhs.strip().split()

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
                i += 1

            data = []
            while i < len(lines):
                row = lines[i].strip()
                if row == "" or row.startswith("#"):
                    break
                data.append([_to_float(x) for x in row.split()])
                i += 1

            blocks.setdefault("vp_table", {})[iproj] = {
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
            i += 1

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
        block: Dict[str, Any] = {"columns": cols if cols else [], "data": arr}

        # Convenience indexing for tables where first column is iproj
        if cols and cols[0].lower() == "iproj" and arr.shape[0] > 0:
            by_iproj: Dict[int, Dict[str, float]] = {}
            for r in arr:
                ip = int(r[0])
                by_iproj[ip] = {
                    cols[j]: r[j]
                    for j in range(
                        min(len(cols), len(r)),
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
        maxmom: int = 4,
        average: bool = False,
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
            "Give the maximum number of projected moments": str(int(maxmom)),
            "Give the number of projected moments": str(int(maxmom)),
            # Always answer output file prompt to avoid blocking.
            "Output file": outname,
        }

        # Choose iwhat values based on average flag
        if average:
            iwhat_intr = 2
            iwhat_proj = 3
        else:
            iwhat_intr = 0
            iwhat_proj = 1

        phase = {"step": 0}  # 0 intrinsic, 1 projected

        def respond(line: str) -> Optional[str]:
            for key, val in answers.items():
                if key in line:
                    return val

            # iwhat prompt
            if "Calculate intrinsic (0) or projected (1)" in line:
                return (
                    str(iwhat_intr)
                    if phase["step"] == 0
                    else str(
                        iwhat_proj,
                    )
                )

            # theta prompt (intrinsic)
            if "Give angle theta in the meridional plane" in line:
                return _fmt(theta)

            # xi prompt (projected)
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
                if phase["step"] == 0:
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

            # ---------------------------------------------------------
            # Parse output: stable preference depends on whether
            # caller requested a file
            # ---------------------------------------------------------

            # (A) If caller explicitly requested a file,
            # prefer parsing that file
            if persist_file and out_path.exists():
                raw = out_path.read_text(encoding="utf-8", errors="replace")
                blocks = parse_scalefree_output(raw)
                return ScaleFreeResult(
                    blocks=blocks,
                    raw_text=raw,
                    output_path=out_path,
                    stdout=stdout_text,
                    stderr=stderr_text,
                )

            # (B) Otherwise, prefer structured blocks from stdout (Option A)
            raw_stdout_struct = _extract_structured_from_stdout(stdout_text)
            if raw_stdout_struct.strip():
                blocks = parse_scalefree_output(raw_stdout_struct)
                if blocks:
                    return ScaleFreeResult(
                        blocks=blocks,
                        raw_text=raw_stdout_struct,
                        output_path=(
                            None
                            if delete_after
                            else (out_path if persist_file else None)
                        ),
                        stdout=stdout_text,
                        stderr=stderr_text,
                    )

            # (C) Fallback: if a file exists (temp or explicit), parse it
            if out_path.exists():
                raw = out_path.read_text(encoding="utf-8", errors="replace")
                blocks = parse_scalefree_output(raw)
                return ScaleFreeResult(
                    blocks=blocks,
                    raw_text=raw,
                    output_path=out_path if persist_file else None,
                    stdout=stdout_text,
                    stderr=stderr_text,
                )

            # (D) Legacy fallback flag (kept, but typically redundant now)
            if parse_stdout_fallback:
                blocks = parse_scalefree_output(raw_stdout_struct)
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
        usevp=usevp,
        verbose_vp=verbose_vp,
        output_path=output_path,
        timeout_s=timeout_s,
        parse_stdout_fallback=parse_stdout_fallback,
        debug_prompts=debug_prompts,
    )
