"""
scalefree.vmoments

Prompt-driven Python wrapper for the ScaleFree Fortran executable.

Key properties
--------------
- Users pass numeric values (floats/ints/bools), not strings.
- We do NOT parse result printouts.
- We DO respond to Fortran prompts (interactive control)
to ensure correct stdin order.
- Results are read from the structured ASCII output file produced
by the modified Fortran code.
- Fallback: if the file is not produced, we can optionally
parse ONLY the structured
  '# kind=...' blocks from stdout (disabled by default).

This wrapper is robust to reordering of prompts
(including an "Output file" prompt),
because we match prompt substrings rather than
relying on a hard-coded input order.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Union, List

import subprocess
import numpy as np
import re


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
                f"Unknown potential='{potential}'."
                + " Use 'kepler'/'logarithmic' or an int 1/2."
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
        "Could not interpret 'potential'."
        + "Provide int 1/2,"
        + " string, "
        + "callable->int, "
        "or an object with ipot/code/fortran_id."
    )


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
# Parser for structured output file
# ---------------------------------------------------------------------


def parse_scalefree_output(text: str) -> Dict[str, Any]:
    """
    Parse the structured ASCII output produced by the modified Fortran code.

    Recognizes:
      - "# kind=XYZ" blocks with optional "# columns: ..." line
      - "# vp_table iproj X" blocks with optional "# columns: ..." line

    Returns a dict of blocks:
      blocks[kind] = {"columns": [...],
      "data": np.ndarray, "by_iproj": {...} (optional)}
      blocks["vp_table"][iproj] = {"columns": [...], "data": np.ndarray}
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
        block = {"columns": cols if cols else [], "data": arr}

        # Convenience indexing for tables where first column is iproj
        if cols and cols[0].lower() == "iproj" and arr.shape[0] > 0:
            by_iproj = {}
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
# Prompt-driven runner
# ---------------------------------------------------------------------


class ScaleFreeRunner:
    """
    Runs the ScaleFree Fortran executable and
    parses its structured output file.
    """

    def __init__(
        self,
        exe_path: Union[str, Path],
        workdir: Optional[Union[str, Path]] = None,
    ):
        self.exe_path = Path(exe_path).expanduser().resolve()
        if not self.exe_path.exists():
            raise FileNotFoundError(f"Executable not found: {self.exe_path}")
        self.workdir = (
            Path(
                workdir,
            )
            .expanduser()
            .resolve()
            if workdir
            else self.exe_path.parent
        )

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
        maxmom: int = 30,
        # max order/number of moments (depends on your Fortran prompts)
        average: bool = False,
        # if True, use iwhat=2 and 3
        usevp: bool = False,
        # currently only influences whether we answer certain prompts
        verbose_vp: int = 0,
        output_path: Optional[Union[str, Path]] = "out.txt",
        timeout_s: int = 120,
        parse_stdout_fallback: bool = False,
        debug_prompts: bool = False,
    ) -> ScaleFreeResult:
        """
        Drives the interactive Fortran executable via prompt matching.

        The run sequence is:
          - intrinsic (iwhat=0 or 2)
          - then projected (iwhat=1 or 3)
        in a single session, using the built-in
        "Calculate something else?" prompt.

        Returns ScaleFreeResult with parsed blocks.
        """
        ipot = _potential_code(potential)

        # Keep output filename short to avoid Fortran CHARACTER truncation
        if output_path is None:
            outname = "out.txt"
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

        # Map prompt fragments to responses
        answers = {
            "Kepler (1) or Logarithmic (2)": str(ipot),
            "Power-law slope gamma": _fmt(gamma),
            "Intrinsic axial ratio q": _fmt(q),
            "Case I (1) or Case II (2) DF": str(int(df)),
            "Case I (1) or Case II (2)": str(int(df)),
            # tolerate variant prompt
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
            "Output file": outname,
            # matches "Output file for results ..." variations
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
            # Direct prompt matches
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

            # Some codes explicitly ask whether to compute/use VPs;
            # if present, answer from usevp.
            if (
                "Calculate VPs" in line
                or "Use VPs" in line
                or "VP" in line
                and "?" in line
            ):
                return "1" if usevp else "0"

            # continue? ("Calculate something else for this model?")
            if "Calculate something else for this model" in line:
                if phase["step"] == 0:
                    phase["step"] = 1
                    return "1"
                return "0"

            return None

        # Run Fortran interactively
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

        if rc != 0 or "STOP Wrong answer" in stderr_text:
            raise RuntimeError(
                "Fortran execution failed.\n\n"
                f"Return code: {rc}\n"
                f"STDERR:\n{stderr_text}\n\n"
                f"STDOUT (first 2000 chars):\n{stdout_text[:2000]}\n"
            )

        # Prefer parsing the file output
        if out_path.exists():
            raw = out_path.read_text(encoding="utf-8", errors="replace")
            blocks = parse_scalefree_output(raw)
            return ScaleFreeResult(
                blocks=blocks,
                raw_text=raw,
                output_path=out_path,
                stdout=stdout_text,
                stderr=stderr_text,
            )

        # Optional fallback: parse structured blocks from stdout
        if parse_stdout_fallback:
            structured = "\n".join(
                [
                    ln
                    for ln in stdout_text.splitlines()
                    if ln.strip().startswith("#")
                    or re.match(r"^\s*[0-9\.\-\+EeDd]+\s+", ln)
                ]
            )
            raw = structured
            blocks = parse_scalefree_output(raw)
            return ScaleFreeResult(
                blocks=blocks,
                raw_text=raw,
                output_path=None,
                stdout=stdout_text,
                stderr=stderr_text,
            )

        raise RuntimeError(
            "Fortran returned success but output file was not found.\n"
            f"Expected: {out_path}\n"
            "If the Fortran code does not prompt"
            + " for an output filename (or writes elsewhere),\n"
            "update the prompt match key in answers['Output file']"
            + " or enable parse_stdout_fallback=True."
        )


# ---------------------------------------------------------------------
# Simple functional API (what most users will call)
# ---------------------------------------------------------------------


def vprofile(
    *,
    exe_path: Union[str, Path],
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
    output_path: Optional[Union[str, Path]] = "out.txt",
    timeout_s: int = 120,
    parse_stdout_fallback: bool = False,
    debug_prompts: bool = False,
) -> ScaleFreeResult:
    """
    Convenience function that instantiates a runner and executes vprofile.

    This keeps the public API simple:
        import scalefree
        res = scalefree.vprofile(exe_path=..., ...)
    """
    runner = ScaleFreeRunner(exe_path=exe_path)
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
