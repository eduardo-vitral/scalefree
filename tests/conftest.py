from __future__ import annotations

from pathlib import Path
import pytest

from scalefree import ScaleFreeRunner


def repo_root() -> Path:
    # tests/ is at repo_root/tests/
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def runner() -> ScaleFreeRunner:
    """
    Use the package's backend resolution (env var / cached exe / auto-compile).
    This matches the intended user experience after `pip install scalefree`.
    """
    try:
        return ScaleFreeRunner()
    except Exception as e:
        pytest.skip(f"ScaleFree backend not available: {e}")


@pytest.fixture(scope="session")
def ref_dir() -> Path:
    d = repo_root() / "tests" / "data"
    if not d.exists():
        pytest.skip("Missing tests/data directory with reference outputs.")
    return d
