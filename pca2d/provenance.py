"""Which code made a run: the git commit, and whether it had local changes.

Written into every run's resolved configuration, its fit.npz, the headers of
its corrected spectra and its LBL template, so that a result can be traced to
the code that produced it and made again.
"""
import functools
import os
import subprocess

from . import __version__

HERE = os.path.dirname(os.path.abspath(__file__))


def _git(*args):
    try:
        r = subprocess.run(["git", "-C", HERE, *args], capture_output=True,
                           text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout.strip() if r.returncode == 0 else None


@functools.lru_cache(maxsize=1)
def code_version():
    """{"commit": the full hash, None outside a git checkout; "dirty": True
    when a tracked file differs from that commit; "version": the package's}.

    Asked once per process: the code that started a run is the code that made
    it, and a correct stage writing hundreds of files asks for every one."""
    commit = _git("rev-parse", "HEAD")
    dirty = None
    if commit:
        # tracked files only: a scratch file beside the code changes nothing
        status = _git("status", "--porcelain", "--untracked-files=no")
        dirty = bool(status) if status is not None else None
    return {"commit": commit, "dirty": dirty, "version": __version__}


def stamp(version=None):
    """The short form, for a FITS card: 12 hex digits, '+' when a tracked file
    had changed, 'unknown' outside a git checkout."""
    v = version or code_version()
    if not v["commit"]:
        return "unknown"
    return v["commit"][:12] + ("+" if v["dirty"] else "")
