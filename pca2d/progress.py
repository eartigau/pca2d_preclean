"""Telling the user what is happening, for the hour it takes.

A run reads three hundred spectra, builds a five gigabyte cube, iterates a
decomposition sixteen times and writes three hundred corrected files. Left to
itself it prints nothing for forty minutes, and there is no way to tell a
working run from a wedged one. Every loop that can take more than a few seconds
goes through here.

Two things, deliberately kept apart:

`bar` is a tqdm with `leave=False`, so a finished loop erases its own bar and
the scrollback holds the narration rather than a graveyard of completed
progress bars. It writes to stderr and disables itself when stderr is not a
terminal, so a log file stays readable.

`log` is that narration: one line per thing that happened, timestamped
`YYMMDD HH:MM:SS.SS | message` and coloured by role. Green for progress, blue
for a number worth remembering, orange for something skipped, red for why the
run is stopping. `stage` bookends a phase and reports how long it took, because
"still going" and "this one is slow" look identical without it.
"""

from __future__ import annotations

import contextlib
import sys
import time

from .logger import log

try:
    from tqdm.auto import tqdm as _tqdm
except ImportError:                                           # pragma: no cover
    _tqdm = None


def bar(iterable=None, total=None, desc="", unit="it", **kwargs):
    """A progress bar that leaves nothing behind when it finishes.

    Falls back to the plain iterable if tqdm is not installed, so nothing here
    is a hard dependency of the science.
    """
    if _tqdm is None:                                         # pragma: no cover
        return iterable if iterable is not None else _Null()
    return _tqdm(iterable, total=total, desc=desc, unit=unit, leave=False,
                 file=sys.stderr, dynamic_ncols=True,
                 disable=not sys.stderr.isatty(), **kwargs)


class _Null:                                                  # pragma: no cover
    """What `bar` returns without tqdm, so callers need no special case."""

    def update(self, *a, **k):
        pass

    def close(self):
        pass

    def set_postfix_str(self, *a, **k):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


@contextlib.contextmanager
def stage(title, level="info"):
    """Announce a phase, then say how long it took.

    The closing line is the useful one: it turns "the run has been going for
    forty minutes" into "the fit took seventeen of them", which is the
    difference between waiting and debugging.
    """
    log(title, level)
    started = time.time()
    try:
        yield
    finally:
        seconds = time.time() - started
        log("%s: %s" % (title, human(seconds)), "value")


def human(seconds):
    """A duration a person reads at a glance."""
    seconds = float(seconds)
    if seconds < 90:
        return "%.1f s" % seconds
    if seconds < 5400:
        return "%d min %02d s" % (int(seconds // 60), int(seconds % 60))
    return "%d h %02d min" % (int(seconds // 3600), int((seconds % 3600) // 60))


def size(nbytes):
    """Bytes, in the unit a person would have used."""
    nbytes = float(nbytes)
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if abs(nbytes) < 1024 or unit == "TB":
            return "%.1f %s" % (nbytes, unit)
        nbytes /= 1024
    return "%.1f TB" % nbytes                                 # pragma: no cover
