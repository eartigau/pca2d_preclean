"""Timestamped, colour-coded console logging.

Convention (shared across projects):
    YYMMDD HH:MM:SS.SS | message

Colours encode the semantic role of the line:
    green  -> general info / progress narration
    blue   -> a reported numeric value (count, median, RMS, ...)
    orange -> something skipped, or a non-critical issue that does not stop the run
    red    -> explains why the run is stopping
"""

from __future__ import annotations

import sys
from datetime import datetime

_COLOURS = {
    "info": "\033[32m",    # green
    "value": "\033[34m",   # blue
    "warn": "\033[33m",    # orange/yellow
    "error": "\033[31m",   # red
}
_RESET = "\033[0m"

_USE_COLOUR = sys.stdout.isatty()


def set_colour(enabled: bool) -> None:
    """Force colour output on or off (auto-detected from the tty by default)."""
    global _USE_COLOUR
    _USE_COLOUR = enabled


def log(message: str, level: str = "info") -> None:
    """Print one timestamped, colour-coded line.

    level is one of 'info', 'value', 'warn', 'error'.
    """
    now = datetime.now()
    stamp = now.strftime("%y%m%d %H:%M:%S.") + "%02d" % (now.microsecond // 10000)
    text = "%s | %s" % (stamp, message)
    if _USE_COLOUR:
        text = _COLOURS.get(level, _COLOURS["info"]) + text + _RESET
    print(text, flush=True)
