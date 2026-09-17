"""A loop with no bar to draw still has to say where it has got to.

A run whose output is a file, or a window's pipe, printed "cube footprint: 310
rows x 577002 pixels = 2.86 GB" and then nothing at all while it read three
hundred spectra off a shared disk: the bar that would have shown it moving was
disabled for not having a terminal, and nothing took its place.
"""
import os
import sys

from pca2d import progress


def narrated(monkeypatch, capsys, items, every=0.0, **kwargs):
    monkeypatch.setattr(progress.sys.stderr, "isatty", lambda: False, raising=False)
    monkeypatch.delenv("PCA2D_COLOUR", raising=False)
    loop = progress.bar(items, desc="reading spectra", unit="file", **kwargs)
    assert isinstance(loop, progress._Narrator), "no terminal, so no bar to draw"
    loop.every = every
    seen = list(loop)
    return seen, capsys.readouterr().out


def test_a_loop_with_no_bar_says_how_far_it_has_got(monkeypatch, capsys):
    seen, out = narrated(monkeypatch, capsys, list(range(4)))
    assert seen == [0, 1, 2, 3], "it is still the loop it wraps"
    lines = [line for line in out.splitlines() if "reading spectra" in line]
    assert len(lines) == 4, "one line per step, at every=0"
    assert "4/4 files" in lines[-1]
    assert "file/s" in lines[-1], "how fast, so a slow disk is visible as one"


def test_it_says_nothing_faster_than_asked(monkeypatch, capsys):
    seen, out = narrated(monkeypatch, capsys, list(range(50)), every=3600.0)
    assert seen == list(range(50))
    assert "reading spectra" not in out, \
        "half a minute apart, or a log nobody reads"


def test_a_length_it_cannot_know_is_counted_not_estimated(monkeypatch, capsys):
    _, out = narrated(monkeypatch, capsys, iter(range(3)))
    last = [line for line in out.splitlines() if "reading spectra" in line][-1]
    assert "3 files" in last and "/" not in last.split("files")[0].split("| ")[-1]
    assert "left" not in last, "no total, so no honest estimate of what is left"


def test_the_window_that_redraws_a_line_gets_the_bar_itself(monkeypatch):
    """The window reads the run through a pipe and redraws a carriage-returned
    line in place (gui.cut_output), so a bar there is a bar."""
    monkeypatch.setattr(progress.sys.stderr, "isatty", lambda: False, raising=False)
    monkeypatch.setenv("PCA2D_COLOUR", "1")
    assert progress._drawn()
    loop = progress.bar([1, 2], desc="reading spectra", unit="file")
    assert not isinstance(loop, progress._Narrator)
