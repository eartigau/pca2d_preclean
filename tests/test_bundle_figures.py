"""Clearing the report's figures before writing them again.

On an exFAT disk, macOS writes a `._name` beside every file and removes it
with the file itself: a joint run died on the sidecar of a figure it had just
removed (2026-09-17)."""

import os

from pca2d.figures.bundle import clear_figures


def test_the_sidecars_and_the_missing_are_not_a_crash(tmp_path):
    folder = tmp_path / "figures"
    folder.mkdir()
    for name in ("01-one.pdf", "._01-one.pdf", "02-two.pdf", "notes.txt",
                 ".DS_Store"):
        (folder / name).write_bytes(b"x")
    assert clear_figures(str(folder)) == 2, "the two figures, nothing else"
    left = sorted(os.listdir(folder))
    assert left == [".DS_Store", "._01-one.pdf", "notes.txt"]
    assert clear_figures(str(folder)) == 0
    assert clear_figures(str(tmp_path / "nowhere")) == 0


def test_a_figure_that_vanishes_under_it_is_not_a_crash(tmp_path, monkeypatch):
    folder = tmp_path / "figures"
    folder.mkdir()
    (folder / "01-one.pdf").write_bytes(b"x")
    (folder / "02-two.pdf").write_bytes(b"x")
    real = os.remove

    def flaky(path):
        if path.endswith("01-one.pdf"):
            real(path)
            raise FileNotFoundError(path)     # as a sidecar's removal does
        return real(path)

    monkeypatch.setattr(os, "remove", flaky)
    assert clear_figures(str(folder)) == 1
    assert os.listdir(folder) == []
