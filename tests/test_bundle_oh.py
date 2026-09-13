"""The bundle draws the OH figure only on spectra that carry an OH model."""
import numpy as np
from astropy.io import fits

from pca2d.figures.bundle import has_oh_model


def _write(path, names):
    """A t.fits with these extensions, EXTNAME in the case the DRS writes it:
    astropy upper-cases a name given to the HDU, not one set on its header."""
    hdus = [fits.PrimaryHDU()]
    for name in names:
        hdu = fits.ImageHDU(np.zeros((2, 4)))
        hdu.header["EXTNAME"] = name
        hdus.append(hdu)
    fits.HDUList(hdus).writeto(path)


def test_spirou_files_have_an_oh_model(tmp_path):
    _write(tmp_path / "2400000o_pp_e2dsff_tcorr_ABt.fits",
           ["FluxAB", "WaveAB", "OHLine"])
    assert has_oh_model(str(tmp_path))


def test_nirps_files_do_not(tmp_path):
    _write(tmp_path / "NIRPS.2023-01-20T08:42:08.941t.fits", ["FluxA", "WaveA"])
    assert not has_oh_model(str(tmp_path))


def test_no_spectra_is_nothing_to_draw_rather_than_a_figure_that_fails(tmp_path):
    """A joint run's object is a name, not a folder, and this figure needs real
    spectra: every joint run had it on the bundle's failure page."""
    assert not has_oh_model(str(tmp_path))
    assert not has_oh_model(str(tmp_path / "PROXIMA+GJ1+GJ3090"))


def test_a_figure_that_skipped_a_window_is_reported_even_though_it_succeeded():
    """A page that is not there is not an error anybody notices: three windows
    of the joint run were skipped for "too few rows" and nobody saw it."""
    from pca2d.figures.bundle import SKIPPED
    said = ["260913 09:34:48.13 |   1667.0-1672.0 nm: too few rows, skipped",
            "   1220.0-1222.0 nm: outside the grid, skipped",
            "   oh_residual.py skipped: these spectra carry no OHLine extension",
            "output.windows: 2450:5 is centred at 2450.0 nm, ... so it is not drawn",
            "   no spectra in data/X, so nothing to draw the airglow from"]
    for line in said:
        assert SKIPPED.search(line), line
    for line in ("260913 10:05:50.30 |   1199.3-1201.3 nm: given 0.0698",
                 "wrote outputs/x.pdf: 6 pages from 8 figures",
                 "drawing 8 windows and binding everything into one PDF"):
        assert not SKIPPED.search(line), line


def test_the_runner_collects_them(monkeypatch):
    import subprocess

    from pca2d.figures import bundle

    class Done:
        returncode = 0
        stdout = ("260913 09:34:48.13 |   1667.0-1672.0 nm: too few rows, skipped\n"
                  "260913 09:34:49.00 |   1199.3-1201.3 nm: given 0.07\n")
        stderr = ""

    monkeypatch.setattr(subprocess, "run", lambda *a, **k: Done())
    failures = []
    assert bundle.run(["python", "sequence.py"], failures) is True
    assert len(failures) == 1, failures
    assert "too few rows" in failures[0][1]
    assert failures[0][2] is True, "a skip is marked as one, not as a failure"


def test_a_real_failure_is_not_marked_as_a_skip(monkeypatch):
    import subprocess

    from pca2d.figures import bundle

    class Broke:
        returncode = 1
        stdout = ""
        stderr = "Traceback ... IndexError"

    monkeypatch.setattr(subprocess, "run", lambda *a, **k: Broke())
    failures = []
    assert bundle.run(["python", "x.py"], failures) is False
    assert failures[0][2] is False and "IndexError" in failures[0][1]
