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
