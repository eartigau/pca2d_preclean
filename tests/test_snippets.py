"""The raw flux around each window is kept once, and the figures read it back.

Drawing the figures used to reopen every spectrum and resample each onto the
whole grid, to keep a few thousand columns of it; and the pages were computed
on the whole grid too. These pin the replacements to what they replace: a
snippet reads back exactly, a figure with a snippet never opens a file, a
block resampled alone is the same block resampled inside a wider grid, and a
page computed on its block is the page computed on the whole grid.
"""

import os

import numpy as np
import pytest
from astropy.table import Table

from pca2d import cache
from pca2d.grids import window_block

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_a_snippet_reads_back_exactly(tmp_path):
    by_file = {"0001t.fits": np.arange(20, dtype=np.float32).reshape(2, 10),
               "0002t.fits": -np.ones((2, 10), dtype=np.float32)}
    cache.write_snippet(str(tmp_path), 100, 110, by_file)
    back = cache.read_snippet(str(tmp_path), 100, 110)
    assert set(back) == set(by_file)
    for name in by_file:
        assert np.array_equal(back[name], by_file[name])
    assert cache.read_snippet(str(tmp_path), 100, 111) is None, "wrong block served"


def test_with_a_snippet_no_spectrum_is_opened(tmp_path, monkeypatch):
    from pca2d import plotting, tfits

    def refuse(*a, **k):
        raise AssertionError("a spectrum was opened although the snippet was there")

    monkeypatch.setattr(tfits, "read_tfits", refuse)
    even = np.full((2, 5), 1.0, dtype=np.float32)
    even[1] = 2.0
    cache.write_snippet(str(tmp_path), 0, 5, {"a.fits": even})
    out = plotting.raw_log_flux_window(str(tmp_path), None, ["a.fits", "a.fits"],
                                       [0, 1], np.arange(5.0), 0, 5)
    assert np.array_equal(out[0], even[0]) and np.array_equal(out[1], even[1])


SPECTRA = os.path.join(REPO, "data", "TOI2120")


def _a_readable_spectrum():
    """Whether this checkout HAS one, which a folder of links does not settle:
    data/ here is links onto a shared disk, and a disk that is not mounted
    leaves the names behind and the files unreachable."""
    import glob

    for path in sorted(glob.glob(os.path.join(SPECTRA, "*t.fits")))[:1]:
        return os.path.exists(path)
    return False


@pytest.mark.skipif(not _a_readable_spectrum(),
                    reason="needs a real t.fits that can be read")
def test_a_block_resampled_alone_is_the_block_resampled_in_a_wider_grid():
    """What makes a snippet from the build equal to the figure's own read."""
    import glob

    from pca2d.config import load_config
    from pca2d.grids import magic_grid
    from pca2d.tfits import raw_block, read_tfits

    path = sorted(glob.glob(os.path.join(SPECTRA, "*t.fits")))[0]
    config = load_config(os.path.join(REPO, "config.yaml"), instrument="SPIROU")
    dom = config["domain"]
    grid = magic_grid(dom["wave0"], dom["dv"], dom["wave_min"], dom["wave_max"])
    a0, b0 = window_block(grid, 1669.5, 5.0, dom["dv"])
    payload = read_tfits(path)
    alone = raw_block(payload, grid, a0, b0, config)
    wide = raw_block(payload, grid, a0 - 400, b0 + 400, config)[:, 400:-400]
    inner = slice(5, -5)            # drop_isolated may differ at the very ends
    assert np.array_equal(np.isnan(alone[:, inner]), np.isnan(wide[:, inner]))
    np.testing.assert_allclose(alone[:, inner], wide[:, inner], rtol=1e-6,
                               equal_nan=True)


# ------------------------------------------ a page on its block, exactly ----
N_ROWS, N_PIX, DV = 8, 1400, 0.5


@pytest.fixture
def cube_and_fit(tmp_path):
    r = np.random.default_rng(12)
    path = tmp_path / "cube"
    path.mkdir()
    grid = 1500.0 * np.exp(np.arange(N_PIX) * DV / 299792.458)
    np.save(path / "grid.npy", grid)
    np.save(path / "data.npy", (0.05 * r.normal(size=(N_ROWS, N_PIX))).astype(np.float32))
    np.save(path / "sigma.npy", np.full((N_ROWS, N_PIX), 0.02, dtype=np.float32))
    meta = Table()
    meta["filename"] = ["%03dt.fits" % (i // 2) for i in range(N_ROWS)]
    meta["snr_band"] = np.full(N_ROWS, 100.0)
    meta["berv"] = np.linspace(-20.0, 20.0, N_ROWS)
    meta.write(path / "meta.fits")
    smooth = lambda k: np.array([np.convolve(r.normal(size=N_PIX),
                                             np.exp(-0.5 * (np.arange(-15, 16) / 3.0) ** 2),
                                             mode="same") for _ in range(k)])
    fit = {"P": smooth(2), "Q": smooth(3), "a": r.normal(size=(N_ROWS, 2)),
           "b": r.normal(size=(N_ROWS, 3)), "berv": np.asarray(meta["berv"])}
    return str(path), fit, grid


def test_a_page_on_its_block_is_the_page_on_the_whole_grid(cube_and_fit):
    from pca2d.figures.sequence import window_arrays
    from pca2d.plotting import live_mask
    from pca2d.twoframe import LanczosShifter, load_cube, star_model

    cube, fit, grid = cube_and_fit
    delta = -fit["berv"] / DV
    means, group = np.zeros((1, N_PIX)), np.zeros(N_ROWS, dtype=int)
    centre, width = float(grid[700]), 0.4

    got = window_arrays(cube, fit, means, group, grid, DV, delta, centre, width)

    # the full-grid computation it replaces
    _, data, w0, _ = load_cube(cube)
    shifter = LanczosShifter(N_PIX, a=8, max_shift=int(np.ceil(np.abs(delta).max())) + 2)
    star = star_model(fit["P"], fit["a"], shifter, delta, N_ROWS, N_PIX)
    earth = fit["b"] @ fit["Q"]
    alive = live_mask(shifter.rows(w0, -delta))
    reference = {"given": data, "model": star + earth, "corrected": data - earth,
                 "nostar": data - star, "resid": data - star - earth}
    inside = np.where((grid >= centre - width / 2) & (grid <= centre + width / 2))[0]
    for name, arr in reference.items():
        want = shifter.rows(arr, -delta)
        want[~alive] = np.nan
        np.testing.assert_allclose(got["home"][name][:, inside - got["a0"]],
                                   want[:, inside].astype(np.float32),
                                   rtol=1e-5, atol=1e-6, equal_nan=True,
                                   err_msg="panel %s differs inside the window" % name)
