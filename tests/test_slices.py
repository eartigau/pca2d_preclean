"""A figure of one window computed on that window alone is the same figure.

The figures used to load the whole cube and carry every model across the whole
grid to draw a few thousand columns: 4.4 GB read for 79 MB used. They now read
a block of columns around each window. That is only acceptable if nothing
changes inside the window, and these tests hold it to that exactly.
"""

import numpy as np
import pytest
from astropy.table import Table

from pca2d.grids import MAX_BERV_KMS, parse_window, window_block
from pca2d.twoframe import LanczosShifter, load_cube, star_model

N_ROWS, N_PIX, DV = 10, 1500, 0.5


@pytest.fixture
def cube(tmp_path):
    r = np.random.default_rng(9)
    path = tmp_path / "cube"
    path.mkdir()
    np.save(path / "grid.npy", 1500.0 * np.exp(np.arange(N_PIX) * DV / 299792.458))
    data = 0.05 * r.normal(size=(N_ROWS, N_PIX))
    sigma = np.full((N_ROWS, N_PIX), 0.02)
    sigma[:, 700:720] = np.nan                       # a hole
    data[:, 900:905] = -2.0                          # below ln_clip_low
    trans = np.clip(1.0 - 0.8 * np.exp(-0.5 * ((np.arange(N_PIX) - 1100) / 30.0) ** 2),
                    0, 1) * np.ones((N_ROWS, 1))
    np.save(path / "data.npy", data.astype(np.float32))
    np.save(path / "sigma.npy", sigma.astype(np.float32))
    np.save(path / "trans.npy", trans.astype(np.float32))
    meta = Table()
    meta["filename"] = ["%03dt.fits" % (i // 2) for i in range(N_ROWS)]
    meta["snr_band"] = np.full(N_ROWS, 100.0)
    meta["berv"] = np.linspace(-20, 20, N_ROWS)
    meta.write(path / "meta.fits")
    return str(path)


@pytest.mark.parametrize("cols", [np.arange(600, 1200), np.arange(10, 300),
                                  np.arange(N_PIX - 250, N_PIX)])
def test_a_slice_of_the_cube_is_the_cube_sliced(cube, cols):
    """Including at the grid's own ends, whose EDGE columns stay zeroed."""
    grid, data, w, _ = load_cube(cube)
    g_s, d_s, w_s, _ = load_cube(cube, columns=cols)
    assert np.array_equal(g_s, grid[cols])
    assert np.array_equal(w_s, w[:, cols]), "weights differ on the slice"
    assert np.array_equal(d_s, data[:, cols])


def test_the_block_covers_the_largest_possible_shift():
    grid = 1500.0 * np.exp(np.arange(N_PIX) * DV / 299792.458)
    centre = float(grid[750])                         # the middle of 1500-1503.75 nm
    a0, b0 = window_block(grid, *parse_window("%.4f:0.5" % centre), DV)
    inside = np.where((grid >= centre - 0.25) & (grid <= centre + 0.25))[0]
    assert inside[0] - a0 >= MAX_BERV_KMS / DV + 8
    assert b0 - 1 - inside[-1] >= MAX_BERV_KMS / DV + 8
    assert window_block(grid, 9999.0, 2.0, DV) is None


def test_the_star_model_on_the_block_is_exact_inside_the_window():
    """The carry is local, so a block with the margin needs nothing outside it."""
    r = np.random.default_rng(2)
    x = np.arange(N_PIX)
    P = np.array([np.convolve(r.normal(size=N_PIX), np.exp(-0.5 * (np.arange(-20, 21) / 4.0) ** 2),
                              mode="same")])
    a = r.normal(size=(N_ROWS, 1))
    delta = r.uniform(-MAX_BERV_KMS / DV, MAX_BERV_KMS / DV, size=N_ROWS)
    full = star_model(P, a, LanczosShifter(N_PIX, a=8), delta, N_ROWS, N_PIX)

    grid = 1500.0 * np.exp(x * DV / 299792.458)
    a0, b0 = window_block(grid, float(grid[750]), 1.0, DV)
    block = star_model(P[:, a0:b0], a, LanczosShifter(b0 - a0, a=8), delta,
                       N_ROWS, b0 - a0)
    inside = np.where((grid >= grid[750] - 0.5) & (grid <= grid[750] + 0.5))[0]
    np.testing.assert_allclose(block[:, inside - a0], full[:, inside],
                               rtol=1e-10, atol=1e-12)
