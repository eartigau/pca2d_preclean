"""The whole fit, on a cube small enough to be a test.

Twice in one afternoon a real run died AFTER the fit: once because four figure
scripts called log() without importing it, once because an extra error column
became an eighth error bar and the plotting could not broadcast it. Both cost
the sixteen minutes of fit that came before them, and neither could have been
caught by a unit test of the solver, because both were in what happens once the
solver is done.

So this runs main() from end to end on twelve rows and two thousand samples: it
builds a cube on disk, fits it, and asserts that every product exists and that
the coefficient table says what it should. It takes about a second, which is
the price of never learning this way again.
"""

import os

import numpy as np
import pytest
from astropy.table import Table

from pca2d.twoframe import main

N_ROWS, N_PIX = 12, 2000
DV = 0.5


@pytest.fixture
def cube(tmp_path):
    """A small two-frame cube on disk, in the shape read_cube_files expects."""
    r = np.random.default_rng(4)
    x = np.arange(N_PIX, dtype=float)
    grid = 1500.0 * np.exp(np.arange(N_PIX) * DV / 299792.458)

    star = -sum(np.exp(-0.5 * ((x - c) / 5.0) ** 2) for c in (300, 900, 1500))
    earth = -sum(np.exp(-0.5 * ((x - c) / 4.0) ** 2) for c in (600, 1200))
    berv = np.linspace(-12.0, 12.0, N_ROWS // 2).repeat(2)
    data = np.zeros((N_ROWS, N_PIX))
    for i in range(N_ROWS):
        shift = int(round(-berv[i] / DV))
        data[i] = np.roll(star, shift) * (1 + 0.05 * r.normal())
        data[i] += earth * (1.0 + 0.3 * r.normal())
    data += 0.01 * r.normal(size=data.shape)
    data -= data.mean(axis=1, keepdims=True)

    path = tmp_path / "cube"
    path.mkdir()
    np.save(path / "grid.npy", grid)
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", np.full_like(data, 0.01))
    meta = Table()
    meta["filename"] = ["%04dt.fits" % i for i in range(N_ROWS)]
    meta["bjd"] = 2459000.0 + np.arange(N_ROWS) * 1.7
    meta["berv"] = berv
    meta["airmass"] = np.full(N_ROWS, 1.2)
    meta["snr_band"] = np.full(N_ROWS, 100.0)
    meta["exposure"] = np.arange(N_ROWS) // 2
    meta["parity"] = np.arange(N_ROWS) % 2
    meta.write(path / "meta.fits", overwrite=True)
    return str(path)


def run(cube, outdir, *extra):
    main(["--cube", cube, "--outdir", str(outdir), "--iters", "3",
          "-k", "2", "-j", "2", "--max-mad", "0", *extra])


def test_the_fit_writes_every_product(cube, tmp_path):
    out = tmp_path / "fit"
    run(cube, out)
    for name in ("coefficients.csv", "fit.npz", "twoframe_components.fits",
                 "coefficients_vs_time.pdf", "variance.pdf", "variance.csv",
                 "components.pdf", "correlations.pdf", "correlations.csv"):
        assert os.path.exists(out / name), "%s was not written" % name


def test_the_table_carries_one_velocity_per_row(cube, tmp_path):
    out = tmp_path / "fit"
    run(cube, out)
    table = Table.read(out / "coefficients.csv", format="csv")
    for column in ("alpha", "vrad_fit", "evrad_fit"):
        assert column in table.colnames, "%s is missing" % column
    assert len(table) == N_ROWS
    # tied by exposure: the two parities of one exposure share their velocity
    v = np.asarray(table["vrad_fit"], float)
    assert np.allclose(v[0::2], v[1::2]), (
        "the two parities of an exposure were given different velocities")


def test_the_same_fit_without_the_term_writes_no_velocity(cube, tmp_path):
    out = tmp_path / "plain"
    run(cube, out, "--no-velocity-term")
    table = Table.read(out / "coefficients.csv", format="csv")
    assert "vrad_fit" not in table.colnames
    assert os.path.exists(out / "twoframe_components.fits")


def test_a_count_of_rows_is_reported_as_a_count_of_exposures():
    """Two rows per exposure, one per order parity: 632 rows are 316 spectra."""
    from pca2d.twoframe import count_exposures

    meta = Table()
    meta["filename"] = ["%04dt.fits" % (i // 2) for i in range(8)]   # 4 exposures
    assert count_exposures(meta) == 4
    dropped = np.array([True, True, False, False, False, False, False, True])
    assert count_exposures(meta, dropped) == 2, (
        "rows 0-1 are one exposure and row 7 is half of another")
    assert count_exposures(meta, ~dropped) == 3
