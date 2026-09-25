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
    # both order parities of an exposure come from one file, as in a real
    # cube: a figure's count of exposures is a count of these names
    meta["filename"] = ["%04dt.fits" % (i // 2) for i in range(N_ROWS)]
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
    run(cube, out, "--velocity-term")
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


def test_a_figure_counts_exposures_not_rows():
    """The correlation matrix said 632 exposures for 316: two rows each."""
    from pca2d.twoframe import exposures_label
    names = ["a.fits", "a.fits", "b.fits", "b.fits", "c.fits", "c.fits"]
    keep = np.array([True, True, True, False, False, False])
    assert exposures_label(names, keep) == "2 exposures"
    assert exposures_label(names, np.ones(6, dtype=bool)) == "3 exposures"
    assert exposures_label(None, keep) == "3 rows", "no names, so it says rows"


def test_the_fit_s_correlation_matrix_is_titled_in_exposures(cube, tmp_path):
    """The title the fit writes on correlations.pdf, on twelve rows of six exposures."""
    from unittest import mock
    from pca2d import plotting
    with mock.patch.object(plotting, "plot_correlations") as drawn:
        run(cube, tmp_path / "fit")
    titles = [call.kwargs.get("title") for call in drawn.call_args_list]
    assert titles and all(t == "6 exposures" for t in titles), titles


def test_the_star_template_reads_the_cube_in_blocks_and_gets_the_same_numbers(cube, tmp_path):
    """A block of columns at a time, as the lbl stage reads it beside a running
    fit, and every column comes out as it does from the whole cube."""
    from pca2d.lbltemplate import star_coverage
    run(cube, tmp_path / "fit")
    fit = np.load(tmp_path / "fit" / "fit.npz")
    whole = star_coverage(cube, fit, block=None)
    blocks = star_coverage(cube, fit, block=300)
    assert np.array_equal(whole[0], blocks[0]) and np.array_equal(whole[3], blocks[3])
    for name, x, y in zip(("count", "weight", "", "residual"), whole[1:], blocks[1:]):
        if name:
            np.testing.assert_allclose(x, y, rtol=1e-6, atol=1e-10, equal_nan=True,
                                       err_msg=name)
    assert np.isfinite(whole[4]).any(), "the residual means exist somewhere"


def test_the_residual_clip_finds_a_spike_and_reads_the_cube_in_blocks(cube, tmp_path):
    """A spike the fit never saw, fifty times the cube's noise, is cut and its
    neighbours are not; and a block of columns at a time cuts exactly what the
    whole cube does. The toy cube's lines are an ln f of -1 deep against a
    noise of 0.01, so three sweeps leave model errors on their flanks that the
    clip cuts too; the test looks at the spike and the samples around it."""
    from pca2d.outliers import residual_outliers
    run(cube, tmp_path / "fit")
    fit = np.load(tmp_path / "fit" / "fit.npz")
    data = np.load(os.path.join(cube, "data.npy"))
    data[4, 1000] += 0.5
    np.save(os.path.join(cube, "data.npy"), data)
    whole, _ = residual_outliers(cube, fit, 8.0, 151, block=None)
    blocks, report = residual_outliers(cube, fit, 8.0, 151, block=300)
    assert whole.keys() == blocks.keys()
    assert report["n_columns"] == 0, "no column test was asked for"
    for name in whole:
        assert np.array_equal(whole[name], blocks[name]), name
    row = whole["0002t.fits"][0]                  # row 4: exposure 2, even orders
    assert row[1000]
    assert not row[985:1000].any() and not row[1001:1016].any(), (
        "the noise beside the spike is left alone")


def test_a_column_every_exposure_disagrees_on_goes_from_all_of_them(cube, tmp_path):
    """The same cube with one observer column given four times its noise in
    every exposure: the column leaves every file, and a spike in one exposure
    still leaves only that one."""
    from pca2d.outliers import residual_outliers
    run(cube, tmp_path / "fit")
    fit = np.load(tmp_path / "fit" / "fit.npz")
    data = np.load(os.path.join(cube, "data.npy"))
    rng = np.random.default_rng(11)
    # 1800 is clear of every line of the toy cube: a deep one is clipped out
    # of the fit's weights, so nothing is measured there to disagree about
    data[:, 1800] += rng.normal(scale=0.04, size=data.shape[0])   # noise is 0.01
    data[4, 1000] += 0.5
    np.save(os.path.join(cube, "data.npy"), data)
    masks, report = residual_outliers(cube, fit, 3.0, 151, block=None,
                                      column_frac=0.10, column_chi2=1.5,
                                      column_min_rows=2)
    assert report["columns"][:, 1800].all(), (
        "frac %s, chi2 %s" % (report["frac"][:, 1800], report["chi2"][:, 1800]))
    assert all(m[0][1800] and m[1][1800] for m in masks.values()), \
        "every exposure loses it, whatever its own residual there was"
    assert masks["0002t.fits"][0][1000], "and the one-exposure spike is still cut"
    # Six exposures per parity is too few for the rule to tell a bad column
    # from one bad exposure: a single clipped exposure is already 17% of them.
    # That discrimination is tested on 120 rows in test_column_outliers.py;
    # what this cube shows is the guard that refuses to judge on too few.
    _, careful = residual_outliers(cube, fit, 3.0, 151, block=None,
                                  column_frac=0.10, column_chi2=1.5,
                                  column_min_rows=8)
    assert careful["n_columns"] == 0, "six exposures of a parity give no verdict"


def test_the_last_iterate_can_be_kept_rather_than_the_best(cube, tmp_path):
    """--keep last and --patience, an experiment's options; the defaults keep
    the lowest chi2 and stop after two worse sweeps, as before."""
    out = tmp_path / "last"
    run(cube, out, "--keep", "last", "--patience", "99")
    assert int(np.load(out / "fit.npz")["best_iter"]) == 2, "the last of three sweeps"


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
