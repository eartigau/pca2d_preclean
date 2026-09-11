"""--mean star: one star spectrum per order parity, nothing in the observer frame.

The two parities see the star's lines differently where the orders overlap
(Proxima, 1201.1 nm, 2026-09-11), and an observer-frame mean cannot hold a
difference that moves with the star. This mode takes a star-frame median per
parity out before the fit and no observer mean at all; the correction then
divides out the observer block alone, which is still exactly panel 3.
"""
import numpy as np
import pytest
from astropy.table import Table

from pca2d.reconstruct import correction_on_grid, load_model, order_correction
from pca2d.twoframe import fit_means, load_cube, main, mean_rows, row_parity

N_ROWS, N_PIX, DV = 48, 2000, 0.5
DEPTH = {0: -0.30, 1: -0.24}          # the even orders see the lines deeper
LINES = (300, 900, 1500)


@pytest.fixture(scope="module")
def fitted(tmp_path_factory):
    r = np.random.default_rng(7)
    x = np.arange(N_PIX, dtype=float)
    grid = 1500.0 * np.exp(np.arange(N_PIX) * DV / 299792.458)
    shape = sum(np.exp(-0.5 * ((x - c) / 5.0) ** 2) for c in LINES)
    earth = -0.2 * sum(np.exp(-0.5 * ((x - c) / 4.0) ** 2) for c in (600, 1200))
    berv = np.linspace(-12.0, 12.0, N_ROWS // 2).repeat(2)
    data = np.zeros((N_ROWS, N_PIX))
    for i in range(N_ROWS):
        data[i] = np.roll(DEPTH[i % 2] * shape, int(round(-berv[i] / DV)))
        data[i] += earth * (1.0 + 0.3 * r.normal())
    data += 0.005 * r.normal(size=data.shape)
    path = tmp_path_factory.mktemp("cube")
    np.save(path / "grid.npy", grid)
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", np.full_like(data, 0.005))
    meta = Table()
    meta["filename"] = ["%04dt.fits" % (i // 2) for i in range(N_ROWS)]
    meta["bjd"] = 2459000.0 + np.arange(N_ROWS) * 1.7
    meta["berv"] = berv
    meta["airmass"] = np.full(N_ROWS, 1.2)
    meta["snr_band"] = np.full(N_ROWS, 100.0)
    meta["exposure"] = np.arange(N_ROWS) // 2
    meta["parity"] = np.arange(N_ROWS) % 2
    meta.write(path / "meta.fits", overwrite=True)
    out = tmp_path_factory.mktemp("fit")
    main(["--cube", str(path), "--outdir", str(out), "--iters", "3",
          "-k", "1", "-j", "2", "--max-mad", "0", "--mean", "star"])
    return str(path), out


def test_the_fit_keeps_a_star_spectrum_per_parity_and_no_observer_mean(fitted):
    _, out = fitted
    fit = np.load(out / "fit.npz")
    assert "templates" in fit.files
    assert not np.any(fit["means"]), "an observer-frame mean was taken out"
    templates = np.asarray(fit["templates"])
    assert templates.shape[0] == 2 and np.any(templates[0]) and np.any(templates[1])


def test_each_parity_s_spectrum_has_that_parity_s_depth(fitted):
    """At BERV 0 the lines sit at LINES, so that is where the star frame has them."""
    _, out = fitted
    templates = np.asarray(np.load(out / "fit.npz")["templates"])
    for parity in (0, 1):
        got = templates[parity][list(LINES)]
        assert np.allclose(got, DEPTH[parity], atol=0.02), (parity, got)


def test_the_correction_is_the_observer_block_alone_and_panel_3(fitted):
    cube, out = fitted
    fit = np.load(out / "fit.npz")
    model = load_model(str(out / "twoframe_components.fits"))
    assert model["templates"] is not None
    _, _, _, meta = load_cube(cube, columns=np.arange(1))
    n, grid = len(meta), model["grid"]
    parity = row_parity(meta, n)
    means, group = fit_means(fit, meta, n, grid.size)
    # sequence.window_arrays, panel 3: data - offset - earth, the offset zero
    subtracted = np.asarray(mean_rows(means, group, n)) + fit["b"] @ fit["Q"]
    coeffs = {str(row["filename"]): row for row in model["coeffs"]}
    inner = slice(100, -100)
    for i in range(n):
        row = coeffs[str(meta["filename"][i])]
        correction, _, j = correction_on_grid(model, row, 0, model["n_earth"])
        values, live = order_correction(model, correction, j, int(parity[i]),
                                        grid[inner])
        assert live.mean() > 0.9, "no support to correct without an observer mean"
        assert np.allclose(values[live], subtracted[i][inner][live],
                           rtol=0, atol=1e-9), i


@pytest.fixture(scope="module")
def smoothed(fitted, tmp_path_factory):
    """The same cube, with the star side smoothed to one resolution element of
    5 samples (R = 120 000 at 0.5 km/s)."""
    cube, _ = fitted
    out = tmp_path_factory.mktemp("fit_smoothed")
    main(["--cube", cube, "--outdir", str(out), "--iters", "3", "-k", "1",
          "-j", "2", "--max-mad", "0", "--mean", "star",
          "--star-resolution", "120000"])
    return out


def _roughness(v):
    """rms of the second difference where the vector has support: its structure
    at the scale of one sample."""
    v = np.asarray(v, dtype=float)
    live = v != 0
    d2 = np.diff(v, 2)[live[1:-1] & live[:-2] & live[2:]]
    return float(np.sqrt(np.mean(d2 ** 2)))


def test_the_star_side_is_smoothed_and_the_observer_side_is_not(fitted, smoothed):
    _, raw = fitted
    before, after = np.load(raw / "fit.npz"), np.load(smoothed / "fit.npz")
    for key in ("templates", "P"):
        for row_b, row_a in zip(np.atleast_2d(before[key]), np.atleast_2d(after[key])):
            assert _roughness(row_a) < 0.5 * _roughness(row_b), key
    rough_before = np.mean([_roughness(r) for r in before["Q"]])
    rough_after = np.mean([_roughness(r) for r in after["Q"]])
    assert rough_after > 0.5 * rough_before, "the observer components were smoothed"
    assert float(after["star_resolution"]) == 120000.0
