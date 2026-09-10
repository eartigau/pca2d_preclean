"""A corrected file holds what panel 3 of the sequence figure shows.

Panel 3, "minus the OBSERVER block: what a corrected file holds", is the cube's
data minus the fit's per-parity mean minus the observer block, with every
sample the fit gave no weight to hidden. On 2026-09-10 the files were found to
keep the parity mean the panel removes (slanted stripes at 1267 nm on TOI2120),
and the panel was made the reference for the corrected spectra. This pins the
two together: what the correction divides out of an order is exactly what
panel 3 subtracts from that order's parity, and the samples it blanks are
exactly the ones the fit did not weight.
"""

import numpy as np
import pytest
from astropy.table import Table

from pca2d.plotting import live_mask
from pca2d.reconstruct import (correction_on_grid, fit_weights_mask, load_model,
                               order_correction)
from pca2d.twoframe import fit_means, load_cube, main, mean_rows, row_parity

N_ROWS, N_PIX, DV = 12, 2000, 0.5


@pytest.fixture(scope="module")
def fitted(tmp_path_factory):
    r = np.random.default_rng(4)
    x = np.arange(N_PIX, dtype=float)
    grid = 1500.0 * np.exp(np.arange(N_PIX) * DV / 299792.458)
    star = -0.3 * sum(np.exp(-0.5 * ((x - c) / 5.0) ** 2) for c in (300, 900, 1500))
    earth = -0.2 * sum(np.exp(-0.5 * ((x - c) / 4.0) ** 2) for c in (600, 1200))
    berv = np.linspace(-12.0, 12.0, N_ROWS // 2).repeat(2)
    data = np.zeros((N_ROWS, N_PIX))
    for i in range(N_ROWS):
        data[i] = np.roll(star, int(round(-berv[i] / DV)))
        data[i] += earth * (1.0 + 0.3 * r.normal()) + 0.02 * (i % 2)
    data += 0.01 * r.normal(size=data.shape)
    data[3, 700:720] = np.nan              # a hole the fit cannot weight
    path = tmp_path_factory.mktemp("cube")
    np.save(path / "grid.npy", grid)
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", np.full_like(data, 0.01))
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
          "-k", "1", "-j", "2", "--max-mad", "0"])
    return str(path), out


def test_what_is_divided_out_is_what_panel_3_subtracts(fitted):
    cube, out = fitted
    fit = np.load(out / "fit.npz")
    model = load_model(str(out / "twoframe_components.fits"))
    _, _, _, meta = load_cube(cube, columns=np.arange(1))
    n, grid = len(meta), model["grid"]
    parity = row_parity(meta, n)
    means, group = fit_means(fit, meta, n, grid.size)
    # sequence.window_arrays, panel 3: data - offset - earth
    subtracted = np.asarray(mean_rows(means, group, n)) + fit["b"] @ fit["Q"]
    coeffs = {str(row["filename"]): row for row in model["coeffs"]}
    inner = slice(100, -100)
    for i in range(n):
        row = coeffs[str(meta["filename"][i])]
        correction, _, j = correction_on_grid(model, row, 0, model["n_earth"])
        values, live = order_correction(model, correction, j, int(parity[i]),
                                        grid[inner])
        assert live.any()
        assert np.allclose(values[live], subtracted[i][inner][live],
                           rtol=0, atol=1e-9), i


def test_the_samples_blanked_are_the_samples_the_fit_did_not_weight(fitted):
    cube, _ = fitted
    alive = fit_weights_mask(cube)
    _, _, w0, meta = load_cube(cube)
    parity = row_parity(meta, len(meta))
    for i in range(len(meta)):
        want = live_mask(w0[i:i + 1])[0]
        assert np.array_equal(alive[str(meta["filename"][i])][int(parity[i]) % 2],
                              want), i
    assert not alive["0001t.fits"][1][700:720].any(), "the hole must be blanked"
