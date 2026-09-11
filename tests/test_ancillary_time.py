"""Every quantity of the correlation matrix drawn against time, each exposure
marked by what the fit did with it."""
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.table import Table

from pca2d.twoframe import main

N_EXP, N_PIX, DV = 12, 1500, 0.5


@pytest.fixture(scope="module")
def fitted(tmp_path_factory):
    r = np.random.default_rng(7)
    n = 2 * N_EXP
    x = np.arange(N_PIX, dtype=float)
    star = -0.3 * sum(np.exp(-0.5 * ((x - c) / 5.0) ** 2) for c in (300, 800, 1200))
    earth = -0.2 * sum(np.exp(-0.5 * ((x - c) / 4.0) ** 2) for c in (500, 1000))
    berv = np.linspace(-12.0, 12.0, N_EXP).repeat(2)
    data = np.zeros((n, N_PIX))
    for i in range(n):
        data[i] = np.roll(star, int(round(-berv[i] / DV))) + earth * (1 + 0.3 * r.normal())
    data += 0.01 * r.normal(size=data.shape)
    path = tmp_path_factory.mktemp("cube")
    np.save(path / "grid.npy", 1500.0 * np.exp(np.arange(N_PIX) * DV / 299792.458))
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", np.full_like(data, 0.01))
    step = np.arange(n) // 2
    meta = Table()
    meta["filename"] = ["%04dt.fits" % k for k in step]
    meta["bjd"] = 2459000.0 + step * 1.7
    meta["mjdmid"] = meta["bjd"] - 2400000.5
    meta["berv"] = berv
    meta["airmass"] = 1.0 + 0.05 * step
    snr = np.full(n, 100.0)
    snr[6:8] = 20.0                    # exposure 3, below half the median
    meta["snr_band"] = snr
    meta["seeing"] = 0.8 + 0.02 * step
    meta["exposure"] = step
    meta["parity"] = np.arange(n) % 2
    meta.write(path / "meta.fits", overwrite=True)
    out = tmp_path_factory.mktemp("fit")
    main(["--cube", str(path), "--outdir", str(out), "--iters", "2",
          "-k", "1", "-j", "1", "--max-mad", "0"])
    return str(path), out


def test_the_quantities_are_the_matrix_s_and_every_exposure_is_there(fitted):
    from pca2d.figures.ancillary_time import exposure_table
    cube, out = fitted
    fit = np.load(out / "fit.npz")
    meta = Table.read(os.path.join(cube, "meta.fits"))
    rjd, values, status = exposure_table(meta, fit)
    assert list(values) == [str(v) for v in fit["anc_labels"] if str(v) != "time"]
    assert len(rjd) == N_EXP, "one point per exposure, not per parity row"
    assert list(np.flatnonzero(status == "not fitted")) == [3], \
        "the exposure the SNR cut dropped is shown, and marked"
    assert np.allclose(values["BERV"], np.linspace(-12.0, 12.0, N_EXP))
    assert np.allclose(rjd, 59000.0 + np.arange(N_EXP) * 1.7)


def test_the_page_draws_as_the_bundle_runs_it(fitted, tmp_path):
    import pca2d.figures.ancillary_time as at
    cube, out = fitted
    pdf = tmp_path / "ancillary.pdf"
    r = subprocess.run([sys.executable, at.__file__, "--cube", cube,
                        "--fit", str(out / "fit.npz"), "--title", "TEST",
                        "--out", str(pdf)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    assert pdf.exists()
