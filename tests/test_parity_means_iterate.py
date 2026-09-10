"""The static part of the model: one mean per parity in each frame, iterated.

A cube built from known pieces: a star with one line width per parity (the
even and the odd orders' resolutions), a static observer-frame pattern per
parity, and a telluric pattern whose depth varies from exposure to exposure
around a non-zero average. The fit must hand the star to the star-frame means
and the static sky to the observer-frame ones, and keep the star smeared over
the BERV out of the latter: estimated once, before the fit, as it used to be,
the observer mean held that smeared star and the observer block spent itself
cancelling it.
"""

import numpy as np
import pytest
from astropy.table import Table

from pca2d.grids import pixel_shift
from pca2d.twoframe import LanczosShifter, main

N_EXP, N_PIX, DV = 24, 2400, 0.5


def _lines(x, centres, width, depth):
    return -depth * sum(np.exp(-0.5 * ((x - c) / width) ** 2) for c in centres)


def _corr(a, b, cut=slice(150, -150)):
    a, b = a[cut] - a[cut].mean(), b[cut] - b[cut].mean()
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


@pytest.fixture(scope="module")
def truth(tmp_path_factory):
    r = np.random.default_rng(7)
    x = np.arange(N_PIX, dtype=float)
    grid = 1500.0 * np.exp(x * DV / 299792.458)
    star = [_lines(x, (300, 900, 1500, 2050), width, 0.4) for width in (4.0, 5.5)]
    sky = [_lines(x, (600, 1200, 1800), width, 0.3) for width in (3.0, 4.0)]
    tell = _lines(x, (700, 1350, 2150), 3.5, 0.2)
    berv = np.linspace(-15.0, 15.0, N_EXP)
    depth = 1.0 + 0.4 * r.normal(size=N_EXP)
    # the star carried into each exposure exactly as the fit carries it
    shifter = LanczosShifter(N_PIX, a=8, max_shift=40)
    placed = [shifter.rows(np.tile(s, (N_EXP, 1)), -pixel_shift(berv, DV))
              for s in star]
    data = np.array([placed[p][e] + sky[p] + depth[e] * tell
                     for e in range(N_EXP) for p in (0, 1)])
    data += 0.003 * r.normal(size=data.shape)
    path = tmp_path_factory.mktemp("cube")
    np.save(path / "grid.npy", grid)
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", np.full_like(data, 0.003))
    meta = Table()
    meta["filename"] = ["%04dt.fits" % e for e in range(N_EXP) for _ in (0, 1)]
    meta["bjd"] = 2459000.0 + np.repeat(np.arange(N_EXP), 2) * 1.3
    meta["berv"] = np.repeat(berv, 2)
    meta["airmass"] = np.full(2 * N_EXP, 1.2)
    meta["snr_band"] = np.full(2 * N_EXP, 100.0)
    meta["exposure"] = np.repeat(np.arange(N_EXP), 2)
    meta["parity"] = np.tile([0, 1], N_EXP)
    meta.write(path / "meta.fits", overwrite=True)
    smeared = [placed[p].mean(axis=0) for p in (0, 1)]
    return {"cube": str(path), "star": star, "sky": sky, "tell": tell,
            "depth": depth, "smeared": smeared}


def _fit(truth, out, *extra):
    main(["--cube", truth["cube"], "--outdir", str(out), "--iters", "8",
          "-k", "1", "-j", "1", "--max-mad", "0", *extra])
    return np.load(out / "fit.npz")


@pytest.fixture(scope="module")
def iterated(truth, tmp_path_factory):
    out = tmp_path_factory.mktemp("iterate")
    return _fit(truth, out, "--mean", "iterate"), out


def test_each_frame_gets_its_own_static_content(truth, iterated):
    fit, _ = iterated
    templates, means = fit["templates"], fit["means"]
    for p in (0, 1):
        static_sky = truth["sky"][p] + truth["depth"].mean() * truth["tell"]
        assert _corr(templates[p], truth["star"][p]) > 0.95, p
        assert _corr(means[p], static_sky) > 0.95, p


def test_the_smeared_star_stays_out_of_the_observer_mean(truth, iterated):
    fit, _ = iterated
    for p in (0, 1):
        assert abs(_corr(fit["means"][p], truth["smeared"][p])) < 0.15, p


def test_each_parity_keeps_its_own_resolution(truth, iterated):
    """The odd orders' lines are wider here, as a lower resolution makes them."""
    fit, _ = iterated
    templates = fit["templates"]
    for p in (0, 1):
        mine = _corr(templates[p], truth["star"][p])
        other = _corr(templates[p], truth["star"][1 - p])
        assert mine > other, (p, mine, other)


def test_the_observer_component_follows_the_telluric_depth(truth, iterated):
    """And not the BERV: when it follows the BERV it is describing the star."""
    fit, _ = iterated
    b = fit["b"][0::2, 0]                 # one row per exposure
    assert abs(np.corrcoef(b, truth["depth"])[0, 1]) > 0.95


def test_the_components_hold_no_static_content(iterated):
    fit, _ = iterated
    b = fit["b"][:, 0]
    assert abs(b.mean()) < 0.1 * b.std(), (b.mean(), b.std())


def test_the_components_file_carries_a_template_per_parity(iterated):
    from pca2d.reconstruct import load_model, template_for

    fit, out = iterated
    model = load_model(str(out / "twoframe_components.fits"))
    assert set(model["templates"]) == {"even", "odd"}
    for p in (0, 1):
        assert np.allclose(template_for(model, p), fit["templates"][p])
    assert np.allclose(model["means"]["even"], fit["means"][0])


def test_the_older_one_shot_mode_still_runs(truth, tmp_path):
    fit = _fit(truth, tmp_path / "offset", "--mean", "offset")
    assert "templates" not in fit.files
    assert str(fit["mean_mode"]) == "offset"
