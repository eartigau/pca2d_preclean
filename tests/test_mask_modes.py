"""One set of lines for every epoch: correct --mask common, and the fit with no
star component that the refit used to break on."""
import numpy as np
import pytest
from astropy.table import Table

from pca2d.reconstruct import fit_weights_mask, star_support

N_ROWS, N_PIX = 8, 600


@pytest.fixture(scope="module")
def cube(tmp_path_factory):
    r = np.random.default_rng(3)
    path = tmp_path_factory.mktemp("cube")
    data = 0.01 * r.normal(size=(N_ROWS, N_PIX))
    sigma = np.full((N_ROWS, N_PIX), 0.01)
    # each exposure loses a different block: its own holes, nobody else's
    for i in range(N_ROWS):
        sigma[i, 100 + 40 * i:130 + 40 * i] = np.nan
    np.save(path / "grid.npy", 1500.0 * np.exp(np.arange(N_PIX) * 0.5 / 299792.458))
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", sigma)
    meta = Table()
    meta["filename"] = ["%04dt.fits" % (i // 2) for i in range(N_ROWS)]
    meta["bjd"] = 2459000.0 + np.arange(N_ROWS) * 1.5
    meta["berv"] = np.linspace(-10.0, 10.0, N_ROWS // 2).repeat(2)
    meta["airmass"] = np.full(N_ROWS, 1.2)
    meta["snr_band"] = np.full(N_ROWS, 100.0)
    meta["exposure"] = np.arange(N_ROWS) // 2
    meta["parity"] = np.arange(N_ROWS) % 2
    meta.write(path / "meta.fits", overwrite=True)
    return str(path)


def test_every_exposure_blanks_its_own_by_default(cube):
    own = fit_weights_mask(cube)
    masks = [m[0] for m in own.values()]
    assert not all(np.array_equal(masks[0], m) for m in masks[1:]), \
        "the exposures are meant to differ here"


def test_common_gives_every_exposure_the_same_set_of_lines(cube):
    own = fit_weights_mask(cube)
    shared = fit_weights_mask(cube, "common")
    first = next(iter(shared.values()))
    for name, mask in shared.items():
        assert np.array_equal(mask, first), "%s carries a different set" % name
        assert not (mask & ~own[name]).any(), "a sample nobody weighted is kept"
    for parity in (0, 1):
        want = np.ones(first.shape[1], dtype=bool)
        for mask in own.values():
            want &= mask[parity]
        assert np.array_equal(first[parity], want)
    assert first.sum() < sum(m.sum() for m in own.values()) / len(own), \
        "a common mask blanks more of each exposure than its own does"


def test_the_mask_option_decides_over_the_refit_s_own():
    """The refit computes a mask from the exposure's own weights. With `common`
    the option wins anyway, or it does nothing at all on a nightly-stacked
    target, where every file goes through that path (2026-09-12)."""
    from pca2d.reconstruct import mask_for
    own = np.array([[True, False], [True, True]])
    shared = np.array([[False, False], [True, False]])
    by_file = {"0001t.fits": shared, "0002t.fits": shared}
    assert mask_for("none", own, by_file, "0001t.fits") is None
    assert np.array_equal(mask_for("common", own, by_file, "0001t.fits"), shared)
    assert np.array_equal(mask_for("exposure", own, by_file, "0001t.fits"), own)
    assert np.array_equal(mask_for("exposure", None, by_file, "0001t.fits"), shared)
    assert mask_for("exposure", None, None, "0001t.fits") is None


def test_a_fit_with_no_star_component_solves_for_the_observer_alone():
    """n_star 0: the star basis is empty, and carrying it broke the refit of
    every exposure. The observer amplitudes must still come out right."""
    from pca2d.twoframe import LanczosShifter, joint_coeffs
    m, rows = 400, 2
    r = np.random.default_rng(0)
    Q = r.normal(size=(3, m))
    truth = np.array([2.0, -1.0, 0.5])
    data = np.tile(truth @ Q, (rows, 1)) + 0.01 * r.normal(size=(rows, m))
    shifter = LanczosShifter(m, a=8, max_shift=4)
    a, b, _, _ = joint_coeffs(data, np.ones((rows, m)), np.zeros((0, m)), Q,
                              shifter, np.zeros(rows),
                              exposure=np.zeros(rows, dtype=int))
    assert a.shape == (rows, 0)
    assert np.allclose(b[0], truth, atol=0.01)


def test_a_fit_with_no_star_component_has_no_support_to_guard():
    """n_star 0: P is empty, and np.any over its first axis gives a scalar,
    which broke gap_guard in the refit on 2026-09-12."""
    assert star_support({"P": np.zeros((0,))}) is None
    assert star_support({"P": np.zeros((0, 5))}) is None
    live = star_support({"P": np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 2.0]])})
    assert live.shape == (3,)
    assert list(live) == [False, True, True]
