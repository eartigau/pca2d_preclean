"""The BERV bias, fitted by MCMC: that it finds one that is there, says
nothing more than a limit about one that is not, and that the sampler under
it samples what it is given."""

import numpy as np
import pytest

from pca2d import bervbias as bb


def a_campaign(amp, sigma, n=150, noise=5.0, seed=2):
    r = np.random.default_rng(seed)
    berv = r.uniform(-25, 25, n)
    e = np.full(n, noise)
    v = 30.0 + bb.shape(berv, amp, sigma) + r.normal(0, noise, n)
    return berv, v, e


def test_the_shape_is_odd_and_peaks_at_sigma():
    x = np.linspace(-30, 30, 6001)
    y = bb.shape(x, 2.0, 7.0)
    assert np.allclose(y, -y[::-1])
    assert x[np.argmax(y)] == pytest.approx(7.0, abs=0.02)
    assert y.max() == pytest.approx(bb.peak(2.0, 7.0))


def test_the_stretch_move_samples_a_known_gaussian():
    rng = np.random.default_rng(0)
    mean, width = np.array([1.0, -2.0]), np.array([0.5, 3.0])

    def log_prob(theta):
        return -0.5 * np.sum(((theta - mean) / width) ** 2, axis=1)

    chain, acceptance = bb.stretch(log_prob, rng.normal(size=(32, 2)), 3000,
                                   rng)
    flat = chain[1000:].reshape(-1, 2)
    assert np.allclose(flat.mean(axis=0), mean, atol=0.15)
    assert np.allclose(flat.std(axis=0), width, rtol=0.1)
    assert 0.2 < acceptance < 0.9


def test_a_bias_that_is_there_is_found_with_its_width():
    berv, v, e = a_campaign(amp=-10.0, sigma=6.0)
    fitted = bb.fit(berv, v, e, seed=3)
    assert fitted["detected"]
    assert fitted["sigma"][0] == pytest.approx(6.0, abs=1.5)
    assert fitted["peak"][0] == pytest.approx(bb.peak(-10.0, 6.0), abs=6.0)
    assert fitted["c"][0] == pytest.approx(30.0, abs=3.0)
    assert fitted["jitter"][0] < 4.0, "no scatter was added past the errors"
    assert "sigma)" in bb.summary(fitted)


def test_a_bias_that_is_not_there_is_an_upper_limit():
    berv, v, e = a_campaign(amp=0.0, sigma=6.0)
    fitted = bb.fit(berv, v, e, seed=4)
    assert not fitted["detected"]
    assert fitted["upper"] < 10.0
    assert bb.summary(fitted).startswith("none detected")


def test_scatter_past_the_error_bars_goes_into_the_jitter():
    berv, v, e = a_campaign(amp=0.0, sigma=6.0, noise=5.0)
    v = v + np.random.default_rng(5).normal(0, 20.0, v.size)
    fitted = bb.fit(berv, v, e, seed=6)
    assert fitted["jitter"][0] == pytest.approx(20.0, rel=0.25)
    assert not fitted["detected"], "scatter is not a bias"


def test_too_little_to_fit_is_said():
    assert bb.fit([1, 2, 3], [1, 2, 3], [1, 1, 1]) is None
    assert bb.summary(None) == "not fitted"


def test_the_envelope_holds_the_truth():
    berv, v, e = a_campaign(amp=-10.0, sigma=6.0)
    fitted = bb.fit(berv, v, e, seed=7)
    grid = np.linspace(-20, 20, 81)
    lo, mid, hi = bb.envelope(fitted, grid)
    truth = bb.shape(grid, -10.0, 6.0)
    inside = (truth >= lo - 2 * (hi - lo)) & (truth <= hi + 2 * (hi - lo))
    assert inside.mean() > 0.9


def test_amp_has_a_flat_prior():
    """With no BERV there is nothing for amp to change in the likelihood, so
    two amplitudes differ by the prior alone, and a flat prior makes them
    equal, on either sign, up to the bound and not beyond it."""
    berv = np.zeros(20)
    v = np.zeros(20)
    e = np.full(20, 2.0)
    theta = np.array([[0.0, np.log(5.0), 0.0, 0.0],
                      [1.0, np.log(5.0), 0.0, 0.0],
                      [3.0, np.log(5.0), 0.0, 0.0],
                      [-3.0, np.log(5.0), 0.0, 0.0],
                      [0.99 * bb.AMP_MAX, np.log(5.0), 0.0, 0.0],
                      [1.01 * bb.AMP_MAX, np.log(5.0), 0.0, 0.0]])
    lp = bb.log_probability(theta, berv, v, e)
    assert np.all(np.isfinite(lp[:5]))
    assert np.allclose(lp[:5], lp[0], rtol=0, atol=1e-12)
    assert lp[5] == -np.inf


def test_both_signs_are_explored_and_found():
    """A bias pulls either way: a positive one is found positive, a negative
    one negative, whichever side the walkers started on."""
    for amp, expect in ((10.0, 1.0), (-10.0, 0.0)):
        berv, v, e = a_campaign(amp=amp, sigma=6.0, seed=11)
        fitted = bb.fit(berv, v, e, seed=12)
        assert fitted["detected"]
        assert fitted["p_positive"] == pytest.approx(expect, abs=0.01)
        assert np.sign(fitted["peak"][0]) == np.sign(amp)
    berv, v, e = a_campaign(amp=0.0, sigma=6.0, seed=13)
    nothing = bb.fit(berv, v, e, seed=14)
    assert 0.05 < nothing["p_positive"] < 0.95, "no bias, no preferred sign"
