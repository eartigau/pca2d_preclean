"""Measuring the correction's amplitudes in the metric that decides velocities.

What a contaminant does to a radial velocity is its overlap with the DERIVATIVE
of the star, since that is what LBL projects a residual on. A contaminant flat
where the star has structure moves no line, however large it is in flux. These
pin the weighting that says so, and the derivative it is built on.
"""

from __future__ import annotations

import numpy as np

from pca2d.reconstruct import VELOCITY_FLOOR, star_derivative, velocity_weights
from pca2d.resolution import smooth


def a_line(x, centre, width=8.0, depth=1.0):
    return -depth * np.exp(-0.5 * ((x - centre) / width) ** 2)


def test_the_savgol_derivative_is_the_derivative():
    """deriv=1 of the template-making filter, against the analytic one. Never
    np.gradient: on noisy data the difference is not cosmetic."""
    x = np.arange(600.0)
    y = a_line(x, 300.0)
    got = smooth(y, 8.0, deriv=1)
    want = (x - 300.0) / 8.0 ** 2 * np.exp(-0.5 * ((x - 300.0) / 8.0) ** 2)
    core = slice(250, 350)
    assert np.max(np.abs(got[core] - want[core])) < 0.05 * np.max(np.abs(want))
    assert abs(got[300]) < 1e-9, "zero at the centre of a symmetric line"
    assert np.max(np.abs(got[:200])) < 1e-12 * np.max(np.abs(want)), \
        "and nothing at all where the line is not"


def test_the_weights_follow_the_star_and_never_invent_a_sample():
    x = np.arange(400.0)
    star = np.vstack([a_line(x, 200.0), a_line(x, 200.0)])
    w = np.ones_like(star)
    w[:, :20] = 0.0                       # what the fit gave no weight
    out = velocity_weights(w, star_derivative(star, 8.0))

    assert np.all(out[:, :20] == 0.0), "a sample with no weight keeps none"
    assert out.shape == w.shape
    wing = int(200 + 8)                   # where the derivative peaks
    assert out[0, wing] > 5 * out[0, 350], "the line wing outweighs the flat"
    assert out[0, 350] > 0, "but the flat is not silenced either"
    assert abs(out[0, 200] - VELOCITY_FLOOR) < 1e-3, \
        "the centre of a line carries no velocity information: the floor"
    kept = out[w > 0].sum() / w[w > 0].sum()
    assert abs(kept - (1.0 + VELOCITY_FLOOR)) < 0.05, \
        "and the weights keep the total they had, so the errors keep theirs"


def test_a_star_with_no_structure_leaves_the_weights_alone():
    """Nothing to say about velocity, so it says nothing rather than zero."""
    w = np.ones((2, 50))
    assert np.array_equal(velocity_weights(w, np.zeros((2, 50))), w)
    assert np.array_equal(velocity_weights(np.zeros((2, 50)),
                                           np.ones((2, 50))), np.zeros((2, 50)))


def test_the_metric_leaves_less_velocity_behind():
    """The point of the whole thing, measured as the thing it is for.

    Two observer components: one a line beside the star's, which moves it, and
    one a broad bump, which does not. What matters is not how well each
    amplitude is recovered but how much VELOCITY the correction leaves behind,
    which is the residual's projection on the star's derivative. Measured over
    sixty noise realisations, in both metrics.
    """
    rng = np.random.default_rng(11)
    x = np.arange(1000.0)
    star = np.vstack([a_line(x, 500.0, 8.0)] * 2)
    sharp = a_line(x, 503.0, 8.0, depth=0.05)          # moves the line
    broad = -0.05 * np.exp(-0.5 * ((x - 500.0) / 160.0) ** 2)   # does not
    basis = np.vstack([sharp, broad])
    truth = np.array([1.0, 1.0])

    deriv = star_derivative(star, 8.0)[0]
    scale = float(deriv @ deriv)
    left = {}
    for name, weights in (("flux", np.ones_like(star)),
                          ("velocity", velocity_weights(np.ones_like(star),
                                                        star_derivative(star, 8.0)))):
        shifts = []
        for _ in range(60):
            data = truth @ basis + rng.normal(0.0, 0.02, size=x.size)
            wr = weights[0]
            normal = (basis * wr) @ basis.T
            fitted = np.linalg.solve(normal, (basis * wr) @ data)
            residual = data - fitted @ basis
            shifts.append(float(residual @ deriv) / scale)
        left[name] = float(np.std(shifts))

    assert left["velocity"] < left["flux"], \
        "the metric that decides velocities leaves less velocity behind"
    assert left["velocity"] < 0.8 * left["flux"], \
        "and by enough to be worth measuring on real spectra (%.3g vs %.3g)" \
        % (left["velocity"], left["flux"])
