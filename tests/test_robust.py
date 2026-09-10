"""Clipping, error propagation and the outlier cut.

These are the parts that decide which data is used, so a silent change here
moves every number downstream without touching a single equation.
"""

from __future__ import annotations

import numpy as np
import pytest

from pca2d.twoframe import (LanczosShifter, _nanmedian_over_rows,
                                 clip_weights, coefficient_errors,
                                 joint_coeffs, mad_outliers)


def test_nanmedian_matches_numpy(rng):
    """The bottleneck path is 6.3x faster and must be bit-identical to numpy.

    It is chosen at import time depending on whether bottleneck is installed,
    which is exactly the kind of difference that would otherwise show up as two
    machines disagreeing on the third decimal.
    """
    z = rng.normal(size=(40, 300))
    z[z > 2.0] = np.nan
    z[:, 17] = np.nan                       # a column that is entirely missing
    got = _nanmedian_over_rows(z)
    with np.errstate(invalid="ignore"):
        want = np.nanmedian(z, axis=0)
    assert np.array_equal(np.nan_to_num(got, nan=-999),
                          np.nan_to_num(want, nan=-999))


def test_clip_weights_only_removes(problem, rng):
    """Soft clipping may lower a weight; it may never raise one.

    The clip is applied to the ORIGINAL weights each iteration and never
    compounded, so a pixel cannot be down-weighted twice for the same reason.
    A weight that came back larger would mean the clip is being applied to its
    own output.
    """
    w0 = problem["w"]
    resid = rng.normal(size=w0.shape) * 0.1
    resid[3, 40:60] = 50.0                  # an obvious outlier stretch
    # returns (weights, hit fraction): the fraction of live samples beyond the
    # clip, which the run log reports so a suddenly noisy night is visible
    w, hit = clip_weights(w0, resid, clip=3.0)
    assert 0.0 <= hit <= 1.0
    assert w.shape == w0.shape
    assert np.all(w <= w0 + 1e-12), "clipping increased a weight"
    assert np.all(w >= 0.0), "clipping produced a negative weight"
    assert w[3, 40:60].sum() < w0[3, 40:60].sum(), "the outlier was not clipped"


def test_clip_never_hard_rejects(problem, rng):
    """A 6-sigma sample keeps a quarter of its weight, it is not thrown away.

    The factor outside the clip is (clip/|z|)^2, which is what you get by
    inflating an outlier's variance until it IS a clip-sigma sample. A hard cut
    on residual would sculpt the very variability the PCA exists to find, so a
    flare must survive at reduced weight rather than vanish.
    """
    w0 = np.ones((60, 40))
    resid = rng.normal(size=w0.shape) * 0.0
    resid += rng.normal(size=w0.shape)      # unit scatter
    resid[0, :] = 6.0                       # one row at 6 sigma everywhere
    w, _ = clip_weights(w0, resid, clip=3.0, min_spectra=10)
    assert np.all(w[0] > 0), "a 6 sigma row was hard-rejected"
    assert np.all(w[0] < w0[0]), "a 6 sigma row was not down-weighted at all"


def test_errors_are_positive_and_scale(problem):
    shifter = LanczosShifter(problem["data"].shape[1], threads=1)
    a, b, _alpha, _ = joint_coeffs(problem["data"], problem["w"], problem["P"],
                           problem["Q"], shifter, problem["delta"], chunk=8)
    formal, scaled, chi2_red, _vel = coefficient_errors(
        problem["data"], problem["w"], problem["P"], problem["Q"], shifter,
        problem["delta"], a, b, chunk=8)
    live = problem["w"].sum(axis=1) > 0
    assert np.all(formal[live] > 0), "a constrained row has a zero error bar"
    assert np.all(np.isfinite(chi2_red[live]))
    # the scaled error is the formal one times sqrt(reduced chi2), by definition
    assert np.allclose(scaled[live], formal[live]
                       * np.sqrt(chi2_red[live])[:, None], rtol=1e-9)


def test_mad_outliers_threshold_is_in_sigmas():
    """The cut is |c - median| / (1.4826 MAD) > threshold, i.e. in sigmas.

    Pinned deliberately. The option is called --max-sigma because of this, and
    the value was NOT changed when it was renamed from --max-mad: a run made
    before the rename must stay comparable with one made after.
    """
    n = 200
    c = np.zeros((n, 1))
    c[:, 0] = np.arange(n, dtype=float) - n / 2.0
    parity = np.zeros(n, dtype=int)
    centre = np.median(c[:, 0])
    mad_plain = np.median(np.abs(c[:, 0] - centre))
    # a point sitting at exactly 4 robust sigmas
    c[0, 0] = centre + 4.0 * 1.4826 * mad_plain
    assert mad_outliers(c, parity, threshold=3.0)[0], "4 sigma survived a 3 sigma cut"
    assert not mad_outliers(c, parity, threshold=5.0)[0], "4 sigma was cut at 5 sigma"


def test_mad_outliers_is_per_parity():
    """Two parities at different levels must not flag each other.

    Pooling them makes each component bimodal, inflates the MAD by the
    separation of the two modes, and hides every real outlier inside it.
    """
    n = 100
    c = np.zeros((n, 1))
    parity = np.tile([0, 1], n // 2)
    r = np.random.default_rng(3)
    c[parity == 0, 0] = r.normal(scale=1.0, size=(parity == 0).sum())
    c[parity == 1, 0] = 1000.0 + r.normal(scale=1.0, size=(parity == 1).sum())
    c[2, 0] = 50.0                          # a genuine outlier within parity 0
    flag = mad_outliers(c, parity, threshold=5.0)
    assert flag[2], "a real outlier was hidden by the parity separation"
    assert flag.sum() < n // 4, "the parity offset was mistaken for outliers"


def test_zero_threshold_disables():
    c = np.random.default_rng(0).normal(size=(50, 3))
    assert not mad_outliers(c, np.zeros(50, int), threshold=0.0).any()
