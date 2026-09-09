"""A valid sample with a gap close on both sides is not a measurement.

It passed the quality cuts, but its neighbourhood did not: its high pass came
from a filter window mostly filled by interpolation, and a shift into another
frame drew it from samples that are not there. On an image those are the
speckles of colour inside a band of yellow, and they read as data.

The three cases below are the ones the rule was specified with, verbatim.
"""

import numpy as np

from pca2d.preprocess import drop_isolated

NAN = np.nan


def run(values, window=3):
    v = np.asarray(values, dtype=float)
    return np.where(drop_isolated(np.isfinite(v), window), v, np.nan)


def same(got, want):
    got, want = np.asarray(got, float), np.asarray(want, float)
    return np.array_equal(np.isnan(got), np.isnan(want)) and \
        np.allclose(got[~np.isnan(got)], want[~np.isnan(want)])


def test_a_single_sample_between_two_gaps_goes():
    assert same(run([1, NAN, 5, NAN, 2]), [1, NAN, NAN, NAN, 2])


def test_a_pair_between_two_gaps_goes_as_well():
    assert same(run([1, NAN, 3, 5, NAN, 2]),
                [1, NAN, NAN, NAN, NAN, 2])


def test_a_run_of_six_is_long_enough_to_stay():
    v = [1, NAN, 5, 5, 5, 5, 3, 5, NAN, 2]
    assert same(run(v), v)


def test_the_ends_of_the_array_are_not_gaps():
    # the 1 has a gap on its right and nothing on its left, so it stays; the
    # rule needs a gap on BOTH sides
    assert same(run([1, NAN, 9, 9, 9, 9, 9]), [1, NAN, 9, 9, 9, 9, 9])


def test_it_acts_along_the_last_axis_of_a_stack():
    v = np.array([[1, NAN, 5, NAN, 2],
                  [1, 2, 3, 4, 5]], dtype=float)
    got = run(v)
    assert same(got[0], [1, NAN, NAN, NAN, 2])
    assert same(got[1], [1, 2, 3, 4, 5])


def test_a_window_of_one_only_kills_a_lone_sample():
    assert same(run([1, NAN, 3, 5, NAN, 2], window=1),
                [1, NAN, 3, 5, NAN, 2])
    assert same(run([1, NAN, 5, NAN, 2], window=1), [1, NAN, NAN, NAN, 2])


def test_nothing_valid_is_created():
    rng = np.random.default_rng(0)
    v = rng.normal(size=(4, 200))
    v[rng.random(v.shape) < 0.2] = np.nan
    before = np.isfinite(v)
    after = drop_isolated(before, 3)
    assert not (after & ~before).any()
