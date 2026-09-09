"""Splitting the observer-frame mean into what to remove and what to keep.

The mean over the exposures, taken per order parity, holds two different
things. The even-minus-odd half-difference is instrumental, the order overlaps
not being at the same resolution, and it has to come out before the fit or it
takes the whole of PC1 (NOTES.md 13.9). The part the two parities share is
whatever is static in the OBSERVER'S frame, which is the atmosphere, and it has
to stay in: the correction removes components and never means, so anything
subtracted here is never removed from a spectrum again.

These pin the arithmetic of that split, which is the whole of the change.
"""

import numpy as np

from pca2d.twoframe import mean_rows, parity_means, subtract_means


def build(n_pixels=64, seed=3):
    rng = np.random.default_rng(seed)
    shared = rng.normal(size=n_pixels)          # the atmosphere, both parities
    offset = rng.normal(size=n_pixels)          # instrumental, equal and opposite
    parity = np.array([0, 1, 0, 1, 0, 1])
    data = np.empty((parity.size, n_pixels))
    for i, p in enumerate(parity):
        data[i] = shared + (offset if p == 0 else -offset)
    return data, np.ones_like(data), parity, shared, offset


def test_the_offset_half_is_equal_and_opposite():
    data, w, parity, shared, offset = build()
    means, _groups = parity_means(data, w, parity)
    common = means.mean(axis=0)
    half = means - common
    assert np.allclose(common, shared)
    assert np.allclose(half[0], offset)
    assert np.allclose(half[0], -half[1])


def test_removing_only_the_offset_leaves_the_shared_part_in_the_data():
    data, w, parity, shared, _offset = build()
    means, groups = parity_means(data, w, parity)
    group = np.searchsorted(groups, parity)
    subtract_means(data, means - means.mean(axis=0), group)
    # every row now carries the shared part and nothing else: the two parities
    # have been brought onto one another, and the atmosphere is still there
    for row in data:
        assert np.allclose(row, shared)


def test_removing_the_whole_mean_takes_the_shared_part_with_it():
    data, w, parity, _shared, _offset = build()
    means, groups = parity_means(data, w, parity)
    group = np.searchsorted(groups, parity)
    subtract_means(data, means, group)
    # nothing left at all: this is what put the atmosphere beyond the model
    assert np.allclose(data, 0.0)


def test_the_two_modes_differ_by_exactly_the_shared_part():
    data, w, parity, shared, _offset = build()
    means, groups = parity_means(data, w, parity)
    group = np.searchsorted(groups, parity)
    full = mean_rows(means, group, parity.size)
    offset_only = mean_rows(means - means.mean(axis=0), group, parity.size)
    assert np.allclose(full - offset_only, shared[None, :])
