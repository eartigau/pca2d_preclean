"""Two things a masked OH core must not do: bleed into the samples beside it,
and open a hole in the star's frame that is only a hole in the observer's.

The second made the vertical lines of the report's panels 2, 4 and 5: the gap
guard took the observer frame's support for the star frame's, cut a hole into
every row at the star column equal to the OH line's observer column, and the
star basis grew a spike at the hole's edge.
"""

import numpy as np

from pca2d.preprocess import erode_edges
from pca2d.twoframe import LanczosShifter, gap_guard, star_support


def test_the_last_valid_pixel_beside_a_gap_goes_and_the_ends_are_not_gaps():
    good = np.array([1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1], dtype=bool)
    assert erode_edges(good, 1).tolist() == [1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1]
    assert erode_edges(good, 0).tolist() == good.tolist()
    assert erode_edges(np.ones(5, dtype=bool), 2).all(), "no gap, nothing to erode"


def _rows(n=24, m=400, hole=200):
    w = np.ones((n, m))
    w[:, hole] = 0.0                    # an OH core, at one observer column in every row
    delta = np.linspace(-15.0, 15.0, n)
    return w, delta


def test_an_observer_hole_is_not_a_star_hole():
    w, delta = _rows()
    shifter = LanczosShifter(w.shape[1], a=8, max_shift=17)
    live = star_support(w, delta, shifter)
    assert live[150:250].all(), "the exposures at other BERVs see that star column"
    guarded = gap_guard(w.copy(), delta, 8, verbose=False, live=live)
    old = gap_guard(w.copy(), delta, 8, verbose=False)
    near = slice(150, 250)
    dropped_new = int(np.sum(guarded[:, near] == 0)) - w.shape[0]    # the core itself
    dropped_old = int(np.sum(old[:, near] == 0)) - w.shape[0]
    assert dropped_new == 0, "nothing beside the core is taken from any row"
    assert dropped_old >= 15 * w.shape[0], "the old rule cut a 2a hole into every row"


def test_a_true_star_hole_is_still_guarded():
    w, delta = _rows()
    w[:, 280:330] = 0.0                 # wider than any shift: a hole in both frames
    shifter = LanczosShifter(w.shape[1], a=8, max_shift=17)
    live = star_support(w, delta, shifter)
    assert not live[300:310].any()
    guarded = gap_guard(w.copy(), delta, 8, verbose=False, live=live)
    assert (guarded[:, 270:280] == 0).any(), "samples whose taps reach it are dropped"
