"""The card that says whether an observer-frame mean was divided out.

A corrected file's headers are read by whoever was not there when the run was
made, so a card that says True when nothing happened is worse than no card.
"""
import numpy as np

from pca2d.reconstruct import mean_divided

STAR = {"means": {"even": np.zeros(8), "odd": np.zeros(8)}}      # mean: star
OFFSET = {"means": {"even": np.zeros(8), "odd": np.r_[np.zeros(7), 0.01]}}


def test_the_nominal_keeps_no_observer_mean_so_the_card_is_false():
    assert mean_divided(STAR, 3) is False, \
        "mean: star leaves the observer mean at zero everywhere"


def test_a_fit_that_keeps_one_says_so():
    assert mean_divided(OFFSET, 3) is True
    assert mean_divided(OFFSET, 1) is True, "any component carries the mean out"


def test_a_file_that_divides_no_observer_component_divides_no_mean():
    assert mean_divided(OFFSET, 0) is False, \
        "t_k-0 keeps the observer block, mean included"
    assert mean_divided(STAR, 0) is False


def test_a_model_without_means_at_all_is_not_a_crash():
    assert mean_divided({}, 3) is False
    assert mean_divided({"means": {}}, 3) is False
