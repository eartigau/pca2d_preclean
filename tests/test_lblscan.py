"""The scan report's grouping and its numbers.

A scan moves one count at a time, so which run belongs to which scan is read
from the tags; and every number the report prints is the one the site prints,
from the same three functions.
"""

import numpy as np
import pytest

from pca2d import lblscan


def test_a_tag_is_star_then_observer_components():
    assert lblscan.parse_tag("2-3v") == (2, 3)
    assert lblscan.parse_tag("10-12") == (10, 12)
    with pytest.raises(SystemExit):
        lblscan.parse_tag("2x3")


def test_the_default_run_is_in_both_scans():
    tags = ["1-3v", "2-3v", "3-3v", "2-2v", "2-4v", "2-5v", "2-6v", "2-7v"]
    star, observer = lblscan.scans(tags, star_default=2, observer_default=3)
    assert star == ["1-3v", "2-3v", "3-3v"]
    assert observer == ["2-2v", "2-3v", "2-4v", "2-5v", "2-6v", "2-7v"]
    assert lblscan.scans(tags) == (star, observer), "the defaults are what most tags share"


def test_the_numbers_are_about_the_median_and_the_nights_are_weighted():
    t = np.array([100.1, 100.2, 101.1, 101.2])
    v = np.array([10.0, 12.0, 20.0, 22.0])
    e = np.array([1.0, 1.0, 1.0, 1.0])
    s = lblscan.velocity_stats(t, v, e)
    assert s["n"] == 4 and s["nights"] == 2
    assert s["rms"] == pytest.approx(np.std(v - np.median(v)))
    assert s["nightly_rms"] == pytest.approx(5.0), "nightly means 11 and 21"
    assert s["median_error"] == 1.0
