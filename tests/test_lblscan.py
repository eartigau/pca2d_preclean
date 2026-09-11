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


def test_the_scans_are_read_from_the_tags():
    tags = ["2-7v", "1-3v", "2-3v", "3-3v", "2-2v", "1-2v", "1-4v", "2-4v"]
    groups = lblscan.scan_groups(tags)
    assert ("observer", 1, ["1-2v", "1-3v", "1-4v"]) in groups
    assert ("observer", 2, ["2-2v", "2-3v", "2-4v", "2-7v"]) in groups
    assert ("star", 3, ["1-3v", "2-3v", "3-3v"]) in groups
    assert not any(kind == "star" and fixed in (2, 4) for kind, fixed, _ in groups), (
        "two runs at one observer count are already in the observer scans")
    assert [g[0] for g in groups] == ["observer", "observer", "star"]


def _fake_run(tag, rms, seed):
    r = np.random.default_rng(seed)
    t = np.sort(59900 + r.uniform(0, 400, 60))
    v = r.normal(0, rms, t.size)
    e = np.full(t.size, 7.0)
    n_star, n_earth = lblscan.parse_tag(tag)
    return {"tag": tag, "n_star": n_star, "n_earth": n_earth, "t": t, "v": v, "e": e,
            "stats": lblscan.velocity_stats(t, v, e)}


def test_the_pages_of_a_partial_matrix_are_drawn():
    """The counts, the matrix and one page per scan, on a grid with holes in it."""
    import matplotlib.pyplot as plt
    tags = ["1-2v", "1-3v", "1-4v", "2-2v", "2-3v", "2-4v", "3-3v"]
    runs = [_fake_run(t, 20 + 5 * i, i) for i, t in enumerate(tags)]
    delivered = _fake_run("0-0", 47, 99)
    by_tag = {r["tag"]: r for r in runs}
    groups = [(kind, fixed, [by_tag[t] for t in members])
              for kind, fixed, members in lblscan.scan_groups(tags)]
    assert len(groups) == 3, "two observer scans and one star scan"
    figures = [lblscan.metrics_page(delivered, groups), lblscan.matrix_page(delivered, runs)]
    figures += [lblscan.series_page(delivered, g[2], lblscan.group_name(g[0], g[1]))
                for g in groups]
    for fig in figures:
        assert fig is not None
        plt.close(fig)
    assert lblscan.matrix_page(delivered, runs[:3]) is None, "one star count is not a matrix"


def test_the_numbers_are_about_the_median_and_the_nights_are_weighted():
    t = np.array([100.1, 100.2, 101.1, 101.2])
    v = np.array([10.0, 12.0, 20.0, 22.0])
    e = np.array([1.0, 1.0, 1.0, 1.0])
    s = lblscan.velocity_stats(t, v, e)
    assert s["n"] == 4 and s["nights"] == 2
    assert s["rms"] == pytest.approx(np.std(v - np.median(v)))
    assert s["nightly_rms"] == pytest.approx(5.0), "nightly means 11 and 21"
    assert s["median_error"] == 1.0
