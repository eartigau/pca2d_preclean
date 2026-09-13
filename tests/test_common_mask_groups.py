"""The common mask is taken within a campaign and within an order parity.

Two ways it was taken over too much, both found on the joint run of 2026-09-13:

  the figure asked a sample to be alive in EVERY row, including rows of the other
  order parity, which never cover that wavelength at all since orders n and n+2
  do not overlap. The intersection was empty everywhere the parities do not meet,
  so the three H-band windows vanished from the report and the J-band ones kept
  only the 15-20% of columns two orders reach.

  the correction intersected over every object of a joint cube, so a sample lost
  by one night of the faintest star was blanked in every exposure of the
  brightest: 9.9% of the H band NaN became 17.1%.
"""
import numpy as np

from pca2d.plotting import live_mask


def parity_weights(n_rows=8, n_cols=40):
    """Weights where each parity covers its own half of the grid, as orders do.

    Row parity alternates; even rows carry the first half, odd rows the second,
    and the two meet over four columns in the middle, which is the overlap of
    two consecutive orders.
    """
    w = np.zeros((n_rows, n_cols), dtype=np.float32)
    half = n_cols // 2
    w[0::2, :half + 2] = 1.0
    w[1::2, half - 2:] = 1.0
    return w, np.arange(n_rows) % 2


def test_over_every_row_the_intersection_is_only_where_orders_overlap():
    """What the figure used to do, kept here as the thing to avoid."""
    w, _parity = parity_weights()
    shared = live_mask(w).all(axis=0)
    assert shared.sum() == 4, "only the four columns both parities reach"


def test_within_a_parity_the_whole_coverage_survives():
    w, parity = parity_weights()
    live = live_mask(w)
    per_group = {g: live[parity == g].all(axis=0) for g in np.unique(parity)}
    assert per_group[0].sum() == 22 and per_group[1].sum() == 22
    assert not (per_group[0] & ~live[0]).any(), "nothing claimed alive that is not"
    # every row keeps what its own parity covers, which is the point
    for i, g in enumerate(parity):
        assert (per_group[g] <= live[i]).all()


def test_one_bad_exposure_costs_its_own_group_and_no_other():
    w, parity = parity_weights()
    w[2, 5:9] = 0.0                       # one even row loses four samples
    live = live_mask(w)
    per_group = {g: live[parity == g].all(axis=0) for g in np.unique(parity)}
    assert per_group[0][5:9].sum() == 0, "its own parity pays for it"
    assert per_group[1].sum() == 22, "the other parity does not"


def test_the_correction_intersects_within_an_object(tmp_path):
    """A joint cube holds several campaigns, and each keeps its own set of
    lines: this is the data, not a figure."""
    from astropy.table import Table

    import pca2d.reconstruct as rec

    n_cols = 12
    w = np.ones((6, n_cols), dtype=np.float32)
    w[0, 3] = 0.0        # PROXIMA loses column 3
    w[4, 7] = 0.0        # GJ3090 loses column 7
    meta = Table({"filename": ["p1", "p1", "p2", "p2", "g1", "g1"],
                  "object": ["PROXIMA", "PROXIMA", "PROXIMA", "PROXIMA",
                             "GJ3090", "GJ3090"],
                  "parity": [0, 1, 0, 1, 0, 1]})

    class Fake:
        @staticmethod
        def load_cube(cube, dtype=None):
            return None, None, w, meta

        @staticmethod
        def row_parity(meta, n):
            return np.asarray(meta["parity"])

    old = rec._bcd
    rec._bcd = Fake
    try:
        out = rec.fit_weights_mask("cube", mode="common")
    finally:
        rec._bcd = old
    assert not out["p1"][0][3], "Proxima's own loss blanks Proxima"
    assert out["g1"][0][3], "and does NOT blank GJ 3090"
    assert not out["g1"][0][7] and out["p1"][0][7], "nor the other way round"
    assert (out["p1"] == out["p2"]).all(), "one set of lines per campaign"
