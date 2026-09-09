"""Which order a window belongs to, where two of them overlap.

Consecutive echelle orders share their ends, and the two parities are separate
cube rows, so an overlapping window is drawn twice unless something picks. The
rule is the order whose middle the window sits closest to, in units of that
order's own half-width; the parity follows from that order's index, and an
off-by-one there sends every overlapping window to the wrong side.
"""

import numpy as np

from pca2d.plotting import parity_from_bounds


# three orders, each 100 nm wide, overlapping their neighbours by 20 nm
LO = np.array([1000.0, 1080.0, 1160.0])
HI = np.array([1100.0, 1180.0, 1260.0])


def test_outside_every_order_declines_to_choose():
    parity, offset, n = parity_from_bounds(900.0, LO, HI)
    assert parity is None and n == 0 and not np.isfinite(offset)


def test_a_window_in_one_order_only_reports_that_it_had_no_choice():
    parity, offset, n = parity_from_bounds(1050.0, LO, HI)
    assert (parity, n) == (0, 1)
    assert offset == 0.0                      # dead centre of order 0


def test_an_overlap_keeps_the_order_whose_centre_is_nearer():
    # 1085 lies in order 0 (edge, offset 0.70) and order 1 (offset 0.90)
    parity, offset, n = parity_from_bounds(1085.0, LO, HI)
    assert n == 2
    assert parity == 0                        # order 0, so even
    assert np.isclose(offset, 0.70)

    # 1095 has moved past order 1's centre-ward side: 0.90 against 0.70
    parity, offset, n = parity_from_bounds(1095.0, LO, HI)
    assert (n, parity) == (2, 1)              # order 1, so odd
    assert np.isclose(offset, 0.70)


def test_the_parity_is_the_order_index_not_its_rank_among_the_covering_ones():
    # only orders 1 and 2 reach 1175, at 0.90 and 0.70 of a half-width; the
    # winner is order 2, parity 0, and a rule that numbered the covering pair
    # 0 and 1 rather than using the true index would answer 1
    parity, offset, n = parity_from_bounds(1175.0, LO, HI)
    assert (n, parity) == (2, 0)
    assert np.isclose(offset, 0.70)


def test_orders_of_unequal_width_compare_by_their_own_half_widths():
    lo = np.array([1000.0, 1090.0])
    hi = np.array([1100.0, 1110.0])           # a 100 nm order and a 20 nm one
    # 1095 is 45 nm from the wide order's centre out of 50, and 5 nm from the
    # narrow one's out of 10: 0.90 against 0.50, so the narrow order wins
    parity, offset, n = parity_from_bounds(1095.0, lo, hi)
    assert (n, parity) == (2, 1)
    assert np.isclose(offset, 0.50)
