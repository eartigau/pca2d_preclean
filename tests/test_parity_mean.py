"""A corrected file loses what panel 3 of the sequence figure shows removed.

Panel 3 is the data minus the parity mean the fit subtracted before solving,
minus the observer block. The file used to lose only the observer block, and
kept the mean: a static observer-frame pattern, slanted stripes in the star's
frame at 1267 nm on TOI-2120 (2026-09-10).
"""

import numpy as np

from pca2d.reconstruct import order_correction


def _model(n=400):
    grid = np.exp(np.linspace(np.log(1260.0), np.log(1270.0), n))
    x = np.linspace(0.0, 8.0 * np.pi, n)
    even, odd = 0.03 * np.sin(x), 0.02 * np.cos(x)
    even[:20] = odd[:20] = 0.0                 # no basis support at the start
    return {"grid": grid, "means": {"even": even, "odd": odd}}, grid


def test_the_parity_mean_comes_out_with_the_observer_block():
    model, grid = _model()
    wave = grid[50:350:3] * (1 + 1e-7)
    correction = np.full(grid.size, 0.01)
    for order, name in ((18, "even"), (19, "odd")):
        values, live = order_correction(model, correction, 3, order, wave)
        want = np.interp(wave, grid, correction + model["means"][name])
        assert live.all()
        assert np.allclose(values, want, atol=2e-4), name


def test_a_file_that_removes_no_observer_component_keeps_the_mean():
    model, grid = _model()
    wave = grid[50:350:3]
    correction = np.full(grid.size, 0.01)
    values, _ = order_correction(model, correction, 0, 18, wave)
    assert np.allclose(values, 0.01, atol=1e-9)


def test_nothing_is_applied_where_the_basis_has_no_support():
    model, grid = _model()
    wave = grid[2:60]
    values, live = order_correction(model, np.zeros(grid.size), 3, 18, wave)
    assert not live[wave < grid[20]].any()
    assert live[wave > grid[22]].all()


def test_an_order_that_misses_the_grid_is_skipped():
    model, grid = _model()
    assert order_correction(model, np.zeros(grid.size), 3, 18,
                            np.full(10, np.nan)) is None


def test_a_pixel_is_blanked_where_the_fit_gave_no_weight():
    """As panel 3 hides it: both grid samples around a pixel must be weighted."""
    from pca2d.reconstruct import weighted_on_pixels

    grid = np.linspace(1000.0, 1010.0, 101)
    alive = np.ones(grid.size, dtype=bool)
    alive[40:45] = False
    wave = np.array([1001.05, 1003.95, 1004.05, 1004.25, 1005.05, 999.0, 1011.0])
    got = weighted_on_pixels(alive, wave, grid)
    # the hole is samples 40-44, 1004.0-1004.4 nm: 1003.95 touches sample 40,
    # 1004.05 and 1004.25 are inside it, 1005.05 sits between 50 and 51, and
    # outside the grid there is no panel 3
    assert got.tolist() == [True, False, False, False, True, True, True]
