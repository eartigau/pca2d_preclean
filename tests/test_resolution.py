"""The star's spectra are smoothed to the resolution exactly as LBL smooths its templates."""
import numpy as np
import pytest

from pca2d.resolution import fwhm_samples, smooth, smooth_rows


def test_one_resolution_element_in_samples_is_lbl_s():
    assert fwhm_samples(70000, 0.5) == 9       # SPIRou
    assert fwhm_samples(80000, 0.5) == 7       # NIRPS


def _spectrum(n=3000, seed=3):
    r = np.random.default_rng(seed)
    x = np.arange(n, dtype=float)
    lines = -0.3 * sum(np.exp(-0.5 * ((x - c) / 6.0) ** 2) for c in range(100, n, 170))
    return lines + 0.01 * r.normal(size=n)


def test_it_is_lbl_s_filter_where_nothing_is_missing():
    lbl_math = pytest.importorskip("lbl.core.math")
    y = _spectrum()
    ours = smooth(y, 9)
    theirs = lbl_math.gaussian_weighted_savgol(y, 9, polyorder=3)
    inner = slice(20, -20)                     # LBL reflects at the ends, we do not
    assert np.allclose(ours[inner], theirs[inner], rtol=0, atol=1e-10)


def test_it_is_lbl_s_refit_beside_a_gap():
    lbl_math = pytest.importorskip("lbl.core.math")
    y = _spectrum()
    y[1000:1040] = 0.0                         # no support: zero here, NaN for LBL
    ours = smooth(y, 9)
    theirs = lbl_math.gaussian_weighted_savgol(np.where(y == 0, np.nan, y), 9,
                                               polyorder=3)
    both = np.isfinite(theirs) & (ours != 0)
    both[:20] = both[-20:] = False
    assert np.allclose(ours[both], theirs[both], rtol=0, atol=1e-10)
    assert not np.any(ours[1000:1040]), "the smoothing grew into the gap"


def test_a_pixel_spike_goes_and_a_resolved_line_stays():
    x = np.arange(2000, dtype=float)
    sigma = 2 * 9 / 2.3548                     # two resolution elements wide
    line = -0.3 * np.exp(-0.5 * ((x - 700) / sigma) ** 2)
    spike = np.zeros(2000)
    spike[1400] = 0.3
    out = smooth(line + spike + 1e-9, 9)       # no exact zero: all of it has support
    assert abs(out[1400]) < 0.1, "a one-pixel spike survived: %.3f" % out[1400]
    assert abs(out[700] - line[700]) < 0.05 * 0.3, "the line lost its depth"


def test_rows_come_back_orthonormal():
    basis = np.array([_spectrum(seed=s) for s in (1, 2)])
    out = smooth_rows(basis, 9)
    assert np.allclose(out @ out.T, np.eye(2), atol=1e-12)
