"""The residual clip: a running robust sigma the width of the high pass, and
NaN beyond nsig of it.

The statistic must describe the noise where it is (a noisier stretch is not a
stretch of outliers), must not be pulled by the outlier it is judging, and must
leave alone what is already NaN.
"""

import numpy as np

from pca2d import outliers


def test_the_running_sigma_follows_the_noise_where_it_is():
    r = np.random.default_rng(3)
    x = r.normal(0, 1, 6000)
    x[3000:] *= 3.0
    med, sig = outliers.running_stats(x, 151)
    assert abs(np.median(sig[300:2700]) - 1.0) < 0.1
    assert abs(np.median(sig[3300:5700]) - 3.0) < 0.3
    assert abs(np.median(med)) < 0.1


def test_a_spike_is_cut_and_the_noise_around_it_is_not():
    r = np.random.default_rng(5)
    x = r.normal(0, 1, 4000)
    x[[500, 1700, 3100]] = [9.0, -12.0, 8.0]
    x[2000:2040] = np.nan
    clipped, cut = outliers.clip(x, 151, 6.0)
    assert cut[[500, 1700, 3100]].all()
    assert cut.sum() == 3, "Gaussian noise never reaches 6 sigma in 4000 samples"
    assert np.isnan(clipped[[500, 1700, 3100]]).all()
    assert not cut[2000:2040].any() and np.isnan(clipped[2000:2040]).all()


def test_rows_are_clipped_each_against_its_own_noise():
    r = np.random.default_rng(7)
    x = r.normal(0, 1, (3, 2000)) * np.array([[1.0], [5.0], [0.2]])
    x[:, 1000] = 4.0
    _, cut = outliers.clip(x, 151, 3.5)
    assert cut[0, 1000] and not cut[1, 1000] and cut[2, 1000], (
        "4 is 4 sigma in the first row, under one in the second, 20 in the third")
