"""The measures that do not reward eating astrophysics.

Measured on TOI-1452, whose binary companion gives 412 m/s of curvature over one
season: the correction absorbed 45% of that curve, the raw rms fell from 100 to
68 m/s, and the scatter inside a night and the error per exposure both got worse.
So `rms` alone cannot be what a correction is judged by.
"""
import numpy as np

from pca2d.lblscan import amplitude_at, signal_stats


def series(n=400, drift=0.0, K=0.0, period=10.0, noise=1.0, seed=4,
           per_night=4):
    """Exposures grouped into nights, with a drift, a sinusoid and noise."""
    rng = np.random.default_rng(seed)
    nights = np.repeat(np.arange(n // per_night), per_night)
    t = 60000.0 + nights * 3.0 + rng.uniform(0.1, 0.4, size=nights.size)
    v = drift * (t - t.mean()) / np.ptp(t)
    v = v + K * np.sin(2 * np.pi * t / period) + rng.normal(0, noise, t.size)
    return t, v, np.full(t.size, noise)


def test_a_sinusoid_is_recovered_with_its_amplitude():
    t, v, e = series(K=20.0, period=17.3, noise=1.0)
    K, sK = amplitude_at(t, v, e, 17.3)
    assert abs(K - 20.0) < 0.5, K
    assert sK < 0.5
    off, _ = amplitude_at(t, v, e, 4.1)
    assert off < 2.0, "nothing at a period the series does not have"


def test_the_drift_is_seen_and_taken_out():
    t, v, e = series(drift=100.0, noise=1.0)
    s = signal_stats(t, v, e)
    assert 90.0 < s["drift"] < 110.0
    assert s["residual"] < 1.3, "once the drift is gone, only the noise"
    assert s["rms"] > 25.0, "the rms is dominated by the drift"
    assert abs(s["in_night"] - 1.0) < 0.4, "a night cannot feel a drift"


def test_eating_the_drift_flatters_the_rms_and_nothing_else():
    """The TOI-1452 case, in miniature: half the drift absorbed, the rms halves,
    and nothing that measures precision improves."""
    t, v, e = series(drift=100.0, noise=1.0)
    before = signal_stats(t, v, e)
    eaten = v - 0.5 * (v - np.median(v))          # 50% of everything removed
    after = signal_stats(t, eaten, e)
    assert after["rms"] < 0.6 * before["rms"], "the rms is delighted"
    assert after["drift"] < 0.6 * before["drift"], "because the drift went"
    assert after["in_night"] < before["in_night"], \
        "a proportional cut shrinks the noise too, so in_night alone is not proof"
    assert after["residual"] < before["residual"]
    # what IS proof: the drift fell, and that is astrophysics
    assert before["drift"] > 50.0 and after["drift"] < before["drift"] * 0.75


def test_a_known_period_is_reported_for_every_series():
    t, v, e = series(K=5.0, period=241.8, noise=1.0)
    s = signal_stats(t, v, e, periods=(241.8, 11.06))
    assert len(s["amplitudes"]) == 2
    assert abs(s["amplitudes"][0][0] - 5.0) < 1.0
    assert s["amplitudes"][1][0] < 1.5, "nothing at the other period"
    assert "rms" in s and "robust" in s, "it still carries the usual statistics"


def test_a_short_series_does_not_raise():
    t, v, e = series(n=4, per_night=1, noise=1.0)
    s = signal_stats(t, v, e, degree=2)
    assert np.isfinite(s["rms"])
    assert np.isnan(s["in_night"]), "no night has three exposures"
