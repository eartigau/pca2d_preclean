"""A column the exposures agree on goes whole, one bad exposure does not.

The residual divided by its own local sigma is a z, so at a fixed OBSERVER
column the exposures can be compared with each other: a column too many of
them are clipped in, whose survivors still carry more variance than noise,
belongs to no exposure's star and is dropped from all of them (asked for on
2026-09-25). Both conditions, because either alone throws away good data: a
cosmic ray shower in one exposure meets neither, and a column where a tenth of
the exposures is genuinely bad but the rest is clean meets only the first.
"""

import numpy as np
import pytest

from pca2d.outliers import expected_chi2, running_stats


def _columns(resid, nsig, window, frac_max, chi2_max, min_rows=10):
    """The column verdict of outliers.residual_outliers, on one parity's rows.

    The same arithmetic, without a cube: the module's own loop reads one from
    disk, and what is under test here is the criterion.
    """
    med, sig = running_stats(resid, window)
    with np.errstate(invalid="ignore"):
        z = (resid - med) / sig
    measured = np.isfinite(z)
    cut = measured & (np.abs(z) > nsig)
    survivor = np.where(measured & ~cut, z, np.nan)
    seen = measured.sum(axis=0)
    frac = np.divide(cut.sum(axis=0), seen, out=np.zeros(resid.shape[1]),
                     where=seen > 0)
    with np.errstate(invalid="ignore"):
        chi2 = np.nanmean(survivor ** 2, axis=0)
    return (seen >= min_rows) & (frac > frac_max) & (chi2 > chi2_max), frac, chi2


def test_the_survivors_of_a_clip_carry_less_than_one():
    """The threshold is measured against this, not against 1: the clip takes
    the tails with it."""
    assert expected_chi2(3.0) == pytest.approx(0.9733, abs=5e-4)
    assert expected_chi2(2.5) == pytest.approx(0.9113, abs=5e-4)
    assert expected_chi2(10.0) == pytest.approx(1.0, abs=1e-6)
    # and the default sits well above the noise it has to ignore
    from pca2d.config import DEFAULTS
    assert DEFAULTS["correct"]["column_chi2"] > expected_chi2(3.0) + 0.4


def test_noise_alone_loses_no_column():
    rng = np.random.default_rng(4)
    resid = rng.normal(size=(120, 900))
    bad, frac, chi2 = _columns(resid, 3.0, 201, 0.10, 1.5)
    assert not bad.any(), "%d columns of pure noise dropped" % bad.sum()
    assert frac.max() < 0.10 and chi2.max() < 1.5


def test_a_column_with_twice_the_noise_goes_whole():
    rng = np.random.default_rng(5)
    resid = rng.normal(size=(120, 900))
    resid[:, 400] *= 2.5
    bad, frac, chi2 = _columns(resid, 3.0, 201, 0.10, 1.5)
    assert bad[400], "frac %.3f, chi2 %.2f" % (frac[400], chi2[400])
    assert bad.sum() == 1, "and nothing else"


def test_one_bad_exposure_does_not_cost_the_column():
    """A cosmic ray in a single exposure is what the per-sample clip is for."""
    rng = np.random.default_rng(6)
    resid = rng.normal(size=(120, 900))
    resid[7, 400] = 40.0
    bad, frac, chi2 = _columns(resid, 3.0, 201, 0.10, 1.5)
    assert not bad.any()
    assert frac[400] == pytest.approx(1 / 120, abs=1e-6), "it was still clipped"


def test_a_tenth_of_the_exposures_bad_but_clean_survivors_keeps_the_column():
    """The chi2 half of the rule: 15 exposures of 120 ruined at one column, the
    other 105 as good as their neighbours. The column stays, those 15 samples
    go."""
    rng = np.random.default_rng(7)
    resid = rng.normal(size=(120, 900))
    resid[:15, 400] = 30.0
    bad, frac, chi2 = _columns(resid, 3.0, 201, 0.10, 1.5)
    assert frac[400] > 0.10, "more than a tenth was clipped"
    assert chi2[400] < 1.5, "but what is left is noise: %.2f" % chi2[400]
    assert not bad[400]


def test_both_thresholds_reach_the_settings_of_a_run():
    """The correct stage is a separate process: a setting that travels under no
    flag is a setting the run ignores."""
    import argparse

    from pca2d.cli import SETTING_FLAGS, add_setting_flags, apply_setting_flags, clip_args
    from pca2d.config import load_config

    known = {path: flag for flag, path, _k, _h in SETTING_FLAGS}
    for path in ("correct.nsig_cut", "correct.column_frac", "correct.column_chi2"):
        assert path in known, path
    args = add_setting_flags(argparse.ArgumentParser()).parse_args(
        ["--nsig-cut", "3", "--column-frac", "0.1", "--column-chi2", "1.5"])
    cfg = apply_setting_flags(load_config("config.yaml", instrument="SPIROU"), args)
    assert cfg["correct"]["nsig_cut"] == 3.0
    passed = clip_args(cfg)
    assert passed[:2] == ["--nsig-cut", "3.0"]
    assert "--column-frac" in passed and "--column-chi2" in passed
    # and nothing about the columns when the clip itself is off
    cfg["correct"]["nsig_cut"] = None
    assert clip_args(cfg) == []
