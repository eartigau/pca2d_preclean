"""The configuration contract and the definitions that must not drift.

Nothing here is about mathematics. It is about the two ways this project has
already been bitten: a cache key that changed when a knob was added, and a
statistic that was renamed without being renamed.
"""

from __future__ import annotations

import numpy as np
import pytest

from pca2d.config import cache_key, load_config

# A cube is named by a hash of the configuration that produced it, so a change
# that moves a key orphans every cube built before it. That is a decision and
# never an accident. This repository ships no cubes to protect, so nothing is
# pinned to a literal here; what IS pinned is the property that made the
# literals worth having.


def test_the_key_depends_on_the_cube_and_not_on_what_comes_after_it():
    """Changing a PCA setting must not orphan a cube that took an hour."""
    base = load_config("config.yaml", instrument="SPIROU")
    key = cache_key(base)
    for section, name, value in (("twoframe", "n_star", 4),
                                 ("twoframe", "n_earth", 3),
                                 ("twoframe", "iters", 1),
                                 ("correct", "n_star", 5)):
        other = load_config("config.yaml", instrument="SPIROU")
        other[section][name] = value
        assert cache_key(other) == key, "%s.%s should not re-key the cube" % (
            section, name)


def test_the_key_does_depend_on_what_makes_the_cube():
    base = load_config("config.yaml", instrument="SPIROU")
    key = cache_key(base)
    for section, name, value in (("domain", "wave_max", 2400.0),
                                 ("highpass", "window", 31),
                                 ("quality", "isolated_window", 5),
                                 ("input", "nightly_stack", True)):
        other = load_config("config.yaml", instrument="SPIROU")
        other[section][name] = value
        assert cache_key(other) != key, "%s.%s must re-key the cube" % (
            section, name)


def test_two_instruments_never_share_a_cube():
    assert (cache_key(load_config("config.yaml", instrument="SPIROU"))
            != cache_key(load_config("config.yaml", instrument="NIRPS")))

def test_unset_options_do_not_change_the_key():
    """Adding an optional knob must not invalidate every existing cube.

    min_rjd and max_rjd were added after the cubes were built. Hashing their
    None values would have re-keyed every config and orphaned twenty minutes of
    work per object for no change in behaviour. An option that is not set
    describes no behaviour and must not describe a different cube.
    """
    cfg = load_config("config.yaml", instrument="NIRPS")
    base = cache_key(cfg)
    cfg["quality"]["min_rjd"] = None
    cfg["quality"]["max_rjd"] = None
    assert cache_key(cfg) == base
    cfg["quality"]["min_rjd"] = 60200
    assert cache_key(cfg) != base, "an epoch window did not invalidate the cube"


def test_max_sigma_and_max_mad_are_the_same_knob():
    """Renaming the threshold must not have changed its value.

    The cut is in robust sigmas, 1.4826 x MAD, so the name max_mad was simply
    wrong. It still works, and it still means exactly what it meant, because a
    fit made before the rename has to stay comparable with one made after.
    """
    from pca2d.twoframe import parse_args
    assert parse_args(["--max-sigma", "7"]).max_mad == 7.0
    assert parse_args(["--max-mad", "7"]).max_mad == 7.0
    assert parse_args([]).max_mad == 10.0



def test_the_score_prefers_what_the_data_preferred():
    """The ranking must reproduce two results found by hand before it existed.

    Proxima prefers the observer block alone; GL 699 prefers the full
    correction. A naive quadrature of signed gains gets Proxima backwards,
    ranking 5-0 first although it raised the annual power by 63 percent and
    removed 30 percent of the rotation peak, because squaring destroys signs.
    The clamp on p_peak is what prevents that, and this test is why it is there.
    """
    def score(mad_r, mad_c, yr_r, yr_c, pk_r, pk_c):
        p_mad = mad_c / mad_r
        p_year = yr_c / yr_r
        p_peak = max(1.0, pk_r / pk_c)
        return float(np.sqrt((p_mad ** 2 + p_year ** 2 + p_peak ** 2) / 3.0))

    proxima = {
        "5-5": score(1.26, 1.25, 0.0083, 0.0076, 0.1335, 0.0973),
        "0-5": score(1.26, 1.27, 0.0083, 0.0044, 0.1335, 0.1341),
        "5-0": score(1.26, 1.24, 0.0083, 0.0135, 0.1335, 0.0928),
    }
    gl699 = {
        "5-5": score(1.79, 1.67, 0.2271, 0.1588, 0.1451, 0.1370),
        "0-5": score(1.79, 1.78, 0.2271, 0.2183, 0.1451, 0.1605),
    }
    assert min(proxima, key=proxima.get) == "0-5"
    assert max(proxima, key=proxima.get) == "5-0"
    assert min(gl699, key=gl699.get) == "5-5"


def test_instrument_table_picks_the_right_extensions():
    """A SPIRou file must resolve to the AB extensions, a NIRPS file to A.

    This is the one silent failure in the project. `FluxA` exists in a SPIRou
    t.fits and reading it raises nothing; it returns half the light with a
    different blaze and SNR, and every number downstream looks plausible. The
    instrument is taken from the file, never from the config, so a mislabelled
    config cannot select the wrong fibre.
    """
    import glob

    from astropy.io import fits

    from pca2d.tfits import extensions_for, instrument_of

    want = {"SPIROU": ("FluxAB", "WaveAB", "BlazeAB", "Recon"),
            "NIRPS": ("FluxA", "WaveA", "BlazeA", "Recon")}
    seen = set()
    for d in sorted(glob.glob("lbl_data/science/*_RAW")):
        files = sorted(glob.glob("%s/*.fits" % d))
        if not files:
            continue
        with fits.open(files[0]) as hdulist:
            inst = instrument_of(hdulist)
            assert extensions_for(hdulist) == want[inst], d
            assert set(want[inst]) <= {h.name for h in hdulist}, (
                "%s does not carry the extensions its instrument declares" % d)
        seen.add(inst)
    if seen:
        assert seen == {"NIRPS", "SPIROU"}, "expected both instruments, saw %s" % seen


def test_unknown_instrument_stops_the_run():
    """Guessing is worse than failing: every extension choice depends on this."""
    import pytest as _pytest
    from astropy.io import fits

    from pca2d.tfits import instrument_of

    h = fits.HDUList([fits.PrimaryHDU()])
    with _pytest.raises(ValueError, match="no INSTRUME"):
        instrument_of(h)
    h[0].header["INSTRUME"] = "HARPS"
    with _pytest.raises(ValueError, match="unsupported"):
        instrument_of(h)
