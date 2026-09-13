"""Keeping the nights a list names, which is not a range.

A joint fit assumes the atmosphere belongs to the night, so it only means
something for stars observed on the SAME nights, and the nights two campaigns
share are scattered rather than contiguous: GJ 1 and GJ 3090 share 47 of 147 and
99, over two seasons.
"""
import numpy as np

from pca2d.config import DEFAULTS, cache_key, load_config
from pca2d.cube import _quality_reason


def meta(rjd, snr=100.0):
    return {"bjd": 2400000.0 + rjd, "snr_band": snr, "median_flux": 1.0,
            "nan_fraction": 0.0, "berv": 10.0}


WAVE = np.linspace(1500.0, 1600.0, 64)


def rejected(rjd, quality):
    """Why this exposure would be dropped, or None."""
    return _quality_reason(meta(rjd), WAVE, quality)


QUALITY = dict(DEFAULTS["quality"])


def test_a_night_in_the_list_is_kept_and_one_outside_is_not():
    quality = dict(QUALITY, nights=[60053, 60060, 60162])
    assert rejected(60053.31, quality) is None
    assert rejected(60060.98, quality) is None, "late in the night"
    assert rejected(60162.04, quality) is None, "after midnight UT"
    said = rejected(60061.5, quality)
    assert said and "not one of the 3 nights" in said


def test_the_list_is_not_a_range():
    """What a min/max window cannot express, which is the whole point."""
    quality = dict(QUALITY, nights=[60000, 60500])
    assert rejected(60000.2, quality) is None
    assert rejected(60500.2, quality) is None
    assert rejected(60250.2, quality) is not None, \
        "between the two, and not asked for"


def test_no_list_keeps_every_night():
    for value in (None, [], False):
        quality = dict(QUALITY, nights=value)
        assert rejected(59000.5, quality) is None
        assert rejected(62000.5, quality) is None


def test_the_nights_are_part_of_the_cube_key_but_only_when_asked_for(tmp_path):
    """They change which exposures the cube holds, so a cube built with them is
    a different cube; unset, they must not orphan the cubes already built."""
    base = load_config(None)
    plain = cache_key(base)
    base["quality"]["nights"] = None
    assert cache_key(base) == plain, "unset changes no key"
    base["quality"]["nights"] = [60053, 60060]
    keyed = cache_key(base)
    assert keyed != plain
    base["quality"]["nights"] = [60053, 60061]
    assert cache_key(base) != keyed, "a different set of nights, a different cube"


def test_it_survives_a_variant_file(tmp_path):
    """This is how the experiment is written down: a variant naming the nights."""
    import yaml
    path = tmp_path / "shared.yaml"
    nights = [60145 + 3 * i for i in range(8)]
    path.write_text(yaml.safe_dump({"quality": {"nights": nights}}))
    variant = yaml.safe_load(path.read_text())
    config = load_config(None, variant=variant)
    assert config["quality"]["nights"] == nights
    quality = config["quality"]
    assert rejected(nights[0] + 0.4, quality) is None
    assert rejected(nights[0] + 1.4, quality) is not None
