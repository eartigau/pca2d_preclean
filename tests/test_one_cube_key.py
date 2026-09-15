"""One run, one cube key.

Every stage is handed the resolved configuration the run wrote and loads it
again; the plan, meanwhile, hashed the configuration it had in memory. If those
two disagree about one hashed value, the run has two cube keys, and on
2026-09-15 it had exactly that:

    260915 16:39:17.44 | cube cache key: cache/cube_tfits_3e01c018993b
    260915 16:39:18.29 | stage fit
    there is no cube at cache/cube_tfits_4e793f8df25f

The cube stage loaded a cube; the fit opened a path nothing had ever written.
--dv from the window changed the grid step after load_config had already
derived highpass.window from highpass.width_kms at the file's step, and the
window is hashed while the width is not.

So: a configuration a run writes must be a fixed point of reading it back.
"""

import argparse

import yaml

from pca2d.cli import SETTING_FLAGS, apply_setting_flags
from pca2d.config import cache_key, load_config

CONFIG = """\
general:
  input:
    directory: data
    format: tfits
  domain:
    wave_min: 1000.0
    wave_max: 2400.0
    dv: 0.5
    smart_dv: false
  highpass:
    mode: log_sub
    width_kms: 100.0
    window: 201
    polyorder: 2
"""


def flags(**given):
    """An args namespace with every setting flag unset but the ones named."""
    args = argparse.Namespace(**{f[2:].replace("-", "_"): None
                                 for f, _p, _k, _h in SETTING_FLAGS})
    for name, value in given.items():
        setattr(args, name, value)
    return args


def written(config, tmp_path, name="resolved_config.yaml"):
    """The config as a run writes it, read back as a stage reads it."""
    path = tmp_path / name
    with open(path, "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return load_config(str(path))


def test_the_grid_step_from_the_window_does_not_make_a_second_key(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(CONFIG)
    config = load_config(str(path))
    assert config["highpass"]["window"] == 201        # 100 km/s at dv = 0.5

    apply_setting_flags(config, flags(dv=1.0))
    assert config["highpass"]["window"] == 101        # 100 km/s at dv = 1.0
    assert cache_key(config) == cache_key(written(config, tmp_path))


def test_the_high_pass_from_the_window_does_not_either(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(CONFIG)
    config = load_config(str(path))

    apply_setting_flags(config, flags(high_pass=50.0))
    assert config["highpass"]["window"] == 101        # 50 km/s at dv = 0.5
    assert cache_key(config) == cache_key(written(config, tmp_path))


def test_a_configuration_nothing_changed_is_already_a_fixed_point(tmp_path):
    """The guard above must not be the only thing holding it: a plain run,
    with no flag at all, hashes the same cube its stages will open."""
    path = tmp_path / "config.yaml"
    path.write_text(CONFIG)
    config = load_config(str(path))
    assert cache_key(config) == cache_key(written(config, tmp_path))


def test_a_flag_that_derives_nothing_leaves_the_key_alone(tmp_path):
    """--iters is downstream of the cube, so it must not re-key it."""
    path = tmp_path / "config.yaml"
    path.write_text(CONFIG)
    config = load_config(str(path))
    before = cache_key(config)
    apply_setting_flags(config, flags(iters=3))
    assert cache_key(config) == before == cache_key(written(config, tmp_path))
