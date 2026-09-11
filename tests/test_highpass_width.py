"""The high pass is a width in km/s; its window in samples follows the grid.

And the velocity term is off unless asked for (2026-09-11).
"""
import yaml

from pca2d.config import DEFAULTS, cache_key, highpass_samples, load_config


def test_a_hundred_km_s_is_201_samples_at_half_a_km_s():
    assert highpass_samples(100, 0.5) == 201


def test_the_same_width_on_a_coarser_grid_is_fewer_samples():
    assert highpass_samples(100, 1.37) == 73
    assert highpass_samples(100, 2.0) == 51          # 50, made odd


def test_the_window_is_never_shorter_than_the_polynomial_allows():
    assert highpass_samples(1, 0.5, polyorder=2) == 5


def test_the_defaults_filter_over_100_km_s():
    cfg = load_config(None)
    assert cfg["highpass"]["width_kms"] == 100.0
    assert cfg["highpass"]["window"] == 201


def _write(tmp_path, block):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({"general": block}))
    return str(path)


def test_a_width_in_the_file_sets_the_window(tmp_path):
    cfg = load_config(_write(tmp_path, {"highpass": {"width_kms": 50}}))
    assert cfg["highpass"]["window"] == 101


def test_a_config_in_samples_keeps_its_window_and_its_cube(tmp_path):
    """A run saved before 2026-09-11 says `window: 151` and nothing in km/s."""
    cfg = load_config(_write(tmp_path, {"highpass": {"window": 151}}))
    assert cfg["highpass"]["window"] == 151
    assert cfg["highpass"]["width_kms"] is None
    # width_kms is not in the key, so the key is the one that cube always had
    without = dict(cfg, highpass={k: v for k, v in cfg["highpass"].items()
                                  if k != "width_kms"})
    assert cache_key(cfg) == cache_key(without)


def test_resolving_one_config_leaves_the_defaults_alone(tmp_path):
    load_config(_write(tmp_path, {"highpass": {"window": 151}}))
    assert DEFAULTS["highpass"]["width_kms"] == 100.0
    assert load_config(None)["highpass"]["window"] == 201


def test_the_velocity_term_is_off_unless_asked_for():
    assert DEFAULTS["twoframe"]["velocity_term"] is False
    assert load_config(None)["twoframe"]["velocity_term"] is False
