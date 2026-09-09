"""Where a run reads and writes, and how the grid step is chosen.

Both of these have already gone wrong in ways nothing downstream notices. The
object name was appended to the input root twice, so every default run looked
for `data/TOI-2120/TOI-2120`; and `domain.dv` is hashed into the cube cache
key, so a step that wobbles by a hair between two reductions silently orphans a
cube that took twenty minutes to build.

No spectrum is read here: `load_config(None, ...)` needs no file at all, and
the arithmetic of the smart step is a pure function of a measured number.
"""

import os

import pytest
import yaml

from pca2d.config import (SMART_DV_FRACTION, cache_key, load_config,
                          smart_dv_from_step, spectra_dir)


# --------------------------------------------------------------- the roots ---
def test_the_object_is_joined_to_the_input_root_exactly_once():
    cfg = load_config(None, object_name="TOI-2120")
    assert cfg["input"]["directory"] == "data", "the root must stay a root"
    assert spectra_dir(cfg) == os.path.join("data", "TOI-2120")


def test_both_roots_can_be_overridden_for_one_run():
    cfg = load_config(None, object_name="TOI-2120", data_dir="/mnt/spirou",
                      out_dir="/scratch/run2")
    assert spectra_dir(cfg) == os.path.join("/mnt/spirou", "TOI-2120")
    assert cfg["output"]["directory"] == "/scratch/run2"


def test_a_saved_config_resolves_to_the_same_folder(tmp_path):
    """A run writes its resolved config out; loading it back must not drift."""
    cfg = load_config(None, object_name="TOI-2120")
    path = tmp_path / "resolved_config.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    again = load_config(str(path), object_name="TOI-2120")
    assert spectra_dir(again) == spectra_dir(cfg)


def test_a_config_written_before_the_split_still_points_at_its_spectra():
    """`directory` used to hold the object's folder itself, not the root."""
    cfg = {"input": {"directory": "data/TOI-2120", "object": "TOI-2120"}}
    assert spectra_dir(cfg) == "data/TOI-2120"


def test_an_object_block_may_move_the_root(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({
        "general": {"input": {"directory": "data"}},
        "objects": {"TOI-2120": {"input": {"directory": "/archive"}}},
    }))
    cfg = load_config(str(path), object_name="TOI-2120")
    assert spectra_dir(cfg) == os.path.join("/archive", "TOI-2120")
    # ... and the command line still outranks it
    cfg = load_config(str(path), object_name="TOI-2120", data_dir="/mnt")
    assert spectra_dir(cfg) == os.path.join("/mnt", "TOI-2120")


def test_nothing_a_run_writes_lands_under_the_input_root(tmp_path):
    """The corrected spectra used to be written into the input directory."""
    from pca2d.cli import parse_args, resolve

    data, out = tmp_path / "data", tmp_path / "outputs"
    (data / "TOI-2120").mkdir(parents=True)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({"general": {"input": {"pattern": "*t.fits"}}}))

    plan = resolve(parse_args(["--object", "TOI-2120", "--config", str(path),
                               "--data-dir", str(data), "--out-dir", str(out)]))
    assert plan["directory"] == str(data / "TOI-2120")
    assert plan["outdir"].startswith(str(out))
    assert plan["corrdir"].startswith(plan["outdir"])
    assert not plan["corrdir"].startswith(str(data))


# ------------------------------------------------------------- the step ------
def test_the_smart_step_is_a_fraction_of_the_finest_pixel():
    # SPIRou samples at about 2.27 km/s, so 70% of it is 1.589, floored to 1.58
    assert smart_dv_from_step(2.27) == pytest.approx(1.58)
    assert SMART_DV_FRACTION == 0.7


def test_the_smart_step_rounds_down_and_never_up():
    for step in (2.27, 2.2701, 2.2749, 3.0, 1.4285):
        dv = smart_dv_from_step(step)
        assert dv <= SMART_DV_FRACTION * step, "the grid must stay finer"
        assert step / dv > 1.0, "the grid must not be coarser than a pixel"
        assert round(dv * 100) == dv * 100, "not rounded to 10 m/s"


def test_a_wavelength_solution_that_moved_by_a_hair_keeps_its_cube():
    """The whole point of rounding: 1 m/s of drift must not re-key the cube."""
    base = load_config(None, object_name="TOI-2120")
    keys = set()
    for step in (2.2700, 2.2703, 2.2699):
        cfg = dict(base, domain=dict(base["domain"],
                                     dv=smart_dv_from_step(step)))
        keys.add(cache_key(cfg))
    assert len(keys) == 1


def test_the_fraction_may_be_given_and_must_be_a_fraction():
    assert smart_dv_from_step(2.27, 0.5) == pytest.approx(1.13)
    for bad in (0, 1, 1.5, -0.3):
        with pytest.raises(ValueError):
            smart_dv_from_step(2.27, bad)


def test_what_the_config_asks_for_is_not_used_when_smart_dv_is_on(tmp_path,
                                                                  monkeypatch):
    """The one thing a reader of the config has to be sure of.

    Whatever `domain.dv` says, smart_dv overwrites it: the same run, the same
    step and the same cube, whether the file asks for 0.5, for 4.0, or for
    nothing at all.
    """
    from pca2d import config as _config

    monkeypatch.setattr(_config, "measure_pixel_dv", lambda *a, **k: 2.27)
    (tmp_path / "data" / "TOI-2120").mkdir(parents=True)

    steps, keys = set(), set()
    for asked in (0.5, 4.0, None):
        path = tmp_path / ("config_%s.yaml" % asked)
        path.write_text(yaml.safe_dump(
            {"general": {"domain": {"smart_dv": True, "dv": asked}}}))
        cfg = load_config(str(path), object_name="TOI-2120",
                          data_dir=str(tmp_path / "data"))
        steps.add(cfg["domain"]["dv"])
        keys.add(cache_key(cfg))
    assert steps == {1.58}, "domain.dv leaked into the resolved step"
    assert len(keys) == 1, "a dead domain.dv re-keyed the cube"


def test_the_step_is_hashed_into_the_cube_cache_key():
    base = load_config(None, object_name="TOI-2120")
    other = dict(base, domain=dict(base["domain"], dv=1.58))
    assert cache_key(other) != cache_key(base)
