"""A fresh clone has no `data/`, and `./check.sh` is the first thing anybody
runs on one. So nothing in the suite may need a spectrum to exist.

This file pins the two halves of that: a configuration resolves for a named
object whose folder is nowhere, and a run still says plainly that the folder is
missing rather than going on without it.
"""

from __future__ import annotations

import types

import pytest

from pca2d.config import cache_key, load_config


def test_a_configuration_resolves_with_no_spectra_anywhere(tmp_path):
    """The check is the whole point of the file: this used to raise SystemExit
    and took nine tests of the suite with it on any clone without data."""
    config = load_config("config.yaml", object_name="PROXIMA",
                         data_dir=str(tmp_path / "no-such-root"))
    assert config["input"]["object"] == "PROXIMA"
    # the object's own block still merged, which is what these tests read
    assert config["target"]["planets"], "the object block is merged regardless"
    # and the key is computable, so a test may compare two configurations
    assert cache_key(config)


def test_a_run_still_refuses_a_missing_folder(tmp_path):
    """Tolerating it in load_config must not make a RUN quietly proceed: the
    clear message moved to the two places that are about to read spectra."""
    from pca2d.cli import joint_members, resolve

    args = types.SimpleNamespace(
        object="PROXIMA", objects=["PROXIMA", "GJ1"], config="config.yaml",
        data_dir=str(tmp_path / "no-such-root"), out_dir=str(tmp_path / "out"),
        instrument="NIRPS", n_star=0, n_earth=None, windows=None,
        rebuild_cube=False, run_lbl=False, variant=None, name=None,
        min_rjd=None, max_rjd=None)
    with pytest.raises(SystemExit):
        resolve(args)
    with pytest.raises(SystemExit):
        joint_members(args, None)
