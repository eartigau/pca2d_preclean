"""Where LBL writes, and how the spectra get there.

Until 2026-09-16 LBL's tree was `lbl` beside wherever the run started: the
window starts a run in its configuration's folder, so SMETHELLS_20's report was
on the data disk and its velocities in another clone's `lbl/`. And the window's
LBL page travelled under no flag at all, so a copy asked for there was never
made. The tree is now under the output root unless it is named, every LBL
setting has a flag, and the link or copy question is asked of the disk.
"""

import argparse
import os

import pytest

from pca2d import lbl as splbl
from pca2d.cli import SETTING_FLAGS, add_setting_flags, apply_setting_flags
from pca2d.config import lbl_directory


def flags(argv):
    parser = add_setting_flags(argparse.ArgumentParser())
    return parser.parse_args(argv)


def test_the_tree_is_under_the_output_root_unless_it_is_named(tmp_path):
    config = {"output": {"directory": str(tmp_path / "out" / "_a_name")},
              "lbl": {"directory": None}}
    # the root as it was before a run name added its level
    assert lbl_directory(config, str(tmp_path / "out")) == \
        str(tmp_path / "out" / "lbl")
    assert lbl_directory(config) == str(tmp_path / "out" / "_a_name" / "lbl")
    config["lbl"]["directory"] = str(tmp_path / "elsewhere")
    assert lbl_directory(config, str(tmp_path / "out")) == \
        str(tmp_path / "elsewhere")


def test_a_relative_tree_is_written_in_full(tmp_path, monkeypatch):
    """Every stage reads it back from the resolved config, and a relative path
    is another folder from another working directory."""
    monkeypatch.chdir(tmp_path)
    assert lbl_directory({"lbl": {"directory": "lbl"}}) == str(tmp_path / "lbl")


def test_every_lbl_setting_of_the_window_has_a_flag():
    from pca2d.gui import LBL_FLAGS, OPTIONS_LBL, OPTIONS_PATHS
    known = {flag: path for flag, path, _kind, _help in SETTING_FLAGS}
    for key, path, _kind in OPTIONS_LBL + OPTIONS_PATHS:
        assert key in LBL_FLAGS, "%s is shown and would travel nowhere" % key
        assert known.get(LBL_FLAGS[key]) == path, key
    assert known["--lbl-dir"] == "lbl.directory"


def test_the_flags_land_in_the_configuration():
    config = {"lbl": {"link": "symlink", "steps": ["template"], "run": True},
              "highpass": {"window": 201, "width_kms": None, "mode": "log_sub",
                           "polyorder": 2},
              "domain": {"dv": 0.5}}
    apply_setting_flags(config, flags(["--lbl-link", "copy", "--lbl-steps",
                                       "compute,compile", "--lbl-run", "false",
                                       "--lbl-dir", "/somewhere/lbl"]))
    assert config["lbl"] == {"link": "copy", "steps": ["compute", "compile"],
                             "run": False, "directory": "/somewhere/lbl"}
    with pytest.raises(SystemExit):
        flags(["--lbl-link", "hardlink"])


def test_copy_asked_is_a_copy_and_a_disk_that_links_gets_links(tmp_path):
    assert splbl.link_mode(str(tmp_path), "copy") == ("copy", None)
    assert splbl.link_mode(str(tmp_path), "symlink") == ("symlink", None)
    assert os.listdir(str(tmp_path / "science")) == [], "the probe cleans up"


def test_a_disk_that_refuses_a_link_gets_copies_and_says_why(tmp_path,
                                                             monkeypatch):
    def refuse(*_args, **_kw):
        raise OSError(45, "Operation not supported")

    monkeypatch.setattr(os, "symlink", refuse)
    mode, why = splbl.link_mode(str(tmp_path), "symlink")
    assert mode == "copy"
    assert "cannot hold a symbolic link" in why and "twice" in why
    assert os.listdir(str(tmp_path / "science")) == []


def test_copies_are_files_and_links_are_links(tmp_path):
    source = tmp_path / "data"
    source.mkdir()
    spectrum = source / "a_t.fits"
    spectrum.write_bytes(b"x" * 100)
    linked, kept, strangers = splbl.link_spectra([str(spectrum)],
                                                 str(tmp_path / "L"), "symlink")
    assert (linked, kept, strangers) == (1, 0, [])
    assert os.path.islink(str(tmp_path / "L" / "a_t.fits"))
    splbl.link_spectra([str(spectrum)], str(tmp_path / "C"), "copy")
    copied = tmp_path / "C" / "a_t.fits"
    assert not os.path.islink(str(copied)) and copied.read_bytes() == b"x" * 100


def test_the_window_always_names_the_tree_and_only_the_changed_settings():
    from pca2d.gui import build_command

    state = {"objects": ["TOI756"], "lbl_dir": "/disk/out/lbl",
             "lbl_link": "copy", "run": True, "lbl_steps": "template, mask"}
    config = {"lbl": {"link": "symlink", "run": True,
                      "steps": ["template", "mask"]}}
    argv = build_command(state, config)
    assert argv[argv.index("--lbl-dir") + 1] == "/disk/out/lbl"
    assert argv[argv.index("--lbl-link") + 1] == "copy"
    assert "--lbl-run" not in argv and "--lbl-steps" not in argv, \
        "what the configuration already says is not repeated"
    # nothing to compare with: every setting shown is said
    everything = build_command(state)
    assert "--lbl-run" in everything
    assert everything[everything.index("--lbl-steps") + 1] == "template,mask"


def test_the_proposed_tree_follows_the_output_root_and_a_typed_one_does_not():
    from pca2d.gui import App, lbl_proposal

    class Var:
        def __init__(self, value=""):
            self.value = value

        def get(self):
            return self.value

        def set(self, value):
            self.value = value

    w = App.__new__(App)
    w.vars = {"out_dir": Var("/disk/a"), "lbl_dir": Var("")}
    w._proposed_lbl = ""
    w._follow_lbl_dir()
    assert w.vars["lbl_dir"].get() == lbl_proposal("/disk/a") == "/disk/a/lbl"
    w.vars["out_dir"].set("/disk/b")
    w._follow_lbl_dir()
    assert w.vars["lbl_dir"].get() == "/disk/b/lbl"
    # an existing tree, typed: left alone from then on
    w.vars["lbl_dir"].set("/Volumes/irrisor/pca2d_preclean/lbl")
    w.vars["out_dir"].set("/disk/c")
    w._follow_lbl_dir()
    assert w.vars["lbl_dir"].get() == "/Volumes/irrisor/pca2d_preclean/lbl"
    assert lbl_proposal("") == "", "no root, no proposal"
