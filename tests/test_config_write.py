"""Saving a setting as the default must not cost the comments.

Every value in config.yaml is followed by the measurement that chose it, which
is most of what the file is worth, and a yaml.safe_dump of a parsed document
deletes all of it. So the writer edits lines.
"""
import os
import shutil

import pytest
import yaml

from pca2d import config

#: the repository's own config.yaml, found from this file rather than from the
#: working directory: the suite is also run from the runtime environment's bin
CONFIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "config.yaml")


def test_a_value_changes_and_its_comment_stays(tmp_path):
    path = str(tmp_path / "config.yaml")
    shutil.copy(CONFIG, path)
    before = open(path).read()
    written = config.update_file(path, {"correct.mask": "exposure",
                                        "twoframe.n_star": 0,
                                        "lbl.run": False})
    assert sorted(written) == ["correct.mask", "lbl.run", "twoframe.n_star"]
    body = yaml.safe_load(open(path))["general"]
    assert body["correct"]["mask"] == "exposure"
    assert body["twoframe"]["n_star"] == 0
    assert body["lbl"]["run"] is False
    after = open(path).read()
    assert "which samples a corrected file blanks" in after, \
        "the comment that explains the value it replaced is still there"
    assert after.count("#") >= before.count("#") - 1
    assert len(after.split("\n")) == len(before.split("\n"))


def test_a_window_is_written_as_text_so_a_reader_gets_it_back(tmp_path):
    """1220:4 unquoted is the base-60 number 73204 (yaml-windows-as-text)."""
    path = str(tmp_path / "config.yaml")
    shutil.copy(CONFIG, path)
    config.update_file(path, {"highpass.width_kms": 150.0,
                              "input.pattern": "*t.fits"})
    body = yaml.safe_load(open(path))["general"]
    assert body["highpass"]["width_kms"] == 150.0
    assert body["input"]["pattern"] == "*t.fits"
    assert config.as_yaml_value("1220:4") == '"1220:4"'
    assert config.as_yaml_value("auto") == "auto"
    assert config.as_yaml_value(True) == "true"
    assert config.as_yaml_value(None) == "null"
    assert config.as_yaml_value([1, 2]) == "[1, 2]"


def test_a_key_that_is_not_there_is_added_under_its_section(tmp_path):
    path = str(tmp_path / "config.yaml")
    shutil.copy(CONFIG, path)
    config.update_file(path, {"twoframe.star_smooth": 7})
    assert yaml.safe_load(open(path))["general"]["twoframe"]["star_smooth"] == 7


def test_a_section_that_is_not_there_is_refused_rather_than_guessed(tmp_path):
    path = str(tmp_path / "config.yaml")
    shutil.copy(CONFIG, path)
    with pytest.raises(SystemExit):
        config.update_file(path, {"nonesuch.key": 1})
    assert open(path).read() == open(CONFIG).read(), "nothing written"
