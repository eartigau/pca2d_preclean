"""The fit finds the spectra a cube was built from, so it can read their
headers (the water column of the correlation matrix)."""
import os

from pca2d.plotting import source_directory


def _cube(tmp_path, text):
    cube = tmp_path / "cube"
    cube.mkdir()
    (cube / "cube_config.yaml").write_text(text)
    return str(cube)


def test_the_object_s_folder_under_the_root(tmp_path):
    """input.directory is the root; the spectra are in root/object. On
    2026-09-11 the root was returned, no header was found and the water row
    had silently gone from every correlation matrix."""
    cube = _cube(tmp_path, "input:\n  directory: data\n  object: GJ1\n")
    assert source_directory(cube) == os.path.join("data", "GJ1")


def test_a_config_from_before_the_split_keeps_its_folder(tmp_path):
    cube = _cube(tmp_path, "input:\n  directory: data/GJ1\n  object: GJ1\n")
    assert source_directory(cube) == "data/GJ1"


def test_no_config_no_folder(tmp_path):
    assert source_directory(str(tmp_path / "nothing")) is None
    assert source_directory(_cube(tmp_path, "output: {}\n")) is None
