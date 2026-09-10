"""The cache is one labelled folder, and emptying it touches only what is ours."""

import os

from pca2d import cache


def config_for(path):
    return {"output": {"cache_directory": str(path)}}


def populate(root):
    (root / "cube_tfits_abc123").mkdir(parents=True)
    (root / "cube_tfits_abc123" / "data.npy").write_bytes(b"x" * 1000)
    (root / "cube_tfits_abc123" / "snippets").mkdir()
    (root / "cube_tfits_abc123" / "snippets" / "w.npz").write_bytes(b"y" * 10)
    (root / "telluric_abc123.npz").write_bytes(b"z" * 100)
    (root / "someone_elses_notes.txt").write_text("keep me")


def test_the_folder_says_what_it_is(tmp_path):
    root = tmp_path / "cache"
    cache.ensure(config_for(root))
    text = (root / "README.txt").read_text()
    assert "pca2d-preclean" in text and "--clean-cache" in text


def test_a_dry_run_removes_nothing(tmp_path):
    root = tmp_path / "cache"
    populate(root)
    count, nbytes = cache.clean(config_for(root), dry_run=True)
    assert count == 2 and nbytes == 1110
    assert (root / "cube_tfits_abc123" / "data.npy").exists()
    assert (root / "telluric_abc123.npz").exists()


def test_cleaning_removes_only_what_this_package_wrote(tmp_path):
    root = tmp_path / "cache"
    populate(root)
    cache.ensure(config_for(root))
    count, _ = cache.clean(config_for(root))
    assert count == 2
    assert not (root / "cube_tfits_abc123").exists()
    assert not (root / "telluric_abc123.npz").exists()
    assert (root / "someone_elses_notes.txt").read_text() == "keep me"
    assert (root / "README.txt").exists(), "the label goes with the folder"


def test_a_link_in_the_cache_is_removed_as_a_link(tmp_path):
    root = tmp_path / "cache"
    root.mkdir()
    precious = tmp_path / "precious"
    precious.mkdir()
    (precious / "result.fits").write_bytes(b"r")
    os.symlink(precious, root / "cube_tfits_linked")
    cache.clean(config_for(root))
    assert (precious / "result.fits").exists(), "a symlink was followed"


def test_an_absent_cache_is_not_an_error(tmp_path):
    assert cache.clean(config_for(tmp_path / "nowhere")) == (0, 0)


def test_the_command_line_offers_it():
    from pca2d.cli import parse_args

    assert parse_args(["--object", "X", "--clean-cache"]).clean_cache is True
    assert parse_args(["--object", "X"]).clean_cache is False
