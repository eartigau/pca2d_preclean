"""A run's folder is a link to the external disk, made before the run writes.

Everything the run writes then lands out there directly. A folder that already
exists for real is moved first, and the links it holds that already point to
their own place out there (what the first move of the FITS left) are taken as
already moved rather than refused. If the disk is not mounted the run stops
rather than fill the internal disk again.
"""

import os

import numpy as np
import pytest
from astropy.io import fits

from pca2d import storage


@pytest.fixture
def tree(tmp_path, monkeypatch):
    repo, disk = tmp_path / "repo", tmp_path / "disk"
    repo.mkdir()
    disk.mkdir()
    monkeypatch.chdir(repo)
    return repo, disk, {"output": {"fits_directory": str(disk / "pca2d")}}


def test_a_new_run_folder_is_one_link_and_its_files_land_outside(tree):
    repo, disk, config = tree
    storage.check(config)
    storage.link_dir(config, "outputs/X/1-3v")
    os.makedirs("outputs/X/1-3v/corrected")
    (repo / "outputs/X/1-3v/twoframe_components.fits").write_bytes(b"basis")
    (repo / "outputs/X/1-3v/corrected/a.fits").write_bytes(b"spectrum")
    assert os.path.islink(repo / "outputs/X/1-3v")
    assert not os.path.islink(repo / "outputs/X")
    out = disk / "pca2d/outputs/X/1-3v"
    assert (out / "twoframe_components.fits").read_bytes() == b"basis"
    assert (out / "corrected/a.fits").read_bytes() == b"spectrum"


def test_astropy_overwriting_inside_a_linked_folder_stays_outside(tree):
    # why whole folders: overwrite removes the file and writes a new one,
    # which, inside a linked folder, is a new one on the external disk
    repo, disk, config = tree
    storage.check(config)
    storage.link_dir(config, "outputs/X/2-7v")
    path = "outputs/X/2-7v/twoframe_components.fits"
    fits.PrimaryHDU(np.zeros(3)).writeto(path)
    fits.PrimaryHDU(np.ones(5)).writeto(path, overwrite=True)
    assert fits.getdata(disk / "pca2d" / path).shape == (5,)


def test_linking_twice_changes_nothing(tree):
    _, _, config = tree
    storage.check(config)
    first = storage.link_dir(config, "lbl/lblrv")
    assert storage.link_dir(config, "lbl/lblrv") == first


def test_an_old_run_folder_is_moved_and_its_links_out_there_recognised(tree):
    repo, disk, config = tree
    storage.check(config)
    # the state the first move left: FITS out there, links to them here
    out = disk / "pca2d/outputs/X/2-7"
    (out / "corrected").mkdir(parents=True)
    (out / "corrected/a.fits").write_bytes(b"spectrum")
    (out / "twoframe_components.fits").write_bytes(b"basis")
    run = repo / "outputs/X/2-7"
    (run / "figs").mkdir(parents=True)
    os.symlink(out / "corrected", run / "corrected")
    os.symlink(out / "twoframe_components.fits", run / "twoframe_components.fits")
    (run / "fit.npz").write_bytes(b"fit")
    (run / "figs/one.pdf").write_bytes(b"pdf")
    storage.link_dir(config, "outputs/X/2-7")
    assert os.path.islink(run)
    for rel, content in [("corrected/a.fits", b"spectrum"),
                         ("twoframe_components.fits", b"basis"),
                         ("fit.npz", b"fit"), ("figs/one.pdf", b"pdf")]:
        assert (run / rel).read_bytes() == content, rel
        assert (out / rel).read_bytes() == content, rel
    assert not os.path.islink(out / "corrected"), "a link went onto the disk"
    assert not os.path.lexists(str(run) + ".__moving__")


def test_a_link_inside_that_points_elsewhere_stops_and_moves_nothing(tree):
    repo, disk, config = tree
    storage.check(config)
    run = repo / "outputs/X/2-7"
    run.mkdir(parents=True)
    (run / "fit.npz").write_bytes(b"fit")
    os.symlink(repo, run / "stray")
    with pytest.raises(SystemExit):
        storage.link_dir(config, "outputs/X/2-7")
    assert not os.path.islink(run) and (run / "fit.npz").read_bytes() == b"fit"


def test_a_link_to_somewhere_else_is_not_overruled(tree):
    repo, disk, config = tree
    storage.check(config)
    elsewhere = disk / "elsewhere"
    elsewhere.mkdir()
    (repo / "lbl").mkdir()
    os.symlink(elsewhere, repo / "lbl/models")
    with pytest.raises(SystemExit):
        storage.link_dir(config, "lbl/models")


def test_nothing_is_linked_inside_a_folder_already_out_there(tree):
    repo, disk, config = tree
    storage.check(config)
    (disk / "pca2d/outputs").mkdir(parents=True)
    os.symlink(disk / "pca2d/outputs", repo / "outputs")
    assert storage.link_dir(config, "outputs/X/1-3") == "outputs/X/1-3"
    assert not os.path.lexists(disk / "pca2d/outputs/X/1-3")


def test_an_unmounted_disk_stops_the_run_and_a_dry_run_only_warns(tmp_path):
    config = {"output": {"fits_directory": str(tmp_path / "Volumes/irrisor/pca2d")}}
    assert storage.check(config, dry_run=True) is None
    with pytest.raises(SystemExit):
        storage.check(config)
    assert not (tmp_path / "Volumes").exists(), "the missing parent was created"


def test_the_directory_itself_is_created_when_its_disk_is_there(tmp_path):
    config = {"output": {"fits_directory": str(tmp_path / "pca2d")}}
    assert storage.check(config, dry_run=True)
    assert not (tmp_path / "pca2d").exists(), "a dry run created it"
    assert storage.check(config) and (tmp_path / "pca2d").is_dir()


def test_without_a_fits_directory_everything_stays_local(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config = {"output": {"fits_directory": None}}
    assert storage.check(config) is None
    assert storage.link_dir(config, "outputs/X/1-3") == "outputs/X/1-3"
    assert not os.path.lexists(tmp_path / "outputs")


def test_every_lbl_folder_but_science_is_linked():
    assert "science" not in storage.LBL_FOLDERS
    assert {"lblrv", "lblrdb", "lblreftable", "templates", "masks",
            "calib"} <= set(storage.LBL_FOLDERS)
