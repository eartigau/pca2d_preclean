"""What pca2d left on the disks, and what the Nettoyage tab may delete.

The dangerous part of a cleanup button is not the deleting, it is the counting:
this tree is full of symlinks, half of them made by storage.link_dir to put a
folder on another disk (those bytes ARE pca2d's) and half of them made by
lbl.link_spectra to point at the observatory's spectra (those are not, and
deleting them frees nothing). Getting the two the wrong way round either offers
to free 25000 raw files or reports a 22 GB tree as empty.
"""

import os

import pytest

from pca2d.housekeeping import (EXPENSIVE, REBUILDABLE, RESULTS, SCRATCH,
                                empty, human, purge, report, survey, totals,
                                tree_size)


def fill(path, n, size):
    os.makedirs(path, exist_ok=True)
    for k in range(n):
        with open(os.path.join(path, "f%d.bin" % k), "wb") as fh:
            fh.write(b"\0" * size)


def test_a_link_inside_a_folder_is_not_counted_as_its_own(tmp_path):
    """lbl/science is 25000 links to the observatory's spectra. Counting those
    would put the raw data on a cleanup button."""
    theirs = tmp_path / "observatory"
    fill(str(theirs), 4, 1000)
    science = tmp_path / "lbl" / "science"
    science.mkdir(parents=True)
    for k in range(4):
        (science / ("f%d.fits" % k)).symlink_to(theirs / ("f%d.bin" % k))

    size, files = tree_size(str(science))
    assert size == 0, "links weigh nothing, because deleting them frees nothing"
    assert files == 4, "but they are still counted, so the line is not blank"
    assert all((theirs / ("f%d.bin" % k)).exists() for k in range(4))


def test_a_folder_that_is_itself_a_link_is_counted(tmp_path):
    """storage.link_dir puts lbl/lblrv on the external disk and leaves a link
    behind. Those 22 GB are pca2d's own, and read as 0 B until 2026-09-15."""
    external = tmp_path / "irrisor" / "lblrv"
    fill(str(external), 5, 2000)
    tree = tmp_path / "lbl"
    tree.mkdir()
    (tree / "lblrv").symlink_to(external, target_is_directory=True)

    size, files = tree_size(str(tree / "lblrv"))
    assert size == 5 * 2000 and files == 5


def test_the_spill_is_not_counted_twice(tmp_path):
    """twoframe.spill_dir puts it inside the cache, so a tab that counted both
    would claim more than the disk holds."""
    cache = tmp_path / "cache"
    fill(str(cache / "cube_tfits_abc"), 3, 1000)
    fill(str(cache / "spill"), 2, 5000)
    items = {it["name"]: it for it in survey(
        {"output": {"cache_directory": str(cache)}}, str(tmp_path / "config.yaml"))}

    assert items["cube cache"]["bytes"] == 3000
    assert items["fit spill files"]["bytes"] == 10000
    assert totals(list(items.values()))[0] == 13000


def test_results_are_counted_and_never_offered(tmp_path):
    fill(str(tmp_path / "out" / "GJ1" / "0-7"), 2, 4000)
    fill(str(tmp_path / "cache"), 1, 1000)
    items = survey({"output": {"cache_directory": str(tmp_path / "cache"),
                               "directory": str(tmp_path / "out")}},
                   str(tmp_path / "config.yaml"))
    by = {it["name"]: it for it in items}
    assert by["reports and corrected spectra"]["bytes"] == 8000
    assert by["reports and corrected spectra"]["removable"] is False
    assert by["cube cache"]["removable"] is True
    total, free = totals(items)
    assert total == 9000 and free == 1000

    freed, gone = purge(items)
    assert freed == 1000, "the results are not in it"
    assert (tmp_path / "out" / "GJ1" / "0-7" / "f0.bin").exists()
    assert not (tmp_path / "cache" / "f0.bin").exists()


def test_a_kept_item_handed_to_purge_is_still_kept(tmp_path):
    """purge trusts the survey's own flag, not its caller's list."""
    fill(str(tmp_path / "results"), 2, 3000)
    item = dict(name="results", path=str(tmp_path / "results"), kind=RESULTS,
                what="", bytes=6000, files=2, removable=False)
    freed, gone = purge([item])
    assert (freed, gone) == (0, [])
    assert (tmp_path / "results" / "f0.bin").exists()


def test_purging_empties_a_folder_and_leaves_it_there(tmp_path):
    """The next run expects the folder; and when it is a link to another disk,
    deleting the link would free nothing and break the layout."""
    external = tmp_path / "irrisor" / "plots"
    fill(str(external), 3, 1000)
    tree = tmp_path / "lbl"
    tree.mkdir()
    link = tree / "plots"
    link.symlink_to(external, target_is_directory=True)

    empty(str(link))
    assert link.is_symlink() and link.is_dir(), "the link is still a folder"
    assert external.exists() and not list(external.iterdir()), "emptied, not removed"


def test_the_dry_run_deletes_nothing_but_says_the_same_number(tmp_path):
    fill(str(tmp_path / "cache"), 4, 2500)
    items = survey({"output": {"cache_directory": str(tmp_path / "cache")}},
                   str(tmp_path / "config.yaml"))
    freed, gone = purge(items, dry_run=True)
    assert freed == 10000 and len(gone) == 1
    assert (tmp_path / "cache" / "f0.bin").exists()
    assert purge(items)[0] == 10000


def test_relative_paths_hang_off_the_config_not_the_working_directory(tmp_path):
    """The window launches a run with its cwd at the config's folder, so that
    is where `cache` and `lbl` are; resolving them against the shell's cwd is
    how a second clone became the LBL tree of a run (2026-09-15)."""
    fill(str(tmp_path / "cache"), 2, 1500)
    items = {it["name"]: it for it in survey({}, str(tmp_path / "config.yaml"))}
    assert items["cube cache"]["path"] == str(tmp_path / "cache")
    assert items["cube cache"]["bytes"] == 3000


def test_the_report_adds_up_and_names_what_is_kept(tmp_path):
    fill(str(tmp_path / "cache"), 1, 2_000_000)
    fill(str(tmp_path / "out"), 1, 1_000_000)
    items = survey({"output": {"cache_directory": str(tmp_path / "cache"),
                               "directory": str(tmp_path / "out")}},
                   str(tmp_path / "config.yaml"))
    text = report(items)
    assert "in all" in text and "can be freed" in text
    assert "(kept)" in text
    assert "3.0 MB" in text and "2.0 MB" in text


@pytest.mark.parametrize("n,said", [(0, "0 B"), (999, "999 B"), (1500, "1.5 kB"),
                                    (2_500_000, "2.5 MB"), (3_400_000_000, "3.4 GB"),
                                    (1.2e12, "1.2 TB")])
def test_sizes_read_the_way_a_person_says_them(n, said):
    assert human(n) == said


def test_a_missing_tree_is_a_zero_and_not_a_crash(tmp_path):
    items = survey({"output": {"cache_directory": str(tmp_path / "nothing"),
                               "directory": str(tmp_path / "neither")}},
                   str(tmp_path / "config.yaml"))
    assert all(it["bytes"] == 0 for it in items)
    assert purge(items) == (0, [])
    assert tree_size(None) == (0, 0)
    assert tree_size(str(tmp_path / "missing")) == (0, 0)


def test_every_lbl_folder_is_classified(tmp_path):
    from pca2d.storage import LBL_FOLDERS

    tree = tmp_path / "lbl"
    for name in LBL_FOLDERS:
        fill(str(tree / name), 1, 100)
    items = {it["name"]: it for it in survey(
        {"lbl": {"directory": str(tree)}}, str(tmp_path / "config.yaml"))}
    for name in LBL_FOLDERS:
        assert "LBL " + name in items, "a folder pca2d makes and nobody counts"
    assert items["LBL lblrdb"]["kind"] == RESULTS, "the velocities stay"
    assert items["LBL lblrv"]["kind"] == EXPENSIVE
    assert items["LBL log"]["kind"] == REBUILDABLE
    assert items["fit spill files"]["kind"] == SCRATCH
