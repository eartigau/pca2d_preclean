"""A fit too big for memory must give what a fit in memory gives.

The arrays go onto a mapped file instead of into RAM, and the kernel pages them:
the same numbers, in the same order, through the same code. If that is not true
then the escape hatch is a second pipeline, and two pipelines that disagree are
worse than one that refuses.
"""

from __future__ import annotations

import numpy as np
import pytest

from pca2d import twoframe
from pca2d.twoframe import main, memory_needed, spill_dir, workspace

from test_fit_end_to_end import cube            # noqa: F401  (the fixture)


def fitted(path):
    """What a run left: the two bases and the coefficients."""
    blob = np.load(path, allow_pickle=True)
    return {k: np.asarray(blob[k]) for k in ("P", "Q", "a", "b")}


def test_a_spilled_fit_is_the_same_fit(cube, tmp_path):
    """The whole point, end to end, on the same cube twice."""
    in_ram = tmp_path / "ram"
    on_disk = tmp_path / "disk"
    main(["--cube", cube, "--outdir", str(in_ram), "--iters", "3",
          "-k", "2", "-j", "2", "--max-mad", "0"])

    # the same run, with every array of (rows, columns) mapped onto a file
    kept = twoframe.spill_dir
    try:
        twoframe.spill_dir = lambda *a, **k: str(tmp_path / "spill")
        main(["--cube", cube, "--outdir", str(on_disk), "--iters", "3",
              "-k", "2", "-j", "2", "--max-mad", "0"])
    finally:
        twoframe.spill_dir = kept

    ram, disk = fitted(str(in_ram / "fit.npz")), fitted(str(on_disk / "fit.npz"))
    for name in ("P", "Q", "a", "b"):
        assert ram[name].shape == disk[name].shape, name
        # a basis vector is defined up to its sign, and nothing here fixes one
        for i in range(ram[name].shape[0] if ram[name].ndim > 1 else 1):
            pair = (ram[name], disk[name])
            same = np.allclose(pair[0], pair[1], rtol=1e-8, atol=1e-10)
            flipped = np.allclose(pair[0], -pair[1], rtol=1e-8, atol=1e-10)
            assert same or flipped, "%s differs between memory and file" % name
            break


def test_a_workspace_on_disk_behaves_as_an_array(tmp_path):
    ram = workspace((4, 9), np.float64)
    mapped = workspace((4, 9), np.float64, spill=str(tmp_path))
    assert ram.shape == mapped.shape == (4, 9)
    assert not np.any(ram) and not np.any(mapped), "both start at zero"
    mapped[1, 2] = 3.5
    assert mapped[1, 2] == 3.5 and mapped.sum() == 3.5
    assert isinstance(mapped, np.memmap) and mapped.filename is not None


def test_the_choice_is_half_the_machine_and_is_said(tmp_path, capsys):
    """Memory when it fits, a file when it does not, and never silence."""
    from pca2d.machine import total_ram_bytes

    if not total_ram_bytes():
        pytest.skip("this machine will not say how much memory it has")
    cube_path = str(tmp_path / "cube_tfits_x")

    assert spill_dir(10, 1000, "float32", cube_path, str(tmp_path)) is None
    small = capsys.readouterr().out
    assert "held in memory" in small

    huge = int(10 * total_ram_bytes() / (2 * 4 + 8) / 1000)
    where = spill_dir(huge, 1000, "float32", cube_path, str(tmp_path))
    said = capsys.readouterr().out
    assert where and where.endswith("spill")
    assert "mapped file" in said and "coadd nights" in said, \
        "it says what it is doing and what would avoid it"
    assert "%.1f" % memory_needed(huge, 1000, "float32") in said
