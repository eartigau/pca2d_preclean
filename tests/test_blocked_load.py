"""Reading a cube in column blocks must give what reading it whole gave.

The peak of a run is while the cube is read, not while it is fitted: all at
once, four arrays of its size are live (data, sigma, transmission, weights),
and on 2026-09-15 a ten-object joint cube asked 38 GB of a 17 GB machine and
was killed by the kernel with nothing said. Block by block the peak is the data
and the weights, which is half of that. The arithmetic per column is unchanged,
so the two must agree to the bit.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from pca2d import twoframe


def a_cube(tmp_path, rows=7, cols=400, seed=3):
    """A cube on disk, with everything the reader has to cope with in it."""
    rng = np.random.default_rng(seed)
    data = rng.normal(0.0, 0.1, size=(rows, cols)).astype(np.float32)
    sigma = np.abs(rng.normal(0.05, 0.01, size=(rows, cols))).astype(np.float32)
    trans = np.clip(rng.uniform(0.2, 1.0, size=(rows, cols)), 0, 1).astype(np.float32)
    data[0, 10:20] = np.nan                 # samples that were never written
    sigma[1, 30:40] = 0.0                   # and samples with no error bar
    sigma[2, 50] = np.nan
    data[3, 60:65] = -3.0                   # deep enough to be clipped away
    path = tmp_path / "cube_tfits_test"
    path.mkdir()
    np.save(path / "grid.npy", np.linspace(1000.0, 1001.0, cols))
    np.save(path / "data.npy", data)
    np.save(path / "sigma.npy", sigma)
    np.save(path / "trans.npy", trans)
    Table({"snr_band": np.full(rows, 100.0),
           "berv": np.linspace(-10, 10, rows),
           "object": ["X"] * rows,
           "filename": ["%04dt.fits" % i for i in range(rows)]}
          ).write(path / "meta.fits")
    return str(path)


def read_whole(path, **kw):
    """What load_cube did before the blocks, written out here as the reference."""
    grid, data, sigma, trans, meta = twoframe.read_cube_files(
        path, np.float32, None, lazy=False)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(np.isfinite(sigma) & (sigma > 0), 1.0 / sigma ** 2, 0.0)
    w[~np.isfinite(data)] = 0.0
    w[data < kw.get("ln_clip_low", -0.5)] = 0.0
    if trans is not None:
        ramp = kw.get("ramp_zero", 0.5)
        w *= np.clip((trans - ramp) / (1.0 - ramp), 0.0, 1.0)
    w[:, :twoframe.EDGE] = 0.0
    w[:, -twoframe.EDGE:] = 0.0
    data = np.where(w > 0, data, 0.0)
    return grid, data, w


def test_block_by_block_is_the_same_cube(tmp_path):
    path = a_cube(tmp_path)
    grid, data, w, _meta = twoframe.load_cube(path, min_snr_frac=0,
                                              dtype=np.float32)
    g0, d0, w0 = read_whole(path)
    assert np.array_equal(grid, g0)
    assert np.array_equal(w, w0), "the weights, to the bit"
    assert np.array_equal(data, d0), "and the data, to the bit"
    assert np.isfinite(w).all() and (w >= 0).all()


def test_a_block_smaller_than_the_cube_changes_nothing(tmp_path):
    """The block size is a memory choice, never a numerical one."""
    path = a_cube(tmp_path, rows=5, cols=311)
    _g, d0, w0, _m = twoframe.load_cube(path, min_snr_frac=0, dtype=np.float32)
    kept = twoframe._LOAD_BLOCK_BYTES
    try:
        twoframe._LOAD_BLOCK_BYTES = 64      # a few columns at a time
        _g, d1, w1, _m = twoframe.load_cube(path, min_snr_frac=0,
                                            dtype=np.float32)
    finally:
        twoframe._LOAD_BLOCK_BYTES = kept
    assert np.array_equal(w0, w1) and np.array_equal(d0, d1)


def test_dropping_rows_still_gives_what_it_gave(tmp_path):
    """The relative SNR cut runs before the weights are built, on both paths."""
    path = a_cube(tmp_path, rows=9)
    table = Table.read(path + "/meta.fits")
    table["snr_band"] = [100.0, 100.0, 5.0, 100.0, 100.0, 4.0, 100.0, 100.0, 100.0]
    table.write(path + "/meta.fits", overwrite=True)
    _g, data, w, meta = twoframe.load_cube(path, min_snr_frac=0.5,
                                           dtype=np.float32)
    assert len(meta) == 7 and data.shape[0] == 7 and w.shape[0] == 7


def test_the_fit_weighs_itself_against_half_the_machine():
    """A run that cannot fit must say so in a second, not be killed in forty
    minutes. Half the machine, never more: the other half is the operating
    system, the file cache this cube is read through, and LBL if it is running
    beside it."""
    import pytest

    from pca2d.machine import DEFAULT_FRACTION, total_ram_bytes
    from pca2d.twoframe import check_memory, memory_needed

    assert DEFAULT_FRACTION == 0.5
    total = total_ram_bytes()
    if not total:
        pytest.skip("this machine will not say how much memory it has")
    budget = total * DEFAULT_FRACTION

    # one that fits in a tenth of the budget, and one ten times the machine
    small = int(budget * 0.1 / (2 * 1.1 * 4 * 1000))
    check_memory(max(1, small), 1000, "float32")
    huge = int(10 * total / (2 * 1.1 * 4 * 1000))
    with pytest.raises(MemoryError) as raised:
        check_memory(huge, 1000, "float32")
    said = str(raised.value)
    assert "coadd nights" in said and "narrow the domain" in said, \
        "it says what to do about it, not only that it will not"
    assert "%.1f" % memory_needed(huge, 1000, "float32") in said, \
        "and how much it would have taken"
