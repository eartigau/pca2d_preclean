"""The fit's steps taken a block at a time give the whole-array numbers.

On 2026-09-16 a three-object joint fit (2.8 GB of float32 data, mapped onto a
file because it would not fit in memory) still held 11 GB in memory and put
the machine 21 of 22 GB into swap: a copy of the weights, the data copied back
out of its file, and in every sweep the clip's float64 array of the cube's size
with two more copies of it for the medians. Those steps now work a block of
rows or columns at a time and store the clipped weights as what changed.

A change of that kind is only acceptable if it changes nothing else, so each
step is held here against the code it replaced, copied verbatim, with blocks
made small enough that every boundary is crossed.
"""

from __future__ import annotations

import os
import subprocess
import sys
import warnings

import numpy as np
import pytest

from pca2d import twoframe
from pca2d.twoframe import (ClippedWeights, LanczosShifter, clip_weights,
                            column_sums, parity_means, star_frame_template,
                            subtract_means, sweep_spill, weighted_power,
                            weighted_power_rows, workspace, zero_unweighted)


# ---------------------------------------------------------------- the old code
def _old_nanmedian(z):
    return twoframe._bn.nanmedian(np.ascontiguousarray(z.T), axis=1) \
        if twoframe._bn is not None else np.nanmedian(z, axis=0)


def old_clip_weights(w0, residual, clip=3.0, min_spectra=20):
    good = w0 > 0
    z = np.sqrt(w0, dtype=np.float64)
    z *= residual
    z[~good] = np.nan
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        centre = _old_nanmedian(z)
        z -= centre[None, :]
        np.abs(z, out=z)
        scale = 1.4826 * _old_nanmedian(z)
        enough = np.sum(good, axis=0) >= min_spectra
        fallback = np.nanmedian(scale[enough]) if np.any(enough) else 1.0
    scale = np.where(enough & np.isfinite(scale) & (scale > 0), scale, fallback)
    if not np.isfinite(fallback) or fallback <= 0:
        scale = np.ones_like(scale)
    with np.errstate(invalid="ignore", divide="ignore"):
        z /= scale[None, :]
    hit = float(np.count_nonzero(z[good] > clip) / max(np.count_nonzero(good), 1))
    with np.errstate(invalid="ignore", divide="ignore"):
        np.maximum(z, clip, out=z)
        z **= -2.0
        z *= clip ** 2
    z[~np.isfinite(z)] = 1.0
    z *= w0
    return z, hit


def old_hierarchical_median(home, berv, bin_kms, min_entries):
    labels = np.floor((berv - np.nanmin(berv)) / bin_kms).astype(int)
    per_bin = []
    for value in np.unique(labels):
        rows = labels == value
        if rows.sum() < min_entries:
            continue
        per_bin.append(np.nanmedian(home[rows], axis=0))
    if len(per_bin) < 2:
        return np.nanmedian(home, axis=0)
    return np.nanmedian(np.vstack(per_bin), axis=0)


def old_star_frame_template(data, w, shifter, delta, min_spectra=20, berv=None,
                            berv_bin=None, berv_min_entries=3):
    home = shifter.rows(data, -delta)
    home_w = shifter.rows(w, -delta)
    home[home_w <= 0] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if berv_bin and berv is not None:
            template = old_hierarchical_median(
                home, np.asarray(berv, dtype=float), float(berv_bin) / 1000.0,
                int(berv_min_entries))
        else:
            template = np.nanmedian(home, axis=0)
    count = np.sum(np.isfinite(home), axis=0)
    template = np.where(np.isfinite(template), template, 0.0)
    template[count < min_spectra] = 0.0
    return template


def old_parity_means(data, w, parity):
    groups = np.unique(parity)
    means = np.zeros((groups.size, data.shape[1]))
    for i, value in enumerate(groups):
        rows = parity == value
        total = w[rows].sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            means[i] = np.where(total > 0,
                                (w[rows] * data[rows]).sum(axis=0)
                                / np.where(total > 0, total, 1.0), 0.0)
    return means, groups


def old_subtract_means(data, means, group):
    if means.shape[0] == 1:
        data -= means[0][None, :]
        return data
    for i in range(means.shape[0]):
        rows = group == i
        data[rows] -= means[i][None, :]
    return data


# ------------------------------------------------------------------- fixtures
@pytest.fixture
def small_blocks(monkeypatch):
    """Blocks of a few hundred samples, so every loop crosses boundaries."""
    monkeypatch.setattr(twoframe, "_BLOCK_SAMPLES", 700)


def a_cube(rows=60, cols=400, dtype=np.float32, seed=5):
    r = np.random.default_rng(seed)
    data = r.normal(size=(rows, cols)).astype(dtype)
    w = (r.random((rows, cols)) * 4 + 0.5).astype(dtype)
    w[r.random((rows, cols)) < 0.1] = 0.0          # holes
    w[:, 7] = 0.0                                  # a dead column
    w[:3, 11] = 0.0
    w[:, 12][r.random(rows) < 0.8] = 0.0           # a poorly covered column
    data[w == 0] = 0.0
    return data, w


def identical(a, b):
    return np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True)


# ---------------------------------------------------------------------- tests
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_the_clip_is_the_old_clip(small_blocks, dtype):
    data, w0 = a_cube(dtype=dtype)
    resid = np.random.default_rng(1).normal(size=w0.shape)
    resid[4, 30:60] = 40.0                          # something to clip
    want, want_hit = old_clip_weights(w0, resid)
    got, hit = clip_weights(w0, resid)
    assert isinstance(got, np.ndarray) and identical(got, want)
    assert hit == want_hit and hit > 0


def test_sparse_clipped_weights_read_as_the_dense_ones(small_blocks):
    data, w0 = a_cube()
    resid = np.random.default_rng(2).normal(size=w0.shape)
    resid[9, 100:180] = -30.0
    want, want_hit = old_clip_weights(w0, resid)
    got, hit = clip_weights(w0, resid, sparse=True)
    assert isinstance(got, ClippedWeights)
    assert hit == want_hit
    assert got.shape == want.shape and len(got) == want.shape[0]
    # by row, by block, from the end, and whole
    for n in (0, 9, 10, want.shape[0] - 1):
        assert identical(got[n], want[n])
    assert identical(got[-1], want[-1])
    for start, stop in ((0, 7), (5, 23), (40, 60), (58, 99)):
        assert identical(got[start:stop], want[start:stop])
    assert identical(np.asarray(got), want)
    assert identical(got[:, 5], want[:, 5])        # anything else: assembled
    # and it is small: only what the clip changed is stored
    assert got.nbytes < 0.5 * want.nbytes


def test_a_clip_that_is_not_exactly_one_inside_stays_dense(small_blocks):
    """At some thresholds (clip/clip)^2 rounds off 1; sparse would be larger."""
    data, w0 = a_cube(dtype=np.float64)
    resid = np.random.default_rng(3).normal(size=w0.shape)
    clip = 3.5
    unit = np.maximum(np.full(1, clip), clip) ** -2.0 * clip ** 2
    got, _ = clip_weights(w0, resid, clip=clip, sparse=True)
    want, _ = old_clip_weights(w0, resid, clip=clip)
    assert identical(got, want)
    if unit[0] != 1.0:
        assert isinstance(got, np.ndarray)


@pytest.mark.parametrize("berv_bin", [None, 3000.0])
def test_the_star_spectrum_is_the_old_one(small_blocks, tmp_path, berv_bin):
    data, w = a_cube(rows=48, cols=300)
    r = np.random.default_rng(6)
    delta = r.uniform(-6, 6, size=data.shape[0])
    berv = r.uniform(-20, 20, size=data.shape[0])
    parity = np.arange(data.shape[0]) % 2
    shifter = LanczosShifter(data.shape[1], a=4, max_shift=10)
    for value in (0, 1):
        rows = parity == value
        want = old_star_frame_template(data[rows], w[rows], shifter,
                                       delta[rows], min_spectra=5,
                                       berv=berv[rows], berv_bin=berv_bin)
        for spill in (None, str(tmp_path)):
            got = star_frame_template(data, w, shifter, delta, min_spectra=5,
                                      berv=berv, berv_bin=berv_bin, rows=rows,
                                      spill=spill)
            assert identical(got, want), (value, spill)
    assert not os.listdir(tmp_path), "the carried copy is not left behind"
    # without `rows`, the whole cube as before
    assert identical(
        star_frame_template(data, w, shifter, delta, min_spectra=5),
        old_star_frame_template(data, w, shifter, delta, min_spectra=5))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_the_parity_means_are_the_old_ones(dtype):
    data, w = a_cube(dtype=dtype)
    parity = np.random.default_rng(7).integers(0, 3, size=data.shape[0])
    got, groups = parity_means(data, w, parity)
    want, want_groups = old_parity_means(data, w, parity)
    assert identical(groups, want_groups) and identical(got, want)


def test_subtracting_the_means_is_the_old_subtraction(small_blocks):
    data, w = a_cube()
    r = np.random.default_rng(8)
    group = r.integers(0, 3, size=data.shape[0])
    means = r.normal(size=(3, data.shape[1]))
    for m, g in ((means, group), (means[:1], np.zeros_like(group))):
        want = old_subtract_means(data.copy(), m, g)
        assert identical(subtract_means(data.copy(), m, g), want)
    # zero means leave the data alone, without reading it
    before = data.copy()
    assert identical(subtract_means(data, np.zeros((3, data.shape[1])), group),
                     before)


def test_zeroing_and_scoring_are_the_old_ones(small_blocks):
    data, w = a_cube()
    data[w == 0] = 5.0                              # something to zero
    want = data.copy()
    want[w <= 0] = 0.0
    assert identical(zero_unweighted(data.copy(), w), want)
    assert identical(np.where(w > 0, data, 0.0), want)
    rows = weighted_power_rows(data, w)
    assert identical(rows, np.sum(w * data ** 2, axis=1, dtype=np.float64))
    total = 0.0
    for s in range(0, data.shape[0], 16):
        total += float(np.sum(w[s:s + 16] * data[s:s + 16] ** 2,
                              dtype=np.float64))
    assert weighted_power(data, w, 16) == total


def test_column_sums_are_numpys():
    data, w = a_cube()
    rows = np.flatnonzero(np.arange(data.shape[0]) % 3 == 1)
    count, total = column_sums(
        [(lambda r: (w[r] > 0).astype(np.int64), len(w)),
         (lambda r: w[r], len(w))], rows=rows)
    assert identical(count, (w[rows] > 0).sum(axis=0))
    assert total.dtype == np.float32 and identical(total, w[rows].sum(axis=0))


def test_mapped_files_of_ended_fits_are_removed(tmp_path, capsys):
    host = twoframe._host()
    done = subprocess.run([sys.executable, "-c", "import os; print(os.getpid())"],
                          capture_output=True, text=True, check=True)
    dead = int(done.stdout)
    if twoframe._alive(dead):
        pytest.skip("the finished process's number was taken again")
    names = {
        "ended": "pca2d-spill-work-%s-%d.dat" % (host, dead),
        "this run": "pca2d-spill-work-%s-%d.dat" % (host, os.getpid()),
        "another machine": "pca2d-spill-work-elsewhere-%d.dat" % dead,
        "not ours": "notes-%d.dat" % dead,
    }
    for name in names.values():
        (tmp_path / name).write_bytes(b"\0" * 1000)
    assert sweep_spill(str(tmp_path)) == 1000
    left = set(os.listdir(tmp_path))
    assert names["ended"] not in left
    assert {names["this run"], names["another machine"],
            names["not ours"]} <= left
    assert "fits that had ended" in capsys.readouterr().out
    assert sweep_spill(str(tmp_path / "absent")) == 0


def test_a_workspace_the_disk_cannot_hold_is_refused(tmp_path, monkeypatch):
    import shutil
    from collections import namedtuple

    usage = namedtuple("usage", "total used free")
    monkeypatch.setattr(shutil, "disk_usage", lambda path: usage(100, 90, 10))
    with pytest.raises(SystemExit, match="GB on disk"):
        workspace((10, 10), np.float64, spill=str(tmp_path), tag="work")
    monkeypatch.undo()
    mapped = workspace((10, 10), np.float64, spill=str(tmp_path), tag="work")
    assert twoframe._host() in os.path.basename(mapped.filename)
    assert not np.any(mapped)
