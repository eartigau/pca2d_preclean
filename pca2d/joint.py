"""One observer basis for several stars.

The atmosphere and the instrument belong to the night, not to the target: two
stars observed in the same campaign share them. Only the star differs, and
with it the BERV coverage and the systemic velocity, so fitting several
objects together buys a cleaner observer basis. What is common to all of them
cannot be a star, and what follows one star cannot be common.

Each object keeps its own star spectrum, per order parity; the observer
components are shared by construction, because there is one basis and every
row of every object sees it.

The joint cube is the objects' own cubes, concatenated row for row on the one
grid they share, with three labels rewritten:

    parity     2 * (the object's index) + the order parity. The star spectra
               are estimated per group, so each object gets its own, while
               `parity % 2` still gives the order parity, which means the same
               thing for every object and is what the rest of the package
               reads.
    exposure   offset per object, so the two rows of an exposure stay tied and
               no two objects share an id
    object     the object's name, which the metadata already carries

Nothing else changes: the fit, the figures and the correction read a cube.
"""
from __future__ import annotations

import os
import shutil

import numpy as np
from astropy.table import Table, vstack

from .logger import log

#: ids are offset by this per object, far above any campaign's exposure count
EXPOSURE_STRIDE = 1000000


def joint_name(objects):
    """PROXIMA+GJ1+GJ3090: the folder and the cache key a joint run is under."""
    return "+".join(str(o) for o in objects)


def cube_path(cache_dir, fmt, key, objects):
    """Where the joint cube lives, keyed by the configuration and the objects."""
    import hashlib
    stamp = hashlib.sha1(joint_name(objects).encode()).hexdigest()[:6]
    return os.path.join(cache_dir, "cube_%s_%s_j%s" % (fmt, key, stamp))


def build(cubes, objects, path, config_file=None, chunk=64):
    """Write the joint cube of `cubes`, in that order, at `path`.

    The rows are copied a chunk at a time through memory maps: three arrays of
    a thousand rows over half a million columns do not belong in memory at
    once. Returns the number of rows written.
    """
    grids = [np.load(os.path.join(c, "grid.npy")) for c in cubes]
    for name, grid in zip(objects, grids):
        if grid.shape != grids[0].shape or not np.allclose(grid, grids[0]):
            raise SystemExit("%s was built on another grid: a joint fit needs"
                             " one domain, one step and one high pass for every"
                             " object" % name)
    metas = [Table.read(os.path.join(c, "meta.fits")) for c in cubes]
    rows = [len(m) for m in metas]
    total = int(sum(rows))
    has_trans = all(os.path.exists(os.path.join(c, "trans.npy")) for c in cubes)

    tmp = path + ".tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp, exist_ok=True)
    np.save(os.path.join(tmp, "grid.npy"), grids[0])
    n_pixels = grids[0].size
    names = ["data.npy", "sigma.npy"] + (["trans.npy"] if has_trans else [])
    out = {name: np.lib.format.open_memmap(
        os.path.join(tmp, name), mode="w+", dtype=np.float32,
        shape=(total, n_pixels)) for name in names}
    start = 0
    for cube, name, n in zip(cubes, objects, rows):
        for what in names:
            source = np.load(os.path.join(cube, what), mmap_mode="r")
            for a in range(0, n, chunk):
                b = min(a + chunk, n)
                out[what][start + a:start + b] = source[a:b]
            del source
        log("  %-10s %4d rows from %s" % (name, n, os.path.basename(cube)))
        start += n
    for handle in out.values():
        handle.flush()
    del out

    stacked = []
    for k, (meta, name) in enumerate(zip(metas, objects)):
        meta = meta.copy()
        parity = (np.asarray(meta["parity"], dtype=int)
                  if "parity" in meta.colnames else np.zeros(len(meta), dtype=int))
        meta["parity"] = 2 * k + (parity % 2)
        exposure = (np.asarray(meta["exposure"], dtype=int)
                    if "exposure" in meta.colnames else np.arange(len(meta)))
        meta["exposure"] = exposure + k * EXPOSURE_STRIDE
        meta["object"] = [str(name)] * len(meta)
        stacked.append(meta)
    vstack(stacked, join_type="outer").write(os.path.join(tmp, "meta.fits"),
                                             overwrite=True)
    if config_file and os.path.exists(config_file):
        shutil.copy(config_file, os.path.join(tmp, "cube_config.yaml"))
    shutil.rmtree(path, ignore_errors=True)
    os.rename(tmp, path)
    log("joint cube %s: %d rows, %d objects, %d columns"
        % (path, total, len(objects), n_pixels), "value")
    return total


def star_group(objects, name):
    """The first star-spectrum group of an object: 2 * its index.

    The correction of that object's files adds the order parity to it, so the
    group is 2 * index + parity, exactly what the cube's labels say.
    """
    return 2 * list(objects).index(name)
