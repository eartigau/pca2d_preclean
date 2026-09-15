#!/usr/bin/env python
"""What pca2d has left on the disks, and which of it can go.

    python -m pca2d.housekeeping --config config.yaml           # just look
    python -m pca2d.housekeeping --config config.yaml --purge   # and delete

A run leaves four kinds of thing, and only the reader knows which is which:

    SCRATCH      the memmap spill of a fit that did not fit in memory. Nothing
                 reads it once the fit has ended; it is only still there when a
                 run was interrupted.
    REBUILDABLE  the cube cache, __pycache__, LBL's plots and logs. Made again
                 from what is still on disk, at a cost in time and nothing else.
    EXPENSIVE    fits, LBL's per-line tables and templates. Also rebuildable,
                 but by re-running hours of work.
    RESULTS      the corrected spectra, the reports, the rdb. What the whole
                 thing was for.

Only the first two are offered for deletion, and each item says what deleting
it would cost. Nothing here follows a symlink: lbl/science is a tree of links
to the spectra, and a cleanup that walked into it would count the raw data as
its own and then delete it.
"""

from __future__ import annotations

import argparse
import os
import shutil

from .logger import log

#: what is offered for deletion, and what is only counted
SCRATCH, REBUILDABLE, EXPENSIVE, RESULTS = "scratch", "rebuildable", "expensive", "results"
#: the kinds the purge button touches, in the order it would remove them
REMOVABLE = (SCRATCH, REBUILDABLE)


def human(n):
    """A byte count as the window shows it."""
    if n is None:
        return "?"
    for unit, scale in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= scale:
            return "%.1f %s" % (n / scale, unit)
    return "%d B" % n


def tree_size(path, exclude=()):
    """(bytes, files) under `path`, never following a symlink.

    A symlinked file counts as zero bytes, which is what deleting it would
    actually free: lbl/science is thousands of links to spectra that belong to
    somebody else, and a cleanup that counted those as its own would offer to
    free the raw data.

    `exclude` are paths counted by another item: the fit spill sits inside the
    cube cache, and a tab whose lines add up to more than the disk holds is a
    tab nobody believes.
    """
    if not path or not os.path.exists(path):
        return 0, 0
    if os.path.islink(path) and os.path.isdir(path):
        # the item's OWN path being a link means the whole folder was moved to
        # another disk (storage.link_dir): those bytes are pca2d's and get
        # counted. A link found INSIDE a folder, below, is somebody else's file
        path = os.path.realpath(path)
    if os.path.islink(path):
        return 0, 1
    if os.path.isfile(path):
        try:
            return os.path.getsize(path), 1
        except OSError:
            return 0, 0
    skip = {os.path.abspath(p) for p in exclude}
    total, count = 0, 0
    for here, dirs, names in os.walk(path, followlinks=False):
        dirs[:] = [d for d in dirs
                   if not os.path.islink(os.path.join(here, d))
                   and os.path.abspath(os.path.join(here, d)) not in skip]
        for name in names:
            full = os.path.join(here, name)
            if os.path.islink(full):
                count += 1
                continue
            try:
                total += os.path.getsize(full)
            except OSError:
                continue
            count += 1
    return total, count


class Item(dict):
    """One line of the tab: a name, a path, what it is, and what it costs."""

    def __init__(self, name, path, kind, what, exclude=(), key=None):
        size, files = tree_size(path, exclude)
        super().__init__(name=name, path=path, kind=kind, what=what,
                         # a stable slug, because `name` is what the window
                         # PRINTS and the window prints it in two languages
                         key=key or name, bytes=size, files=files,
                         removable=kind in REMOVABLE)

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key) from None


def work_dir(config_path):
    """Where the relative paths of a config resolve: beside the config file."""
    if not config_path:
        return os.getcwd()
    return os.path.dirname(os.path.abspath(config_path)) or os.getcwd()


def under(base, path):
    """`path` as the config means it: relative ones hang off the config's own
    folder, which is how a run launched from the window resolves them."""
    if not path:
        return None
    return path if os.path.isabs(path) else os.path.join(base, path)


def survey(config=None, config_path=None, out_root=None, package=None):
    """Every place pca2d puts bytes, measured, biggest kind first.

    `config` is a resolved or raw configuration (only output.cache_directory,
    output.directory and lbl.directory are read, all optional); `out_root`
    overrides output.directory, which is what the window passes because the
    window's field is the truth while it is open.
    """
    config = config or {}
    base = work_dir(config_path)
    out = config.get("output") or {}
    lbl = config.get("lbl") or {}
    items = []

    cache = under(base, out.get("cache_directory") or "cache")
    # the spill sits beside the cube it belonged to (twoframe.spill_dir), so it
    # is INSIDE the cache and has to come out of the cache's own number
    spill = os.path.join(cache, "spill") if cache else None
    items.append(Item("cube cache", cache, REBUILDABLE,
                      "the cubes, read back instead of re-reading the spectra."
                      " Rebuilt on the next run, minutes per campaign",
                      exclude=[spill] if spill else (), key="cache"))
    items.append(Item("fit spill files", spill, SCRATCH,
                      "the mapped scratch of a fit too big for memory. Nothing"
                      " reads it once the fit has ended", key="spill"))

    root = out_root or out.get("directory")
    if root:
        items.append(Item("reports and corrected spectra", root, RESULTS,
                          "what the runs produced. Not offered for deletion",
                          key="results"))

    tree = under(base, lbl.get("directory") or "lbl")
    if tree and os.path.isdir(tree):
        for folder, kind, what in (
                ("science", REBUILDABLE,
                 "links to the spectra LBL measures. Remade by the lbl stage;"
                 " symlinks, so this is rarely more than a few MB"),
                ("plots", REBUILDABLE, "LBL's own figures, one per exposure"),
                ("log", REBUILDABLE, "LBL's logs"),
                ("lblrv", EXPENSIVE,
                 "LBL's per-line velocity tables, one per exposure. Remade only"
                 " by running LBL again"),
                ("templates", EXPENSIVE, "the templates LBL built"),
                ("masks", EXPENSIVE, "the line masks LBL built"),
                ("models", EXPENSIVE, "LBL's models"),
                ("calib", EXPENSIVE, "LBL's calibrations"),
                ("lblreftable", EXPENSIVE, "LBL's reference tables"),
                ("lblrdb", RESULTS, "the velocities. Not offered for deletion")):
            items.append(Item("LBL " + folder, os.path.join(tree, folder),
                              kind, what, key="lbl_" + folder))

    if package:
        items.append(Item("__pycache__", None, REBUILDABLE, "", key="pycache"))
        items[-1].update(_pycache(package))

    items.sort(key=lambda it: (-it["bytes"], it["name"]))
    return items


def _pycache(package):
    """Every __pycache__ under the package, as one item's worth of numbers."""
    total, files, paths = 0, 0, []
    for here, dirs, _names in os.walk(package, followlinks=False):
        for d in list(dirs):
            if d == "__pycache__":
                full = os.path.join(here, d)
                size, n = tree_size(full)
                total += size
                files += n
                paths.append(full)
                dirs.remove(d)
    return {"bytes": total, "files": files, "path": paths or None,
            "what": "compiled Python, remade the next time it is imported"}


def totals(items):
    """(all bytes, bytes that can be freed) of a survey."""
    return (sum(it["bytes"] for it in items),
            sum(it["bytes"] for it in items if it["removable"]))


def empty(path):
    """Remove everything INSIDE a folder, keeping the folder itself.

    Not rmtree: half these folders are symlinks onto another disk made before
    the run started (storage.link_dir), and a cleanup that deleted the link
    would free nothing and break the layout, while one that followed it and
    deleted the target would leave a link pointing at nothing. Emptying works
    the same either way, and what the next run expects to find is a folder.
    """
    for name in os.listdir(path):
        full = os.path.join(path, name)
        try:
            if os.path.isdir(full) and not os.path.islink(full):
                shutil.rmtree(full)
            else:
                os.remove(full)
        except OSError as exc:
            log("could not remove %s (%s)" % (full, exc), "warn")


def purge(items, dry_run=False):
    """Delete the removable items of a survey. Returns (bytes freed, paths).

    Only items the survey itself marked removable are touched, so a caller
    cannot turn a results folder into a purge by passing it here.
    """
    freed, gone = 0, []
    for it in items:
        if not it.get("removable"):
            continue
        paths = it["path"] if isinstance(it["path"], list) else [it["path"]]
        for path in paths:
            if not path or not os.path.exists(path):
                continue
            size, _files = tree_size(path)
            if not dry_run:
                try:
                    if os.path.isdir(path):
                        empty(path)
                    else:
                        os.remove(path)
                except OSError as exc:
                    log("could not remove %s (%s)" % (path, exc), "warn")
                    continue
            freed += size
            gone.append(path)
    return freed, gone


def report(items):
    """The survey as the lines the tab and the command line both show."""
    width = max([len(it["name"]) for it in items] + [12])
    out = []
    for it in items:
        out.append("%-*s  %9s  %7s files  %-12s %s"
                   % (width, it["name"], human(it["bytes"]), it["files"],
                      it["kind"], "" if it["removable"] else "(kept)"))
    total, free = totals(items)
    out.append("%-*s  %9s" % (width, "in all", human(total)))
    out.append("%-*s  %9s" % (width, "can be freed", human(free)))
    return "\n".join(out)


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="python -m pca2d.housekeeping",
                                description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--out-dir", default=None,
                   help="override output.directory")
    p.add_argument("--purge", action="store_true",
                   help="delete the scratch and rebuildable items")
    p.add_argument("--dry-run", action="store_true",
                   help="with --purge: say what would go, delete nothing")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    from .config import load_config
    config = load_config(args.config) if os.path.exists(args.config) else {}
    items = survey(config, args.config, args.out_dir,
                   package=os.path.dirname(os.path.abspath(__file__)))
    log("what pca2d has left on these disks:", "info")
    for line in report(items).splitlines():
        log("  " + line, "value")
    if not args.purge:
        return 0
    freed, gone = purge(items, dry_run=args.dry_run)
    log("%s %s from %d places"
        % ("would free" if args.dry_run else "freed", human(freed), len(gone)),
        "info")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
