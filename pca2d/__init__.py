"""pca2d-preclean: the window, and the LBL wrap file it runs.

The two-frame precleaning itself lives in LBL since 2026-10-19, as the module
lbl.pca2d, which lbl_wrap calls when its run parameters carry a PCA2D dict
(the pca2d-module branch of https://github.com/njcuk9999/lbl). This package
keeps what is not the algorithm:

    gui            the window: objects, settings, the run's log
    gui_es, gui_pt its Spanish and Portuguese
    scan           what a data root holds, read once and remembered
    housekeeping   what the runs leave on the disks, and clearing it
    wrapfile       the LBL wrap script a run is, with its PCA2D dict

Every other pca2d.<module> IS lbl.pca2d.<module>: the same module object
under a second name, so that `python -m pca2d.cli`, the scripts under tools/,
docs/make_figures.py and the tests run on the one copy of the code, the one
LBL runs. Nothing is copied back here, so the two cannot drift apart.

Documentation: https://eartigau.github.io/pca2d_preclean/
"""
from __future__ import annotations

import importlib
import importlib.abc
import importlib.util
import os
import sys

__version__ = "0.4.0"

#: the modules that are this package's own; every other one is LBL's
OWN = frozenset({"gui", "gui_es", "gui_pt", "scan", "housekeeping", "wrapfile"})
#: where the algorithm is
ENGINE = "lbl.pca2d"
#: what to say when it is not there
MISSING = ("pca2d needs LBL with its pca2d module: the pca2d-module branch of"
           " https://github.com/njcuk9999/lbl, installed in this environment"
           " (git clone -b pca2d-module https://github.com/njcuk9999/lbl.git;"
           " pip install -e ./lbl). See the README, 'Installing it'.")


def _not_lbl(entry: str) -> bool:
    """Whether a sys.path entry holds an lbl/ folder that is not a package.

    An LBL data tree is called lbl/ too (calib, masks, science...), and a
    folder of that name ahead of the installed LBL on the path makes `import
    lbl` an empty namespace package: an editable LBL is found through a hook
    that comes after the path, so the data folder wins. This repository has
    one at its root, and `python -m` puts the working directory first."""
    folder = os.path.join(entry or os.getcwd(), "lbl")
    return os.path.isdir(folder) and not os.path.exists(
        os.path.join(folder, "__init__.py"))


def engine():
    """lbl.pca2d, imported past any lbl/ data folder on the path."""
    try:
        return importlib.import_module(ENGINE)
    except ImportError as first:
        found = sys.modules.get("lbl")
        if found is None or getattr(found, "__file__", None) is not None:
            raise ImportError(MISSING) from first
    # `lbl` resolved to a data folder: forget it, and look again without the
    # entries that hold one. The real package, once imported, stays in
    # sys.modules, so the path can be put back as it was
    for name in [n for n in sys.modules if n == "lbl" or n.startswith("lbl.")]:
        del sys.modules[name]
    saved = list(sys.path)
    sys.path[:] = [p for p in sys.path if not _not_lbl(p)]
    try:
        return importlib.import_module(ENGINE)
    except ImportError as second:
        raise ImportError(MISSING) from second
    finally:
        sys.path[:] = saved


class _Alias(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """pca2d.<name> as lbl.pca2d.<name>: the same module object, never a copy.

    First on sys.meta_path, ahead of an editable install of an older
    pca2d-preclean that would otherwise hand out its own copy of a module
    this package no longer has."""

    prefix = __name__ + "."

    def _target(self, name):
        sub = name[len(self.prefix):]
        if not name.startswith(self.prefix) or sub.split(".")[0] in OWN:
            return None
        return ENGINE + "." + sub

    def find_spec(self, name, path=None, target=None):
        real_name = self._target(name)
        if real_name is None:
            return None
        engine()
        real = importlib.util.find_spec(real_name)
        if real is None:
            return None
        return importlib.util.spec_from_loader(
            name, self, origin=real.origin,
            is_package=real.submodule_search_locations is not None)

    def create_module(self, spec):
        return importlib.import_module(self._target(spec.name))

    def exec_module(self, module):
        """Nothing: the module ran once, under its own name in LBL."""

    def get_code(self, name):
        """What `python -m pca2d.cli` runs: the code of lbl/pca2d/cli.py."""
        real = importlib.util.find_spec(self._target(name))
        return real.loader.get_code(real.name)


if not any(isinstance(f, _Alias) for f in sys.meta_path):
    sys.meta_path.insert(0, _Alias())
