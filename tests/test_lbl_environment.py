"""Which LBL measures: lbl.environment, and the interpreter it names.

The speed branch of LBL lives in a conda environment of its own (lbl-rapide),
the default since 2026-09-18, and `current` is the LBL installed beside this
package. None of this needs either one: the environments are folders made
here, and the LBL they import is a stand-in.
"""

import os
import stat
import sys

import pytest

from pca2d import lbl as splbl
from pca2d import lbltemplate
from pca2d.config import load_config


@pytest.fixture(autouse=True)
def _fresh_probe():
    splbl.probe.cache_clear()
    yield
    splbl.probe.cache_clear()


def _python(path, body="exec %s \"$@\"\n" % sys.executable):
    """An executable standing in for an environment's python."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def _nowhere(monkeypatch, tmp_path, this="pca2d-preclean"):
    """This process in <tmp>/envs/<this>, and no conda to ask."""
    monkeypatch.setattr(sys, "prefix", str(tmp_path / "envs" / this))
    monkeypatch.delenv("CONDA_EXE", raising=False)
    monkeypatch.setattr(splbl.shutil, "which", lambda _name: None)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))


def test_the_default_is_the_speed_branch():
    cfg = load_config("config.yaml", instrument="SPIROU")
    assert cfg["lbl"]["environment"] == splbl.DEFAULT_ENVIRONMENT == "lbl-rapide"
    assert splbl.asked_environment({}) == "lbl-rapide"
    assert "test-speed-260918-110104" in splbl.FAST_RECIPE


def test_current_is_this_interpreter():
    python, what = splbl.interpreter({"lbl": {"environment": "current"}})
    assert python == sys.executable and "pca2d-preclean runs in" in what


def test_an_environment_is_found_by_name_beside_this_one(tmp_path, monkeypatch):
    _nowhere(monkeypatch, tmp_path)
    fast = _python(tmp_path / "envs" / "lbl-rapide" / "bin" / "python")
    python, what = splbl.interpreter({"lbl": {"environment": "lbl-rapide"}})
    assert python == str(fast) and "lbl-rapide" in what


def test_an_environment_that_is_not_there_says_so(tmp_path, monkeypatch):
    _nowhere(monkeypatch, tmp_path)
    code = splbl.chosen({"lbl": {"environment": "lbl-rapide"}})
    assert code["python"] is None and not code["ok"]
    assert "not on this machine" in splbl.say_chosen(code)
    fix = splbl.how_to_get(code)
    assert "conda create -n lbl-rapide python=3.12" in fix
    assert "--lbl-env current" in fix


def test_a_path_is_a_python_or_an_environments_folder(tmp_path):
    fake = _python(tmp_path / "env" / "bin" / "python")
    for asked in (str(fake), str(tmp_path / "env")):
        python, _ = splbl.interpreter({"lbl": {"environment": asked}})
        assert python == str(fake), asked
    python, what = splbl.interpreter(
        {"lbl": {"environment": str(tmp_path / "nothing")}})
    assert python is None and "not a python" in what


def test_the_probe_names_the_lbl_that_python_imports(tmp_path):
    package = tmp_path / "site" / "lbl"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("__version__ = '9.9.9'\n")
    fake = _python(tmp_path / "env" / "bin" / "python",
                   "PYTHONPATH=%s exec %s \"$@\"\n"
                   % (tmp_path / "site", sys.executable))
    ok, version, where = splbl.probe(str(fake))
    assert ok and version == "9.9.9"
    assert os.path.realpath(where) == os.path.realpath(package)


def test_another_environment_does_not_get_this_ones_path(tmp_path, monkeypatch):
    """The window puts the repository first on PYTHONPATH, and its lbl/ folder
    then shadowed an LBL installed in place in lbl-rapide."""
    package = tmp_path / "site" / "lbl"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("__version__ = '9.9.9'\n")
    fake = _python(tmp_path / "env" / "bin" / "python",
                   "if [ -n \"$PYTHONPATH\" ]; then echo \"PYTHONPATH reached"
                   " it: $PYTHONPATH\" >&2; exit 1; fi\n"
                   "PYTHONPATH=%s exec %s \"$@\"\n"
                   % (tmp_path / "site", sys.executable))
    monkeypatch.setenv("PYTHONPATH", str(tmp_path))
    ok, version, _ = splbl.probe(str(fake))
    assert ok and version == "9.9.9", version
    assert "PYTHONPATH" not in splbl.environment_for(str(fake))
    assert splbl.environment_for(sys.executable)["PYTHONPATH"] == str(tmp_path), (
        "this interpreter's own runs keep it, as they always did")


def test_a_python_without_lbl_is_not_ok(tmp_path):
    fake = _python(tmp_path / "env" / "bin" / "python",
                   "exec %s -S -c 'raise ImportError(\"no module named lbl\")'\n"
                   % sys.executable)
    ok, why, where = splbl.probe(str(fake))
    assert not ok and where is None and "lbl" in why


def test_the_runner_is_written_for_its_python(tmp_path):
    cfg = load_config("config.yaml", instrument="SPIROU")
    runs = [("AFTER", "comment", splbl.runparams(
        cfg, "lbl", "SPIROU", "CADC", ["TOI-2120_PCA2D_2-3v"], "lbl_config.yaml"))]
    strpca = dict(fit="fit.npz", template="t.fits", mask="m.fits",
                  models_dir="models", prefix="TOI-2120_PCA2D_2-3v", run="2-3v")
    path = tmp_path / "run_lbl.py"
    splbl.write_runner(str(path), runs, strpca, "TOI-2120", "2-3v",
                       "lbl_config.yaml", python="/envs/lbl-rapide/bin/python")
    text = path.read_text()
    assert text.splitlines()[0] == "#!/envs/lbl-rapide/bin/python"
    assert "    /envs/lbl-rapide/bin/python %s" % path in text
    # pca2d found by an environment that does not hold it, and only once LBL
    # is imported: the repository's lbl/ folder would otherwise be `lbl`
    root = os.path.dirname(os.path.dirname(os.path.abspath(splbl.__file__)))
    append = text.index("sys.path.append(%r)" % root)
    assert text.index("from lbl.recipes import lbl_wrap") < append
    assert append < text.index("from pca2d.lbltemplate import strpca_from")

    splbl.write_runner(str(path), runs, None, "TOI-2120", "2-3v",
                       "lbl_config.yaml")
    text = path.read_text()
    assert text.splitlines()[0] == "#!/usr/bin/env python"
    assert "sys.path" not in text, "no STRPCA, no pca2d to find"


def test_the_speed_branch_still_divides_in_place(tmp_path):
    """Its numba kernel keeps main's division, so the STRPCA3+ warning must
    stay when that is the LBL that runs, though general.py has no alias left."""
    (tmp_path / "science").mkdir()
    (tmp_path / "core").mkdir()
    (tmp_path / "science" / "general.py").write_text("x = 1\n")
    kernel = tmp_path / "core" / "fastmath.py"
    kernel.write_text("for i in range(nseg):\n    x[i] = 0\n")
    assert not lbltemplate.resproj_divides_in_place(str(tmp_path))
    kernel.write_text("for ikey in range(nproj):\n"
                      "    for i in range(nseg):\n"
                      "        diff_seg[i] /= bn_seg[i]\n")
    assert lbltemplate.resproj_divides_in_place(str(tmp_path))
    (tmp_path / "science" / "general.py").write_text(
        "frac_diff_seg = diff_seg\nfrac_diff_seg /= (b * n)\n")
    kernel.unlink()
    assert lbltemplate.resproj_divides_in_place(str(tmp_path))


def test_the_window_passes_the_old_lbl_only_when_it_is_chosen():
    from pca2d.gui import OPTIONS_LBL, build_command
    cfg = load_config("config.yaml", instrument="SPIROU")
    kind = dict((key, kind) for key, _path, kind in OPTIONS_LBL)["lbl_env"]
    assert kind[0] == cfg["lbl"]["environment"], "the default is listed first"
    assert "current" in kind, "and the old LBL is offered"
    state = {"objects": ["TOI-2120"], "lbl_env": "lbl-rapide"}
    assert "--lbl-env" not in build_command(state, cfg)
    state["lbl_env"] = "current"
    argv = build_command(state, cfg)
    assert argv[argv.index("--lbl-env") + 1] == "current"


def test_the_flag_sets_the_key():
    import argparse

    from pca2d.cli import add_setting_flags, apply_setting_flags
    args = add_setting_flags(argparse.ArgumentParser()).parse_args(
        ["--lbl-env", "current"])
    cfg = apply_setting_flags(load_config("config.yaml", instrument="SPIROU"),
                              args)
    assert cfg["lbl"]["environment"] == "current"
