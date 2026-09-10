"""run_lbl.py, run against a stand-in for LBL that records what it is handed.

The runner is the one place the order of LBL's steps is decided, and for STRPCA
the order is the point: the tables are written in the rest frame the mask step
measures, so the corrected object's mask has to exist before they are written
and they have to exist before its velocities are measured. The delivered
spectra never get them.
"""

import sys
import types

import numpy as np

from pca2d import lbl as splbl
from pca2d import lbltemplate
from pca2d.config import load_config

AFTER = "TOI-2120_PCA2D_2-3v"


def _runner(tmp_path, strpca=None):
    cfg = load_config("config.yaml", instrument="SPIROU")
    runs = [(name, "comment", splbl.runparams(cfg, "lbl", "SPIROU", "CADC", [obj],
                                              "lbl_config.yaml"))
            for name, obj in (("BEFORE", "TOI-2120"), ("AFTER", AFTER))]
    path = tmp_path / "run_lbl.py"
    splbl.write_runner(str(path), runs, strpca, "TOI-2120", "2-3v", "lbl_config.yaml")
    return path.read_text()


def _execute(text, monkeypatch):
    """Run a runner with lbl_wrap.main replaced by a recorder."""
    calls = []
    recipes = types.ModuleType("lbl.recipes")
    recipes.lbl_wrap = types.SimpleNamespace(main=lambda p: calls.append(dict(p)))
    monkeypatch.setitem(sys.modules, "lbl", types.ModuleType("lbl"))
    monkeypatch.setitem(sys.modules, "lbl.recipes", recipes)
    monkeypatch.setattr(lbltemplate, "strpca_from",
                        lambda **kw: {"STRPCA2": "STRPCA2_%s.fits" % kw["prefix"]})
    exec(compile(text, "run_lbl.py", "exec"), {"__name__": "__main__"})
    return calls


STRPCA = dict(fit="fit.npz", template="star_template.fits", mask="mask.fits",
              models_dir="lbl/models", prefix=AFTER, run="2-3v")


def test_one_dict_per_object_and_no_tables_without_a_second_component(tmp_path, monkeypatch):
    calls = _execute(_runner(tmp_path), monkeypatch)
    assert [c["OBJECT_SCIENCE"] for c in calls] == [["TOI-2120"], [AFTER]]
    assert not any("RESPROJ_TABLES" in c for c in calls)
    assert all(c["RUN_LBL_COMPUTE"] for c in calls)


def test_strpca_is_written_between_the_mask_and_the_velocities(tmp_path, monkeypatch):
    calls = _execute(_runner(tmp_path, STRPCA), monkeypatch)
    before, mask, after = calls
    assert before["OBJECT_SCIENCE"] == ["TOI-2120"] and "RESPROJ_TABLES" not in before, (
        "the delivered spectra were measured without the tables")
    assert mask["OBJECT_SCIENCE"] == [AFTER] and mask["RUN_LBL_MASK"]
    assert not mask["RUN_LBL_COMPUTE"] and not mask["RUN_LBL_COMPILE"]
    assert "RESPROJ_TABLES" not in mask
    assert after["RESPROJ_TABLES"] == {"STRPCA2": "STRPCA2_%s.fits" % AFTER}
    assert after["RUN_LBL_COMPUTE"] and after["RUN_LBL_COMPILE"]


def test_the_runner_says_what_each_dict_is(tmp_path):
    text = _runner(tmp_path, STRPCA)
    assert "BEFORE = dict(" in text and "AFTER = dict(" in text
    assert "STRPCA = dict(" in text and "    prefix=%r," % AFTER in text
    assert "INSTRUMENT='SPIROU'" in text, "LBL's keys stay upper case"


def test_the_number_of_star_components_is_read_from_the_fit(tmp_path):
    assert splbl.fit_components(str(tmp_path)) == 0
    np.savez(tmp_path / "fit.npz", P=np.zeros((3, 10)))
    assert splbl.fit_components(str(tmp_path)) == 3
