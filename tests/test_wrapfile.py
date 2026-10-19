"""The window's wrap file: LBL's runparams with a PCA2D dict (pca2d.wrapfile).

The window no longer runs the pipeline; it writes an LBL wrap script and
runs it with LBL's python, which corrects the spectra (lbl.recipes.lbl_pca2d)
before measuring them. What is checked: what the PCA2D dict carries, where
the script and the products go, and that the delivered object never gets one.
"""
import os

import numpy as np
import pytest
from astropy.io import fits

CONFIG = """
general:
  input: {directory: data, pattern: "*t.fits"}
  output: {directory: outputs}
  lbl: {steps: [template, mask, compute, compile], run: true, before: true}
instruments:
  NIRPS:
    lbl: {instrument: NIRPS_HE, data_source: CADC}
objects:
  STAR:
    input: {object_header: STAR_HDR}
    quality: {drs_modes: [HE]}
    output: {cache_directory: /elsewhere}
"""


@pytest.fixture
def setup(tmp_path, monkeypatch):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(CONFIG)
    for name in ("STAR", "OTHER"):
        folder = tmp_path / "data" / name
        folder.mkdir(parents=True)
        hdu = fits.PrimaryHDU(np.zeros(2))
        hdu.header["INSTRUME"] = "NIRPS"
        hdu.header["OBJTEMP"] = 3300.0
        hdu.writeto(str(folder / "a1t.fits"))
    from pca2d import lbl as _lbl
    monkeypatch.setattr(_lbl, "chosen", lambda config: dict(
        asked="lbl-rapide", python="/fake/python", what="test", ok=True,
        detail="0.67", where="/fake/lbl"))
    return tmp_path, str(cfg)


def _state(cfg, **more):
    state = dict(objects=["STAR"], config=cfg, n_star="1", n_earth="3",
                 stage_figures=True, stage_lbl=True)
    state.update(more)
    return state


def _exec(text):
    """The runparams dicts the script defines, without running LBL."""
    names = {}
    exec(compile(text.split('if __name__ == "__main__":')[0]
                 .replace("from lbl.recipes import lbl_wrap", ""),
                 "run.py", "exec"), names)
    return names


def test_the_corrected_dict_carries_the_counts_and_the_object_block(setup):
    from pca2d import wrapfile
    _, cfg = setup
    made = wrapfile.write(_state(cfg), None, dry=True)
    after = _exec(made["text"])["RPARAMS"]
    assert after["PCA2D"]["n_star"] == 1 and after["PCA2D"]["n_earth"] == 3
    assert after["PCA2D"]["figures"] is True
    # the object's block travels, less what LBL decides
    assert after["PCA2D"]["quality"] == {"drs_modes": ["HE"]}
    assert "output" not in after["PCA2D"]
    assert "object_header" not in (after["PCA2D"].get("input") or {})
    assert after["INSTRUMENT"] == "NIRPS_HE" and after["OBJECT_TEFF"] == [3300.0]


def test_one_dict_and_lbl_measures_the_delivered_object_itself(setup):
    """No second dict for the delivered spectra since PCA2D compare=True:
    the wrap file is the one anybody would write."""
    from pca2d import wrapfile
    _, cfg = setup
    names = _exec(wrapfile.write(_state(cfg), None, dry=True)["text"])
    assert "BEFORE" not in names and "AFTER" not in names
    assert names["RPARAMS"]["OBJECT_SCIENCE"] == ["STAR"]
    assert names["RPARAMS"]["PCA2D"]["compare"] is True
    names = _exec(wrapfile.write(_state(cfg, lbl_before=False), None,
                                 dry=True)["text"])
    assert names["RPARAMS"]["PCA2D"]["compare"] is False


def test_the_script_and_the_report_are_in_lbl_s_tree(setup):
    from pca2d import wrapfile
    root, cfg = setup
    made = wrapfile.write(_state(cfg, run_name="trial"), None, dry=True)
    data_dir = os.path.join(str(root), "outputs", "lbl")
    assert made["data_dir"] == data_dir
    assert made["script"] == os.path.join(data_dir, "run_STAR_pca2d_trial.py")
    assert made["report_dirs"] == [os.path.join(data_dir, "pca2d", "_trial",
                                                "STAR", "1-3")]
    assert _exec(made["text"])["RPARAMS"]["PCA2D"]["name"] == "trial"
    assert wrapfile.script_path(_state(cfg, run_name="trial")) == made["script"]


def test_several_objects_are_one_joint_correction(setup):
    from pca2d import wrapfile
    _, cfg = setup
    made = wrapfile.write(_state(cfg, objects=["STAR", "OTHER"]), None, dry=True)
    after = _exec(made["text"])["RPARAMS"]
    assert after["PCA2D"]["joint"] is True
    assert after["OBJECT_SCIENCE"] == ["STAR", "OTHER"]
    assert made["script"].endswith("run_STAR+OTHER_pca2d.py")


def test_without_the_lbl_stage_only_the_correction_is_made(setup):
    from pca2d import wrapfile
    _, cfg = setup
    names = _exec(wrapfile.write(_state(cfg, stage_lbl=False), None,
                                 dry=True)["text"])
    assert names["RPARAMS"]["PCA2D"]["compare"] is False
    assert not any(names["RPARAMS"]["RUN_LBL_%s" % s] for s in
                   ("TEMPLATE", "MASK", "COMPUTE", "COMPILE"))


def test_a_dry_run_writes_and_links_nothing(setup):
    from pca2d import wrapfile
    root, cfg = setup
    made = wrapfile.write(_state(cfg), None, dry=True)
    assert not os.path.exists(made["script"])
    assert not os.path.exists(os.path.join(made["data_dir"], "science"))


def test_a_real_one_links_the_spectra_and_writes_the_script(setup):
    from pca2d import wrapfile
    root, cfg = setup
    made = wrapfile.write(_state(cfg), None)
    assert os.path.exists(made["script"])
    assert os.listdir(os.path.join(made["data_dir"], "science", "STAR")) == ["a1t.fits"]
    assert open(made["script"]).read() == made["text"]
