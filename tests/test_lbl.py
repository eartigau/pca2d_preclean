"""What is handed to LBL, and in LBL's own spelling.

None of this reads a spectrum, and none of it needs LBL to be installed except
the two tests that ask LBL to validate what we wrote, which skip when it is
not there. That pairing is deliberate: the shape of the handoff is pinned
everywhere, and where the environment has LBL, LBL itself is the judge of
whether our config and our runparams are things it will accept.
"""

import os

import pytest
import yaml

from pca2d import lbl as splbl
from pca2d.config import load_config


def spirou(**over):
    cfg = load_config("config.yaml", instrument="SPIROU")
    cfg["input"]["object"] = "TOI-2120"
    cfg["lbl"].update(over)
    return cfg


# ------------------------------------------------------------ the objects ---
def test_both_objects_go_in_before_and_after():
    cfg = spirou()
    before, after = splbl.object_names(cfg, "TOI-2120", "2-7")
    assert before == "TOI-2120", "the delivered spectra keep the plain name"
    assert after.startswith("TOI-2120") and after != before


def test_the_corrected_object_is_named_for_the_run_that_made_it():
    """Two component counts must not share one LBL science folder.

    LBL globs the folder, so a 2-7 and a 3-5 correction landing in the same one
    would be measured as a single series, and nothing would say so.
    """
    cfg = spirou()
    _, a = splbl.object_names(cfg, "TOI-2120", "2-7")
    _, b = splbl.object_names(cfg, "TOI-2120", "3-5")
    assert a != b
    assert a == "TOI-2120_PCA2D_2-7" and b == "TOI-2120_PCA2D_3-5"


def test_a_suffix_without_the_tag_is_still_allowed():
    cfg = spirou(suffix="_corr")
    assert splbl.object_names(cfg, "TOI-2120", "2-7")[1] == "TOI-2120_corr"


# ----------------------------------------------------------- the profile ----
def test_the_profile_comes_from_the_instrument_block():
    assert splbl.profile(spirou())[:2] == ("SPIROU", "APERO")
    nirps = load_config("config.yaml", instrument="NIRPS")
    assert splbl.profile(nirps)[:2] == ("NIRPS_HE", "APERO")


def test_a_spectrograph_with_no_profile_stops_rather_than_guesses():
    """Choosing this wrongly is not an error, it is other velocities."""
    cfg = spirou(instrument=None, data_source=None)
    cfg["input"]["instrument"] = "MYSPECTROGRAPH"
    with pytest.raises(SystemExit):
        splbl.profile(cfg)


# --------------------------------------------------------- the runparams ----
def required(params):
    return all(key in params for key in
               ("INSTRUMENT", "DATA_DIR", "DATA_SOURCE", "DATA_TYPES",
                "OBJECT_SCIENCE", "OBJECT_COMPARISON", "OBJECT_TEFF"))


def build(cfg, objects=("TOI-2120", "TOI-2120_PCA2D_2-7")):
    return splbl.runparams(cfg, "lbl", "SPIROU", "APERO", list(objects),
                           "lbl_config.yaml")


def test_the_runparams_are_spelled_the_way_lbl_reads_them():
    """lbl_wrap reads runparams['INSTRUMENT'], not ['instrument']."""
    params = build(spirou())
    assert all(key.isupper() for key in params)
    assert required(params)


def test_every_object_is_measured_as_science_and_against_itself():
    params = build(spirou())
    assert params["DATA_TYPES"] == ["SCIENCE", "SCIENCE"]
    assert params["OBJECT_COMPARISON"] == params["OBJECT_SCIENCE"], (
        "each object builds its own template; one shared template would be a"
        " ruler that fits neither")


def test_a_named_template_overrides_that_for_both():
    params = build(spirou(template="TOI-2120"))
    assert params["OBJECT_COMPARISON"] == ["TOI-2120", "TOI-2120"]


def test_the_steps_asked_for_are_the_ones_turned_on():
    params = build(spirou(steps=["compute", "compile"]))
    assert params["RUN_LBL_COMPUTE"] and params["RUN_LBL_COMPILE"]
    assert not params["RUN_LBL_TEMPLATE"] and not params["RUN_LBL_MASK"]


def test_an_unknown_step_stops_the_run():
    with pytest.raises(SystemExit):
        build(spirou(steps=["compute", "reduce"]))


# ------------------------------------------------------------ the staging ---
def test_linking_is_idempotent_and_names_what_is_not_ours(tmp_path):
    source = tmp_path / "spectra"
    source.mkdir()
    files = []
    for i in range(3):
        path = source / ("281117%dt.fits" % i)
        path.write_text("")
        files.append(str(path))
    target = str(tmp_path / "lbl" / "science" / "TOI-2120")

    linked, kept, strangers = splbl.link_spectra(files, target)
    assert (linked, kept, strangers) == (3, 0, [])
    # a second run must move nothing: LBL globs the folder, and what is already
    # right is already right
    linked, kept, strangers = splbl.link_spectra(files, target)
    assert (linked, kept, strangers) == (0, 3, [])

    (tmp_path / "lbl" / "science" / "TOI-2120" / "someone_elses.fits").write_text("")
    _, _, strangers = splbl.link_spectra(files, target)
    assert strangers == ["someone_elses.fits"], (
        "a file LBL would measure with the rest went unmentioned")


def test_the_glob_written_out_finds_the_corrected_files_too():
    """They are <stem>t_<M>-<N>.fits, so '*t.fits' would find none of them."""
    import fnmatch

    document = splbl.config_document(spirou(), "lbl", "SPIROU", "APERO")
    pattern = document["INPUT_FILE"]
    assert fnmatch.fnmatch("2811170t.fits", pattern)
    assert fnmatch.fnmatch("2811170t_2-7.fits", pattern)


# -------------------------------------------------- and LBL's own verdict ---
def test_lbl_reads_the_config_we_write_and_picks_the_right_class(tmp_path):
    pytest.importorskip("lbl", reason="LBL is in environment.yml")
    from lbl.instruments import select

    document = splbl.config_document(spirou(), str(tmp_path), "SPIROU", "APERO")
    path = tmp_path / "lbl_config.yaml"
    path.write_text(yaml.safe_dump(document, sort_keys=False))
    keys = ["INSTRUMENT", "DATA_DIR", "DATA_SOURCE", "DATA_TYPE", "INPUT_FILE",
            "SKIP_DONE"]
    args = select.parse_args(keys, dict(config_file=str(path)), "test",
                             parse=False)
    assert args["INSTRUMENT"] == "SPIROU"
    assert type(select.load_instrument(args, plogger=None)).__name__ == "Spirou"


def test_lbl_accepts_every_key_the_runner_hands_it():
    pytest.importorskip("lbl", reason="LBL is in environment.yml")
    from lbl.resources import lbl_misc

    params = build(spirou())
    for key in ("INSTRUMENT", "DATA_DIR", "DATA_SOURCE", "DATA_TYPES",
                "OBJECT_SCIENCE", "OBJECT_COMPARISON", "OBJECT_TEFF"):
        lbl_misc.check_runparams(params, key)
