"""What is handed to LBL, and in LBL's own spelling.

None of this reads a spectrum, and none of it needs LBL to be installed except
the two tests that ask LBL to validate what we wrote, which skip when it is
not there. That pairing is deliberate: the shape of the handoff is pinned
everywhere, and where the environment has LBL, LBL itself is the judge of
whether our config and our runparams are things it will accept.
"""

import os
import re

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
    """CADC, not APERO, even for an APERO reduction.

    The pairing picks the class that reads the files, and only LBL's CADC
    classes read the named fibre extensions where a t.fits keeps its
    wavelength solution. Getting this wrong cost a run.
    """
    assert splbl.profile(spirou())[:2] == ("SPIROU", "CADC")
    nirps = load_config("config.yaml", instrument="NIRPS")
    assert splbl.profile(nirps)[:2] == ("NIRPS_HE", "CADC")


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


# --------------------------------------------------------------- the Teff ---
def test_the_teff_is_read_from_the_spectra_by_default(monkeypatch):
    """APERO writes it as OBJTEMP; LBL stops without it and never looks."""
    monkeypatch.setattr(splbl, "teff_from_header", lambda path: (3179.0, "OBJTEMP"))
    teff, whence = splbl.resolve_teff(spirou(), ["2811170t.fits"])
    assert teff == 3179.0 and "OBJTEMP" in whence


def test_a_number_in_the_config_beats_the_header(monkeypatch):
    monkeypatch.setattr(splbl, "teff_from_header", lambda path: (3179.0, "OBJTEMP"))
    teff, whence = splbl.resolve_teff(spirou(teff=3400), ["2811170t.fits"])
    assert teff == 3400.0 and "config" in whence


def test_a_header_without_a_teff_says_so_rather_than_guessing(monkeypatch):
    monkeypatch.setattr(splbl, "teff_from_header", lambda path: (None, None))
    teff, whence = splbl.resolve_teff(spirou(), ["2811170t.fits"])
    assert teff is None and whence == "nowhere"


def test_the_resolved_teff_is_what_lbl_is_handed(monkeypatch):
    monkeypatch.setattr(splbl, "teff_from_header", lambda path: (3179.0, "OBJTEMP"))
    cfg = spirou()
    teff, _ = splbl.resolve_teff(cfg, ["2811170t.fits"])
    params = splbl.runparams(cfg, "lbl", "SPIROU", "APERO",
                             ["TOI2120", "TOI2120_PCA2D_2-7"], "c.yaml", teff)
    assert params["OBJECT_TEFF"] == [3179.0, 3179.0]
    assert splbl.config_document(cfg, "lbl", "SPIROU", "APERO",
                                 teff)["OBJECT_TEFF"] == 3179.0


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


# ------------------------------------- the coefficients a corrected file keeps
def test_the_corrected_header_names_every_component_the_fit_has():
    """Not only the ones divided out.

    correct.n_star is 0 by default, so the star coefficients are exactly what
    stays in the flux, and nothing downstream carries them: an rdb knows only
    what LBL measured. The header is the only place they can be lined up with
    the exposure they belong to.
    """
    from pca2d.reconstruct import coefficient_cards

    model = {"n_star": 2, "n_earth": 7}
    row = {"a1": -0.5, "a2": 0.25, **{"b%d" % (i + 1): 0.1 * i for i in range(7)}}
    cards = coefficient_cards(model, row, k=0, j=7)

    keys = [key for key, _, _ in cards]
    assert keys[:4] == ["PCASTR_N", "PCASTR_D", "PCASTR01", "PCASTR02"]
    assert keys[4:] == ["PCAOBS_N", "PCAOBS_D"] + [
        "PCAOBS0%d" % (i + 1) for i in range(7)]
    assert dict((k, v) for k, v, _ in cards)["PCASTR01"] == -0.5


def test_the_fit_velocity_is_written_from_a_real_table_row(tmp_path):
    """`"vrad_fit" in row` asks a FITS row about its values, not its columns.

    The tests above hand coefficient_cards a dict, where `in` means keys, so
    they passed while every corrected file went out without PCASTR_V.
    """
    from astropy.io import fits
    from astropy.table import Table
    from pca2d.reconstruct import coefficient_cards

    table = Table({"a1": [0.5], "b1": [0.1], "vrad_fit": [12.5]})
    table.write(tmp_path / "coeffs.fits")
    model = {"n_star": 1, "n_earth": 1}
    for row in (fits.getdata(tmp_path / "coeffs.fits", 1)[0], table[0]):
        values = {key: value for key, value, _ in coefficient_cards(model, row, 0, 1)}
        assert values["PCASTR_V"] == 12.5
    without = Table({"a1": [0.5], "b1": [0.1]})[0]
    assert "PCASTR_V" not in {key for key, _, _ in coefficient_cards(model, without, 0, 1)}


def test_the_counts_say_how_many_cards_follow():
    """And are not the same number as how many were divided out."""
    from pca2d.reconstruct import coefficient_cards

    model = {"n_star": 2, "n_earth": 7}
    row = {"a1": 0.0, "a2": 0.0, **{"b%d" % (i + 1): 0.0 for i in range(7)}}
    cards = coefficient_cards(model, row, k=0, j=7)   # the default correction
    values = {key: value for key, value, _ in cards}
    assert values["PCASTR_N"] == 2, "two star amplitudes are written"
    assert values["PCAOBS_N"] == 7
    assert values["PCASTR_D"] == 0, "and none of them was divided out"
    assert values["PCAOBS_D"] == 7
    listed = [key for key in values if key.startswith("PCASTR") and key[-2:].isdigit()]
    assert len(listed) == values["PCASTR_N"], "the count must match the cards"


def test_the_count_never_promises_a_card_that_was_dropped():
    """Past 99 the extra cards are not written, so the count must not name them."""
    from pca2d.reconstruct import coefficient_cards

    model = {"n_star": 1, "n_earth": 120}
    row = {"a1": 0.0}
    row.update({"b%d" % (i + 1): 0.0 for i in range(120)})
    cards = coefficient_cards(model, row, 0, 120)
    values = {key: value for key, value, _ in cards}
    written = sum(1 for key, _, _ in cards
                  if key.startswith("PCAOBS") and key[-2:].isdigit())
    assert values["PCAOBS_N"] == written == 99


def test_no_keyword_is_longer_than_fits_allows():
    """PCASTR001 would be nine characters and astropy would write HIERARCH."""
    from pca2d.reconstruct import coefficient_cards

    model = {"n_star": 12, "n_earth": 40}
    row = {"a%d" % (i + 1): 0.0 for i in range(12)}
    row.update({"b%d" % (i + 1): 0.0 for i in range(40)})
    for key, _, comment in coefficient_cards(model, row, 0, 40):
        assert len(key) <= 8, key
        assert len(comment) <= 47, "%s: comment would be truncated" % key


def test_a_card_says_whether_that_component_was_removed():
    from pca2d.reconstruct import coefficient_cards

    model = {"n_star": 2, "n_earth": 3}
    row = {"a1": 1.0, "a2": 2.0, "b1": 3.0, "b2": 4.0, "b3": 5.0}
    comments = {key: c for key, _, c in coefficient_cards(model, row, k=1, j=2)}
    assert "divided out" in comments["PCASTR01"]
    assert "left in the flux" in comments["PCASTR02"]
    assert "divided out" in comments["PCAOBS02"]
    assert "left in the flux" in comments["PCAOBS03"]


def test_the_header_spells_the_two_frames_one_way_only():
    """STR for the star, OBS for the observer, in every keyword of the family.

    The counts were PCA2NSTA and PCA2NEAR while the amplitudes were PCASTR01
    and PCAOBS01, so one header called the same two blocks by four names.
    """
    import inspect

    from pca2d import reconstruct

    source = inspect.getsource(reconstruct.correct_file)
    written = set(re.findall(r'head\["(PCA[A-Z0-9_]+)"\]', source))
    assert written, "no header keywords found; did correct_file change?"
    model = {"n_star": 2, "n_earth": 3}
    row = {"a1": 0.0, "a2": 0.0, "b1": 0.0, "b2": 0.0, "b3": 0.0}
    written |= {key for key, _, _ in
                reconstruct.coefficient_cards(model, row, 1, 2)}
    for key in written:
        assert len(key) <= 8, key
        assert "EAR" not in key and "STA" not in key, (
            "%s does not use the STR / OBS spelling the amplitude cards use" % key)


# ---------------------------------------- the name `lbl` must mean LBL -------
def test_nothing_in_the_package_edits_sys_path():
    """A package directory on the path turns every module in it into a top-level
    name. pca2d/lbl.py then shadows LBL, and pca2d/io.py the standard io.

    The CALLS, read from the syntax tree, and not the words: run_lbl.py, which
    pca2d/lbl.py writes as text for another process, appends the repository
    to that process's path once LBL is imported there (lbl.write_runner), and
    a string is not this package editing its own path."""
    import ast
    import pathlib

    def edits_path(node):
        return (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr in ("insert", "append", "extend")
                and isinstance(node.func.value, ast.Attribute)
                and node.func.value.attr == "path"
                and isinstance(node.func.value.value, ast.Name)
                and node.func.value.value.id == "sys")

    root = pathlib.Path(__file__).resolve().parent.parent / "pca2d"
    offenders = [str(p.relative_to(root.parent)) for p in root.rglob("*.py")
                 if any(edits_path(node)
                        for node in ast.walk(ast.parse(p.read_text())))]
    assert not offenders, "edits sys.path: %s" % ", ".join(offenders)


def test_lbl_is_still_found_after_the_figures_stage():
    """The order of a full run: the figures stage imports bundle in-process,
    and the lbl stage comes after it. It used to find pca2d/lbl.py instead."""
    pytest.importorskip("lbl", reason="LBL is in environment.yml")
    import pathlib
    import subprocess
    import sys

    code = ("import pca2d.figures.bundle\n"
            "from pca2d import lbl as stage\n"
            "ok, why = stage.available()\n"
            "assert ok, why\n")
    repo = pathlib.Path(__file__).resolve().parent.parent
    result = subprocess.run([sys.executable, "-c", code], cwd=repo,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-600:]


def test_an_appledouble_file_is_not_a_corrected_spectrum(tmp_path):
    """exFAT puts `._name.fits` beside files; it holds no spectrum."""
    (tmp_path / "2811170t_0-7.fits").write_bytes(b"x")
    (tmp_path / "._2811170t_0-7.fits").write_bytes(b"y")
    (tmp_path / "notes.txt").write_text("z")
    found = [os.path.basename(p) for p in splbl.corrected_files(str(tmp_path))]
    assert found == ["2811170t_0-7.fits"]
    assert splbl.corrected_files(str(tmp_path / "absent")) == []
