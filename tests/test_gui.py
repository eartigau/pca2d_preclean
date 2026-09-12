"""The window's own logic: what it finds, what it would run, what it exports.

The widgets are not tested here, since a test machine has no screen; what is
tested is everything the window decides before it draws anything.
"""
import os

from pca2d.gui import OPTIONS, build_command, objects_in, variant_yaml


def test_the_objects_of_a_data_root_are_the_folders_with_spectra(tmp_path):
    for name, n in (("PROXIMA", 3), ("GJ1", 2), ("EMPTY", 0)):
        folder = tmp_path / name
        folder.mkdir()
        for i in range(n):
            (folder / ("%04dt.fits" % i)).write_text("")
    (tmp_path / "notes.txt").write_text("")
    assert objects_in(str(tmp_path)) == [("GJ1", 2), ("PROXIMA", 3)]
    assert objects_in(str(tmp_path / "nothing")) == []


def test_one_object_and_several_are_different_commands():
    one = build_command({"objects": ["PROXIMA"], "config": "config.yaml"})
    assert one[:5] == ["pca2d-preclean", "--object", "PROXIMA", "--config",
                       "config.yaml"]
    several = build_command({"objects": ["PROXIMA", "GJ1", "GJ3090"]})
    assert several[1:3] == ["--objects", "PROXIMA,GJ1,GJ3090"], \
        "several objects are fitted together against one observer basis"


def test_the_command_carries_the_counts_the_stages_and_the_variant():
    argv = build_command({"objects": ["TOI4552"], "n_star": 0, "n_earth": 3,
                          "variant": "k0", "stage_cube": True, "stage_fit": True,
                          "dry_run": True})
    assert "--n-star" in argv and argv[argv.index("--n-star") + 1] == "0"
    assert argv[argv.index("--stages") + 1] == "cube,fit"
    assert argv[argv.index("--variant") + 1] == "k0"
    assert argv[-1] == "--dry-run"
    every = build_command({"objects": ["TOI4552"],
                           **{"stage_" + s: True for s in
                              ("cube", "fit", "figures", "correct", "lbl")}})
    assert "--stages" not in every, "every stage is the default, so it is not said"
    assert "--variant" not in build_command({"objects": ["X"], "variant": "(none)"})


def test_the_export_writes_what_was_changed_and_nothing_else():
    defaults = {"twoframe": {"n_star": 1, "mean": "star"},
                "correct": {"shrink": True, "mask": "common"},
                "domain": {"dv": 0.5}}
    same = variant_yaml({"n_star": "1", "mean": "star", "shrink": True,
                         "mask": "common", "dv": "0.5"}, defaults)
    assert same == {}, "a variant that repeats the nominal says nothing"
    changed = variant_yaml({"n_star": "0", "mask": "exposure", "dv": "0.5"},
                           defaults)
    assert changed == {"twoframe": {"n_star": 0}, "correct": {"mask": "exposure"}}


def test_every_option_names_a_real_configuration_key():
    from pca2d.config import DEFAULTS
    for _key, path, _kind in OPTIONS:
        section, name = path.split(".")
        assert section in DEFAULTS, path
        assert name in DEFAULTS[section], path


def test_every_item_explains_itself_in_both_languages():
    """The window is for people meeting the pipeline: every item says what the
    choice implies, in whichever of the two languages is on."""
    from pca2d.gui import EN, FR, STAGES, text
    keys = ["help_data_dir", "help_config", "help_out_dir", "help_objects",
            "help_variant", "help_command", "help_log", "help_rescan",
            "help_lang", "help_run_button", "help_stop_button",
            "help_dry_button", "help_export_button", "help_savelog_button",
            "help_openout_button"]
    keys += ["help_" + key for key, _p, _k in OPTIONS]
    keys += ["help_stage_" + stage for stage in STAGES]
    for key in keys:
        for name, table in (("en", EN), ("fr", FR)):
            assert key in table, (name, key)
            assert len(table[key]) > 40, ("too short to explain anything",
                                          name, key)
    for key, _path, _kind in OPTIONS:          # the labels, which are short
        for name, table in (("en", EN), ("fr", FR)):
            assert table.get("opt_" + key), (name, key)
    assert set(EN) == set(FR), "the two languages say the same things"
    assert text("fr", "run") == "Lancer" and text("en", "run") == "Run"
    assert text("de", "run") == "Run", "an unknown language falls back"
    assert text("en", "nothing at all") == "nothing at all"
