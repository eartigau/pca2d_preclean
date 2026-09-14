"""The window's own logic: what it finds, what it would run, what it exports.

The widgets are not tested here, since a test machine has no screen; what is
tested is everything the window decides before it draws anything.
"""
import os

from pca2d.gui import (ALL_OPTIONS, OPTIONS, OPTIONS_LBL, build_command,
                       instruments_of, objects_in, variant_yaml)


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


def test_the_command_carries_the_counts_and_the_stages():
    argv = build_command({"objects": ["TOI4552"], "n_star": 0, "n_earth": 3,
                          "stage_cube": True, "stage_fit": True,
                          "dry_run": True})
    assert "--n-star" in argv and argv[argv.index("--n-star") + 1] == "0"
    assert argv[argv.index("--stages") + 1] == "cube,fit"
    assert argv[-1] == "--dry-run"
    every = build_command({"objects": ["TOI4552"],
                           **{"stage_" + s: True for s in
                              ("cube", "fit", "figures", "correct", "lbl")}})
    assert "--stages" not in every, "every stage is the default, so it is not said"


def test_the_window_does_not_run_a_variant_any_more():
    """The parameters converged, so a second set of settings offered beside the
    settings is a window that contradicts itself. `--variant` stays on the
    command line, which is what variants/README.md reproduces from."""
    assert "--variant" not in build_command({"objects": ["X"], "variant": "k0"})


def test_the_export_writes_what_was_changed_and_nothing_else():
    defaults = {"twoframe": {"n_star": 1, "n_earth": 3},
                "correct": {"shrink": True}, "domain": {"dv": 0.5}}
    same = variant_yaml({"n_star": "1", "n_earth": "3", "shrink": True,
                         "dv": "0.5"}, defaults)
    assert same == {}, "a variant that repeats the nominal says nothing"
    changed = variant_yaml({"n_star": "0", "shrink": False, "dv": "0.5"},
                           defaults)
    assert changed == {"twoframe": {"n_star": 0}, "correct": {"shrink": False}}
    settled = variant_yaml({"mean": "offset", "mask": "exposure"}, defaults)
    assert settled == {}, \
        "what the window no longer offers cannot be exported from it either"


def test_the_lbl_window_exports_its_own_block(tmp_path):
    """The LBL settings are settings like the others: what differs is written."""
    defaults = {"lbl": {"run": False, "star_template": False, "suffix": "_P_{tag}",
                        "steps": ["template", "mask", "compute", "compile"],
                        "template": None, "link": "symlink"}}
    same = variant_yaml({"lbl_star_template": False, "lbl_suffix": "_P_{tag}",
                         "lbl_steps": "template, mask, compute, compile",
                         "lbl_link": "symlink", "lbl_template": ""}, defaults)
    assert same == {}, "a window that repeats the configuration says nothing"
    changed = variant_yaml({"lbl_star_template": True, "lbl_steps": "compute,compile",
                            "lbl_link": "copy"}, defaults)
    assert changed == {"lbl": {"star_template": True, "link": "copy",
                               "steps": ["compute", "compile"]}}


def test_two_instruments_are_seen_before_a_run_is_started():
    rows = {"PROXIMA": {"instrument": "NIRPS"}, "GJ1": {"instrument": "NIRPS"},
            "TOI1452": {"instrument": "SPIROU"}, "NEW": {"instrument": "?"}}
    assert instruments_of(rows, ["PROXIMA", "GJ1"]) == {"NIRPS"}
    assert instruments_of(rows, ["PROXIMA", "TOI1452"]) == {"NIRPS", "SPIROU"}
    assert instruments_of(rows, ["PROXIMA", "NEW"]) == {"NIRPS"}, \
        "an instrument that could not be read is not evidence of a second one"
    assert instruments_of(rows, ["NOBODY"]) == set()


def test_every_option_names_a_real_configuration_key():
    from pca2d.config import DEFAULTS
    for _key, path, _kind in ALL_OPTIONS:
        section, name = path.split(".")
        assert section in DEFAULTS, path
        assert name in DEFAULTS[section], path
    assert len({key for key, _p, _k in ALL_OPTIONS}) == len(ALL_OPTIONS), \
        "one widget per setting: two would disagree"
    assert len({path for _k, path, _kind in ALL_OPTIONS}) == len(ALL_OPTIONS)


def test_every_item_explains_itself_in_both_languages():
    """The window is for people meeting the pipeline: every item says what the
    choice implies, in whichever of the two languages is on."""
    from pca2d.gui import EN, FR, STAGES, text
    keys = ["help_data_dir", "help_config", "help_out_dir", "help_objects",
            "help_command", "help_log", "help_rescan", "help_quit_button",
            "help_lang", "help_run_button", "help_stop_button",
            "help_dry_button", "help_export_button", "help_savelog_button",
            "help_openout_button", "help_savedefaults_button",
            "help_lblwin_button", "help_all_button", "help_check",
            "help_col_snr", "help_col_exptime", "help_col_mag",
            "help_berv"]
    keys += ["help_" + key for key, _p, _k in ALL_OPTIONS]
    keys += ["help_stage_" + stage for stage in STAGES]
    for key in keys:
        for name, table in (("en", EN), ("fr", FR)):
            assert key in table, (name, key)
            assert len(table[key]) > 40, ("too short to explain anything",
                                          name, key)
    for key, _path, _kind in ALL_OPTIONS:      # the labels, which are short
        for name, table in (("en", EN), ("fr", FR)):
            assert table.get("opt_" + key), (name, key)
    assert set(EN) == set(FR), "the two languages say the same things"
    assert text("fr", "run") == "Lancer" and text("en", "run") == "Run"
    assert text("de", "run") == "Run", "an unknown language falls back"
    assert text("en", "nothing at all") == "nothing at all"


def speaking(lang):
    """The window's text machinery without a screen: it needs only a language."""
    from pca2d.gui import App
    window = App.__new__(App)
    window.lang = lang
    return lambda key, *args: App._line(window, key, *args)


def test_what_the_window_says_is_timestamped_and_in_its_language():
    """`YYMMDD HH:MM:SS.SS | message`, the convention every log here follows."""
    import re
    english, french = speaking("en"), speaking("fr")
    line = english("log_scan_done", "science_ab12.json", 3, 10, 0, 1)
    assert re.match(r"^\d{6} \d\d:\d\d:\d\d\.\d\d \| ", line), line
    assert line.endswith("\n") and "science_ab12.json" in line
    assert "3 objects, 10 spectra read" in line
    assert "3 objets, 10 spectres lus" in french("log_scan_done", "x", 3, 10, 0, 1)


def test_a_message_whose_placeholders_drifted_does_not_stop_the_window():
    line = speaking("en")("log_scan_done", "only one argument")
    assert "only one argument" in line, "said badly rather than not at all"


def test_a_path_is_shown_as_it_will_be_read(tmp_path):
    """`data` and `config.yaml` say nothing about WHERE: the same two words mean
    a different folder from a different working directory."""
    import os

    from pca2d.gui import absolute
    assert absolute("data") == os.path.join(os.getcwd(), "data")
    assert absolute("~") == os.path.expanduser("~")
    assert absolute("/Volumes/irrisor/x") == "/Volumes/irrisor/x"
    assert absolute("") == "" and absolute(None) == ""


def test_the_output_is_proposed_beside_the_data():
    """Corrected spectra are a copy of the campaign: they belong on the disk the
    campaign is already on, not on whatever disk the window started from."""
    from pca2d.gui import corrected_dir
    assert corrected_dir("/Volumes/irrisor/pca2d_preclean/science") == \
        "/Volumes/irrisor/pca2d_preclean/corrected"
    assert corrected_dir("/Volumes/irrisor/data/") == "/Volumes/irrisor/corrected"
    assert corrected_dir("/Volumes/irrisor/corrected") == \
        "/Volumes/irrisor/corrected", "not nested inside itself"
    assert corrected_dir("") == ""


def test_an_instrument_not_read_yet_is_not_painted_as_a_third_one():
    """A tint says "this instrument". A row whose instrument the scan has not
    reached must not look like a third one, and it did: the fallback palette
    painted it the pale yellow this project keeps for a missing sample."""
    from pca2d.gui import App, INSTRUMENT_TINT, UNREAD

    window = App.__new__(App)
    window._tints = {}

    class Tree:
        def __init__(self):
            self.tags = {}

        def tag_configure(self, tag, background):
            self.tags[tag] = background

    window.tree = Tree()
    assert window._tint("NIRPS") == "inst_NIRPS"
    assert window.tree.tags["inst_NIRPS"] == INSTRUMENT_TINT["NIRPS"]
    assert window._tint("?") == "", "unknown gets no tag at all"
    assert window._tint(UNREAD) == "" and window._tint("") == ""
    assert len(window.tree.tags) == 1, "and no colour was invented for it"


def test_the_column_says_not_read_yet_rather_than_unknown():
    from pca2d.gui import App, UNREAD

    values = App._values(App.__new__(App),
                         {"files": 500, "snr": None, "exptime": None,
                          "mag": None, "instrument": None}, approximate=True)
    assert values[-1] == UNREAD, "the scan has not got there"
    read = App._values(App.__new__(App),
                       {"files": 3, "snr": 100.0, "exptime": 60.0,
                        "mag": None, "instrument": "?"})
    assert read[-1] == "?", "read, and it declared nothing"
    known = App._values(App.__new__(App),
                        {"files": 3, "snr": 100.0, "exptime": 60.0,
                         "mag": None, "instrument": "SPIROU"})
    assert known[-1] == "SPIROU"


def rows_for_sorting():
    return [
        {"object": "PROXIMA", "instrument": "NIRPS", "files": 782, "snr": 182.0,
         "exptime": 201.0, "mag": 5.36},
        {"object": "GJ1", "instrument": "NIRPS", "files": 292, "snr": 163.0,
         "exptime": 178.0, "mag": 5.33},
        {"object": "TOI2120", "instrument": "SPIROU", "files": 321, "snr": 33.0,
         "exptime": 903.0, "mag": 10.45},
        {"object": "GL699_SPIROU", "instrument": None, "files": 535, "snr": None,
         "exptime": None, "mag": None},
    ]


def order(window, rows):
    from pca2d.gui import App
    return [r["object"] for r in App._sorted(window, rows)]


def window_with(column=None, reverse=False):
    from pca2d.gui import App
    w = App.__new__(App)
    w.sort_column, w.sort_reverse = column, reverse
    return w


def test_the_default_order_groups_by_instrument_then_name():
    rows = rows_for_sorting()
    assert order(window_with(), rows) == ["GJ1", "PROXIMA", "TOI2120",
                                          "GL699_SPIROU"], \
        "a run is one instrument, so that is the order the list is for"


def test_a_column_sorts_and_a_second_click_reverses():
    rows = rows_for_sorting()
    assert order(window_with("mag"), rows)[:3] == ["GJ1", "PROXIMA", "TOI2120"]
    assert order(window_with("mag", True), rows)[:3] == ["TOI2120", "PROXIMA",
                                                         "GJ1"]
    assert order(window_with("snr"), rows)[:3] == ["TOI2120", "GJ1", "PROXIMA"]
    assert order(window_with("files"), rows) == ["GJ1", "TOI2120",
                                                 "GL699_SPIROU", "PROXIMA"], \
        "the file count IS known for the one still being scanned"
    assert order(window_with("#0"), rows)[0] == "GJ1", "by name"
    assert order(window_with("#0", True), rows)[0] == "TOI2120"


def test_what_is_not_known_sorts_last_either_way():
    """A blank is not a small number, and a target the scan has not reached
    should not head the list because of it."""
    rows = rows_for_sorting()
    for column in ("mag", "snr", "exptime", "instrument"):
        for reverse in (False, True):
            assert order(window_with(column, reverse), rows)[-1] == \
                "GL699_SPIROU", (column, reverse)


def test_the_command_is_found_beside_the_interpreter_that_runs_the_window():
    """The window is started by an entry point in an environment's bin, and its
    sibling is the command, whatever PATH the window inherited. Started from
    another shell it reported "No such file or directory: 'pca2d-preclean'"
    although the command was installed all along."""
    import os
    import sys

    from pca2d.gui import preclean_argv

    argv = preclean_argv()
    assert argv, "this environment has the package, so something must work"
    if len(argv) == 1:
        assert os.path.isabs(argv[0]) and os.access(argv[0], os.X_OK)
    else:
        assert argv[:2] == [sys.executable, "-m"], argv
        assert argv[2] == "pca2d.cli"


def test_the_shown_command_stays_the_one_to_paste_in_a_terminal():
    """What is run may be an absolute path or `python -m`; what is SHOWN is
    the plain command, since that is what somebody copies elsewhere."""
    from pca2d.gui import build_command
    argv = build_command({"objects": ["GL699_SPIROU"], "n_star": 0})
    assert argv[0] == "pca2d-preclean"


def test_the_run_name_and_the_dates_reach_the_command():
    """Two runs of the same targets at different settings must not write into
    one folder nor under one LBL object name."""
    from pca2d.gui import build_command
    argv = build_command({"objects": ["GL699_SPIROU"], "n_star": 0,
                          "run_name": " saison1 ", "min_rjd": "58383",
                          "max_rjd": "58700"})
    assert argv[argv.index("--name") + 1] == "saison1", "trimmed"
    assert argv[argv.index("--min-rjd") + 1] == "58383"
    assert argv[argv.index("--max-rjd") + 1] == "58700"

    plain = build_command({"objects": ["GL699_SPIROU"], "run_name": "  ",
                           "min_rjd": "", "max_rjd": None})
    assert "--name" not in plain and "--min-rjd" not in plain, \
        "empty is the nominal path, not a run called nothing"


def test_the_command_line_names_a_run_the_same_way_the_window_does():
    """The window shows `--name X`; the run must put it where the window says."""
    import types

    from pca2d.cli import name_run, window_label

    config = {"output": {"directory": "outputs"}, "lbl": {"suffix": "_PCA2D_{tag}"}}
    args = types.SimpleNamespace(name="saison1", min_rjd=None, max_rjd=None)
    assert name_run(config, args) == "saison1"
    assert config["output"]["directory"] == "outputs/_saison1"
    assert config["lbl"]["suffix"] == "_PCA2D_{tag}_saison1"

    # no name: the dates name it themselves
    config = {"output": {"directory": "outputs"}, "lbl": {"suffix": "_PCA2D_{tag}"}}
    args = types.SimpleNamespace(name=None, min_rjd=58383.0, max_rjd=58700.0)
    assert window_label(args) == "rjd58383-58700"
    assert name_run(config, args) == "rjd58383-58700"
    assert config["output"]["directory"] == "outputs/_rjd58383-58700"

    # neither: the nominal path, untouched
    config = {"output": {"directory": "outputs"}, "lbl": {"suffix": "_PCA2D_{tag}"}}
    args = types.SimpleNamespace(name=None, min_rjd=None, max_rjd=None)
    assert name_run(config, args) is None
    assert config["output"]["directory"] == "outputs"

    # a name with a slash in it cannot climb out of the output root
    config = {"output": {"directory": "outputs"}, "lbl": {"suffix": "_P"}}
    args = types.SimpleNamespace(name="../../etc", min_rjd=None, max_rjd=None)
    assert "/" not in name_run(config, args)


def test_a_finished_scan_redraws_the_panels_it_filled():
    """The coverage histogram and the timeline are drawn from what has been
    READ. Nothing redrew them when a scan ended, so the panel a twenty-minute
    run is decided on kept the bars of part of the campaign, and the banner
    that says it is still reading, until a tick changed or the window was
    resized. Seen on a complete index: 8885 spectra read, "UNDER CONSTRUCTION"
    still on the panel."""
    from pca2d.gui import App

    window = App.__new__(App)
    drawn = []
    window.scanning = True
    window._fill = lambda rows: None
    window._state = lambda: None
    window._say = lambda *a, **k: None
    window._draw_berv = lambda: drawn.append("berv")
    window._draw_time = lambda: drawn.append("time")
    index = {"version": 1, "objects": {"GL205": {"files": {}}}}
    window._scanned("root", index, {"read": 0, "kept": 3, "gone": 0}, "i.json")
    assert window.scanning is False
    assert drawn == ["berv", "time"], "both panels, once the scan is over"


def test_a_campaign_finished_mid_scan_redraws_them_only_if_it_is_shown():
    """The bars grow campaign by campaign, and a campaign nobody ticked is in
    neither panel: redrawing for it, every ten spectra, would be hundreds of
    redraws of two canvases for nothing."""
    from pca2d.gui import App

    window = App.__new__(App)
    drawn = []
    window.index = {"version": 1, "objects": {"GL205": {"files": {}}}}
    window.rows = {}
    window.names = {}
    window._values = lambda row, approximate=False: ()
    window._draw_berv = lambda: drawn.append("berv")
    window._draw_time = lambda: drawn.append("time")

    window.picked = lambda: ["GL205"]
    window._one("GL205", 300, 635)               # still reading it
    assert drawn == [], "a partial campaign is redrawn by the next one, not now"
    window._one("GL205", 635, 635)               # done
    assert drawn == ["berv", "time"]

    drawn.clear()
    window.picked = lambda: []                   # ticked by nobody
    window._one("GL205", 635, 635)
    assert drawn == [], "not in either panel, so not redrawn for it"


def test_a_warning_over_the_bars_gets_its_own_ground():
    """"UNDER CONSTRUCTION: still reading" was written straight onto the
    histogram and read against the bars, worst where the campaign is densest.
    The box is measured from the text that was drawn, never guessed from a
    character count, and it goes UNDER it."""
    from pca2d.gui import chip

    class Canvas:
        def __init__(self):
            self.calls = []

        def create_text(self, x, y, **kw):
            self.calls.append(("text", x, y, kw))
            return "t1"

        def bbox(self, item):
            assert item == "t1", "the box is measured from the text itself"
            return (40, 10, 160, 24)

        def create_rectangle(self, *xy, **kw):
            self.calls.append(("rect", xy, kw))
            return "r1"

        def tag_lower(self, below, above):
            self.calls.append(("lower", below, above))

    canvas = Canvas()
    item, box = chip(canvas, 100, 17, "still reading", ink="#b26a00",
                     fill="#fdf1dd", outline="#e0a94a",
                     font=("Helvetica", 10, "bold"), pad=5)
    kinds = [c[0] for c in canvas.calls]
    assert kinds == ["text", "rect", "lower"], "text, then its box, then under"
    assert canvas.calls[1][1] == (35, 7, 165, 27), \
        "the bbox, 5 px either side and 3 above and below"
    assert canvas.calls[1][2]["fill"] == "#fdf1dd"
    assert canvas.calls[2][1:] == ("r1", "t1"), "the box below the text"
    assert (item, box) == ("t1", "r1")


def test_a_fresh_window_opens_on_its_own_config_and_no_roots():
    """The window remembered the last paths, and a first run had none to
    remember: it opened on `data` and `config.yaml` resolved against whatever
    folder it was started from, which on another machine is somebody else's
    layout or nothing at all. A fresh window opens on the config.yaml that came
    with the code that is running, and on NO roots: where the spectra are and
    where the copies go are choices about somebody's disks."""
    from pca2d.gui import installed_config, opening_paths

    fresh = opening_paths({})
    assert fresh["data_dir"] == "", "no data root is invented"
    assert fresh["out_dir"] == "", "and no output root either"
    assert fresh["config"] == installed_config()
    assert os.path.isfile(fresh["config"]), "the one beside the package"
    assert os.path.basename(fresh["config"]) == "config.yaml"

    kept = opening_paths({"data_dir": "~/spectra", "config": "/tmp/mine.yaml",
                          "out_dir": "/tmp/out"})
    assert kept["data_dir"] == os.path.expanduser("~/spectra"), "shown in full"
    assert kept["config"] == "/tmp/mine.yaml" and kept["out_dir"] == "/tmp/out"


def test_the_status_line_says_what_is_missing_rather_than_idle():
    """An empty table with "idle" under it says nothing about what to do."""
    from pca2d.gui import App

    window = App.__new__(App)
    window.lang = "en"
    window.scanning = False
    window.proc = None
    said = []

    class Label:
        def configure(self, text):
            said.append(text)

    window.status = Label()

    class Var:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

    window.vars = {"data_dir": Var("   ")}
    window._state()
    assert said[-1] == window.t("pick_root")
    window.vars["data_dir"] = Var("/data")
    window._state()
    assert said[-1] == window.t("idle")
    window.scanning = True
    window._state()
    assert len(said) == 2, "a scan's progress is left alone"


def test_the_command_box_says_what_is_missing_rather_than_half_a_command():
    """With nothing ticked the box read `pca2d-preclean --config ... --n-star 1`,
    which is not a command anybody can paste, and it made a window with no data
    root look ready to run."""
    from pca2d.gui import command_line

    hint = "tick a target"
    assert command_line({"objects": [], "config": "c.yaml"}, hint) == hint
    assert command_line({"config": "c.yaml"}, hint) == hint
    ready = command_line({"objects": ["PROXIMA"], "config": "c.yaml"}, hint)
    assert ready.startswith("pca2d-preclean --object PROXIMA")


def test_a_reduction_is_named_by_its_targets_and_a_hash_of_the_rest():
    """An empty field says a reduction needs no name, and then two of them land
    in one folder and under one LBL object, where LBL measures the mixture. A
    folder name has to say WHAT was reduced and whether it was reduced like the
    one beside it: the targets, then six characters of a hash of everything
    that makes it a different result."""
    from pca2d.gui import suggested_run_name

    base = {"objects": ["GL205", "GL48"], "n_star": "0", "n_earth": "3",
            "dv": "0.5"}
    name = suggested_run_name(base)
    assert name.startswith("GL205+GL48_"), "the targets are readable"
    assert len(name.rsplit("_", 1)[1]) == 6
    assert suggested_run_name(dict(base)) == name, "the same run, the same name"
    for key, value in (("n_earth", "4"), ("n_star", "1"), ("dv", "0.25"),
                       ("shrink", False), ("min_rjd", "58661"),
                       ("width_kms", "75")):
        assert suggested_run_name(dict(base, **{key: value})) != name, \
            "%s changes the result, so it changes the name" % key
    assert suggested_run_name(dict(base, objects=["GL48", "GL205"])) == name, \
        "the same pair in another order is the same reduction"
    assert suggested_run_name(dict(base, objects=["GL205"])) != name

    many = suggested_run_name({"objects": ["A", "B", "C", "D", "E"]})
    assert many.startswith("A+B+3_"), "past three targets, how many is enough"
    assert suggested_run_name({}).startswith("run_"), "nothing ticked yet"
    assert "/" not in suggested_run_name({"objects": ["a/b"]}), \
        "it becomes a folder name"


def test_the_settings_fill_three_columns_however_many_there_are():
    """Four rows per column was written for eleven settings; two of them have
    since been settled and taken out of the window, which left the ninth alone
    in a column of its own."""
    from pca2d.gui import OPTIONS

    per_column = -(-len(OPTIONS) // 3)
    columns = {}
    for i in range(len(OPTIONS)):
        columns.setdefault(i // per_column, []).append(i)
    assert len(columns) <= 3, "three columns, never a fourth"
    assert max(len(c) for c in columns.values()) - \
        min(len(c) for c in columns.values()) <= 1, "and evenly filled"


def test_a_campaign_being_read_shows_how_much_of_it_is_in():
    """A row whose numbers are still an estimate says so with a tilde, which
    says nothing about how close it is. A slice of a disc beside the name does:
    five steps, and nothing at all once every file is in, so a finished list is
    not a column of symbols."""
    from pca2d.gui import PIE, pie_glyph

    assert pie_glyph(0, 100) == PIE[0], "read nothing yet, an empty disc"
    assert pie_glyph(40, 100) == PIE[2]
    assert pie_glyph(90, 100) == PIE[-1], "nearly there, a full one"
    assert pie_glyph(100, 100) == "", "and gone once it is all in"
    assert pie_glyph(120, 100) == "", "a folder that lost files is not partial"
    assert pie_glyph(3, 0) == "" and pie_glyph(None, None) == ""
    assert pie_glyph("x", 100) == "", "a row that knows nothing shows nothing"
    assert len({pie_glyph(k, 100) for k in (0, 25, 45, 65, 85)}) == 5, \
        "the five steps are five different marks"


def test_quitting_asks_only_while_a_run_would_go_with_it():
    """The settings are written at every change, so leaving loses nothing of the
    window. A run is a subprocess of it and goes when it does, which is worth a
    question; with nothing running, a question would be noise."""
    from pca2d.gui import App

    window = App.__new__(App)
    closed = []
    window._close = lambda: closed.append(True)

    window.proc = None
    window.quit_window(confirm=lambda: (_ for _ in ()).throw(
        AssertionError("nothing is running, so nothing is asked")))
    assert closed == [True]

    window.proc = object()
    window.quit_window(confirm=lambda: False)
    assert closed == [True], "asked, answered no, still here"
    window.quit_window(confirm=lambda: True)
    assert closed == [True, True]


def test_the_disc_of_a_campaign_being_read_has_twelve_steps():
    """Written as a character the disc is the size of the name beside it, so it
    is drawn, and a drawing can have any number of steps. None of them once
    every file is in: a disc that never goes away is decoration."""
    from pca2d.gui import pie_step

    assert pie_step(0, 100) == 0, "nothing read yet is an empty dial"
    assert pie_step(50, 100) == 6
    assert pie_step(99, 100) == 11, "the last step is still not the whole disc"
    assert pie_step(100, 100) is None, "and it goes when they are all in"
    assert pie_step(120, 100) is None
    assert pie_step(3, 0) is None and pie_step(None, None) is None
    assert pie_step("x", 100) is None
    assert len({pie_step(k, 100) for k in range(0, 100, 4)}) == 12, \
        "every step of the twelve is reachable"


def test_a_campaign_copied_in_while_the_window_is_open_is_noticed():
    """Spectra are copied into the data root while the window sits there, and a
    list that only changes when somebody presses Rescan is quietly wrong."""
    from pca2d.gui import folder_news

    shown = {"GL205": 635, "GL48": 989}
    same = folder_news({"GL205": 635, "GL48": 989}, shown,
                       seen={"GL205": 635, "GL48": 989})
    assert same == ([], [], []), "nothing moved, so nothing is read again"

    added, gone, grown = folder_news({"GL205": 635, "GL48": 989, "TOI4552": 119},
                                     shown, seen=None)
    assert (added, gone, grown) == (["TOI4552"], [], []), "a new folder"

    added, gone, grown = folder_news({"GL205": 635}, shown, seen=shown)
    assert (added, gone, grown) == ([], ["GL48"], []), "and one that went"

    added, gone, grown = folder_news({"GL205": 700, "GL48": 989}, shown,
                                     seen=shown)
    assert (added, gone, grown) == ([], [], ["GL205"]), "spectra copied in"


def test_the_count_is_read_against_the_last_look_never_against_the_list():
    """The list counts what the INDEX holds. A spectrum the scan could not read
    is missing from it for good, so comparing the two would ask for a rescan
    every ten seconds for ever."""
    from pca2d.gui import folder_news

    on_disk = {"GL205": 635}
    in_list = {"GL205": 634}                 # one file the scan could not read
    assert folder_news(on_disk, in_list, seen=on_disk) == ([], [], []), \
        "the same disk as a moment ago is nothing to do"
    assert folder_news(on_disk, in_list, seen={"GL205": 600}) == \
        ([], [], ["GL205"]), "but a disk that changed is"
