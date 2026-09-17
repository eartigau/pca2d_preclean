"""The window's own logic: what it finds, what it would run, what it exports.

The widgets are not tested here, since a test machine has no screen; what is
tested is everything the window decides before it draws anything.
"""
import os
import sys

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
            "help_openout_button", "help_openpdf_button",
            "help_savedefaults_button",
            "help_lblwin_button", "help_all_button", "help_check",
            "help_col_snr", "help_col_exptime", "help_col_mag",
            "help_berv", "help_runs", "help_runs_rescan",
            "help_runs_open_pdf", "help_runs_open_folder"]
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


def test_the_command_that_is_spawned_is_this_package():
    """The window used to look for the pca2d-preclean script beside its own
    interpreter, which is right until two of them are installed. It runs its
    own module now, through its own python: the two are one version, and the
    window never depends on a PATH it inherited from a shell."""
    import os
    import sys

    from pca2d.gui import preclean_argv

    argv = preclean_argv()
    assert argv, "this environment has the package, so something must work"
    if len(argv) == 1:                       # the fallback, without the package
        assert os.path.isabs(argv[0]) and os.access(argv[0], os.X_OK)
    else:
        assert argv[0] == sys.executable
        assert argv[-2:] == ["-m", "pca2d.cli"], argv


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


def test_unticking_a_stage_unticks_what_comes_after_it():
    """The stages happen in an order and depend on each other in that order. A
    run that says it will correct spectra it is not fitting is a run that fails
    twenty minutes in."""
    from pca2d.gui import follow_stages

    all_on = {"cube": True, "fit": True, "figures": True, "correct": True,
              "lbl": True}
    off = dict(all_on, fit=False)
    after = follow_stages(off, "fit")
    assert after["correct"] is False and after["lbl"] is False
    assert after["cube"] is True, "what came before is untouched"
    assert after["figures"] is True, "figures is not in the chain"

    off = dict(all_on, cube=False)
    assert follow_stages(off, "cube") == {"cube": False, "fit": False,
                                          "figures": True, "correct": False,
                                          "lbl": False}


def test_there_is_no_lbl_without_the_correction_it_measures():
    from pca2d.gui import follow_stages

    state = {"cube": False, "fit": False, "figures": False, "correct": False,
             "lbl": True}
    assert follow_stages(state, "lbl")["correct"] is True
    assert follow_stages(state, "lbl")["fit"] is False, \
        "the fit is another matter: it is always redone when it is asked for"

    state = {"cube": False, "fit": False, "figures": False, "correct": True,
             "lbl": False}
    assert follow_stages(state, "correct") == state, \
        "correcting again with the fit that is there is a thing to want"


def test_a_disk_that_is_not_there_is_not_offered(tmp_path):
    """output.fits_directory is a path on ONE machine. A configuration that
    carries one hands every fresh install a disk it has never heard of, and a
    run that finds it missing stops. The window empties the field instead, and
    says so."""
    from pca2d.gui import opening_paths

    here = str(tmp_path)
    gone = str(tmp_path / "not-mounted")

    assert opening_paths({}, gone)["fits_dir"] == "", "blanked, whatever named it"
    assert opening_paths({"fits_dir": gone}, here)["fits_dir"] == "", \
        "including one this window itself remembered"
    assert opening_paths({}, here)["fits_dir"] == here, "a real one stands"
    assert opening_paths({"fits_dir": here}, gone)["fits_dir"] == here, \
        "and what the window kept beats what the configuration says"
    assert opening_paths({})["fits_dir"] == "", "nothing anywhere, nothing shown"


def test_an_empty_products_field_is_a_decision_the_command_carries():
    """`--fits-dir ""` is not a command anybody can paste: a shell drops the
    empty word and argparse then asks for the argument it was promised."""
    from pca2d.gui import build_command

    empty = build_command({"objects": ["X"], "fits_dir": ""})
    assert "--no-fits-dir" in empty and "--fits-dir" not in empty
    assert "" not in empty, "no empty word in a command meant to be pasted"

    named = build_command({"objects": ["X"], "fits_dir": "/mnt/disk"})
    assert named[named.index("--fits-dir") + 1] == "/mnt/disk"
    assert "--no-fits-dir" not in named

    silent = build_command({"objects": ["X"]})
    assert "--fits-dir" not in silent and "--no-fits-dir" not in silent, \
        "a window that never had the field says nothing about it"


def test_the_window_runs_the_code_it_is_itself():
    """A window that had just written --no-fits-dir launched a pipeline that
    had never heard of it (2026-09-15): two editable installs of the package,
    and the run starts in the config file's folder, so the copy sitting there
    won. What is spawned is now this interpreter and this package, said twice:
    as a module, and as the first entry of PYTHONPATH."""
    import pca2d
    from pca2d.gui import package_home, preclean_argv

    argv = preclean_argv()
    assert argv[0] == sys.executable, "the python running the window"
    assert argv[-2:] == ["-m", "pca2d.cli"], "its own module, not a script"
    if sys.version_info >= (3, 11):
        assert "-P" in argv, \
            "without it, `-m` puts the working directory ahead of PYTHONPATH" \
            " and a second clone sitting there is what runs"

    home = package_home()
    assert os.path.isdir(os.path.join(home, "pca2d"))
    assert os.path.dirname(os.path.abspath(pca2d.__file__)) == \
        os.path.join(home, "pca2d"), "the package the window imported"


def test_every_setting_the_window_shows_reaches_the_run():
    """The trap this closes: only --n-star and --n-earth used to travel, so the
    high pass, the shrinkage, the sweeps, the grid step, the coadding and the
    velocity term were shown, changed, and then ignored by the run, which used
    the configuration's own values and said nothing."""
    from pca2d.cli import SETTING_FLAGS
    from pca2d.gui import OPTIONS, build_command

    state = {"objects": ["X"], "n_star": 0, "n_earth": 3, "weight": "velocity",
             "width_kms": 100.0, "shrink": True, "iters": 16, "dv": 0.5,
             "nightly_stack": "auto", "velocity_term": False, "mask": "common"}
    argv = build_command(state)
    carried = {flag for flag in argv if flag.startswith("--")}

    flags = {path: flag for flag, path, _k, _h in SETTING_FLAGS}
    for key, path, _kind in OPTIONS:
        if path in ("twoframe.n_star", "twoframe.n_earth"):
            continue                        # those two have always travelled
        assert path in flags, \
            "%s is on the settings page with no flag to travel under" % path
        assert flags[path] in carried, "%s never reaches the run" % key

    assert argv[argv.index("--weight") + 1] == "velocity"
    assert argv[argv.index("--shrink") + 1] == "true", "a bool as the CLI wants"
    assert argv[argv.index("--velocity-term") + 1] == "false"


def test_a_setting_a_run_was_asked_for_lands_in_its_configuration():
    """And the other half: the flag has to change what the run resolves."""
    import types

    from pca2d.cli import apply_setting_flags, parse_args

    args = parse_args(["--object", "X", "--weight", "velocity",
                       "--high-pass", "75", "--shrink", "false", "--dv", "0.25"])
    config = {"correct": {"weight": "flux", "shrink": True},
              "highpass": {"width_kms": 100.0}, "domain": {"dv": 0.5}}
    apply_setting_flags(config, args)
    assert config["correct"] == {"weight": "velocity", "shrink": False}
    assert config["highpass"]["width_kms"] == 75.0
    assert config["domain"]["dv"] == 0.25

    untouched = {"correct": {"weight": "flux"}}
    apply_setting_flags(untouched, parse_args(["--object", "X"]))
    assert untouched == {"correct": {"weight": "flux"}}, \
        "a flag nobody passed changes nothing"


def test_a_two_valued_setting_is_two_buttons_holding_its_own_words():
    """A menu of two words is a menu too many, and a tick box can only name one
    of the two. Two buttons on one variable: choosing one releases the other by
    construction, each says what it IS (F, (dF/dv)^2), and the variable carries
    the config's own word so nothing downstream translates anything back."""
    from pca2d.gui import OPTIONS

    kinds = {key: kind for key, _path, kind in OPTIONS}
    kind = kinds["weight"]
    assert kind[0] == "radio"
    assert [value for value, _caption in kind[1:]] == ["flux", "velocity"]
    assert [caption for _value, caption in kind[1:]] == ["F", "(dF/dv)\u00b2"], \
        "the buttons say the quantity they weigh by, with a real superscript"
    for key, one in kinds.items():
        if isinstance(one, tuple) and one and one[0] == "radio":
            assert len(one) >= 3, "%s: a choice needs at least two buttons" % key
            for pair in one[1:]:
                assert len(pair) == 2, "%s: (value, what the button says)" % key


def test_every_window_setting_is_on_the_report_page():
    """config.WINDOW_SETTINGS is what the report's summary page prints, and the
    window is what sets them. They have to be the same list in both directions:
    a setting added to the window and not here would change a run and never
    appear on the run's own front page, which is how a (dF/dv)^2 run and a flux
    run read identically on 2026-09-15."""
    from pca2d.config import WINDOW_SETTINGS
    from pca2d.gui import window_settings

    # the LBL folder is a path field since 2026-09-16, and still a setting
    window = window_settings()
    listed = [path for path, _what in WINDOW_SETTINGS]
    assert sorted(window) == sorted(listed), (
        "only in the window: %s; only on the page: %s"
        % (sorted(set(window) - set(listed)), sorted(set(listed) - set(window))))
    assert window == listed, "same order, so the page reads like the window"
    assert all(what for _path, what in WINDOW_SETTINGS), "one line each"


def test_a_setting_the_config_does_not_carry_is_said_out_loud():
    from pca2d.config import setting_value

    config = {"correct": {"weight": "velocity"}, "twoframe": {"n_star": 0}}
    assert setting_value(config, "correct.weight") == "velocity"
    assert setting_value(config, "twoframe.n_star") == 0
    assert setting_value(config, "lbl.run") == "(not set)"
    assert setting_value(config, "correct.weight.deeper") == "(not set)"


# ---- the cleanup page ---------------------------------------------------
def test_the_window_has_a_cleanup_page():
    """The disks fill with cubes and nothing in the window ever said so."""
    import inspect
    import re

    import pca2d.gui as gui
    from pca2d.gui import EN, FR
    source = inspect.getsource(gui.App.__init__)
    tabs = re.search(r'for key in \(([^)]*)\)', source).group(1)
    assert '"tab_clean"' in tabs, "the page is in the notebook"
    assert EN["tab_clean"].strip() == "cleanup"
    assert FR["tab_clean"].strip() == "nettoyage"


def test_the_purge_button_asks_first_and_takes_no_for_an_answer(monkeypatch,
                                                                tmp_path):
    from pca2d.gui import App
    from pca2d.housekeeping import survey

    (tmp_path / "cache").mkdir()
    (tmp_path / "cache" / "big.npy").write_bytes(b"\0" * 5000)
    items = survey({"output": {"cache_directory": str(tmp_path / "cache")}},
                   str(tmp_path / "config.yaml"))

    w = App.__new__(App)
    w.lang, w.clean_items = "en", items
    said = []
    w._say = lambda key, *a, **k: said.append(key)
    w.measure_disks = lambda: said.append("measured")
    asked = []

    def refuse(_self, _title, body):
        asked.append(body)
        return False

    # the question is this window's own now, with a broom where macOS drew the
    # application's icon, so it is the method that is stood in for
    monkeypatch.setattr(App, "_ask_delete", refuse)

    App.purge_disks(w)
    assert asked and "5.0 kB" in asked[0], "it says how much before it asks"
    assert "cube cache" in asked[0], "and what it is about to empty"
    assert (tmp_path / "cache" / "big.npy").exists(), "no means no"
    assert "measured" not in said


def test_the_purge_button_empties_only_what_can_go(monkeypatch, tmp_path):
    from pca2d.gui import App
    from pca2d.housekeeping import survey

    for name, size in (("cache", 4000), ("out", 9000)):
        (tmp_path / name).mkdir()
        (tmp_path / name / "f.bin").write_bytes(b"\0" * size)
    items = survey({"output": {"cache_directory": str(tmp_path / "cache"),
                               "directory": str(tmp_path / "out")}},
                   str(tmp_path / "config.yaml"))

    w = App.__new__(App)
    w.lang, w.clean_items = "fr", items
    said = []
    w._say = lambda key, *a, **k: said.append((key, a))
    w.measure_disks = lambda: said.append(("remeasured", ()))
    monkeypatch.setattr(App, "_ask_delete", lambda *_a: True)

    App.purge_disks(w)
    assert not (tmp_path / "cache" / "f.bin").exists()
    assert (tmp_path / "out" / "f.bin").exists(), "the results are results"
    assert (tmp_path / "cache").is_dir(), "the folder the next run expects"
    assert said[0][0] == "clean_freed" and said[0][1][0] == "4.0 kB"
    assert ("remeasured", ()) in said, "the numbers are read again after"


def test_nothing_to_free_says_so_rather_than_asking(monkeypatch, tmp_path):
    from pca2d.gui import App

    w = App.__new__(App)
    w.lang, w.clean_items = "en", []
    said = []
    w._say = lambda key, *a, **k: said.append(key)
    monkeypatch.setattr(App, "_ask_delete",
                        lambda *_a: pytest.fail("it should not ask"))
    App.purge_disks(w)                       # nothing measured yet
    w.clean_items = [{"name": "results", "kind": "results", "bytes": 10,
                      "removable": False, "path": None}]
    App.purge_disks(w)                       # measured, nothing free in it
    assert said == ["clean_none", "clean_nothing"]


def test_every_row_can_be_deleted_by_picking_it_and_the_question_is_priced():
    """The button is not greyed and the results are not refused: what separates
    them from the cache is what it costs to have them back, so that is what the
    question says, in the terms of the most expensive kind picked."""
    from pca2d.gui import App

    w = App.__new__(App)
    w.lang = "en"
    w.clean_items = [
        {"name": "cube cache", "kind": "rebuildable", "bytes": 3_000_000,
         "removable": True, "path": "cache"},
        {"name": "LBL templates", "kind": "expensive", "bytes": 2_000_000,
         "removable": False, "path": "lbl/templates"},
        {"name": "reports", "kind": "results", "bytes": 1_000_000,
         "removable": False, "path": "out"},
    ]
    asked, deleted, said = [], [], []
    w._say = lambda key, *a, **k: said.append(key)
    w._ask_delete = lambda title, question: asked.append(question) or True
    w._delete = lambda going: deleted.append([it["name"] for it in going])
    w.clean_tree = type("T", (), {"selection": staticmethod(lambda: ())})()

    App.purge_selected(w)
    assert said == ["clean_pick"] and not asked, "nothing picked, nothing asked"

    w.clean_tree.selection = staticmethod(lambda: ("0",))
    App.purge_selected(w)
    assert "made again from what stays here" in asked[-1]
    assert deleted[-1] == ["cube cache"]

    w.clean_tree.selection = staticmethod(lambda: ("0", "1"))
    App.purge_selected(w)
    assert "Hours, not minutes" in asked[-1], "the dearest kind sets the price"

    w.clean_tree.selection = staticmethod(lambda: ("1", "2"))
    App.purge_selected(w)
    assert "IS A RESULT" in asked[-1]
    assert deleted[-1] == ["LBL templates", "reports"]
    # biggest first, whatever order they were picked in
    assert asked[-1].index("LBL templates") < asked[-1].index("reports")


def test_every_line_of_the_cleanup_list_explains_itself_in_both_languages(
        tmp_path):
    """The list is 14 folders and the window is bilingual; the explanation is
    the only part that says what deleting one would cost."""
    import os

    from pca2d.gui import EN, FR
    from pca2d.housekeeping import survey

    # a tree of its own: the LBL lines are listed only for a tree that exists,
    # and `lbl` beside wherever the tests run is there in one checkout only
    (tmp_path / "lbl").mkdir()
    items = survey({"lbl": {"directory": str(tmp_path / "lbl")}},
                   str(tmp_path / "config.yaml"), out_root="out",
                   package=os.path.dirname(os.path.abspath(
                       __import__("pca2d").__file__)))
    assert len(items) >= 13
    for item in items:
        key = "clean_" + item["key"]
        assert key in EN and key in FR, "%s says nothing" % item["name"]
        assert len(EN[key]) > 10 and len(FR[key]) > 10


def test_the_window_knows_which_pdf_is_this_run_s(tmp_path):
    """Open compil PDF opens a file, not the folder above it, so it has to
    name the run exactly as cli.resolve does: the name, the object or the
    joint set, and the tag WITH the v of the velocity term."""
    from pca2d.gui import report_pdf, run_folder

    state = {"objects": ["TOI2120"], "n_star": 0, "n_earth": 3,
             "velocity_term": False, "run_name": ""}
    assert run_folder(state, str(tmp_path)) == str(tmp_path / "TOI2120" / "0-3")
    assert report_pdf(state, str(tmp_path)).endswith("0-3/TOI2120_0-3.pdf")

    # the velocity term is another run, and another folder
    assert report_pdf({**state, "velocity_term": True},
                      str(tmp_path)).endswith("0-3v/TOI2120_0-3v.pdf")
    # a named run keeps its own root, joint runs their own subfolder
    assert run_folder({**state, "run_name": "test1"}, str(tmp_path)) == \
        str(tmp_path / "_test1" / "TOI2120" / "0-3")
    assert run_folder({**state, "objects": ["PROXIMA", "GJ1"]},
                      str(tmp_path)) == str(tmp_path / "joint" / "PROXIMA+GJ1"
                                            / "0-3")
    # and a joint run's PDF carries the joint name, as its folder does
    assert report_pdf({**state, "objects": ["PROXIMA", "GJ1"]},
                      str(tmp_path)) == str(tmp_path / "joint" / "PROXIMA+GJ1"
                                            / "0-3" / "PROXIMA+GJ1_0-3.pdf")
    assert run_folder({"objects": []}, str(tmp_path)) is None


def test_the_proposed_name_follows_the_targets_and_a_typed_one_does_not():
    """The name is the folder, and a folder still reading GJ1 over a run of
    TOI2120 is worse than no name. So the proposal follows the ticks; a name
    somebody typed is a decision and stays."""
    from pca2d.gui import App, suggested_run_name

    class Var:
        def __init__(self, value=""):
            self.value = value

        def get(self):
            return self.value

        def set(self, value):
            self.value = value

    w = App.__new__(App)
    picked = ["GJ1"]
    w.vars = {"run_name": Var()}
    w.state = lambda: {"objects": list(picked), "n_star": 0, "n_earth": 3}
    first = suggested_run_name(w.state())
    w._proposed = first
    w.vars["run_name"].set(first)

    picked[:] = ["TOI2120"]
    w._follow_name()
    assert w.vars["run_name"].get() == suggested_run_name(w.state()) != first

    # the settings are in the name's hash, so they move it too
    w.state = lambda: {"objects": list(picked), "n_star": 2, "n_earth": 3}
    w._follow_name()
    assert w.vars["run_name"].get() == suggested_run_name(w.state())

    # typed by hand: left alone from then on, whatever is ticked
    w.vars["run_name"].set("la bonne")
    picked[:] = ["PROXIMA"]
    w._follow_name()
    assert w.vars["run_name"].get() == "la bonne"

    # emptied on purpose: still a decision
    w.vars["run_name"].set("")
    w._follow_name()
    assert w.vars["run_name"].get() == ""

    # Auto puts it back under the window's care
    App._auto_name(w)
    assert w.vars["run_name"].get() == suggested_run_name(w.state())
    picked[:] = ["GJ1", "PROXIMA"]
    w._follow_name()
    assert w.vars["run_name"].get() == suggested_run_name(w.state())


def test_the_report_of_the_run_that_was_launched_is_still_the_report(tmp_path):
    """The folder is named from the settings, and the proposed name follows the
    ticks. A run launched as GL406_3b7955 finished into that folder while the
    window had moved on to SMETHELLS_20_5b5342, so the button computed a path
    nothing had ever written and said there was no report, twenty seconds after
    the run had said where it had put one. What was launched is on disk, and is
    what the button means."""
    from pca2d.gui import App

    written = tmp_path / "_GL406_3b7955" / "SMETHELLS_20" / "0-3"
    written.mkdir(parents=True)
    pdf = written / "SMETHELLS_20_0-3.pdf"
    pdf.write_text("")

    window = App.__new__(App)
    said, opened = [], []
    window.state = lambda: {"objects": ["SMETHELLS_20"], "n_star": 0,
                            "n_earth": 3, "run_name": "SMETHELLS_20_5b5342"}
    window._out_root = lambda: str(tmp_path)
    window._say = lambda key, *a, **k: said.append(key)
    window._open = lambda path: opened.append(path)

    window.launched_report = str(pdf)
    window.open_report()
    assert opened == [str(pdf)], "the run that was launched wrote this one"
    assert said == ["log_report_launched"], "and the window says it is that one"

    # nothing launched from this window: the warning is still the warning
    del opened[:], said[:]
    window.launched_report = None
    window.open_report()
    assert opened == [] and said == ["log_no_report"]

    # the settings still name a report that exists: that one, with no remark
    del opened[:], said[:]
    here = tmp_path / "_SMETHELLS_20_5b5342" / "SMETHELLS_20" / "0-3"
    here.mkdir(parents=True)
    (here / "SMETHELLS_20_0-3.pdf").write_text("")
    window.launched_report = str(pdf)
    window.open_report()
    assert opened == [str(here / "SMETHELLS_20_0-3.pdf")] and said == []


def test_a_bar_is_cut_where_it_redraws_itself_and_not_at_every_update():
    """tqdm writes a bar, a carriage return, the bar again. Read as text,
    Python turns each of those returns into a newline, and LBL's bars (drawn
    into a pipe, unlike this package's own) arrived as a column:
    ' 33%|...| 42/127', ' 35%|...| 45/127', one finished line per update."""
    from pca2d.gui import cut_output

    segments, rest = cut_output(" 33%|x| 42/127\r 35%|xx| 45/127\r")
    assert segments == [(" 33%|x| 42/127", True), (" 35%|xx| 45/127", True)]
    assert rest == "", "a frame is handed on as it is read, not one late"

    # a real line ends with a newline and is not transient
    segments, rest = cut_output("260915 18:01:04.81 | iter 0\nhalf a li")
    assert segments == [("260915 18:01:04.81 | iter 0", False)]
    assert rest == "half a li", "what has no terminator waits for the next chunk"

    # and the chunk boundary falls wherever the pipe filled
    first, rest = cut_output(" 33%|x| 4")
    second, rest2 = cut_output(rest + "2/127\rdone\n")
    assert first == [] and rest == " 33%|x| 4"
    assert second == [(" 33%|x| 42/127", True), ("done", False)] and rest2 == ""


def test_the_reader_marks_what_the_next_line_has_to_replace():
    """The window's own contract: a line handed over with a carriage return in
    front of it replaces the line written last (App._write). So every segment
    that followed a bar frame carries one, and nothing else does."""
    import queue as _queue

    from pca2d.gui import App

    class Stream:
        def __init__(self, chunks):
            self.chunks = list(chunks)

        def read1(self, _n):
            return self.chunks.pop(0) if self.chunks else b""

    class Proc:
        def __init__(self, stream):
            self.stdout = stream

        def wait(self):
            return 0

    window = App.__new__(App)
    window.lines = _queue.Queue()
    window._line = lambda key, *a: "ended"
    window.proc = Proc(Stream([b" 10%|x| 1/10\r 20%|xx| 2/10\r",
                               b" 30%|xxx| 3/10\rwriting the RDB\nnext\n"]))
    window._reader()

    got = []
    while not window.lines.empty():
        item = window.lines.get()
        if isinstance(item, str):
            got.append(item)
    assert got == [" 10%|x| 1/10\n",
                   "\r 20%|xx| 2/10\n",
                   "\r 30%|xxx| 3/10\n",
                   "\rwriting the RDB\n",
                   "next\n"], \
        "each bar frame replaces the one before it; the line after it too"


def test_the_scan_listing_and_the_clean_listing_are_two_methods():
    """Both were called _listed, so the second was the only one there was: the
    queue handed the scan's {name: count} to the housekeeping formatter and a
    launch died in it on `it["bytes"]`, with a name in `it`."""
    from pca2d.gui import App

    window = App.__new__(App)
    filled = []
    window.index = {"version": 1, "objects": {}}
    window._fill = lambda rows: filled.append(rows)
    window._listed({"GL205": 3, "PROXIMA": 7})
    assert [row["object"] for row in filled[0]] == ["GL205", "PROXIMA"]
    assert [row["files"] for row in filled[0]] == [3, 7], \
        "the folder listing is the truth for the count"

    lines = App._as_lines([{"bytes": 2048, "name": "cube"},
                           {"bytes": 4096, "name": "corrected"}])
    assert lines.splitlines()[0].endswith("corrected"), "biggest first"


def test_the_cleaning_question_wears_a_broom_and_not_the_system_icon():
    """macOS draws the application's icon in a messagebox, and this one is a
    python in a conda environment: the question before deleting six gigabytes
    came up under a generic folder. Same answer as the quit question, which
    has carried its own face since it was written: a Toplevel, and the emoji
    drawn where the system icon was."""
    import pytest

    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except Exception:                      # a machine with no screen
        pytest.skip("no display to open a window on")
    root.withdraw()
    from tkinter import ttk

    from pca2d.gui import App

    window = App.__new__(App)
    window.tk, window.ttk, window.root = tk, ttk, root
    window.t = lambda key: {"clean_go": "On y va", "clean_keep": "Garder"}[key]

    seen = {}

    def answer(which):
        """Find the live dialog, read what it shows, press one of its buttons."""
        tries = []
        def act():
            shown = [w for w in root.winfo_children()
                     if isinstance(w, tk.Toplevel)]
            tries.append(1)
            if len(tries) < 100 and not (
                    shown and shown[0].winfo_viewable()
                    and int(shown[0].attributes("-topmost"))):
                # not shown and raised yet: showing it handles events on
                # macOS, and this can come in the middle of it
                root.after(20, act)
                return
            win = shown[0]
            labels, buttons = [], []
            def walk(widget):
                for child in widget.winfo_children():
                    if isinstance(child, ttk.Label):
                        labels.append(child.cget("text"))
                    if isinstance(child, ttk.Button):
                        buttons.append(child)
                    walk(child)
            walk(win)
            seen["labels"] = labels
            seen["buttons"] = [b.cget("text") for b in buttons]
            # kept above the window it interrupts: on macOS it had dropped
            # behind it, holding the grab, and nothing answered clicks
            seen["topmost"] = bool(int(win.attributes("-topmost")))
            buttons[{"keep": 0, "go": 1}[which]].invoke()
        return act

    root.after(50, answer("go"))
    assert window._ask_delete("effacer", "6.0 GB. On y va ?") is True
    assert "\U0001F9F9" in seen["labels"], "the broom, where the icon was"
    assert "6.0 GB. On y va ?" in seen["labels"]
    assert seen["buttons"] == ["Garder", "On y va"], \
        "the answers named, not Yes and No"
    assert seen["topmost"]

    root.after(50, answer("keep"))
    assert window._ask_delete("effacer", "6.0 GB. On y va ?") is False
    root.destroy()


def test_no_tooltip_opens_while_a_question_waits():
    """A tooltip is a new window, and one opening over the cleaning question
    is what could send the question behind the main window on macOS."""
    import pytest

    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except Exception:                      # a machine with no screen
        pytest.skip("no display to open a window on")
    from tkinter import ttk

    from pca2d.gui import Tip

    app = type("A", (), {"t": staticmethod(lambda key: key), "fonts": {}})()
    button = ttk.Button(root, text="purge")
    button.pack()
    tip = Tip(app, button, "help_clean_purge")
    question = tk.Toplevel(root)
    root.update()
    question.grab_set()
    tip.show()
    assert tip.window is None, "held by a question: no tooltip"
    question.grab_release()
    tip.show()
    assert tip.window is not None, "and one again once it is answered"
    tip.leave()
    root.destroy()


def test_the_runs_tab_lists_what_was_run_and_finds_its_pdf(tmp_path):
    """The page reads the runs themselves: a row each, newest first, and the
    picked one's every option under it."""
    import queue

    import pytest

    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except Exception:                      # a machine with no screen
        pytest.skip("no display to open a window on")
    root.withdraw()
    from tkinter import ttk

    from pca2d import runs as runs_module
    from pca2d.gui import EN, App

    import yaml
    where = tmp_path / "GL406" / "0-7"
    (where / "corrected").mkdir(parents=True)
    (where / "resolved_config.yaml").write_text(yaml.safe_dump(
        {"input": {"object": "GL406"}, "twoframe": {"n_star": 0, "n_earth": 7},
         "provenance": {"command": "pca2d-preclean --object GL406",
                        "started": "2026-09-17T06:53:09"}}))
    (where / "GL406_0-7.pdf").write_bytes(b"%PDF-1.4\n")

    window = App.__new__(App)
    window.tk, window.ttk, window.root = tk, ttk, root
    window.t = lambda key: EN.get(key, key)
    window._register = lambda *a, **k: None
    window._tip = lambda *a, **k: None
    window.lines = queue.Queue()
    window._say = lambda *a, **k: None
    opened = []
    window._open = opened.append
    App._build_runs(window, ttk.Frame(root))
    App._runs_found(window, str(tmp_path), runs_module.listed(str(tmp_path)))
    rows = window.runs_tree.get_children()
    assert len(rows) == 1
    values = window.runs_tree.item(rows[0], "values")
    assert values[0] == "2026-09-17 06:53:09" and values[1] == "GL406"
    assert values[2] == "0-7" and values[4] == "yes"
    assert "1 run(s)" in window.runs_totals.cget("text")

    window.runs_tree.selection_set(rows[0])
    App._run_selected(window)
    options = [window.runs_options.item(i, "values")
               for i in window.runs_options.get_children()]
    assert ("twoframe.n_earth", "7") in options
    assert "--object GL406" in window.runs_command.cget("text")
    App.open_run_pdf(window)
    assert opened and opened[0].endswith("GL406_0-7.pdf")
    App.open_run_folder(window)
    assert opened[1] == str(where)
    root.destroy()
