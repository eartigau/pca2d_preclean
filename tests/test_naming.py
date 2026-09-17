"""What a run is called: its targets, then a hash of the command it is.

The window and the command line must agree, and the hash must move with
everything the command says and with nothing else."""

import types

from pca2d import naming
from pca2d.cli import name_run
from pca2d.gui import build_command, suggested_run_name


def line(**changed):
    argv = ["--object", "GL406", "--config", "/r/config.yaml",
            "--data-dir", "/r/data/NIRPS", "--out-dir", "/r/outputs/NIRPS",
            "--lbl-dir", "/r/lbl/NIRPS", "--n-star", "0", "--n-earth", "7",
            "--weight", "velocity", "--high-pass", "100"]
    for flag, value in changed.items():
        flag = "--" + flag.replace("_", "-")
        argv[argv.index(flag) + 1] = value
    return argv


def test_the_name_moves_with_everything_the_command_says():
    base = naming.run_name(line())
    assert base.startswith("GL406_") and len(base) == len("GL406_") + 6
    assert naming.run_name(line()) == base, "the same command, the same name"
    for flag, value in (("data_dir", "/r/data/SPIROU"),
                        ("out_dir", "/r/outputs/SPIROU"),
                        ("lbl_dir", "/r/lbl/SPIROU"), ("n_earth", "3"),
                        ("weight", "flux"), ("high_pass", "50")):
        assert naming.run_name(line(**{flag: value})) != base, flag
    # the instrument is in the command through its directories, so the same
    # target on two instruments is two runs
    assert naming.run_name(line(data_dir="/r/data/SPIROU")).startswith("GL406_")


def test_the_name_ignores_itself_and_the_targets_order():
    base = naming.run_name(line())
    for extra in (["--name", "auto"], ["--name=auto"], ["--name", "anything"]):
        assert naming.run_name(line() + extra) == base, extra
    joint = ["--objects", "B,A", "--n-star", "0"]
    assert naming.run_name(joint) == naming.run_name(["--objects", "A,B",
                                                      "--n-star", "0"])
    assert naming.run_name(joint).startswith("A+B_")
    assert naming.run_name(["--objects", "A,B,C,D,E"]).startswith("A+B+3_")
    assert naming.run_name([]).startswith("run_")


def test_auto_names_the_run_on_the_command_line(monkeypatch):
    """`--name auto` is the same machinery, fed sys.argv."""
    argv = line() + ["--name", "auto"]
    monkeypatch.setattr("sys.argv", ["pca2d-preclean"] + argv)
    config = {"output": {"directory": "outputs"},
              "lbl": {"suffix": "_PCA2D_{tag}"}}
    args = types.SimpleNamespace(name="auto", object="GL406", objects=None,
                                 min_rjd=None, max_rjd=None)
    label = name_run(config, args)
    assert label == naming.run_name(argv, ["GL406"])
    assert config["output"]["directory"] == "outputs/_" + label
    assert config["lbl"]["suffix"] == "_PCA2D_{tag}_" + label, \
        "its own LBL object too, or two scenarios are one"


def test_the_window_and_the_command_line_agree():
    state = {"objects": ["GL205", "GL48"], "n_star": 0, "n_earth": 7,
             "config": "/r/config.yaml", "out_dir": "/r/outputs"}
    name = suggested_run_name(state)
    argv = build_command(dict(state, objects=sorted(state["objects"]),
                              run_name=""))
    assert name == naming.run_name(argv, sorted(state["objects"]))
    assert name.startswith("GL205+GL48_")
