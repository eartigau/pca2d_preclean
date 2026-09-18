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


def test_a_resume_is_the_same_run():
    """What a run resumed with --stages lbl, after it died in LBL, must land
    in the folder it left: the flags that say how much is run this time are
    out of the hash (TOI2120 earth3 on rali, 2026-09-18)."""
    first = line()
    base = naming.run_hash(first)
    for extra in (["--stages", "lbl"], ["--stages=figures,correct,lbl"],
                  ["--dry-run"], ["--rebuild-cube"], ["--clean-cache"],
                  ["--lbl-before", "false"], ["--lbl-before=true"],
                  ["--stages", "lbl", "--lbl-before", "false", "--name", "x"]):
        assert naming.run_hash(first + extra) == base, extra
    # the program's own name is not the run either: window and command line
    assert naming.run_hash(["pca2d-preclean"] + first) == base
    assert naming.run_hash(["/usr/bin/python", "-m"] + first) != base, \
        "a flag of python's own is still a flag"
    # what DOES change the result still changes the name
    assert naming.run_hash(first + ["--lbl-after", "false"]) != base
    assert naming.run_hash(first + ["--weight", "flux"]) != base


def test_a_report_name_carries_the_hash_and_an_old_one_does_not():
    assert naming.report_name("GL406", "0-7", "abc123") == "GL406_0-7_abc123"
    assert naming.report_name("GL406", "0-7") == "GL406_0-7"
    assert naming.label_with_hash("earth3", "abc123") == "earth3_abc123"
    assert naming.label_with_hash("earth3_abc123", "abc123") == "earth3_abc123"


def a_fit(root, label, star="GL406", tag="0-7", run_hash=None, age=0.0):
    """A run folder holding a fit, as the fit stage leaves it."""
    import os
    import time

    import yaml
    folder = os.path.join(root, ("_" + label) if label else "", star, tag)
    os.makedirs(folder, exist_ok=True)
    for name in ("fit.npz", "twoframe_components.fits"):
        with open(os.path.join(folder, name), "wb") as handle:
            handle.write(b"x")
    provenance = {"command": "pca2d-preclean --object %s" % star}
    if run_hash:
        provenance["run_hash"] = run_hash
    with open(os.path.join(folder, "resolved_config.yaml"), "w") as handle:
        yaml.safe_dump({"input": {"object": star},
                        "provenance": provenance}, handle)
    when = time.time() - age
    os.utime(os.path.join(folder, "fit.npz"), (when, when))
    return folder


def test_reuse_fit_reaches_a_named_run(tmp_path):
    """reuse_fit used to look under the unnamed root only, so a fit made as
    `--name nominal` was out of reach, and with every name carrying a hash
    every fit would have been."""
    from pca2d.cli import reused_fit

    root = str(tmp_path)
    named = a_fit(root, "nominal_ab12cd", run_hash="ab12cd")
    assert reused_fit(root, "nominal", "GL406", "0-7") == named
    assert reused_fit(root, "ab12cd", "GL406", "0-7") == named, "by its hash"
    assert reused_fit(root, "nominal_ab12cd", "GL406", "0-7") == named
    assert reused_fit(root, named, "GL406", "0-7") == named, "by its folder"
    # the unnamed nominal, when there is one, is what `nominal` has meant
    unnamed = a_fit(root, "", run_hash="ffffff")
    assert reused_fit(root, "nominal", "GL406", "0-7") == unnamed
    # a scenario made before the hash, and the same one made again after it:
    # the newest fit
    old = a_fit(root, "earth3", age=3600)
    new = a_fit(root, "earth3_0a0b0c", run_hash="0a0b0c")
    assert reused_fit(root, "earth3", "GL406", "0-7") == new
    import os
    os.utime(os.path.join(new, "fit.npz"), (1, 1))
    assert reused_fit(root, "earth3", "GL406", "0-7") == old


def test_two_scenarios_given_the_same_name_do_not_share_a_folder(
        tmp_path, monkeypatch):
    """A folder of that name made by ANOTHER command keeps its own; this run
    takes the name with its hash. One made by this command, or before any
    hash was recorded, is this run."""
    import types

    from pca2d.cli import name_run

    root = str(tmp_path)
    argv = line()
    monkeypatch.setattr("sys.argv", ["pca2d-preclean"] + argv)
    stamp = naming.run_hash(argv)

    def named():
        config = {"output": {"directory": root},
                  "lbl": {"suffix": "_PCA2D_{tag}"}}
        args = types.SimpleNamespace(name="earth3", object="GL406",
                                     objects=None, min_rjd=None, max_rjd=None)
        return name_run(config, args)

    a_fit(root, "earth3", run_hash="999999")          # another command's
    assert named() == "earth3_" + stamp
    import shutil
    shutil.rmtree(tmp_path / "_earth3")
    a_fit(root, "earth3", run_hash=stamp)             # this command's
    assert named() == "earth3"
    shutil.rmtree(tmp_path / "_earth3")
    a_fit(root, "earth3")                             # before any hash
    assert named() == "earth3"


def test_old_runs_are_read_with_a_hash_of_their_own_command(tmp_path):
    from pca2d import runs

    folder = a_fit(str(tmp_path), "earth3")
    with open(folder + "/GL406_0-7.pdf", "wb") as handle:
        handle.write(b"%PDF")
    found = runs.listed(str(tmp_path))[0]
    assert found["hash"] == naming.run_hash(
        ["pca2d-preclean", "--object", "GL406"]), \
        "a run with no recorded hash gets the one its command gives"
    assert found["report"].endswith("GL406_0-7.pdf"), "its old PDF name"
    stamped = a_fit(str(tmp_path), "earth1_abcdef", run_hash="abcdef")
    with open(stamped + "/GL406_0-7_abcdef.pdf", "wb") as handle:
        handle.write(b"%PDF")
    new = [r for r in runs.listed(str(tmp_path)) if r["hash"] == "abcdef"][0]
    assert new["report"].endswith("GL406_0-7_abcdef.pdf")
