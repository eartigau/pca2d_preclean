"""A cube a stage is about to open is checked before the first stage runs.

On 2026-09-15 a user emptied their cache, started a run with the cube stage
unticked because it had always been cached before, and got this:

    260915 15:45:05.09 | stage fit: 2.5 s
    FileNotFoundError: cache/cube_tfits_4e793f8df25f

The fit stage had already announced itself and timed itself before anything
told them the cube was gone. There is exactly one thing to do about a missing
cube, which is to build it; the only choice is whether that is said at the top
or found out from a traceback in np.load.
"""

import os
import shutil

import pytest

from pca2d.cli import (CUBE_FILES, check_cubes, cube_ready, cubes_of,
                       missing_cubes)


def cube(path, complete=True):
    os.makedirs(path, exist_ok=True)
    for name in CUBE_FILES if complete else CUBE_FILES[1:]:
        open(os.path.join(path, name), "wb").close()
    return str(path)


def plan_for(cube_path, obj="GJ1", members=None):
    plan = {"cube": cube_path, "config": {"input": {"object": obj}}}
    if members:
        plan["members"] = members
    return plan


def test_a_cube_that_is_there_changes_nothing(tmp_path):
    plan = plan_for(cube(tmp_path / "cube_tfits_abc"))
    assert missing_cubes(plan) == []
    wanted = ["cube", "fit", "figures", "correct", "lbl"]
    assert check_cubes(plan, wanted) == wanted


def test_a_missing_cube_puts_the_cube_stage_back(tmp_path):
    plan = plan_for(str(tmp_path / "cube_tfits_gone"))
    assert [why for _w, _p, why in missing_cubes(plan)] == ["not there"]
    assert check_cubes(plan, ["fit", "correct"]) == ["cube", "fit", "correct"]
    assert check_cubes(plan, ["fit"])[0] == "cube"


def test_a_half_emptied_cube_counts_as_missing(tmp_path):
    """A cache emptied while a run was not looking can leave the folder and
    take data.npy with it, and np.load's error for that is three stages away."""
    path = cube(tmp_path / "cube_tfits_part", complete=False)
    plan = plan_for(path)
    gone = missing_cubes(plan)
    assert len(gone) == 1 and "data.npy" in gone[0][2]
    assert check_cubes(plan, ["fit"]) == ["cube", "fit"]


def test_the_cube_stage_already_asked_for_is_left_alone(tmp_path):
    plan = plan_for(str(tmp_path / "nothing"))
    wanted = ["cube", "fit"]
    assert check_cubes(plan, wanted) == wanted, "no second cube stage"


def test_a_run_that_opens_no_cube_is_not_made_to_build_one(tmp_path):
    """lbl works from the corrected files, so a re-run of it alone must not be
    turned into a few minutes of reading every spectrum."""
    plan = plan_for(str(tmp_path / "nothing"))
    assert check_cubes(plan, ["lbl"]) == ["lbl"]


def test_every_cube_of_a_joint_run_is_checked(tmp_path):
    """A joint run builds one cube per member and then the cube their rows
    share; a member's missing is the same surprise one stage later."""
    members = [{"object": "PROXIMA", "cube": cube(tmp_path / "c_proxima")},
               {"object": "GJ1", "cube": str(tmp_path / "c_gj1_gone")}]
    plan = plan_for(cube(tmp_path / "c_joint"), "PROXIMA+GJ1", members)

    assert [w for w, _p in cubes_of(plan)] == \
        ["PROXIMA, its own", "GJ1, its own", "PROXIMA+GJ1"]
    gone = missing_cubes(plan)
    assert len(gone) == 1 and gone[0][0] == "GJ1, its own"
    assert check_cubes(plan, ["fit"]) == ["cube", "fit"]


def test_what_is_wrong_is_said_before_any_stage_runs(tmp_path, capsys):
    plan = plan_for(str(tmp_path / "cube_tfits_4e793f8df25f"))
    check_cubes(plan, ["fit", "correct"])
    out = capsys.readouterr().out + capsys.readouterr().err
    assert "cube_tfits_4e793f8df25f" in out, "the path, so it can be looked for"
    assert "not there" in out
    assert "fit, correct" in out, "and which stages were going to open it"


# ---- and everything else a stage opens ----------------------------------
from pca2d.cli import FIT_FILES, check_inputs, missing_fit   # noqa: E402


def full_plan(tmp_path, cube_ok=True, fit_ok=True, fitdir=None):
    out = tmp_path / "out" / "0-3"
    out.mkdir(parents=True, exist_ok=True)
    if fit_ok:
        for name in FIT_FILES.values():
            open(out / name, "wb").close()
    path = cube(tmp_path / "cube_tfits_abc") if cube_ok \
        else str(tmp_path / "cube_tfits_gone")
    plan = plan_for(path)
    plan["outdir"] = str(out)
    if fitdir:
        plan["fitdir"] = fitdir
    return plan


def test_correcting_without_a_fit_runs_the_fit_rather_than_crashing(tmp_path):
    """reconstruct.py raises a bare FileNotFoundError on a missing
    twoframe_components.fits: the cube trap, one stage further along."""
    plan = full_plan(tmp_path, fit_ok=False)
    gone = missing_fit(plan, ["correct"])
    assert [g[1] for g in gone] == ["twoframe_components.fits"]
    assert check_inputs(plan, ["correct"]) == ["fit", "correct"]


def test_drawing_without_a_fit_runs_the_fit_too(tmp_path):
    """A missing fit.npz does not crash the bundle, it quietly puts sequence.py
    on the report's "what did not build" page: a report with no river plot in
    it, which is exactly what went wrong on SMETHELLS_20."""
    plan = full_plan(tmp_path, fit_ok=False)
    assert [g[1] for g in missing_fit(plan, ["figures"])] == ["fit.npz"]
    assert check_inputs(plan, ["figures"]) == ["fit", "figures"]


def test_a_fit_put_back_gets_its_cube_checked_too(tmp_path):
    """Backwards along the chain, or the fit that was just added would open a
    cube nobody looked for."""
    plan = full_plan(tmp_path, cube_ok=False, fit_ok=False)
    assert check_inputs(plan, ["correct"]) == ["cube", "fit", "correct"]


def test_a_fit_that_is_there_changes_nothing(tmp_path):
    plan = full_plan(tmp_path)
    assert missing_fit(plan, ["figures", "correct"]) == []
    assert check_inputs(plan, ["figures", "correct"]) == ["figures", "correct"]


def test_the_fit_stage_already_asked_for_is_left_alone(tmp_path):
    plan = full_plan(tmp_path, fit_ok=False)
    assert check_inputs(plan, ["fit", "correct"]) == ["fit", "correct"], \
        "no second fit stage"


def test_a_variant_reusing_another_fit_is_left_to_check_reused_fit(tmp_path):
    """Its fit is somewhere else, and putting a fit stage back here would make
    a variant rebuild the very fit it exists to reuse."""
    plan = full_plan(tmp_path, fit_ok=False, fitdir=str(tmp_path / "base"))
    assert missing_fit(plan, ["correct"]) == []
    assert check_inputs(plan, ["correct"]) == ["correct"]


def test_lbl_alone_is_never_turned_into_a_fit(tmp_path):
    """The one thing lbl reads from an earlier stage is the corrected spectra,
    and lbl.prepare already says in words that it found none."""
    plan = full_plan(tmp_path, cube_ok=False, fit_ok=False)
    assert check_inputs(plan, ["lbl"]) == ["lbl"]


# ---- reusing something from disk means checking it first ----------------
def test_one_place_decides_whether_a_cube_on_disk_is_a_cube():
    """Every caller about to reuse a cube asks cube_ready, and no caller asks
    os.path.isdir on its own: "is the folder there" is a different question,
    and answering it instead is how a half-emptied cache reports a cache hit
    and then dies in np.load."""
    import inspect

    from pca2d import cli
    for func in (cli.run_cube, cli.run_joint_cube):
        body = inspect.getsource(func)
        assert "cube_ready" in body, "%s reuses without asking" % func.__name__


def test_a_half_emptied_cube_is_rebuilt_not_reused(tmp_path, monkeypatch):
    from pca2d import cli
    from pca2d.cli import cube_ready

    path = cube(tmp_path / "cube_tfits_part", complete=False)
    assert cube_ready(path) and "data.npy" in cube_ready(path)
    assert cube_ready(cube(tmp_path / "cube_tfits_whole")) is None

    built = []
    monkeypatch.setattr("pca2d.build.main", lambda argv: built.append(argv))
    plan = {"cube": path, "files": ["a.fits"],
            "written_config": str(tmp_path / "c.yaml"),
            "config": {"output": {"use_cache": True},
                       "domain": {"wave_min": 965.0, "wave_max": 1950.0,
                                  "dv": 0.5}}}
    cli.run_cube(plan)
    assert built, "an incomplete cube is built again, not handed to the fit"


def test_a_complete_cube_is_still_reused(tmp_path, monkeypatch):
    from pca2d import cli

    built = []
    monkeypatch.setattr("pca2d.build.main", lambda argv: built.append(argv))
    plan = {"cube": cube(tmp_path / "cube_tfits_whole"), "files": ["a.fits"],
            "written_config": str(tmp_path / "c.yaml"),
            "config": {"output": {"use_cache": True},
                       "domain": {"wave_min": 965.0, "wave_max": 1950.0,
                                  "dv": 0.5}}}
    cli.run_cube(plan)
    assert built == [], "the whole point of the cache"


def test_the_fit_stage_itself_says_a_missing_cube_in_words(tmp_path):
    """The other half of the same trap: twoframe is also run on its own with
    --cube, and then check_cubes has not looked at anything. What used to come
    out of np.load was a bare FileNotFoundError, three frames in."""
    from pca2d.twoframe import cube_shape

    with pytest.raises(SystemExit) as gone:
        cube_shape(str(tmp_path / "cube_tfits_4e793f8df25f"))
    assert "cube_tfits_4e793f8df25f" in str(gone.value)
    assert "--stages cube,fit" in str(gone.value) and "window" in str(gone.value)

    half = tmp_path / "cube_tfits_part"
    half.mkdir()
    with pytest.raises(SystemExit) as empty:
        cube_shape(str(half))
    assert "data.npy" in str(empty.value)


def test_a_cache_emptied_mid_run_is_built_again_not_raised(tmp_path, monkeypatch):
    """The other half of the trap, and the one that cost 20 minutes on
    2026-09-15: check_inputs looked before anything ran, and everything was
    there. A joint run then built one member's cube for ten minutes while the
    cleanup page emptied the cache, and the joint build opened a grid.npy that
    had been there when the loop passed it."""
    from pca2d import cli

    built = []

    def fake_build(argv):
        # what the real one does, as far as this test is concerned: it writes
        # the cube its config names, and takes long enough for anything to
        # happen meanwhile
        import yaml
        config = yaml.safe_load(open(argv[0]))
        name = config["input"]["object"]
        built.append(name)
        cube(str(tmp_path / ("cube_" + name)))
        if name == "GL48":                 # the purge, while the run is going
            shutil.rmtree(str(tmp_path / "cube_GL205"))

    plan = {"config": {"output": {"use_cache": True, "reuse_cache": True},
                       "input": {"object": "GL205+GL48"}},
            "cube": str(tmp_path / "cube_joint"),
            "outdir": str(tmp_path), "written_config": str(tmp_path / "c.yaml"),
            "members": [{"object": name, "files": [],
                         "cube": str(tmp_path / ("cube_" + name)),
                         "config": {"input": {"object": name}}}
                        for name in ("GL205", "GL48")]}
    cube(str(tmp_path / "cube_GL205"))     # already built by an earlier run

    joined = {}
    monkeypatch.setattr("pca2d.joint.build",
                        lambda cubes, objects, path, config=None: joined.update(
                            cubes=list(cubes), objects=list(objects)))
    cli.run_joint_cube(plan, fake_build)

    assert built == ["GL48", "GL205"], \
        "GL48 was missing and built; GL205 was there, then was not, and was" \
        " built again rather than opened"
    assert joined["objects"] == ["GL205", "GL48"]
    assert all(cube_ready(c) is None for c in joined["cubes"])


def test_the_stage_about_to_read_a_cube_builds_it_rather_than_dying(tmp_path,
                                                                    monkeypatch):
    """Between two stages a cache can be emptied, and the second of them would
    open what the first one left. Whatever deleted it, there is one thing to do
    about a missing cube."""
    from pca2d import cli

    plan = {"config": {"output": {"use_cache": True, "reuse_cache": True},
                       "input": {"object": "TOI2120"}},
            "cube": str(tmp_path / "cube_tfits_gone")}
    calls = []
    monkeypatch.setattr(cli, "run_cube",
                        lambda p: calls.append(p) or cube(p["cube"]))
    cli.rebuild_missing_cubes(plan, "stage fit")
    assert len(calls) == 1 and cube_ready(plan["cube"]) is None

    cli.rebuild_missing_cubes(plan, "stage figures")
    assert len(calls) == 1, "a cube that is there is not built a second time"


def test_a_rebuild_that_produces_nothing_stops_with_a_reason(tmp_path,
                                                             monkeypatch):
    from pca2d import cli

    plan = {"config": {"output": {"use_cache": False, "reuse_cache": True},
                       "input": {"object": "TOI2120"}},
            "cube": str(tmp_path / "cube_tfits_gone")}
    monkeypatch.setattr(cli, "run_cube", lambda p: None)
    said = []
    monkeypatch.setattr(cli, "log", lambda text, level="info": said.append(
        (level, text)))
    with pytest.raises(SystemExit):
        cli.rebuild_missing_cubes(plan, "stage fit")
    assert any(level == "error" and "use_cache" in text for level, text in said)


def test_the_joint_build_names_the_object_whose_cube_is_gone(tmp_path):
    """Defence in depth: cli builds what is missing before coming here, so
    reaching this means it could not, and np.load's own error names a grid.npy
    and no object at all."""
    from pca2d.joint import build

    with pytest.raises(SystemExit) as gone:
        build([str(tmp_path / "cube_a"), str(tmp_path / "cube_b")],
              ["GL205", "GL48"], str(tmp_path / "joint"))
    assert "GL205" in str(gone.value) and "GL48" in str(gone.value)
    assert "cube stage" in str(gone.value)
