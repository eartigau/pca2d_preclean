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

import pytest

from pca2d.cli import CUBE_FILES, check_cubes, cubes_of, missing_cubes


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
