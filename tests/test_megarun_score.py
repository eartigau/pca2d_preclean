"""The sweep's table: one row per star per run, with what each run was given
and what it did to the robust sigma. No network, no cube: the velocities are
written as rdb files, as LBL leaves them."""

import os
import sys

import numpy as np
import yaml
from astropy.table import Table

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "megarun"))
import score                                                    # noqa: E402


def rdb(path, t, v, e):
    Table({"rjd": t, "vrad": v, "svrad": e}).write(path, format="ascii.rdb",
                                                   overwrite=True)


def a_run(root, star, tag, scatter_before, scatter_after, n_earth=7,
          weight="velocity"):
    """A run folder as the stages leave it, with LBL's velocities beside it."""
    folder = os.path.join(root, "outputs", star, tag)
    tree = os.path.join(root, "lbl", "lblrdb")
    os.makedirs(folder, exist_ok=True)
    os.makedirs(tree, exist_ok=True)
    config = {"input": {"object": star, "directory": root},
              "twoframe": {"n_star": 0, "n_earth": n_earth,
                           "velocity_term": False},
              "correct": {"weight": weight, "shrink": True},
              "highpass": {"width_kms": 100.0},
              "lbl": {"directory": os.path.join(root, "lbl"),
                      "suffix": "_PCA2D_{tag}"},
              "provenance": {"started": "2026-09-17T12:00:00",
                             "command": "pca2d-preclean --object %s" % star}}
    with open(os.path.join(folder, "resolved_config.yaml"), "w") as handle:
        yaml.safe_dump(config, handle)
    r = np.random.default_rng(3)
    t = np.sort(60000 + r.uniform(0, 300, 60))
    for name, scatter in ((star, scatter_before),
                          ("%s_PCA2D_%s" % (star, tag), scatter_after)):
        rdb(os.path.join(tree, "lbl_%s_%s.rdb" % (name, name)), t,
            1000 + r.normal(0, scatter, t.size), np.full(t.size, 1.0))
    return folder


def test_each_run_is_scored_by_what_it_did_to_the_robust_sigma(tmp_path):
    root = str(tmp_path)
    a_run(root, "GL406", "0-7", 10.0, 5.0)
    a_run(root, "GL725B", "0-3", 10.0, 12.0, n_earth=3, weight="flux")
    rows = score.table([os.path.join(root, "outputs")])
    assert [row["star"] for row in rows] == ["GL406", "GL725B"], "best first"
    good, bad = rows
    assert good["gain_robust"] > 1.5 and bad["gain_robust"] < 1.0
    assert good["n_earth"] == 7 and bad["n_earth"] == 3
    assert good["weight"] == "velocity" and bad["weight"] == "flux"
    assert good["high_pass_kms"] == 100.0 and good["n"] == 60
    assert good["report"] == "", "no PDF written in this test"

    csv_path = os.path.join(root, "results.csv")
    status = os.path.join(root, "status.md")
    score.write_csv(rows, csv_path)
    score.write_status(rows, status, [root], waiting=[
        {"objects": "TOI782", "name": "TOI782", "started": "2026-09-17"}])
    text = open(status).read()
    assert "| GL406 |" in text and "**" in text
    assert "TOI782" in text and "no velocities yet" in text
    assert open(csv_path).readline().startswith("star,run,hash,started")


def test_a_joint_run_is_scored_star_by_star(tmp_path):
    root = str(tmp_path)
    folder = a_run(root, "A", "0-7", 10.0, 6.0)
    # a joint run leaves one cube_config per member beside the run
    tree = os.path.join(root, "lbl", "lblrdb")
    for star in ("A", "B"):
        with open(os.path.join(folder, "cube_config_%s.yaml" % star),
                  "w") as handle:
            yaml.safe_dump({"target": {}}, handle)
        r = np.random.default_rng(5)
        t = np.sort(60000 + r.uniform(0, 300, 60))
        for name, scatter in ((star, 8.0), ("%s_PCA2D_0-7" % star, 4.0)):
            rdb(os.path.join(tree, "lbl_%s_%s.rdb" % (name, name)), t,
                1000 + r.normal(0, scatter, t.size), np.full(t.size, 1.0))
    rows = score.table([os.path.join(root, "outputs")])
    assert sorted(row["star"] for row in rows) == ["A", "B"]
    assert all(row["joint"] for row in rows)
