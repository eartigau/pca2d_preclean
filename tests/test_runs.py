"""The runs an output root holds, read from what they wrote down."""

import os

import yaml

from pca2d import runs


def a_run(root, folder, object_name, tag, started=None, pdf=True):
    where = os.path.join(root, folder, tag)
    os.makedirs(os.path.join(where, "corrected"), exist_ok=True)
    provenance = {"command": "pca2d-preclean --object %s" % object_name,
                  "commit": "0123456789abcdef", "dirty": False}
    if started:
        provenance["started"] = started
    with open(os.path.join(where, "resolved_config.yaml"), "w") as handle:
        yaml.safe_dump({"input": {"object": object_name},
                        "twoframe": {"n_star": 0, "n_earth": 7},
                        "correct": {"shrink": True},
                        "provenance": provenance}, handle)
    # a run's own folders hold no run, and are not walked into
    with open(os.path.join(where, "corrected", "resolved_config.yaml"),
              "w") as handle:
        handle.write("{}\n")
    if pdf:
        with open(os.path.join(where, "%s_%s.pdf" % (object_name, tag)),
                  "wb") as handle:
            handle.write(b"%PDF-1.4\n")
    return where


def test_every_run_is_found_newest_first(tmp_path):
    root = str(tmp_path)
    old = a_run(root, "GL406", "GL406", "0-7", started="2026-09-01T10:00:00")
    new = a_run(root, os.path.join("_named", "joint", "A+B"), "A+B", "0-3",
                started="2026-09-17T12:00:00", pdf=False)
    found = runs.listed(root)
    assert [r["folder"] for r in found] == [new, old], "newest first"
    assert [r["objects"] for r in found] == ["A+B", "GL406"]
    assert found[1]["report"].endswith("GL406_0-7.pdf")
    assert found[0]["report"] is None, "no PDF written yet"
    assert found[0]["tag"] == "0-3" and found[0]["components"] == "0 + 7"
    assert "pca2d-preclean" in found[0]["command"]
    assert runs.listed(str(tmp_path / "nowhere")) == []


def test_a_run_without_a_start_time_is_dated_by_its_configuration(tmp_path):
    where = a_run(str(tmp_path), "X", "X", "0-3")
    os.utime(os.path.join(where, "resolved_config.yaml"), (1e9, 1e9))
    found = runs.listed(str(tmp_path))[0]
    assert found["started"].startswith("2001-09-")


def test_the_settings_of_a_run_are_the_ones_the_window_can_set(tmp_path):
    a_run(str(tmp_path), "X", "X", "0-3", started="2026-09-17T12:00:00")
    rows = dict(runs.settings(runs.listed(str(tmp_path))[0]))
    assert rows["targets"] == "X" and rows["started"].startswith("2026-09-17")
    assert rows["twoframe.n_earth"] == "7"
    assert rows["correct.shrink"] == "True"
    assert rows["compilation PDF"].endswith(".pdf")
    assert rows["code"].startswith("0123456789")
